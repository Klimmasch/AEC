classdef Model < handle
    properties
        scmodel_Large;      %SparseCoding class (downsampled version)
        scmodel_Small;      %SparseCoding class
        rlmodel;            %ReinforcementLearning class

        focalLength;        %focal length [px]
        baseline;           %interocular distance
        objDistMin;         %object distance to eyes [m]
        objDistMax;
        muscleInitMin;      %minimal initial muscle innervation
        muscleInitMax;      %maximal --"--
        interval;           %period of eye stimulus change
        desiredAngleMin;    %min/max desired vergence angle
        desiredAngleMax;

        learnedFile;        %previously learned model
        textureFile;        %config file containing texture stimulus list
        trainTime;          %# training iterations
        simulatedTime;     %how long did the model take to be learned (min)

        sparseCodingType;   %type of sparse coding
        stopSC;
        SCInterval;

        lambdaMuscleFB;     %factor of muscle activity feedback to RL feature vector
        lambdaMet;          %factor of metCosts for reward function
        lambdaRec;          %reconstruction error factor
        lambdaV;            %value networks input->output weights factor
        lambdaP1;           %policy networks input->hidden weights factor
        lambdaP2;           %policy networks hidden->output weights factor

        % Model data history
        recerr_hist;        %reconstruction error [coarse scale, fine scale]
        disp_hist;          %disparity
        vergerr_hist;       %vergence error
        verge_actual;       %actual vergence angle
        verge_desired;      %desired vergence angle
        Z;                  %object depth
        fixZ;               %fixation depth
        g_hist;             %natural gradient change
        td_hist;            %temporal difference (td) error
        feature_hist;       %feature vector
        cmd_hist;           %vergence commands
        relCmd_hist;        %relativ changes in motor commands
        l12_weights;        %L1/L2, i.e. sum abs, sum pow2 weights of actor and critic
        reward_hist;        %reward function
        metCost_hist;       %metabolic costs
        variance_hist;      %exploratory variance of actor
        savePath;           %where all the data are stored

        % Model results at testing procedure
        disZtest;
        fixZtest;
        vergErrTest;
        responseResults;
    end

    methods
        function obj = Model(PARAM)
            obj.learnedFile = PARAM{1}{1};
            obj.textureFile = PARAM{1}{2};
            obj.trainTime = PARAM{1}{3};
            obj.sparseCodingType = PARAM{1}{4};
            obj.focalLength = PARAM{1}{5};
            obj.baseline = PARAM{1}{6};
            obj.objDistMin = PARAM{1}{7};
            obj.objDistMax = PARAM{1}{8};
            obj.muscleInitMin = PARAM{1}{9};
            obj.muscleInitMax = PARAM{1}{10};
            obj.interval = PARAM{1}{11};
            obj.lambdaMuscleFB = PARAM{1}{12};
            obj.lambdaMet = PARAM{1}{13};
            obj.lambdaRec = PARAM{1}{14};
            obj.lambdaV = PARAM{1}{15};
            obj.lambdaP1 = PARAM{1}{16};
            obj.lambdaP2 = PARAM{1}{17};

            % single eye
            obj.desiredAngleMin = atand(obj.baseline / (2 * obj.objDistMax));
            obj.desiredAngleMax = atand(obj.baseline / (2 * obj.objDistMin));

            % Discrete or continuous policy
            if (PARAM{3}{14})
                obj.rlmodel = ReinforcementLearningCont(PARAM{3});
            else
                obj.rlmodel = ReinforcementLearning(PARAM{3});
            end

            % Create 2 sparse coding objects with same parameters
            if (strcmp(obj.sparseCodingType, 'homeo'))
                obj.scmodel_Large = SparseCodingHomeo(PARAM{2}{1}); %coarse scale
                obj.scmodel_Small = SparseCodingHomeo(PARAM{2}{2}); %fine scale
            else
                obj.scmodel_Large = SparseCoding(PARAM{2}{1});      %coarse scale
                obj.scmodel_Small = SparseCoding(PARAM{2}{2});      %fine scale
            end
            obj.stopSC = obj.trainTime;
            obj.SCInterval = 1;

            obj.recerr_hist = zeros(obj.trainTime, 2);
            obj.disp_hist = zeros(obj.trainTime, 1);
            obj.vergerr_hist = zeros(obj.trainTime, 1);
            obj.verge_actual = zeros(obj.trainTime, 1);
            obj.verge_desired = zeros(obj.trainTime, 1);
            obj.Z = zeros(obj.trainTime, 1);
            obj.fixZ = zeros(obj.trainTime, 1);

            obj.g_hist = zeros(obj.trainTime, 1);
            obj.td_hist = zeros(obj.trainTime, 1);
            % obj.feature_hist = zeros(obj.trainTime, obj.rlmodel.inputDim);
            obj.cmd_hist = zeros(obj.trainTime, 2);
            obj.relCmd_hist = zeros(obj.trainTime, 1);
            obj.l12_weights = zeros(obj.trainTime, 4);
            obj.reward_hist = zeros(obj.trainTime, 1);
            obj.metCost_hist = zeros(obj.trainTime, 1);
            obj.variance_hist = zeros(obj.trainTime, 1);

            obj.disZtest = [];
            obj.fixZtest = [];
            obj.vergErrTest = [];
            obj.responseResults = struct();
        end

        %%% Make a (deep) copy of a handle object
        function new = copy(this)
            % Instantiate new object of the same class
            new = feval(class(this));

            % Copy all non-hidden properties
            p = properties(this);
            % Copy hidden properties which don't show up in the result
            % of properties(this)
            % p = fieldnames(struct(this));
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end

        %%% Generate Feature Vector and Reward
        % Feature is concatenation of small and large scale feature vector
        % Reward is sum of separate sc errors (small and large scale)
        % Error is sum of errors
        function [feature, reward, Error, Error_L, Error_S] = generateFR(this, Images)
           %LARGE SCALE (downsampled version)
           imagesLarge = Images{1};

            imPind = find(sum(imagesLarge .^ 2)); %find non-zero patches (columns)
            if (isempty(imPind))
                feature_L = zeros(this.scmodel_Large.Basis_num, 1);
                reward_L = this.rlmodel.J;
                return;
            end
            [Coef_L, Error_L] = this.scmodel_Large.sparseEncode(imagesLarge);

            %cost_L = sum(Error_L.^2);        %Error for each patch
            %Error_L = mean(cost_L(imPind));  %can be changed with new rec error formula  (error relative to original patch)
            Error_L = sum(sum(Error_L .^ 2)) / sum(sum(imagesLarge .^ 2));

            reward_L = -Error_L;
            feature_L = mean(Coef_L(:, imPind) .^ 2, 2); %feature vector for large scale (*5), average activation, Eq. 3.1

            %SMALL SCALE
            imagesSmall = Images{2};

            imPind = find(sum(imagesSmall .^ 2)); %find non-zero patches
            if (isempty(imPind))
                feature_S = zeros(this.scmodel_Small.Basis_num, 1);
                reward_S = this.rlmodel.J;
                return;
            end
            [Coef_S, Error_S] = this.scmodel_Small.sparseEncode(imagesSmall);

            %cost_S = sum(Error_S.^2);          %Error for each patch
            %Error_S = mean(cost_S(imPind));    %can be changed with new rec error formula (error relative to original patch)
            Error_S = sum(sum(Error_S .^ 2)) / sum(sum(imagesSmall .^ 2));

            reward_S = -Error_S;
            feature_S = mean(Coef_S(:, imPind) .^ 2, 2); %feature vector for small scale

            %compute total rec error and reward
            Error = Error_L + Error_S;          %total reconstruction error
            reward = reward_L + reward_S;       %sum rewards
            feature = [feature_L; feature_S];   %join feature vectors
        end

        %% Plotting everything and save graphs
        function allPlotSave(this)
            windowSize = 125;
            if (this.trainTime < windowSize * this.interval)
                windowSize = round(this.trainTime / this.interval / 5);
            end
            % only take the last value before the image/texture is changed
            ind = this.interval : this.interval : this.trainTime;

            %% Simple Moving Average Vergence Error
            vergerr = filter(ones(1, windowSize) / windowSize, 1, abs(this.vergerr_hist(ind)));

            figure;
            hold on;
            grid on;
            % Raw vergance error values
            plot(this.interval : this.interval : size(this.vergerr_hist), abs(this.vergerr_hist(ind)), ...
                 'color', [1, 0.549, 0], 'LineWidth', 0.8);

            % Simple Moving Average
            plot((windowSize + 1) * this.interval : this.interval : size(this.vergerr_hist), vergerr(windowSize + 1 : end), ...
                 'color', 'b', 'LineWidth', 1.3);

            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Vergence Error [deg]', 'FontSize', 12);
            title('Moving Average of the Vergence Error');
            legend('|verg_{err}|', 'SMA(|verg_{err}|)');
            plotpath = sprintf('%s/mvngAvgVergErr', this.savePath);
            saveas(gcf, plotpath, 'png');

            %% Root Mean Squared Error
            vergerr = this.vergerr_hist(ind);
            rmse = zeros(length(1 : windowSize : length(vergerr) - mod(length(vergerr), 2)), 1); %cut if odd length
            k = 1 : windowSize : length(vergerr) - mod(length(vergerr), 2);
            for i = 1 : length(rmse)
                rmse(i) = mean(vergerr(k(i) : k(i) + windowSize - 1) .^ 2);
            end
            rmse = sqrt(rmse);

            figure;
            hold on;
            grid on;
            plot(windowSize * this.interval : windowSize * this.interval : (length(vergerr) - mod(length(vergerr), 2)) * this.interval, ...
                 rmse, 'LineWidth', 1.3);
            axis([-inf, inf, 0, inf]);
            xlabel(sprintf('Iteration # (windowSize=%d)', windowSize * this.interval), 'FontSize', 12);
            ylabel('RMSE Vergence Error [deg]', 'FontSize', 12);
            title('RMSE of the Vergence Error');
            plotpath = sprintf('%s/rmseVergErr', this.savePath);
            saveas(gcf, plotpath, 'png');

            %% Reconstruction Error
            recerr_L = filter(ones(1, windowSize) / windowSize, 1, this.recerr_hist(ind, 1));
            recerr_S = filter(ones(1, windowSize) / windowSize, 1, this.recerr_hist(ind, 2));

            figure;
            hold on;
            grid on;
            plot((windowSize + 1) * this.interval : this.interval : size(this.recerr_hist), recerr_L(windowSize + 1 : end), 'r', 'LineWidth', 1.3);
            plot((windowSize + 1) * this.interval : this.interval : size(this.recerr_hist), recerr_S(windowSize + 1 : end), 'b', 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Reconstruction Error [AU]', 'FontSize', 12);
            legend('Coarse', 'Fine');
            plotpath = sprintf('%s/recErr', this.savePath);
            saveas(gcf, plotpath, 'png');

            %% Vergence angle
            figure;
            hold on;
            grid on;
            plot(this.verge_desired, 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.3);
            plot(this.verge_actual, 'b', 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Angle [deg]', 'FontSize', 12);
            legend('desired', 'actual');
            title('Vergence');
            plotpath = sprintf('%s/vergenceAngle', this.savePath);
            saveas(gcf, plotpath, 'png');

            %% Muscel graphs
            figure;
            hold on;
            grid on;
            subplot(3, 1, 1);
            plot(this.cmd_hist(:, 2), 'color', [rand, rand, rand], 'LineWidth', 1.3);
            ylabel('Value', 'FontSize', 12);
            title('Absolute Muscle Commands');

            subplot(3, 1, 2);
            plot(this.relCmd_hist, 'color', [rand, rand, rand], 'LineWidth', 1.3);
            ylabel('Value', 'FontSize', 12);
            title('Relative Muscle Commands');

            subplot(3, 1, 3);
            plot(this.metCost_hist, 'color', [rand, rand, rand], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            title('Metabolic Costs');

            plotpath = sprintf('%s/muscleGraphs', this.savePath);
            saveas(gcf, plotpath, 'png');

            %% Weights
            % L1 norm
            figure;
            hold on;
            grid on;
            plot(this.l12_weights(:, 1), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
            plot(this.l12_weights(:, 3), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
            plot(this.l12_weights(:, 5), 'color', [1, 0.5098, 0.1961], 'LineWidth', 1.3);
            % plot(this.l12_weights(:, 7), 'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('\Sigma \midweights\mid', 'FontSize', 12);
            % legend('w_{Vji}', 'w_{Pji}', 'w_{Pkj}', 'w_{Pnji}', 'Location', 'best');
            legend('w_{Vji}', 'w_{Pji}', 'w_{Pkj}', 'Location', 'best');
            title('Model weights (L1)')
            plotpath = sprintf('%s/weightsL1', this.savePath);
            saveas(gcf, plotpath, 'png');

            % L2 norm
            % figure;
            % hold on;
            % grid on;
            % plot(this.l12_weights(:, 2), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
            % plot(this.l12_weights(:, 4), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
            % % plot(this.l12_weights(:, 6), 'color', [1, 0.5098, 0.1961], 'LineWidth', 1.3);
            % % plot(this.l12_weights(:, 8), 'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            % xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            % ylabel('\Sigma weights^{2}', 'FontSize', 12);
            % % legend('w_{Vji}', 'w_{Pji}', 'w_{Pkj}', 'w_{Pnji}', 'Location', 'best');
            % legend('w_{Vji}', 'w_{Pki}', 'Location', 'best');
            % title('Model weights (L2)')
            % plotpath = sprintf('%s/weightsL2', this.savePath);
            % saveas(gcf, plotpath, 'png');

            %% Reward
            figure;
            hold on;
            grid on;
            % r = [- this.lambdaMet * this.metCost_hist, ...
            %      - this.lambdaP2 * this.l12_weights(:, 5), ...
            %      - this.lambdaP1 * this.l12_weights(:, 3), ...
            %      - this.lambdaV * this.l12_weights(:, 1), ...
            %      - this.lambdaRec * (this.recerr_hist(:, 1) + this.recerr_hist(:, 2))];
            r = [- this.lambdaMet * this.metCost_hist, ...
                 - this.lambdaRec * (this.recerr_hist(:, 1) + this.recerr_hist(:, 2))];
            handle = area(r, 'LineStyle','none');
            handle(1).FaceColor = [1, 0.25, 0];
            handle(2).FaceColor = [1, 0.549, 0];
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            % l = legend('\lambdametCost', '\lambdaL1(w_{Pkj})', '\lambdaL1(w_{Pji})', '\lambdaL1(w_{Vji})', '\lambdaRecErr');
            l = legend('\lambdametCost', '\lambdaRecErr');
            if(version('-release') == '2015b')
                l.Location = 'southwest';
            end
            % title('Reward composition (L1)');
            title('Reward composition');
            plotpath = sprintf('%s/rewardComp', this.savePath);
            saveas(gcf, plotpath, 'png');

            %% delta_MF(Vergence_error)
            %
            % desired_angle_min @ 2m = 2 * atand(baseline / (2 * 2)) = 1.6042
            % desired_angle_max @ 0.5m = 2 * atand(baseline / (2 * 0.5)) = 6.4104
            % actual_distance_min = (baseline / 2) / tand(results_deg(11,1)*2 / 2) = 0.0389 [m]
            % actual_distance_max = (baseline / 2) / tand(results_deg(1,1)*2 / 2) = 3.2219 [m]
            % angle_min @ actual_distance_max = results_deg(1,1) * 2 = 0.9958 deg
            % angle_max @ actual_distance_min = results_deg(11,1) * 2 = 71.5164 deg
            % verg_err_min = desired_angle_min - angle_max = 1.6042 - 71.5164 = -69.9122
            % Verg_err_max = desired_angle_max - angle_min = 6.4104 - 0.9958 = 5.4146

            degrees = load('Degrees.mat');
            % angleMin = degrees.results_deg(1, 1);
            % angleMax = degrees.results_deg(11, 1);
            % vergErrMin = this.desiredAngleMin - angleMax;
            % vergErrMax = this.desiredAngleMax - angleMin;

            resolution = 10001;
            approx = spline(1:11, degrees.results_deg(:, 1));

            xValPos = ppval(approx, 1:0.001:11)';
            yValPos = linspace(0, 1, resolution)';

            xValNeg = flipud(ppval(approx, 1:0.001:11)' * -1);
            yValNeg = linspace(-1, 0, resolution)';

            % calculate muscle function :=  mf(vergence_angle) = muscle force [single muuscle]
            mf = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
            dmf = diff(mf(1:2, 1)); % delta in angle
            indZero = find(mf(:, 2) == 0); % MF == 0_index
            indMaxFix = find(mf(:, 1) <= this.desiredAngleMin + dmf & mf(:, 1) >= this.desiredAngleMin - dmf); % MF(desiredAngleMin)_index
            indMinFix = find(mf(:, 1) <= this.desiredAngleMax + dmf & mf(:, 1) >= this.desiredAngleMax - dmf); % MF(desiredAngleMax)_index

            % perfect_response := [max_fixation_x, max_fixation_y, min_fixation_x, min_fixation_y]
            % x = vergenceError, y = deltaMuscelForce
            perfectResponseMaxFix = [(mf(indMaxFix, 1) - flipud(mf(indMaxFix : end, 1))) * 2, ...
                                     (mf(indMaxFix, 2) - flipud(mf(indMaxFix : end, 2))); ...
                                     (mf(indMaxFix, 1) - flipud(mf(indZero : indMaxFix - 1, 1))) * 2, ...
                                     (mf(indMaxFix, 2) - flipud(mf(indZero : indMaxFix - 1, 2)))];

            perfectResponseMinFix = [(mf(indMinFix, 1) - flipud(mf(indMinFix : end, 1))) * 2, ...
                                     (mf(indMinFix, 2) - flipud(mf(indMinFix : end, 2))); ...
                                     (mf(indMinFix, 1) - flipud(mf(indZero : indMinFix - 1, 1))) * 2, ...
                                     (mf(indMinFix, 2) - flipud(mf(indZero : indMinFix - 1, 2)))];

            perfectResponse = [perfectResponseMaxFix, perfectResponseMinFix];
            actualResponse = [this.vergerr_hist, this.relCmd_hist];

            % observation Window, i.e. plot statistics over last #obsWin iterations
            obsWin = 1000;
            if (size(actualResponse, 1) <= obsWin)
                obsWin = size(actualResponse, 1) - 1;
            end
            nVal = 20; % #bins of statistics

            tmpRsp = sortrows(actualResponse(end - obsWin : end, :));
            deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nVal;
            % tmp = nVal x 3 = [index_x = vergence_error angle, mean_muscle_force, std_muscle_force]
            tmp = zeros(nVal, 3);

            for i = 1:nVal
                tmp(i, 1) = (tmpRsp(1, 1) - deltaVergErr / 2) + i * deltaVergErr;
                tmp(i, 2) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                             & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
                tmp(i, 3) = std(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                            & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
            end
            actualResponseStat = tmp;
            actualResponseStat(isnan(actualResponseStat(:, 2)), :) = []; % drop NaN elements

            figure;
            hold on;
            grid on;
            % perfect response to vergence error
            plot(perfectResponse(:, 1), perfectResponse(:, 2), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
            plot(perfectResponse(:, 3), perfectResponse(:, 4), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
            % error bars of actual response of the model
            errorbar(actualResponseStat(:, 1), actualResponseStat(:, 2), actualResponseStat(:, 3),'color', [1, 0.5098, 0.1961], 'LineWidth', 0.9);
            % actual response of the model
            plot(actualResponseStat(:, 1), actualResponseStat(:, 2),'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            l = legend('perfect (fixDist_{max})', 'perfect (fixDist_{min})', 'actual');
            if(version('-release') == '2015b')
                l.FontSize = 7;
                l.Orientation = 'horizontal';
                l.Location = 'southoutside';
            end
            % adjust axis to actual response ranges + std deviation
            xmin = min(actualResponseStat(:, 1)) * 1.1;
            xmax = max(actualResponseStat(:, 1)) * 1.1;
            if (xmin >= xmax)
                xmin = -5.5;
                xmax = 5.5;
            end
            ymin = min([-0.1, (min(actualResponseStat(:, 2)) - max(actualResponseStat(:, 3))) * 1.1]);
            ymax = max([max(perfectResponse(:, 4)) * 1.1, (max(actualResponseStat(:, 2)) + max(actualResponseStat(:, 3))) * 1.1]);
            if (ymin >= ymax)
                ymin = -0.1;
                ymax = 0.1;
            end
            plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
            plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
            axis([xmin, xmax, ymin, ymax]);
            xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
            ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
            title(strcat('\Delta MF(verg_{err}) response after ', sprintf(' %d iterations', size(actualResponse, 1) - obsWin)));
            plotpath = sprintf('%s/deltaMFasFktVerErr', this.savePath);
            saveas(gcf, plotpath, 'png');
        end

        %% plot & save delta_MF(Vergence_error)
        % observation Window = obsWin, i.e. plot statistics over last #obsWin iterations
        function deltaMFplotObsWin(this, obsWin)
            % desired_angle_min @ 2m = 2 * atand(baseline / (2 * 2)) = 1.6042
            % desired_angle_max @ 0.5m = 2 * atand(baseline / (2 * 0.5)) = 6.4104
            % actual_distance_min = (baseline / 2) / tand(results_deg(11,1)*2 / 2) = 0.0389 [m]
            % actual_distance_max = (baseline / 2) / tand(results_deg(1,1)*2 / 2) = 3.2219 [m]
            % angle_min @ actual_distance_max = results_deg(1,1) * 2 = 0.9958 deg
            % angle_max @ actual_distance_min = results_deg(11,1) * 2 = 71.5164 deg
            % verg_err_min = desired_angle_min - angle_max = 1.6042 - 71.5164 = -69.9122
            % Verg_err_max = desired_angle_max - angle_min = 6.4104 - 0.9958 = 5.4146

            degrees = load('Degrees.mat');
            % angleMin = degrees.results_deg(1, 1);
            % angleMax = degrees.results_deg(11, 1);
            % vergErrMin = this.desiredAngleMin - angleMax;
            % vergErrMax = this.desiredAngleMax - angleMin;

            resolution = 10001;
            approx = spline(1:11, degrees.results_deg(:, 1));

            xValPos = ppval(approx, 1:0.001:11)';
            yValPos = linspace(0, 1, resolution)';

            xValNeg = flipud(ppval(approx, 1:0.001:11)' * -1);
            yValNeg = linspace(-1, 0, resolution)';

            % calculate muscle function :=  mf(vergence_angle) = muscle force [single muuscle]
            mf = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
            dmf = diff(mf(1:2, 1)); % delta in angle
            indZero = find(mf(:, 2) == 0); % MF == 0_index
            indMaxFix = find(mf(:, 1) <= this.desiredAngleMin + dmf & mf(:, 1) >= this.desiredAngleMin - dmf); % MF(desiredAngleMin)_index
            indMinFix = find(mf(:, 1) <= this.desiredAngleMax + dmf & mf(:, 1) >= this.desiredAngleMax - dmf); % MF(desiredAngleMax)_index

            % perfect_response := [max_fixation_x, max_fixation_y, min_fixation_x, min_fixation_y]
            % x = vergenceError, y = deltaMuscelForce
            perfectResponseMaxFix = [(mf(indMaxFix, 1) - flipud(mf(indMaxFix : end, 1))) * 2, ...
                                     (mf(indMaxFix, 2) - flipud(mf(indMaxFix : end, 2))); ...
                                     (mf(indMaxFix, 1) - flipud(mf(indZero : indMaxFix - 1, 1))) * 2, ...
                                     (mf(indMaxFix, 2) - flipud(mf(indZero : indMaxFix - 1, 2)))];

            perfectResponseMinFix = [(mf(indMinFix, 1) - flipud(mf(indMinFix : end, 1))) * 2, ...
                                     (mf(indMinFix, 2) - flipud(mf(indMinFix : end, 2))); ...
                                     (mf(indMinFix, 1) - flipud(mf(indZero : indMinFix - 1, 1))) * 2, ...
                                     (mf(indMinFix, 2) - flipud(mf(indZero : indMinFix - 1, 2)))];

            perfectResponse = [perfectResponseMaxFix, perfectResponseMinFix];
            actualResponse = [this.vergerr_hist, this.relCmd_hist];

            if (size(actualResponse, 1) <= obsWin)
                obsWin = size(actualResponse, 1) - 1;
            end
            nVal = 20; % #bins of statistics

            tmpRsp = sortrows(actualResponse(end - obsWin : end, :));
            deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nVal;
            % tmp = nVal x 3 = [index_x = vergence_error angle, mean_muscle_force, std_muscle_force]
            tmp = zeros(nVal, 3);

            for i = 1:nVal
                tmp(i, 1) = (tmpRsp(1, 1) - deltaVergErr / 2) + i * deltaVergErr;
                tmp(i, 2) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                             & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
                tmp(i, 3) = std(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                            & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
            end
            actualResponseStat = tmp;
            actualResponseStat(isnan(actualResponseStat(:, 2)), :) = []; % drop NaN elements

            figure;
            hold on;
            grid on;
            % perfect response to vergence error
            plot(perfectResponse(:, 1), perfectResponse(:, 2), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
            plot(perfectResponse(:, 3), perfectResponse(:, 4), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
            % error bars of actual response of the model
            errorbar(actualResponseStat(:, 1), actualResponseStat(:, 2), actualResponseStat(:, 3),'color', [1, 0.5098, 0.1961], 'LineWidth', 0.9);
            % actual response of the model
            plot(actualResponseStat(:, 1), actualResponseStat(:, 2),'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            l = legend('perfect (fixDist_{max})', 'perfect (fixDist_{min})', 'actual');
            if(version('-release') == '2015b')
                l.FontSize = 7;
                l.Orientation = 'horizontal';
                l.Location = 'southoutside';
            end
            % adjust axis to actual response ranges + std deviation
            xmin = min(actualResponseStat(:, 1)) * 1.1;
            xmax = max(actualResponseStat(:, 1)) * 1.1;
            if (xmin >= xmax)
                xmin = -5.5;
                xmax = 5.5;
            end
            ymin = min([-0.1, (min(actualResponseStat(:, 2)) - max(actualResponseStat(:, 3))) * 1.1]);
            ymax = max([max(perfectResponse(:, 4)) * 1.1, (max(actualResponseStat(:, 2)) + max(actualResponseStat(:, 3))) * 1.1]);
            if (ymin >= ymax)
                ymin = -0.1;
                ymax = 0.1;
            end
            plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
            plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
            axis([xmin, xmax, ymin, ymax]);
            xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
            ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
            title(strcat('\Delta MF(verg_{err}) response after ', sprintf(' %d iterations', size(actualResponse, 1) - obsWin)));
            plotpath = sprintf('%s/deltaMFasFktVerErrObsWin', this.savePath);
            saveas(gcf, plotpath, 'png');
        end

        %% plot & save delta_MF(Vergence_error)
        % startIter, endIter = plot statistics over #startIter : endIter iterations
        function deltaMFplotStartEnd(this, startIter, endIter)
            % desired_angle_min @ 2m = 2 * atand(baseline / (2 * 2)) = 1.6042
            % desired_angle_max @ 0.5m = 2 * atand(baseline / (2 * 0.5)) = 6.4104
            % actual_distance_min = (baseline / 2) / tand(results_deg(11,1)*2 / 2) = 0.0389 [m]
            % actual_distance_max = (baseline / 2) / tand(results_deg(1,1)*2 / 2) = 3.2219 [m]
            % angle_min @ actual_distance_max = results_deg(1,1) * 2 = 0.9958 deg
            % angle_max @ actual_distance_min = results_deg(11,1) * 2 = 71.5164 deg
            % verg_err_min = desired_angle_min - angle_max = 1.6042 - 71.5164 = -69.9122
            % Verg_err_max = desired_angle_max - angle_min = 6.4104 - 0.9958 = 5.4146

            degrees = load('Degrees.mat');
            % angleMin = degrees.results_deg(1, 1);
            % angleMax = degrees.results_deg(11, 1);
            % vergErrMin = this.desiredAngleMin - angleMax;
            % vergErrMax = this.desiredAngleMax - angleMin;

            resolution = 10001;
            approx = spline(1:11, degrees.results_deg(:, 1));

            xValPos = ppval(approx, 1:0.001:11)';
            yValPos = linspace(0, 1, resolution)';

            xValNeg = flipud(ppval(approx, 1:0.001:11)' * -1);
            yValNeg = linspace(-1, 0, resolution)';

            % calculate muscle function :=  mf(vergence_angle) = muscle force [single muuscle]
            mf = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
            dmf = diff(mf(1:2, 1)); % delta in angle
            indZero = find(mf(:, 2) == 0); % MF == 0_index
            indMaxFix = find(mf(:, 1) <= this.desiredAngleMin + dmf & mf(:, 1) >= this.desiredAngleMin - dmf); % MF(desiredAngleMin)_index
            indMinFix = find(mf(:, 1) <= this.desiredAngleMax + dmf & mf(:, 1) >= this.desiredAngleMax - dmf); % MF(desiredAngleMax)_index

            % perfect_response := [max_fixation_x, max_fixation_y, min_fixation_x, min_fixation_y]
            % x = vergenceError, y = deltaMuscelForce
            perfectResponseMaxFix = [(mf(indMaxFix, 1) - flipud(mf(indMaxFix : end, 1))) * 2, ...
                                     (mf(indMaxFix, 2) - flipud(mf(indMaxFix : end, 2))); ...
                                     (mf(indMaxFix, 1) - flipud(mf(indZero : indMaxFix - 1, 1))) * 2, ...
                                     (mf(indMaxFix, 2) - flipud(mf(indZero : indMaxFix - 1, 2)))];

            perfectResponseMinFix = [(mf(indMinFix, 1) - flipud(mf(indMinFix : end, 1))) * 2, ...
                                     (mf(indMinFix, 2) - flipud(mf(indMinFix : end, 2))); ...
                                     (mf(indMinFix, 1) - flipud(mf(indZero : indMinFix - 1, 1))) * 2, ...
                                     (mf(indMinFix, 2) - flipud(mf(indZero : indMinFix - 1, 2)))];

            perfectResponse = [perfectResponseMaxFix, perfectResponseMinFix];
            actualResponse = [this.vergerr_hist, this.relCmd_hist];

            if (endIter > size(actualResponse, 1))
                endIter = size(actualResponse, 1);
            end
            nVal = 20; % #bins of statistics

            tmpRsp = sortrows(actualResponse(startIter : endIter, :));
            deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nVal;
            % tmp = nVal x 3 = [index_x = vergence_error angle, mean_muscle_force, std_muscle_force]
            tmp = zeros(nVal, 3);

            for i = 1:nVal
                tmp(i, 1) = (tmpRsp(1, 1) - deltaVergErr / 2) + i * deltaVergErr;
                tmp(i, 2) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                             & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
                tmp(i, 3) = std(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                            & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
            end
            actualResponseStat = tmp;
            actualResponseStat(isnan(actualResponseStat(:, 2)), :) = []; % drop NaN elements

            figure;
            hold on;
            grid on;
            % perfect response to vergence error
            plot(perfectResponse(:, 1), perfectResponse(:, 2), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
            plot(perfectResponse(:, 3), perfectResponse(:, 4), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
            % error bars of actual response of the model
            errorbar(actualResponseStat(:, 1), actualResponseStat(:, 2), actualResponseStat(:, 3),'color', [1, 0.5098, 0.1961], 'LineWidth', 0.9);
            % actual response of the model
            plot(actualResponseStat(:, 1), actualResponseStat(:, 2),'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            l = legend('perfect (fixDist_{max})', 'perfect (fixDist_{min})', 'actual');
            if(version('-release') == '2015b')
                l.FontSize = 7;
                l.Orientation = 'horizontal';
                l.Location = 'southoutside';
            end
            % adjust axis to actual response ranges + std deviation
            xmin = min(actualResponseStat(:, 1)) * 1.1;
            xmax = max(actualResponseStat(:, 1)) * 1.1;
            if (xmin >= xmax)
                xmin = -5.5;
                xmax = 5.5;
            end
            ymin = min([-0.1, (min(actualResponseStat(:, 2)) - max(actualResponseStat(:, 3))) * 1.1]);
            ymax = max([max(perfectResponse(:, 4)) * 1.1, (max(actualResponseStat(:, 2)) + max(actualResponseStat(:, 3))) * 1.1]);
            if (ymin >= ymax)
                ymin = -0.1;
                ymax = 0.1;
            end
            plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
            plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
            axis([xmin, xmax, ymin, ymax]);
            xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
            ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
            title(strcat('\Delta MF(verg_{err}) at iteration ', sprintf(' %d through %d', startIter, endIter)));
            plotpath = sprintf('%s/deltaMFasFktVerErrStartEnd', this.savePath);
            saveas(gcf, plotpath, 'png');
        end
    end
end

classdef Model < handle
    properties
        scmodel_Large;      %SparseCoding class (downsampled version)
        scmodel_Small;      %SparseCoding class
        rlmodel;            %ReinforcementLearning class

        f;                  %focal length [px]
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
        simulated_time;     %how long did the model take to be learned (min)

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
        savePath;           %where all the data are stored
    end

    methods
        function obj = Model(PARAM)
            obj.learnedFile = PARAM{1}{1};
            obj.textureFile = PARAM{1}{2};
            obj.trainTime = PARAM{1}{3};
            obj.sparseCodingType = PARAM{1}{4};
            obj.f = PARAM{1}{5};
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
            obj.feature_hist = zeros(obj.trainTime, obj.rlmodel.S0);
            obj.cmd_hist = zeros(obj.trainTime, 2);
            obj.relCmd_hist = zeros(obj.trainTime, 1);
            obj.l12_weights = zeros(obj.trainTime, 4);
            obj.reward_hist = zeros(obj.trainTime, 1);
            obj.metCost_hist = zeros(obj.trainTime, 1);
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

        %% Plotting errors
        function errPlot(this)
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
            legend('verg err', 'SMA(verg err)');

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
            xlabel(sprintf('Iteration # (windowSize=%d)', windowSize * this.interval), 'FontSize', 12);
            ylabel('RMSE Vergence Error [deg]', 'FontSize', 12);
            title('RMSE of the Vergence Error');

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
            legend('verg err', 'SMA(verg err)');
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
            % xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            title('Absolute Muscle Commands');

            subplot(3, 1, 2);
            plot(this.relCmd_hist, 'color', [rand, rand, rand], 'LineWidth', 1.3);
            % xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
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
            figure;
            hold on;
            grid on;
            plot(this.l12_weights(:, 1), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
            plot(this.l12_weights(:, 3), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
            % plot(this.l12_weights(:, 5), 'color', [1, 0.5098, 0.1961], 'LineWidth', 1.3);
            % plot(this.l12_weights(:, 7), 'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('\Sigma \midweights\mid', 'FontSize', 12);
            % legend('w_{Vji}', 'w_{Pji}', 'w_{Pkj}', 'w_{Pnji}', 'Location', 'best');
            legend('w_{Vji}', 'w_{Pki}', 'Location', 'best');
            title('Model weights (L1)')
            plotpath = sprintf('%s/weightsL1', this.savePath);
            saveas(gcf, plotpath, 'png');

            figure;
            hold on;
            grid on;
            plot(this.l12_weights(:, 2), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
            plot(this.l12_weights(:, 4), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
            % plot(this.l12_weights(:, 6), 'color', [1, 0.5098, 0.1961], 'LineWidth', 1.3);
            % plot(this.l12_weights(:, 8), 'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('\Sigma weights^{2}', 'FontSize', 12);
            % legend('w_{Vji}', 'w_{Pji}', 'w_{Pkj}', 'w_{Pnji}', 'Location', 'best');
            legend('w_{Vji}', 'w_{Pki}', 'Location', 'best');
            title('Model weights (L2)')
            plotpath = sprintf('%s/weightsL2', this.savePath);
            saveas(gcf, plotpath, 'png');

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
            area(r, 'LineStyle','none');
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            % l = legend('\lambdametCost', '\lambdaL1(w_{Pkj})', '\lambdaL1(w_{Pji})', '\lambdaL1(w_{Vji})', '\lambdaRecErr');
            l = legend('\lambdametCost', '\lambdaRecErr');
            % l.FontSize = 7;
            l.Location = 'southwest';
            title('Reward composition (L1)');
            plotpath = sprintf('%s/rewardCompL1', this.savePath);
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
            angleMin = degrees.results_deg(1, 1);
            angleMax = degrees.results_deg(11, 1);
            vergErrMin = this.desiredAngleMin - angleMax;
            vergErrMax = this.desiredAngleMax - angleMin;

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
            l.FontSize = 7;
            l.Orientation = 'horizontal';
            l.Location = 'southoutside';
            % axis
            xmin = min(actualResponseStat(:, 1)) * 1.1;
            xmax = max(actualResponseStat(:, 1)) * 1.1;
            ymin = min([-0.1, (min(actualResponseStat(:, 2)) - max(actualResponseStat(:, 3))) * 1.2]);
            ymax = max([max(perfectResponse(:, 4)) * 1.2, (max(actualResponseStat(:, 2)) + max(actualResponseStat(:, 3))) * 1.2]);
            plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
            plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
            axis([xmin, xmax, ymin, ymax]);
            xlabel(sprintf('Vergence Error [deg], bin size = %.3g deg', deltaVergErr), 'FontSize', 12);
            ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
            title(strcat('\Delta MF(Vergence_{error}) response after ', sprintf(' %d iterations', size(actualResponse, 1) - obsWin)));
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
            angleMin = degrees.results_deg(1, 1);
            angleMax = degrees.results_deg(11, 1);
            vergErrMin = this.desiredAngleMin - angleMax;
            vergErrMax = this.desiredAngleMax - angleMin;

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
            l.FontSize = 7;
            l.Orientation = 'horizontal';
            l.Location = 'southoutside';
            % axis
            xmin = min(actualResponseStat(:, 1)) * 1.1;
            xmax = max(actualResponseStat(:, 1)) * 1.1;
            ymin = min([-0.1, (min(actualResponseStat(:, 2)) - max(actualResponseStat(:, 3))) * 1.2]);
            ymax = max([max(perfectResponse(:, 4)) * 1.2, (max(actualResponseStat(:, 2)) + max(actualResponseStat(:, 3))) * 1.2]);
            plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
            plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
            axis([xmin, xmax, ymin, ymax]);
            xlabel(sprintf('Vergence Error [deg], bin size = %.3g deg', deltaVergErr), 'FontSize', 12);
            ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
            title(strcat('\Delta MF(Vergence_{error}) response after ', sprintf(' %d iterations', size(actualResponse, 1) - obsWin)));
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
            angleMin = degrees.results_deg(1, 1);
            angleMax = degrees.results_deg(11, 1);
            vergErrMin = this.desiredAngleMin - angleMax;
            vergErrMax = this.desiredAngleMax - angleMin;

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

            %%% generating fixed distances
            % [vergErrs, relCmds] = generateRelCmds(this, [0.5,1,2], [-5:1:5], 10);
            % vergRange = [-5 : 0.5 : 5];
            % [vergErrs, relCmds] = generateRelCmds(this, [0.5, 2], vergRange, 10)
            % actualResponse = [vergErrs, relCmds];

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
            l.FontSize = 7;
            l.Orientation = 'horizontal';
            l.Location = 'southoutside';
            % axis
            xmin = min(actualResponseStat(:, 1)) * 1.1;
            xmax = max(actualResponseStat(:, 1)) * 1.1;
            ymin = min([-0.1, (min(actualResponseStat(:, 2)) - max(actualResponseStat(:, 3))) * 1.2]);
            ymax = max([max(perfectResponse(:, 4)) * 1.2, (max(actualResponseStat(:, 2)) + max(actualResponseStat(:, 3))) * 1.2]);
            plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
            plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
            axis([xmin, xmax, ymin, ymax]);
            xlabel(sprintf('Vergence Error [deg], bin size = %.3g deg', deltaVergErr), 'FontSize', 12);
            ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
            title(strcat('\Delta MF(Vergence_{error}) betweem iteration ', sprintf(' %d and %d', startIter, endIter)));
            plotpath = sprintf('%s/deltaMFasFktVerErrStartEnd', this.savePath);
            saveas(gcf, plotpath, 'png');
        end

        %% plot & save delta_MF(Vergence_error)
        % objRange = range of object distances being tested
        % vergRange = range of vergences being tested
        % repeat = #repetitions of testing procedure
        function deltaMFplotGenDist(this, objRange, vergRange, repeat)
            % desired_angle_min @ 2m = 2 * atand(baseline / (2 * 2)) = 1.6042
            % desired_angle_max @ 0.5m = 2 * atand(baseline / (2 * 0.5)) = 6.4104
            % actual_distance_min = (baseline / 2) / tand(results_deg(11,1)*2 / 2) = 0.0389 [m]
            % actual_distance_max = (baseline / 2) / tand(results_deg(1,1)*2 / 2) = 3.2219 [m]
            % angle_min @ actual_distance_max = results_deg(1,1) * 2 = 0.9958 deg
            % angle_max @ actual_distance_min = results_deg(11,1) * 2 = 71.5164 deg
            % verg_err_min = desired_angle_min - angle_max = 1.6042 - 71.5164 = -69.9122
            % Verg_err_max = desired_angle_max - angle_min = 6.4104 - 0.9958 = 5.4146

            degrees = load('Degrees.mat');
            angleMin = degrees.results_deg(1, 1);
            angleMax = degrees.results_deg(11, 1);
            vergErrMin = this.desiredAngleMin - angleMax;
            vergErrMax = this.desiredAngleMax - angleMin;

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

            %%% generating fixed distances
            % actualResponse = [vergErrs, relCmds];
            responseResults = generateRelCmds(this, objRange, vergRange, repeat);
            actualResponse = [responseResults.vergErrs, responseResults.relCmds];
            nVal = size(vergRange, 2); % #bins of statistics

            tmpRsp = sortrows(actualResponse);
            deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nVal;
            % tmp = [index_x = vergence_error angle, mean_muscle_force, std_muscle_force]
            tmp = zeros(nVal, 3);

            for i = 1:nVal
                tmp(i, 1) = vergRange(i);
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
            l.FontSize = 7;
            l.Orientation = 'horizontal';
            l.Location = 'southoutside';
            % axis
            xmin = min(actualResponseStat(:, 1)) * 1.1;
            xmax = max(actualResponseStat(:, 1)) * 1.1;
            ymin = min([-0.1, (min(actualResponseStat(:, 2)) - max(actualResponseStat(:, 3))) * 1.2]);
            ymax = max([max(perfectResponse(:, 4)) * 1.2, (max(actualResponseStat(:, 2)) + max(actualResponseStat(:, 3))) * 1.2]);
            plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
            plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
            axis([xmin, xmax, ymin, ymax]);
            xlabel(sprintf('Vergence Error [deg], bin size = %.3g deg', deltaVergErr), 'FontSize', 12);
            ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
            title('\Delta MF(Vergence_{error}) response at Testing procedure');
            if ~isempty(this.savePath)
                plotpath = sprintf('%s/deltaMFasFktVerErrGenDist', this.savePath);
                saveas(gcf, plotpath, 'png');
            end

        end

        %% plot & save reconstructionErr(Vergence_error)
        % objRange = range of object distances being tested
        % vergRange = range of vergences being tested
        % repeat = #repetitions of testing procedure
        function recErrPlotGenDist(this, objRange, vergRange, repeat)
            %plotting the resonstruction error of basis functions over
            %different disparities
            responseResults = generateRelCmds(this, objRange, vergRange, repeat);
            recResponse = [responseResults.vergErrs, responseResults.recErrs];
            recResponseLarge = [responseResults.vergErrs, responseResults.recErrsLarge];
            recResponseSmall = [responseResults.vergErrs, responseResults.recErrsSmall];
            nVal = size(vergRange, 2); % #bins of statistics

            %calculate mean and std of reconstruction error
            tmpRsp = sortrows(recResponse);
            deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nVal;
            % tmp = [index_x = vergence_error angle, mean_recError, std_recError]
            tmp = zeros(nVal, 3);

            for i = 1:nVal
                tmp(i, 1) = vergRange(i);
                tmp(i, 2) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                             & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
                tmp(i, 3) = std(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                            & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
            end
            recErrs = tmp;
            recErrs(isnan(recErrs(:, 2)), :) = []; % drop NaN elements

            %calculate mean and std of large scale reconstruction error
            tmpRsp = sortrows(recResponseLarge);
            deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nVal;
            % tmp = [index_x = vergence_error angle, mean_recError, std_recError]
            tmp = zeros(nVal, 3);

            for i = 1:nVal
                tmp(i, 1) = vergRange(i);
                tmp(i, 2) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                             & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
                tmp(i, 3) = std(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                            & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
            end
            recErrsLarge = tmp;
            recErrsLarge(isnan(recErrsLarge(:, 2)), :) = []; % drop NaN elements

            %calculate mean and std of small scale reconstruction error
            tmpRsp = sortrows(recResponseSmall);
            deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nVal;
            % tmp = [index_x = vergence_error angle, mean_recError, std_recError]
            tmp = zeros(nVal, 3);

            for i = 1:nVal
                tmp(i, 1) = vergRange(i);
                tmp(i, 2) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                             & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
                tmp(i, 3) = std(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                            & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
            end
            recErrsSmall = tmp;
            recErrsSmall(isnan(recErrsSmall(:, 2)), :) = []; % drop NaN elements

            figure;
            hold on;
            grid on;
            % error bars of reconstruction of the model
            errorbar(recErrs(:, 1), recErrs(:, 2), recErrs(:, 3), 'LineWidth', 0.9); %'color', [1, 0.5098, 0.1961],
            errorbar(recErrsLarge(:, 1), recErrsLarge(:, 2), recErrsLarge(:, 3), 'LineWidth', 0.9);%'color', [1, 0.5098, 0.1961],
            errorbar(recErrsSmall(:, 1), recErrsSmall(:, 2), recErrsSmall(:, 3), 'LineWidth', 0.9);%'color', [1, 0.5098, 0.1961],
            % reconstruction of the model
%             plot(recErrs(:, 1), recErrs(:, 2),'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            l = legend('reconstruction Error', 'small scale recErr', 'large scale recErr');
            l.FontSize = 7;
            l.Orientation = 'horizontal';
            l.Location = 'southoutside';
            % axis
%             xmin = -10;
%             xmax = -xmin;
%             ymin = -0.1;
%             ymax = -ymin;
%             plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
%             plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
%             axis([xmin, xmax, ymin, ymax]);
            xlabel(sprintf('Vergence Error [deg], bin size = %.3g deg', deltaVergErr), 'FontSize', 12);
            ylabel('resonstruction Error', 'FontSize', 12);
            title(sprintf('Reconstruction Error over different disparities\nobject distances: [%s]',num2str(objRange)));

            if ~ isempty(this.savePath)
                plotpath = sprintf('%s/recErrVsVerErrGenDist', this.savePath);
                saveas(gcf, plotpath, 'png');
            end
        end

        % DEPRECATED
        %% Train the whole system
        %
        % trainSC indicates whether the sparse coding part will be
        % trained
        % trainRL indicates whether the reinforcement learning part will
        % be trained
        % debugmode indicates whether some intermedia should be recorded;
        function traindebug(this, trainSC, trainRL)
            command = [];
            Current_View = this.envmodel.generate(1, command);
            [~, reward] = this.generateFR(Current_View);
            this.rlmodel.J = reward;
            Errortmp = [];  %temp variable to compute rec error (Luca)

            for t = 1:this.T_train
                Current_View = this.envmodel.generate(~mod(t - 1, this.interval), command);

                if (trainSC && ~mod(t - 1, this.SCInterval) && (t<this.stopSC))
                    %train 2 sparse coding models
                    this.scmodel_Large.stepTrain(Current_View{1});
                    this.scmodel_Small.stepTrain(Current_View{2});
                end

                [feature, reward, Error] = this.generateFR(Current_View); %Added Error output (Luca)

                %save Rec Error during training (Luca)
                Errortmp = [Errortmp Error];

                if (trainRL)
                    command = this.rlmodel.stepTrain(feature, reward, mod(t - 1, this.interval));
                else
                    command = this.rlmodel.Action(randi(this.rlmodel.Action_num));
                end

                if (~mod(t, ceil(this.T_train/100)))
                    disp([num2str(t * 100/this.T_train) '% is finished']);

                    %**added by Luca - Display feature and weights value***
                    %every 1/100 of the total # of iterations
                    if (trainRL)
                        %save norm of weights
                        win = norm(this.rlmodel.Weights{2});                %norm of input weights
                        wout = sqrt(sum(this.rlmodel.Weights{1}.^2, 2));    %norm of output weights

                        w = [win; mean(wout)];
                        this.weights_hist = [this.weights_hist w];

                        sprintf('Norm of Input weights = %.3g \nMean of Norm of Output weights = %.3g \nMin of Output weights = %.3g \nMax of Output weights = %.3g \nInput mean = %3.g \nInput variance = %.3g ', ...
                                mean(win), mean(wout), min(wout), max(wout), mean(this.rlmodel.X), var(this.rlmodel.X))

                        %save TD error and norm of natural gradient
                        this.td_hist = [this.td_hist this.rlmodel.td];
                        this.g_hist = [this.g_hist norm(this.rlmodel.g)];

                        %display mean and variance of input vector

                        %display current view
                        % figure(1);
                        % clf
                        % subplot(1, 2, 1);
                        % display_network(Current_View(1:100, :), 0);
                        % subplot(1, 2, 2);
                        % display_network(Current_View(101:end, :), 0);

                        %display mean squared disparity
                        ds = this.envmodel.Disphist;
                        ds = sqrt(mean (ds (end-this.T_train/100:end).^2 ) );
                        disp(['rms disparity = ' num2str(ds)])
                    end
                end
            end
        end

        % DEPRICATED
        %% Test the greedy policy of RL model to
        function testRL(this)
            for (t = 1:this.trainTime)
                Current_View = this.envmodel.generate_test(t);      %generate image
                feature = this.generateFR(Current_View);            %feature
                command = this.rlmodel.softmaxAct(feature);         %generate command (and saves action prob in pol_hist)

                this.feature_hist = [this.feature_hist feature];    %save feature vector activation

                if (~mod(t, ceil(this.trainTime/100)))
                    disp([num2str(t*100/this.trainTime) '% is finished']);
                end
            end
        end

        % DEPRICATED
        %% Display following variables
        % - Histogram of disparities
        % - Histogram of actions
        % - Histogram of vergence values
        % - Moving average of Vergence Error
        % - Basis functions (To be included)
        %
        % (if variable test is set to 1 display also):
        % - Desired and fixated depth
        % - Desired and actual depth
        % - Desired and actual vergence
        function testmodel(this, test)
            %some vars (need to be included in the model properties in
            %future versions)
            baseline = 0.0698;  %interocular distance (baseline)
            f = 257.34;         %focal length [px]
            vergeMax = 20;
            windowSize = this.trainTime / 100;

            %HISTOGRAMS
            figure
            if (~isempty(this.Z))
                subplot(311);
                hold on
                %show disparity histogram through training
                thetaZmin = 2 * atan(baseline / (2 * min(this.Z)));   %verge angle for obj at min distance
                thetaZmax = 2 * atan(baseline / (2 * max(this.Z)));   %verge angle for obj at max distance
                dmax = 2 * f * tan(thetaZmin / 2);
                dmin = 2 * f * tan((thetaZmax - vergeMax * pi / 180) / 2);
                hist(this.disp_hist, dmin:2:dmax)
                title('Histogram of Disparities')
            end

            subplot(312);
            hold on
            hist(this.rlmodel.Action_hist, [min(this.rlmodel.Action):max(this.rlmodel.Action)]);
            title('Histogram of Actions')
            subplot(313);
            hold on
            hist(this.verge_actual, 0:0.5:vergeMax);
            title('Histogram of vergence values')

            %MOVING AVERAGE OF VERGENCE ERROR
            ind = 1:this.trainTime;
            i = mod(ind, this.interval);
            ind = find(i); %to exclude trials where vergence is resetted

            %Absolute Mean Error (AME)
            vergerr = filter(ones(1, windowSize) / windowSize, 1, abs(this.vergerr_hist(ind))); %AME
            vergerr = vergerr(windowSize:end); %discard first values (mean not reliable)

            figure;
            hold on;
            plot(vergerr, 'LineWidth', 2);
            xlabel('Trial #', 'FontSize', 14);
            ylabel('Vergence Error [deg]', 'FontSize', 14);
            title('Moving Average of the vergence error')
            set(gca, 'XTick', 0:windowSize * 10:length(vergerr) + windowSize, 'FontSize', 14);
            grid on

            % AME over last trial before object switch
            % windowSize = 200;
            % indB = this.interval - 1:this.interval:ind(end);
            % vergerr = filter(ones(1, windowSize) / windowSize, 1, abs(this.vergerr_hist(indB))); %AME
            % vergerr = vergerr(windowSize:end); %discard first values (mean not reliable)
            % figure;
            % hold on;
            % plot(vergerr, 'LineWidth', 2);
            % xlabel('Iteration #', 'FontSize', 14);
            % ylabel('Vergence Error [deg]', 'FontSize', 14);
            % title('Moving Average of the vergence error (before stimulus change)')
            % set(gca, 'XTick', 0:windowSize * 10:length(vergerr) + windowSize, 'FontSize', 14);
            % grid on

            %WEIGHTS NORM
            t = 0:this.trainTime / 10:this.trainTime;
            wPol = this.rlmodel.Weights_hist{1};
            wVal = this.rlmodel.Weights_hist{2};

            %norm of the weights
            try
                wPol = reshape(wPol, [size(wPol, 1) * size(wPol, 2), 1, size(wPol, 3)]);  %make a single vector of weights
                normwPol = sqrt(sum(wPol.^2, 1));                                       %norm of policy weights
                normwPol = reshape(normwPol, [size(normwPol, 3) 1]);

                normwVal = sqrt(sum(wVal.^2, 2));                                       %norm of value weights
                normwVal = reshape(normwVal, [size(normwVal, 3) 1]);

                figure;
                subplot(1, 2, 1);
                plot(t, normwVal, 'o-');
                title('Value Net: |v(t)|')
                subplot(1, 2, 2);
                plot(t, normwPol, 'o-');
                title('Policy Net: |w(t)|')

                %value network
                for i = 2:size(wVal, 3)
                    deltaD0(i - 1) = (wVal(1, :, i) * wVal(1, :, 1)') / (normwVal(i) * normwVal(1));    %direction change with respect to start
                    deltaN0(i - 1) = normwVal(i) / normwVal(1);                                     %norm change with respect to start

                    deltaD(i - 1) = (wVal(1, :, i - 1) * wVal(1, :, i)') / (normwVal(i - 1) * normwVal(1)); %direction change trial by trial
                    deltaN(i - 1) = normwVal(i) / normwVal(i - 1);                                    %norm change trial by trial
                end

                figure
                subplot(2, 2, 1);
                hold on;
                plot(t(2:end), deltaD0, 'o-');
                title('VALUE NETWORK: w0 x wt')
                subplot(2, 2, 2);
                hold on;
                plot(t(2:end), deltaN0, 'o-');
                title('|wt|/|w0|')
                subplot(2, 2, 3);
                hold on;
                plot(t(2:end), deltaD, 'o-');
                title('w(t+1) x w(t)')
                subplot(2, 2, 4);
                hold on;
                plot(t(2:end), deltaN, 'o-');
                title('|w(t+1)|/|w(t)|')

                %policy network
                %norm of the weights at each iteration

                for i = 2:size(wPol, 3)
                    deltaD0(i - 1) = (wPol(1, :, i) * wPol(1, :, 1)') / (normwPol(i) * normwPol(1));    %direction change with respect to start
                    deltaN0(i - 1) =  normwPol(i) / normwPol(1);                                    %norm change with respect to start

                    deltaD(i - 1) = (wPol(1, :, i - 1) * wPol(1, :, i)') / (normwPol(i - 1) * normwPol(1)); %direction change trial by trial
                    deltaN(i - 1) = normwPol(i) / normwPol(i - 1);                                    %norm change trial by trial
                end

                figure
                subplot(2, 2, 1);
                hold on;
                plot(t(2:end), deltaD0, 'o-');
                title('POLICY NETWORK: w0 x wt')
                subplot(2, 2, 2);
                hold on;
                plot(t(2:end), deltaN0, 'o-');
                title('|wt|/|w0|')
                subplot(2, 2, 3);
                hold on;
                plot(t(2:end), deltaD, 'o-');
                title('w(t+1) x w(t)')
                subplot(2, 2, 4);
                hold on;
                plot(t(2:end), deltaN, 'o-');
                title('|w(t+1)|/|w(t)|')
            catch
                display('error in norm of weigths')
            end

            %RECERROR
            recerr_L = filter(ones(1, windowSize) / windowSize, 1, this.recerr_hist(:, 1));
            recerr_S = filter(ones(1, windowSize) / windowSize, 1, this.recerr_hist(:, 2));
            recerr_L = recerr_L(windowSize:end);
            recerr_S = recerr_S(windowSize:end);
            figure;
            hold on
            plot(recerr_L, 'r', 'LineWidth', 2);
            xlabel('Trial #', 'FontSize', 14);
            ylabel('Reconstruction Error', 'FontSize', 14);
            plot(recerr_S, 'b', 'LineWidth', 2);
            % set(gca, 'XTick', 0:windowSize:t, 'FontSize', 14);
            % grid on;
            legend('Coarse', 'Fine')
            %plot sum of recerrors
            figure;
            hold on
            plot(recerr_L + recerr_S, 'b', 'LineWidth', 2);
            xlabel('Trial #', 'FontSize', 14);
            ylabel('Reconstruction Error', 'FontSize', 14);

            %BASES FUNCTIONS
            %display the basis according to the binocularity
            r = 16;
            c = 18; %how to arrange the basis in rows and cols
            Nscales = 2;                             %Number of Sparse coding models
            len = size(this.scmodel_Small.Basis_hist, 3);  % # of trials saved
            basisTrack = cell(Nscales, len);           %variable to store all the saved basis

            %copy into the new var
            for j = 1:len
                basisTrack{1, j} = this.scmodel_Small.Basis_hist(:, :, j);
            end

            % h = figure;
            % scrsz = get(0, 'ScreenSize');
            % set(h, 'Position', [scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);

            %loop for Large and Small scale sparse coding models (if
            %available)

            for (s = 1:Nscales)
                %sort basis according to LEFT energy norm
                endBasis = basisTrack{s, end}(1:end / 2, :);
                leftEnergy = abs(sum(endBasis.^2) - 0.5);
                [~, I] = sort(leftEnergy);

                % subplot(1, 2, s);
                [di, num] = size(basisTrack{s, 1});

                fun1 = @(blc_struct) padarray(padarray(reshape(permute(padarray(reshape(blc_struct.data, sqrt(di / 2), ...
                        sqrt(di / 2), 2), [1, 1], 'pre'), [1, 3, 2]), (sqrt(di / 2) + 1) * 2, ...
                        sqrt(di / 2) + 1), [1, 1], 'post') - 1, [1 1], 'pre') + 1;

                % UNCOMMENT THIS PART TO SEE MOVIE OF DEVELOPING BASES
                % for j = 1:len
                %     A = basisTrack{s, j}(:, I);
                %     % A = basisTrack{s, j}(:, I(ceil(linspace(1, 288, 200))));
                %     % A = basisTrack{s, j}(:,[13, 25, 50, 150, 160, 19, 93, 195, 132, 116]);
                %     % B = reshape(A, di*sqrt(num/2), sqrt(num/2)*2);
                %     B = reshape(A, di*r, c);
                %     B = B/max(max(abs(B))) + 0.5;
                %     C = padarray(padarray(blockproc(B,[di, 1], fun1)-1,[1 1],'post')+1,[2, 2]);
                %     imshow(C);
                %     title(num2str(this.T_train*0.1*(j-1)));
                %     drawnow;
                %     pause(.5);
                % end

                %Display Histogram of Norm or Energy of Left (non blurred) Bases
                leftNorm = abs(sum(endBasis.^2));
                figure(10);
                subplot(1, 2, s);
                hold on;
                hist(leftNorm, 0:0.1:1);
                xlabel('Basis norm');
                ylabel('# of bases');
            end

            if (test)
                %compute mean and std dev of vergence error before stimulus
                %change
                ind = this.interval - 1:this.interval:this.T_train - 1;
                merr = mean(abs(this.vergerr_hist(ind)));
                serr = std(abs(this.vergerr_hist(ind)));

                %display desired vergence and vergence error
                figure;
                subplot(211);
                hold on
                plot(this.verge_actual, 'r.-', 'LineWidth', 2, 'MarkerSize', 18);
                xlabel('Iteration #', 'FontSize', 14);
                ylabel('Vergence [deg]', 'FontSize', 14);
                plot(this.vergerr_hist + this.verge_actual, 'b.-', 'LineWidth', 2, 'MarkerSize', 18);
                xlabel('Iteration #', 'FontSize', 14);
                grid on
                legend('Actual Vergence', 'Desired Vergence')
                set(gca, 'XTick', 0:this.interval:t, 'FontSize', 18);
                grid on
                subplot(212);
                hold on;
                plot(this.vergerr_hist, '.-', 'LineWidth', 2);
                xlabel('Iteration #', 'FontSize', 14);
                ylabel('Vergence Error [deg]', 'FontSize', 14);
                set(gca, 'XTick', 0:this.interval:t, 'FontSize', 18);
                grid on
                title(['mean error = ' num2str(merr) 'deg   \sigma = ' num2str(serr)]);

                %DESIRED vs ACTUAL FIXATION DEPTH
                figure;
                hold on
                plot(this.Z, 'r.-', 'LineWidth', 2, 'MarkerSize', 18);
                xlabel('Iteration #', 'FontSize', 14);
                ylabel('Depth [m]', 'FontSize', 14);
                plot(this.fixZ, 'b.-', 'LineWidth', 2, 'MarkerSize', 18);
                xlabel('Iteration #', 'FontSize', 14);
                grid on
                legend('Object Depth', 'Fixation Depth')

                %CONTROL ERROR (shows the command and the actual vergence)
                %display actual vergence and vergence command at each time step
                figure;
                subplot(2, 1, 1);
                hold on
                plot(this.verge_actual(1:end), 'r.-', 'LineWidth', 2, 'MarkerSize', 18);
                xlabel('Trial #', 'FontSize', 14);
                ylabel('Vergence angle [deg]', 'FontSize', 14);
                plot(this.verge_command(1:end), 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
                grid on;
                legend('Actual', 'Vergence Command')
                subplot(2, 1, 2);
                hold on
                plot(this.verge_actual - this.verge_command, 'g-.', 'LineWidth', 2)
                set(gca, 'XTick', 0:this.interval:this.T_train, 'FontSize', 14);
                legend('Control Error')
            end
        end
    end
end

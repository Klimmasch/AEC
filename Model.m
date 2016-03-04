classdef Model < handle
    properties
        scmodel_Large;      %SparseCoding class (downsampled version)
        scmodel_Small;      %SparseCoding class
        rlmodel;            %ReinforcementLearning class

        interval;           %period to change a new environment for the eye

        learnedFile;        %file with config
        textureFile;        %texture file that was used
        trainTime;          %train
        simulated_time;     %how long did the model take to be learned (min)
        sparseCodingType;   %type of sparse coding

        stopSC = 0;
        SCInterval = 1;

        lambdaMuscleFB;     %factor of musscle activity feedback to RL feature vector
        lambdaMet;          %proportion of recErr and metCosts for reward function
        lambdaRec;          %reconstruction error factor
        lambdaV;            %value networks input->output weights factor | L1 norm 0.1 | L2 norm 0.04
        lambdaP1;           %policy networks input->hidden weights factor | L1 norm 0.1 | L2 norm 3.2
        lambdaP2;           %policy networks hidden->output weights factor | L1 norm 1.2 | L2 norm 430.8

        %model data
        recerr_hist;        %history of rec error
        disp_hist;          %history of disparity
        vergerr_hist;       %history of vergence error
        verge_actual;       %actual vergence angle
        verge_desired;      %desired vergence angle
        Z;                  %object depth
        fixZ;               %depth of fixation
        g_hist;             %history of nat gradient change
        td_hist;            %history of td error
        feature_hist;       %history of feature vector
        cmd_hist;           %history of vergence commands
        AC_norm_weights;    %history of norm of the weights of the actor and critic
        l12_weights;        %history of L1/L2, i.e. sum abs, sum pow2 weights of the actor and critic
        relCmd_hist;        %relativ changes in motor commands
        reward_hist;        %reward function
        metCost_hist;       %metabolic costs
    end

    methods
        function obj = Model(PARAM)
            obj.learnedFile = PARAM{1}{1};
            obj.textureFile = PARAM{1}{2};
            obj.trainTime = PARAM{1}{3};
            obj.sparseCodingType = PARAM{1}{4};
            obj.interval = PARAM{1}{5};
            obj.lambdaMuscleFB = PARAM{1}{6};
            obj.lambdaMet = PARAM{1}{7};
            obj.lambdaRec = PARAM{1}{8};
            obj.lambdaV = PARAM{1}{9};
            obj.lambdaP1 = PARAM{1}{10};
            obj.lambdaP2 = PARAM{1}{11};

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

            obj.recerr_hist = zeros(obj.trainTime, 2);     %coarse and fine scale error
            obj.disp_hist = zeros(obj.trainTime, 1);       %saved disparity values
            obj.vergerr_hist = zeros(obj.trainTime, 1);    %saved vergence values
            obj.verge_actual = zeros(obj.trainTime, 1);    %saved vergence values
            obj.verge_desired = zeros(obj.trainTime, 1);   %saved vergence values
            obj.Z = zeros(obj.trainTime, 1);               %saved object depth
            obj.fixZ = zeros(obj.trainTime, 1);            %saved fixation depth

            obj.g_hist = zeros(obj.trainTime, 1);          %history of nat gradient change
            obj.td_hist = zeros(obj.trainTime, 1);         %history of td error
            % obj.feature_hist = zeros(obj.trainTime, obj.rlmodel.S0);    %history of feature vector ##!
            obj.cmd_hist = zeros(obj.trainTime, 2);        %history of vergence commands
            obj.relCmd_hist = zeros(obj.trainTime, 1);
            obj.AC_norm_weights = zeros(obj.trainTime, 7);
            obj.l12_weights = zeros(obj.trainTime, 8);
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

            imPind = find(sum(imagesLarge .^ 2));  %find non-zero patches (columns)
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
            % feature_L = mean(Coef_L(:, imPind).^2, 2) * 1; %feature vector for large scale (*5), average activation, Eq. 3.1
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
            % feature_S = mean(Coef_S(:, imPind).^2, 2)*1;  %feature vector for small scale
            feature_S = mean(Coef_S(:, imPind) .^ 2, 2);  %feature vector for small scale

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

        %% Plotting errors and save graphs
        function allPlotSave(this, savePath)
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
            plotpath = sprintf('%s/mvngAvgVergErr', savePath);
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
            xlabel(sprintf('Iteration # (windowSize=%d)', windowSize * this.interval), 'FontSize', 12);
            ylabel('RMSE Vergence Error [deg]', 'FontSize', 12);
            title('RMSE of the Vergence Error');
            plotpath = sprintf('%s/rmseVergErr', savePath);
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
            plotpath = sprintf('%s/recErr', savePath);
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
            plotpath = sprintf('%s/vergenceAngle', savePath);
            saveas(gcf, plotpath, 'png');

            %% Muscel graphs
            figure;
            hold on;
            grid on;
            subplot(3, 1, 1);
            plot(this.cmd_hist(:, 2), 'color', [rand, rand, rand], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            title('Absolute Muscel Commands');

            subplot(3, 1, 2);
            plot(this.relCmd_hist, 'color', [rand, rand, rand], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            title('Relative Muscel Commands');

            subplot(3, 1, 3);
            plot(this.metCost_hist, 'color', [rand, rand, rand], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            title('Metabolic Costs');

            plotpath = sprintf('%s/muscelGraphs', savePath);
            saveas(gcf, plotpath, 'png');

            %% Weights
            figure;
            hold on;
            grid on;
            plot(this.l12_weights(:, 1), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
            plot(this.l12_weights(:, 3), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
            plot(this.l12_weights(:, 5), 'color', [1, 0.5098, 0.1961], 'LineWidth', 1.3);
            plot(this.l12_weights(:, 7), 'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('\Sigma \midweights\mid', 'FontSize', 12);
            legend('w_{Vji}', 'w_{Pji}', 'w_{Pkj}', 'w_{Pnji}');
            title('Model weights (L1)')
            plotpath = sprintf('%s/weightsL1', savePath);
            saveas(gcf, plotpath, 'png');

            figure;
            hold on;
            grid on;
            plot(this.l12_weights(:, 2), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
            plot(this.l12_weights(:, 4), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
            plot(this.l12_weights(:, 6), 'color', [1, 0.5098, 0.1961], 'LineWidth', 1.3);
            plot(this.l12_weights(:, 8), 'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('\Sigma weights^{2}', 'FontSize', 12);
            legend('w_{Vji}', 'w_{Pji}', 'w_{Pkj}', 'w_{Pnji}');
            title('Model weights (L2)')
            plotpath = sprintf('%s/weightsL2', savePath);
            saveas(gcf, plotpath, 'png');

            %% Reward
            figure;
            hold on;
            grid on;
            r = [- this.lambdaMet * this.metCost_hist, ...
                 - this.lambdaP2 * this.l12_weights(:, 5), ...
                 - this.lambdaP1 * this.l12_weights(:, 3), ...
                 - this.lambdaV * this.l12_weights(:, 1), ...
                 - this.lambdaRec * (this.recerr_hist(:, 1) + this.recerr_hist(:, 2))];
            area(r, 'LineStyle','none');
            % TODO: function is delayed by 10 iteration steps
            % plot(this.reward_hist, 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            % legend('\lambdametCost', '\lambdaL1(w_{Pkj})', '\lambdaL1(w_{Pji})', '\lambdaL1(w_{Vji})', '\lambdaRecErr', 'Reward');
            legend('\lambdametCost', '\lambdaL1(w_{Pkj})', '\lambdaL1(w_{Pji})', '\lambdaL1(w_{Vji})', '\lambdaRecErr');
            title('Reward composition (L1)');
            plotpath = sprintf('%s/rewardCompL1', savePath);
            saveas(gcf, plotpath, 'png');
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

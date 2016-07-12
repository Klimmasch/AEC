classdef Model < handle
    properties
        scModel;            %SparseCoding cell array containing respective classes
        rlModel;            %ReinforcementLearning class

        focalLength;        %focal length [px]
        baseline;           %interocular distance
        objDistMin;         %object distance to eyes [m]
        objDistMax;
        muscleInitMin;      %minimal initial muscle innervation
        muscleInitMax;      %maximal --"--
        interval;           %period of eye stimulus change
        desiredAngleMin;    %min/max desired vergence angle
        desiredAngleMax;
        fixDistMin;
        fixDistMax;
        vergAngleMin;
        vergAngleMax;
        vergAngleFixMin;
        vergAngleFixMax;

        textureFile;        %config file containing texture stimulus list
        trainTime;          %number of training (intended) iterations
        trainedUntil;       %how long did the training actually proceed?
        simulatedTime;      %how long did the model take to be learned (min)

        sparseCodingType;   %type of sparse coding

        lambdaMuscleFB;     %factor of muscle activity feedback to RL feature vector
        lambdaRec;          %reconstruction error factor
        lambdaMet;          %factor of metCosts for reward function

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
        weight_hist;        %L1/L2, i.e. sum abs, sum pow2 weights of actor and critic
        reward_hist;        %reward function
        metCost_hist;       %metabolic costs
        variance_hist;      %exploratory variance of actor
        savePath;           %where all the data are stored
        notes;              %is there any special things about this model to note?

        % Model results at testing procedure
        disZtest;
        fixZtest;
        vergErrTest;
        responseResults;
        testResult;
        testResult2;
        testResult3;
        testResult4;
        testResult5;

        % Image processing
        patchSize;
        pxFieldOfView
        dsRatio;
        columnInd;
        stride;
        patchesLeft;
        patchesRight;
        overlap;
        cutout;
    end

    methods
        function obj = Model(PARAM)
            obj.textureFile = PARAM{1}{1};
            obj.trainTime = PARAM{1}{2};
            obj.sparseCodingType = PARAM{1}{3};
            obj.focalLength = PARAM{1}{4};
            obj.baseline = PARAM{1}{5};
            obj.objDistMin = PARAM{1}{6};
            obj.objDistMax = PARAM{1}{7};
            obj.muscleInitMin = PARAM{1}{8};
            obj.muscleInitMax = PARAM{1}{9};
            obj.interval = PARAM{1}{10};
            obj.lambdaMuscleFB = PARAM{1}{11};
            obj.lambdaRec = PARAM{1}{12};
            obj.lambdaMet = PARAM{1}{13};
            obj.fixDistMin = PARAM{1}{18};
            obj.fixDistMax = PARAM{1}{19};

            % single eye
            obj.desiredAngleMin = atand(obj.baseline / (2 * obj.objDistMax));
            obj.desiredAngleMax = atand(obj.baseline / (2 * obj.objDistMin));

            obj.vergAngleMin = 2 * atand(obj.baseline / (2 * obj.fixDistMax));
            obj.vergAngleMax = 2 * atand(obj.baseline / (2 * obj.fixDistMin));

            % obj.vergAngleFixMin = 2 * atand(obj.baseline / (2 * 2));
            obj.vergAngleFixMin = 2 * atand(obj.baseline / (2 * obj.objDistMax));
            obj.vergAngleFixMax = 2 * atand(obj.baseline / (2 * obj.objDistMin));

            %%% Create RL models
            % Discrete or continuous policy
            if (PARAM{3}{11} == 1)
                obj.rlModel = ReinforcementLearningCont(PARAM{3});
            else
                obj.rlModel = ReinforcementLearning(PARAM{3});
            end

            obj.recerr_hist = zeros(obj.trainTime, length(PARAM{2}{1})); % recerr_hist = t x #SC_scales
            obj.disp_hist = zeros(obj.trainTime, 1);
            obj.vergerr_hist = zeros(obj.trainTime, 1);
            obj.verge_actual = zeros(obj.trainTime, 1);
            obj.verge_desired = zeros(obj.trainTime, 1);
            obj.Z = zeros(obj.trainTime, 1);
            obj.fixZ = zeros(obj.trainTime, 1);

            obj.g_hist = zeros(obj.trainTime, 1);
            obj.td_hist = zeros(obj.trainTime, 1);
            % obj.feature_hist = zeros(obj.trainTime, PARAM{3}{9}(1));
            obj.cmd_hist = zeros(obj.trainTime, 2);
            obj.relCmd_hist = zeros(obj.trainTime, PARAM{3}{9}(3)); % relCmd_hist = t x output_dim
            % obj.weight_hist = zeros(obj.trainTime, 4);
            obj.weight_hist = zeros(obj.trainTime, 6); % for also traking change in weights
            obj.reward_hist = zeros(obj.trainTime, 1);
            obj.metCost_hist = zeros(obj.trainTime, 1);
            obj.variance_hist = zeros(obj.trainTime, 1);

            obj.disZtest = [];
            obj.fixZtest = [];
            obj.vergErrTest = [];
            obj.responseResults = struct();
            obj.testResult = [];
            obj.testResult2 = [];
            obj.testResult3 = [];
            obj.testResult4 = [];
            obj.testResult5 = [];
            obj.simulatedTime = 0;
            obj.trainedUntil = 0;
            obj.notes = '';

            %%% Generate image processing constants
            obj.patchSize = PARAM{1}{14};
            obj.pxFieldOfView = PARAM{1}{15};
            obj.dsRatio = PARAM{1}{16};
            obj.stride = PARAM{1}{17};
            obj.overlap = PARAM{1}{20};
            obj.cutout = PARAM{1}{21};

            % Prepare index matrix for image patches
            obj.prepareColumnInd();
            % cut out central region
            if (obj.cutout == 1)
                obj.prepareCutout();
            end

            %%% Create SC models
            obj.scModel = {};
            PARAM{2}{end + 1} = [cellfun('length', obj.columnInd)]; % append image batch size's 2nd dimensions
            if (obj.sparseCodingType == 0)
                for i = 1 : length(PARAM{2}{1})
                    % pick respective parameters from PARAM cell array for ith SC model constructor
                    obj.scModel{end + 1} = SparseCoding2(cellfun(@(x) x(i), PARAM{2}, 'UniformOutput', true));
                end
            else
                %TODO: update depricated SparseCodingHomeo class
                sprintf('SparseCodingHomeo class is DEPRICATED and therefore currently not supported!')
                return;
                % obj.scModel_Large = SparseCodingHomeo(PARAM{2}{1}); %coarse scale
                % obj.scModel_Small = SparseCodingHomeo(PARAM{2}{2}); %fine scale
            end

            % Intermediate patch matricies
            obj.patchesLeft = cell(1, length(obj.scModel));
            obj.patchesRight = cell(1, length(obj.scModel));
            for i = 1 : length(obj.scModel)
                obj.patchesLeft{i} = zeros(obj.patchSize ^ 2, length(obj.columnInd{i}));
                obj.patchesRight{i} = zeros(obj.patchSize ^ 2, length(obj.columnInd{i}));
            end
        end

        %%% Copy constructor
        % Make a (deep) copy of a handle object
        function new = copy(this)
            % Instantiate new object of the same class
            new = feval(class(this));

            % Copy all non-hidden properties
            p = properties(this);
            % Copy hidden properties which don't show up in the result
            % of properties(this)
            % p = fieldnames(struct(this));
            for i = 1 : length(p)
                new.(p{i}) = this.(p{i});
            end
        end

        % Generates index matrix for image patches
        % Filled image batch, i.e. all patches for the respective scale are used
        function prepareColumnInd(this)
            % index matrix
            this.columnInd = {};
            % (start) index of left upper corner of respective patches
            npc = this.pxFieldOfView - this.patchSize + 1;
            % for #scales
            for i = 1 : length(this.dsRatio)
                k = 1;
                % number of patches per column/row
                l = length(1 : this.stride(i) : npc(i));
                % each scale has l * l patches
                this.columnInd{end + 1} = zeros(1, l ^ 2);
                % calculate index of left upper corner of respective patches
                for j = 1 : this.stride(i) : npc(i)
                    this.columnInd{i}((k - 1) * l + 1 : k * l) = (j - 1) * npc(i) + 1 : this.stride(i) : j * npc(i);
                    k = k + 1;
                end
            end
        end

        % Cuts out smaller scales from larger scales' fields of view
        function prepareCutout(this)
            % most inner layer does not need a cutout
            % therefore n-1 cuts per n layers
            for i = 1 : (length(this.dsRatio) - 1)
                % Calculate cutout area, upsample to orig and downsample to current (substract more if needed)
                foveaSmallerAdjusted = ceil((this.pxFieldOfView(i + 1) - 2 * this.overlap(i)) * (this.dsRatio(i + 1) / this.dsRatio(i)));

                % Calculate offset of coarse scale [px]
                offsetCoarse = (this.pxFieldOfView(i) - foveaSmallerAdjusted) / 2;

                % # stride applications =^ offset width
                offsetWidth = ceil(((offsetCoarse - this.patchSize) / this.stride(i)) + 1);

                % number of patches per row/column
                nprc = floor((this.pxFieldOfView(i) - this.patchSize) / this.stride(i) + 1);

                % remaining patches vector [indices]
                remainder = [];

                % full rows + first left border part
                startPart = 1 : (offsetWidth * (nprc + 1));
                % remainder(1 : length(startPart)) = startPart;
                remainder = [remainder, startPart];

                % right border + next row's left border patches
                for j = 1 : (nprc - 2 * offsetWidth - 1)
                    borderPart = (nprc * (offsetWidth + j) - offsetWidth + 1) : (nprc * (offsetWidth + j) + offsetWidth);
                    % remainder(length(startPart) + length(borderPart) * (j - 1) + 1 : length(startPart) + length(borderPart) * j) = borderPart;
                    remainder = [remainder, borderPart];
                end

                % right border + remaining full rows
                endPart = (nprc * (offsetWidth + j + 1) - offsetWidth + 1) : nprc ^ 2;
                % remainder(length(startPart) + length(borderPart) * j + 1 : end) = endPart;
                remainder = [remainder, endPart];

                % Cutting columnInd to eradicate the inner patches in the outer layer
                this.columnInd{i} = this.columnInd{i}(:, remainder);
            end
        end

        %%% Patch generation
        % img:      image to be processed
        % scScale:  SC scale index elem {coarse, ..., fine}
        % eyePos:   eye position index elem {1 := left, 2 := right}
        function preprocessImage(this, img, scScale, eyePos)
            % down scale image
            for k = 1 : log2(this.dsRatio(scScale))
                img = impyramid(img, 'reduce');
            end

            % convert to double
            img = double(img);

            % cut fovea in the center
            [h, w, ~] = size(img);
            img = img(fix(h / 2 + 1 - this.pxFieldOfView(scScale) / 2) : fix(h / 2 + this.pxFieldOfView(scScale) / 2), ...
                      fix(w / 2 + 1 - this.pxFieldOfView(scScale) / 2) : fix(w / 2 + this.pxFieldOfView(scScale) / 2));

            % cut patches and store them as col vectors
            patches = im2col(img, [this.patchSize, this.patchSize], 'sliding'); % slide window of 1 px

            % take patches by application of respective strides (8 px)
            patches = patches(:, this.columnInd{scScale});

            % pre-processing steps (0 mean, unit norm)
            patches = patches - repmat(mean(patches), [size(patches, 1), 1]);   % 0 mean
            normp = sqrt(sum(patches .^ 2));                                      % patches norm

            % normalize patches to norm 1
            normp(normp == 0) = eps;                                            % regularizer
            patches = patches ./ repmat(normp, [size(patches, 1), 1]);          % normalized patches

            if (eyePos == 1)
                this.patchesLeft{scScale} = patches;
            else
                this.patchesRight{scScale} = patches;
            end
        end

        %%% Generate Feature Vector and Reward
        function [totalFeature, totalReward, errorArray] = generateFR(this, imageBatch)
            errorArray = zeros(1, length(this.scModel));
            rewardArray = zeros(1, length(this.scModel));
            featureArray = cell(1, 2);

            for i = 1 : length(this.scModel)
                tmpImages = imageBatch{i};
                imPind = find(sum(tmpImages .^ 2)); %find non-zero patches (columns)
                if (isempty(imPind))
                    % try
                    %     feature_L = zeros(this.scModel_Large.Basis_num, 1);
                    %     reward_L = this.rlModel.J; %TODO: not supported for continuous rlModel
                    % catch
                    % end
                    sprintf('All values in the extracted patches are zero. Check if the image rendering is all right!')
                    return;
                end
                this.scModel{i}.sparseEncode(tmpImages);
                errorArray(i) = sum(sum(this.scModel{i}.currentError .^ 2)) / sum(sum(tmpImages .^ 2));
                rewardArray(i) = -errorArray(i);
                % feature vector, i.e. average activation, Eq. 3.1
                featureArray{i} = mean(this.scModel{i}.currentCoef(:, imPind) .^ 2, 2);
            end

            % Compute total reconstruction error and reward
            % totalError = sum(errorArray);             % total reconstruction error
            totalReward = sum(rewardArray);             % sum rewards
            totalFeature = vertcat(featureArray{:});    % join feature vectors
        end

        %% Plotting everything and save graphs
        function allPlotSave(this)
            % windowSize = 125;
            windowSize = 1000;
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
            title('Moving Average of Vergence Error');
            legend('|verg_{err}|', 'SMA(|verg_{err}|)');
            plotpath = sprintf('%s/mvngAvgVergErrFull', this.savePath);
            saveas(gcf, plotpath, 'png');

            figure;
            hold on;
            grid on;
            % Simple Moving Average
            plot((windowSize + 1) * this.interval : this.interval : size(this.vergerr_hist), vergerr(windowSize + 1 : end), ...
                 'LineWidth', 1.3);

            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('SMA(|verg_{err}|) [deg]', 'FontSize', 12);
            title('Moving Average of Vergence Error');
            plotpath = sprintf('%s/mvngAvgVergErr', this.savePath);
            saveas(gcf, plotpath, 'png');

            %% Root Mean Squared Error
            % windowSize = 125;
            windowSize = 250;
            vergerr = this.vergerr_hist(ind);
            rmse = zeros(length(1 : windowSize : length(vergerr) - mod(length(vergerr), 2)), 1); %cut if odd length
            k = 1 : windowSize : length(vergerr) - mod(length(vergerr), 2);
            for i = 1 : length(rmse)
                try
                    rmse(i) = sqrt(mean(vergerr(k(i) : k(i) + windowSize - 1) .^ 2));
                catch
                    % if windowsize > values in vergerr
                    rmse(i) = 0;
                end
            end

            try
                figure;
                hold on;
                grid on;
                plot(windowSize * this.interval : windowSize * this.interval : (length(vergerr) - mod(length(vergerr), 2)) * this.interval, ...
                     rmse, 'LineWidth', 1.3);
                axis([-inf, inf, 0, inf]);
                xlabel(sprintf('Iteration # (windowSize=%d)', windowSize * this.interval), 'FontSize', 12);
                ylabel('RMSE(verg_{err}) [deg]', 'FontSize', 12);
                title('RMSE of Vergence Error');
                plotpath = sprintf('%s/rmseVergErr', this.savePath);
                saveas(gcf, plotpath, 'png');
            catch
                % if windowsize > values in vergerr
                sprintf('Warning: windowSize >= vergerr')
            end

            %% Reconstruction Error
            try
                figure;
                hold on;
                grid on;
                handleArray = zeros(1, length(this.scModel));
                for i = 1 : length(this.scModel)
                    tmpError = filter(ones(1, windowSize) / windowSize, 1, this.recerr_hist(ind, i));
                    handleArray(i) = plot((windowSize + 1) * this.interval : this.interval : size(this.recerr_hist, 1), tmpError(windowSize + 1 : end), ...
                                          'color', [rand, rand, rand], 'LineWidth', 1.3);
                end
                xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
                ylabel('Reconstruction Error [AU]', 'FontSize', 12);
                captions = cell(1, length(handleArray));
                for i = 1 : length(handleArray)
                    captions{i} = strcat('scale', sprintf(' %d', i));
                end
                l = legend(handleArray, captions);
                plotpath = sprintf('%s/recErr', this.savePath);
                saveas(gcf, plotpath, 'png');
            catch
                % if windowsize > values in recerr_hist
                sprintf('Warning: windowSize >= recerr_hist')
            end

            %% Vergence angle
            figure;
            hold on;
            grid on;
            if length(this.verge_desired) >= 200
                plot(this.verge_desired(end-200:end), 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.3);
                plot(this.verge_actual(end-200:end), 'b', 'LineWidth', 1.3);
            else
                plot(this.verge_desired, 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.3);
                plot(this.verge_actual, 'b', 'LineWidth', 1.3);
            end
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Angle [deg]', 'FontSize', 12);
            legend('desired', 'actual');
            title('Vergence at last 200 steps of training');
            plotpath = sprintf('%s/vergenceAngle', this.savePath);
            saveas(gcf, plotpath, 'png');

            %% Muscel graphs
            if (this.rlModel.continuous == 1)
                % Lateral Rectus
                if (this.rlModel.CActor.output_dim == 2)
                    figure;
                    hold on;
                    grid on;
                    subplot(3, 1, 1);
                    plot(this.cmd_hist(ind, 1), 'color', [rand, rand, rand], 'LineWidth', 1.3);
                    ylabel('Value', 'FontSize', 12);
                    title('Total Muscle Commands (lateral rectus)');

                    subplot(3, 1, 2);
                    plot(this.relCmd_hist(ind, 1), 'color', [rand, rand, rand], 'LineWidth', 1.3);
                    ylabel('Value', 'FontSize', 12);
                    title('\Delta Muscle Commands (lateral rectus)');

                    subplot(3, 1, 3);
                    plot(this.metCost_hist(ind), 'color', [rand, rand, rand], 'LineWidth', 1.3);
                    xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
                    ylabel('Value', 'FontSize', 12);
                    title('Metabolic Costs');

                    plotpath = sprintf('%s/muscleGraphsLateralRectus', this.savePath);
                    saveas(gcf, plotpath, 'png');
                end

                % Medial Rectus
                figure;
                hold on;
                grid on;
                subplot(3, 1, 1);
                plot(this.cmd_hist(ind, 2), 'color', [rand, rand, rand], 'LineWidth', 1.3);
                ylabel('Value', 'FontSize', 12);
                title('Total Muscle Commands (medial rectus)');

                subplot(3, 1, 2);
                if (this.rlModel.CActor.output_dim == 2)
                    plot(this.relCmd_hist(ind, 2), 'color', [rand, rand, rand], 'LineWidth', 1.3);
                else
                    plot(this.relCmd_hist(ind), 'color', [rand, rand, rand], 'LineWidth', 1.3);
                end
                ylabel('Value', 'FontSize', 12);
                title('\Delta Muscle Commands (medial rectus)');

                subplot(3, 1, 3);
                plot(this.metCost_hist(ind), 'color', [rand, rand, rand], 'LineWidth', 1.3);
                xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
                ylabel('Value', 'FontSize', 12);
                title('Metabolic Costs');

                plotpath = sprintf('%s/muscleGraphsMedialRectus', this.savePath);
                saveas(gcf, plotpath, 'png');

                if (this.rlModel.CActor.output_dim == 2)
                    % Muscle correlation check
                    % extract last x per cent of iterations
                    pcOffset = 0.01;
                    tmpOffset = ceil(size(this.cmd_hist, 1) * pcOffset);
                    % Total
                    figure;
                    hold on;
                    scatter(this.cmd_hist(end - tmpOffset : end, 1), this.cmd_hist(end - tmpOffset : end, 2), 5,'MarkerFaceColor',[0, 0.7, 0.7]);
                    corrl = corr(this.cmd_hist(end - tmpOffset : end, 1), this.cmd_hist(end - tmpOffset : end, 2));
                    xlabel('Lateral rectus [%]', 'FontSize', 12);
                    ylabel('Medial rectus [%]', 'FontSize', 12);
                    title(strcat('Total Muscle Commands (training)', sprintf('\nCorrelation = %1.2e at last %d iterations', corrl, tmpOffset)));
                    plotpath = sprintf('%s/muscleGraphsScatterTotalTraining', this.savePath);
                    saveas(gcf, plotpath, 'png');

                    % Delta
                    figure;
                    hold on;
                    scatter(this.relCmd_hist(end - tmpOffset : end, 1), this.relCmd_hist(end - tmpOffset : end, 2), 5,'MarkerFaceColor',[0, 0.7, 0.7]);
                    corrl = corr(this.relCmd_hist(end - tmpOffset : end, 1), this.relCmd_hist(end - tmpOffset : end, 2));
                    xlabel('Lateral rectus [%]', 'FontSize', 12);
                    ylabel('Medial rectus [%]', 'FontSize', 12);
                    title(strcat('\Delta Muscle Commands (training)', sprintf('\nCorrelation = %1.2e at last %d iterations', corrl, tmpOffset)));
                    plotpath = sprintf('%s/muscleGraphsScatterDeltaTraining', this.savePath);
                    saveas(gcf, plotpath, 'png');
                end
            end

            %% Weights
            if (this.rlModel.continuous == 1)
                % L1 norm
                figure;
                hold on;
                grid on;
                subplot(3, 1, 1);
                plot(this.weight_hist(:, 1), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
                ylabel('\Sigma \midweights\mid', 'FontSize', 12);
                title('w_{Vji}');
                subplot(3, 1, 2);
                plot(this.weight_hist(:, 2), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
                ylabel('\Sigma \midweights\mid', 'FontSize', 12);
                title('w_{Pji}');
                subplot(3, 1, 3);
                plot(this.weight_hist(:, 3), 'color', [1, 0.5098, 0.1961], 'LineWidth', 1.3);
                xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
                ylabel('\Sigma \midweights\mid', 'FontSize', 12);
                title('w_{Pkj}');
                plotpath = sprintf('%s/weightsDevelopment', this.savePath);
                saveas(gcf, plotpath, 'png');
            else
                % L1 norm
                figure;
                hold on;
                grid on;
                subplot(2, 1, 1);
                plot(this.weight_hist(:, 1), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
                ylabel('\Sigma \midweights\mid', 'FontSize', 12);
                title('w_{Vji}');
                subplot(2, 1, 2);
                plot(this.weight_hist(:, 2), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
                ylabel('\Sigma \midweights\mid', 'FontSize', 12);
                xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
                title('w_{Pji}');
                plotpath = sprintf('%s/weightsDevelopment', this.savePath);
                saveas(gcf, plotpath, 'png');
            end

            %% Reward composition
            if (this.rlModel.continuous == 1)
                figure;
                hold on;
                grid on;
                r = [- this.lambdaMet * this.metCost_hist, ...
                     - this.lambdaRec * sum(this.recerr_hist, 2)];
                handle = area(r, 'LineStyle','none');
                xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
                ylabel('Value', 'FontSize', 12);
                % % l = legend('\lambdametCost', '\lambdaL1(w_{Pkj})', '\lambdaL1(w_{Pji})', '\lambdaL1(w_{Vji})', '\lambdaRecErr');
                l = legend('\lambdametCost', '\lambdaRecErr');
                if(version('-release') == '2015b')
                    handle(1).FaceColor = [1, 0.25, 0];
                    handle(2).FaceColor = [1, 0.549, 0];
                    l.Location = 'southwest';
                end
                % % title('Reward composition (L1)');
                title('Reward composition');
                plotpath = sprintf('%s/rewardComp', this.savePath);
                saveas(gcf, plotpath, 'png');

                %% Total reward
                % figure;
                % hold on;
                % grid on;
                % plot(this.reward_hist, 'color', [1, 0.25, 0]);
                % xlabel('Iteration #', 'FontSize', 12);
                % ylabel('Value', 'FontSize', 12);
                % title('Reward');
                % plotpath = sprintf('%s/rewardTotal', this.savePath);
                % saveas(gcf, plotpath, 'png');
            else
                figure;
                hold on;
                grid on;
                plot(-sum(this.recerr_hist, 2), 'color', [1, 0.549, 0]);
                xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
                ylabel('Value', 'FontSize', 12);
                title('Reward');
                plotpath = sprintf('%s/rewardHistory', this.savePath);
                saveas(gcf, plotpath, 'png');
            end

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

            % actualResponse = [this.vergerr_hist, this.relCmd_hist];

            % % observation Window, i.e. plot statistics over last #obsWin iterations
            % obsWin = 1000;
            % if (size(actualResponse, 1) <= obsWin)
            %     obsWin = size(actualResponse, 1) - 1;
            % end
            % nVal = 20; % #bins of statistics

            % actualResponse = sortrows(actualResponse(end - obsWin + 1 : end, :));
            % if (all(actualResponse == 0)) % skip plot if empty or still training
            %     return;
            % end
            % deltaVergErr = (abs(actualResponse(1, 1)) + abs(actualResponse(end, 1))) / nVal;
            % % actualResponseStat = nVal x 3 = [index_x = vergence_error angle, mean_muscle_force, std_muscle_force]
            % actualResponseStat = zeros(nVal, 3);

            % for i = 1:nVal
            %     actualResponseStat(i, 1) = (actualResponse(1, 1) - deltaVergErr / 2) + i * deltaVergErr;
            %     actualResponseStat(i, 2) = mean(actualResponse(find(actualResponse(:, 1) >= actualResponse(1, 1) + (i - 1) * deltaVergErr ...
            %                                     & actualResponse(:, 1) <= actualResponse(1, 1) + i * deltaVergErr), 2));
            %     actualResponseStat(i, 3) = std(actualResponse(find(actualResponse(:, 1) >= actualResponse(1, 1) + (i - 1) * deltaVergErr ...
            %                                    & actualResponse(:, 1) <= actualResponse(1, 1) + i * deltaVergErr), 2));
            % end
            % actualResponseStat(isnan(actualResponseStat(:, 2)), :) = []; % drop NaN elements

            % degrees = load('Degrees.mat');
            % resolution = 10001;
            % approx = spline(1:11, degrees.results_deg(:, 1));

            % xValPos = ppval(approx, 1:0.001:11)';
            % yValPos = linspace(0, 1, resolution)';

            % xValNeg = flipud(ppval(approx, 1:0.001:11)' * -1);
            % yValNeg = linspace(-1, 0, resolution)';

            % % calculate muscle function :=  mf(vergence_angle) = muscle force [single muuscle]
            % mf = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
            % dmf = diff(mf(1:2, 1)); % delta in angle
            % indZero = find(mf(:, 2) == 0); % MF == 0_index
            % indMaxFix = find(mf(:, 1) <= this.desiredAngleMin + dmf & mf(:, 1) >= this.desiredAngleMin - dmf); % MF(desiredAngleMin)_index
            % indMinFix = find(mf(:, 1) <= this.desiredAngleMax + dmf & mf(:, 1) >= this.desiredAngleMax - dmf); % MF(desiredAngleMax)_index

            % % perfect_response := [max_fixation_x, max_fixation_y, min_fixation_x, min_fixation_y]
            % % x = vergenceError, y = deltaMuscelForce
            % perfectResponseMaxFix = [(mf(indMaxFix, 1) - flipud(mf(indMaxFix : end, 1))) * 2, ...
            %                          (mf(indMaxFix, 2) - flipud(mf(indMaxFix : end, 2))); ...
            %                          (mf(indMaxFix, 1) - flipud(mf(indZero : indMaxFix - 1, 1))) * 2, ...
            %                          (mf(indMaxFix, 2) - flipud(mf(indZero : indMaxFix - 1, 2)))];

            % perfectResponseMinFix = [(mf(indMinFix, 1) - flipud(mf(indMinFix : end, 1))) * 2, ...
            %                          (mf(indMinFix, 2) - flipud(mf(indMinFix : end, 2))); ...
            %                          (mf(indMinFix, 1) - flipud(mf(indZero : indMinFix - 1, 1))) * 2, ...
            %                          (mf(indMinFix, 2) - flipud(mf(indZero : indMinFix - 1, 2)))];

            % perfectResponse = [perfectResponseMaxFix, perfectResponseMinFix];

            % figure;
            % hold on;
            % grid on;
            % % perfect response to vergence error
            % plot(perfectResponse(:, 1), perfectResponse(:, 2), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
            % plot(perfectResponse(:, 3), perfectResponse(:, 4), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
            % % error bars of actual response of the model
            % errorbar(actualResponseStat(:, 1), actualResponseStat(:, 2), actualResponseStat(:, 3),'color', [1, 0.5098, 0.1961], 'LineWidth', 0.9);
            % % actual response of the model
            % plot(actualResponseStat(:, 1), actualResponseStat(:, 2),'color', [1, 0.0784, 0], 'LineWidth', 1.3);
            % l = legend('perfect (fixDist_{max})', 'perfect (fixDist_{min})', 'actual');
            % if(version('-release') == '2015b')
            %     l.FontSize = 7;
            %     l.Orientation = 'horizontal';
            %     l.Location = 'southoutside';
            % end
            % % adjust axis to actual response ranges + std deviation
            % xmin = min(actualResponseStat(:, 1)) * 1.1;
            % xmax = max(actualResponseStat(:, 1)) * 1.1;
            % if (xmin >= xmax)
            %     xmin = -5.5;
            %     xmax = 5.5;
            % end
            % ymin = min([-0.1, (min(actualResponseStat(:, 2)) - max(actualResponseStat(:, 3))) * 1.1]);
            % ymax = max([max(perfectResponse(:, 4)) * 1.1, (max(actualResponseStat(:, 2)) + max(actualResponseStat(:, 3))) * 1.1]);
            % if (ymin >= ymax)
            %     ymin = -0.1;
            %     ymax = 0.1;
            % end
            % plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
            % plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
            % axis([xmin, xmax, ymin, ymax]);
            % xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
            % ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
            % title(strcat('\Delta MF(verg_{err}) response after ', sprintf(' %d iterations', size(this.vergerr_hist, 1) - obsWin)));
            % plotpath = sprintf('%s/deltaMFasFktVerErr', this.savePath);
            % saveas(gcf, plotpath, 'png');
        end

        %% plot & save delta_MF(Vergence_error)
        % observation Window = obsWin, i.e. plot statistics over last #obsWin iterations
        % see delta_MF(Vergence_error) @ allPlotSave() for remarks regarding calculations
        function deltaMFplotObsWin(this, obsWin)
            actualResponse = [this.vergerr_hist, this.relCmd_hist];

            if (size(actualResponse, 1) <= obsWin)
                obsWin = size(actualResponse, 1) - 1;
            end
            nVal = 20; % #bins of statistics

            actualResponse = sortrows(actualResponse(end - obsWin + 1 : end, :));
            if (all(actualResponse == 0)) % skip plot if empty or still training
                sprintf('No model response within observation window.')
                return;
            end
            deltaVergErr = (abs(actualResponse(1, 1)) + abs(actualResponse(end, 1))) / nVal;
            % actualResponseStat = nVal x 3 = [index_x = vergence_error angle, mean_muscle_force, std_muscle_force]
            actualResponseStat = zeros(nVal, 3);

            for i = 1:nVal
                actualResponseStat(i, 1) = (actualResponse(1, 1) - deltaVergErr / 2) + i * deltaVergErr;
                actualResponseStat(i, 2) = mean(actualResponse(find(actualResponse(:, 1) >= actualResponse(1, 1) + (i - 1) * deltaVergErr ...
                                                & actualResponse(:, 1) <= actualResponse(1, 1) + i * deltaVergErr), 2));
                actualResponseStat(i, 3) = std(actualResponse(find(actualResponse(:, 1) >= actualResponse(1, 1) + (i - 1) * deltaVergErr ...
                                               & actualResponse(:, 1) <= actualResponse(1, 1) + i * deltaVergErr), 2));
            end
            actualResponseStat(isnan(actualResponseStat(:, 2)), :) = []; % drop NaN elements

            degrees = load('Degrees.mat');
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
            title(strcat('\Delta MF(verg_{err}) response after ', sprintf(' %d iterations', size(this.vergerr_hist, 1) - obsWin)));
            plotpath = sprintf('%s/deltaMFasFktVerErrObsWin', this.savePath);
            saveas(gcf, plotpath, 'png');
        end

        %% plot & save delta_MF(Vergence_error)
        % startIter, endIter = plot statistics over #startIter : endIter iterations
        % see delta_MF(Vergence_error) @ allPlotSave() for remarks regarding calculations
        function deltaMFplotStartEnd(this, startIter, endIter)
            actualResponse = [this.vergerr_hist, this.relCmd_hist];

            if (endIter > size(actualResponse, 1))
                endIter = size(actualResponse, 1);
            end
            nVal = 20; % #bins of statistics

            actualResponse = sortrows(actualResponse(startIter : endIter, :));
            if (all(actualResponse == 0)) % skip plot if empty or still training
                sprintf('No model response within observation window.')
                return;
            end
            deltaVergErr = (abs(actualResponse(1, 1)) + abs(actualResponse(end, 1))) / nVal;
            % actualResponseStat = nVal x 3 = [index_x = vergence_error angle, mean_muscle_force, std_muscle_force]
            actualResponseStat = zeros(nVal, 3);

            for i = 1:nVal
                actualResponseStat(i, 1) = (actualResponse(1, 1) - deltaVergErr / 2) + i * deltaVergErr;
                actualResponseStat(i, 2) = mean(actualResponse(find(actualResponse(:, 1) >= actualResponse(1, 1) + (i - 1) * deltaVergErr ...
                                                & actualResponse(:, 1) <= actualResponse(1, 1) + i * deltaVergErr), 2));
                actualResponseStat(i, 3) = std(actualResponse(find(actualResponse(:, 1) >= actualResponse(1, 1) + (i - 1) * deltaVergErr ...
                                               & actualResponse(:, 1) <= actualResponse(1, 1) + i * deltaVergErr), 2));
            end
            actualResponseStat(isnan(actualResponseStat(:, 2)), :) = []; % drop NaN elements

            degrees = load('Degrees.mat');
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

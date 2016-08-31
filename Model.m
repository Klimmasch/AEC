classdef Model < handle
    properties
        scModel;            % SparseCoding cell array containing respective classes
        rlModel;            % ReinforcementLearning class

        focalLength;        % focal length [px]
        baseline;           % interocular distance
        objDistMin;         % object distance to eyes [m]
        objDistMax;
        muscleInitMin;      % minimal initial muscle innervation
        muscleInitMax;      % maximal --"--
        interval;           % period of eye stimulus change
        desiredAngleMin;    % min/max desired vergence angle
        desiredAngleMax;
        fixDistMin;
        fixDistMax;
        vergAngleMin;
        vergAngleMax;
        vergAngleFixMin;
        vergAngleFixMax;

        textureFile;        % config file containing texture stimulus list
        trainTime;          % number of training (intended) iterations
        trainedUntil;       % how long did the training actually proceed?
        simulatedTime;      % how long did the model take to be learned (min)

        sparseCodingType;   % type of sparse coding

        lambdaMuscleFB;     % factor of muscle activity feedback to RL feature vector
        lambdaRec;          % reconstruction error factor
        lambdaMet;          % factor of metCosts for reward function

        % Model data history
        recerr_hist;        % reconstruction error [coarse scale, fine scale]
        disp_hist;          % disparity
        vergerr_hist;       % vergence error
        verge_actual;       % actual vergence angle
        verge_desired;      % desired vergence angle
        Z;                  % object depth
        fixZ;               % fixation depth
        g_hist;             % natural gradient change
        td_hist;            % temporal difference (td) error
        feature_hist;       % feature vector
        cmd_hist;           % vergence commands
        relCmd_hist;        % relativ changes in motor commands
        weight_hist;        % L1/L2, i.e. sum abs, sum pow2 weights of actor and critic
        reward_hist;        % reward function
        metCost_hist;       % metabolic costs
        variance_hist;      % exploratory variance of actor
        savePath;           % where all the data are stored
        notes;              % is there any special things about this model to note?

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

        % Providing relevant methods
        degrees;        % contains a tabular that maps muscle activity to angles
        degreesIncRes;  % contains a part of degrees with increased resolution
        degDiff;        % difference between entries in degreesIncRes
        metCosts;       % contains a tabular that maps muscle activity to metabolic costs
        mfunctionMR;    % contains a function approx. of the mapping from activity to angle
                        % for the medial rectus with activation of lateral rectus equals 0
        mfunctionLR;    % analogous for the lateral rectus
        dmf;            % delta in muscle force (for mfunctions) dmf2
        dAngleMR;       % delta in angle (for mfunctions) dmf
        dAngleLR;       % delta in muscle force (for mfunctionLR) dmf3
        angleMin;       % minimal and maximal angle that can be reached by one-dimensional muscle commands
        angleMax;

        imgRawLeft;     % images that are updated by refreshImages
        imgRawRight;
        imgGrayLeft;
        imgGrayRight;
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
            % obj.feature_hist = zeros(obj.trainTime, 1);%zeros(obj.trainTime, PARAM{3}{9}(1));
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

            %%% Preparing Functions
            % getMF functions
            obj.degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
            obj.metCosts = load('MetabolicCosts.mat');

            resFactor = 10; % factor by which the resolution of the tabular should be increased
            % resFactor = 13; % results in comparable amount of entries as the mfunctions
            %resFac, size of tabular: 5,33 | 6,65 | 10,1025 | 13,8193 | x,(2^x)+1
            % between 0 and 0.1 mus. act. the resulution is increased
            obj.degreesIncRes = interp2(obj.degrees.results_deg(1:2, 1:2), resFactor);
            obj.degDiff = max(max(diff(obj.degreesIncRes)));    % distance between the entries

            % muscle function :=  mf(vergence_angle) = muscle force [single muscle]
            resolution = 100001;
            approx = spline(1 : 11, obj.degrees.results_deg(:, 1));

            xValPos = ppval(approx, 1 : 0.0001 : 11)';
            yValPos = linspace(0, 1, resolution)';

            % xValNeg = flipud(ppval(approx, 1 : 0.0001 : 11)' * -1);
            % yValNeg = linspace(-1, 0, resolution)';

            % mfunction = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
            obj.mfunctionMR = [xValPos, yValPos];
            obj.mfunctionMR(:, 1) = obj.mfunctionMR(:, 1) * 2;  % angle for two eyes
            obj.dAngleMR = abs(diff(obj.mfunctionMR(1 : 2, 1)));   % delta in angle
            obj.dmf = diff(obj.mfunctionMR(1 : 2, 2));       % delta in muscle force

            approx = spline(1 : 11, obj.degrees.results_deg(1, :));
            xValPos = ppval(approx, 1 : 0.0001 : 11)';
            yValPos = linspace(0, 1, resolution)';

            obj.mfunctionLR = [xValPos, yValPos];
            obj.mfunctionLR(:, 1) = obj.mfunctionLR(:, 1) * 2;    % angle for two eyes
            obj.dAngleLR = abs(diff(obj.mfunctionLR(1 : 2, 1)));     % delta in angle

            obj.angleMin = min(obj.mfunctionLR(obj.mfunctionLR(:, 1) > 0));
            obj.angleMax = obj.mfunctionMR(end, 1);

            % refreshImages
            obj.imgRawLeft = uint8(zeros(240, 320, 3));
            obj.imgRawRight = uint8(zeros(240, 320, 3));
            obj.imgGrayLeft = uint8(zeros(240, 320, 3));
            obj.imgGrayRight = uint8(zeros(240, 320, 3));
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
            % pixel index of left upper corner of last patch
            ilp = this.pxFieldOfView - this.patchSize + 1;

            % for #scales
            for i = 1 : length(this.dsRatio)
                k = 1;
                % number of patches per row/column
                nprc = length(1 : this.stride(i) : ilp(i));

                % each scale has nprc * nprc patches
                this.columnInd{end + 1} = zeros(1, nprc ^ 2);

                % calculate indices of patches with respect to stride
                for j = 1 : this.stride(i) : ilp(i)
                    this.columnInd{i}((k - 1) * nprc + 1 : k * nprc) = (j - 1) * ilp(i) + 1 : this.stride(i) : j * ilp(i);
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
        function preprocessImage(this, leftBool, scScale, eyePos)

            if leftBool
                img = this.imgGrayLeft;
            else
                img = this.imgGrayRight;
            end

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
            normp = sqrt(sum(patches .^ 2));                                    % patches norm

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

        %% calculating muscle inits and angles

        % maps {vergenceAngle} -> {muscleForce} for one muscle
        function mf = getMF(this, vergAngle)
            % look up index of vergAngle
            if (vergAngle >= this.degrees.results_deg(1, 1) * 2)
                indVergAngle = find(this.mfunctionMR(:, 1) <= vergAngle + this.dAngleMR & this.mfunctionMR(:, 1) >= vergAngle - this.dAngleMR);
                mf = this.mfunction(indVergAngle, 2);
                mf = mf(ceil(length(mf) / 2));
            else
                mf = 0;
            end
        end

        % Calculates muscle force for two muscles, whereas the activation
        % of on muscle is always 0.
        function [mf, angleInit] = getMF2(this, objDist, desVergErr)
            % correct vergence angle for given object distance
            angleCorrect = 2 * atand(this.baseline / (2 * objDist));
            % desired init angle for given vergence error [deg]
            angleInit = angleCorrect - desVergErr;
            % look up index of angleInit
            % if objDist not fixateable with medial rectus, use lateral rectus
            if (angleInit >= this.mfunctionMR(1, 1))
                indAngleInit = find(this.mfunctionMR(:, 1) <= angleInit + this.dAngleMR & this.mfunctionMR(:, 1) >= angleInit - this.dAngleMR);
                mf = this.mfunctionMR(indAngleInit, 2);
                mf = [0; mf(ceil(length(mf) / 2))];
            else
                indAngleInit = find(this.mfunctionLR(:, 1) <= angleInit + this.dAngleLR & this.mfunctionLR(:, 1) >= angleInit - this.dAngleLR);
                mf = this.mfunctionLR(indAngleInit, 2);
                mf = [mf(ceil(length(mf) / 2)); 0];
            end
        end

        % same as above, but now the returned muscle activations are
        % greater 0.
        % == get muscle force equally distributed over object distance.
        function [mfLR, mfMR, angleInit] = getMFedood(this, objDist, desVergErr, useBounds)
            angleCorrect = 2 * atand(this.baseline / (2 * objDist));
            angleInit = angleCorrect - desVergErr;
            % angleInit is the angle for both eyes, but degreesIncRes only
            % contains angles for one eye
            [xi, yi] = find(this.degreesIncRes <= (angleInit/2) + this.degDiff & this.degreesIncRes >= (angleInit/2) - this.degDiff);

            if ~useBounds
                i = randi(length(xi));

                % now transform indizes to muscle activities
                numInds = length(this.degreesIncRes);
                scaleFac = 0.1 / numInds;
                mfMR = xi(i) * scaleFac;
                mfLR = yi(i) * scaleFac;
            else
                mfMR = 1;
                mfLR = 1;
                while mfMR > 0.1 || mfLR > 0.015 %fist part should be unnecessary
                    i = randi(length(xi));

                    % now transform indizes to muscle activities
                    numInds = length(this.degreesIncRes);
                    scaleFac = 0.1 / numInds;
                    mfMR = xi(i) * scaleFac;
                    mfLR = yi(i) * scaleFac;
                end
            end
        end

        % mapping muscle activity to angle (one eye)
        function angle = getAngle(this, command)
            cmd = (command * 10) + 1;                                       % scale commands to table entries
            angle = interp2(this.degrees.results_deg, cmd(1), cmd(2), 'spline'); % interpolate in table by cubic splines
        end

        % depricated method: slightly faster variant, also a bit less precise,
        % only works for the medial rectus
        function angle = getAngle2(this, command)
            angleIndex = find(this.mfunctionMR(:, 2) <= command(2) + this.dmf & this.mfunctionMR(:, 2) >= command(2) - this.dmf);
            angle = this.mfunctionMR(angleIndex, 1);
            angle = angle(ceil(length(angle) / 2));
        end

        % calculation of metabolic costs
        function tmpMetCost = getMetCost(this, command)
            cmd = (command * 10) + 1;                                           % scale commands to table entries
            tmpMetCost = interp2(this.metCosts.results, cmd(1), cmd(2), 'spline');   % interpolate in table by cubic splines
        end

        %%% Helper function for calculating {objDist} -> {maxVergErr}
        function vergErrMax = getVergErrMax(this, objDist)
            % correct vergence angle for given object distance
            angleCorrect = 2 * atand(this.baseline / (2 * objDist));
            vergErrMax = angleCorrect - this.angleMin;
        end

        %% Updating images during simulation
        %%% Generates two new images for both eyes
        % simulator:    a renderer instance
        % texture:      file path of texture input
        % eyeAngle:     angle of single eye (rotation from offspring)
        % objDist:      distance of stimulus
        % scaleImSize:  scaling factor of stimulus plane [m]
        function refreshImages(this, simulator, texture, eyeAngle, objDist, scaleImSize)
            simulator.add_texture(1, texture);
            simulator.set_params(1, eyeAngle, objDist, 0, scaleImSize); % scaling of obj plane size

            result1 = simulator.generate_left();
            result2 = simulator.generate_right();

            this.imgRawLeft = permute(reshape(result1, ...
                                         [size(this.imgRawLeft, 3), ...
                                          size(this.imgRawLeft, 2), ...
                                          size(this.imgRawLeft, 1)]), ...
                                         [3, 2, 1]);

            this.imgRawRight = permute(reshape(result2, ...
                                          [size(this.imgRawRight, 3), ...
                                           size(this.imgRawRight, 2), ...
                                           size(this.imgRawRight, 1)]), ...
                                          [3, 2, 1]);

            % convert images to gray scale
            this.imgGrayLeft = 0.2989 * this.imgRawLeft(:, :, 1) + 0.5870 * this.imgRawLeft(:, :, 2) + 0.1140 * this.imgRawLeft(:, :, 3);
            this.imgGrayRight = 0.2989 * this.imgRawRight(:, :, 1) + 0.5870 * this.imgRawRight(:, :, 2) + 0.1140 * this.imgRawRight(:, :, 3);
        end

        %%% Generates two new images for both eyes for experimental renderer
        % simulator:    a renderer instance
        % textureNumber:    index of stimulus in memory buffer
        % eyeAngle:         angle of single eye (rotation from offspring)
        % objDist:          distance of stimulus
        % scaleImSize:  scaling factor of stimulus plane [m]
        function refreshImagesNew(this, simulator, textureNumber, eyeAngle, objDist, scaleImSize)
            simulator.set_params(textureNumber, eyeAngle, objDist, 0, scaleImSize); % scaling of obj plane size

            result1 = simulator.generate_left();
            result2 = simulator.generate_right();

            this.imgRawLeft = permute(reshape(result1, ...
                                         [size(this.imgRawLeft, 3), ...
                                          size(this.imgRawLeft, 2), ...
                                          size(this.imgRawLeft, 1)]), ...
                                         [3, 2, 1]);

            this.imgRawRight = permute(reshape(result2, ...
                                          [size(this.imgRawRight, 3), ...
                                           size(this.imgRawRight, 2), ...
                                           size(this.imgRawRight, 1)]), ...
                                          [3, 2, 1]);

            % convert images to gray scale
            this.imgGrayLeft = 0.2989 * this.imgRawLeft(:, :, 1) + 0.5870 * this.imgRawLeft(:, :, 2) + 0.1140 * this.imgRawLeft(:, :, 3);
            this.imgGrayRight = 0.2989 * this.imgRawRight(:, :, 1) + 0.5870 * this.imgRawRight(:, :, 2) + 0.1140 * this.imgRawRight(:, :, 3);
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

            %% Root Mean Squared Vergence Error
            % windowSize = 125;
            % windowSize = 250;
            windowSize = 1000;
            if (this.trainTime < windowSize * this.interval)
                windowSize = round(this.trainTime / this.interval / 5);
            end

            vergerr = this.vergerr_hist(ind);
            rmse = zeros(length(1 : windowSize : length(vergerr) - mod(length(vergerr), 2)), 1); % cut if odd length
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

            %% Vergence angle / fixation distance
            obsWin = 99; %249; % #last iterations to plot
            figure;
            hold on;
            grid on;
            if (length(this.verge_desired) >= obsWin)
                % plot(this.verge_desired(end - obsWin : end), 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.8);
                % plot(this.verge_actual(end - obsWin : end), 'color', [0, 0.6863, 1.0000], 'LineWidth', 1.3);
                plot(this.Z(this.trainedUntil - obsWin : this.trainedUntil), 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.8);
                plot(this.fixZ(this.trainedUntil - obsWin : this.trainedUntil), 'color', [0, 0.6863, 1.0000], 'LineWidth', 1.3);
            else
                % plot(this.verge_desired, 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.8);
                % plot(this.verge_actual, 'color', [0, 0.6863, 1.0000], 'LineWidth', 1.3);
                plot(this.Z, 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.8);
                plot(this.fixZ, 'color', [0, 0.6863, 1.0000], 'LineWidth', 1.3);
            end

            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            % ylabel('Angle [deg]', 'FontSize', 12);
            ylabel('Object Distance [m]', 'FontSize', 12);
            ylim([this.objDistMin - 1, this.objDistMax + 1]);
            legend('desired', 'actual');
            title(sprintf('Vergence at last %d steps of training', obsWin + 1));
            % plotpath = sprintf('%s/vergenceAngle', this.savePath);
            plotpath = sprintf('%s/fixationDistTraining', this.savePath);
            saveas(gcf, plotpath, 'png');

            %% Muscel graphs
            if (this.rlModel.continuous == 1)
                % windowSize = 1000;
                windowSize = ceil(this.trainTime * 0.002);
                windowSize2 = ceil(this.trainTime * 0.02);
                windowSize3 = ceil(this.trainTime * 0.01);
                if (this.trainTime < windowSize * this.interval)
                    windowSize = round(this.trainTime / this.interval / 5);
                end
                cmd_hist_sma = filter(ones(1, windowSize) / windowSize, 1, this.cmd_hist(ind, 1));
                relCmd_hist_sma = filter(ones(1, windowSize2) / windowSize2, 1, this.relCmd_hist(ind, 1));
                metCost_hist_sma = filter(ones(1, windowSize3) / windowSize3, 1, this.metCost_hist(ind));

                % Two eye muscles
                if (this.rlModel.CActor.output_dim == 2)
                    xVal = [1 : length(ind)];
                    cmd_hist_sma = [cmd_hist_sma, filter(ones(1, windowSize) / windowSize, 1, this.cmd_hist(ind, 2))];
                    relCmd_hist_sma = [relCmd_hist_sma, filter(ones(1, windowSize2) / windowSize2, 1, this.relCmd_hist(ind, 2))];
                    figHandle = figure('OuterPosition', [100, 100, 768, 1024]);

                    % Temporal Rectus
                    % Total muscle commands
                    subplot(3, 2, 1);
                    hold on;
                    grid on;
                    tmpSTD = movingstd(cmd_hist_sma(:, 1), windowSize, 'backward');
                    [hl1, hp1] = boundedline(xVal, ...
                                             cmd_hist_sma(:, 1), ...
                                             tmpSTD, ...
                                             'alpha');

                    hl1.Color = [rand, rand, rand];
                    hp1.FaceColor = hl1.Color;
                    axis([windowSize * 2, length(cmd_hist_sma(:, 1)), ...
                          min(cmd_hist_sma(windowSize * 2 : end, 1) - tmpSTD(windowSize * 2 : end)) * 0.9, ...
                          max(cmd_hist_sma(windowSize * 2 : end, 1) + tmpSTD(windowSize * 2 : end)) * 1.1]);

                    xlabel('Iteration * interval^{-1}', 'FontSize', 8);
                    ylabel(sprintf('Total Muscle\nCommands [%%]'), 'FontSize', 12);
                    title('Lateral Rectus', 'fontweight','normal');

                    % Delta muscle commands
                    subplot(3, 2, 3);
                    hold on;
                    grid on;
                    tmpSTD = movingstd(relCmd_hist_sma(:, 1), windowSize2, 'backward');
                    [hl2, hp2] = boundedline(xVal, ...
                                            relCmd_hist_sma(:, 1), ...
                                            tmpSTD, ...
                                            'alpha');

                    hl2.Color = [rand, rand, rand];
                    hp2.FaceColor = hl2.Color;
                    axis([windowSize2 * 2, length(relCmd_hist_sma(:, 1)), ...
                          min(relCmd_hist_sma(windowSize2 * 2 : end, 1) - tmpSTD(windowSize2 * 2 : end)) * 1.1, ...
                          max(relCmd_hist_sma(windowSize2 * 2 : end, 1) + tmpSTD(windowSize2 * 2 : end)) * 1.1]);

                    xlabel('Iteration * interval^{-1}', 'FontSize', 8);
                    ylabel(strcat('\Delta', sprintf('Muscle\nCommands [%%]')), 'FontSize', 12);

                    % Medial Rectus
                    % Total muscle commands
                    subplot(3, 2, 2);
                    hold on;
                    grid on;
                    tmpSTD = movingstd(cmd_hist_sma(:, 2), windowSize, 'backward');
                    [hl3, hp3] = boundedline(xVal, ...
                                            cmd_hist_sma(:, 2), ...
                                            tmpSTD, ...
                                            'alpha');

                    hl3.Color = [rand, rand, rand];
                    hp3.FaceColor = hl3.Color;
                    axis([windowSize * 2, length(cmd_hist_sma(:, 2)), ...
                          min(cmd_hist_sma(windowSize * 2 : end, 2) - tmpSTD(windowSize * 2 : end)) * 0.9, ...
                          max(cmd_hist_sma(windowSize * 2 : end, 2) + tmpSTD(windowSize * 2 : end)) * 1.1]);

                    xlabel('Iteration * interval^{-1}', 'FontSize', 8);
                    % ylabel('Value', 'FontSize', 12);
                    % set(gca,'yaxislocation','right');
                    title('Medial rectus', 'fontweight','normal');

                    % Delta muscle commands
                    subplot(3, 2, 4);
                    hold on;
                    grid on;
                    tmpSTD = movingstd(relCmd_hist_sma(:, 2), windowSize2, 'backward');
                    [hl4, hp4] = boundedline(xVal, ...
                                            relCmd_hist_sma(:, 2), ...
                                            tmpSTD, ...
                                            'alpha');

                    hl4.Color = [rand, rand, rand];
                    hp4.FaceColor = hl4.Color;
                    axis([windowSize2 * 2, length(relCmd_hist_sma(:, 2)), ...
                          min(relCmd_hist_sma(windowSize2 * 2 : end, 2) - tmpSTD(windowSize2 * 2 : end)) * 1.1, ...
                          max(relCmd_hist_sma(windowSize2 * 2 : end, 2) + tmpSTD(windowSize2 * 2 : end)) * 1.1]);

                    xlabel('Iteration * interval^{-1}', 'FontSize', 8);
                    % ylabel('Value', 'FontSize', 12);
                    % set(gca,'yaxislocation','right');

                    % Metabolic costs
                    subplot(3, 2, 5 : 6);
                    hold on;
                    grid on;
                    tmpSTD = movingstd(metCost_hist_sma, windowSize3, 'backward');
                    [hl5, hp5] = boundedline(xVal, ...
                                            metCost_hist_sma, ...
                                            tmpSTD, ...
                                            'alpha');

                    hl5.Color = [rand, rand, rand];
                    hp5.FaceColor = hl5.Color;
                    axis([windowSize3 * 2, length(metCost_hist_sma), ...
                          min(metCost_hist_sma(windowSize3 * 2 : end) - tmpSTD(windowSize3 * 2 : end)) * 0.9, ...
                          max(metCost_hist_sma(windowSize3 * 2 : end) + tmpSTD(windowSize3 * 2 : end)) * 1.1]);

                    xlabel(sprintf('Iteration * interval^{-1} # (interval=%d)', this.interval), 'FontSize', 8);
                    ylabel('Value', 'FontSize', 12);
                    title('Metabolic Costs', 'FontSize', 14, 'FontWeight','normal');

                    % Subplot overall title
                    suptitle('Muscle Activities');

                    set(figHandle,'PaperPositionMode','auto'); % keep aspect ratio
                    plotpath = sprintf('%s/muscleGraphs', this.savePath);
                    saveas(gcf, plotpath, 'png');
                else
                    % Medial Rectus
                    figure;

                    % Total muscle commands
                    subplot(3, 1, 1);
                    hold on;
                    grid on;
                    tmpSTD = movingstd(cmd_hist_sma(:, 2), windowSize, 'backward');
                    [hl3, hp3] = boundedline(xVal, ...
                                            cmd_hist_sma(:, 2), ...
                                            tmpSTD, ...
                                            'alpha');

                    hl3.Color = [rand, rand, rand];
                    hp3.FaceColor = hl3.Color;
                    axis([windowSize * 2, length(cmd_hist_sma(:, 2)), ...
                          min(cmd_hist_sma(windowSize * 2 : end, 2) - tmpSTD(windowSize * 2 : end)) * 0.9, ...
                          max(cmd_hist_sma(windowSize * 2 : end, 2) + tmpSTD(windowSize * 2 : end)) * 1.1]);

                    xlabel('Iteration * interval^{-1}', 'FontSize', 8);
                    ylabel('Value', 'FontSize', 12);
                    title('Total Muscle Commands', 'fontweight','normal');

                    % Delta muscle commands
                    subplot(3, 1, 2);
                    hold on;
                    grid on;
                    tmpSTD = movingstd(relCmd_hist_sma(:, 2), windowSize2, 'backward');
                    [hl4, hp4] = boundedline(xVal, ...
                                            relCmd_hist_sma(:, 2), ...
                                            tmpSTD, ...
                                            'alpha');

                    hl4.Color = [rand, rand, rand];
                    hp4.FaceColor = hl4.Color;
                    axis([windowSize2 * 2, length(relCmd_hist_sma(:, 2)), ...
                          min(relCmd_hist_sma(windowSize2 * 2 : end, 2) - tmpSTD(windowSize2 * 2 : end)) * 1.1, ...
                          max(relCmd_hist_sma(windowSize2 * 2 : end, 2) + tmpSTD(windowSize2 * 2 : end)) * 1.1]);

                    xlabel('Iteration * interval^{-1}', 'FontSize', 8);
                    ylabel('Value', 'FontSize', 12);
                    title('\Delta Muscle Commands');

                    % Metabolic costs
                    subplot(3, 1, 3);
                    hold on;
                    grid on;
                    tmpSTD = movingstd(metCost_hist_sma, windowSize3, 'backward');
                    [hl5, hp5] = boundedline(xVal, ...
                                            metCost_hist_sma, ...
                                            tmpSTD, ...
                                            'alpha');

                    hl5.Color = [rand, rand, rand];
                    hp5.FaceColor = hl5.Color;
                    axis([windowSize3 * 2, length(metCost_hist_sma), ...
                          min(metCost_hist_sma(windowSize3 * 2 : end) - tmpSTD(windowSize3 * 2 : end)) * 0.9, ...
                          max(metCost_hist_sma(windowSize3 * 2 : end) + tmpSTD(windowSize3 * 2 : end)) * 1.1]);

                    xlabel(sprintf('Iteration * interval^{-1} # (interval=%d)', this.interval), 'FontSize', 8);
                    ylabel('Value', 'FontSize', 12);
                    title('Metabolic Costs', 'FontSize', 14, 'FontWeight','normal');

                    % Subplot overall title
                    suptitle('Muscle Activities (medial rectus)');

                    plotpath = sprintf('%s/muscleGraphsMedialRectus', this.savePath);
                    saveas(gcf, plotpath, 'png');
                end

                if (this.rlModel.CActor.output_dim == 2)
                    % Muscle correlation check
                    % extract last x per cent of iterations
                    pcOffset = 0.01;
                    tmpOffset = ceil(size(this.cmd_hist, 1) * pcOffset);
                    tmpEnd = this.trainedUntil;
                    if tmpOffset > tmpEnd
                        tmpOffset = tmpEnd;
                    end
                    % Total
                    figure;
                    hold on;
                    scatter(this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 1), this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 2), 5,'MarkerFaceColor',[0, 0.7, 0.7]);
                    corrl = corr(this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 1), this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 2));
                    xlabel('Lateral rectus [%]', 'FontSize', 12);
                    ylabel('Medial rectus [%]', 'FontSize', 12);
                    title(strcat('Total Muscle Commands (training)', sprintf('\nCorrelation = %1.2e at last %d iterations', corrl, tmpOffset)));
                    plotpath = sprintf('%s/muscleGraphsScatterTotalTraining', this.savePath);
                    saveas(gcf, plotpath, 'png');

                    % Delta
                    figure;
                    hold on;
                    scatter(this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 1), this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 2), 5,'MarkerFaceColor',[0, 0.7, 0.7]);
                    corrl = corr(this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 1), this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 2));
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
        end

        % Generates anaglyphs of the large and small scale fovea and
        % one of the two unpreprocessed gray scale images
        function generateAnaglyphs(this, identifier, markScales, infos, imgNumber)

            numberScales = length(this.scModel);

            % defining colors in the image: (from larges to smallest scale)
            scalingColors = {'blue', 'red', 'green'};
            textColor = 'yellow';

            scaleImages = cell(numberScales);

            for scale = 1 : numberScales
                % Downsampling Large
                imgLeft = this.imgGrayLeft(:);
                imgLeft = reshape(imgLeft, size(this.imgGrayLeft));
                imgRight = this.imgGrayRight(:);
                imgRight = reshape(imgRight, size(this.imgGrayRight));

                for ds = 1 : log2(model.dsRatio(scale))
                    imgLeft = impyramid(imgLeft, 'reduce');
                    imgRight = impyramid(imgRight, 'reduce');
                end

                % cut fovea in the center
                [h, w, ~] = size(imgLeft);
                cutOutVertical = [fix(h / 2 + 1 - this.pxFieldOfView(scale) / 2), fix(h / 2 + this.pxFieldOfView(scale) / 2)];
                cutOutHorizontal = [fix(w / 2 + 1 - this.pxFieldOfView(scale) / 2), fix(w / 2 + this.pxFieldOfView(scale) / 2)];


                imgLeft = imgLeft(cutOutVertical(1) : cutOutVertical(2) , cutOutHorizontal(1) : cutOutHorizontal(2));
                imgRight = imgRight(cutOutVertical(1) : cutOutVertical(2) , cutOutHorizontal(1) : cutOutHorizontal(2));

                % create an anaglyph of the two pictures, scale it up and save it
                anaglyph = imfuse(imgLeft, imgRight, 'falsecolor');

                % [hAna, vAna, ~] = size(anaglyph);
                if (markScales)
                    anaglyph = insertShape(anaglyph, 'rectangle', [model.pxFieldOfView(scale) + 1 - model.patchSize, 1, model.patchSize, model.patchSize], 'color', scalingColors(scale));
                end

                scaleImages{scale} = anaglyph;
            end

            % imwrite(anaglyph, [savePath '/anaglyph.png']);
            anaglyph = imfuse(this.imgGrayLeft, this.imgGrayRight, 'falsecolor');
            [h, w, ~] = size(anaglyph);

            if (markScales)
                for scale = 1:numberScales
                    % this creates a rectangle inside the image for each scale and
                    % an examplary patch
                    scaleSizeOrigImg = model.pxFieldOfView(scale) * model.dsRatio(scale);
                    patchSizeOrigImg = model.patchSize * model.dsRatio(scale);
                    windowStart = [fix(w / 2 + 1 - scaleSizeOrigImg / 2), fix(h / 2 + 1 - scaleSizeOrigImg / 2)];

                    % indexMatrix = [x, y, scaleWindow, scaleWindow; x, y, patchSize, patchSize]
                    indexMatrix = [windowStart(1), windowStart(2), scaleSizeOrigImg, scaleSizeOrigImg; windowStart(1), windowStart(2), patchSizeOrigImg, patchSizeOrigImg];
                    anaglyph = insertShape(anaglyph, 'rectangle', indexMatrix, 'Color', scalingColors(scale)); % whole region of scale
                end
            end
            % anaglyph = imfuse(leftGray, rightGray, 'falsecolor');

            % imwrite(anaglyph, sprintf('%s/anaglyph%.3d.png', savePath, identifier));

            % imwrite(imresize(anaglyphL, 20), [savePath '/anaglyphLargeScale.png']);
            % imwrite(anaglyphL, sprintf('%s/anaglyphLargeScale%.3d.png', savePath, identifier));
            % largeScaleView = imfuse(imgLeftL, imgRightL, 'montage');
            % imwrite(imresize(largeScaleView, 20), sprintf('%s/LargeScaleMontage%.3d.png', savePath, identifier));

            % imwrite(imresize(anaglyphS, 8), [savePath '/anaglyphSmallScale.png']);

            % imwrite(anaglyphS, sprintf('%s/anaglyphSmallScale%.3d.png', savePath, identifier));
            % smallScaleView = imfuse(imgLeftL, imgRightL, 'montage');
            % imwrite(imresize(smallScaleView, 8), sprintf('%s/smallScaleMontage%.3d.png', savePath, identifier));

            % xRange = [0, 320]; yRange = [0, 240];
            %% infos = {iteration, tmpObjRange, vseRange(vseIndex), angleDes - angleNew, command', relativeCommand', reward, recErrorArray};
            xPos = [10, 220, 10, 10, 260, 10, 10, 240]; %display in headful modus: [10, 200, 10, 10, 260, 10, 10]
            yPos = [10, 10, 230, 220, 230, 190, 200, 220];
            % imName = strsplit(currentTexture, '/');           % stable renderer
            % imName = strsplit(texture{currentTexture}, '/');  % experimental renderer
            imName = strsplit(imgNumber, '/');                  % still necessary / functional?
            imName = imName{end};
            insert = {sprintf('Image: \t%s', imName), ...
                        sprintf('Object distance: \t%.2f', infos{2}), ...
                        sprintf('Vergence Error:          \t%.3f', infos{4}), ...
                        sprintf('Start Vergence Error: \t%.3f', infos{3}), ...
                        sprintf('Iteration: \t%d', infos{1}), ...
                        sprintf('Muscle Activation:    \t%.3f,  \t%.3f', infos{5}(1), infos{5}(2)), ...
                        sprintf('Relative Command: \t%.3f,  \t%.3f', infos{6}(1), infos{6}(2)), ...
                        sprintf('Reward: \t%.3f', infos{7})};

            fig = figure('Position', [0, 0, 960, 640]); % create a figure with a specific size
            % subplot(2,3,[1,2,4,5]);
            % title(sprintf('Object distance: %d,\titeration: %d,\tvergence error: %d', infos(1), infos(2), infos(3)));
            subplot('Position', [0.05, 0.12, 0.6, 0.8]);
            imshow(anaglyph);
            text(xPos, yPos, insert, 'color', textColor); % these numbers resemble white, but white will appear black after saving ;P
            % text(10, 20, num2str(size(anaglyph)), 'color', 'yellow'); %[1-eps, 1-eps, 1-eps]
            % title('Anaglyph');
            % subplot(2,3,3);
            % positionVector = [left, bottom, width, height], all between 0 and 1
            sizeScaleImg = 0.6/numberScales;
            for sp = 1:numberScales
                if (sp == 1)
                    subplot('Position', [0.7, 0.9 - sizeScaleImg, sizeScaleImg, sizeScaleImg])
                else
                    subplot('Position', [0.7, 0.9 - ((sizeScaleImg + 0.05) * sp), sizeScaleImg, sizeScaleImg])
                end
                imshow(scaleImages{sp});
                descr = {sprintf('Scale %d', sp), sprintf('Reconstruction Error: %.3f', infos{8}(sp))};
                % todo: the fist two values in text have to be adapted, especially
                % for more than 2 scales, also the size of the letters,
                text(0.03, 0.1, descr, 'color', textColor, 'Units', 'normalized')
            end
            saveas(fig, sprintf('%s/anaglyph.png', this.savePath), 'png');
            fig.delete();
        end

    end
end

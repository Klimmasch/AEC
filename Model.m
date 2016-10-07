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
        testAt;             % at which iteration steps online testing is performed

        sparseCodingType;   % type of sparse coding

        lambdaMuscleFB;     % factor of muscle activity feedback to RL feature vector
        lambdaRec;          % reconstruction error factor
        lambdaMet;          % factor of metCosts for reward function
        reward_mean;        % tracks the exponentially weighted average of rewards
        reward_variance;    % and the variance for normalizing the reward (if desired)

        % Model data history
        recerr_hist;        % reconstruction error [coarse scale, fine scale]
        % disp_hist;          % disparity
        vergerr_hist;       % vergence error
        % verge_actual;       % actual vergence angle
        % verge_desired;      % desired vergence angle
        % Z;                  % object depth
        % fixZ;               % fixation depth
        td_hist;            % temporal difference (td) error
        feature_hist;       % feature vector
        cmd_hist;           % vergence commands
        relCmd_hist;        % relativ changes in motor commands
        weight_hist;        % L1/L2, i.e. sum abs, sum pow2 weights of actor and critic
        % reward_hist;        % reward function
        metCost_hist;       % metabolic costs
        variance_hist;      % exploratory variance of actor
        savePath;           % where all the data are stored
        notes;              % is there any special things about this model to note?

        % Model results at testing procedure
        responseResults;
        testResult;
        testResult2;
        testResult3;
        testResult4;
        testResult5;
        testResult6;
        testResult7;
        testHist;           % history of testing performance over traintime

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
        metCostsIncRes; % contains a part of metCosts with increased resolution
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
        scaleFacMR;
        scaleFacLR;

        imgRawLeft;     % images that are updated by refreshImages
        imgRawRight;
        imgGrayLeft;
        imgGrayRight;
    end

    methods
        function obj = Model(PARAM)
            obj.textureFile = PARAM{1}{1};
            obj.trainTime = PARAM{1}{2};
            obj.testAt = PARAM{1}{3};
            obj.sparseCodingType = PARAM{1}{4};
            obj.focalLength = PARAM{1}{5};
            obj.baseline = PARAM{1}{6};
            obj.objDistMin = PARAM{1}{7};
            obj.objDistMax = PARAM{1}{8};
            obj.muscleInitMin = PARAM{1}{9};
            obj.muscleInitMax = PARAM{1}{10};
            obj.interval = PARAM{1}{11};
            obj.lambdaMuscleFB = PARAM{1}{12};
            obj.lambdaRec = PARAM{1}{13};
            obj.lambdaMet = PARAM{1}{14};
            obj.fixDistMin = PARAM{1}{19};
            obj.fixDistMax = PARAM{1}{20};

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

            % obj.recerr_hist = zeros(obj.trainTime, length(PARAM{2}{1})); % recerr_hist = t x #SC_scales
            obj.recerr_hist = zeros(obj.trainTime / obj.interval, length(PARAM{2}{1})); % recerr_hist = t x #SC_scales
            % obj.disp_hist = zeros(obj.trainTime, 1);
            obj.vergerr_hist = zeros(obj.trainTime, 1);
            % obj.verge_actual = zeros(obj.trainTime, 1);
            % obj.verge_desired = zeros(obj.trainTime, 1);
            % obj.Z = zeros(obj.trainTime, 1);
            % obj.fixZ = zeros(obj.trainTime, 1);

            obj.td_hist = zeros(obj.trainTime, 1);
            % obj.feature_hist = zeros(obj.trainTime, 1);%zeros(obj.trainTime, PARAM{3}{9}(1));
            obj.cmd_hist = zeros(obj.trainTime, 2);
            obj.relCmd_hist = zeros(obj.trainTime, PARAM{3}{9}(3)); % relCmd_hist = t x output_dim
            % obj.weight_hist = zeros(obj.trainTime, 4);
            obj.weight_hist = zeros(obj.trainTime, 6); % for also traking change in weights
            % obj.reward_hist = zeros(obj.trainTime, 1);
            obj.metCost_hist = zeros(obj.trainTime, 1);
            obj.variance_hist = zeros(obj.trainTime, 1);

            obj.responseResults = struct();
            obj.testResult = [];
            obj.testResult2 = [];
            obj.testResult3 = [];
            obj.testResult4 = [];
            obj.testResult5 = [];
            obj.testResult6 = [];
            obj.testResult7 = [];
            % rmse(vergerr), mean(abs(vergErr)), std(abs(vergErr)), rmse(deltaMetCost), mean(abs(deltaMetCost)), std(abs(deltaMetCost))
            obj.testHist = zeros(length(obj.testAt), 6);
            obj.testHist(1, :) = [1.1593, -0.074, 7.4197, 1.0736, 0.9747, 3.3798]; % insert modelAt0 entry
            obj.simulatedTime = 0;
            obj.trainedUntil = 0;
            obj.notes = '';

            obj.reward_mean = 0;
            obj.reward_variance = 1;

            %%% Generate image processing constants
            obj.patchSize = PARAM{1}{15};
            obj.pxFieldOfView = PARAM{1}{16};
            obj.dsRatio = PARAM{1}{17};
            obj.stride = PARAM{1}{18};
            obj.overlap = PARAM{1}{21};
            obj.cutout = PARAM{1}{22};

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

            % between 0 and 0.2/0.1 mus. act. the resulution is increased
            usedRows = 3;
            usedCols = 2;
            obj.degreesIncRes = interp2(obj.degrees.results_deg(1 : usedRows, 1 : usedCols), resFactor);
            obj.degDiff = max(max(diff(obj.degreesIncRes)));                                                % distance between the entries
            obj.scaleFacMR = ((usedRows - 1) / 10) / size(obj.degreesIncRes, 1);                            % table scaling factors for backwards
            obj.scaleFacLR = ((usedCols - 1) / 10) / size(obj.degreesIncRes, 2);                            % calculation of table_index -> muscle inervation

            % increased resolution of metCosts table
            obj.metCostsIncRes = interp2(obj.metCosts.results(1 : usedRows, 1 : usedCols), resFactor);

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
        % scScale:  SC scale index elem {coarse := 1, 2, ..., fine = #}
        % eyePos:   eye position index elem {1 := left, 2 := right}
        function preprocessImage(this, scScale, eyePos)
            if (eyePos == 1)
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
                    sprintf('All values in the extracted patches are zero. Check if the image rendering is alright!')
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

        % Calculates muscle force for medial rectus
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

        % Calculates muscle force for two muscles
        % the activation of on muscle is always 0.
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

        % Calculates muscle force for two muscles
        % drawn from all permitted mf(l, m) ^= f(objDist, desVergErr), where l, m >= 0
        % == get muscle force equally distributed over object distance.
        function [command, angleInit] = getMFedood(this, objDist, desVergErr)
            angleCorrect = 2 * atand(this.baseline / (2 * objDist));
            angleInit = angleCorrect - desVergErr;

            % angleInit is the angle for both eyes, but degreesIncRes only
            % contains angles for one eye
            [xi, yi] = find(this.degreesIncRes <= (angleInit / 2) + this.degDiff & this.degreesIncRes >= (angleInit / 2) - this.degDiff);

            i = randi(length(xi));

            % transform indizes to muscle activities
            mfMR = xi(i) * this.scaleFacMR;
            mfLR = yi(i) * this.scaleFacLR;

            command = [mfLR; mfMR];
        end

        % Maps {objDist, desVergErr} -> {medialRectusActivations, lateralRectusActivations},
        % i.e. calculates all muscle activity cominations corresponding to specified {objDist, desVergErr}
        function [mfLR, mfMR] = getAnglePoints(this, objDist, desVergErr)
            angleCorrect = 2 * atand(this.baseline / (2 * objDist));
            angleInit = angleCorrect - desVergErr;

            % angleInit is the angle for both eyes, but degreesIncRes only
            % contains angles for one eye
            [xi, yi] = find(this.degreesIncRes <= (angleInit / 2) + this.degDiff & this.degreesIncRes >= (angleInit / 2) - this.degDiff);

            % transform indizes to muscle activities
            mfMR = xi .* this.scaleFacMR;
            mfLR = yi .* this.scaleFacLR;
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

        %% Plot all gathered performance data and save graphs
        %  @param level:    # of plot elem range [1, 7]
        function allPlotSave(this, level)

            % only take the last value before the image/texture is changed
            ind = this.interval : this.interval : this.trainTime;
            windowSize3 = ceil(this.trainTime * 0.01);
            metCost_hist_sma = filter(ones(1, windowSize3) / windowSize3, 1, this.metCost_hist(ind));

            windowSize = 1000;

            %% Vergence error
            if (~isempty(find(level == 1)))
                % windowSize = 125;
                % windowSize = 1000;
                if (this.trainTime < windowSize * this.interval)
                    windowSize = round(this.trainTime / this.interval / 5);
                end

                %% Simple Moving Average Vergence Error
                vergerr = filter(ones(1, windowSize) / windowSize, 1, abs(this.vergerr_hist(ind)));

                figure;
                hold on;
                grid on;
                plot((windowSize + 1) * this.interval : this.interval : size(this.vergerr_hist), vergerr(windowSize + 1 : end), ...
                     'LineWidth', 1.3);

                xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
                ylabel('SMA(|verg_{err}|) [deg]', 'FontSize', 12);
                title('Moving Average of Vergence Error');
                plotpath = sprintf('%s/vergErrSMA', this.savePath);
                saveas(gcf, plotpath, 'png');

                %% Root Mean Squared Vergence Error
                % windowSize = 125;
                % windowSize = 250;
                % windowSize = 1000;
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
                    plotpath = sprintf('%s/vergErrRMSE', this.savePath);
                    saveas(gcf, plotpath, 'png');
                catch
                    % if windowsize > values in vergerr
                    sprintf('Warning: windowSize >= vergerr')
                end
            end

            %% Reconstruction Error
            if (~isempty(find(level == 2)))
                try
                    figure;
                    hold on;
                    grid on;
                    handleArray = zeros(1, length(this.scModel));
                    for i = 1 : length(this.scModel)
                        % tmpError = filter(ones(1, windowSize) / windowSize, 1, this.recerr_hist(ind, i));
                        % handleArray(i) = plot((windowSize + 1) * this.interval : this.interval : size(this.recerr_hist, 1), ...
                        %                       tmpError(windowSize + 1 : end), ...
                        %                       'color', [rand, rand, rand], 'LineWidth', 1.3);
                        tmpError = filter(ones(1, windowSize) / windowSize, 1, this.recerr_hist(:, i));
                        handleArray(i) = plot((windowSize + 1) * this.interval : size(this.recerr_hist, 1), ...
                                              tmpError(windowSize + 1 : end), ...
                                              'color', [rand, rand, rand], 'LineWidth', 1.3);
                    end
                    xlabel(sprintf('Episode # (interval=%d)', this.interval), 'FontSize', 12);
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
            end

            % %% Vergence angle / fixation distance
            % if (~isempty(find(level == 3)))
            %     obsWin = 99; %249; % #last iterations to plot
            %     figure;
            %     hold on;
            %     grid on;
            %     if (length(this.verge_desired) >= obsWin)
            %         % plot(this.verge_desired(end - obsWin : end), 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.8);
            %         % plot(this.verge_actual(end - obsWin : end), 'color', [0, 0.6863, 1.0000], 'LineWidth', 1.3);
            %         plot(this.Z(this.trainedUntil - obsWin : this.trainedUntil), 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.8);
            %         plot(this.fixZ(this.trainedUntil - obsWin : this.trainedUntil), 'color', [0, 0.6863, 1.0000], 'LineWidth', 1.3);
            %     else
            %         % plot(this.verge_desired, 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.8);
            %         % plot(this.verge_actual, 'color', [0, 0.6863, 1.0000], 'LineWidth', 1.3);
            %         plot(this.Z, 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.8);
            %         plot(this.fixZ, 'color', [0, 0.6863, 1.0000], 'LineWidth', 1.3);
            %     end

            %     xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            %     % ylabel('Angle [deg]', 'FontSize', 12);
            %     ylabel('Object Distance [m]', 'FontSize', 12);
            %     ylim([this.objDistMin - 1, this.objDistMax + 1]);
            %     legend('desired', 'actual');
            %     title(sprintf('Vergence at last %d steps of training', obsWin + 1));
            %     % plotpath = sprintf('%s/vergenceAngle', this.savePath);
            %     plotpath = sprintf('%s/fixationDistTraining', this.savePath);
            %     saveas(gcf, plotpath, 'png');
            % end

            %% Muscel graphs
            if (~isempty(find(level == 4)))
                if (this.rlModel.continuous == 1)
                    ind2 = 1 : 25 : this.trainTime;
                    % windowSize = 1000;
                    windowSize = ceil(this.trainTime * 0.005);
                    windowSize2 = ceil(this.trainTime * 0.01);
                    if (this.trainTime < windowSize * this.interval)
                        windowSize = round(this.trainTime / this.interval / 5);
                    end
                    cmd_hist_sma = filter(ones(1, windowSize) / windowSize, 1, this.cmd_hist(ind2, 1));
                    relCmd_hist_sma = filter(ones(1, windowSize2) / windowSize2, 1, this.relCmd_hist(ind2, 1));

                    % Two eye muscles
                    if (this.rlModel.CActor.output_dim == 2)
                        xVal = [1 : length(cmd_hist_sma)];
                        cmd_hist_sma = [cmd_hist_sma, filter(ones(1, windowSize) / windowSize, 1, this.cmd_hist(ind2, 2))];
                        relCmd_hist_sma = [relCmd_hist_sma, filter(ones(1, windowSize2) / windowSize2, 1, this.relCmd_hist(ind2, 2))];
                        figHandle = figure('OuterPosition', [100, 100, 768, 1024]); % static figure resolution/size

                        % Lateral Rectus
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

                        xlabel('Iteration #', 'FontSize', 8);
                        ylabel(sprintf('Total Muscle\nCommands [%%]'), 'FontSize', 12);
                        title('Lateral Rectus', 'fontweight','normal');

                        % Delta muscle commands
                        xVal = [1 : length(relCmd_hist_sma)];
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

                        xlabel('Iteration #', 'FontSize', 8);
                        ylabel(strcat('\Delta', sprintf('Muscle\nCommands [%%]')), 'FontSize', 12);

                        % Medial Rectus
                        % Total muscle commands
                        xVal = [1 : length(cmd_hist_sma)];
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

                        xlabel('Iteration #', 'FontSize', 8);
                        % ylabel('Value', 'FontSize', 12);
                        % set(gca,'yaxislocation','right');
                        title('Medial rectus', 'fontweight','normal');

                        % Delta muscle commands
                        xVal = [1 : length(relCmd_hist_sma)];
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

                        xlabel('Iteration #', 'FontSize', 8);
                        % ylabel('Value', 'FontSize', 12);
                        % set(gca,'yaxislocation','right');

                        % Metabolic costs
                        xVal = [1 : length(metCost_hist_sma)];
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
                              min(metCost_hist_sma(windowSize3 * 2 : end) - tmpSTD(windowSize3 * 2 : end)) * 0.95, ...
                              max(metCost_hist_sma(windowSize3 * 2 : end) + tmpSTD(windowSize3 * 2 : end)) * 1.05]);

                        xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 8);
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
                        xVal = [1 : length(cmd_hist_sma)];
                        subplot(3, 1, 1);
                        hold on;
                        grid on;
                        tmpSTD = movingstd(cmd_hist_sma, windowSize, 'backward');
                        [hl3, hp3] = boundedline(xVal, ...
                                                cmd_hist_sma, ...
                                                tmpSTD, ...
                                                'alpha');

                        hl3.Color = [rand, rand, rand];
                        hp3.FaceColor = hl3.Color;
                        axis([windowSize * 2, length(cmd_hist_sma), ...
                              min(cmd_hist_sma(windowSize * 2 : end) - tmpSTD(windowSize * 2 : end)) * 0.9, ...
                              max(cmd_hist_sma(windowSize * 2 : end) + tmpSTD(windowSize * 2 : end)) * 1.1]);

                        xlabel('Iteration * interval^{-1}', 'FontSize', 8);
                        ylabel('Value', 'FontSize', 12);
                        title('Total Muscle Commands', 'fontweight','normal');

                        % Delta muscle commands
                        xVal = [1 : length(relCmd_hist_sma)];
                        subplot(3, 1, 2);
                        hold on;
                        grid on;
                        tmpSTD = movingstd(relCmd_hist_sma, windowSize2, 'backward');
                        [hl4, hp4] = boundedline(xVal, ...
                                                relCmd_hist_sma, ...
                                                tmpSTD, ...
                                                'alpha');

                        hl4.Color = [rand, rand, rand];
                        hp4.FaceColor = hl4.Color;
                        axis([windowSize2 * 2, length(relCmd_hist_sma), ...
                              min(relCmd_hist_sma(windowSize2 * 2 : end) - tmpSTD(windowSize2 * 2 : end)) * 1.1, ...
                              max(relCmd_hist_sma(windowSize2 * 2 : end) + tmpSTD(windowSize2 * 2 : end)) * 1.1]);

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
                        % scatter(this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 1), this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 2), 5,'MarkerFaceColor',[0, 0.7, 0.7]);
                        histHandle = hist3(this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, :), [40, 40]);

                        corrl = corr(this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 1), this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 2));
                        xb = linspace(min(this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 1)), max(this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 1)), size(histHandle, 1));
                        yb = linspace(min(this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 2)), max(this.cmd_hist(tmpEnd - tmpOffset : tmpEnd, 2)), size(histHandle, 1));

                        xlabel('Lateral rectus [%]', 'FontSize', 12);
                        ylabel('Medial rectus [%]', 'FontSize', 12);
                        title(strcat('Total Muscle Commands (training)', sprintf('\nCorrelation = %1.2e at last %d iterations', corrl, tmpOffset)));

                        pcHandle = pcolor(xb, yb, histHandle);
                        axis([0, xb(end), 0, yb(end)]);
                        shading interp;
                        set(pcHandle, 'EdgeColor', 'none');

                        colormap(createCM(1));
                        cb = colorbar();
                        cb.Label.String = '# Occurences';

                        plotpath = sprintf('%s/muscleTotalCmdTraining', this.savePath);
                        saveas(gcf, plotpath, 'png');

                        % Delta
                        figure;
                        hold on;
                        % scatter(this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 1), this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 2), 5,'MarkerFaceColor',[0, 0.7, 0.7]);
                        histHandle = hist3(this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, :), [40, 40]);

                        corrl = corr(this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 1), this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 2));
                        xb = linspace(min(this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 1)), max(this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 1)), size(histHandle, 1));
                        yb = linspace(min(this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 2)), max(this.relCmd_hist(tmpEnd - tmpOffset : tmpEnd, 2)), size(histHandle, 1));

                        xlabel('Lateral rectus [%]', 'FontSize', 12);
                        ylabel('Medial rectus [%]', 'FontSize', 12);
                        title(strcat('\Delta Muscle Commands (training)', sprintf('\nCorrelation = %1.2e at last %d iterations', corrl, tmpOffset)));

                        pcHandle = pcolor(xb, yb, histHandle);
                        axis([xb(1), xb(end), yb(1), yb(end)]);
                        shading interp;
                        set(pcHandle, 'EdgeColor', 'none');

                        colormap(createCM(1));
                        cb = colorbar();
                        cb.Label.String = '# Occurences';

                        plotpath = sprintf('%s/muscleDeltaCmdTraining', this.savePath);
                        saveas(gcf, plotpath, 'png');
                    end
                end
            end

            %% Weights
            if (~isempty(find(level == 5)))
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
            end

            %% Reward composition
            if (~isempty(find(level == 6)))
                if (this.rlModel.continuous == 1)
                    windowSize = ceil(this.trainTime * 0.01);

                    % backward compatibility
                    if (length(this.recerr_hist) == length(this.trainTime))
                        tmp = this.recerr_hist(ind, :);
                        this.recerr_hist = zeros(this.trainTime / this.interval, size(this.recerr_hist, 2));
                        this.recerr_hist = tmp;
                    end

                    recerr_hist_sma = filter(ones(1, windowSize) / windowSize, 1, sum(this.recerr_hist, 2));

                    figure;
                    hold on;
                    grid on;
                    r = [- this.lambdaMet * metCost_hist_sma, ...
                         - this.lambdaRec * recerr_hist_sma];
                    handle = area(r, 'LineStyle','none');
                    xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
                    ylabel('Value', 'FontSize', 12);
                    l = legend('\lambdaMetabolic_{cost}', 'Reconstruction_{error}');
                    handle(1).FaceColor = [1, 0.25, 0];
                    handle(2).FaceColor = [1, 0.549, 0];
                    axis([windowSize, inf, -inf, 0]);
                    % l.Location = 'southwest';
                    l.Location = 'best';
                    tmpstr = '\lambda = ';
                    title(strcat(sprintf('Reward composition (SMA)\n%s %6.4f, ~%4.1f%% Metabolic Costs',tmpstr, this.lambdaMet, (this.lambdaMet / 0.1622) * 100)));
                    plotpath = sprintf('%s/rewardComp', this.savePath);
                    saveas(gcf, plotpath, 'png');
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

            %% Testing performance as a function of traintime
            if (~isempty(find(level == 7)))
                if (isfield(this, 'testHist'))
                    % RMSE vergErr
                    figure;
                    subplot(2, 2, 1);
                    hold on;
                    grid on;
                    plot(this.testAt, this.testHist(:, 1), 'x-', 'LineWidth', 1.3);

                    xlabel('Traintime', 'FontSize', 12);
                    ylabel('RMSE(verg_{err}) [deg]', 'FontSize', 12);

                    % mean, std vergErr
                    subplot(2, 2, 2);
                    hold on;
                    grid on;
                    [hl, hp] = boundedline(this.testAt, this.testHist(:, 2), this.testHist(:, 3), 'alpha');

                    hl.Marker = 'x';
                    hl.MarkerSize = 4;

                    hl.Color = [rand, rand, rand];
                    hp.FaceColor = hl.Color;
                    hl.LineWidth = 1.6;

                    xlabel('Traintime', 'FontSize', 12);
                    ylabel('|verg_{err}| [deg]', 'FontSize', 12);

                    % RMSE deltaMetCost
                    subplot(2, 2, 3);
                    hold on;
                    grid on;
                    plot(this.testAt, this.testHist(:, 4), 'x-', 'LineWidth', 1.3);

                    xlabel('Traintime', 'FontSize', 12);
                    ylabel('RMSE(|\Deltamc|)', 'FontSize', 12);

                    % mean, std deltaMetCost
                    subplot(2, 2, 4);
                    hold on;
                    grid on;
                    [hl, hp] = boundedline(this.testAt, this.testHist(:, 5), this.testHist(:, 6), 'alpha');

                    hl.Marker = 'x';
                    hl.MarkerSize = 4;

                    hl.Color = [rand, rand, rand];
                    hp.FaceColor = hl.Color;
                    hl.LineWidth = 1.6;

                    xlabel('Traintime', 'FontSize', 12);
                    ylabel('|\Deltamc| = |mc_{actual} - mc_{desired}|', 'FontSize', 12);

                    % Subplot overall title
                    suptitle('Performance vs. Traintime');

                    plotpath = sprintf('%s/performanceVsTraintime', this.savePath);
                    saveas(gcf, plotpath, 'png');
                else
                    % TODO: add backward compatibility support
                    % sprintf('Model is too old for allPlotSave(7), backward compatibility is not supported yet!')
                    % return;

                    % try to update the model
                    try
                        % generate new model attributes
                        clone = this.copy();
                        tmpStr = GetFullPath(clone.savePath)
                        tmpIndex = strfind(tmpStr, 'model');
                        if (length(tmpIndex) == 1)
                            clone.savePath = strcat('/home/aecgroup/aecdata/Results/', tmpStr(tmpIndex(1) : end));
                        else
                            clone.savePath = strcat('/home/aecgroup/aecdata/Results/', tmpStr(tmpIndex(1) : tmpIndex(2) - 2));
                        end

                        % Get a list of all files and folders in this folder
                        files = dir(clone.savePath);
                        % Get a logical vector that tells which is a directory
                        dirFlags = [files.isdir];
                        % Extract only those that are directories.
                        subFolders = files(dirFlags);

                        % Extract all relevant data
                        for k = 1 : length(subFolders)
                            % load intermediate model objects
                            tmpModel = load(strcat(clone.savePath, '/', subFolders(k).name));

                            % get test iteration number
                            clone.testAt(k + 1) = str2num(subFolders(k).name(9 : end));

                            % generate test history
                            % assumes testInterval = 20
                            clone.testHist(k + 1, :) = [sqrt(mean(tmpModel.testResult3(:, 20) .^ 2)), ...
                                                        mean(abs(tmpModel.testResult3(:, 20) .^ 2)), ...
                                                        std(abs(tmpModel.testResult3(:, 20) .^ 2)), ...
                                                        sqrt(mean(tmpModel.testResult7(:, 20) .^ 2)), ...
                                                        mean(abs(tmpModel.testResult7(:, 20) .^ 2)), ...
                                                        std(abs(tmpModel.testResult7(:, 20) .^ 2))];
                        end

                        % delete(model);
                        % clear model;
                        % model = clone;

                        % overwrite old model
                        save(strcat(clone.savePath, '/model'), 'model');

                        % delete(clone);
                        % clear clone;

                        % try to plot again
                        clone.allPlotSave(7);
                    catch
                        sprintf('Error: One or more file operations failed. Check path strings.')
                        return;
                    end
                end
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

                for ds = 1 : log2(this.dsRatio(scale))
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
                    anaglyph = insertShape(anaglyph, 'rectangle', [this.pxFieldOfView(scale) + 1 - this.patchSize, 1, this.patchSize, this.patchSize], 'color', scalingColors(scale));
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
                    scaleSizeOrigImg = this.pxFieldOfView(scale) * this.dsRatio(scale);
                    patchSizeOrigImg = this.patchSize * this.dsRatio(scale);
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

        %%% This method creates a trajectory from the given paramters and plots it
        %%% on the plane of object depth.
        %%% Note that this script is intended solely for continuous models.
        % @param objDist                the object distance
        % @param startVergErr           the vergence error to muscles start with
        % @param initMethod             either 'simple' or 'random'
        % @param numIters               number of iterations that are executed
        % @param stimuliIndices         an array of indizes from the texture files
        % @param simulator              either a simulator object or [] for a new one
        % @param titleStr               string identifier that is used for the title and the saved image
        % @param savePlot               true or false if the resulting plots should be saved
        %%%
        %%TODO:
        % enable multiple fixation dists in one plot with same init values
        % idea: instead of contourf, just plot single lines that correspond to
        % spec. obj. dists
        function h = plotTrajectory(this, objDist, startVergErr, initMethod, numIters, stimuliIndices, simulator, directory, titleStr, savePlot)
            % simulator check
            if (isempty(simulator))
                sprintf('An initialized simulator is necessary to continue.\nPlease execute simulator = prepareSimulator();')
                return;
            end

            %%% Saturation function that keeps motor commands in [0, 1]
            %   corresponding to the muscelActivity/metabolicCost tables
            function [cmd] = checkCmd(cmd)
                i0 = cmd < 0;
                cmd(i0) = 0;
                i1 = cmd > 1;
                cmd(i1) = 1;
            end

            % preperation
            rng(22);

            if strcmp(initMethod, 'advanced')
                initMethod = uint8(0);

            elseif strcmp(initMethod, 'fixed')
                initMethod = uint8(1);
                % hand-picked inits for muscles, used in initMethod 'random'
                cmdInit = [[0.03; 0.16], [0.05; 0.12], [0.07; 0.12], [0.04; 0.08], [0.06; 0.06], [0.08; 0.06]];

            elseif strcmp(initMethod, 'simple')
                initMethod = uint8(2);

            else
                sprintf('Muscle initialization method %s not supported.', initMethod)
                return;
            end

            plotAnaglyphs = false;
            nStimuli = length(stimuliIndices);
            trajectory = zeros(length(objDist), length(startVergErr), nStimuli, numIters + 1, 2);

            %% main loop
            if (plotAnaglyphs == true)
                figure;
                figIter = 1;
            end

            for odIndex = 1 : length(objDist)
                angleDes = 2 * atand(this.baseline / (2 * objDist(odIndex)));

                for stimIter = 1 : nStimuli
                    % currentTexture = stimuliIndices(stimIter);
                    currentTexture = randi(40); % hack

                    for vergErrIndex = 1 : length(startVergErr)
                        % muscle init
                        if (initMethod == 0)
                            try
                                % catch negative/divergent vergence angle
                                [command, angleNew] = this.getMFedood(objDist(odIndex), min(startVergErr(vergErrIndex), this.getVergErrMax(objDist(odIndex))));
                            catch
                                % catch non-existing variables error, occuring in non-up-to-date models
                                try
                                    clone = this.copy();
                                    delete(this);
                                    clear this;
                                    this = clone;
                                    [command, angleNew] = this.getMFedood(objDist(odIndex), startVergErr(vergErrIndex));
                                    delete(clone);
                                    clear clone;
                                catch
                                    % catch when new model property isn't present in Model class yet
                                    sprintf('Error: One or more new model properties (variables) are not present in Model.m class yet!')
                                    return;
                                end
                            end
                        elseif (initMethod == 1)
                            command = cmdInit(:, stimIter);
                            angleNew = this.getAngle(command);
                        elseif (initMethod == 2)
                            [command, angleNew] = this.getMF2(objDist(odIndex), startVergErr(vergErrIndex));
                        end
                        trajectory(odIndex, vergErrIndex, stimIter, 1, :) = command;

                        for iter = 1 : numIters
                            this.refreshImagesNew(simulator, currentTexture, angleNew / 2, objDist(odIndex), 3);

                            % show anaglyphs for quit performance check
                            if (plotAnaglyphs && ((iter == 1) || (iter == numIters)))
                                subplot(length(objDist) * length(startVergErr) * nStimuli, 2, figIter);
                                % imshow(stereoAnaglyph(this.imgGrayLeft, this.imgGrayRight))
                                imshow(imfuse(this.imgGrayLeft, this.imgGrayRight, 'falsecolor'))
                                if (iter == 1)
                                    title(sprintf('fix. depth = %1.1fm (%.3f°)\nverg. error = %.3f', ...
                                          (this.baseline / 2) / tand(angleNew / 2), angleNew, angleDes - angleNew));
                                end
                                figIter = figIter + 1;
                            end

                            for i = 1 : length(this.scModel)
                                this.preprocessImage(i, 1);
                                this.preprocessImage(i, 2);
                                currentView{i} = vertcat(this.patchesLeft{i}, this.patchesRight{i});
                            end

                            [bfFeature, ~, ~] = this.generateFR(currentView);              % encode image patches
                            feature = [bfFeature; command * this.lambdaMuscleFB];          % append muscle activities to feature vector
                            relativeCommand = this.rlModel.act(feature);                   % generate change in muscle activity
                            command = checkCmd(command + relativeCommand);                 % calculate new muscle activities
                            angleNew = this.getAngle(command) * 2;                         % transform into angle

                            trajectory(odIndex, vergErrIndex, stimIter, iter + 1, :) = command;

                            if (plotAnaglyphs && (iter == numIters))
                                title(sprintf('fix. depth = %1.1fm (%.3f°)\nverg. error = %.3f', ...
                                              (this.baseline / 2) / tand(angleNew / 2), angleNew, angleDes - angleNew));
                            end
                        end
                    end
                end
            end

            %% Plot results
            h = figure();
            hold on;
            if (isempty(titleStr))
                title('Object Fixation Trajectories');
            else
                title(sprintf('Object Fixation Trajectories\nmodel trained for #%s iterations', titleStr));
            end

            % pcHandle = pcolor(this.degreesIncRes); % use vergence degree as color dimension (background)
            pcHandle = pcolor(this.metCostsIncRes);  % use metabolic costs as color dimension (background)
            % shading interp;
            set(pcHandle, 'EdgeColor', 'none');

            % colormap(createCM(1));
            cb = colorbar();
            % cb.Label.String = 'vergence degree'; % use vergence degree as color dimension (background)
            cb.Label.String = 'metabolic costs';   % use metabolic costs as color dimension (background)

            ax = gca;
            set(ax, 'Layer','top'); % bring axis to the front

            ax.XTick = linspace(1, size(this.degreesIncRes, 2), 11);
            ax.YTick = linspace(1, size(this.degreesIncRes, 1), 11);

            ax.XTickLabel = strsplit(num2str(linspace(0, 10, 11)));
            ax.YTickLabel = strsplit(num2str(linspace(0, 20, 11)));

            axis([1, size(this.degreesIncRes, 2), 1, size(this.degreesIncRes, 1)]);

            for odIndex = 1 : length(objDist)
                % draw +1 pixel offset in respect to desired vergence distance
                [lateralDes, medialDes] = this.getAnglePoints(objDist(odIndex), 0.22);
                plot(lateralDes ./ this.scaleFacLR, medialDes ./ this.scaleFacMR, ...
                     'color', [0, 0.5882, 0], 'LineStyle', ':', 'LineWidth', 1.8);

                % draw -1 pixel offset in respect to desired vergence distance
                [lateralDes, medialDes] = this.getAnglePoints(objDist(odIndex), -0.22);
                plot(lateralDes ./ this.scaleFacLR, medialDes ./ this.scaleFacMR, ...
                     'color', [0, 0.5882, 0], 'LineStyle', ':', 'LineWidth', 1.8);

                % draw a line of points into the plane that represent the desired vergence
                [lateralDes, medialDes] = this.getAnglePoints(objDist(odIndex), 0);
                plot(lateralDes ./ this.scaleFacLR, medialDes ./ this.scaleFacMR, ...
                     'color', [0.6510, 1.0000, 0.6588], 'LineWidth', 1.8);

                % add corresponding distance value to desired vergence graph
                text(lateralDes(end - ceil(length(lateralDes) / 10)) / this.scaleFacLR, ...
                     medialDes(end - ceil(length(medialDes) / 10)) / this.scaleFacMR, ...
                     sprintf('%3.1fm', objDist(odIndex)));
            end

            % draw trajectories
            for odIndex = 1 : length(objDist)
                % tmpRed = [rand, 0, 0];
                % lightOrange = [1, 201 / 255, 41 / 255]
                % orange = [1, 94 / 255, 41 / 255],
                for stim = 1 : length(stimuliIndices)
                    for vergErrIndex = 1 : length(startVergErr)

                        % first plot all
                        plot(reshape(trajectory(odIndex, vergErrIndex, stim, :, 1), [numIters + 1, 1]) ./ this.scaleFacLR + 1, ...
                             reshape(trajectory(odIndex, vergErrIndex, stim, :, 2), [numIters + 1, 1]) ./ this.scaleFacMR + 1, ...
                             '.-', 'color', [1, 201 / 255, 41 / 255],  'LineWidth', 1.5, 'MarkerSize', 7);

                        % plot iter 1-interval in differen color if numIters >= model.interval
                        if (numIters >= this.interval)
                            plot(reshape(trajectory(odIndex, vergErrIndex, stim, 1 : this.interval, 1), [this.interval, 1]) ./ this.scaleFacLR + 1, ...
                                 reshape(trajectory(odIndex, vergErrIndex, stim, 1 : this.interval, 2), [this.interval, 1]) ./ this.scaleFacMR + 1, ...
                                 '.-', 'color', 'r', 'LineWidth', 1.5, 'MarkerSize', 7); % vielleicht ne andere Farbe ...
                        end

                        % plot init point
                        plot(trajectory(odIndex, vergErrIndex, stim, 1, 1) / this.scaleFacLR + 1, ...
                             trajectory(odIndex, vergErrIndex, stim, 1, 2) / this.scaleFacMR + 1, ...
                            'r.', 'MarkerSize', 15);
                            % 'color', tmpRed, 'MarkerSize', 15); 'color', [41 / 255, 1, 94 / 255],

                        % plot destination point
                        plot(trajectory(odIndex, vergErrIndex, stim, end, 1) / this.scaleFacLR + 1, ...
                             trajectory(odIndex, vergErrIndex, stim, end, 2) / this.scaleFacMR + 1, ...
                             'g.',  'MarkerSize', 15);
                    end
                end
            end

            xlabel('lateral rectus activation [%]');
            ylabel('medial rectus activation [%]');

            if (savePlot == 1)
                if (~isempty(directory))
                    plotpath = sprintf('%s/muscleActivityTrajectory.png', directory);
                else
                    plotpath = 'muscleActivityTrajectory.png';
                end
                saveas(h, plotpath, 'png');
            end
        end

        %%% This methods displays the current binocular basis functions of the model.
        %%% If the history of basis functions was saved, their development is displayed.
        function displayBasis(this, savePlot)
            % r = 16; c = 18; %how to arrange the basis in rows and cols
            r = 20;

            numScales = length(this.scModel);
            len = size(this.scModel{1}.basisHist, 3);  %# of trials saved
            basisTrack = cell(numScales, len);          %variable to store all the saved basis


            for scale = 1:numScales
                if len == 1
                    basisTrack{scale, 1} = this.scModel{scale}.basis;
                else
                    for j = 1:len
                        basisTrack{scale, j} = this.scModel{scale}.basisHist(:, :, j);
                    end
                end
            end

            h = figure(1);
            scrsz = get(0,'ScreenSize');
            set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);

            % loop over scales
            for s = 1:numScales
                % sort basis according to left energy norm
                endBasis = basisTrack{s,end}(1:end/2,:);
                leftEnergy = abs(sum(endBasis.^2)-0.5);
                [~,I] = sort(leftEnergy);

                subplot(1,numScales,s);
                [di,num] = size(basisTrack{s,1});

                fun1 = @(blc_struct) padarray(padarray(reshape(permute(padarray(reshape(blc_struct.data, sqrt(di / 2), ...
                         sqrt(di / 2), 2), [1, 1], 'pre'), [1, 3, 2]), (sqrt(di / 2) + 1) * 2, sqrt(di / 2) + 1), [1, 1], ...
                         'post') - 1, [1 1], 'pre') + 1;

                for j = 1:len
                    A = basisTrack{s,j}(:,I);
                    % B = reshape(A,di*r,c);
                    B = reshape(A,di * r,num / r); %hotfix!
                    B = B/max(max(abs(B))) + 0.5;
                    C = padarray(padarray(blockproc(B,[di,1],fun1)-1,[1 1],'post')+1,[2,2]);
                    imshow(C);
                    title(sprintf('Basis functions at\n%d%% of training\n(scale %d)', ceil((j / len) * 100), s))
                    % title(num2str(this.trainTime*0.1*(j-1)));
                    drawnow; pause(.001);
                end
            end

            if savePlot
                saveas(h, sprintf('%s/basisFunctions.png', this.savePath), 'png');
            end
        end
    end
end

%% Not overlapping Patch generation
% function patchesNoOv = preprocessImageNoOv(img, fovea, downSampling, patchSize)
%     img = .2989 * img(:,:,1) + .5870 * img(:,:,2) + .1140 * img(:,:,3);
%     for i = 1:log2(downSampling)
%         img = impyramid(img, 'reduce');
%     end

%     % convert to double
%     img = double(img);

%     % cut fovea in the center
%     [h, w, ~] = size(img);
%     img = img(fix(h / 2 + 1 - fovea / 2) : fix(h / 2 + fovea / 2), ...
%               fix(w / 2 + 1 - fovea / 2) : fix(w / 2 + fovea / 2));

%     % cut patches and store them as col vectors
%     % no overlapping patches (for display)
%     patchesNoOv = im2col(img, [patchSize patchSize], 'distinct');
% end

%% Generation of random vergence angles according to truncated Laplace distribution
% function l = truncLaplacian(diversity, range)
%     % see wikipedia for the generation of random numbers according to the
%     % LaPlace distribution via the inversion method
%     r = rand;

%     switch r < 0.5
%         case 1
%             l = 1 / diversity * log(2 * r);
%         case 0
%             l = -1 / diversity * log(2 * (1 - r));
%     end

%     if (abs(l) > range)
%         l = 0;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function trains a model with predifined basis functions as given in
% the model file. These are not trained, but just passively generate a
% feature vector.

function TrainWithBF(trainTime,randomizationSeed,description,pathToBFModel)
    modelBF = load(sprintf('%s/model.mat', pathToBFModel), 'model');
    rng(randomizationSeed);
    learnedFile = '';
    textureFile = 'Textures_New.mat';
    sparseCodingType = 'nonhomeo';

    % Plotting flag
    % Whether figures should be generated, saved and plotted
    % plotIt:   0 = no plot
    %           1 = plot
    plotIt = uint8(1);

    % Save model and conduct testing every saveInterval training iterations (+1)
    saveInterval = 1000;
    if (trainTime < saveInterval)
        saveInterval = trainTime;
    end

    % Instantiate and initiate model and test_data objects
    model = config(learnedFile, textureFile, trainTime, sparseCodingType);

    if (trainTime <= model.interval)
        sprintf('trainTime[%d] must be > model.interval[%d]', trainTime, model.interval)
        return;
    elseif (~exist(fullfile(cd, 'checkEnvironment'), 'file'))
        sprintf('Rendering binary \"checkEnvironment\" not present in current dir\n\"%s\"', cd)
        return;
    end

    % File management
    savePath = sprintf('model_%s_%i_%i_%i_%s_%i_%s', ...
                        datestr(now), ...
                        trainTime, ...
                        sparseCodingType, ...
                        randomizationSeed, ...
                        description);
    % folder = '~/projects/RESULTS/';
    folder = './results/';
    mkdir(folder, savePath);
    savePath = strcat(folder, savePath);

    % Image process variables
    patchSize = 8;

    dsRatioL = modelBF.model.scmodel_Large.Dsratio; %downsampling ratio (Large scale) | original 8
    dsRatioS = modelBF.model.scmodel_Small.Dsratio; %downsampling ratio (Small scale) | original 2

    % fovea = [128 128];
    foveaL = patchSize + patchSize^2 / 2^log2(dsRatioL); %fovea size (Large scale) | 16
    foveaS = patchSize + patchSize^2 / 2^log2(dsRatioS); %fovea size (Small scale) | 40

    stOvL = patchSize / dsRatioL; %steps of overlap in the ds image | 1
    stOvS = patchSize / dsRatioS; %steps of overlap in the ds image | 4

    ncL = foveaL - patchSize + 1; %number of patches per column (slide of 1 px) | 9
    ncS = foveaS - patchSize + 1; %number of patches per column (slide of 1 px) | 33

    % Prepare index matricies for image patches
    columnIndL = [];
    for kc = 1:stOvL:ncL
        tmpInd = (kc - 1) * ncL + 1 : stOvL : kc * ncL;
        columnIndL = [columnIndL tmpInd];
    end
    columnIndS = [];
    for kc = 1:stOvS:ncS
        tmpInd = (kc - 1) * ncS + 1 : stOvS : kc * ncS;
        columnIndS = [columnIndS tmpInd];
    end

    % Camera parameters
    % offset = 0;       %vertical offset between left and right (0 in the Simulator!!!)
    f = 257.34;         %focal length [px]
    baseline = 0.056;   %interocular distance (baseline)

    % Textures
    texturePath = sprintf('config/%s', textureFile);
    texture = load(texturePath);
    texture = texture.texture;
    nTextures = length(texture);
    % currentTexture = texture{1}; %choose first texture as initial

    % Object distance to eyes [m]
    objDistMin = 0.5;
    objDistMax = 2;
    % objDist = objDistMax; %object init position

    degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
    metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

    command = [0, 0];
    t = 0;

    %%% Helper function that maps muscle activities to resulting angle
    function [angle] = getAngle(command)
        cmd = (command * 10) + 1;                               % scale commands to table entries
        angle = interp2(degrees.results_deg, cmd(1), cmd(2));   % interpolate in tabular
    end

    %%% Helper function that maps muscle activities to resulting metabolic costs
    function [tmpMetCost] = getMetCost(command)
        cmd = (command * 10) + 1;                               % scale commands to table entries
        tmpMetCost = interp2(metCosts.results, cmd(1), cmd(2)); % interpolate in tabular
    end

    %%% Main execution loop
    tic % start time count
    for iter1 = 1 : (model.trainTime / model.interval)

        % pick random texture every #interval times
        currentTexture = texture{(randi(nTextures, 1))};

        % random depth
        objDist = objDistMin + (objDistMax - objDistMin) * rand(1, 1);

        % reset muscle activities to random values
        % initialization for muscle in between borders of desired actvity
        % i.e. min and max stimulus distance
        command(2) = 0.00807 + (0.07186 - 0.00807) * rand(1,1);
        % command(2) = 0.1 * rand(1, 1); % random policy

        angleNew = getAngle(command) * 2;
        [status, res] = system(sprintf('./checkEnvironment %s %d %d left.png right.png', ...
                                       currentTexture, objDist, angleNew));

        % Abort execution if error occured
        if (status)
            sprintf('Error in checkEnvironment:\n%s', res)
            return;
        end

        for iter2 = 1 : model.interval
            t = t + 1;
            % Read input images and convert to gray scale
            imgRawLeft = imread('left.png');
            imgRawRight = imread('right.png');
            imgGrayLeft = .2989 * imgRawLeft(:,:,1) + .5870 * imgRawLeft(:,:,2) + .1140 * imgRawLeft(:,:,3);
            imgGrayRight = .2989 * imgRawRight(:,:,1) + .5870 * imgRawRight(:,:,2) + .1140 * imgRawRight(:,:,3);

            % anaglyph = stereoAnaglyph(imgGrayLeft, imgGrayRight);
            % imwrite(anaglyph, 'anaglyph.png');

            % Image patch generation: left{small scale, large scale}, right{small scale, large scale}
            [patchesLeftSmall] = preprocessImage(imgGrayLeft, foveaS, dsRatioS, patchSize, columnIndS);
            [patchesLeftLarge] = preprocessImage(imgGrayLeft, foveaL, dsRatioL, patchSize, columnIndL);
            [patchesRightSmall] = preprocessImage(imgGrayRight, foveaS, dsRatioS, patchSize, columnIndS);
            [patchesRightLarge] = preprocessImage(imgGrayRight, foveaL, dsRatioL, patchSize, columnIndL);

            % Image patches matrix (input to model)
            currentView = {[patchesLeftLarge; patchesRightLarge] [patchesLeftSmall; patchesRightSmall]};

            % Generate input feature vector from current images
            [feature, reward, errorTotal, errorLarge, errorSmall] = modelBF.model.generateFR(currentView);

            % % incorporationg the current muscle activity into feature vector
            % % and scaling it to the value
            % % range of BF activations
            % feature = [feature; command' * model.lambdaMuscleFB];

            %%% Feedback
        % Absolute command feedback # concatination
        feature = [feature; command(2) * model.lambdaMuscleFB];
        % Relative command feedback # concatination
        % if (iter2 > 1)
        %     feature = [feature; model.relCmd_hist(t-1) * model.lambdaMuscleFB];
        % else
        %     feature = [feature; 0];
        % end

        %% Absolute command feedback # additive
        % feature = feature + command(2) * model.lambdaMuscleFB;
        %% Absolute command feedback # multiplicative
        % feature = feature * (command(2) * model.lambdaMuscleFB);
        %% Relative command feedback # additive
        % if (iter2 > 1)
        %     feature = feature + model.relCmd_hist(t - 1) * model.lambdaMuscleFB;
        % end
        %% Relative command feedback # multiplicative
        % if (iter2 > 1)
        %     feature = feature * model.relCmd_hist(t - 1) * model.lambdaMuscleFB;
        % end

        %%% Calculate metabolic costs
        metCost = getMetCost(command) * 2;

        %%% Calculate reward function
        % rewardFunction = (model.lambdaMet * reward) + ((1 - model.lambdaMet) * - metCost);

        %%% Weight L1 regularization
        rewardFunction = model.lambdaRec * reward ...
                         - model.lambdaMet * metCost ...
                         - model.lambdaV * (sum(sum(abs(model.rlmodel.CCritic.v_ji)))) ...
                         - model.lambdaP1 * (sum(sum(abs(model.rlmodel.CActor.wp_ji)))) ...
                         - model.lambdaP2 * (sum(sum(abs(model.rlmodel.CActor.wp_kj))));

        %%% Weight L2 regularization
        % rewardFunction = model.lambdaRec * reward ...
        %                  - model.lambdaMet * metCost ...
        %                  - model.lambdaV * (sum(sum(model.rlmodel.CCritic.v_ji .^ 2))) ...
        %                  - model.lambdaP1 * (sum(sum(model.rlmodel.CActor.wp_ji .^ 2))) ...
        %                  - model.lambdaP2 * (sum(sum(model.rlmodel.CActor.wp_kj .^ 2)));

        %%% Learning
        % Sparse coding models
            %%% Learning
            % Sparse coding models
%             model.scmodel_Large.stepTrain(currentView{1});
%             model.scmodel_Small.stepTrain(currentView{2});
            % RL model
            relativeCommand = model.rlmodel.stepTrain(feature, rewardFunction, (iter2 > 1));

            % add the change in muscle Activities to current ones
            command(2) = command(2) + relativeCommand';
            command = checkCmd(command);            %restrain motor commands to [0,1]
            angleNew = getAngle(command) * 2;       %resulting angle is used for both eyes

            % generate new view (two pictures) with new vergence angle
            [status, res] = system(sprintf('./checkEnvironment %s %d %d left.png right.png', ...
                                           currentTexture, objDist, angleNew));

            % Abort execution if error occured
            if (status)
                sprintf('Error in checkEnvironment:\n%s', res)
                return;
            end

            %%%%%%%%%%%%%%%% TRACK ALL PARAMETERS %%%%%%%%%%%%%%%%%%

            %Compute desired vergence command, disparity and vergence error
            fixDepth = (0.5 * baseline) / tand(angleNew / 2);
            angledes = 2 * atan(baseline / (2 * objDist));                  %desired vergence [rad]
            anglerr = angledes * 180 / pi - angleNew;                       %vergence error [deg]
            disparity = 2 * f * tan((angledes - angleNew * pi / 180) / 2);  %current disp [px]

            %save them
            model.Z(t) = objDist;
            model.fixZ(t) = fixDepth;
            model.disp_hist(t) = disparity;
            model.vergerr_hist(t) = anglerr;
            model.recerr_hist(t, :) = [errorLarge; errorSmall];
            model.verge_actual(t) = angleNew;
            model.relCmd_hist(t) = relativeCommand;
            model.cmd_hist(t, :) = command;
            model.reward_hist(t) = rewardFunction;
            % model.feature_hist(t, :) = feature;
            model.metCost_hist(t) = metCost;
            model.td_hist(t) = model.rlmodel.CCritic.delta;
            model.g_hist(t) = model.rlmodel.CActor.params(7);
            model.weight_hist(t, 1) = model.rlmodel.CCritic.params(1);
            model.weight_hist(t, 2) = model.rlmodel.CCritic.params(2);
            model.weight_hist(t, 3) = model.rlmodel.CActor.params(1);
            model.weight_hist(t, 4) = model.rlmodel.CActor.params(2);
            model.weight_hist(t, 5) = model.rlmodel.CActor.params(3);
            model.weight_hist(t, 6) = model.rlmodel.CActor.params(4);
            model.weight_hist(t, 7) = model.rlmodel.CActor.params(5);
            model.weight_hist(t, 8) = model.rlmodel.CActor.params(6);
        end

        sprintf('Training Iteration = %d\nCommand = [%.3g,\t%.3g]\tCurrent Vergence = %.3g\nRec Error = %.3g\tVergence Error =\n[%.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g]', ...
            t, command(1), command(2), angleNew, errorTotal, model.vergerr_hist(t - model.interval + 1 : t))

        % Display per cent completed of training and save model
        if (~mod(t, saveInterval))
            sprintf('%g%% is finished', (t / model.trainTime * 100))
            save(strcat(savePath, '/model'), 'model');

            %save Basis
            model.scmodel_Large.saveBasis;
            model.scmodel_Small.saveBasis;

            %save Weights
            % model.rlmodel.saveWeights; %save policy and value net weights
        end
    end
    elapsedTime = toc;

    % Total simulation time
    model.simulatedTime = elapsedTime / 60;
    sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
            elapsedTime / 3600, elapsedTime / 60, elapsedTime, trainTime / elapsedTime)

    % Plot results
    if (plotIt)
        % model.errPlot();
        model.allPlotSave(savePath);
    end

    % Save results data
    save(strcat(savePath, '/model'), 'model');

    %%% Testing procedure
    %%% TODO: implement function or script call here
end

%%% Saturation function that keeps motor commands in [0, 1]
%   corresponding to the muscelActivity/metabolicCost tables
function [cmd] = checkCmd(cmd)
    i0 = cmd < 0;
    cmd(i0) = 0;
    i1 = cmd > 1;
    cmd(i1) = 1;
end

%%% Helper functions for image preprocessing
%% Patch generation
function [patches] = preprocessImage(img, fovea, downSampling, patchSize, columnIndicies)
    % img = .2989 * img(:,:,1) + .5870 * img(:,:,2) + .1140 * img(:,:,3);
    for i = 1:log2(downSampling)
        img = impyramid(img, 'reduce');
    end

    % convert to double
    img = double(img);

    % cut fovea in the center
    [h, w, ~] = size(img);
    img = img(fix(h / 2 + 1 - fovea / 2) : fix(h / 2 + fovea / 2), ...
              fix(w / 2 + 1 - fovea / 2) : fix(w / 2 + fovea / 2));

    % cut patches and store them as col vectors
    patches = im2col(img, [patchSize patchSize], 'sliding');            %slide window of 1 px

    % take patches at steps of s (8 px)
    patches = patches(:, columnIndicies);                               %81 patches

    % pre-processing steps (0 mean, unit norm)
    patches = patches - repmat(mean(patches), [size(patches, 1) 1]);    %0 mean
    normp = sqrt(sum(patches.^2));                                      %patches norm

    % normalize patches to norm 1
    normp(normp == 0) = eps;                                            %regularizer
    patches = patches ./ repmat(normp, [size(patches, 1) 1]);           %normalized patches
end

%% Not overlapping Patch generation
function patchesNoOv = preprocessImageNoOv(img, fovea, downSampling, patchSize)
    img = .2989 * img(:,:,1) + .5870 * img(:,:,2) + .1140 * img(:,:,3);
    for i = 1:log2(downSampling)
        img = impyramid(img, 'reduce');
    end

    % convert to double
    img = double(img);

    % cut fovea in the center
    [h, w, ~] = size(img);
    img = img(fix(h / 2 + 1 - fovea / 2) : fix(h / 2 + fovea / 2), ...
              fix(w / 2 + 1 - fovea / 2) : fix(w / 2 + fovea / 2));

    % cut patches and store them as col vectors
    % no overlapping patches (for display)
    patchesNoOv = im2col(img, [patchSize patchSize], 'distinct');
end

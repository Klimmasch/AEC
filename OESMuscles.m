%%% Main script for launching experimental procedure
% @param trainTime           training time in number of iterations
% @param randomizationSeed   randomization seed
% @param fileDescription     description of approach used as file name
%%%
function OESMuscles(trainTime, randomizationSeed, fileDescription)
    rng(randomizationSeed);

    % useLearnedFile: [use previously learned policy specified in learnedFile,
    %                  whether the training shall be completed, or retrained with same/new parameters]
    useLearnedFile = [0, 0];
    learnedFile = '';
    % learnedFile = '/home/klimmasch/projects/results/model_13-Apr-2016_14:04:55_100_nonhomeo_2_testTrainOn/model.mat';
    % learnedFile = '/home/lelais/Documents/MATLAB/results/model_18-Apr-2016_18:27:54_200000_nonhomeo_1_CACLAVar_NewHiddenUpdate_init00017-01-004_alpha10_var-5/model.mat';

    %%% Stimulus declaration
    textureFile = 'Textures_vanHaterenTrain.mat';   % vanHateren database
    % textureFile = 'Textures_celine.mat';          % Celine's images

    % Sparse coding approach
    % sparseCodingType: 0 = non-homeostatic
    %                   1 = homeostatic
    sparseCodingType = uint8(0);
    sparseCodingTypeName = cellstr(['nonhomeo'; 'homeo   ']);

    % Plotting flag
    % Whether figures should be generated and saved
    % plotIt: [training, test]
    %            0 = don't do it
    %            1 = do it
    plotIt = [uint8(1), uint8(1)];

    % Testing flag
    % Whether the testing procedure shall be executed after training
    % testIt:   0 = don't do it
    %           1 = do it
    testIt = uint8(1);

    % Load model from file or instantiate and initiate new model object
    if (useLearnedFile(1) == 1)
        if isempty(learnedFile)
            sprintf('could not open the learned file! %s', learnedFile)
            return;
        else
            model = load(learnedFile, 'model');
            model = model.model;
        end
    else
        model = config(textureFile, trainTime, sparseCodingType);
    end

    % safety check for plotting functions
    if (trainTime <= model.interval)
        sprintf('trainTime[%d] must be > model.interval[%d]', trainTime, model.interval)
        return;
    end

    % File management: either complete training with existing folder etc.,
    % or create a new one
    if ((useLearnedFile(1) == 1) && (useLearnedFile(2) == 1))
        timeToTrain = model.trainTime - model.trainedUntil;
    else
        modelName = sprintf('model_%s_%i_%s_%i_%s', ...
                            datestr(now, 'dd-mmm-yyyy_HH:MM:SS'), ...
                            trainTime, ...
                            sparseCodingTypeName{sparseCodingType + 1}, ...
                            randomizationSeed, ...
                            fileDescription);
        folder = '../results/';
        mkdir(folder, modelName);
        model.savePath = strcat(folder, modelName);

        % backup all used files
        copyfile(strcat(mfilename, '.m'), model.savePath);
        copyfile('config.m', model.savePath);
        copyfile(strcat(class(model), '.m'), model.savePath);
        copyfile(strcat(class(model.rlmodel), '.m'), model.savePath);
        copyfile(strcat(class(model.rlmodel.CCritic), '.m'), model.savePath);
        copyfile(strcat(class(model.rlmodel.CActor), '.m'), model.savePath);

        timeToTrain = model.trainTime;
    end
    % additional notes/infromation to this model/approach
    model.notes = [model.notes fileDescription];

    % Save model every #saveInterval training iterations
    saveInterval = ceil(model.trainTime / 5);

    % Textures
    texture = load(sprintf('config/%s', textureFile));
    texture = texture.texture;
    nTextures = length(texture);

    degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
    metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

    % muscle function :=  mf(vergence_angle) = muscle force [single muscle]
    resolution = 10001;
    approx = spline(1:11, degrees.results_deg(:, 1));

    xValPos = ppval(approx, 1:0.001:11)';
    yValPos = linspace(0, 1, resolution)';

    xValNeg = flipud(ppval(approx, 1:0.001:11)' * -1);
    yValNeg = linspace(-1, 0, resolution)';

    mfunction = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
    mfunction(:, 1) = mfunction(:, 1) * 2;  % angle for two eyes
    dmf = diff(mfunction(1 : 2, 1));        % delta in angle
    dmf2 = diff(mfunction(1 : 2, 2));       % delta in mf

    %%% New renderer
    simulator = OpenEyeSim('create');
    simulator.initRenderer();
    % simulator.reinitRenderer(); % for debugging

    imgRawLeft = uint8(zeros(240, 320, 3));
    imgRawRight = uint8(zeros(240, 320, 3));

    % Generates two new images for both eyes
    function refreshImages(texture, vergAngle, objDist)
        simulator.add_texture(1, texture);
        simulator.set_params(1, vergAngle, objDist);

        result1 = simulator.generate_left();
        result2 = simulator.generate_right();

        k = 1;
        l = 1;
        for i = 1 : 3 : length(result1)
            imgRawLeft(k, l, 1) = result1(i);
            imgRawLeft(k, l, 2) = result1(i + 1);
            imgRawLeft(k, l, 3) = result1(i + 2);

            imgRawRight(k, l, 1) = result2(i);
            imgRawRight(k, l, 2) = result2(i + 1);
            imgRawRight(k, l, 3) = result2(i + 2);

            l = l + 1;
            if (l > 320)
                l = 1;
                k = k + 1;
            end
        end
    end

    %%% Helper function that maps {vergenceAngle} -> {muscleForce}
    function mf = getMF(vergAngle)
        % look up index of vergAngle
        indVergAngle = find(mfunction(:, 1) <= vergAngle + dmf & mfunction(:, 1) >= vergAngle - dmf);
        mf = mfunction(indVergAngle, 2);
        mf = mf(ceil(length(mf) / 2));
    end

    %%% Helper function that maps muscle activities to resulting angle
    function [angle] = getAngle(command)
        cmd = (command * 10) + 1;                               % scale commands to table entries
        angle = interp2(degrees.results_deg, cmd(1), cmd(2));   % interpolate in tabular
    end

    function angle = getAngle2(command)
        angleIndex = find(mfunction(:, 2) <= command(2) + dmf2 & mfunction(:, 2) >= command(2) - dmf2);
        angle = mfunction(angleIndex, 1);
        angle = angle(ceil(length(angle) / 2));
    end

    %%% Helper function that maps muscle activities to resulting metabolic costs
    function [tmpMetCost] = getMetCost(command)
        cmd = (command * 10) + 1;                               % scale commands to table entries
        tmpMetCost = interp2(metCosts.results, cmd(1), cmd(2)); % interpolate in tabular
    end

    %%% Saturation function that keeps motor commands in [0, 1]
    %   corresponding to the muscelActivity/metabolicCost tables
    function [cmd] = checkCmd(cmd)
        i0 = cmd < 0;
        cmd(i0) = 0;
        i1 = cmd > 1;
        cmd(i1) = 1;
    end

    %%% Helper function for image preprocessing
    %% Patch generation
    function [patches] = preprocessImage(img, imgSize)
        % img = .2989 * img(:,:,1) + .5870 * img(:,:,2) + .1140 * img(:,:,3);

        if (imgSize == 0)
            % small scale
            for i = 1:log2(model.dsRatioS)
                img = impyramid(img, 'reduce');
            end

            % convert to double
            img = double(img);

            % cut fovea in the center
            [h, w, ~] = size(img);
            img = img(fix(h / 2 + 1 - model.foveaS / 2) : fix(h / 2 + model.foveaS / 2), ...
                      fix(w / 2 + 1 - model.foveaS / 2) : fix(w / 2 + model.foveaS / 2));
        else
            % large scale
            for i = 1:log2(model.dsRatioL)
                img = impyramid(img, 'reduce');
            end

            % convert to double
            img = double(img);

            % cut fovea in the center
            [h, w, ~] = size(img);
            img = img(fix(h / 2 + 1 - model.foveaL / 2) : fix(h / 2 + model.foveaL / 2), ...
                      fix(w / 2 + 1 - model.foveaL / 2) : fix(w / 2 + model.foveaL / 2));
        end

        % cut patches and store them as col vectors
        patches = im2col(img, [model.patchSize model.patchSize], 'sliding');            %slide window of 1 px

        % take patches at steps of s (8 px)
        if (imgSize == 0)
            patches = patches(:, model.columnIndS);                               %81 patches
        else
            patches = patches(:, model.columnIndL);                               %81 patches
        end

        % pre-processing steps (0 mean, unit norm)
        patches = patches - repmat(mean(patches), [size(patches, 1) 1]);    %0 mean
        normp = sqrt(sum(patches.^2));                                      %patches norm

        % normalize patches to norm 1
        normp(normp == 0) = eps;                                            %regularizer
        patches = patches ./ repmat(normp, [size(patches, 1) 1]);           %normalized patches
    end

    %%% Main execution loop
    t = model.trainedUntil; % this is zero in new initiated model
    command = [0, 0];
    % rewardFunction_prev = 0;
    tic; % start time count
    for iter1 = 1 : (timeToTrain / model.interval)
        % pick random texture every #interval times
        currentTexture = texture{(randi(nTextures, 1))};

        % random depth
        objDist = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);

        % reset muscle activities to random values
        % initialization for muscle in between borders of desired actvity
        % i.e. min and max stimulus distance
        command(1) = 0; % single muscle
        % command(1) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1, 1); % two muscles
        % command(2) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1, 1);
        % command(2) = getMF(model.vergAngleMin + (model.vergAngleMax - model.vergAngleMin) * rand(1, 1));
        initDist = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);
        command(2) = getMF(2 * atand(model.baseline / (2 * initDist)));

        % testing input distribution
        % nSamples = 10000;
        % commands = zeros(nSamples,1);
        % angles = zeros(nSamples, 1);
        % dists = zeros(nSamples, 1);
        % for i = 1:10000
        %     initDist = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);
        %     initAngle = atand(model.baseline / (2 * initDist));
        %     commands(i) = getMF(initAngle*2);
        %     angles(i) = getAngle([0, commands(i)]);
        %     dists(i) = (model.baseline/ (2 * tand(angles(i))));
        % end
        % figure; histogram(commands); title('commands');
        % figure; histogram(angles); title('angles');
        % figure; histogram(dists); title('distances');

        % angleNew = getAngle(command) * 2;
        angleNew = getAngle2(command);

        for iter2 = 1 : model.interval
            t = t + 1;

            % update stimuli
            refreshImages(currentTexture, angleNew / 2, objDist);

            % convert images to gray scale
            imgGrayLeft = .2989 * imgRawLeft(:,:,1) + .5870 * imgRawLeft(:,:,2) + .1140 * imgRawLeft(:,:,3);
            imgGrayRight = .2989 * imgRawRight(:,:,1) + .5870 * imgRawRight(:,:,2) + .1140 * imgRawRight(:,:,3);

            % Generate & save the anaglyph picture
            % anaglyph = stereoAnaglyph(imgGrayLeft, imgGrayRight); % only for matlab 2015 or newer
            % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [model.savePath '/anaglyph.png']); %this one works for all tested matlab
            % more advanced functions that generated the anaglyphs of the foveal views
            % generateAnaglyphs(imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS, model.savePath);

            % Image patch generation: left{small scale, large scale}, right{small scale, large scale}
            [patchesLeftSmall] = preprocessImage(imgGrayLeft, 0);
            [patchesLeftLarge] = preprocessImage(imgGrayLeft, 1);
            [patchesRightSmall] = preprocessImage(imgGrayRight, 0);
            [patchesRightLarge] = preprocessImage(imgGrayRight, 1);

            % Image patches matrix (input to model)
            currentView = {[patchesLeftLarge; patchesRightLarge] [patchesLeftSmall; patchesRightSmall]};

            % Generate input feature vector from current images
            [feature, reward, ~, errorLarge, errorSmall] = model.generateFR(currentView);

            %%% Feedback
            % Absolute command feedback # concatination
            feature = [feature; command(2) * model.lambdaMuscleFB];

            %%% Calculate metabolic costs
            metCost = getMetCost(command) * 2;

            %%% Calculate reward function
            %% Standard reward
            rewardFunction = model.lambdaRec * reward - model.lambdaMet * metCost;
            % rewardFunction = (model.lambdaMet * reward) + ((1 - model.lambdaMet) * - metCost);

            %% Delta reward
            % rewardFunctionReal = model.lambdaRec * reward - model.lambdaMet * metCost;
            % rewardFunction = rewardFunctionReal - rewardFunction_prev;

            % Stasis punishment, i.e. punish non-movement of eyes
            % if (abs(rewardFunctionReal - rewardFunction_prev) < 1e-5)
            %     rewardFunction = rewardFunctionReal - rewardFunction_prev - 1e-5;
            % end

            % rewardFunction_prev = rewardFunctionReal;

            %%% Weight L1 regularization
            % rewardFunction = model.lambdaRec * reward ...
            %                  - model.lambdaMet * metCost ...
            %                  - model.lambdaV * (sum(sum(abs(model.rlmodel.CCritic.v_ji)))) ...
            %                  - model.lambdaP1 * (sum(sum(abs(model.rlmodel.CActor.wp_ji)))) ...
            %                  - model.lambdaP2 * (sum(sum(abs(model.rlmodel.CActor.wp_kj))));

            %%% Weight L2 regularization
            % rewardFunction = model.lambdaRec * reward ...
            %                  - model.lambdaMet * metCost ...
            %                  - model.lambdaV * (sum(sum(model.rlmodel.CCritic.v_ji .^ 2))) ...
            %                  - model.lambdaP1 * (sum(sum(model.rlmodel.CActor.wp_ji .^ 2))) ...
            %                  - model.lambdaP2 * (sum(sum(model.rlmodel.CActor.wp_kj .^ 2)));

            %%% Learning
            % Sparse coding models
            model.scmodel_Large.stepTrain(currentView{1});
            model.scmodel_Small.stepTrain(currentView{2});

            % RL model
            % Variance decay, i.e. reduction of actor's output perturbation
            % if ((model.rlmodel.continuous == 1) && (model.rlmodel.CActor.varDec > 0))
            %     model.rlmodel.CActor.variance = model.rlmodel.CActor.varianceRange(1) * 2 ^ (-t / model.rlmodel.CActor.varDec);
            % end

            relativeCommand = model.rlmodel.stepTrain(feature, rewardFunction, (iter2 > 1));

            % add the change in muscle Activities to current ones
            % command = command + relativeCommand';     %two muscles
            command(1) = 0;
            command(2) = command(2) + relativeCommand;  %one muscle
            command = checkCmd(command);                %restrain motor commands to [0,1]

            % angleNew = getAngle(command) * 2; %resulting angle is used for both eyes
            angleNew = getAngle2(command);

            %%%%%%%%%%%%%%%% TRACK ALL PARAMETERS %%%%%%%%%%%%%%%%%%

            % compute desired vergence command, disparity and vergence error
            fixDepth = (model.baseline / 2) / tand(angleNew / 2);   %fixation depth [m]
            angleDes = 2 * atand(model.baseline / (2 * objDist));   %desired vergence [deg]
            anglerr = angleDes - angleNew;                          %vergence error [deg]
            disparity = 2 * model.focalLength * tand(anglerr / 2);  %current disp [px]

            % save state
            model.Z(t) = objDist;
            model.fixZ(t) = fixDepth;
            model.disp_hist(t) = disparity;
            model.vergerr_hist(t) = anglerr;
            model.recerr_hist(t, :) = [errorLarge; errorSmall];
            model.verge_actual(t) = angleNew;
            model.verge_desired(t) = angleDes;
            model.relCmd_hist(t) = relativeCommand;
            model.cmd_hist(t, :) = command;
            model.reward_hist(t) = rewardFunction;
            % model.feature_hist(t, :) = feature;
            model.metCost_hist(t) = metCost;
            model.td_hist(t) = model.rlmodel.CCritic.delta;
            model.weight_hist(t, 1) = model.rlmodel.CCritic.params(1);
            model.weight_hist(t, 2) = model.rlmodel.CActor.params(1);
            model.weight_hist(t, 3) = model.rlmodel.CActor.params(2);
            model.weight_hist(t, 4) = model.rlmodel.CActor.params(3);
            model.variance_hist(t) = model.rlmodel.CActor.variance;

            model.trainedUntil = t;
        end

        sprintf('Training Iteration = %d\nAbs Command =\t[%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f]\nRel Command = \t[%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f]\nVer Error =\t[%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f]', ...
                t, model.cmd_hist(t - model.interval + 1 : t, 2), model.relCmd_hist(t - model.interval + 1 : t), model.vergerr_hist(t - model.interval + 1 : t))

        % Display per cent completed of training and save model
        if (~mod(t, saveInterval))
            sprintf('%g%% is finished', (t / timeToTrain * 100))
            save(strcat(model.savePath, '/model'), 'model');

            % save Basis
            model.scmodel_Large.saveBasis();
            model.scmodel_Small.saveBasis();
        end
    end
    elapsedTime = toc;

    % Total simulation time
    model.simulatedTime = elapsedTime / 60;
    sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
            elapsedTime / 3600, elapsedTime / 60, elapsedTime, trainTime / elapsedTime)

    % store simulated time
    % if useLearnedFile(2)
    %     model.trainTime = model.trainTime + model.trainedUntil;
    % end
    save(strcat(model.savePath, '/model'), 'model');

    % plot results
    if (plotIt(1) == 1)
        model.allPlotSave();
    end

    %%% Testing procedure
    if (testIt)
        % testModel(model, randomizationSeed, objRange, vergRange, repeat, randStimuli, randObjRange, plotIt, saveTestResults)
        % testModel(model, randomizationSeed, [0.5, 1, 1.5, 2], [-3 : 0.2 : 3], [50, 50], 0, 1, plotIt(2), 1);
        % testModel2(model, nStim, plotIt, saveTestResults)
        testModel2(model, 50, plotIt(2), 1, 1);
    end
end

%this function generates anaglyphs of the large and small scale fovea and
%one of the two unpreprocessed gray scale images
function generateAnaglyphs(leftGray, rightGray, dsRatioL, dsRatioS, foveaL, foveaS, savePath)
    anaglyph = imfuse(leftGray, rightGray, 'falsecolor');
    imwrite(anaglyph, [savePath '/anaglyph.png']);

    %Downsampling Large
    imgLeftL = leftGray(:);
    imgLeftL = reshape(imgLeftL, size(leftGray));
    imgRightL = rightGray(:);
    imgRightL = reshape(imgRightL, size(rightGray));
    for i = 1:log2(dsRatioL)
        imgLeftL = impyramid(imgLeftL, 'reduce');
        imgRightL = impyramid(imgRightL, 'reduce');
    end

    % cut fovea in the center
    [h, w, ~] = size(imgLeftL);
    imgLeftL = imgLeftL(fix(h / 2 + 1 - foveaL / 2) : fix(h / 2 + foveaL / 2), ...
              fix(w / 2 + 1 - foveaL / 2) : fix(w / 2 + foveaL / 2));
    imgRightL = imgRightL(fix(h / 2 + 1 - foveaL / 2) : fix(h / 2 + foveaL / 2), ...
              fix(w / 2 + 1 - foveaL / 2) : fix(w / 2 + foveaL / 2));

    %create an anaglyph of the two pictures, scale it up and save it
    anaglyphL = imfuse(imgLeftL, imgRightL, 'falsecolor');
    imwrite(imresize(anaglyphL, 20), [savePath '/anaglyphLargeScale.png']);
    largeScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(largeScaleView, 20), [savePath '/LargeScaleMontage.png']);

    %Downsampling Small
    imgLeftS = leftGray(:);
    imgLeftS = reshape(imgLeftS, size(leftGray));
    imgRightS = rightGray(:);
    imgRightS = reshape(imgRightS, size(rightGray));
    for i = 1:log2(dsRatioS)
        imgLeftS = impyramid(imgLeftS, 'reduce');
        imgRightS = impyramid(imgRightS, 'reduce');
    end

    % cut fovea in the center
    [h, w, ~] = size(imgLeftS);
    imgLeftS = imgLeftS(fix(h / 2 + 1 - foveaS / 2) : fix(h / 2 + foveaS / 2), ...
              fix(w / 2 + 1 - foveaS / 2) : fix(w / 2 + foveaS / 2));
    imgRightS = imgRightS(fix(h / 2 + 1 - foveaS / 2) : fix(h / 2 + foveaS / 2), ...
              fix(w / 2 + 1 - foveaS / 2) : fix(w / 2 + foveaS / 2));

    %create an anaglyph of the two pictures, scale it up and save it
    anaglyphS = imfuse(imgLeftS, imgRightS, 'falsecolor');
    imwrite(imresize(anaglyphS, 8), [savePath '/anaglyphSmallScale.png']);
    smallScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(smallScaleView, 8), [savePath '/smallScaleMontage.png']);
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

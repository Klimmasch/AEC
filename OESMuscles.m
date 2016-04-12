%%% Main script for launching experimental procedure
%@param trainTime           training time in number of iterations
%@param randomizationSeed   randomization seed
%@param fileDescription     description of approach used as file name
%
% learnedFile:              previously learned policy and sparse coding model
% textureFile:              texture settings files
% sparseCodingType:         type of sparse coding approach
%%%
function OESMuscles(trainTime, randomizationSeed, fileDescription)

rng(randomizationSeed);
learnedFile = '';
% textureFile = 'Textures_celine.mat';
textureFile = 'Textures_vanHaterenTrain.mat';
sparseCodingType = 'nonhomeo';

% Plotting and saving flag
% Whether figures should be generated, saved and plotted
% additionally the relevant scripts are backed up
% plotNsave: [training, test]
%            0 = don't do it
%            1 = do it
plotNsave = [uint8(1), uint8(1)];

% Testing flag
% Whether the testing procedure shall be executed after training
% testIt:   0 = don't do it
%           1 = do it
testIt = uint8(1);

% Save model every #saveInterval training iterations
saveInterval = 1000;
if (trainTime < saveInterval)
    saveInterval = trainTime;
end

% Instantiate and initiate model object
model = config(learnedFile, textureFile, trainTime, sparseCodingType);
%%%%%%%%%%%%%
% n = 10000;
% commands = zeros(n,2);
% angles = zeros(n,1);
% degrees = load('Degrees.mat');
% for i = 1:n
%     commands(i,2) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1,1);
%     angles(i) = getAngle(commands(i,:))*2;
% end
% histogram(angles); 
% min(angles) % = 1.6045
% max(angles) % = 6.4094
%%%%%%%%%%%%%%%%%
if (trainTime <= model.interval)
    sprintf('trainTime[%d] must be > model.interval[%d]', trainTime, model.interval)
    return;
elseif (~exist(fullfile(cd, 'checkEnvironment'), 'file'))
    sprintf('Rendering binary \"checkEnvironment\" not present in current dir\n\"%s\"', cd)
    return;
end

% File management
modelName = sprintf('model_%s_%i_%s_%i_%s', ...
                    datestr(now, 'dd-mmm-yyyy_HH:MM:SS'), ...
                    trainTime, ...
                    sparseCodingType, ...
                    randomizationSeed, ...
                    fileDescription);
% folder = '~/projects/RESULTS/';
% folder = './results/';
folder = '../results/';
mkdir(folder, modelName);
model.savePath = strcat(folder, modelName);

% Image process variables
patchSize = 8;

dsRatioL = model.scmodel_Large.Dsratio; %downsampling ratio (Large scale) | original 8
dsRatioS = model.scmodel_Small.Dsratio; %downsampling ratio (Small scale) | original 2

% fovea = [128 128];
foveaL = patchSize + patchSize ^ 2 / 2 ^ log2(dsRatioL); %fovea size (Large scale) | 16
foveaS = patchSize + patchSize ^ 2 / 2 ^ log2(dsRatioS); %fovea size (Small scale) | 40

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

% Textures
texturePath = sprintf('config/%s', textureFile);
texture = load(texturePath);
texture = texture.texture;
nTextures = length(texture);

degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

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
t = 0;
command = [0, 0];
rewardFunction_prev = 0;
tic; % start time count
for iter1 = 1 : (model.trainTime / model.interval)

    % pick random texture every #interval times
    currentTexture = texture{(randi(nTextures, 1))};

    % random depth
    objDist = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);

    % reset muscle activities to random values
    % initialization for muscle in between borders of desired actvity
    % i.e. min and max stimulus distance
    command(1) = 0; % single muscle
    % command(1) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1, 1); % two muscles
    command(2) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1, 1);

    angleNew = getAngle(command) * 2;

    [status, res] = system(sprintf('./checkEnvironment %s %d %d left.png right.png', ...
                                   currentTexture, objDist, angleNew));

    % abort execution if error occured
    if (status)
        sprintf('Error in checkEnvironment:\n%s', res)
        return;
    end

    for iter2 = 1 : model.interval
        t = t + 1;
        % read input images and convert to gray scale
        imgRawLeft = imread('left.png');
        imgRawRight = imread('right.png');
        imgGrayLeft = .2989 * imgRawLeft(:,:,1) + .5870 * imgRawLeft(:,:,2) + .1140 * imgRawLeft(:,:,3);
        imgGrayRight = .2989 * imgRawRight(:,:,1) + .5870 * imgRawRight(:,:,2) + .1140 * imgRawRight(:,:,3);

        % Generate & save the anaglyph picture
        % anaglyph = stereoAnaglyph(imgGrayLeft, imgGrayRight); % only for matlab 2015 or newer
        % generateAnaglyphs(imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS);

        % Image patch generation: left{small scale, large scale}, right{small scale, large scale}
        [patchesLeftSmall] = preprocessImage(imgGrayLeft, foveaS, dsRatioS, patchSize, columnIndS);
        [patchesLeftLarge] = preprocessImage(imgGrayLeft, foveaL, dsRatioL, patchSize, columnIndL);
        [patchesRightSmall] = preprocessImage(imgGrayRight, foveaS, dsRatioS, patchSize, columnIndS);
        [patchesRightLarge] = preprocessImage(imgGrayRight, foveaL, dsRatioL, patchSize, columnIndL);

        % Image patches matrix (input to model)
        currentView = {[patchesLeftLarge; patchesRightLarge] [patchesLeftSmall; patchesRightSmall]};

        % Generate input feature vector from current images
        [feature, reward, ~, errorLarge, errorSmall] = model.generateFR(currentView);

        %%% Feedback
        % Absolute command feedback # concatination
        if (model.rlmodel.continuous == 1)
            feature = [feature; command(2) * model.lambdaMuscleFB];
        end

        %%% Calculate metabolic costs
        metCost = getMetCost(command) * 2;

        %%% Calculate reward function
        %% Standard reward
        % rewardFunction = model.lambdaRec * reward - model.lambdaMet * metCost;
        % rewardFunction = (model.lambdaMet * reward) + ((1 - model.lambdaMet) * - metCost);

        %% Delta reward
        rewardFunctionReal = model.lambdaRec * reward - model.lambdaMet * metCost;
        rewardFunction = rewardFunctionReal - rewardFunction_prev;

        % Stasis punishment, i.e. punish non-movement of eyes
        % if (abs(rewardFunctionReal - rewardFunction_prev) < 1e-5)
        %     rewardFunction = rewardFunctionReal - rewardFunction_prev - 1e-5;
        % end

        rewardFunction_prev = rewardFunctionReal;

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
        if (model.rlmodel.continuous == 1)
            model.rlmodel.CActor.variance = model.rlmodel.CActor.varianceRange(1) * 2 ^ (-t / model.rlmodel.CActor.varDec);
        end

        relativeCommand = model.rlmodel.stepTrain(feature, rewardFunction, (iter2 > 1));

        % add the change in muscle Activities to current ones
        % command = command + relativeCommand';     %two muscels
        command(1) = 0;
        command(2) = command(2) + relativeCommand;  %one muscel
        command = checkCmd(command);                %restrain motor commands to [0,1]

        if (model.rlmodel.continuous == 1)
            angleNew = getAngle(command) * 2; %resulting angle is used for both eyes
        else
            angleNew = angleNew + relativeCommand;
            if (angleNew > 71.5 || angleNew < 0.99) % analogous to checkCmd
                command = [0, 0];
                command(2) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1,1);
                angleNew = getAngle(command) * 2;
            end
        end

        % generate new view (two pictures) with new vergence angle
        [status, res] = system(sprintf('./checkEnvironment %s %d %d left.png right.png', ...
                                       currentTexture, objDist, angleNew));

        % abort execution if error occured
        if (status)
            sprintf('Error in checkEnvironment:\n%s', res)
            return;
        end

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
        if (model.rlmodel.continuous == 1)
            model.td_hist(t) = model.rlmodel.CCritic.delta;
            % model.g_hist(t) = model.rlmodel.CActor.params(7);
            model.weight_hist(t, 1) = model.rlmodel.CCritic.params(1);
            model.weight_hist(t, 2) = model.rlmodel.CActor.params(1);
            if (model.rlmodel.rlFlavour(2) >= 4)
                model.weight_hist(t, 3) = model.rlmodel.CActor.params(2);
                if (model.rlmodel.rlFlavour(2) == 5)
                    model.weight_hist(t, 4) = model.rlmodel.CActor.params(3);
                end
            end
            model.variance_hist(t) = model.rlmodel.CActor.variance;
        end
    end

    sprintf('Training Iteration = %d\nAbs Command =\t[%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f]\nRel Command = \t[%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f]\nVer Error =\t[%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f]', ...
            t, model.cmd_hist(t - model.interval + 1 : t, 2), model.relCmd_hist(t - model.interval + 1 : t), model.vergerr_hist(t - model.interval + 1 : t))

    % Display per cent completed of training and save model
    if (~mod(t, saveInterval))
        sprintf('%g%% is finished', (t / model.trainTime * 100))
        save(strcat(model.savePath, '/model'), 'model');

        % save Basis
        model.scmodel_Large.saveBasis();
        model.scmodel_Small.saveBasis();

        % TODO: DEPRECATED
        % save Weights
        % model.rlmodel.saveWeights();
    end
end
elapsedTime = toc;

% Total simulation time
model.simulatedTime = elapsedTime / 60;
sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
        elapsedTime / 3600, elapsedTime / 60, elapsedTime, trainTime / elapsedTime)

% Backup scripts & plot results
if (plotNsave(1) == 1)
    save(strcat(model.savePath, '/model'), 'model'); % storing simulated time

    copyfile('OESMuscles.m', model.savePath);
    copyfile('config.m', model.savePath);
    copyfile('Model.m', model.savePath);

    if (model.rlmodel.continuous == 1)
        copyfile('ReinforcementLearningCont.m', model.savePath);
    else
        copyfile('ReinforcementLearning.m', model.savePath);
    end

    switch model.rlmodel.rlFlavour(1)
        case 0
            %% Chong's implementation
            copyfile('CCriticG.m', model.savePath);
        case 1
            %% CRG
            copyfile('CRGCritic.m', model.savePath);
        case 2
            %% CACLA
            copyfile('CACLACritic.m', model.savePath);
    end

    switch model.rlmodel.rlFlavour(2)
        case 0
            %% Chong's implementation
            copyfile('CActorG.m', model.savePath);
        case 1
            %% CRG
            copyfile('CRGActor.m', model.savePath);
        case 2
            %% CACLA linear
            copyfile('CACLAActorLin.m', model.savePath);
        case 3
            %% CACLAVar linear
            copyfile('CACLAVarActorLin.m', model.savePath);
        case 4
            %% CACLA
            copyfile('CACLAActor.m', model.savePath);
        case 5
            %% CACLAVar
            copyfile('CACLAVarActor.m', model.savePath);
    end
    model.allPlotSave();
end

%%% Testing procedure
if (testIt)
    % testModel(model, randomizationSeed, objRange, vergRange, repeat, randStimuli, randObjRange, plotIt, saveTestResults)
    testModel(model, randomizationSeed, [0.5, 1, 1.5, 2], [-3 : 0.2 : 3], [50, 50], 0, 1, plotNsave(2), 1);
end

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

%% Generation of random vergence angles according to truncated Laplace distribution
function l = truncLaplacian(diversity, range)
    % see wikipedia for the generation of random numbers according to the
    % LaPlace distribution via the inversion method
    r = rand;

    switch r < 0.5
        case 1
            l = 1 / diversity * log(2 * r);
        case 0
            l = -1 / diversity * log(2 * (1 - r));
    end

    if (abs(l) > range)
        l = 0;
    end
end

%this function generates anaglyphs of the large and small scale fovea and
%one of the two unpreprocessed gray scale images
% TODO: adjust the sizes of the montage view
function generateAnaglyphs(leftGray, rightGray, dsRatioL, dsRatioS, foveaL, foveaS)
    anaglyph = imfuse(leftGray, rightGray, 'falsecolor');
    imwrite(anaglyph, 'anaglyph.png');

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
    imwrite(imresize(anaglyphL, 20), 'anaglyphLargeScale.png');
    largeScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(largeScaleView, 20), 'LargeScaleMontage.png');

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
    imwrite(imresize(anaglyphS, 16), 'anaglyphSmallScale.png');
    smallScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(smallScaleView, 8), 'smallScaleMontage.png');
end

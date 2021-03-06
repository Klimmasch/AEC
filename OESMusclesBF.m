%%% Script for running application with pre-learned basis functions
%@param trainTime           training time in number of iterations
%@param randomizationSeed   randomization seed
%
% learnedFile:              file with policy and sparse coding to test if not empty
% textureFile:              texture settings files
% sparseCodingType:         type of sparse coding approach
%
% Attention: Function is DEPICATED
%%%
function OESMusclesBF(trainTime, randomizationSeed, fileDescription)

rng(randomizationSeed);
learnedFile = '';
textureFile = 'Textures_vanHaterenTrain.mat';
sparseCodingType = 'nonhomeo';

% Plotting and saving flag
% Whether figures should be generated, saved and plotted
% additionally some relevant scripts are backed up
% plotNsave:    0 = don't do it
%               1 = do it
plotNsave = uint8(1);

% Testing flag
% Whether the testing procedure shall be executed in the end
% testIt:   0 = don't do it
%           1 = do it
testIt = uint8(1);

% Save model and conduct testing every saveInterval training iterations (+1)
saveInterval = 1000;
if (trainTime < saveInterval)
    saveInterval = trainTime;
end

% Instantiate and initiate model and test_data objects
model = config(learnedFile, textureFile, trainTime, sparseCodingType);
modelBF = load('/home/klimmasch/projects/results/model_12-Apr-2016_17:41:14_100000_nonhomeo_1_Bestdiscrete_highRes_cluster/model.mat');
model.scModel_Small.Basis = modelBF.model.scModel_Small.Basis;
model.scModel_Large.Basis = modelBF.model.scModel_Large.Basis;

if (trainTime <= model.interval)
    sprintf('trainTime[%d] must be > model.interval[%d]', trainTime, model.interval)
    return;
elseif (~exist(fullfile(cd, 'checkEnvironment'), 'file'))
    sprintf('Rendering binary \"checkEnvironment\" not present in current dir\n\"%s\"', cd)
    return;
end

% File management
modelName = sprintf('model_%s_%i_%i_%i_%s_%i_%s', ...
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
copyfile('config.m', model.savePath);
copyfile('OESMusclesBF.m', model.savePath);
copyfile('ReinforcementLearningCont.m', model.savePath);
copyfile('CRGActor.m', model.savePath);
copyfile('CRGCritic.m', model.savePath);
copyfile('Model.m', model.savePath);

model.notes = [model.notes fileDescription]; %just and idea to store some more information

% Image process variables
patchSize = 8;

dsRatioL = model.scModel_Large.Dsratio; %downsampling ratio (Large scale) | original 8
dsRatioS = model.scModel_Small.Dsratio; %downsampling ratio (Small scale) | original 2

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

% Textures
texture = load(sprintf('config/%s', textureFile));
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

%%% New renderer
simulator = OpenEyeSim('create');
simulator.initRenderer();
% simulator.reinitRenderer();

imgRawLeft = uint8(zeros(240, 320, 3));
imgRawRight = uint8(zeros(240, 320, 3));

function [imLeft, imRight] = refreshImages(texture, vergAngle, objDist)
    simulator.add_texture(1, texture);
    simulator.set_params(1, vergAngle, objDist); %2-angle 3-distance

    result = simulator.generate_left;
    result2 = simulator.generate_right;

    imLeft=uint8(zeros(240, 320, 3));
    k=1;l=1;
    for i = 1:3:length(result)
        imLeft(k,l,1) = result(i);
        imLeft(k,l,2) = result(i+1);
        imLeft(k,l,3) = result(i+2);

        l=l+1;
        if (l>320)
            l=1;
            k=k+1;
        end
    end
%     imLeft = COLOR;     %320x240 image

    imRight=uint8(zeros(240, 320, 3));
    k=1;l=1;
    for i = 1:3:length(result2)
        imRight(k,l,1) = result2(i);
        imRight(k,l,2) = result2(i+1);
        imRight(k,l,3) = result2(i+2);

        l=l+1;
        if (l>320)
            l=1;
            k=k+1;
        end
    end
%     imRight = COLOR2;     %320x240 image
end


%%% Main execution loop
t = 0;
% rewardFunction_prev = -5;
tic; % start time count
for iter1 = 1 : (model.trainTime / model.interval)

    % pick random texture every #interval times
    currentTexture = texture{(randi(nTextures, 1))};

    % random depth
    objDist = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);

    % reset muscle activities to random values
    % initialization for muscle in between borders of desired actvity
    % i.e. min and max stimulus distance
    command = [0, 0];
    command(2) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1,1); %only for one muscle
    % command(1) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1,1); %two muscles
    % command(2) = command(1);
    % command(2) = 0.1 * rand(1, 1); % random policy

    angleNew = getAngle(command) * 2;

    % Training of basis functions:
    % b = 1;      % diversity or variance of laplacian dist.
    % range = 5;  % maximum vergence is 5 degree
    % angleDes = 2 * atand(model.baseline / (2 * objDist));
    % angleNew = angleDes + truncLaplacian(b, range);

%     [status, res] = system(sprintf('./checkEnvironment %s %d %d %s/left.png %s/right.png', ...
%                                    currentTexture, objDist, angleNew, model.savePath, model.savePath));
%
%     % abort execution if error occured
%     if (status)
%         sprintf('Error in checkEnvironment:\n%s', res)
%         return;
%     end

%     [imgRawLeft, imgRawRight] = refreshImages(currentTexture, angleNew, objDist);

    for iter2 = 1 : model.interval
        t = t + 1;

        [imgRawLeft, imgRawRight] = refreshImages(currentTexture, -angleNew/2, objDist);
        % read input images and convert to gray scale
%         imgRawLeft = imread([model.savePath '/left.png']);
%         imgRawRight = imread([model.savePath '/right.png']);
        imgGrayLeft = .2989 * imgRawLeft(:,:,1) + .5870 * imgRawLeft(:,:,2) + .1140 * imgRawLeft(:,:,3);
        imgGrayRight = .2989 * imgRawRight(:,:,1) + .5870 * imgRawRight(:,:,2) + .1140 * imgRawRight(:,:,3);

        % Generate & save the anaglyph picture
        % anaglyph = stereoAnaglyph(imgGrayLeft, imgGrayRight); % only for matlab 2015 or newer
        % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [model.savePath '/anaglyph.png']); %this one works for all tested matlab
        %more advanced functions that generated the anaglyphs of the foveal views
%         generateAnaglyphs(imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS, model.savePath);

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
        if (model.rlModel.continuous == 1)
            feature = [feature; command(2) * model.lambdaMuscleFB];
        end

        %%% Calculate metabolic costs
        metCost = getMetCost(command) * 2;

        %%% Calculate reward function
        %% Standard reward
        rewardFunction = model.lambdaRec * reward - model.lambdaMet * metCost;
        % rewardFunction = (model.lambdaMet * reward) + ((1 - model.lambdaMet) * - metCost);

        %% Delta reward
        % rewardFunctionReal = model.lambdaRec * reward - model.lambdaMet * metCost;
        % rewardFunction = rewardFunctionReal - rewardFunction_prev;

        % counter balance 0 movement by small negative bias
        % rewardFunction = rewardFunctionReal - rewardFunction_prev - 1e-5;
        % rewardFunction_prev = rewardFunctionReal;


        %%% Weight L1 regularization
        % rewardFunction = model.lambdaRec * reward ...
        %                  - model.lambdaMet * metCost ...
        %                  - model.lambdaV * (sum(sum(abs(model.rlModel.CCritic.v_ji)))) ...
        %                  - model.lambdaP1 * (sum(sum(abs(model.rlModel.CActor.wp_ji)))) ...
        %                  - model.lambdaP2 * (sum(sum(abs(model.rlModel.CActor.wp_kj))));

        %%% Weight L2 regularization
        % rewardFunction = model.lambdaRec * reward ...
        %                  - model.lambdaMet * metCost ...
        %                  - model.lambdaV * (sum(sum(model.rlModel.CCritic.v_ji .^ 2))) ...
        %                  - model.lambdaP1 * (sum(sum(model.rlModel.CActor.wp_ji .^ 2))) ...
        %                  - model.lambdaP2 * (sum(sum(model.rlModel.CActor.wp_kj .^ 2)));

        %%% Learning
        % Sparse coding models
        %model.scModel_Large.stepTrain(currentView{1});
        %model.scModel_Small.stepTrain(currentView{2});

        % RL model
        % Variance decay, i.e. reduction of actor's output perturbation
        if ((model.rlModel.continuous == 1) && (model.rlModel.CActor.varDec > 0))
            model.rlModel.CActor.variance = model.rlModel.CActor.varianceRange(1) * 2 ^ (-t / model.rlModel.CActor.varDec);
        end

        relativeCommand = model.rlModel.stepTrain(feature, rewardFunction, (iter2 > 1));

        % add the change in muscle Activities to current ones
        % command = command + relativeCommand';     %two muscels
        command(2) = command(2) + relativeCommand;  %one muscel
        command = checkCmd(command);                %restrain motor commands to [0,1]
        % angleNew = getAngle(command) * 2;           %resulting angle is used for both eyes

        if (model.rlModel.continuous == 1)
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
%         [status, res] = system(sprintf('./checkEnvironment %s %d %d %s/left.png %s/right.png', ...
%                                    currentTexture, objDist, angleNew, model.savePath, model.savePath));
%
%         % abort execution if error occured
%         if (status)
%             sprintf('Error in checkEnvironment:\n%s', res)
%             return;
%         end

%         [imgRawLeft, imgRawRight] = refreshImages(currentTexture, angleNew, objDist);

        %%%%%%%%%%%%%%%% TRACK ALL PARAMETERS %%%%%%%%%%%%%%%%%%

        % compute desired vergence command, disparity and vergence error
        fixDepth = (model.baseline / 2) / tand(angleNew / 2);   %fixation depth [m]
        angleDes = 2 * atand(model.baseline / (2 * objDist));   %desired vergence [deg]
        anglerr = angleDes - angleNew;                          %vergence error [deg]
        disparity = 2 * model.focalLength * tand(anglerr / 2);            %current disp [px]

        % save state
        % model.Z(t) = objDist;
        % model.fixZ(t) = fixDepth;
        % model.disp_hist(t) = disparity;

        model.vergerr_hist(t) = anglerr;
        model.recerr_hist(t, :) = [errorLarge; errorSmall];
        % model.verge_actual(t) = angleNew;
        % model.verge_desired(t) = angleDes;
        model.relCmd_hist(t) = relativeCommand;
        model.cmd_hist(t, :) = command;
        % model.reward_hist(t) = rewardFunction;
        % model.feature_hist(t, :) = feature;
        model.metCost_hist(t) = metCost;
        if (model.rlModel.continuous == 1)
            model.td_hist(t) = model.rlModel.CCritic.delta;
            % model.g_hist(t) = model.rlModel.CActor.params(7);
            model.weight_hist(t, 1) = model.rlModel.CCritic.params(1);
            model.weight_hist(t, 2) = model.rlModel.CActor.params(1);
            if (model.rlModel.rlFlavour(2) >= 4)
                model.weight_hist(t, 3) = model.rlModel.CActor.params(2);
                if ((model.rlModel.rlFlavour(2) == 5) || (model.rlModel.rlFlavour(2) == 7))
                    model.weight_hist(t, 4) = model.rlModel.CActor.params(3);
                end
            end
            model.variance_hist(t) = model.rlModel.CActor.variance;
        end
        model.trainedUntil = t;
    end

    sprintf('Training Iteration = %d\nAbs Command =\t[%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f]\nRel Command = \t[%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f]\nVer Error =\t[%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f]', ...
            t, model.cmd_hist(t - model.interval + 1 : t, 2), model.relCmd_hist(t - model.interval + 1 : t), model.vergerr_hist(t - model.interval + 1 : t))

    % Display per cent completed of training and save model
    if (~mod(t, saveInterval))
        sprintf('%g%% is finished', (t / model.trainTime * 100))
        save(strcat(model.savePath, '/model'), 'model');

        % save Basis
        model.scModel_Large.saveBasis();
        model.scModel_Small.saveBasis();

        % save Weights
        % model.rlModel.saveWeights();
    end
end
elapsedTime = toc;

% Total simulation time
model.simulatedTime = elapsedTime / 60;
sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
        elapsedTime / 3600, elapsedTime / 60, elapsedTime, trainTime / elapsedTime)

% plot results
if (plotNsave(1) == 1)
%     if useLearnedFile(2)
%         model.trainTime = model.trainTime + model.trainedUntil;
%     end
    save(strcat(model.savePath, '/model'), 'model'); % storing simulated time

    if (model.rlModel.continuous == 1)
        copyfile('ReinforcementLearningCont.m', model.savePath);
    else
        copyfile('ReinforcementLearning.m', model.savePath);
    end

    switch model.rlModel.rlFlavour(1)
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

    switch model.rlModel.rlFlavour(2)
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
        case 6
            %% CACLA2
            copyfile('CACLAActor2.m', model.savePath);
        case 7
            %% CACLAVar2
            copyfile('CACLAVarActor2.m', model.savePath);
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
            l = 1/diversity*log(2*r);
        case 0
            l = -1/diversity*log(2*(1-r));
    end

    if abs(l) > range
        l = 0;
    end
end

%this function generates anaglyphs of the large and small scale fovea and
%one of the two unpreprocessed gray scale images
% TODO: adjust the sizes of the montage view
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


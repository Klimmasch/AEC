%%% Script for running application with pre-learned basis functions
%@param trainTime           training time in number of iterations
%@param randomizationSeed   randomization seed
%
% learnedFile:              file with policy and sparse coding to test if not empty
% textureFile:              texture settings files
% sparseCodingType:         type of sparse coding approach
%%%
function OESMusclesBF(trainTime, randomizationSeed, fileDescription)

rng(randomizationSeed);
learnedFile = '';
textureFile = 'Textures_New.mat';
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
% modelBF = load('/tmp/model_11-Mar-2016 20:38:05_100000_110_111_nhomeo_1_trainBF_laplacian_b=1/model.mat');
% model.scmodel_Small.Basis = modelBF.model.scmodel_Small.Basis;
% model.scmodel_Large.Basis = modelBF.model.scmodel_Large.Basis;

if (trainTime <= model.interval)
    sprintf('trainTime[%d] must be > model.interval[%d]', trainTime, model.interval)
    return;
elseif (~exist(fullfile(cd, 'checkEnvironment'), 'file'))
    sprintf('Rendering binary \"checkEnvironment\" not present in current dir\n\"%s\"', cd)
    return;
end

% File management
modelName = sprintf('model_%s_%i_%i_%i_%s_%i_%s', ...
                    datestr(now), ...
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
rewardFunction_prev = -5;
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
    b = 1;      % diversity or variance of laplacian dist.
    range = 5;  % maximum vergence is 5 degree
    angleDes = 2 * atand(model.baseline / (2 * objDist));
    angleNew = angleDes + truncLaplacian(b, range);

    [status, res] = system(sprintf('./checkEnvironment %s %s %d %d left.png right.png %d', ...
                                   currentTexture, currentTexture, objDist, objDist, angleNew));

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
        % anaglyph = stereoAnaglyph(imgGrayLeft, imgGrayRight); % only for
        % matlab 2015 or newer
        anaglyph = imfuse(imgGrayLeft, imgGrayRight, 'falsecolor');
        imwrite(anaglyph, 'anaglyph.png');
        
        % Image patch generation: left{small scale, large scale}, right{small scale, large scale}
        [patchesLeftSmall] = preprocessImage(imgGrayLeft, foveaS, dsRatioS, patchSize, columnIndS);
        [patchesLeftLarge] = preprocessImage(imgGrayLeft, foveaL, dsRatioL, patchSize, columnIndL);
        [patchesRightSmall] = preprocessImage(imgGrayRight, foveaS, dsRatioS, patchSize, columnIndS);
        [patchesRightLarge] = preprocessImage(imgGrayRight, foveaL, dsRatioL, patchSize, columnIndL);

        % Image patches matrix (input to model)
        currentView = {[patchesLeftLarge; patchesRightLarge] [patchesLeftSmall; patchesRightSmall]};

        % Generate input feature vector from current images
        [feature, reward, errorTotal, errorLarge, errorSmall] = model.generateFR(currentView);

        %%% Feedback
        % Absolute command feedback # concatination
        feature = [feature; command(2) * model.lambdaMuscleFB];

        %%% Calculate metabolic costs
        metCost = getMetCost(command) * 2;

        %%% Calculate reward function
        % delta reward
        rewardFunctionReal = model.lambdaRec * reward - model.lambdaMet * metCost;
        % rewardFunction = rewardFunctionReal - rewardFunction_prev;

        % counter balance 0 movement by small negative bias
        rewardFunction = rewardFunctionReal - rewardFunction_prev - 1e-5;
        rewardFunction_prev = rewardFunctionReal;

        % standard reward
        % rewardFunction = model.lambdaRec * reward - model.lambdaMet * metCost;
        % rewardFunction = (model.lambdaMet * reward) + ((1 - model.lambdaMet) * - metCost);

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
        %model.scmodel_Large.stepTrain(currentView{1});
        %model.scmodel_Small.stepTrain(currentView{2});

        % RL model
        % decay of actor's output perturbation
        % variance(t = [1, 100k]) ~= [0.001, 1e-5]
        % model.rlmodel.CActor.variance = 0.001 * 2 ^ (-t / 15000);
        % variance(t = [1, 100k]) ~= [0.001, 1e-4]
        model.rlmodel.CActor.variance = 0.001 * 2 ^ (-t / 30200);
        % variance(t = [1, 100k]) ~= [0.01, 1e-4]
        % model.rlmodel.CActor.variance = 0.01 * 2 ^ (-t / 15100);

        relativeCommand = model.rlmodel.stepTrain(feature, rewardFunction, (iter2 > 1));

        % add the change in muscle Activities to current ones
        % command = command + relativeCommand';     %two muscels
        command(2) = command(2) + relativeCommand;  %one muscel
        command = checkCmd(command);                %restrain motor commands to [0,1]
        angleNew = getAngle(command) * 2;           %resulting angle is used for both eyes

        % in case you want to train basisfunctions tuned to a specific
        % disparity:
%         angleDes = 2 * atand(model.baseline / (2 * objDist)); %desired vergence [deg]
%         angleNew = angleDes + truncLaplacian(b,range);
%         testLP = [];
%         for test = 1:10000
%             testLP = [testLP truncLaplacian(b,range)];
%         end
%         figure;
%         hold on;
%         histogram(testLP);
%         hold off;

        % generate new view (two pictures) with new vergence angle
        [status, res] = system(sprintf('./checkEnvironment %s %s %d %d left.png right.png %d', ...
                               currentTexture, currentTexture, objDist, objDist, angleNew));

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
        disparity = 2 * model.f * tand(anglerr / 2);            %current disp [px]

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
        model.feature_hist(t, :) = feature;
        model.metCost_hist(t) = metCost;
        model.td_hist(t) = model.rlmodel.CCritic.delta;
        % model.g_hist(t) = model.rlmodel.CActor.params(7);
        model.l12_weights(t, 1) = model.rlmodel.CCritic.params(1);
        model.l12_weights(t, 2) = model.rlmodel.CCritic.params(2);
        model.l12_weights(t, 3) = model.rlmodel.CActor.params(1);
        model.l12_weights(t, 4) = model.rlmodel.CActor.params(2);
        % model.l12_weights(t, 5) = model.rlmodel.CActor.params(3);
        % model.l12_weights(t, 6) = model.rlmodel.CActor.params(4);
        % model.l12_weights(t, 7) = model.rlmodel.CActor.params(5);
        % model.l12_weights(t, 8) = model.rlmodel.CActor.params(6);
        % plot(model.td_hist);
        % if (model.feature_hist(end,end) > sum(model.feature_hist(t,1:end-1),2))
        %     figure
        %     hold on;
        %     plot(sum(model.feature_hist(:,1:end-1),2),'r');
        %     plot(model.feature_hist(:,end),'b');
        %     title(sprintf('%g', model.feature_hist(end,end)));
        % end
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

        % save Weights
        % model.rlmodel.saveWeights();
    end
end
elapsedTime = toc;

% Total simulation time
model.simulated_time = elapsedTime / 60;
sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
        elapsedTime / 3600, elapsedTime / 60, elapsedTime, trainTime / elapsedTime)

% Plot results
if (plotNsave)
    % model.errPlot();
    model.allPlotSave();
    copyfile('config.m', model.savePath);
    copyfile('OESMuscles.m', model.savePath);
    copyfile('OESMusclesBF.m', model.savePath);
    copyfile('ReinforcementLearningCont.m', model.savePath);
    copyfile('CRGActor.m', model.savePath);
    copyfile('CRGCritic.m', model.savePath);
    copyfile('Model.m', model.savePath);
end

%%% Testing procedure
if (testIt)
    % TestTrial(model, randomizationSeed, fileDescription);
    model.deltaMFplotGenDist([0.5, 1, 2], [-5:5], 5);
    model.recErrPlotGenDist([0.5, 1, 2], [-5:5], 5);
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
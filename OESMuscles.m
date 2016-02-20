%%% Main script for launching application
%@param trainTime           training time in number of iterations
%@param randomizationSeed   randomization seed
%
% learnedFile:              file with policy and sparse coding to test if not empty
% textureFile:              texture settings files
% sparseCodingType:         type of sparse coding approach
%%%

function OESMuscles(trainTime, randomizationSeed, description)

rng(randomizationSeed);
learnedFile = '';
textureFile = 'Textures_New.mat';
sparseCodingType = 'nonhomeo';

% Plotting flag
% Whether figures should be generated, saved and plotted
% plotIt:   0 = no plot
%           1 = plot
plotIt = uint8(1);

% Learning flag
% learning: 0 = no model is being learned
%           1 = model is being learned
learning = uint8(isempty(learnedFile));

% Predefine testing phases
testingWindow = uint8(0);
testingCounter = 0;

% Save model and conduct testing every saveInterval training iterations (+1)
saveInterval = 10;

% Testing flag
% testing or training, we do this every 10% of learning
% trialPhase: 0 = training
%             1 = testing
trialPhase = uint8(0);

% Instantiate and initiate model and test_data objects
model = config(learnedFile, textureFile, trainTime, sparseCodingType);
modelData = ModelTestdata((trainTime / saveInterval) * testingWindow + 1);

if (trainTime <= model.interval)
    sprintf('trainTime[%d] must be > model.interval[%d]', trainTime, model.interval)
    return;
elseif (~exist(fullfile(cd, 'checkEnvironment'), 'file'))
    sprintf('Rendering binary \"checkEnvironment\" not present in current dir\n\"%s\"', cd)
    return;
end

% File management
% description = 'contin_';
savePath = sprintf('model_%s_%i_%i_%i_%s_%i_%s', ...
                    datestr(now), ...
                    trainTime, ...
                    sparseCodingType, ...
                    randomizationSeed, ...
                    description);
% mkdir('results', savePath);
savePath = strcat('~/projects/RESULTS/', savePath);

% Control Parameters
% vergeMax = 16;
% angleNew = 0; %init new vergence command

% Iteration counter
t = 0;
testIter = 0;

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

% Camera parameters
% offset = 0;       %vertical offset between left and right (0 in the Simulator!!!)
f = 257.34;         %focal length [px]
baseline = 0.056;   %interocular distance (baseline)

% Textures
texturePath = sprintf('config/%s', textureFile);
texture = load(texturePath);
texture = texture.texture;
nTextures = length(texture);
currentTexture = texture{1}; %choose first texture as initial

% Object distance to eyes [m]
objDistMin = 0.5;
objDistMax = 2;
objDist = objDistMax; %object init position

degrees = load('Degrees.mat');    %loads tabular for resulting degrees as 'results_deg'
metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

%%% Helper function that maps muscle activities to resulting angle
function [angle] = getAngle(command)
    cmd = (command * 10) + 1; % scale commands to table entries
    angle = interp2(degrees.results_deg, cmd(1), cmd(2)); % interpolate in tabular
end

%%% Helper function that maps muscle activities to resulting angle
function [angle] = getMetCost(command)
    cmd = (command * 10) + 1; % scale commands to table entries
    angle = interp2(metCosts.results, cmd(1), cmd(2)); % interpolate in tabular
end

%%% Main execution loop
tic %start time count
while (true)
    if (trialPhase == 0)
        t = t + 1;
        if (t > model.trainTime)
            trialPhase = uint8(1);
        end

        % pick random texture every #interval times
        if (~mod(t - 1, model.interval))
            currentTexture = texture{(randi(nTextures, 1))};
            % random depth
            objDist = objDistMin + (objDistMax - objDistMin) * rand(1, 1);
            % reset muscle activities to random values
			command = 0.0075 + (0.06 - 0.0075) * rand(1,1); % initialization for muscle in between borders of desired actvity
            command = [0, command]; %try to learn only one dimension
            angleNew = getAngle(command) * 2;
            [status, res] = system(sprintf('./checkEnvironment %s %s %d %d left.png right.png %d', ...
                                           currentTexture, currentTexture, objDist, objDist, angleNew));

            % Abort execution if error occured
            if (status)
                sprintf('Error in checkEnvironment:\n%s', res)
                return;
            end
        end
    elseif (trialPhase == 1)
        testingCounter = testingCounter + 1;
        testIter = testIter + 1;
        if (testingCounter > testingWindow)
            if (testingWindow > 0)
                save(strcat(savePath, '/modelData'), 'modelData');
            end

            % Exit condition
            if (t == model.trainTime)
                break;
            end

            trialPhase = uint8(0);
            testingCounter = 0;
            t = t + 1;
        end

        % pick random texture every #interval times
        if (~mod(testingCounter - 1, model.interval))
            currentTexture = texture{(randi(nTextures, 1))};
            % random depth
            objDist = objDistMin + (objDistMax - objDistMin) * rand(1, 1);
            % reset muscle activities to random values
% 			command = rand(2,1);
            command = 0.0075 + (0.06 - 0.0075) * rand(1,1); % initialization for muscle in between borders of desired actvity
            command = [0, command]; %try to learn only one dimension
            angleNew = getAngle(command) * 2;
            [status, res] = system(sprintf('./checkEnvironment %s %s %d %d left.png right.png %d', ...
                                           currentTexture, currentTexture, objDist, objDist, angleNew));

            % Abort execution if error occured
            if (status)
                sprintf('Error in checkEnvironment:\n%s', res)
                return;
            end
        end
    end

    % Read input images
    imgRawLeft = imread('left.png');
    imgRawRight = imread('right.png');

    % Image patch generation: left{small scale, large scale}, right{small scale, large scale}
    [patchesLeftSmall] = preprocessImage(imgRawLeft, foveaS, dsRatioS, patchSize, columnIndS);
    [patchesLeftLarge] = preprocessImage(imgRawLeft, foveaL, dsRatioL, patchSize, columnIndL);
    [patchesRightSmall] = preprocessImage(imgRawRight, foveaS, dsRatioS, patchSize, columnIndS);
    [patchesRightLarge] = preprocessImage(imgRawRight, foveaL, dsRatioL, patchSize, columnIndL);

    % image patches matrix (input to model)
    currentView = {[patchesLeftLarge; patchesRightLarge] [patchesLeftSmall; patchesRightSmall]};

    if (trialPhase == 0)
        %generate input feature vector from current images
        [feature, reward, errorTotal, errorLarge, errorSmall] = model.generateFR(currentView);
        feature = [feature; command(2)*0.01]; % incorporationg the current muscle activity into feature vector
                                              % and scaling it to the value
                                              % range of BF activations
        
        %calculate metabolic costs and reward function
        metCost = 2 * getMetCost(command);
        rewardFunction = (model.lambdaMet * reward) + ((1-model.lambdaMet) * - metCost);
        
        if (learning)
            %train 2 sparse coding models
            model.scmodel_Large.stepTrain(currentView{1});
            model.scmodel_Small.stepTrain(currentView{2});
            [relativeCommand, paramsC, paramsA] = model.rlmodel.stepTrain(feature, rewardFunction, mod(t - 1, model.interval));
            %### why calculate relative vergance command after first error measure?
        end
        
        relativeCommand %print for debugging
        
		% add the change in muscle Activities to current ones
		command(2) = command(2) + relativeCommand;
		command = checkCmd(command); %restrain motor commands to [0,1]
		angleNew = getAngle(command) * 2; %resulting angle is used for both eyes
        
        %angleNew = max(0.01, command + angleNew); %command is relative angle - constrain to positive vergence
        %safety on control command - RESET TO a new vergence angle
        %if (angleNew > vergeMax)
        %    angleNew = randi(vergeMax, 1);     %relax the eyes
        %end
        
        %generate new view (two pictures) with new vergence angle
        [status, res] = system(sprintf('./checkEnvironment %s %s %d %d left.png right.png %d', ...
                               currentTexture, currentTexture, objDist, objDist, angleNew));

        % Abort execution if error occured
        if (status)
            sprintf('Error in checkEnvironment:\n%s', res)
            return;
        end
        %%%%%%%%%%%%%%%% TRACK ALL PARAMETERS %%%%%%%%%%%%%%%%%%
        
        %Compute desired vergence command, disparity and vergence error
        fixDepth = (0.5 * baseline) / tand(angleNew / 2);
        angledes = 2 * atan(baseline / (2 * objDist)); %desired vergence [rad]
        anglerr = angledes * 180 / pi - angleNew; %vergence error [deg]
        disparity = 2 * f * tan((angledes - angleNew * pi / 180) / 2); %current disp [px]
        
        %save them
        model.Z(t) = objDist;
        model.fixZ(t) = fixDepth;
        model.disp_hist(t) = disparity;
        model.vergerr_hist(t) = anglerr;
        model.recerr_hist(t, :) = [errorLarge; errorSmall];
        model.verge_actual(t) = angleNew; 
        model.relCmd_hist(t,:) = relativeCommand;
        model.cmd_hist(t,:) = command;
        model.reward_hist(t) = rewardFunction;
        model.feature_hist(t,:) = feature;
        model.medCost_hist(t) = metCost;
        model.td_hist(t) = paramsC(2);
        model.g_hist(t) = paramsA(7);
        model.AC_norm_weights(t,1) = paramsC(1);%norm(v_ji)
        model.AC_norm_weights(t,2) = paramsA(1);%norm(wp_ji)
        model.AC_norm_weights(t,3) = paramsA(2);%norm(dwp)
        model.AC_norm_weights(t,4) = paramsA(3);%norm(wp_kj)
        model.AC_norm_weights(t,5) = paramsA(4);%norm(dvp)
        model.AC_norm_weights(t,6) = paramsA(5);%norm(wn_ij)
        model.AC_norm_weights(t,7) = paramsA(6);%psi' * this.wn_ji
        
        sprintf('Training Iteration = %d\nCommand = %.3g,%.3g\tCurrent Vergence = %.3g\tVergence Error = %.3g\nRec Error = %.3g', ...
                t, command(1), command(2), angleNew, anglerr, errorTotal)
            
        % Display per cent completed of training and save model
        if (~mod(t, model.trainTime / saveInterval) && (learning))
            sprintf('%g%% is finished', (t / model.trainTime * 100))
            save(strcat(savePath, '/model'), 'model');

            %save Basis
            model.scmodel_Large.saveBasis;
            model.scmodel_Small.saveBasis;

            %save Weights
            model.rlmodel.saveWeights; %save policy and value net weights

            trialPhase = uint8(1);
        end

    elseif (trialPhase == 1)
        %Compute Vergence Command
        fixDepth = (0.5 * baseline) / tand(angleNew / 2);
        %save fixaton and object depth
        modelData.Z(testIter) = objDist;
        modelData.fixZ(testIter) = fixDepth;

        %compute current disparity
        angledes = 2 * atan(baseline / (2 * objDist)); %desired vergence [rad]
        disparity = 2 * f * tan((angledes - angleNew * pi / 180) / 2); %current disp [px]
        modelData.disp_hist(testIter) = disparity;

        %compute current vergence error
        anglerr = angledes * 180 / pi - angleNew; %vergence error [deg]
        %save verge error
        modelData.vergerr_hist(testIter) = anglerr;

        %generate input feature vector from current images
        [feature, ~, errorTotal, errorLarge, errorSmall] = model.generateFR(currentView);
        modelData.recerr_hist(testIter, :) = [errorLarge; errorSmall];
        modelData.verge_actual(testIter) = angleNew; %current vergence angle
        feature = [feature; command(2)*0.01];
        if (learning)
            relativeCommand = model.rlmodel.softmaxAct(feature);
            %### why calculate relative vergance command after first error measure?
        end
		command(2) = command(2) + relativeCommand;
		checkCmd(command)
		angleNew = getAngle(command) * 2;
        sprintf('Testing Iteration = %d\nCommand = %.3g,%.3g\tCurrent Vergence = %.3g\tVergence Error = %.3g\nRec Error = %.3g', ...
                testingCounter, command(1), command(2), angleNew, anglerr, errorTotal)

        %angleNew = max(0.01, command + angleNew); %command is relative angle - constrain to positive vergence
        %safety on control command - RESET TO a new vergence angle
        %if (angleNew > vergeMax)
        %    angleNew = randi(vergeMax, 1);     %relax the eyes
        %end
        [status, res] = system(sprintf('./checkEnvironment %s %s %d %d left.png right.png %d', ...
                                       currentTexture, currentTexture, objDist, objDist, angleNew));

        % Abort execution if error occured
        if (status)
            sprintf('Error in checkEnvironment:\n%s', res)
            return;
        end
    end
end
elapsedTime = toc;

% Total simulation time
model.simulated_time = elapsedTime / 60;
sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]', elapsedTime / 3600, elapsedTime / 60, elapsedTime)
sprintf('Frequency = %.4f [iterations/sec]', trainTime / elapsedTime)

% Plot results
if (plotIt)
    % model.errPlot();
    model.errPlotSave(savePath);
end

% Save results data
save(strcat(savePath, '/model'), 'model');
save(strcat(savePath, '/modelData'), 'modelData');

end


%%% Helper function that keeps motor commands in [0, 1]
function [cmd] = checkCmd(cmd)
	i0 = cmd < 0;
	cmd(i0) = 0;
	i1 = cmd > 1;
	cmd(i1) = 1;
end

%%% Helper functions for image preprocessing
%% Patch generation
function [patches] = preprocessImage(img, fovea, downSampling, patchSize, columnIndicies)
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

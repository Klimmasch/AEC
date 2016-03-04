%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes the current model loaded in the workspace and
% simulates a test trial with it:
%       1.000 testing episodes
%       traking of vergence Error, reconstruction Error, muscle forces
%           and metabolic costs
% No learning occures during this trial and the results are saved in a
% specified folder (under ./results/) in the file modelTestData.mat.
function TestTrial(model, randomizationSeed, fileDescription, savePathMo)

% TODO: range of disparity ueberpruefen auf relCmd!!

function TestTrial(model, randomizationSeed, fileDescription)

    numberTrials = 100;
    modelTest = ModelTestData(numberTrials * model.interval, model.interval);
    folder = strcat(savePathMo, './testResults/');
    savePath = sprintf('TestedModel_%s_%s', datestr(now), fileDescription);
    mkdir(folder, savePath);
    savePath = strcat(folder, savePath);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% predefining variables %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f = 257.34;         %focal length [px]
    baseline = 0.056;   %interocular distance (baseline)

    objDistMin = 0.5;   %minimal object distance
    objDistMax = 2;     %maximal object distance

    muscleInitMin = 0.00807;    %minimal initial muscle innervation
    muscleInitMax = 0.07186;    %maximal --"--

    degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
    metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

    % command = [0, 0];   %initialization of muscle commands

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

    % preparing Textures
%     texturePath = sprintf('config/%s', 'Textures_New.mat');
    texture = load('config/Textures_New.mat');
    texture = texture.texture;
    nTextures = length(texture);

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% starting the main loop %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = 1;
    rng(randomizationSeed);
    tic
    for iter1 = 1 : numberTrials

        % pick random texture every #interval times
        currentTexture = texture{(randi(nTextures, 1))};
        % random depth
        objDist = objDistMin + (objDistMax - objDistMin) * rand(1, 1);
        % reset muscle activities to random values
        command = [0, 0];
        command(2) = muscleInitMin + (muscleInitMax - muscleInitMin) * rand(1, 1); %only for one muscle
        angleNew = getAngle(command) * 2;

        %generate two new pictures
        [status, res] = system(sprintf('./checkEnvironment %s %s %d %d left.png right.png %d', ...
                                   currentTexture, currentTexture, objDist, objDist, angleNew));

        % Abort execution if error occured
        if (status)
            sprintf('Error in checkEnvironment:\n%s', res)
            return;
        end

        for iter2 = 1 : model.interval
            % Read input images and convert to gray scale
            imgRawLeft = imread('left.png');
            imgRawRight = imread('right.png');
            imgGrayLeft = .2989 * imgRawLeft(:,:,1) + .5870 * imgRawLeft(:,:,2) + .1140 * imgRawLeft(:,:,3);
            imgGrayRight = .2989 * imgRawRight(:,:,1) + .5870 * imgRawRight(:,:,2) + .1140 * imgRawRight(:,:,3);

            anaglyph = stereoAnaglyph(imgGrayLeft, imgGrayRight);
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
            % feature = [feature; 0; 0]; %just for testing purposes

            %%% Feedback
            % Absolute command feedback # concatination
            feature = [feature; command(2) * model.lambdaMuscleFB];
            % feature = [feature; command' * 0.01]; % just to make it how I trained it before ('ChongsParams')
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


            % generation of motor command without learning and noise
            % [relativeCommand, ~, ~] = model.rlmodel.stepTrain(feature, rewardFunction, 0);
            relativeCommand = model.rlmodel.softmaxAct(feature);

            % command = command + relativeCommand';     %two muscels
            command(2) = command(2) + relativeCommand;  %one muscel
            command = checkCmd(command);                %restrain motor commands to [0,1]
            angleNew = getAngle(command) * 2;           %resulting angle is used for both eyes

            % generate new view (two pictures) with new vergence angle
            [status, res] = system(sprintf('./checkEnvironment %s %s %d %d left.png right.png %d', ...
                                   currentTexture, currentTexture, objDist, objDist, angleNew));

            % Abort execution if error occured
            if (status)
                sprintf('Error in checkEnvironment:\n%s', res)
                return;
            end

            %%%%%%%%%%%%%%%% TRACK ALL PARAMETERS %%%%%%%%%%%%%%%%%%

            %Compute desired vergence command, disparity and vergence error
            fixDepth = (baseline / 2) / tand(angleNew / 2);
            angleDes = 2 * atand(baseline / (2 * objDist)); %desired vergence [deg]
            anglerr = angleDes - angleNew;                  %vergence error [deg]
            disparity = 2 * f * tand(anglerr / 2);          %current disp [px]

            %save them
            modelTest.Z(t) = objDist;
            modelTest.fixZ(t) = fixDepth;
            modelTest.disp_hist(t) = disparity;
            modelTest.vergerr_hist(t) = anglerr;
            modelTest.recerr_hist(t, :) = [errorLarge; errorSmall];
            modelTest.verge_actual(t) = angleNew;
            modelTest.verge_desired(t) = angleDes;
            modelTest.relCmd_hist(t, 2) = relativeCommand;          %one muscle
            % modelTest.relCmd_hist(t, :) = relativeCommand;          %two muscles
            modelTest.cmd_hist(t, :) = command;
            % modelTest.reward_hist(t) = rewardFunction;
            modelTest.metCost_hist(t) = metCost;

            t = t + 1;
        end
        sprintf('Testing Iteration = %d\nCommand = [%.3g,\t%.3g]\tCurrent Vergence = %.3g\nRec Error = %.3g\tVergence Error =\n[%.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g]', ...
            t, command(1), command(2), angleNew, errorTotal, modelTest.vergerr_hist(t - modelTest.interval : t - 1))
    end
    elapsedTime = toc;
    sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
            elapsedTime / 3600, elapsedTime / 60, elapsedTime, t - 1 / elapsedTime)

    % Save and plot results data
    save(strcat(savePath, '/modelTestData'), 'modelTest');
    modelTest.testPlotSave(savePath);
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

%This function goes through all object distances and averages the motor command 
% and reconstruction error for different vergence errors.
%TODO what about muscle innervations
function [avgCmd, avgRecErr] = averageForVergenceErrors(model, objectRange, vergRange)
    avgCmd = [];
    avgVergErr = [];
    focalLength = 257.34;        	%focal length [px]
    baseline = 0.056;   		%interocular distance (baseline)
    for objDist = objectRange(1):0.1:objectRange(2)
        angleDes = 2 * atand(baseline / (2 * objDist));
        for verg = vergRange(1):0.1:vergeRange(2)
            relCmd = model.rlmodel.actHard;
        end
    end
end

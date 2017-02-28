%%% Main script for launching experimental procedure
% @param trainTime           training time in number of iterations
% @param randomizationSeed   randomization seed
% @param clusterCall         whether this is a cluster job 0 (no) 1 (yes)
% @param inputParams         a list composed of 'identifier' followed by value, e.g {'alpha', 1, 'beta', 2, ...}
%                            whereas 'identifier' has to be present in configVar.m
% @param experimentName     description of approach used as file name
%%%

function OES2Muscles(trainTime, randomizationSeed, clusterCall, inputParams, experimentDirName, experimentName)

    % Initialize random number generator with given seed
    rng(randomizationSeed);

    % useLearnedFile(1):    0 = don't do it
    %                       1 = use previously learned policy specified in learnedFile
    % useLearnedFile(2):    0 = retrain with same/new parameters
    %                       1 = complete/continue training
    %
    % If useLearnedFile = [1, 1] and you want to continue training, trainTime must be the overall desired train time,
    % i.e. trained at first 100k iterations with useLearnedFile = [0, 0], then decided to continue training for 100k more
    % iterations, then the overall train time for the model is 200k and you set useLearnedFile = [1, 1] and execute with
    % OES2Muscles(200000, randomizationSeed, clusterCall, inputParams, experimentName)
    useLearnedFile = [0, 0];
    learnedFile = '';
    % useLearnedFile = [1, 1];
    % learnedFile = '/home/aecgroup/aecdata/Results/17-01-11_2000000iter_1_feature_Vector_Normalize/model.mat';

    %%% Stimulus declaration
    % textureFile = 'Textures_mcgillManMadeTrain(jpg).mat';     % McGill man made database
    % textureFile = 'Textures_mcgillFruitsAll.mat';             % McGill fruits database
    % textureFile = 'Textures_mcgillFoliageTrain(jpg).mat';     % McGill foliage database
    % textureFile = 'Textures_vanHaterenTrain.mat';             % vanHateren database
    % textureFile = 'Textures_celine.mat';                      % Celine's images

    % for the new renderer, all textures to be used during training and
    % testing have to be loaded into the buffer at the beginning
    % per convention, the testing images are given in the first entry!!
    % textureFiles = {'mcGillTest2.mat', 'mcGillTest1.mat'}; % test files containing less images
    textureFiles = {'Textures_mcgillManMade40.mat', 'Textures_mcgillManMade100.mat'};

    %%% Execute intermediate test procedure during training
    % testAt = [500000 : 500000 : trainTime];
    testAt = [250000 : 250000 : trainTime];

    %%% Testing flag
    % Whether the testing procedure shall be executed after training
    % testIt:   0 = don't do it
    %           1 = do it
    testIt = uint8(1);

    %%% Amount of test stimuli
    nStimTest = 40;

    % Save model every #saveInterval training iterations
    saveInterval = ceil(trainTime / 20);
    % saveInterval = trainTime;

    % Track the evolution of all basis functions of the respective sparse coders
    trackSCBasisHistory = uint8(0);

    %%% Sparse coding approach
    % sparseCodingType: 0 = non-homeostatic
    %                   1 = homeostatic
    sparseCodingType = uint8(0);
    sparseCodingTypeName = cellstr(['nonhomeo'; 'homeosSC']);

    %%% Plotting flag
    % Whether figures should be generated and saved
    % plotIt: [training, testing]
    %            0 = don't do it
    %            1 = do it
    plotIt = [uint8(1), uint8(1)];

    %%% Whether figures should be closed after generation
    % closeFigures: 0 = don't do it
    %               1 = do it
    closeFigures = uint8(1); % maybe necessary to use 0 when running headless

    % check whether given (cluster) job can/shall be continued
    % if (clusterCall == 1)
        if (isempty(experimentDirName))
            absoluteDir = '/home/aecgroup/aecdata/Results';
        else
            absoluteDir = strcat('/home/aecgroup/aecdata/Results/', experimentDirName);
        end

        fullDir = dir(sprintf('%s/*%s*', absoluteDir, experimentName));
        if length(fullDir) > 1
            fullDir = fullDir(1);   % in case there are more than one folder (different days, etc)
        end
        if (~isempty(fullDir))
            learnedFile = strcat(absoluteDir, '/', fullDir.name, '/model.mat');

            if (exist(learnedFile, 'file') == 2) % indicates being a file (7 == directory)
                useLearnedFile = [1, 1];
            else
                warning('%s folder already exists, but no model.mat file was found.\nThis experiment will be reset.', fullDir.name);
            end
        end
    % end

    % Load model from file or instantiate and initiate new model object
    if (useLearnedFile(1) == 1)
        if (isempty(learnedFile))
            error('learnedFile is an empty string.')
        else
            try
                model = load(learnedFile, 'model');
            catch
                error('Could not open the learned file %s. File seems to be corrupted.', learnedFile);
            end
            model = model.model;
            model.trainTime = trainTime;

            try
                model.savePath = strcat(absoluteDir, '/', fullDir.name);
            catch
                warning('Either absoluteDir or fullDir not defined.\nmodel.savePath remains %s', model.savePath);
            end
        end
    else
        % old static version of config.m
        % model = config(textureFiles, trainTime, testAt, sparseCodingType);

        % for configVar, at first copy values from before ...
        standardParams = {'textureFile', textureFiles, 'trainTime', trainTime, 'testAt', testAt, 'sparseCodingType', sparseCodingType};
        % ... and then add those handled in the function call

        paramVector = [standardParams, inputParams];
        model = configVar(paramVector);
    end

    % check if main script and model are compatible
    if (model.rlModel.continuous == 0)
        error('This training/main script is not compatible with discrete action space models!\nPlease execute OESDiscrete.m instead.');
    elseif (model.rlModel.CActor.output_dim < 2)
        error('This training/main script is not compatible with %d eye muscle models!\nPlease execute OES1Muscle.m instead.', ...
                model.rlModel.CActor.output_dim);
    end

    % safety check for plotting functions
    if (trainTime <= model.interval)
        error('trainTime[%d] must be > model.interval[%d]', trainTime, model.interval);
    end

    % File management: either complete training with existing folder etc., or create a new one
    if ((useLearnedFile(1) == 1) && (useLearnedFile(2) == 1))
        timeToTrain = model.trainTime - model.trainedUntil;
        % cancel experiment if already finished
        if (timeToTrain <= 0)
            warning('timeToTrain = %d, i.e. training finished, therefore training procedure aborted.', timeToTrain);
            if exist(sprintf('%s/modelAt%d/model.mat',model.savePath,model.trainedUntil), 'file') == 2
                return;     % abort simulation if the last testinf phase is also done
            else
                sprintf('... but testing will be started for modelAt%d', model.trainedUntil)
            end
        else
            sprintf('Training procedure will be continued for %d iterations.', timeToTrain)
        end
    else
        if (sparseCodingType > 0)
            modelName = sprintf('%s_%iiter_%s_%i_%s', ...
                                datestr(now, 'yy-mm-dd'), ...
                                trainTime, ...
                                sparseCodingTypeName{sparseCodingType + 1}, ...
                                randomizationSeed, ...
                                experimentName);
        else
            modelName = sprintf('%s_%iiter_%i_%s', ...
                                datestr(now, 'yy-mm-dd'), ...
                                trainTime, ...
                                randomizationSeed, ...
                                experimentName);
        end

        if (~isempty(experimentDirName))
            % folder = sprintf('../results/%s/', experimentDirName);                       % local destination
            folder = sprintf('/home/aecgroup/aecdata/Results/%s/', experimentDirName);     % group folder destination
        else
            % folder = '../results/';                       % local destination
            folder = '/home/aecgroup/aecdata/Results/';     % group folder destination
        end

        mkdir(folder, modelName);
        model.savePath = strcat(folder, modelName);

        % backup all used files
        copyfile(strcat(mfilename, '.m'), model.savePath);
        copyfile('config.m', model.savePath);
        copyfile('configVar.m', model.savePath);
        copyfile(strcat(class(model), '.m'), model.savePath);
        copyfile(strcat(class(model.rlModel), '.m'), model.savePath);
        copyfile(strcat(class(model.rlModel.CCritic), '.m'), model.savePath);
        copyfile(strcat(class(model.rlModel.CActor), '.m'), model.savePath);
        copyfile('results.ods', model.savePath);
        if ((model.rlModel.continuous == 1) && (testIt == 1))
            copyfile('testModelContinuous.m', model.savePath);
        end
        if (clusterCall == 1)
            copyfile('submitClusterJob.sh', model.savePath);
            copyfile('parOES.m', model.savePath);
        end

        timeToTrain = model.trainTime;
        model.inputParams = inputParams;
    end

    % additional notes/information to this model/approach
    model.notes = [model.notes, experimentName];

    %% initialize renderer
    % simulator = OpenEyeSim('create'); % stable renderer
    simulator = OpenEyeSimV5('create'); % experimental version

    simulator.initRenderer();
    % simulator.reinitRenderer(); % for debugging

    % Prepare Textures
    % load all stimuli into memory for experimental renderer
    nTextures = 0;
    tmpTexInd = 1;
    for i = 1 : length(textureFiles)
        texture = load(sprintf('config/%s', textureFiles{i}));
        texture = texture.texture;
        textureCnt = length(texture);

        if (i == 1)
            nTestTextures = textureCnt; % save number of test textures for proper indexing later
        elseif (i == 2)
            nStimTrain = textureCnt;
        end

        nTextures = nTextures + textureCnt;
        textureInd = 1;
        for j = tmpTexInd : nTextures % 140
            simulator.add_texture(j, texture{textureInd});
            textureInd = textureInd + 1;
        end
        tmpTexInd = tmpTexInd + nTextures;
    end
    sprintf('%d textures were added to the buffer from the training and testing sets', nTextures)

    if (nStimTest > nTestTextures)
        warning('%d images were requested as training stimuli, but the renderer only holds %d.', nStimTest, nTestTextures);
        nStimTest = nTestTextures;
    end

    % Image patches cell array (input to model)
    currentView = cell(1, length(model.scModel));

    if clusterCall
        sprintf('Starting the simulation now ...')
    end

    %%% Saturation function that keeps motor commands in [0, 1]
    %   corresponding to the muscelActivity/metabolicCost tables
    function cmd = checkCmd(cmd)
        i0 = cmd < 0;
        cmd(i0) = 0;
        i1 = cmd > 1;
        cmd(i1) = 1;
    end

    %% these parameters can be put in the model:
    % vergErrMax = 0;
    % angleMin = (atand(model.baseline / (2 * model.objDistMax)) * 2) - vergErrMax; %angle for both eyes
    % angleMax = (atand(model.baseline / (2 * model.objDistMin)) * 2) + vergErrMax;
    angleMin = 0;                                                           % corresponds to parallel looking eyes
    angleMax = (atand(model.baseline / (2 * model.objDistMax)) * 2) + 5;    % corr. to maximal obj dist but still in the large scale patches

    %%% Main execution loop
    t = model.trainedUntil; % this is zero in newly initiated model
    command = [0; 0]; %[0.1; 0.1];
    % rewardFunction_prev = 0;
    elapsedTime = 0;
    for iter1 = 1 : (timeToTrain / model.interval)
        % intermediate testing during training
        rngState = rng; % store current state
        if ((testIt == 1) & find(testAt == t)) % have to use single & here, because the last statement is a scalar
            testModelContinuous(model, nStimTest, plotIt(2), 1, 0, simulator, 0, sprintf('modelAt%d', t), [1 : 5]);
            close all;
        end
        rng(rngState); % restore state after testing, to not mess up the experiment

        tic; % start/continue time count

        %% Draw new stimulus
        % currentTexture = texture{(randi(nTextures, 1))};          % stable renderer
        currentTexture = nTestTextures + randi(nStimTrain, 1);      % experimental renderer
        % currentTexture = nStimTest + randi(nStimTrain, 1);

        %% Draw new object depth
        objDist = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);
        angleDes = 2 * atand(model.baseline / (2 * objDist));   % desired vergence [deg]

        %% Initialize muscle activities depending on init method
        if model.initMethod < 3
            if model.initMethod < 1
                % initMethod = 0: completely random muscle commands
                command = rand(2, 1) .* [0.1; 0.2];
                angleNew = model.getAngle(command) * 2;
            else
                % Uniform object fixation distribution
                fixationDist = model.fixDistMin + (model.fixDistMax - model.fixDistMin) * rand(1, 1);

                if model.initMethod == 1
                     % initMethod = 1: Calculate corresponding single muscle activity, i.e. one muscle = 0 activity
                    [command, angleNew] = model.getMF2(fixationDist, 0);
                else
                    % initMethod = 2: Uniform muscle activation distribution for two muscles
                    [command, angleNew] = model.getMFedood(fixationDist, 0);
                end
            end
        else
            % initMethod = 3: add a random vector to the last command
            % the new command should lie within some boundaries, specified by angleMin and angleMax
            dummy = true;
            commandOld = command;
            rangeMax = 0.02;

            while dummy
                [relativeCommand(1, 1), relativeCommand(2, 1)] = pol2cart(rand(1) * pi * 2, rand(1) * rangeMax);
                command = command + relativeCommand;
                command = checkCmd(command);
                angleNew = model.getAngle(command) * 2;

                if (angleNew > angleMin) && (angleNew < angleMax)
                    dummy = false;
                else
                    command = commandOld;
                    rangeMax = rangeMax + 0.005;
                end
            end
        end


        % testing input distribution
        % nSamples = 10000;
        % commands = zeros(nSamples,2);
        % angles = zeros(nSamples, 1);
        % dists = zeros(nSamples, 1);
        % for i = 1:nSamples
        %     initDist = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);
        %     [cmds, angles(i)] = getMF2(initDist, 0);
        %     commands(i, :) = cmds;
        %     initAngle = atand(model.baseline / (2 * initDist));
        %     commands(i) = getMF(initAngle*2);
        %     angles(i) = getAngle([0, commands(i)]);
        %     dists(i) = ((model.baseline / 2)/ (tand(angles(i) / 2)));
        % end
        % figure; histogram(commands); title('commands');
        % figure; histogram(angles); title('angles');
        % figure; histogram(dists); title('distances');

        for iter2 = 1 : model.interval
            % if mod(t, 500) == 0
                % if mod(t, 100000) == 0
                    % investigateMotorSpace(model, simulator, 3, [], 5, 1); % safe the model only every second time
                    % close all;
                % else
                    % investigateMotorSpace(model, simulator, 3, [], 5, 0);
                % end
                % close all;
            % end
            t = t + 1;

            %% Update retinal images
            % refreshImages(currentTexture, angleNew / 2, objDist, 3); % stable renderer
            model.refreshImagesNew(simulator, currentTexture, angleNew / 2, objDist, 3); % experimental renderer

            %% Generate & save the anaglyph picture
            % anaglyph = stereoAnaglyph(imgGrayLeft, imgGrayRight); % only for matlab 2015 or newer
            % imwrite(imfuse(model.imgGrayLeft, model.imgGrayRight, 'falsecolor'), sprintf('%s/anaglyph.png', model.savePath)); %this one works for all tested matlab
            % more advanced functions that generated the anaglyphs of the foveal views
            % generateAnaglyphs(imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS, model.savePath);

            %% Image patch generation
            for i = 1 : length(model.scModel)
                model.preprocessImage(i, 1);
                model.preprocessImage(i, 2);
                currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
            end

            %% Generate basis function feature vector from current images
            [bfFeature, reward, recErrorArray] = model.generateFR(currentView);

            %%% Feedback
            % Generate RL model's input feature vector by
            % basis function feature vector & total muscle command concatination

            if (model.normFeatVect == 0)
                feature = [bfFeature; command * model.lambdaMuscleFB];
            else
                %% Normalized feature vector
                % z-transform raw feature vector (no muscle feedback scaling)
                feature = [bfFeature; command];
                for i = 1 : length(feature)
                    feature(i) = model.onlineNormalize(t, feature(i), i, 1);
                end
                feature = [feature(1 : end - 2); feature(end - 1 : end) * model.lambdaMuscleFB];
            end

            %% bias analysis
            % feature = [feature; 1];

            %%% Calculate metabolic costs
            metCost = model.getMetCost(command) * 2;

            %%% Calculate reward function
            %% Standard reward
            rewardFunction = model.lambdaRec * reward - model.lambdaMet * metCost;

            %% Vergence error reward
            % rewardFunction = -abs(angleDes - angleNew);

            %% Delta reward (part 1)
            % rewardFunctionReal = model.lambdaRec * reward - model.lambdaMet * metCost;
            % rewardFunction = rewardFunctionReal - rewardFunction_prev;

            %% Stasis punishment, i.e. punish non-movement of eyes
            % if (abs(rewardFunctionReal - rewardFunction_prev) < 1e-5)
            %     rewardFunction = rewardFunctionReal - rewardFunction_prev - 1e-5;
            % end

            %% Delta reward (part 2)
            % rewardFunction_prev = rewardFunctionReal;

            %% Normalized reward mean = 0, std = 1
            % rewardFunction = rewardFunction / 100;
            % alpha = model.rlModel.CCritic.gamma; %weighting range (equal to reinforcement running average constant)
            % delta = rewardFunction - model.reward_mean;
            % model.reward_mean = (1 - alpha) * model.reward_mean + (alpha * delta);
            % model.reward_variance = (1 - alpha) * model.reward_variance + (alpha * delta^2);
            % rewardFunction = (rewardFunction - model.reward_mean) / sqrt(model.reward_variance);
            %
            %% norm(reward)
            % rewardFunction = model.onlineNormalize(t, rewardFunction, 3, 1);
            %% norm(recErr)
            % rewardFunction = model.onlineNormalize(t, model.lambdaRec * reward, 1, 1) - model.lambdaMet * metCost;
            %% norm(metCost)
            % rewardFunction = model.lambdaRec * reward - model.onlineNormalize(t, model.lambdaMet * metCost, 2, 1);
            %% norm(recErr) norm(metCost)
            % rewardFunction = model.onlineNormalize(t, model.lambdaRec * reward, 1, 1) + model.onlineNormalize(t, model.lambdaMet * metCost, 2, 1);

            %%% Learning
            %% Sparse coding models
            for i = 1 : length(model.scModel)
                model.scModel{i}.stepTrain();
            end

            % generate delta of muscle activations
            relativeCommand = model.rlModel.stepTrain(feature, rewardFunction, (iter2 > 1));

            % apply the change in muscle activations
            command = command + relativeCommand;    % two muscles
            command = checkCmd(command);            % restrain motor commands to [0, 1]

            % calculate resulting angle which is used for both eyes
            angleNew = model.getAngle(command) * 2;
            % angleNew = getAngle2(command);

            %%% Save statistics
            % compute desired vergence command, disparity and vergence error
            fixDepth = (model.baseline / 2) / tand(angleNew / 2);       % fixation depth [m]
            anglerr = angleDes - angleNew;                              % vergence error [deg]
            % disparity = 2 * model.focalLength * tand(anglerr / 2);    % current disp [px]

            %%% track all movement-related stuff
            model.vergerr_hist(t) = anglerr; % every 10th => adjust displayBasisNEW.m and testModelContinuous.m
            model.relCmd_hist(t, :) = relativeCommand;
            model.cmd_hist(t, :) = command;

            model.td_hist(t) = model.rlModel.CCritic.delta;

            %%% track all reward associated stuff
            model.reward_hist(t) = rewardFunction;
            model.metCost_hist(t) = metCost;
            %  model.recerr_hist(t, :) = recErrorArray; % is tracked every 10th iteration

            %%% track all weight-related stuff
            model.weight_hist(t, 1) = model.rlModel.CCritic.params(1);
            model.weight_hist(t, 2) = model.rlModel.CActor.params(1);
            model.weight_hist(t, 3) = model.rlModel.CActor.params(2);
            model.weight_hist(t, 4) = model.rlModel.CActor.params(3);
            model.weight_hist(t, 5) = model.rlModel.CActor.params(4); % weight change hidden layer
            model.weight_hist(t, 6) = model.rlModel.CActor.params(5); % weight change output layer

            % track all variables that may be decaying
            model.variance_hist(t) = model.rlModel.CActor.variance;
            model.criticLR_hist(t) = model.rlModel.CCritic.alpha_v;
            model.actorLR_hist(t) = model.rlModel.CActor.beta_p;
            model.lambdaMet_hist(t) = model.lambdaMet;

            model.trainedUntil = t;

            %% RL model
            % critic learn rate decay
            if (model.rlModel.criticDecFac > 0)
                %% exponential decay
                % model.rlModel.CCritic.alpha_v = model.rlModel.criticLearningRange(1) * 2 ^ (-t / model.rlModel.criticDecFac);

                %% linear decay
                model.rlModel.CCritic.alpha_v = model.rlModel.CCritic.alpha_v - (model.rlModel.criticDecFac / model.trainTime);
            end

            % actor learn rate decay
            if (model.rlModel.actorDecFac > 0)
                %% exponential decay
                % model.rlModel.CActor.beta_p = model.rlModel.actorLearningRange(1) * 2 ^ (-t / model.rlModel.actorDecFac);

                %% linear decay
                model.rlModel.CActor.beta_p = model.rlModel.CActor.beta_p - (model.rlModel.actorDecFac / model.trainTime);
            end

            % Variance decay, i.e. reduction of actor's output perturbation
            if (model.rlModel.CActor.varDec > 0)
                %% exponential decay
                % model.rlModel.CActor.variance = model.rlModel.CActor.varianceRange(1) * 2 ^ (-t / model.rlModel.CActor.varDec);

                %% linear decay
                model.rlModel.CActor.variance = model.rlModel.CActor.variance - (model.rlModel.CActor.varDec / model.trainTime);
            end

            % metabolic cost decay
            if (model.metCostDec > 0)
                %% exponential decay
                % model.lambdaMet = model.metCostRange(1) * 2 ^ (-t / model.metCostDec);

                %% linear decay
                model.lambdaMet = model.lambdaMet - (model.metCostDec / model.trainTime);
            end

            %% Removed
            % model.verge_actual(t) = angleNew;
            % model.verge_desired(t) = angleDes;
            % model.reward_hist(t) = rewardFunction;
            % model.feature_hist(t, :) = mean(bfFeature);
            % model.Z(t) = objDist;
            % model.fixZ(t) = fixDepth;
            % model.disp_hist(t) = disparity;
        end

        % store every 10th iteration
        model.recerr_hist(t / model.interval, :) = recErrorArray;

        if ((clusterCall == 0 ) && (mod(t, 100) == 0))
            % offers an insight into the models view while it learns
            % imwrite(stereoAanglyph(model.imgGrayLeft, model.imgGrayRight), strcat(model.savePath, '/anaglyph.png'))
            % imwrite(imfuse(model.imgGrayLeft, model.imgGrayRight), strcat(model.savePath, '/anaglyph.png'))

            sprintf('Training Iteration: %d\nObjectDistance:\t%6.2fm\tStart Error:\t%6.3f°\nEnd Fixation:\t%6.2fm\tEnd Error:\t%6.3f°\nMuscle Activations:\t[%.3f, %.3f]\nMean Relative Commands:\t[%.3f, %.3f]', ...
                    t, objDist, model.vergerr_hist(t - model.interval + 1), fixDepth, model.vergerr_hist(t), model.cmd_hist(t, :), ...
                    mean(model.relCmd_hist(t - model.interval + 1 : t, 1)), mean(model.relCmd_hist(t - model.interval + 1 : t, 2)))
        end

        % display per cent completed of training and save model
        if (~mod(t, saveInterval))
            % sprintf('%g%% is finished %d/%d iterations', (t / timeToTrain * 100), t, timeToTrain)
            sprintf('%g%% is finished %d/%d iterations', (t / model.trainTime * 100), t, model.trainTime)
            save(strcat(model.savePath, '/model'), 'model');

            % track basis function history
            if (trackSCBasisHistory == 1)
                for i = 1 : length(model.scModel)
                    model.scModel{i}.saveBasis();
                end
            end
        end

        % if (~mod(t, 10000))
        %     imwrite(imfuse(model.imgGrayLeft, model.imgGrayRight), strcat(model.savePath, sprintf('/anaglyph%d.png', ceil(t / 10000))))
        % end

        elapsedTime = elapsedTime + toc;
    end

    % Total simulation time
    model.simulatedTime = model.simulatedTime + elapsedTime / 60; % [min]
    try
        sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
                elapsedTime / 3600, elapsedTime / 60, elapsedTime, timeToTrain / elapsedTime)
    catch
        warning('elapsedTime = %d, catched erroneous printout.', elapsedTime);
    end

    % store simulated time
    save(strcat(model.savePath, '/model'), 'model');

    %%% Final testing procedure
    if (testIt == 1)
        % testModelContinuous(model, nStim, plotIt, saveTestResults, verbose, simulator, reinitRenderer, experimentDirName, level)
        rngState = rng; % store current state
        testModelContinuous(model, nStimTest, plotIt(2), 1, 0, simulator, 0, sprintf('modelAt%d', t), [1 : 5]);
        rng(rngState); % restore state after testing, to not mess up the experiment

        % print the time again after the line output of the testing script
        try
            sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
                elapsedTime / 3600, elapsedTime / 60, elapsedTime, timeToTrain / elapsedTime)
        catch
            warning('elapsedTime = %d, catched erroneous printout.', elapsedTime);
        end
    end

    % plot results
    if (plotIt(1) == 1)
        if (isempty(testAt))
            model.allPlotSave([1 : 6]); % no test procedure during training -> no testPerformanceVsTraintime plot
        else
            model.allPlotSave([1 : 7]);
        end
    end

    if ((clusterCall == 0) && (closeFigures == 1))
        close all;
    end

    % close the job after completion and release the matlab licence u.u
    if (clusterCall == 0)
        quit
    end
end

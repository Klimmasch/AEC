%%% Main script for launching experimental procedure
% @param trainTime           training time in number of iterations
% @param randomizationSeed   randomization seed
% @param fileDescription     description of approach used as file name
%%%
function OES2Muscles(trainTime, randomizationSeed, fileDescription)
    rng(randomizationSeed);

    % useLearnedFile(1):    0 = don't do it
    %                       1 = use previously learned policy specified in learnedFile
    % useLearnedFile(2):    0 = retrain with same/new parameters
    %                       1 = complete/continue training
    %
    % If useLearnedFile = [1, 1] and you want to continue training, trainTime must be the overall desired train time,
    % i.e. trained at first 100k iterations with useLearnedFile = [0, 0], then decided to continue training for 100k more
    % iterations, then the overall train time for the model is 200k and you set useLearnedFile = [1, 1] and execute with
    % OES2Muscles(200000, randomizationSeed, fileDescription)
    useLearnedFile = [0, 0];
    learnedFile = '';
    % learnedFile = '/home/klimmasch/projects/results/model_05-Jul-2016_21:41:46_1000000_nonhomeo_1_sparseLearning001_finerLS_OD15-6_increasedInit_noMet/model.mat';

    %%% Stimulus declaration
    % textureFile = 'Textures_mcgillManMadeTrain(jpg).mat';       % McGill man made database
    % textureFile = 'Textures_mcgillFruitsAll.mat';             % McGill fruits database
    % textureFile = 'Textures_mcgillFoliageTrain(jpg).mat';     % McGill foliage database
    % textureFile = 'Textures_vanHaterenTrain.mat';             % vanHateren database
    % textureFile = 'Textures_celine.mat';                      % Celine's images

    % for the new renderer, all textures to be used during training and
    % testing have to be loaded into the buffer at the beginning
    % per convention, the testing images are given in the first entry!!
    textureFiles = {'Textures_mcgillManMade40.mat', 'Textures_mcgillManMade100.mat'};

    %%% executing the test procedure during training?
    testAt = [500000 : 500000 : trainTime];

    %%% Sparse coding approach
    % sparseCodingType: 0 = non-homeostatic
    %                   1 = homeostatic
    sparseCodingType = uint8(0);
    % sparseCodingTypeName = cellstr(['nonhomeo'; 'homeo___']);

    %%% Plotting flag
    % Whether figures should be generated and saved
    % plotIt: [training, testing]
    %            0 = don't do it
    %            1 = do it
    plotIt = [uint8(1), uint8(1)];

    %%% Whether figures should be closed after generation
    % closeFigures: 0 = don't do it
    %               1 = do it
    closeFigures = uint8(1);

    %%% Testing flag
    % Whether the testing procedure shall be executed after training
    % testIt:   0 = don't do it
    %           1 = do it
    testIt = uint8(1);

    %%% Amount of test stimuli
    nStimTest = 100;

    % Load model from file or instantiate and initiate new model object
    if (useLearnedFile(1) == 1)
        if isempty(learnedFile)
            sprintf('could not open the learned file! %s', learnedFile)
            return;
        else
            model = load(learnedFile, 'model');
            model = model.model;
            model.trainTime = trainTime;
        end
    else
        model = config(textureFiles(2), trainTime, sparseCodingType);
    end

    % check if main script and model are compatible
    if (model.rlModel.continuous == 0)
        sprintf('Error: This training/main script is not compatible with discrete action space models!\nPlease execute OESDiscrete.m instead.')
        return;
    elseif (model.rlModel.CActor.output_dim < 2)
        sprintf('Error: This training/main script is not compatible with %d eye muscle models!\nPlease execute OES1Muscle.m instead.', ...
                model.rlModel.CActor.output_dim)
        return;
    end

    % safety check for plotting functions
    if (trainTime <= model.interval)
        sprintf('Error: trainTime[%d] must be > model.interval[%d]', trainTime, model.interval)
        return;
    end

    % File management: either complete training with existing folder etc., or create a new one
    if ((useLearnedFile(1) == 1) && (useLearnedFile(2) == 1))
        timeToTrain = model.trainTime - model.trainedUntil;
    else
        % modelName = sprintf('model_%s_%i_%s_%i_%s', ...
        %                     datestr(now, 'dd-mmm-yyyy_HH:MM:SS'), ...
        %                     trainTime, ...
        %                     sparseCodingTypeName{sparseCodingType + 1}, ...
        %                     randomizationSeed, ...
        %                     fileDescription);
        modelName = sprintf('model_%s_%i_%i_%s', ...
                            datestr(now, 'dd-mmm-yyyy_HH:MM:SS'), ...
                            trainTime, ...
                            randomizationSeed, ...
                            fileDescription);
        folder = '../results/';
        mkdir(folder, modelName);
        model.savePath = strcat(folder, modelName);

        % backup all used files
        copyfile(strcat(mfilename, '.m'), model.savePath);
        copyfile('config.m', model.savePath);
        copyfile(strcat(class(model), '.m'), model.savePath);
        copyfile(strcat(class(model.rlModel), '.m'), model.savePath);
        copyfile(strcat(class(model.rlModel.CCritic), '.m'), model.savePath);
        copyfile(strcat(class(model.rlModel.CActor), '.m'), model.savePath);
        if (model.rlModel.continuous == 1)
            copyfile('testModelContinuous.m', model.savePath);
        end
        copyfile('results.ods', model.savePath);

        timeToTrain = model.trainTime;
    end

    % additional notes/information to this model/approach
    model.notes = [model.notes, fileDescription];

    % Save model every #saveInterval training iterations
    % saveInterval = ceil(model.trainTime / 10);
    saveInterval = trainTime;

    % Track the evolution of all basis functions of the respective sparse coders
    trackSCBasisHistory = uint8(0);

    degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
    metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

    % muscle function :=  mf(vergence_angle) = muscle force [single muscle]
    resolution = 100001;
    approx = spline(1 : 11, degrees.results_deg(:, 1));

    xValPos = ppval(approx, 1 : 0.0001 : 11)';
    yValPos = linspace(0, 1, resolution)';

    % xValNeg = flipud(ppval(approx, 1 : 0.0001 : 11)' * -1);
    % yValNeg = linspace(-1, 0, resolution)';

    % mfunction = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
    mfunction = [xValPos, yValPos];
    mfunction(:, 1) = mfunction(:, 1) * 2;  % angle for two eyes
    dmf = abs(diff(mfunction(1 : 2, 1)));   % delta in angle
    dmf2 = diff(mfunction(1 : 2, 2));       % delta in mf

    approx = spline(1 : 11, degrees.results_deg(1, :));
    xValPos = ppval(approx, 1 : 0.0001 : 11)';
    yValPos = linspace(0, 1, resolution)';

    mfunction2 = [xValPos, yValPos];
    mfunction2(:, 1) = mfunction2(:, 1) * 2;    % angle for two eyes
    dmf3 = abs(diff(mfunction2(1 : 2, 1)));     % delta in angle

    %%% New renderer
    % simulator = OpenEyeSim('create'); % stable renderer
    simulator = OpenEyeSimV5('create'); % experimental version

    simulator.initRenderer();
    % simulator.reinitRenderer(); % for debugging

    % Prepare Textures
    % load all stimuli into memory for experimental renderer
    nTextures = 0;
    tmpTexInd = 1;
    for i = 1:length(textureFiles)
        texture = load(sprintf('config/%s', textureFiles{i}));
        texture = texture.texture;
        nTextures = nTextures + length(texture);
        if i == 1
            nTestTextures = nTextures;      % save number of test textures for proper indexing later
        end
        textureInd = 1;
        for j = tmpTexInd : nTextures % 140
            simulator.add_texture(j, texture{textureInd});
            textureInd = textureInd + 1;
        end
        tmpTexInd = tmpTexInd + nTextures;
    end
    sprintf('%d textures were added to the buffer from the training and testing files', nTextures)

    if nStimTest > nTestTextures
        nStimTest = nTestTextures;
        nStimTrain = nTextures - nStimTest;
        sprintf('The file for testing textures only contains %d images, but I will use them all!', nStimTest)
    end

    imgRawLeft = uint8(zeros(240, 320, 3));
    imgRawRight = uint8(zeros(240, 320, 3));
    imgGrayLeft = uint8(zeros(240, 320, 3));
    imgGrayRight = uint8(zeros(240, 320, 3));

    % Image patches cell array (input to model)
    currentView = cell(1, length(model.scModel));

    %%% Generates two new images for both eyes
    % texture:      file path of texture input
    % eyeAngle:     angle of single eye (rotation from offspring)
    % objDist:      distance of stimulus
    % scaleImSize:  scaling factor of stimulus plane [m]
    function refreshImages(texture, eyeAngle, objDist, scaleImSize)
        simulator.add_texture(1, texture);
        simulator.set_params(1, eyeAngle, objDist, 0, scaleImSize); % scaling of obj plane size

        result1 = simulator.generate_left();
        result2 = simulator.generate_right();

        imgRawLeft = permute(reshape(result1, ...
                                     [size(imgRawLeft, 3), ...
                                      size(imgRawLeft, 2), ...
                                      size(imgRawLeft, 1)]), ...
                                     [3, 2, 1]);

        imgRawRight = permute(reshape(result2, ...
                                      [size(imgRawRight, 3), ...
                                       size(imgRawRight, 2), ...
                                       size(imgRawRight, 1)]), ...
                                      [3, 2, 1]);

        % convert images to gray scale
        imgGrayLeft = 0.2989 * imgRawLeft(:, :, 1) + 0.5870 * imgRawLeft(:, :, 2) + 0.1140 * imgRawLeft(:, :, 3);
        imgGrayRight = 0.2989 * imgRawRight(:, :, 1) + 0.5870 * imgRawRight(:, :, 2) + 0.1140 * imgRawRight(:, :, 3);
    end

    %%% Generates two new images for both eyes for experimental renderer
    % textureNumber:    index of stimulus in memory buffer
    % eyeAngle:         angle of single eye (rotation from offspring)
    % objDist:          distance of stimulus
    % scaleImSize:  scaling factor of stimulus plane [m]
    function refreshImagesNew(textureNumber, eyeAngle, objDist, scaleImSize)
        simulator.set_params(textureNumber, eyeAngle, objDist, 0, scaleImSize); % scaling of obj plane size

        result1 = simulator.generate_left();
        result2 = simulator.generate_right();

        imgRawLeft = permute(reshape(result1, ...
                                     [size(imgRawLeft, 3), ...
                                      size(imgRawLeft, 2), ...
                                      size(imgRawLeft, 1)]), ...
                                     [3, 2, 1]);

        imgRawRight = permute(reshape(result2, ...
                                      [size(imgRawRight, 3), ...
                                       size(imgRawRight, 2), ...
                                       size(imgRawRight, 1)]), ...
                                      [3, 2, 1]);

        % convert images to gray scale
        imgGrayLeft = 0.2989 * imgRawLeft(:, :, 1) + 0.5870 * imgRawLeft(:, :, 2) + 0.1140 * imgRawLeft(:, :, 3);
        imgGrayRight = 0.2989 * imgRawRight(:, :, 1) + 0.5870 * imgRawRight(:, :, 2) + 0.1140 * imgRawRight(:, :, 3);
    end

    %%% Helper function that maps {vergenceAngle} -> {muscleForce}
    function mf = getMF(vergAngle)
        % look up index of vergAngle
        if (vergAngle >= degrees.results_deg(1, 1) * 2)
            indVergAngle = find(mfunction(:, 1) <= vergAngle + dmf & mfunction(:, 1) >= vergAngle - dmf);
            mf = mfunction(indVergAngle, 2);
            mf = mf(ceil(length(mf) / 2));
        else
            mf = 0;
        end
    end

    %%% Helper function that maps muscle activities to resulting angle (for
    %%% one eye)
    function angle = getAngle(command)
        cmd = (command * 10) + 1;                                       % scale commands to table entries
        angle = interp2(degrees.results_deg, cmd(1), cmd(2), 'spline'); % interpolate in table by cubic splines
    end

    function angle = getAngle2(command)
        angleIndex = find(mfunction(:, 2) <= command(2) + dmf2 & mfunction(:, 2) >= command(2) - dmf2);
        angle = mfunction(angleIndex, 1);
        angle = angle(ceil(length(angle) / 2));
    end

    %%% Helper function that maps muscle activities to resulting metabolic costs
    function tmpMetCost = getMetCost(command)
        cmd = (command * 10) + 1;                                           % scale commands to table entries
        tmpMetCost = interp2(metCosts.results, cmd(1), cmd(2), 'spline');   % interpolate in table by cubic splines
    end

    %%% Saturation function that keeps motor commands in [0, 1]
    %   corresponding to the muscelActivity/metabolicCost tables
    function cmd = checkCmd(cmd)
        i0 = cmd < 0;
        cmd(i0) = 0;
        i1 = cmd > 1;
        cmd(i1) = 1;
    end

    % Calculates muscle force for two muscles
    function [mf, angleInit] = getMF2(objDist, desVergErr)
        % correct vergence angle for given object distance
        angleCorrect = 2 * atand(model.baseline / (2 * objDist));
        % desired init angle for given vergence error [deg]
        angleInit = angleCorrect - desVergErr;
        % look up index of angleInit
        % if objDist not fixateable with medial rectus, use lateral rectus
        if (angleInit >= mfunction(1, 1))
            indAngleInit = find(mfunction(:, 1) <= angleInit + dmf & mfunction(:, 1) >= angleInit - dmf);
            mf = mfunction(indAngleInit, 2);
            mf = [0; mf(ceil(length(mf) / 2))];
        else
            indAngleInit = find(mfunction2(:, 1) <= angleInit + dmf3 & mfunction2(:, 1) >= angleInit - dmf3);
            mf = mfunction2(indAngleInit, 2);
            mf = [mf(ceil(length(mf) / 2)); 0];
        end
    end

    %%% Main execution loop
    t = model.trainedUntil; % this is zero in newly initiated model
    command = [0; 0];
    % rewardFunction_prev = 0;
    elapsedTime = 0;
    tic; % start time count
    for iter1 = 1 : (timeToTrain / model.interval)
        % pick random texture every #interval times
        % currentTexture = texture{(randi(nTextures, 1))};  % stable renderer
        currentTexture = nTestTextures + randi(nStimTrain, 1);               % experimental renderer

        % random depth
        objDist = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);
        angleDes = 2 * atand(model.baseline / (2 * objDist));   % desired vergence [deg]

        % reset muscle activities to random values
        % initialization for muscle in between borders of desired actvity
        % i.e. min and max stimulus distance
        % command(1) = 0; % single muscle
        % command(1) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1, 1); % two muscles
        % command(2) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1, 1);
        % initDist = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);
        % command(2) = getMF(2 * atand(model.baseline / (2 * initDist)));
        % command(2) = getMF(model.vergAngleMin + (model.vergAngleMax - model.vergAngleMin) * rand(1, 1));
        % command(1) = model.muscleInitMin(1) + (model.muscleInitMax(1) - model.muscleInitMin(1)) * rand(1, 1);
        % command(2) = model.muscleInitMin(2) + (model.muscleInitMax(2) - model.muscleInitMin(2)) * rand(1, 1);
        fixationDist = model.fixDistMin + (model.fixDistMax - model.fixDistMin) * rand(1, 1);
        [command, angleNew] = getMF2(fixationDist, 0);

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

        % angleNew = getAngle(command) * 2;
        % angleNew = getAngle2(command);

        for iter2 = 1 : model.interval
            t = t + 1;

            % update stimuli
            % refreshImages(currentTexture, angleNew / 2, objDist, 3);  % stable renderer
            refreshImagesNew(currentTexture, angleNew / 2, objDist, 3); % experimental renderer

            % Generate & save the anaglyph picture
            % anaglyph = stereoAnaglyph(imgGrayLeft, imgGrayRight); % only for matlab 2015 or newer
            % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [model.savePath '/anaglyph.png']); %this one works for all tested matlab
            % more advanced functions that generated the anaglyphs of the foveal views
            % generateAnaglyphs(imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS, model.savePath);

            % Image patch generation
            for i = 1 : length(model.scModel)
                model.preprocessImage(imgGrayLeft, i, 1);
                model.preprocessImage(imgGrayRight, i, 2);
                currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
            end

            % Generate basis function feature vector from current images
            [bfFeature, reward, recErrorArray] = model.generateFR(currentView);

            %%% Feedback
            % Generate RL model's input feature vector by
            % basis function feature vector & total muscle command concatination
            feature = [bfFeature; command * model.lambdaMuscleFB];

            %%% Calculate metabolic costs
            metCost = getMetCost(command) * 2;

            %%% Calculate reward function
            %% Standard reward
            rewardFunction = model.lambdaRec * reward - model.lambdaMet * metCost;

            %% Vergence error reward
            % rewardFunction = -abs(angleDes - angleNew);

            %% Delta reward
            % rewardFunctionReal = model.lambdaRec * reward - model.lambdaMet * metCost;
            % rewardFunction = rewardFunctionReal - rewardFunction_prev;

            % Stasis punishment, i.e. punish non-movement of eyes
            % if (abs(rewardFunctionReal - rewardFunction_prev) < 1e-5)
            %     rewardFunction = rewardFunctionReal - rewardFunction_prev - 1e-5;
            % end

            % rewardFunction_prev = rewardFunctionReal;

            %%% Learning
            %% Sparse coding models
            for i = 1 : length(model.scModel)
                model.scModel{i}.stepTrain();
            end

            %% RL model
            % Variance decay, i.e. reduction of actor's output perturbation
            % if (model.rlModel.CActor.varDec > 0)
            %     model.rlModel.CActor.variance = model.rlModel.CActor.varianceRange(1) * 2 ^ (-t / model.rlModel.CActor.varDec);
            % end

            relativeCommand = model.rlModel.stepTrain(feature, rewardFunction, (iter2 > 1));

            % add the change in muscle Activities to current ones
            command = command + relativeCommand;    % two muscles
            command = checkCmd(command);            % restrain motor commands to [0, 1]

            % calculate resulting angle which is used for both eyes
            angleNew = getAngle(command) * 2;
            % angleNew = getAngle2(command);

            %%%%%%%%%%%%%%%% TRACK ALL PARAMETERS %%%%%%%%%%%%%%%%%%

            % compute desired vergence command, disparity and vergence error
            fixDepth = (model.baseline / 2) / tand(angleNew / 2);   % fixation depth [m]
            anglerr = angleDes - angleNew;                          % vergence error [deg]
            disparity = 2 * model.focalLength * tand(anglerr / 2);  % current disp [px]

            % save state
            model.Z(t) = objDist;
            model.fixZ(t) = fixDepth;
            model.disp_hist(t) = disparity;
            model.vergerr_hist(t) = anglerr;
            model.recerr_hist(t, :) = recErrorArray;
            model.verge_actual(t) = angleNew;
            model.verge_desired(t) = angleDes;
            model.relCmd_hist(t, :) = relativeCommand;
            model.cmd_hist(t, :) = command;
            model.reward_hist(t) = rewardFunction;
            % model.feature_hist(t) = mean(bfFeature);
            model.metCost_hist(t) = metCost;
            model.td_hist(t) = model.rlModel.CCritic.delta;
            model.weight_hist(t, 1) = model.rlModel.CCritic.params(1);
            model.weight_hist(t, 2) = model.rlModel.CActor.params(1);
            model.weight_hist(t, 3) = model.rlModel.CActor.params(2);
            model.weight_hist(t, 4) = model.rlModel.CActor.params(3);
            model.weight_hist(t, 5) = model.rlModel.CActor.params(4); % weight change hidden layer
            model.weight_hist(t, 6) = model.rlModel.CActor.params(5); % weight change output layer
            model.variance_hist(t) = model.rlModel.CActor.variance;

            model.trainedUntil = t;
        end

        if (mod(t, 100) == 0)
            % infos = {t, objDist, model.vergerr_hist(t - model.interval + 1), angleDes - angleNew, command', relativeCommand', reward, recErrorArray};
            % generateAnaglyphs(t, 1, infos);
            imwrite(stereoAnaglyph(imgGrayLeft, imgGrayRight), strcat(model.savePath, '/anaglyph.png')) % offers an insight into the models view while it learns

            sprintf('Training Iteration: %d\nObjectDistance:\t%6.2fm\tStart Error:\t%6.3f°\nEnd Fixation:\t%6.2fm\tEnd Error:\t%6.3f°\nMuscle Activations:\t[%.3f, %.3f]\nMean Relative Commands:\t[%.3f, %.3f]', ...
                    t, objDist, model.vergerr_hist(t - model.interval + 1), fixDepth, model.vergerr_hist(t), model.cmd_hist(t, :), ...
                    mean(model.relCmd_hist(t - model.interval + 1 : t, 1)), mean(model.relCmd_hist(t - model.interval + 1 : t, 2)))
        end

        elapsedTime = elapsedTime + toc;

        % testing during training!
        if (testIt & find(testAt == t)) % have to use single & here, because the last statement is a scalar
            testModelContinuous(model, nStimTest, plotIt(2), 1, simulator, 0, sprintf('modelAt%d', t));
            close all;
        end

        tic; % skip testing measure

        % Display per cent completed of training and save model
        if (~mod(t, saveInterval))
            sprintf('%g%% is finished', (t / timeToTrain * 100))
            save(strcat(model.savePath, '/model'), 'model');

            % track basis history
            if (trackSCBasisHistory == 1)
                for i = 1 : length(model.scModel)
                    model.scModel{i}.saveBasis();
                end
            end
        end
    end
    elapsedTime = elapsedTime + toc;

    % Total simulation time
    model.simulatedTime = model.simulatedTime + elapsedTime / 60;
    sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
            elapsedTime / 3600, elapsedTime / 60, elapsedTime, timeToTrain / elapsedTime)

    % store simulated time
    save(strcat(model.savePath, '/model'), 'model');

    % plot results
    if (plotIt(1) == 1)
        model.allPlotSave();
    end

    %%% Testing procedure
    if (testIt && isempty(testAt))
        % testModelContinuous(model, nStim, plotIt, saveTestResults, simulatorHandle, reinitRenderer, tempResultsFolderName)
        testModelContinuous(model, nStimTest, plotIt(2), 1, simulator, 0, sprintf('modelAt%d', t));
    end

    if (closeFigures == 1)
        close all;
    end


    % Generates anaglyphs of the large and small scale fovea and
    % one of the two unpreprocessed gray scale images
    % function generateAnaglyphs(leftGray, rightGray, dsRatioL, dsRatioS, foveaL, foveaS, savePath, identifier, markScales, infos)
    function generateAnaglyphs(identifier, markScales, infos)

        numberScales = length(model.scModel);

        % defining colors in the image: (from larges to smallest scale)
        scalingColors = {'blue', 'red', 'green'};
        textColor = 'yellow';

        scaleImages = cell(numberScales);

        for scale = 1 : numberScales
            % Downsampling Large
            imgLeft = imgGrayLeft(:);
            imgLeft = reshape(imgLeft, size(imgGrayLeft));
            imgRight = imgGrayRight(:);
            imgRight = reshape(imgRight, size(imgGrayRight));

            for ds = 1 : log2(model.dsRatio(scale))
                imgLeft = impyramid(imgLeft, 'reduce');
                imgRight = impyramid(imgRight, 'reduce');
            end

            % cut fovea in the center
            [h, w, ~] = size(imgLeft);
            cutOutVertical = [fix(h / 2 + 1 - model.pxFieldOfView(scale) / 2), fix(h / 2 + model.pxFieldOfView(scale) / 2)];
            cutOutHorizontal = [fix(w / 2 + 1 - model.pxFieldOfView(scale) / 2), fix(w / 2 + model.pxFieldOfView(scale) / 2)];


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
        anaglyph = imfuse(imgGrayLeft, imgGrayRight, 'falsecolor');
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
        %% todo: insert reconstruction error for all and for each scale
        %% infos = {iteration, tmpObjRange, vseRange(vseIndex), angleDes - angleNew, command', relativeCommand', reward, recErrorArray};
        xPos = [10, 220, 10, 10, 260, 10, 10, 240]; %display in headful modus: [10, 200, 10, 10, 260, 10, 10]
        yPos = [10, 10, 230, 220, 230, 190, 200, 220];
        % imName = strsplit(currentTexture, '/');           % stable renderer
        imName = strsplit(texture{currentTexture}, '/');    % experimental renderer
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
        saveas(fig, sprintf('%s/anaglyph.png', model.savePath), 'png');
        fig.delete();
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

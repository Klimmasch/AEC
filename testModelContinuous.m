%%% Model testing procedure
% param model:              respective model object to be tested
% param nStim:              # stimuli to be tested, provide 0 if you just want to plot results
% pram plotIt:              whether plots shall be generated
% param saveTestResults:    whether to save the results (not recommended if model is still trained!)
% param verbose:            generation of text output 0 (no) 1 (yes)
% param simulator:          simulator handle, provide [] if there is no simulator handle yet
% param reinitRenderer:     1 if renderer was already initialized
%                           0 if renderer wasn't initialized yet
% param folderName:         subfolder name of this testing instance, default: 'modelAtX', X = completed training # iteration
% param level:              array of level numbers, i.e. which testing/plotting parts shall be executed
%                           elem. [1, 6] (1: vergenceStartErrors, 2: reconstrErrCritic, 3: ObjDistvsFixDist, 4: approachings,
%                           5: trajectory, 6: critics response of muscle plane)
function testModelContinuous(model, nStim, plotIt, saveTestResults, verbose, simulator, reinitRenderer, folderName, level)

    % should the simulation time be measured?
    % not recommended for use during training
    measureTime = false;

    % Vergence error resolution for 2nd testing procedure
    % Needs to be odd number to include vergErr = 0
    test2Resolution = 101;

    %%% Stimulus declaration
    % textureFile = 'Textures_mcgillManMadeTrain(jpg).mat';     % McGill man made database
    % textureFile = 'Textures_mcgillManMadeTest(jpg).mat';
    textureFile = 'Textures_mcgillManMade40.mat';
    % textureFile = 'Textures_mcgillFruitsAll(jpg).mat';        % McGill fruits database
    % textureFile = 'Textures_mcgillFoliageTrain(jpg).mat';     % McGill foliage database
    % textureFile = 'Textures_mcgillFoliageTest(jpg).mat';
    % textureFile = 'Textures_vanHaterenTrain.mat';             % vanHateren database
    % textureFile = 'Textures_vanHaterenTest.mat';
    % textureFile = 'Textures_celine.mat';                      % Celine's images

    % cancel testing procedure
    if ((nStim == 0) && (isempty(model.testResult)))
        error('Model has no testResults!');
    elseif (nStim == 1)
        error('nStim must be != 1');
    end

    % user handling
    if (isempty(level))
        level = 1 : 6;
    elseif (level(1) < 1)
        error('level(1) = %d must be > 0!', level(1));
    end

    % backward compatibility for feat. norm.
    % TODO: remove as soon as it won't be needed
    if ((isempty(model.normFeatVect)) && (isempty(model.desiredStdZT)))
        model.normFeatVect = 0;
        model.desiredStdZT = 1;
    elseif ((isempty(model.normFeatVect)) && (model.desiredStdZT ~= 1))
        model.normFeatVect = 1;
    elseif ((isempty(model.normFeatVect)) && (model.desiredStdZT == 1))
        model.normFeatVect = 0;
    end

    %%% New renderer
    if (isempty(simulator) && (nStim ~= 0))
        % simulator = OpenEyeSim('create'); % stable renderer
        simulator = OpenEyeSimV5('create'); % latest renderer

        if (reinitRenderer == 0)
            simulator.initRenderer();
        else
            % for debugging purposes
            simulator.reinitRenderer();
        end

        % load all stimuli into memory for experimental renderer, if no
        % renderer is provided by the training procedure
        texture = load(sprintf('config/%s', textureFile));
        texture = texture.texture;
        for i = 1 : nStim
            simulator.add_texture(i, texture{i});
        end
        sprintf('%d textures added to the testing simulator', nStim)
    end

    %%% creating a new directory if (folderName ~= '/.')
    if (folderName(1) ~= '/')
        folderName = ['/' folderName];
    end
    imageSavePath = [model.savePath folderName];
    mkdir(imageSavePath);

    if (saveTestResults == 1)
        save(strcat(imageSavePath, '/model'), 'model');
    end

    % fixation interval at testing procedure
    if (isempty(model.testInterval))
        model.testInterval = model.interval * 2;
        % model.testInterval = 200;
    end
    
    strabAng = model.strabAngle; % copy for saving
    model.strabAngle = 0; % remove strabismic angle for testing only
    
    command = [0; 0];
%     objRange = [model.objDistMin : 0.5 : model.objDistMax];
%     if (model.objDistMin == 0.5 && model.objDistMax == 6)
%         objRange = [0.5, 1 : 6];
%     end
    objRange = [0.5, 1 : 6]; % same conditions for all the models during testing

    tmpResult1 = zeros(nStim, model.testInterval + 1);
    tmpResult2 = zeros(nStim, model.testInterval + 1);
    tmpResult3 = zeros(nStim, model.testInterval + 1);

    testResult = zeros(length(objRange), 7, 6 * model.testInterval + 6); % mean and std of vergenceError, deltaMF and critics response for different obj dists and starting pos
    testResult2 = zeros(length(objRange) * 7 * nStim * model.testInterval, 1 + length(model.scModel)); % reconstruction error statistics

    testResult3 = zeros(length(objRange) * 7 * nStim, model.testInterval); % ALL single values
    testResult4 = zeros(length(objRange), test2Resolution, nStim * (2 + length(model.scModel)));
    testResult5 = zeros(length(objRange) * 7 * nStim * model.testInterval, model.rlModel.CActor.output_dim * 2); % correlation between abs muscle activations and deltaMFs
    testResult6 = zeros(model.testInterval * 10, 2);
    testResult7 = zeros(length(objRange) * 7 * nStim, model.testInterval); % ALL single values for metCost
    vergenceAngleApproach = []; % is defined further below after some parameters are set
    metCostsApproach = [];
    musclePlaneResponse = []; % is defined further below after some parameters are set

    % here, the images are safed that start at the maximal vergence errors (neg & pos) and that end up worse than they started
    % this tabular is going to be safed inside the models folder and
    % histograms will be generated
    reallyBadImages = zeros(2, length(objRange), nStim);

    % Image patches cell array (input to model)
    currentView = cell(1, length(model.scModel));

    %%% Saturation function that keeps motor commands in [0, 1]
    %   corresponding to the muscelActivity/metabolicCost tables
    function [cmd] = checkCmd(cmd)
        i0 = cmd < 0;
        cmd(i0) = 0;
        i1 = cmd > 1;
        cmd(i1) = 1;
    end

    tr2Ind = 1;
    tr3Ind = 1;
    tr5Ind = 1;
    tr7Ind = 1;
    if (measureTime == true)
        tic;
    end

    % don't repeat testing procedure if nStim == 0, but just plot the results
    if (nStim > 0)
        sprintf('Testing procedure at iter = %s started...', folderName(9 : end))

        %%% Vergence behaviour for different vergence start errors
        if ((~isempty(find(level == 1))) || (~isempty(find(level == 4))))
            if (verbose == 1)
                sprintf('Level 1/4 Vergence behaviour for different vergence start errors')
            end

            for odIndex = 1 : length(objRange)
                if (verbose == 2)
                    sprintf('Level 1/4 Test iteration = %d/%d', odIndex, size(objRange, 2) * 2)
                end

                % vergence start error
                vergMax = model.getVergErrMax(objRange(odIndex));
                if vergMax > 2
                    vergMax = 2;
                end
                vseRange = [linspace(-2, 0, 4), linspace(0, vergMax, 4)];
                vseRange = [vseRange(1 : 3), vseRange(5 : end)]; % remove one 0
                % vseRange = [-3:3];
                % vseRange = linspace(-1, 1, 7);
                [cmdDesired, angleDes] = model.getMF2(objRange(odIndex), 0);
                metCostDesired = model.getMetCost(cmdDesired) * 2;

                for vseIndex = 1 : length(vseRange)
                    tmpResult1(:, 1) = vseRange(vseIndex);

                    for stimulusIndex = 1 : nStim
                        % currentTexture = texture{stimulusIndex};  % stable renderer
                        currentTexture = stimulusIndex;             % experimental renderer

                        % Calculate corresponding single muscle activity, i.e. one muscle = 0 activity
                        % [command, angleNew] = model.getMF2(objRange(odIndex), vseRange(vseIndex));

                        % Uniform muscle activation distribution for two muscles
                        [command, angleNew] = model.getMFedood(objRange(odIndex), vseRange(vseIndex));
                        randForLeftFilt = rand(1,1);
                        randForRightFilt = rand(1,1);

                        for iter = 2 : model.testInterval + 1
                            % update stimuli
                            % refreshImages(currentTexture, angleNew / 2, objRange(odIndex), 3);                    % stable renderer
                            model.refreshImagesNew(simulator, currentTexture, angleNew / 2, objRange(odIndex), 3, [0,0,0]);  % experimental renderer

                            %% change left and right images to simulate altered rearing conditions
                            % if ~isempty(model.filterLeft)
                                % if randForLeftFilt < model.filterLeftProb
                                    % model.imgGrayLeft = conv2(model.imgGrayLeft, model.filterLeft, 'same');
                                    % sligthly faster version
                                    % model.imgGrayLeft = conv2(model.filterLeft{1}, model.filterLeft{2}, model.imgGrayLeft, 'same');
                                % end
                            % end
                            % if ~isempty(model.filterRight)
                                % if randForRightFilt < model.filterRightProb
                                    % model.imgGrayRight = conv2(model.imgGrayRight, model.filterRight, 'same');
                                    % slightly faster version
                                    % model.imgGrayRight = conv2(model.filterRight{1}, model.filterRight{2}, model.imgGrayRight, 'same');
                                % end
                            % end

                            % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [imageSavePath '/anaglyph.png']);

                            % Image patch generation
                            for i = 1 : length(model.scModel)
                                model.preprocessImage(i, 1);
                                model.preprocessImage(i, 2);
                                currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
                            end

                            % Generate input feature vector from current images
                            [bfFeature, ~, recErrorArray] = model.generateFR(currentView);

                            % Track reconstruction error statistics
                            testResult2(tr2Ind, :) = [angleDes - angleNew, recErrorArray];
                            tr2Ind = tr2Ind + 1;

                            % Absolute command feedback # concatination
                            if (model.normFeatVect == 0)
                                if (model.rlModel.continuous == 1)
                                    if (model.rlModel.CActor.output_dim == 1)
                                        feature = [bfFeature; command(2) * model.lambdaMuscleFB];   % single muscle
                                    else
                                        feature = [bfFeature; command * model.lambdaMuscleFB];      % two muscles
                                    end
                                else
                                    feature = bfFeature;
                                end
                            else
                                if (model.rlModel.continuous == 1)
                                    if (model.rlModel.CActor.output_dim == 1)
                                        % z-transformed raw feature vector (no muscle feedback scaling)
                                        feature = [bfFeature; command(2)];
                                    else
                                        % z-transformed raw feature vector (no muscle feedback scaling)
                                        feature = [bfFeature; command];
                                    end
                                else
                                    feature = bfFeature;
                                end

                                % z-transformed raw feature vector (no muscle feedback scaling)
                                for i = 1 : length(feature)
                                    feature(i) = model.onlineNormalize(model.trainedUntil, feature(i), i, 0);
                                end
                                % post normalization muscle feedback scaling
                                feature = [feature(1 : end - 2); feature(end - 1 : end) * model.lambdaMuscleFB];
                            end

                            %% bias analysis
                            if (model.rlModel.bias > 0)
                                feature = [feature; model.rlModel.bias];
                            end

                            %%% Calculate metabolic costs
                            % metCost = getMetCost(command) * 2;

                            %%% Action
                            relativeCommand = model.rlModel.act(feature);

                            % add the change in muscle activities to current ones
                            if (model.rlModel.continuous == 1)
                                if (model.rlModel.CActor.output_dim == 1)
                                    command(2) = command(2) + relativeCommand;  % single muscle
                                else
                                    command = command + relativeCommand;        % two muscles
                                end
                                command = checkCmd(command);                    % restrain motor commands to [0,1]
                                angleNew = model.getAngle(command) * 2;         % resulting angle is used for both eyes
                            else
                                angleNew = angleNew + relativeCommand;
                                if (angleNew > angleMax || angleNew < angleMin)
                                    angleNew = model.vergAngleFixMin + (model.vergAngleFixMax - model.vergAngleFixMin) * rand(1, 1);
                                end
                            end

                            % track commands for correlation plot
                            % and metabolic costs
                            if (model.rlModel.CActor.output_dim >= 2)
                                testResult5(tr5Ind, :) = [command; relativeCommand]';
                                tr5Ind = tr5Ind + 1;

                                testResult7(tr7Ind, iter - 1) = model.getMetCost(command) * 2 - metCostDesired;
                            end

                            % track bad or redundant stimuli
                            % if (iter == 11)
                            %     if (abs(angleDes - angleNew) > 0.5)
                            %         sprintf('VergErr = %.1f\timage = %s\tstimulusIndex = %d\tobjDist = %.2f', (angleDes - angleNew), currentTexture, stimulusIndex, objRange(odIndex))
                            %     end
                            % end

                            % temporary results
                            tmpResult1(stimulusIndex, iter) = angleDes - angleNew;
                            if (model.rlModel.CActor.output_dim == 2)
                                tmpResult2(stimulusIndex, iter) = relativeCommand(2); %TODO: fix that, extend to 2 muscles!
                            else
                                tmpResult2(stimulusIndex, iter) = relativeCommand; %TODO: fix that, extend to 2 muscles!
                            end
                            tmpResult3(stimulusIndex, iter) = model.rlModel.CCritic.v_ji * feature;

                            % total error measurement
                            testResult3(tr3Ind, iter - 1) = angleDes - angleNew;
                        end
                        tr3Ind = tr3Ind + 1;
                        tr7Ind = tr7Ind + 1;

                        % anaglyph
                        % if (abs(tmpResult1(stimulusIndex, 11)) > 3)
                        %     imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), ...
                        %             sprintf('%s/anaglyph%d_vergerr_%.2f_img%d.png', imageSavePath, tr3Ind, tmpResult1(stimulusIndex, 11), stimulusIndex));
                        % end

                        % analysis of outliers due to 'bad' stimuli
                        if vseIndex == 1        % first vergence error to be tested
                            if (angleDes - angleNew) < vseRange(vseIndex)
                                reallyBadImages(1, odIndex, stimulusIndex) = 1;
                            end
                        elseif vseIndex == 7    % last vergence error to be tested
                            if (angleDes - angleNew) > vseRange(vseIndex)
                                reallyBadImages(2, odIndex, stimulusIndex) = 1;
                            end
                        end
                    end

                    % final results
                    testResult(odIndex, vseIndex, 1 : model.testInterval + 1) = mean(tmpResult1);                         % vergErr
                    testResult(odIndex, vseIndex, model.testInterval + 2 : 2 * model.testInterval + 2) = std(tmpResult1);
                    testResult(odIndex, vseIndex, 2 * model.testInterval + 3 : 3 * model.testInterval + 3) = mean(tmpResult2);  % deltaMF
                    testResult(odIndex, vseIndex, 3 * model.testInterval + 4 : 4 * model.testInterval + 4) = std(tmpResult2);
                    testResult(odIndex, vseIndex, 4 * model.testInterval + 5 : 5 * model.testInterval + 5) = mean(tmpResult3);  % critic's response
                    testResult(odIndex, vseIndex, 5 * model.testInterval + 6 : 6 * model.testInterval + 6) = std(tmpResult3);
                end
            end
            save(strcat(imageSavePath, '/reallyBadImages'), 'reallyBadImages');
        end

        %% Reconstruction error and critic's response additional testing procedure
        if (~isempty(find(level == 2)))
            if (verbose == 1)
                sprintf('Level 2/4 Reconstruction error and critics response additional testing procedure')
            end

            recErrCritVal = zeros(1, nStim * (2 + length(model.scModel)));
            % vergence start error
            vseRange = linspace(-1, 1, test2Resolution);

            for odIndex = 1 : length(objRange)
                if (verbose == 2)
                    sprintf('Level 2/4 Test iteration = %d/%d', odIndex + size(objRange, 2), size(objRange, 2) * 2)
                end

                for vseIndex = 1 : length(vseRange)
                    if (model.getVergErrMax(objRange(odIndex)) < vseRange(vseIndex))
                        continue
                    end

                    % Calculate corresponding single muscle activity, i.e. one muscle = 0 activity
                    % [command, angleNew] = model.getMF2(objRange(odIndex), vseRange(vseIndex));

                    % Uniform muscle activation distribution for two muscles
                    [command, angleNew] = model.getMFedood(objRange(odIndex), vseRange(vseIndex));

                    % randForLeftFilt = rand(1,1);
                    % randForRightFilt = rand(1,1);

                    for stimulusIndex = 1 : nStim
                        % update stimuli
                        % currentTexture = texture{stimulusIndex};                              % stable renderer
                        % refreshImages(currentTexture, angleNew / 2, objRange(odIndex), 3);
                        currentTexture = stimulusIndex;                                         % experimental renderer
                        model.refreshImagesNew(simulator, currentTexture, angleNew / 2, objRange(odIndex), 3, [0,0,0]);

                        %% change left and right images to simulate altered rearing conditions
                        % if ~isempty(model.filterLeft)
                            % if randForLeftFilt < model.filterLeftProb
                                % model.imgGrayLeft = conv2(model.imgGrayLeft, model.filterLeft, 'same');
                                % sligthly faster version
                                % model.imgGrayLeft = conv2(model.filterLeft{1}, model.filterLeft{2}, model.imgGrayLeft, 'same');
                            % end
                        % end
                        % if ~isempty(model.filterRight)
                            % if randForRightFilt < model.filterRightProb
                                % model.imgGrayRight = conv2(model.imgGrayRight, model.filterRight, 'same');
                                % slightly faster version
                                % model.imgGrayRight = conv2(model.filterRight{1}, model.filterRight{2}, model.imgGrayRight, 'same');
                            % end
                        % end

                        % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [imageSavePath '/anaglyph.png']);
                        % generateAnaglyphs(imageSavePath, imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS);

                        % Image patch generation
                        for i = 1 : length(model.scModel)
                            model.preprocessImage(i, 1);
                            model.preprocessImage(i, 2);
                            currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
                        end

                        % Generate input feature vector from current images
                        [bfFeature, ~, recErrorArray] = model.generateFR(currentView);

                        % Absolute command feedback # concatination
                        if (model.normFeatVect == 0)
                            if (model.rlModel.continuous == 1)
                                if (model.rlModel.CActor.output_dim == 1)
                                    feature = [bfFeature; command(2) * model.lambdaMuscleFB];   % single muscle
                                else
                                    feature = [bfFeature; command * model.lambdaMuscleFB];      % two muscles
                                end
                            else
                                feature = bfFeature;
                            end
                        else
                            if (model.rlModel.continuous == 1)
                                if (model.rlModel.CActor.output_dim == 1)
                                    % z-transformed raw feature vector (no muscle feedback scaling)
                                    feature = [bfFeature; command(2)];
                                else
                                    % z-transformed raw feature vector (no muscle feedback scaling)
                                    feature = [bfFeature; command];
                                end
                            else
                                feature = bfFeature;
                            end

                            % z-transformed raw feature vector (no muscle feedback scaling)
                            for i = 1 : length(feature)
                                feature(i) = model.onlineNormalize(model.trainedUntil, feature(i), i, 0);
                            end
                            % post normalization muscle feedback scaling
                            feature = [feature(1 : end - 2); feature(end - 1 : end) * model.lambdaMuscleFB];
                        end

                        %% bias analysis
                        if (model.rlModel.bias > 0)
                            feature = [feature; model.rlModel.bias];
                        end

                        % Track reconstruction error and Critic's response
                        % recErrCritVal(stimulusIndex, :) = [model.rlModel.CCritic.v_ji * feature, sum(recErrorArray), recErrorArray];
                        recErrCritVal(1 + (stimulusIndex - 1) * (2 + length(model.scModel)) : stimulusIndex * (2 + length(model.scModel))) = ...
                            [model.rlModel.CCritic.v_ji * feature, sum(recErrorArray), recErrorArray];
                    end
                    testResult4(odIndex, vseIndex, :) = recErrCritVal;
                end
            end
            testResult4(testResult4 == 0) = NaN;
        end

        %% Object distance vs. Fixation distance
        if (~isempty(find(level == 3)))
            if (verbose == 1)
                sprintf('Level 3/4 Object distance vs. Fixation distance')
            end

            rng(42); % lets search for another seed with more big angles!
            objRange2 = [model.objDistMax, ...
                         model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, (length(testResult6) / model.testInterval) - 2), ...
                         model.objDistMin];

            % let's try this part without muscle-reset between trials!
            [command, angleNew] = model.getMFedood(objRange2(odIndex), vseRange(randi(length(vseRange))));

            tmpcnt = 1;
            for odIndex = 1 : length(testResult6) / model.testInterval
                if (verbose == 2)
                    sprintf('Level 3/4 Test iteration = %d/%d', odIndex, length(testResult6) / model.testInterval)
                end

                % vergence start error
                vergMax = model.getVergErrMax(objRange2(odIndex));
                if vergMax > 2
                    vergMax = 2;
                end
                vseRange = [linspace(-2, 0, 4), linspace(0, vergMax, 4)];
                vseRange = [vseRange(1 : 3), vseRange(5 : end)];

                % Calculate corresponding single muscle activity, i.e. one muscle = 0 activity
                % [command, angleNew] = model.getMF2(objRange2(i), vseRange(randi(length(vseRange))));

                % Uniform muscle activation distribution for two muscles
                % [command, angleNew] = model.getMFedood(objRange2(odIndex), vseRange(randi(length(vseRange))));

                currentTexture = randi(nStim);
                randForLeftFilt = rand(1,1);
                randForRightFilt = rand(1,1);

                for iter = 1 : model.testInterval
                    model.refreshImagesNew(simulator, currentTexture, angleNew / 2, objRange2(odIndex), 3, [0,0,0]);

                    %% change left and right images to simulate altered rearing conditions
                    % if ~isempty(model.filterLeft)
                        % if randForLeftFilt < model.filterLeftProb
                            % model.imgGrayLeft = conv2(model.imgGrayLeft, model.filterLeft, 'same');
                            % sligthly faster version
                            % model.imgGrayLeft = conv2(model.filterLeft{1}, model.filterLeft{2}, model.imgGrayLeft, 'same');
                        % end
                    % end
                    % if ~isempty(model.filterRight)
                        % if randForRightFilt < model.filterRightProb
                            % model.imgGrayRight = conv2(model.imgGrayRight, model.filterRight, 'same');
                            % slightly faster version
                            % model.imgGrayRight = conv2(model.filterRight{1}, model.filterRight{2}, model.imgGrayRight, 'same');
                        % end
                    % end

                    % Image patch generation
                    for i = 1 : length(model.scModel)
                        model.preprocessImage(i, 1);
                        model.preprocessImage(i, 2);
                        currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
                    end

                    % Generate input feature vector from current images
                    [bfFeature, ~, ~] = model.generateFR(currentView);

                    % Absolute command feedback # concatination
                    if (model.normFeatVect == 0)
                        if (model.rlModel.continuous == 1)
                            if (model.rlModel.CActor.output_dim == 1)
                                feature = [bfFeature; command(2) * model.lambdaMuscleFB];   % single muscle
                            else
                                feature = [bfFeature; command * model.lambdaMuscleFB];      % two muscles
                            end
                        else
                            feature = bfFeature;
                        end
                    else
                        if (model.rlModel.continuous == 1)
                            if (model.rlModel.CActor.output_dim == 1)
                                % z-transformed raw feature vector (no muscle feedback scaling)
                                feature = [bfFeature; command(2)];
                            else
                                % z-transformed raw feature vector (no muscle feedback scaling)
                                feature = [bfFeature; command];
                            end
                        else
                            feature = bfFeature;
                        end

                        % z-transformed raw feature vector (no muscle feedback scaling)
                        for i = 1 : length(feature)
                            feature(i) = model.onlineNormalize(model.trainedUntil, feature(i), i, 0);
                        end
                        % post normalization muscle feedback scaling
                        feature = [feature(1 : end - 2); feature(end - 1 : end) * model.lambdaMuscleFB];
                    end

                    %% bias analysis
                    if (model.rlModel.bias > 0)
                        feature = [feature; model.rlModel.bias];
                    end

                    %%% Action
                    relativeCommand = model.rlModel.act(feature);

                    % add the change in muscle activities to current ones
                    if (model.rlModel.continuous == 1)
                        if (model.rlModel.CActor.output_dim == 1)
                            command(2) = command(2) + relativeCommand;  % single muscle
                        else
                            command = command + relativeCommand;        % two muscles
                        end
                        command = checkCmd(command);                    % restrain motor commands to [0,1]
                        angleNew = model.getAngle(command) * 2;         % resulting angle is used for both eyes
                    else
                        angleNew = angleNew + relativeCommand;
                        if (angleNew > angleMax || angleNew < angleMin)
                            angleNew = model.vergAngleFixMin + (model.vergAngleFixMax - model.vergAngleFixMin) * rand(1, 1);
                        end
                    end
                    % testResult6(tmpcnt, :) = [objRange2(odIndex), (model.baseline / 2) / tand(angleNew / 2)];
                    testResult6(tmpcnt, :) = [atand(model.baseline / (2 * objRange2(odIndex))), angleNew / 2];
                    tmpcnt = tmpcnt + 1;
                end
            end
        end

        %% Desired vergence angle and metabolic costs approach [%] vs. iteration
        if (~isempty(find(level == 4)))
            if (verbose == 1)
                sprintf('Level 4/4 Desired vergence angle and metabolic costs approach [%] vs. iteration')
            end

            vergenceAngleApproach = zeros(size(testResult5, 1) / model.testInterval, model.testInterval);
            metCostsApproach = zeros(size(testResult5, 1) / model.testInterval, model.testInterval);
            vergAngle = zeros(1, model.testInterval);
            metCost = zeros(1, model.testInterval);

            % calculate all desired vergence angles and metabolic costs
            angleDes = objRange';
            metCostDesired = objRange';
            for odIndex = 1 : length(objRange)
                [cmdDesired, angleDes(odIndex)] = model.getMF2(objRange(odIndex), 0);
                metCostDesired(odIndex) = model.getMetCost(cmdDesired) * 2;
            end
            angleDes = repelem(angleDes, 7 * nStim, 1);
            metCostDesired = repelem(metCostDesired, 7 * nStim, 1);

            for trial = 1 : size(vergenceAngleApproach, 1)
                if (verbose == 2)
                    sprintf('Level 4/4 Test iteration = %d/%d', trial, size(vergenceAngleApproach, 1))
                end
                % starting point in muscle space
                cmdStart = [testResult5(trial * model.testInterval - model.testInterval + 1, 1) - testResult5(trial * model.testInterval - model.testInterval + 1, 3); ...
                            testResult5(trial * model.testInterval - model.testInterval + 1, 2) - testResult5(trial * model.testInterval - model.testInterval + 1, 4)];

                angleStart = model.getAngle(cmdStart) * 2;
                metCostStart = model.getMetCost(cmdStart) * 2;

                % deltaVergenceAngleDesired = vergAngleDesired - vergAngleStart
                deltaVergAnglDesired = angleDes(trial) - angleStart;

                % deltaMetabolicCostsDesired = metCostDesired - metCostStart
                deltaMetCostDesired = metCostDesired(trial) - metCostStart;

                for iter = 1 : model.testInterval
                    % vergenceAngle(iteration)
                    vergAngle(iter) = model.getAngle([testResult5(trial * model.testInterval - model.testInterval + iter, 1); ...
                                                      testResult5(trial * model.testInterval - model.testInterval + iter, 2)]) * 2;

                    % metabolicCosts(iteration)
                    metCost(iter) = model.getMetCost([testResult5(trial * model.testInterval - model.testInterval + iter, 1); ...
                                                      testResult5(trial * model.testInterval - model.testInterval + iter, 2)]) * 2;
                end

                % vergenceAngleApproach(iteration) = (vergAngle(iteration) - vergAngleStart) / deltaDesired
                vergenceAngleApproach(trial, :) = (vergAngle - angleStart) / deltaVergAnglDesired;

                % metCostsApproach(iteration) = (metCost(iteration) - metCostStart) / deltaDesired
                metCostsApproach(trial, :) = (metCost - metCostStart) / deltaMetCostDesired;
            end
            % scale to [%]
            vergenceAngleApproach = vergenceAngleApproach .* 100;
            metCostsApproach = metCostsApproach .* 100;

            % remove trials with 0° vergence start error
            memberChange = false;
            if (nStim == 0)
                nStim = 40; % catch 'just plot it' case
                memberChange = true;
            end

            odIndex2 = 0 : length(objRange) - 1;         % object distance iterator
            idxStart = (odIndex2 * 7 + 3) * nStim + 1;   % start indicies of respective trials
            idxEnd = (odIndex2 * 7 + 4) * nStim;         % end indicies of respective trials

            for k = flip([1 : length(objRange)])
                vergenceAngleApproach(idxStart(k) : idxEnd(k), :) = [];
            end

            % undo change
            if (memberChange == true)
                nStim = 0;
                memberChange = false;
            end
        end

        %% calculate response over the whole muscle plane for actor, critic and basis functions
        if (~isempty(find(level == 6)))
            if (verbose == 1)
                sprintf('Level 5/6 average responses over muscle plane')
            end
            % in the standard model, model.degreesIncRes(end,end) = model.degrees.results_deg(3,2)
            usedRows = 3;
            usedCols = 2;

            nStimuli = 5; % trade of for saving comp. time
            resolution = 3; % choose 4 for a higher resolution of muscle commands, but increased computational time

            angles = interp2(model.degrees.results_deg(1:usedRows, 1:usedCols), resolution, 'spline');
            scaleFacMR = ((usedRows - 1) / 10) / size(angles, 1);
            scaleFacLR = ((usedCols - 1) / 10) / size(angles, 2);

            if resolution == 4
                objDists = model.baseline ./ (2 * tand(angles(ceil(end/2) : 2 : end, end)));
            else
                objDists = model.baseline ./ (2 * tand(angles(ceil(end/2) : end, end)));
            end
            objDists(objDists < 0.5) = []; % erase values that are not trained

            [h, w] = size(angles);
            nMeasurements = 4;
            musclePlaneResponse = zeros(h, w, nMeasurements);

            for obj = 1 : length(objDists)
                if (verbose == 2)
                    sprintf('Level 5/5 Test iteration = %d/%d', obj, length(objDists))
                end
                objDist = objDists(obj);
                for stim = 1 : nStimuli
                    for x = 1 : h
                        for y = 1 : w
                            %% setting up the muscles and the images
                            med = x * scaleFacMR; % medial rectus activation
                            lat = y * scaleFacLR; % lateral rectus activation

                            command = [lat; med];
                            angleNew = model.getAngle(command);

                            model.refreshImagesNew(simulator, stim, angleNew, objDist, 3, [0,0,0]);

                            %% change left and right images to simulate altered rearing conditions
                            % if ~isempty(model.filterLeft)
                                % if randForLeftFilt < model.filterLeftProb
                                    % model.imgGrayLeft = conv2(model.imgGrayLeft, model.filterLeft, 'same');
                                    % sligthly faster version
                                    % model.imgGrayLeft = conv2(model.filterLeft{1}, model.filterLeft{2}, model.imgGrayLeft, 'same');
                                % end
                            % end
                            % if ~isempty(model.filterRight)
                                % if randForRightFilt < model.filterRightProb
                                    % model.imgGrayRight = conv2(model.imgGrayRight, model.filterRight, 'same');
                                    % slightly faster version
                                    % model.imgGrayRight = conv2(model.filterRight{1}, model.filterRight{2}, model.imgGrayRight, 'same');
                                % end
                            % end

                            for i = 1 : length(model.scModel)
                                model.preprocessImage(i, 1);
                                model.preprocessImage(i, 2);
                                currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
                            end

                            %% feature vector generation, edit normalization if necessary
                            [bfFeature, reward, recErrorArray] = model.generateFR(currentView);

                            if (model.normFeatVect == 0)
                                if (model.rlModel.continuous == 1)
                                    if (model.rlModel.CActor.output_dim == 1)
                                        feature = [bfFeature; command(2) * model.lambdaMuscleFB];   % single muscle
                                    else
                                        feature = [bfFeature; command * model.lambdaMuscleFB];      % two muscles
                                    end
                                else
                                    feature = bfFeature;
                                end
                            else
                                if (model.rlModel.continuous == 1)
                                    if (model.rlModel.CActor.output_dim == 1)
                                        % z-transformed raw feature vector (no muscle feedback scaling)
                                        feature = [bfFeature; command(2)];
                                    else
                                        % z-transformed raw feature vector (no muscle feedback scaling)
                                        feature = [bfFeature; command];
                                    end
                                else
                                    feature = bfFeature;
                                end

                                % z-transformed raw feature vector (no muscle feedback scaling)
                                for i = 1 : length(feature)
                                    feature(i) = model.onlineNormalize(model.trainedUntil, feature(i), i, 0);
                                end
                                % post normalization muscle feedback scaling
                                feature = [feature(1 : end - 2); feature(end - 1 : end) * model.lambdaMuscleFB];
                            end
                            % feature = [bfFeature; command];
                            % for i = 1 : length(feature)
                            %     feature(i) = model.onlineNormalize(t, feature(i), i, 1);
                            % end
                            % feature = [feature(1 : end - 2); feature(end - 1 : end) * model.lambdaMuscleFB];

                            if (model.rlModel.bias > 0)
                                feature = [feature; model.rlModel.bias];
                            end

                            %% extraction of critic and actor response
                            relativeCommand = model.rlModel.act(feature);
                            musclePlaneResponse(x, y, 1) = musclePlaneResponse(x, y, 1) + relativeCommand(1);
                            musclePlaneResponse(x, y, 2) = musclePlaneResponse(x, y, 2) + relativeCommand(2);
                            musclePlaneResponse(x, y, 3) = musclePlaneResponse(x, y, 3) + model.rlModel.CCritic.v_ji * feature;
                            musclePlaneResponse(x, y, 4) = musclePlaneResponse(x, y, 4) + reward;
                        end
                    end
                end
            end
            musclePlaneResponse = musclePlaneResponse ./ (length(objDists) * nStimuli); % forming the average
        end

        if (measureTime == true)
            elapsedTime = toc;
            sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
                    elapsedTime / 3600, elapsedTime / 60, elapsedTime, ...
                    (length(objRange) * length(vseRange) * nStim * model.testInterval + length(objRange) * length(vseRange) * nStim) / elapsedTime)
        end

        % save test results
        try
            model.testResult = testResult;
            model.testResult2 = testResult2;
            model.testResult3 = testResult3;
            model.testResult4 = testResult4;
            model.testResult5 = testResult5;
            model.testResult6 = testResult6;
            model.testResult7 = testResult7;
            model.vergenceAngleApproach = vergenceAngleApproach;
            model.metCostsApproach = metCostsApproach;
            model.musclePlaneResponse = musclePlaneResponse;
            if (saveTestResults == 1)
                model.strabAngle = strabAng;
                save(strcat(imageSavePath, '/model'), 'model');
                model.strabAngle = 0;
            end
        catch
            % catch non-existing variables error, occuring in non-up-to-date models
            try
                clone = model.copy();                                   % create deep copy & handle/pointer of model object
                delete(model);                                          % delete old object instance
                clear model;                                            % delete handle/pointer to deleted model object
                model = clone;                                          % copy new handle/pointer
                model.testResult = testResult;
                model.testResult2 = testResult2;
                model.testResult3 = testResult3;
                model.testResult4 = testResult4;
                model.testResult5 = testResult5;
                model.testResult6 = testResult6;
                model.testResult7 = testResult7;
                model.vergenceAngleApproach = vergenceAngleApproach;
                model.metCostsApproach = metCostsApproach;
                model.musclePlaneResponse = musclePlaneResponse;
                if (saveTestResults == 1)
                    model.strabAngle = strabAng;
                    save(strcat(imageSavePath, '/model'), 'model');
                    model.strabAngle = 0;
                end
                clear clone;                                            % delete obsolete handle/pointer
            catch
                % catch when new model property isn't present in Model class yet
                error('One or more new model properties (variables) are not present in Model.m class yet!');
            end
        end
    end

    sprintf('Testing procedure at iter = %s finished. Graph generation started...', folderName(9 : end))

    %%% Plotting
    if (plotIt == 1)
        % Vergence Error vs. iteration
        if (~isempty(find(level == 1)))
            rng(0);
            for odIndex = 1 : size(objRange, 2)
                figure;
                hold on;
                grid on;
                grid minor;
                for vseIndex = 1 : size(model.testResult, 2)
                    errorbar(0 : model.testInterval, squeeze(model.testResult(odIndex, vseIndex, 1 : model.testInterval + 1)), squeeze(model.testResult(odIndex, vseIndex, model.testInterval + 2 : 2 * model.testInterval + 2)), ...
                             'color', [rand, rand, rand], 'LineWidth', 1.3);
                end
                axis([-1, model.testInterval + 1, -inf, inf]);
                if (nStim > 0)
                    xlabel(sprintf('Iteration step (#stimuli=%d)', nStim), 'FontSize', 12);
                else
                    xlabel('Iteration step', 'FontSize', 12); % TODO: add nstim from model size
                end
                ylabel('Vergence Error [deg]', 'FontSize', 12);
                title(sprintf('Avg Vergence Error over Trial at Testing\nObject Distance = %.2fm', objRange(odIndex)));
                plotpath = sprintf('%s/AvgVergErrOverTrial_objDist[%.2fm].png', imageSavePath, objRange(odIndex));
                saveas(gcf, plotpath, 'png');
            end

            % 3D plot
            % figure;
            % hold on;
            % grid on;
            % [x, y] = meshgrid(0 : model.testInterval, objRange);
            % for vseIndex = 1 : size(model.testResult, 2)
            %     surf(x, y, reshape(model.testResult(:, vseIndex, 1 : 11), [4, 11]));
            % end
            % xlabel('Iteration step', 'FontSize', 12);
            % ylabel('Object distance [m]', 'FontSize', 12);
            % zlabel('Vergence Error [deg]', 'FontSize', 12);
            % title('Avg Vergence Error over Trial at Testing');
            % plotpath = sprintf('%s/AvgVergErrOverTrial3D', imageSavePath);
            % saveas(gcf, plotpath, 'png');

            % deltaMuscleForce vs Vergence Error
            figure;
            hold on;
            grid on;
            % perfect response to vergence error
            % hl1 = plot(perfectResponse(:, 1), perfectResponse(:, 2), 'color', [0.5882, 0.9608, 0], ...
            %            'DisplayName', 'perfect (fixDist_{max})', 'LineWidth', 1.3);
            % hl2 = plot(perfectResponse(:, 3), perfectResponse(:, 4), 'color', [0, 0.5882, 0.9608], ...
            %            'DisplayName', 'perfect (fixDist_{min})', 'LineWidth', 1.3);
            % lineHandles = zeros(1, 2 + length(objRange));
            % lineHandles(1 : 2) = [hl1, hl2];
            lineHandles = zeros(1, length(objRange));

            % actual response
            % xmin = 0;
            % xmax = 0;
            % for odIndex = 1 : length(objRange)
            %     % delta_mf_t+1(vergAngle_t)
            %     % hl3 = errorbar(reshape(reshape(model.testResult(odIndex, :, 1 : model.testInterval), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval]), ...
            %     %                reshape(reshape(model.testResult(odIndex, :, 24 : 33), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval]), ...
            %     %                reshape(reshape(model.testResult(odIndex, :, 35 : 44), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval]), ...
            %     %                'DisplayName', sprintf('%.2fm objDist', objRange(odIndex)), 'Marker', '*', 'MarkerSize', 2.5, ...
            %     %                'color', [rand, rand, rand], 'LineWidth', 0.7, 'LineStyle', 'none');

            %     tmpMat = sortrows([reshape(reshape(model.testResult(odIndex, :, 1 : model.testInterval), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval])', ...
            %                        reshape(reshape(model.testResult(odIndex, :, 2 * model.testInterval + 4 : 3 * model.testInterval + 3), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval])', ...
            %                        reshape(reshape(model.testResult(odIndex, :, 3 * model.testInterval + 5 : 4 * model.testInterval + 4), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval])']);

            %     [hl3, hp] = boundedline(tmpMat(:, 1), tmpMat(:, 2), tmpMat(:, 3), 'alpha');

            %     hl3.DisplayName = sprintf('%.2fm objDist', objRange(odIndex));
            %     hl3.Marker = '*';
            %     hl3.MarkerSize = 2.5;
            %     hl3.Color = [rand, rand, rand];
            %     hp.FaceColor = hl3.Color;
            %     hl3.LineStyle = 'none';
            %     % outlinebounds(hl3, hp);
            %     % lineHandles(odIndex + 2) = hl3;
            %     lineHandles(odIndex) = hl3;

            %     % for axis adjustment
            %     tmp = [min(tmpMat(:, 1)), max(tmpMat(:, 1))];
            %     if (xmin > tmp(1))
            %         xmin = tmp(1);
            %     end
            %     if (xmax < tmp(2))
            %         xmax = tmp(2);
            %     end
            % end
            % l = legend(lineHandles);
            % l.Location = 'southeast';
            % l.Box = 'off';

            % % adjust axis to actual response ranges + offset
            % ymin = -0.1;
            % ymax = 0.1;
            % plot([xmin * 1.1, xmax * 1.1], [0, 0], 'k', 'LineWidth', 0.2);
            % plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.2);
            % axis([xmin * 1.1, xmax * 1.1, ymin, ymax]);
            % xlabel(sprintf('Vergence Error [deg] (#stimuli=%d)', nStim), 'FontSize', 12);
            % ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
            % title('\Delta MF(verg_{err}) response at Testing procedure');
            % % if (~isempty(model.savePath))
            %     plotpath = sprintf('%s/deltaMFasFktVerErr', imageSavePath);
            %     saveas(gcf, plotpath, 'png');
            % % end

            % critic's response
            figure;
            hold on;
            grid on;
            grid minor;
            markers = {'p', '+', 'o', '*', '.', 'x', 's', 'd'};
            xmin = 0;
            xmax = 0;
            lineHandles = zeros(1, length(objRange));
            for odIndex = 1 : length(objRange)
                % delta_mf_t+1(vergAngle_t)
                % errorbar(reshape(reshape(model.testResult(odIndex, :, 1 : model.testInterval), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval]), ...
                %          reshape(reshape(model.testResult(odIndex, :, 46 : 55), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval]), ...
                %          reshape(reshape(model.testResult(odIndex, :, 57 : 66), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval]), ...
                %          'Marker', '*', 'MarkerSize', 2.5, 'LineWidth', 0.9, 'LineStyle', 'none');

                tmpMat = sortrows([reshape(reshape(model.testResult(odIndex, :, 1 : model.testInterval), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval])', ...
                                   reshape(reshape(model.testResult(odIndex, :, 4 * model.testInterval + 6 : 5 * model.testInterval + 5), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval])', ...
                                   reshape(reshape(model.testResult(odIndex, :, 5 * model.testInterval + 7 : 6 * model.testInterval + 6), [size(model.testResult, 2), model.testInterval])', [1, size(model.testResult, 2) * model.testInterval])']);

                [hl, hp] = boundedline(tmpMat(:, 1), tmpMat(:, 2), tmpMat(:, 3), 'alpha');

                hl.DisplayName = sprintf('%.2fm objDist', objRange(odIndex));
                if (odIndex <= length(markers))
                    hl.Marker = markers{odIndex};
                else
                    hl.Marker = markers{randi(length(markers))};
                end
                hl.MarkerSize = 2.5;
                hl.Color = [rand, rand, rand];
                hp.FaceColor = hl.Color;
                hl.LineStyle = 'none';

                lineHandles(odIndex) = hl;

                % for axis adjustment
                tmp = [min(tmpMat(:, 1)), max(tmpMat(:, 1))];
                if (xmin > tmp(1))
                    xmin = tmp(1);
                end
                if (xmax < tmp(2))
                    xmax = tmp(2);
                end
            end

            l = legend(lineHandles);
            l.Box = 'off';

            % adjust axis to actual response ranges + std deviation
            xlim([xmin * 1.1, xmax * 1.1]);
            xlabel(sprintf('Vergence Error [deg] (#stimuli=%d)', nStim), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            title('Critic Value vs. Vergence Error');
            plotpath = sprintf('%s/criticValvsVerErr', imageSavePath);
            saveas(gcf, plotpath, 'png');

            %%% Plot the resonstruction error of basis functions over different disparities
            % nBins = 1000;
            % % calculate mean and std of reconstruction error
            % tmpRsp = sortrows(model.testResult2);
            % deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nBins;
            % % recErrs = nBins x [recErr; total_mean; total_std; scale1_mean; scale1_std; ...]
            % recErrs = zeros(nBins, 1 + 2 * (length(model.scModel) + 1));
            % tmp = zeros(nBins, 3);

            % % total reconstruction error
            % for i = 1 : nBins
            %     tmp(i, 1) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
            %                                  & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 1));

            %     tmp(i, 2) = mean(sum(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
            %                                 & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2 : end), 2));

            %     tmp(i, 3) = std(sum(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
            %                                & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2 : end), 2));
            % end
            % recErrs(:, 1 : 3) = tmp;

            % reconstruction error over different scales
            % tmp = zeros(nBins, 2);
            % k = 4;
            % for i = 2 : length(model.scModel) + 1
            %     for j = 1 : nBins
            %         tmp(j, 1) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (j - 1) * deltaVergErr ...
            %                                     & tmpRsp(:, 1) <= tmpRsp(1, 1) + j * deltaVergErr), i));

            %         tmp(j, 2) = std(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (j - 1) * deltaVergErr ...
            %                                    & tmpRsp(:, 1) <= tmpRsp(1, 1) + j * deltaVergErr), i));
            %     end
            %     recErrs(:, k : k + 1) = tmp;
            %     k = k + 2;
            % end
            % recErrs(isnan(recErrs(:, 2)), :) = []; % drop NaN elements

            % figure;
            % hold on;
            % grid on;
            % grid minor;
            % handleArray = zeros(1, 1 + length(model.scModel));

            % k = 2;
            % for i = 1 : length(model.scModel) + 1
            %     handleArray(i) = errorbar(recErrs(:, 1), recErrs(:, k), recErrs(:, k + 1), 'LineWidth', 0.9);
            %     k = k + 2;
            % end

            % captions = cell(1, length(handleArray));
            % captions{1} = 'Total Error';
            % for i = 2 : length(handleArray)
            %     captions{i} = sprintf('Scale %d Error', i - 1);
            % end
            % l = legend(handleArray, captions);

            % l.FontSize = 7;
            % l.Orientation = 'horizontal';
            % l.Location = 'southoutside';
            % xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
            % ylabel('Resonstruction Error', 'FontSize', 12);
            % title(sprintf('Reconstruction Error over different disparities\nobject distances: [%s]', num2str(objRange)));

            % plotpath = sprintf('%s/recErrVsVergErr[%.2fm,%.2fm].png', imageSavePath, objRange(1), objRange(end));
            % saveas(gcf, plotpath, 'png');

            % Total error
            totelErrorHandle = figure();
            hold on;
            grid on;
            b = boxplot(model.testResult3);

            % remove outliers
            outl = findobj(b,'tag','Outliers');
            set(outl, 'Visible', 'off');

            % rescale axis to whiskers + offset
            upWi = findobj(b, 'tag', 'Upper Whisker');
            lowWi = findobj(b, 'tag', 'Lower Whisker');
            axis([0, model.testInterval + 1, ...
                  min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
                  max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

            if (nStim > 0)
                xlabel(sprintf('Iteration step (#stimuli=%d)', nStim), 'FontSize', 12);
            else
                xlabel('Iteration step', 'FontSize', 12);
            end
            ylabel('Vergence Error [deg]', 'FontSize', 12);
            title(sprintf('Total Vergence Error over Trial at Testing\nMean = %.4f°, Median = %.4f°,\n4*IQR = %.4f, RMSE = %.4f° at %dth step', ...
                          mean(model.testResult3(:, model.testInterval)), median(model.testResult3(:, model.testInterval)), ...
                          iqr(model.testResult3(:, model.testInterval)) * 4, sqrt(mean(model.testResult3(:, model.testInterval) .^ 2)), model.testInterval));

            plotpath = sprintf('%s/totalError', imageSavePath);
            saveas(gcf, plotpath, 'png');

            % remove outliers for later
            % outl = findobj(b,'tag','Outliers');
            % set(outl, 'Visible', 'off');

            %% Check for bias at 0° vergence start error
            if (nStim == 0)
                nStim = size(model.testResult3, 1) / (length(objRange) * 7);
                memberChange = true;
            end

            boxlabels = cell(1, length(objRange));
            tmp = zeros(nStim, length(objRange));
            startInd = nStim * 3 + 1;
            for i = 1 : length(objRange)
                endInd = startInd + nStim - 1;
                tmp(:, i) = model.testResult3(startInd : endInd, model.testInterval);
                startInd = endInd + nStim * 6 + 1;
                boxlabels{i} = num2str(objRange(i));
            end

            figure;
            axArray = zeros(1, 2);
            axArray(1) = subplot(1, 4, [1, 3]);
            hold on;
            grid on;
            boxplot(tmp, 'labels', boxlabels);
            xlabel('Object Distance [m]');
            ylabel('Vergence Error [deg]', 'FontSize', 12);

            axArray(2) = subplot(1, 4, 4);
            hold on;
            grid on;
            boxplot(model.testResult3(:, model.testInterval), 'widths', 0.5, 'labels', 'total');
            xlabel('\forallVerg_{err}, \forallObj_{dist}');

            for i = 1 : length(axArray)
                % remove outliers
                outl = findobj(axArray(i),'tag','Outliers');
                set(outl, 'Visible', 'off');

                % rescale axis to whiskers + offset
                upWi = findobj(axArray(i), 'tag', 'Upper Whisker');
                lowWi = findobj(axArray(i), 'tag', 'Lower Whisker');
                ylim([min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
                      max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);
            end
            % force same y-axis ranges
            linkaxes(fliplr(axArray), 'y');

            % if (nStim > 0)
            %     xlabel(sprintf('Iteration step (#stimuli=%d)', nStim), 'FontSize', 12);
            % else
            %     xlabel('Iteration step', 'FontSize', 12);
            % end
            suptitle(sprintf('Model bias at testing\n 0° vergence start error & total performance\nafter %d iterations for %d stimuli', ...
                             model.testInterval, nStim));
            plotpath = sprintf('%s/ModelBiasAt0VergErr', imageSavePath);
            saveas(gcf, plotpath, 'png');

            %%% Metabolic costs box plot
            mcDeltaHandle = figure();
            hold on;
            grid on;
            b2 = boxplot(model.testResult7);

            % remove outliers
            outl = findobj(b2,'tag','Outliers');
            set(outl, 'Visible', 'off');

            % rescale axis to whiskers + offset
            upWi = findobj(b2, 'tag', 'Upper Whisker');
            lowWi = findobj(b2, 'tag', 'Lower Whisker');
            axis([0, model.testInterval + 1, ...
                  min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
                  max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

            if (nStim > 0)
                xlabel(sprintf('Iteration step (#stimuli=%d)', nStim), 'FontSize', 12);
            else
                xlabel('Iteration step', 'FontSize', 12);
            end
            ylabel('\Deltamc = mc_{actual} - mc_{desired}', 'FontSize', 12);
            title(sprintf('Metabolic Costs over Trial at Testing (trainTime=%g)\nMean = %.4f, Median = %.4f,\n4*IQR = %.4f, RMSE = %.4f at %dth step', ...
                          model.trainedUntil, ...
                          mean(model.testResult7(:, model.testInterval)), median(model.testResult7(:, model.testInterval)), ...
                          iqr(model.testResult7(:, model.testInterval)) * 4, sqrt(mean(model.testResult7(:, model.testInterval) .^ 2)), model.testInterval));

            plotpath = sprintf('%s/totalMetCostsBox', imageSavePath);
            saveas(gcf, plotpath, 'png');

            %%% Metabolic cost delta vs. objDist
            % TODO: plot former and this graph into one subfigure
            figure;
            hold on;
            grid on;
            grid minor;
            ymin = 0;
            ymax = 0;
            lineHandles = zeros(1, length(objRange));
            nEntry = length(model.testResult7) / length(objRange); % #entries per objDist
            for odIndex = 1 : length(objRange)

                tmpMat = [[1 : model.testInterval]', ...
                          [mean(model.testResult7((odIndex - 1) * nEntry + 1 : odIndex * nEntry, :))]', ...
                          [std(model.testResult7((odIndex - 1) * nEntry + 1 : odIndex * nEntry, :))]'];

                [hl, hp] = boundedline(tmpMat(:, 1), tmpMat(:, 2), tmpMat(:, 3), 'alpha');

                hl.DisplayName = sprintf('%.2fm objDist', objRange(odIndex));

                hl.Marker = 'x';
                hl.MarkerSize = 4;

                hl.Color = [rand, rand, rand];
                hp.FaceColor = hl.Color;
                hl.LineWidth = 1.6;

                lineHandles(odIndex) = hl;

                % for axis adjustment
                tmp = [min(tmpMat(:, 2) - tmpMat(:, 3)), max(tmpMat(:, 2) + tmpMat(:, 3))];
                if (ymin > tmp(1))
                    ymin = tmp(1);
                end
                if (ymax < tmp(2))
                    ymax = tmp(2);
                end
            end
            l = legend(lineHandles);
            % l.Location = 'southeast';
            l.Box = 'off';

            % adjust axis to actual response ranges + std deviation
            axis([0, model.testInterval + 1, ymin * 1.1, ymax * 1.1]);
            if (nStim > 0)
                xlabel(sprintf('Iteration step (#stimuli=%d)', nStim), 'FontSize', 12);
            else
                xlabel('Iteration step', 'FontSize', 12);
            end
            ylabel('\Deltamc = mc_{actual} - mc_{desired}', 'FontSize', 12);
            title('\DeltaMC vs. iteration step at various objDist');
            plotpath = sprintf('%s/totalMetCostsObjDist', imageSavePath);
            saveas(gcf, plotpath, 'png');

            if (memberChange == true)
                nStim = 0;
                memberChange = false;
            end

            %%% Muscle correlation check
            % Total
            if (model.rlModel.CActor.output_dim >= 2)
                figure;
                hold on;
                % scatter(model.testResult5(:, 1), model.testResult5(:, 2), 5,'MarkerFaceColor',[0, 0.7, 0.7]);
                histHandle = hist3(model.testResult5(:, 1 : 2), [40, 40]);

                corrl = corr(model.testResult5(:, 1), model.testResult5(:, 2));
                xb = linspace(min(model.testResult5(:, 1)), max(model.testResult5(:, 1)), size(histHandle, 1));
                yb = linspace(min(model.testResult5(:, 2)), max(model.testResult5(:, 2)), size(histHandle, 1));

                pcHandle = pcolor(xb, yb, histHandle);

                axis([0, max([xb(end), 0.01]), 0, max([yb(end), 0.01])]); % take the max between those values in case xb(end) or yb(end) == 0
                shading interp;
                set(pcHandle, 'EdgeColor', 'none');

                colormap(createCM(1));
                cb = colorbar();
                cb.Label.String = '# Occurences';

                xlabel('Lateral rectus [%]', 'FontSize', 12);
                ylabel('Medial rectus [%]', 'FontSize', 12);
                title(strcat('Total Muscle Commands (testing)', sprintf('\nCorrelation = %1.2e', corrl)));
                plotpath = sprintf('%s/muscleTotalCmdTesting', imageSavePath);
                saveas(gcf, plotpath, 'png');

                % Delta
                figure;
                hold on;
                % scatter(model.testResult5(:, 3), model.testResult5(:, 4), 5,'MarkerFaceColor',[0, 0.7, 0.7]);
                histHandle = hist3(model.testResult5(:, 3 : 4), [40, 40]);

                corrl = corr(model.testResult5(:, 3), model.testResult5(:, 4));
                xb = linspace(min(model.testResult5(:, 3)), max(model.testResult5(:, 3)), size(histHandle, 1));
                yb = linspace(min(model.testResult5(:, 4)), max(model.testResult5(:, 4)), size(histHandle, 1));

                pcHandle = pcolor(xb, yb, histHandle);
                try
                    axis([xb(1), xb(end), yb(1), yb(end)]);
                catch
                    sprintf('Axis at delta muscle correlation plot have not been set.')
                end
                shading interp;
                set(pcHandle, 'EdgeColor', 'none');

                colormap(createCM(1));
                cb = colorbar();
                cb.Label.String = '# Occurences';

                xlabel('Lateral rectus [%]', 'FontSize', 12);
                ylabel('Medial rectus [%]', 'FontSize', 12);
                title(strcat('\Delta Muscle Commands (testing)', sprintf('\nCorrelation = %1.2e', corrl)));
                plotpath = sprintf('%s/muscleDeltaCmdTesting', imageSavePath);
                saveas(gcf, plotpath, 'png');
            end
        end

        % critic's response fine resolution
        % average over all objDist and all stimuli
        if (~isempty(find(level == 2)))
            vseRange = linspace(-1, 1, test2Resolution);
            figure;
            hold on;
            grid on;
            grid minor;

            tmpMatrix = permute(model.testResult4(:, :, 1 : 2 + length(model.scModel) : end), [1, 3, 2]);                 % swap 2nd and 3rd dimension
            tmpMatrix = reshape(tmpMatrix, [], size(model.testResult4(:, :, 1 : 2 + length(model.scModel) : end), 2), 1); % concatenate 3rd dimension long 1st dimension of new 2d matrix

            % handle NaNs
            tmpMean = mean(tmpMatrix, 1, 'omitnan');
            % tmpMean(isnan(tmpMean)) = [];

            vseRange = vseRange(1 : length(tmpMean));
            tmpMax = vseRange(tmpMean == max(tmpMean, [], 'omitnan'));

            tmpStd = std(tmpMatrix, 0, 1, 'omitnan');
            % tmpStd(isnan(tmpStd)) = [];

            % [hl, hp] = boundedline(vseRange, tmpMean, tmpStd, 'alpha');

            % hl.Marker = '*';
            % hl.MarkerSize = 2.5;
            % hl.Color = [rand, rand, rand];
            % hp.FaceColor = hl.Color;
            % hl.LineStyle = 'none';
            % outlinebounds(hl3, hp);

            if (nStim > 0)
                xlabel(sprintf('Vergence Error [deg] (#stimuli=%d)', nStim), 'FontSize', 12);
            else
                xlabel('Vergence Error [deg]', 'FontSize', 12);
            end
            ylabel('Value', 'FontSize', 12);
            title(sprintf('Critic Value vs. Vergence Error\nfor all objDist, max@vergErr = %.3f°', tmpMax));
            plotpath = sprintf('%s/criticValvsVerErrFineAllObjDist.png', imageSavePath);
            saveas(gcf, plotpath, 'png');

            % average over all stimuli for each objDist
            % vseRange = linspace(-1, 1, test2Resolution);
            for odIndex = 1 : length(objRange)
                figure;
                hold on;
                grid on;
                grid minor;

                % handle NaNs
                tmpMean = mean(model.testResult4(odIndex, :, 1 : 2 + length(model.scModel) : end), 3, 'omitnan');
                tmpMean(isnan(tmpMean)) = [];

                vseRange = vseRange(1 : length(tmpMean));
                tmpMax = vseRange(tmpMean == max(tmpMean));

                tmpStd = std(model.testResult4(odIndex, :, 1 : 2 + length(model.scModel) : end), 0, 3, 'omitnan');
                tmpStd(isnan(tmpStd)) = [];

                % [hl, hp] = boundedline(vseRange, ...
                %                        tmpMean, ...
                %                        tmpStd, ...
                %                        'alpha');

                % hl.Marker = '*';
                % hl.MarkerSize = 2.5;
                % hl.Color = [rand, rand, rand];
                % hp.FaceColor = hl.Color;
                % hl.LineStyle = 'none';
                % outlinebounds(hl3, hp);

                if (nStim > 0)
                    xlabel(sprintf('Vergence Error [deg] (#stimuli=%d)', nStim), 'FontSize', 12);
                else
                    xlabel('Vergence Error [deg]', 'FontSize', 12);
                end
                ylabel('Value', 'FontSize', 12);
                title(sprintf('Critic Value vs. Vergence Error\nobjDist = %.1fm, max@vergErr = %.3f°', objRange(odIndex), tmpMax));
                plotpath = sprintf('%s/criticValvsVerErrFine[%.1fm].png', imageSavePath, objRange(odIndex));
                saveas(gcf, plotpath, 'png');
            end

            % Reconstruction error fine
            vseRange = linspace(-1, 1, test2Resolution);
            handleArray = zeros(1, 1 + length(model.scModel));
            captions = cell(1, length(handleArray));
            captions{1} = 'Total Error';
            if (length(model.scModel) == 2)
                captions{2} = sprintf('Coarse Scale Error');
                captions{3} = sprintf('Fine Scale Error');
            else
                for i = 2 : length(handleArray)
                    captions{i} = sprintf('Scale %d Error', i - 1);
                end
            end

            % average over all objDist and all stimuli
            figure;
            hold on;
            grid on;
            grid minor;

            tmpMatrix = permute(model.testResult4, [1, 3, 2]);                 % swap 2nd and 3rd dimension
            tmpMatrix = reshape(tmpMatrix, [], size(model.testResult4, 2), 1); % concatenate 3rd dimension long 1st dimension of new 2d matrix

            for i = 2 : length(handleArray) + 1
                % handleArray(i - 1) = errorbar(vseRange, ...
                %                               mean(model.testResult4(odIndex, :, i : 2 + length(model.scModel) : end), 3, 'omitnan'), ...
                %                               std(model.testResult4(odIndex, :, i : 2 + length(model.scModel) : end), 0, 3, 'omitnan'), ...
                %                               'LineWidth', 0.9);
                handleArray(i - 1) = errorbar(vseRange, ...
                                              mean(tmpMatrix(i : 2 + length(model.scModel) : end, :), 1, 'omitnan'), ...
                                              std(tmpMatrix(i : 2 + length(model.scModel) : end, :), 0, 1, 'omitnan'), ...
                                              'LineWidth', 0.9);
            end

            % get min value in total reconstruction error
            tmpMean = mean(tmpMatrix(:, 2 : 2 + length(model.scModel) : end), 1, 'omitnan');
            tmpMin = vseRange(tmpMean == min(tmpMean, [], 'omitnan'));

            l = legend(handleArray, captions);
            l.FontSize = 7;
            l.Orientation = 'horizontal';
            l.Location = 'southoutside';
            xlim([vseRange(1) * 1.1, vseRange(end) * 1.1]);
            % xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
            xlabel('Vergence Error [deg]', 'FontSize', 12);
            ylabel('Resonstruction Error', 'FontSize', 12);
            title(sprintf('Reconstruction Error vs. Vergence Error\nfor all objDist, TotalMin@vergErr = %.4f°', tmpMin));

            plotpath = sprintf('%s/recErrVsVergErrFineAllObjDist.png', imageSavePath);
            saveas(gcf, plotpath, 'png');

            % average over all stimuli for each objDist
            for odIndex = 1 : length(objRange)
                figure;
                hold on;
                grid on;
                grid minor;

                for i = 2 : length(handleArray) + 1
                    handleArray(i - 1) = errorbar(vseRange, ...
                                                  mean(model.testResult4(odIndex, :, i : 2 + length(model.scModel) : end), 3, 'omitnan'), ...
                                                  std(model.testResult4(odIndex, :, i : 2 + length(model.scModel) : end), 0, 3, 'omitnan'), ...
                                                  'LineWidth', 0.9);
                end

                % get min value in total reconstruction error
                tmpMean = mean(model.testResult4(odIndex, :, 2 : 2 + length(model.scModel) : end), 3, 'omitnan');
                tmpMin = vseRange(tmpMean == min(tmpMean));

                l = legend(handleArray, captions);
                l.FontSize = 7;
                l.Orientation = 'horizontal';
                l.Location = 'southoutside';
                xlim([vseRange(1) * 1.1, vseRange(end) * 1.1]);
                % xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
                xlabel('Vergence Error [deg]', 'FontSize', 12);
                ylabel('Resonstruction Error', 'FontSize', 12);
                title(sprintf('Reconstruction Error vs. Vergence Error\nobjDist = %.1fm, TotalMin@vergErr = %.4f°', objRange(odIndex), tmpMin));

                plotpath = sprintf('%s/recErrVsVergErrFine[%.1fm].png', imageSavePath, objRange(odIndex));
                saveas(gcf, plotpath, 'png');
            end
        end

        %% Object distance vs. fixation distance: in degree (one eye)
        if (~isempty(find(level == 3)))
            figure;
            hold on;
            grid on;
            plot(model.testResult6(:, 1), 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.8);
            plot(model.testResult6(:, 2), 'color', [0, 0.6863, 1.0000], 'LineWidth', 1.3);

            xlabel(sprintf('Iteration # (interval=%d)', model.testInterval), 'FontSize', 12);
            ylabel('Angle [deg]', 'FontSize', 12);
            % ylim([0, model.objDistMax + 1]);
            ylim([0, (atand(model.baseline / (2 * model.objDistMin))) + 0.5]); % adjust to one muscle
            legend('desired (ObjDist)', 'actual (FixDist)');
            title('Vergence movements at testing');
            plotpath = sprintf('%s/fixationDistTesting_degree', imageSavePath);
            saveas(gcf, plotpath, 'png');

            %% Object distance vs. fixation distance: in degree: in meter
            figure;
            hold on;
            grid on;
            plot((model.baseline) ./ (tand(model.testResult6(:, 1)) * 2), 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.8); % adopt for two eyes
            plot((model.baseline) ./ (tand(model.testResult6(:, 2)) * 2), 'color', [0, 0.6863, 1.0000], 'LineWidth', 1.3);

            xlabel(sprintf('Iteration # (interval=%d)', model.testInterval), 'FontSize', 12);
            ylabel('Distance [m]', 'FontSize', 12);
            ylim([0, model.objDistMax + 1]);
            % ylim([0, atand(model.baseline / (2 * model.objDistMin)) + 1]);
            legend('desired (ObjDist)', 'actual (FixDist)');
            title('Vergence movements at testing');
            plotpath = sprintf('%s/fixationDistTesting_meter', imageSavePath);
            saveas(gcf, plotpath, 'png');

            %%% Final Figure
            % subplot(2, 2, 1)
        end

        %%% Desired vergence angle approach [%] vs. iteration
        if (((~isempty(find(level == 1))) && (~isempty(find(level == 4)))) ...
        ||  ((~isempty(find(level == 4))) && (~isempty(model.testResult5))))
            % Just for backward compatibility
            % TODO: remove this redundancy as soon as backward compatibility fades off
            if ((nStim == 0) && (isempty(model.vergenceAngleApproach) || isempty(model.metCostsApproach)))
                nStim = 40; % catch 'just plot it' case
                memberChange = true;

                %% Desired vergence angle and metabolic costs approach [%] vs. iteration
                vergenceAngleApproach = zeros(size(model.testResult5, 1) / model.testInterval, model.testInterval);
                metCostsApproach = zeros(size(model.testResult5, 1) / model.testInterval, model.testInterval);
                vergAngle = zeros(1, model.testInterval);
                metCost = zeros(1, model.testInterval);

                % calculate all desired vergence angles and metabolic costs
                angleDes = objRange';
                metCostDesired = objRange';
                for odIndex = 1 : length(objRange)
                    [cmdDesired, angleDes(odIndex)] = model.getMF2(objRange(odIndex), 0);
                    metCostDesired(odIndex) = model.getMetCost(cmdDesired) * 2;
                end
                angleDes = repelem(angleDes, 7 * nStim, 1);
                metCostDesired = repelem(metCostDesired, 7 * nStim, 1);

                for trial = 1 : size(vergenceAngleApproach, 1)
                    % starting point in muscle space
                    cmdStart = [model.testResult5(trial * model.testInterval - model.testInterval + 1, 1) - model.testResult5(trial * model.testInterval - model.testInterval + 1, 3); ...
                                model.testResult5(trial * model.testInterval - model.testInterval + 1, 2) - model.testResult5(trial * model.testInterval - model.testInterval + 1, 4)];

                    angleStart = model.getAngle(cmdStart) * 2;
                    metCostStart = model.getMetCost(cmdStart) * 2;

                    % deltaVergenceAngleDesired = vergAngleDesired - vergAngleStart
                    deltaVergAnglDesired = angleDes(trial) - angleStart;

                    % deltaMetabolicCostsDesired = metCostDesired - metCostStart
                    deltaMetCostDesired = metCostDesired(trial) - metCostStart;

                    for iter = 1 : model.testInterval
                        % vergenceAngle(iteration)
                        vergAngle(iter) = model.getAngle([model.testResult5(trial * model.testInterval - model.testInterval + iter, 1); ...
                                                          model.testResult5(trial * model.testInterval - model.testInterval + iter, 2)]) * 2;

                        % metabolicCosts(iteration)
                        metCost(iter) = model.getMetCost([model.testResult5(trial * model.testInterval - model.testInterval + iter, 1); ...
                                                          model.testResult5(trial * model.testInterval - model.testInterval + iter, 2)]) * 2;
                    end

                    % vergenceAngleApproach(iteration) = (vergAngle(iteration) - vergAngleStart) / deltaDesired
                    vergenceAngleApproach(trial, :) = (vergAngle - angleStart) / deltaVergAnglDesired;

                    % metCostsApproach(iteration) = (metCost(iteration) - metCostStart) / deltaDesired
                    metCostsApproach(trial, :) = (metCost - metCostStart) / deltaMetCostDesired;
                end
                % scale to [%]
                vergenceAngleApproach = vergenceAngleApproach .* 100;
                metCostsApproach = metCostsApproach .* 100;

                % remove trials with 0° vergence start error
                odIndex2 = 0 : length(objRange) - 1;         % object distance iterator
                idxStart = (odIndex2 * 7 + 3) * nStim + 1;   % start indicies of respective trials
                idxEnd = (odIndex2 * 7 + 4) * nStim;         % end indicies of respective trials

                for k = flip([1 : length(objRange)])
                    vergenceAngleApproach(idxStart(k) : idxEnd(k), :) = [];
                end

                % undo change
                if (memberChange == true)
                    nStim = 0;
                    memberChange = false;
                end

                % Backward compatibility
                if (saveTestResults == 1)
                    % catch non-existing variables error, occuring in non-up-to-date models
                    try
                        clone = model.copy();
                        delete(model);
                        clear model;
                        model = clone;
                        model.vergenceAngleApproach = vergenceAngleApproach;
                        model.metCostsApproach = metCostsApproach;
                        
                        model.strabAngle = strabAng;
                        save(strcat(imageSavePath, '/model'), 'model');
                        model.strabAngle = 0;
                        clear clone;
                    catch
                        % catch when new model property isn't present in Model class yet
                        error('One or more new model properties (variables) are not present in Model.m class yet!');
                    end
                end
            end

            % TODO: use model.vergenceAngleApproach as soon as backward compatibility fades off
            figure();
            hold on;
            grid on;
            try
                if (~isempty(model.vergenceAngleApproach))
                    b = boxplot(model.vergenceAngleApproach);
                end
            catch
                b = boxplot(vergenceAngleApproach);
            end

            % remove outliers
            outl = findobj(b,'tag','Outliers');
            set(outl, 'Visible', 'off');

            % rescale axis to whiskers + offset
            upWi = findobj(b, 'tag', 'Upper Whisker');
            lowWi = findobj(b, 'tag', 'Lower Whisker');
            try
                axis([0, model.testInterval + 1, ...
                      min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
                      max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);
            catch
                sprintf('Axis at desired verg err approach have not be set.')
            end

            xlabel('Iteration step', 'FontSize', 12);
            ylabel('Vergence Angle Approach [%]', 'FontSize', 12);
            title(sprintf('Vergence Angle Approach over Trial\nmodel trained for #%d iterations', model.trainedUntil));

            plotpath = sprintf('%s/vergenceAngleApproach', imageSavePath);
            saveas(gcf, plotpath, 'png');

            %%% Desired metabolic costs approach [%] vs. iteration
            % TODO: use model.metCostsApproach as soon as backward compatibility fades off
            figure();
            hold on;
            grid on;

            try
                if (~isempty(model.metCostsApproach))
                    b = boxplot(model.metCostsApproach);
                end
            catch
                b = boxplot(metCostsApproach);
            end

            % remove outliers
            outl = findobj(b,'tag','Outliers');
            set(outl, 'Visible', 'off');

            % rescale axis to whiskers + offset
            upWi = findobj(b, 'tag', 'Upper Whisker');
            lowWi = findobj(b, 'tag', 'Lower Whisker');
            try
                axis([0, model.testInterval + 1, ...
                      min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
                      max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);
            catch
                sprintf('Axis at desired met costs approach have not be set.')
            end

            xlabel('Iteration step', 'FontSize', 12);
            ylabel('Metabolic Costs Approach [%]', 'FontSize', 12);
            title(sprintf('Metabolic Costs Approach over Trial\nmodel trained for #%d iterations', model.trainedUntil));

            plotpath = sprintf('%s/metCostsApproach', imageSavePath);
            saveas(gcf, plotpath, 'png');
        end

        %% Generate muscle activation trajectories
        if (~isempty(find(level == 5)))
            if verbose
                sprintf('generating trajectory ...')
            end
            model.plotTrajectory([0.5, 6], [-2, 0, 2], 'advanced', 200, randi(max(nStim, 40)), simulator, imageSavePath, folderName(9 : end), plotIt);
        end

        %% plot responses over muscle plane
        if (~isempty(find(level == 6)))
            usedRows = 3;
            usedCols = 2;
            scaleFacMR = ((usedRows - 1) / 10) / size(model.musclePlaneResponse, 1);    % table scaling factors for backwards
            scaleFacLR = ((usedCols - 1) / 10) / size(model.musclePlaneResponse, 2);

            % zRange = [-0.5, 0];
            [h, w, ~] = size(model.musclePlaneResponse);
            [x, y] = meshgrid(1 : h, 1:w);
            fig1 = figure;

            s1 = subplot(2, 2, 1); % critic sideView
            % s1 = subplot(2, 2, [1, 3]);
            surf(model.musclePlaneResponse(:, :, 3));
            % title(sprintf('Critics Response at %d iterations', model.trainedUntil));
            title('critic response');
            view([-37.5, 30]);
            axis 'tight';
            % zlim(zRange);
            zlabel('value');
            xLabels = get(s1, 'XTickLabel');
            yLabels = get(s1, 'YTickLabel');
            set(s1, 'XTickLabel', num2str(linspace(0, 0.1, length(xLabels))', '%0.2f'));
            set(s1, 'YTickLabel', num2str(linspace(0, 0.2, length(yLabels))', '%0.2f'));
            xlabel('lr [%]');
            ylabel('mr [%]');

            s2 = subplot(2, 2, 2); % critic top view
            pcolor(model.musclePlaneResponse(:, :, 3));
            hold on;
            title('critic response');
            % view(2);
            axis 'tight';
            xLabels = get(s2, 'XTickLabel');
            yLabels = get(s2, 'YTickLabel');
            set(s2, 'XTickLabel', num2str(linspace(0, 0.1, length(xLabels))', '%0.2f'));
            set(s2, 'YTickLabel', num2str(linspace(0, 0.2, length(yLabels))', '%0.2f'));
            xlabel('lr [%]');
            ylabel('mr [%]');
            colorbar();

            objDists = [7.6434, 2.2482, 1.3110, 0.9218, 0.7087, 0.5742];

            for od = 1 : length(objDists)
                objDist = objDists(od);
                [lr, mr] = model.getAnglePoints(objDist, 0);
                xpos = (lr * (size(model.musclePlaneResponse, 2) - 1) * (10 / (usedCols - 1))) + 1; % convert muscle activities to table entries
                ypos = (mr * (size(model.musclePlaneResponse, 1) - 1) * (10 / (usedRows - 1))) + 1;
                plot(xpos, ypos, 'color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 0.1); %[0, 0.5882, 0]
            end

            s3 = subplot(2, 2, 3);
            surf(model.musclePlaneResponse(:, :, 4));
            % title(sprintf('Critics Response at %d iterations', model.trainedUntil));
            title('reconstruction errror');
            view([-37.5, 30]);
            axis 'tight';
            % zlim(zRange);
            zlabel('value');
            xLabels = get(s3, 'XTickLabel');
            yLabels = get(s3, 'YTickLabel');
            set(s3, 'XTickLabel', num2str(linspace(0, 0.1, length(xLabels))', '%0.2f'));
            set(s3, 'YTickLabel', num2str(linspace(0, 0.2, length(yLabels))', '%0.2f'));
            xlabel('lr [%]');
            ylabel('mr [%]');

            s4 = subplot(2, 2, 4); % basis function response
            pcolor(model.musclePlaneResponse(:, :, 4));
            hold on;
            title('reconstruction errror');
            view(2);
            axis 'tight';
            xLabels = get(s4, 'XTickLabel');
            yLabels = get(s4, 'YTickLabel');
            set(s4, 'XTickLabel', num2str(linspace(0, 0.1, length(xLabels))', '%0.2f'));
            set(s4, 'YTickLabel', num2str(linspace(0, 0.2, length(yLabels))', '%0.2f'));
            xlabel('lr [%]');
            ylabel('mr [%]');
            colorbar();

            objDists = [7.6434, 2.2482, 1.3110, 0.9218, 0.7087, 0.5742];

            for od = 1 : length(objDists)
                objDist = objDists(od);
                [lr, mr] = model.getAnglePoints(objDist, 0);
                xpos = (lr * (size(model.musclePlaneResponse, 2) - 1) * (10 / (usedCols - 1))) + 1; % convert muscle activities to table entries
                ypos = (mr * (size(model.musclePlaneResponse, 1) - 1) * (10 / (usedRows - 1))) + 1;
                plot(xpos, ypos, 'color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', .5); %[0, 0.5882, 0]
            end


            suptitle(sprintf('%d Iterations of Training', model.trainedUntil));
            % saveas(fig1, sprintf('%s/overviewAt%06d.png', folder, index));
            saveas(fig1, sprintf('%s/overviewAt%07d.png', imageSavePath, model.trainedUntil));

            close(fig1);
        end
    end

    %%% Results overview table generation
    resultsFN = strcat(model.savePath, '/results.ods'); % file name
    resultsFID = fopen(resultsFN, 'a');                 % file descriptor

    resultsOverview = cell(1);                          % results value vector
    resultsOverview{end} = '';

    % traintime
    if (model.trainTime >= 1e6)
        resultsOverview{end + 1} = strcat(num2str(model.trainedUntil / 1e6), 'mio');
    elseif (model.trainTime >= 1e3)
        resultsOverview{end + 1} = strcat(num2str(model.trainedUntil / 1e3), 'k');
    else
        resultsOverview{end + 1} = num2str(model.trainedUntil);
    end
    formatSpec = '%s, %s,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % nBasis
    resultsOverview = {};
    formatSpec = {''};
    for k = 1 : length(model.dsRatio)
        resultsOverview{k} = model.scModel{k}.nBasis;
        formatSpec = strcat(formatSpec, {'%.0f '});
    end
    formatSpec = strcat(formatSpec, ',');
    fprintf(resultsFID, formatSpec{1 : end}, resultsOverview{1 : end});

    % lambdas
    resultsOverview = {model.lambdaRec, model.lambdaMuscleFB, model.lambdaMet, (model.lambdaMet / 0.1622) * 100};
    formatSpec = '%.0f, %.4f, %.4f, %.1f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % field of view
    resultsOverview = {model.pxFieldOfView .* model.dsRatio};
    formatSpec = {''};
    for k = 1 : length(model.dsRatio)
        formatSpec = strcat(formatSpec, {'%.0f '});
    end
    formatSpec = strcat(formatSpec, ',');
    fprintf(resultsFID, formatSpec{1 : end}, resultsOverview{1 : end});

    % dsRatio
    resultsOverview = {model.dsRatio};
    fprintf(resultsFID, formatSpec{1 : end}, resultsOverview{1 : end});

    % stride
    resultsOverview = {model.stride / model.patchSize};
    formatSpec = strrep(formatSpec, '0', '1');
    fprintf(resultsFID, formatSpec{1 : end}, resultsOverview{1 : end});

    % learning rates
    resultsOverview = {model.rlModel.CCritic.alpha_v, model.rlModel.CActor.beta_p, model.scModel{1}.eta};
    formatSpec = '%.2f, %.2f, %.2f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % weight init
    resultsOverview = {model.rlModel.weight_range .* ...
                       [model.rlModel.CCritic.input_dim, ...
                       (model.rlModel.CActor.input_dim * model.rlModel.CActor.hidden_dim), ...
                       (model.rlModel.CActor.hidden_dim * model.rlModel.CActor.output_dim)]};
    formatSpec = '%.0f %.0f %.0f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % regularizer
    resultsOverview = {model.rlModel.CActor.regularizer};
    formatSpec = '%.4f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % column fill
    resultsOverview = {''};
    formatSpec = '%s,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % results vergErr
    resultsOverview = {sqrt(mean(model.testResult3(:, model.testInterval) .^ 2)), median(model.testResult3(:, model.testInterval)), iqr(model.testResult3(:, model.testInterval)) * 4};
    formatSpec = '%.4f, %.4f, %.4f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % column fill
    resultsOverview = {''};
    formatSpec = '%s,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % results delta metCost
    resultsOverview = {sqrt(mean(model.testResult7(:, model.testInterval) .^ 2)), median(model.testResult7(:, model.testInterval)), iqr(model.testResult7(:, model.testInterval)) * 4};
    formatSpec = '%.4f, %.4f, %.4f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % directory
    resultsOverview = {'', '', '', GetFullPath(imageSavePath)};
    formatSpec = '%s, %s, %s, %s\n';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % close file
    fclose(resultsFID);

    % save testing performance history over traintime
    try
        model.testHist(find(model.testAt == str2num(folderName(9 : end))), :) = ...
                      [sqrt(mean(model.testResult3(:, model.testInterval) .^ 2)), ...
                       mean(abs(model.testResult3(:, model.testInterval) .^ 2)), ...
                       std(abs(model.testResult3(:, model.testInterval) .^ 2)), ...
                       sqrt(mean(model.testResult7(:, model.testInterval) .^ 2)), ...
                       mean(abs(model.testResult7(:, model.testInterval) .^ 2)), ...
                       std(abs(model.testResult7(:, model.testInterval) .^ 2))];

        if (saveTestResults == 1)
            model.strabAngle = strabAng;
            save(strcat(imageSavePath, '/model'), 'model');
            model.strabAngle = 0;
        end
    catch
        warning('Current model has no \"testHist\" field. Performance history will not be stored.');
    end
    if (saveTestResults == 1)
        model.strabAngle = strabAng;
        save(strcat(imageSavePath, '/model'), 'model');
        model.strabAngle = 0;
    end
    sprintf('Testing procedure at iter = %s finished. Graph generation finished. Model saved.', folderName(9 : end))
end

%this function generates anaglyphs of the large and small scale fovea and
%one of the two unpreprocessed gray scale images
function generateAnaglyphs(imageSavePath, leftGray, rightGray, dsRatioL, dsRatioS, foveaL, foveaS)
    anaglyph = imfuse(leftGray, rightGray, 'falsecolor');
    imwrite(anaglyph,  sprintf('%s/anaglyph.png', imageSavePath));

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
    imwrite(imresize(anaglyphL, 20), sprintf('%s/anaglyphLargeScale.png', imageSavePath));
    largeScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(largeScaleView, 20), sprintf('%s/LargeScaleMontage.png', imageSavePath));

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
    imwrite(imresize(anaglyphS, 8), sprintf('%s/anaglyphSmallScale.png', imageSavePath));
    smallScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(smallScaleView, 8), sprintf('%s/smallScaleMontage.png', imageSavePath));
end

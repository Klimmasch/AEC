%%% Model testing procedure
%@param model               respective model object to be tested
%@param nStim               # stimuli to be tested, provide 0 if you just want to plot results
%@pram plotIt               whether plots shall be generated
%@param saveTestResults     whether to save the results (not recommended if model is still trained!)
%@param simulator           simulator handle, provide [] if there is no simulator handle yet
%@param reinitRenderer      1 if renderer was already initialized
%                           0 if renderer wasn't initialized yet
%%%
function testModelContinuous(model, nStim, plotIt, saveTestResults, simulator, reinitRenderer, folder)

    % should the simulation time be measured?
    % not recommended for use during training

    measureTime = false;

    % Results overview table generation
    resultsFN = strcat(model.savePath, '/results.ods'); % file name
    resultsFID = fopen(resultsFN, 'a');                 % file descriptor

    resultsOverview = cell(1);                          % results value vector
    resultsOverview{end} = '';

    if (model.trainTime >= 1e6)
        resultsOverview{end + 1} = strcat(num2str(model.trainedUntil / 1e6), 'mio');
    elseif (model.trainTime >= 1e3)
        resultsOverview{end + 1} = strcat(num2str(model.trainedUntil / 1e3), 'k');
    else
        resultsOverview{end + 1} = num2str(model.trainedUntil);
    end
    formatSpec = '%s, %s,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    resultsOverview = {model.lambdaRec, model.lambdaMuscleFB, model.lambdaMet};
    formatSpec = '%.0f, %.4f, %.4f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    resultsOverview = {model.pxFieldOfView .* model.dsRatio};
    formatSpec = {''};
    for k = 1 : length(model.dsRatio)
        formatSpec = strcat(formatSpec, {'%.0f '});
    end
    formatSpec = strcat(formatSpec, ',');
    fprintf(resultsFID, formatSpec{1 : end}, resultsOverview{1 : end});

    resultsOverview = {model.dsRatio};
    fprintf(resultsFID, formatSpec{1 : end}, resultsOverview{1 : end});

    resultsOverview = {model.stride / model.patchSize};
    formatSpec = strrep(formatSpec, '0', '1');
    fprintf(resultsFID, formatSpec{1 : end}, resultsOverview{1 : end});

    resultsOverview = {model.rlModel.CCritic.alpha_v, model.rlModel.CActor.beta_p, model.scModel{1}.eta};
    formatSpec = '%.2f, %.2f, %.2f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    resultsOverview = {model.rlModel.weight_range .* ...
                       [model.rlModel.CCritic.input_dim, ...
                       (model.rlModel.CActor.input_dim * model.rlModel.CActor.hidden_dim), ...
                       (model.rlModel.CActor.hidden_dim * model.rlModel.CActor.output_dim)]};
    formatSpec = '%.0f %.0f %.0f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    resultsOverview = {model.rlModel.CActor.regularizer};
    formatSpec = '%.4f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    resultsOverview = {model.muscleInitMax};
    formatSpec = '%.4f %.4f,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % column fill
    resultsOverview = {'', '', '', ''};
    formatSpec = '%s, %s, %s, %s,';
    fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

    % Vergence error resolution for 2nd testing procedure
    % Needs to be odd number to include vergErr = 0
    test2Resolution = 101;

    %%% Stimulus declaration
    % textureFile = 'Textures_mcgillManMadeTrain(jpg).mat';     % McGill man made database
    textureFile = 'Textures_mcgillManMadeTest(jpg).mat';
    % textureFile = 'Textures_mcgillFruitsAll(jpg).mat';        % McGill fruits database
    % textureFile = 'Textures_mcgillFoliageTrain(jpg).mat';     % McGill foliage database
    % textureFile = 'Textures_mcgillFoliageTest(jpg).mat';
    % textureFile = 'Textures_vanHaterenTrain.mat';             % vanHateren database
    % textureFile = 'Textures_vanHaterenTest.mat';
    % textureFile = 'Textures_celine.mat';                      % Celine's images

    % Prepare Textures
    texture = load(sprintf('config/%s', textureFile));
    texture = texture.texture;
    nTextures = length(texture);
    if (nTextures < nStim)
        sprintf('The texture file only contains %d images, but I will use them all!', nTextures)
        nStim = nTextures;
    end

    % cancel testing procedure
    if ((nStim == 0) && (isempty(model.testResult)))
        sprintf('Error: Model has no testResults!')
        return;
    elseif (nStim == 1)
        sprintf('Error: nStim must be != 1')
        return;
    end

    %%% New renderer
    if (isempty(simulator))
        % simulator = OpenEyeSim('create'); % stable renderer
        simulator = OpenEyeSimV5('create'); % experimental renderer

        if (reinitRenderer == 0)
            simulator.initRenderer();
        else
            % for debugging purposes
            simulator.reinitRenderer();
        end
    end

    % load all stimuli into memory for experimental renderer
    for i = 1 : nStim
        simulator.add_texture(i, texture{i});
    end

    %%% creating a new directory if (folder ~= '/.')
    if (folder(1) ~= '/')
        folder = ['/' folder];
    end
    imageSavePath = [model.savePath folder];
    mkdir(imageSavePath);

    % fixation interval at testing procedure
    testInterval = model.interval * 2;

    command = [0; 0];
    objRange = [model.objDistMin : 0.5 : model.objDistMax];

    tmpResult1 = zeros(nStim, testInterval + 1);
    tmpResult2 = zeros(nStim, testInterval + 1);
    tmpResult3 = zeros(nStim, testInterval + 1);

    testResult = zeros(length(objRange), 7, 66); % mean and std of vergenceError, deltaMF and critics response for different obj dists and starting pos
    testResult2 = zeros(length(objRange) * 7 * nStim * testInterval, 1 + length(model.scModel)); % reconstruction error statistics

    testResult3 = zeros(length(objRange) * 7 * nStim, testInterval); % ALL single values
    testResult4 = zeros(length(objRange), test2Resolution, nStim * (2 + length(model.scModel)));
    testResult5 = zeros(length(objRange) * 7 * nStim * testInterval, model.rlModel.CActor.output_dim * 2); % correlation between abs muscle activations and deltaMFs

    realyBadImages = zeros(2, length(objRange), nStim); % here, the images are safed that start at the maximal vergence errors (pos & neg) and that end up worse than they started
    %this tabular is going to be safed inside the models folder and
    %histograms will be generated

    degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
    % metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

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
    indZero = find(mfunction(:, 2) == 0);   % MF == 0_index

    approx = spline(1 : 11, degrees.results_deg(1, :));
    xValPos = ppval(approx, 1 : 0.0001 : 11)';
    yValPos = linspace(0, 1, resolution)';

    % xValNeg = flipud(ppval(approx, 1 : 0.0001 : 11)' * -1);
    % yValNeg = linspace(-1, 0, resolution)';

    mfunction2 = [xValPos, yValPos];
    mfunction2(:, 1) = mfunction2(:, 1) * 2;    % angle for two eyes
    dmf3 = abs(diff(mfunction2(1 : 2, 1)));     % delta in angle
    dmf4 = diff(mfunction2(1 : 2, 2));          % delta in mf
    indZero = find(mfunction2(:, 2) == 0);      % MF == 0_index

    %%% Perfect Response function
    % indMaxFix = find(mfunction(:, 1) <= model.vergAngleFixMin + dmf & mfunction(:, 1) >= model.vergAngleFixMin - dmf); % MF(vergAngleFixMin)_index
    % indMaxFix = indMaxFix(1);
    % indMinFix = find(mfunction(:, 1) <= model.vergAngleFixMax + dmf & mfunction(:, 1) >= model.vergAngleFixMax - dmf); % MF(vergAngleFixMax)_index
    % indMinFix = indMinFix(1);

    % perfect_response := [max_fixation_x, max_fixation_y, min_fixation_x, min_fixation_y]
    % x = vergenceError, y = deltaMuscelForce
    % perfectResponseMaxFix = [(mfunction(indMaxFix, 1) - flipud(mfunction(indMaxFix : end, 1))), ...
    %                          (mfunction(indMaxFix, 2) - flipud(mfunction(indMaxFix : end, 2))); ...
    %                          (mfunction(indMaxFix, 1) - flipud(mfunction(indZero : indMaxFix - 1, 1))), ...
    %                          (mfunction(indMaxFix, 2) - flipud(mfunction(indZero : indMaxFix - 1, 2)))];

    % perfectResponseMinFix = [(mfunction(indMinFix, 1) - flipud(mfunction(indMinFix : end, 1))), ...
    %                          (mfunction(indMinFix, 2) - flipud(mfunction(indMinFix : end, 2))); ...
    %                          (mfunction(indMinFix, 1) - flipud(mfunction(indZero : indMinFix - 1, 1))), ...
    %                          (mfunction(indMinFix, 2) - flipud(mfunction(indZero : indMinFix - 1, 2)))];

    % perfectResponse = [perfectResponseMaxFix, perfectResponseMinFix];

    % minimal and maximal angle that can be reached by one-dimensional muscle commands
    angleMin = min(mfunction2(mfunction2(:, 1) > 0));
    angleMax = mfunction(end, 1);

    % Color images for left & right eye
    imgRawLeft = uint8(zeros(240, 320, 3));
    imgRawRight = uint8(zeros(240, 320, 3));
    imgGrayLeft = uint8(zeros(240, 320, 3));
    imgGrayRight = uint8(zeros(240, 320, 3));

    % Image patches cell array (input to model)
    currentView = cell(1, length(model.scModel));

    %%% Helper function that maps {objDist, desiredVergErr} -> {muscleForce, angleInit}
    function [mf, angleInit] = getMF(objDist, desVergErr)
        % correct vergence angle for given object distance
        angleCorrect = 2 * atand(model.baseline / (2 * objDist));
        % desired init angle for given vergence error [deg]
        angleInit = angleCorrect - desVergErr;
        % look up index of angleInit
        indAngleInit = find(mfunction(:, 1) <= angleInit + dmf & mfunction(:, 1) >= angleInit - dmf);
        mf = mfunction(indAngleInit, 2);
        mf = mf(ceil(length(mf) / 2));
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

    % %%% Helper function for calculating {objDist} -> {minVergErr, maxVergErr}
    % function [vergErrMin, vergErrMax] = getVergErrMinMax(objDist)
    %     % correct vergence angle for given object distance
    %     angleCorrect = 2 * atand(model.baseline / (2 * objDist));
    %     vergErrMin = angleCorrect - angleMax;
    %     vergErrMax = angleCorrect - angleMin;
    % end

    %%% Helper function for calculating {objDist} -> {maxVergErr}
    function vergErrMax = getVergErrMax(objDist)
        % correct vergence angle for given object distance
        angleCorrect = 2 * atand(model.baseline / (2 * objDist));
        vergErrMax = angleCorrect - angleMin;
    end

    %%% Helper function that maps muscle activities to resulting angle
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
    % function tmpMetCost = getMetCost(command)
    %     cmd = (command * 10) + 1;                                           % scale commands to table entries
    %     tmpMetCost = interp2(metCosts.results, cmd(1), cmd(2), 'spline');   % interpolate in table by cubic splines
    % end

    %%% Generates two new images for both eyes
    % texture:      file path of texture input
    % eyeAngle:     angle of single eye (rotation from offspring)
    % objDist:      distance of stimulus
    % scaleImSize:  scaling factor of stimulus plane [m]
    function refreshImages(texture, eyeAngle, objDist, scalingFactor)
        simulator.add_texture(1, texture);
        simulator.set_params(1, eyeAngle, objDist, 0, scalingFactor);

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
    if (measureTime == true)
        tic;
    end
    % don't repeat testing procedure if nStim == 0, but just plot the results
    if (nStim > 0)
        for odIndex = 1 : length(objRange)
            sprintf('Test iteration = %d/%d', odIndex, size(objRange, 2) * 2)

            % vergence start error
            vergMax = getVergErrMax(objRange(odIndex));
            if vergMax > 2
                vergMax = 2;
            end
            vseRange = [linspace(-2, 0, 4), linspace(0, vergMax, 4)];
            vseRange = [vseRange(1 : 3), vseRange(5 : end)];
            % vseRange = [-3:3];
            % vseRange = linspace(-1, 1, 7);
            angleDes = 2 * atand(model.baseline / (2 * objRange(odIndex)));

            for vseIndex = 1 : length(vseRange)
                tmpResult1(:, 1) = vseRange(vseIndex);

                for stimulusIndex = 1 : nStim
                    % currentTexture = texture{stimulusIndex};  % stable renderer
                    currentTexture = stimulusIndex;             % experimental renderer

                    % command(1) = 0;
                    % [command(2), angleNew] = getMF(objRange(odIndex), vseRange(vseIndex));
                    [command, angleNew] = getMF2(objRange(odIndex), vseRange(vseIndex));

                    for iter = 2 : testInterval + 1
                        % update stimuli
                        % refreshImages(currentTexture, angleNew / 2, objRange(odIndex), 3);    % stable renderer
                        refreshImagesNew(currentTexture, angleNew / 2, objRange(odIndex), 3);   % experimental renderer

                        % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [imageSavePath '/anaglyph.png']);
                        % generateAnaglyphs(imageSavePath, imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS);

                        % Image patch generation
                        for i = 1 : length(model.scModel)
                            model.preprocessImage(imgGrayLeft, i, 1);
                            model.preprocessImage(imgGrayRight, i, 2);
                            currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
                        end

                        % Generate input feature vector from current images
                        [bfFeature, ~, recErrorArray] = model.generateFR(currentView);

                        % Track reconstruction error statistics
                        testResult2(tr2Ind, :) = [angleDes - angleNew, recErrorArray];
                        tr2Ind = tr2Ind + 1;

                        % Absolute command feedback # concatination
                        if (model.rlModel.continuous == 1)
                            if (model.rlModel.CActor.output_dim == 1)
                                feature = [bfFeature; command(2) * model.lambdaMuscleFB];   % single muscle
                            else
                                feature = [bfFeature; command * model.lambdaMuscleFB];      % two muscles
                            end
                        else
                            feature = bfFeature;
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
                            angleNew = getAngle(command) * 2;               % resulting angle is used for both eyes
                        else
                            angleNew = angleNew + relativeCommand;
                            if (angleNew > angleMax || angleNew < angleMin)
                                angleNew = model.vergAngleFixMin + (model.vergAngleFixMax - model.vergAngleFixMin) * rand(1, 1);
                            end
                        end

                        % track commands for correlation plot
                        if (model.rlModel.CActor.output_dim >= 2)
                            testResult5(tr5Ind, :) = [command; relativeCommand];
                            tr5Ind = tr5Ind + 1;
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

                    % anaglyph
                    % if (abs(tmpResult1(stimulusIndex, 11)) > 3)
                    %     imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), ...
                    %             sprintf('%s/anaglyph%d_vergerr_%.2f_img%d.png', imageSavePath, tr3Ind, tmpResult1(stimulusIndex, 11), stimulusIndex));
                    % end

                    if vseIndex == 1 %first vergence error to be tested
                        if (angleDes - angleNew) < vseRange(vseIndex)
                            realyBadImages(1, odIndex, stimulusIndex) = 1;
                        end
                    elseif vseIndex == 7 %last vergence error to be tested
                        if (angleDes - angleNew) > vseRange(vseIndex)
                            realyBadImages(2, odIndex, stimulusIndex) = 1;
                        end
                    end
                end

                % final results

                testResult(odIndex, vseIndex, 1 : testInterval + 1) = mean(tmpResult1);                         % vergErr
                testResult(odIndex, vseIndex, testInterval + 2 : 2 * testInterval + 2) = std(tmpResult1);
                testResult(odIndex, vseIndex, 2 * testInterval + 3 : 3 * testInterval + 3) = mean(tmpResult2);  % deltaMF
                testResult(odIndex, vseIndex, 3 * testInterval + 4 : 4 * testInterval + 4) = std(tmpResult2);
                testResult(odIndex, vseIndex, 4 * testInterval + 5 : 5 * testInterval + 5) = mean(tmpResult3);  % critic's response
                testResult(odIndex, vseIndex, 5 * testInterval + 6 : 6 * testInterval + 6) = std(tmpResult3);
                % testResult(odIndex, vseIndex, 1 : 11) = mean(tmpResult1);   % vergErr
                % testResult(odIndex, vseIndex, 12 : 22) = std(tmpResult1);
                % testResult(odIndex, vseIndex, 23 : 33) = mean(tmpResult2);  % deltaMF
                % testResult(odIndex, vseIndex, 34 : 44) = std(tmpResult2);
                % testResult(odIndex, vseIndex, 45 : 55) = mean(tmpResult3);  % critic's response
                % testResult(odIndex, vseIndex, 56 : 66) = std(tmpResult3);
            end
        end
        save(strcat(imageSavePath, '/realyBadImages'), 'realyBadImages');

        %% Reconstruction error and critic's response additional testing procedure
        tmp = zeros(1, nStim * (2 + length(model.scModel)));
        % vergence start error
        vseRange = linspace(-1, 1, test2Resolution);

        for odIndex = 1 : length(objRange)
            sprintf('Test iteration = %d/%d', odIndex + size(objRange, 2), size(objRange, 2) * 2)

            for vseIndex = 1 : length(vseRange)
                if (getVergErrMax(objRange(odIndex)) < vseRange(vseIndex))
                    continue
                end
                % command(1) = 0;
                % [command(2), angleNew] = getMF(objRange(odIndex), vseRange(vseIndex));
                [command, angleNew] = getMF2(objRange(odIndex), vseRange(vseIndex));
                for stimulusIndex = 1 : nStim
                    % update stimuli
                    % currentTexture = texture{stimulusIndex};                              % stable renderer
                    % refreshImages(currentTexture, angleNew / 2, objRange(odIndex), 3);
                    currentTexture = stimulusIndex;                                         % experimental renderer
                    refreshImagesNew(currentTexture, angleNew / 2, objRange(odIndex), 3);

                    % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [imageSavePath '/anaglyph.png']);
                    % generateAnaglyphs(imageSavePath, imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS);

                    % Image patch generation
                    for i = 1 : length(model.scModel)
                        model.preprocessImage(imgGrayLeft, i, 1);
                        model.preprocessImage(imgGrayRight, i, 2);
                        currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
                    end

                    % Generate input feature vector from current images
                    [bfFeature, ~, recErrorArray] = model.generateFR(currentView);

                    % Absolute command feedback # concatination
                    if (model.rlModel.continuous == 1)
                        if (model.rlModel.CActor.output_dim == 1)
                            feature = [bfFeature; command(2) * model.lambdaMuscleFB];   % single muscle
                        else
                            feature = [bfFeature; command * model.lambdaMuscleFB];      % two muscles
                        end
                    else
                        feature = bfFeature;
                    end

                    % Track reconstruction error and Critic's response
                    % tmp(stimulusIndex, :) = [model.rlModel.CCritic.v_ji * feature, sum(recErrorArray), recErrorArray];
                    tmp(1 + (stimulusIndex - 1) * (2 + length(model.scModel)) : stimulusIndex * (2 + length(model.scModel))) = ...
                        [model.rlModel.CCritic.v_ji * feature, sum(recErrorArray), recErrorArray];
                end
                % testResult4(odIndex, vseIndex, :) = mean(tmp);
                % testResult4(odIndex, vseIndex, :) = reshape(tmp', [1, size(tmp, 1) * size(tmp, 2)]);
                testResult4(odIndex, vseIndex, :) = tmp;
            end
        end
        testResult4(testResult4 == 0) = NaN;

        if (measureTime == true)
            elapsedTime = toc;
            sprintf('Time = %.2f [h] = %.2f [min] = %f [sec]\nFrequency = %.4f [iterations/sec]', ...
                    elapsedTime / 3600, elapsedTime / 60, elapsedTime, ...
                    (length(objRange) * length(vseRange) * nStim * testInterval + length(objRange) * length(vseRange) * nStim) / elapsedTime)
        end

        % save test results
        try
            model.testResult = testResult;
            model.testResult2 = testResult2;
            model.testResult3 = testResult3;
            model.testResult4 = testResult4;
            model.testResult5 = testResult5;
            if (saveTestResults == 1)
                save(strcat(imageSavePath, '/model'), 'model');
            end
        catch
            % catch non-existing variables error, occuring in non-up-to-date models
            try
                clone = model.copy();
                delete(model);
                clear model;
                model = clone;
                model.testResult = testResult;
                model.testResult2 = testResult2;
                model.testResult3 = testResult3;
                model.testResult4 = testResult4;
                model.testResult5 = testResult5;
                if (saveTestResults == 1)
                    save(strcat(imageSavePath, '/model'), 'model');
                end
                delete(clone);
                clear clone;
            catch
                % catch when new model property isn't present in Model class yet
                sprintf('Error: One or more new model properties (variables) are not present in Model.m class yet!')
                return;
            end
        end
    end

    %%% Plotting
    if (plotIt == 1)
        % Vergence Error vs. iteration
        rng(0);
        for odIndex = 1 : size(objRange, 2)
            figure;
            hold on;
            grid on;
            grid minor;
            for vseIndex = 1 : size(model.testResult, 2)
                errorbar(0 : testInterval, model.testResult(odIndex, vseIndex, 1 : testInterval + 1), model.testResult(odIndex, vseIndex, testInterval + 2 : 2 * testInterval + 2), ...
                         'color', [rand, rand, rand], 'LineWidth', 1.3);
            end
            axis([-1, testInterval + 1, -inf, inf]);
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
        % [x, y] = meshgrid(0 : testInterval, objRange);
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
        xmin = 0;
        xmax = 0;
        for odIndex = 1 : length(objRange)
            % delta_mf_t+1(vergAngle_t)
            % hl3 = errorbar(reshape(reshape(model.testResult(odIndex, :, 1 : testInterval), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %                reshape(reshape(model.testResult(odIndex, :, 24 : 33), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %                reshape(reshape(model.testResult(odIndex, :, 35 : 44), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %                'DisplayName', sprintf('%.2fm objDist', objRange(odIndex)), 'Marker', '*', 'MarkerSize', 2.5, ...
            %                'color', [rand, rand, rand], 'LineWidth', 0.7, 'LineStyle', 'none');

            tmpMat = sortrows([reshape(reshape(model.testResult(odIndex, :, 1 : testInterval), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])', ...
                               reshape(reshape(model.testResult(odIndex, :, 2 * testInterval + 4 : 3 * testInterval + 3), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])', ...
                               reshape(reshape(model.testResult(odIndex, :, 3 * testInterval + 5 : 4 * testInterval + 4), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])']);

            [hl3, hp] = boundedline(tmpMat(:, 1), tmpMat(:, 2), tmpMat(:, 3), 'alpha');

            hl3.DisplayName = sprintf('%.2fm objDist', objRange(odIndex));
            hl3.Marker = '*';
            hl3.MarkerSize = 2.5;
            hl3.Color = [rand, rand, rand];
            hp.FaceColor = hl3.Color;
            hl3.LineStyle = 'none';
            % outlinebounds(hl3, hp);
            % lineHandles(odIndex + 2) = hl3;
            lineHandles(odIndex) = hl3;

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
        l.Location = 'southeast';
        l.Box = 'off';

        % adjust axis to actual response ranges + offset
        ymin = -0.1;
        ymax = 0.1;
        plot([xmin * 1.1, xmax * 1.1], [0, 0], 'k', 'LineWidth', 0.2);
        plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.2);
        axis([xmin * 1.1, xmax * 1.1, ymin, ymax]);
        xlabel(sprintf('Vergence Error [deg] (#stimuli=%d)', nStim), 'FontSize', 12);
        ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
        title('\Delta MF(verg_{err}) response at Testing procedure');
        % if (~isempty(model.savePath))
            plotpath = sprintf('%s/deltaMFasFktVerErr', imageSavePath);
            saveas(gcf, plotpath, 'png');
        % end

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
            % errorbar(reshape(reshape(model.testResult(odIndex, :, 1 : testInterval), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %          reshape(reshape(model.testResult(odIndex, :, 46 : 55), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %          reshape(reshape(model.testResult(odIndex, :, 57 : 66), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %          'Marker', '*', 'MarkerSize', 2.5, 'LineWidth', 0.9, 'LineStyle', 'none');

            tmpMat = sortrows([reshape(reshape(model.testResult(odIndex, :, 1 : testInterval), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])', ...
                               reshape(reshape(model.testResult(odIndex, :, 4 * testInterval + 6 : 5 * testInterval + 5), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])', ...
                               reshape(reshape(model.testResult(odIndex, :, 5 * testInterval + 7 : 6 * testInterval + 6), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])']);

            [hl, hp] = boundedline(tmpMat(:, 1), tmpMat(:, 2), tmpMat(:, 3), 'alpha');

            hl.DisplayName = sprintf('%.2fm objDist', objRange(odIndex));
            if (odIndex <= length(markers))
                hl.Marker = markers{odIndex};
            else
                hl.Marker = markers{randi(length(markers))};
            end
            hl.MarkerSize = 2.5;
            % hl.Color = [0, 0.5882, 0.9608];
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
        % l.Location = 'southeast';
        l.Box = 'off';

        % adjust axis to actual response ranges + std deviation
        xlim([xmin * 1.1, xmax * 1.1]);
        xlabel(sprintf('Vergence Error [deg] (#stimuli=%d)', nStim), 'FontSize', 12);
        ylabel('Value', 'FontSize', 12);
        title('Critic Value vs. Vergence Error');
        % if (~isempty(model.savePath))
            plotpath = sprintf('%s/criticValvsVerErr', imageSavePath);
            saveas(gcf, plotpath, 'png');
        % end

        % critic's response fine resolution
        vseRange = linspace(-1, 1, test2Resolution);
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

            [hl, hp] = boundedline(vseRange, ...
                                   tmpMean, ...
                                   tmpStd, ...
                                   'alpha');

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
            title(sprintf('Critic Value vs. Vergence Error\nobjDist = %.2fm, max@vergErr = %.3f°', objRange(odIndex), tmpMax));
            % if (~isempty(model.savePath))
                plotpath = sprintf('%s/criticValvsVerErrFine[%.2f].png', imageSavePath, objRange(odIndex));
                saveas(gcf, plotpath, 'png');
            % end
        end

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

        % if(version('-release') == '2015b')
        %     l.FontSize = 7;
        %     l.Orientation = 'horizontal';
        %     l.Location = 'southoutside';
        % end
        % xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
        % ylabel('Resonstruction Error', 'FontSize', 12);
        % title(sprintf('Reconstruction Error over different disparities\nobject distances: [%s]', num2str(objRange)));

        % if (~ isempty(model.savePath))
        %     plotpath = sprintf('%s/recErrVsVergErr[%.2fm,%.2fm].png', imageSavePath, objRange(1), objRange(end));
        %     saveas(gcf, plotpath, 'png');
        % end

        % reconstruction error fine
        vseRange = linspace(-1, 1, test2Resolution);
        handleArray = zeros(1, 1 + length(model.scModel));
        captions = cell(1, length(handleArray));
        captions{1} = 'Total Error';
        for i = 2 : length(handleArray)
            captions{i} = sprintf('Scale %d Error', i - 1);
        end

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
            if(version('-release') == '2015b')
                l.FontSize = 7;
                l.Orientation = 'horizontal';
                l.Location = 'southoutside';
            end
            xlim([vseRange(1) * 1.1, vseRange(end) * 1.1]);
            % xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
            xlabel('Vergence Error [deg]', 'FontSize', 12);
            ylabel('Resonstruction Error', 'FontSize', 12);
            title(sprintf('Reconstruction Error vs. Vergence Error\nobjDist = %.2fm, TotalMin@vergErr = %.4f°', objRange(odIndex), tmpMin));

            % if (~ isempty(model.savePath))
                plotpath = sprintf('%s/recErrVsVergErrFine[%.2fm].png', imageSavePath, objRange(odIndex));
                saveas(gcf, plotpath, 'png');
            % end
        end

        % Total error
        figure;
        hold on;
        grid on;
        b = boxplot(model.testResult3);

        % remove outliers
        % outl = findobj(b,'tag','Outliers');
        % set(outl, 'Visible', 'off');

        % rescale axis to whiskers + offset
        upWi = findobj(b, 'tag', 'Upper Whisker');
        lowWi = findobj(b, 'tag', 'Lower Whisker');
        axis([0, testInterval + 1, ...
              min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
              max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

        if (nStim > 0)
            xlabel(sprintf('Iteration step (#stimuli=%d)', nStim), 'FontSize', 12);
        else
            xlabel('Iteration step', 'FontSize', 12);
        end
        ylabel('Vergence Error [deg]', 'FontSize', 12);
        title(sprintf('Total Vergence Error over Trial at Testing\nMean = %.4f°, Median = %.4f°,\n4*IQR = %.4f, RMSE = %.4f° at %dth step', ...
                      mean(model.testResult3(:, testInterval)), median(model.testResult3(:, testInterval)), iqr(model.testResult3(:, testInterval)) * 4, sqrt(mean(model.testResult3(:, testInterval) .^ 2)), testInterval));
        % if (~isempty(model.savePath))
        plotpath = sprintf('%s/totalError', imageSavePath);
        saveas(gcf, plotpath, 'png');
        % end

        % results vector
        resultsOverview = {sqrt(mean(model.testResult3(:, testInterval) .^ 2)), median(model.testResult3(:, testInterval)), iqr(model.testResult3(:, testInterval)) * 4};
        formatSpec = '%.4f, %.4f, %.4f,';
        fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

        %% Check for bias at 0° vergence start error
        if (nStim == 0)
            nStim = size(model.testResult3, 1) / (length(objRange) * 7);
        end

        boxlabels = cell(1, length(objRange));
        tmp = zeros(nStim, length(objRange));
        startInd = nStim * 3 + 1;
        for i = 1 : length(objRange)
            endInd = startInd + nStim - 1;
            tmp(:, i) = model.testResult3(startInd : endInd, testInterval);
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
        boxplot(model.testResult3(:, testInterval), 'widths', 0.5, 'labels', 'total');
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
        suptitle(sprintf('Model bias with 0° vergence start error and total\nafter %d iterations for %d stimuli', ...
                         testInterval, nStim));
        % if (~isempty(model.savePath))
            plotpath = sprintf('%s/ModelBiasAt0VergErr', imageSavePath);
            saveas(gcf, plotpath, 'png');
        % end

        %%% Muscle correlation check
        % Total
        if (model.rlModel.CActor.output_dim >= 2)
            figure;
            hold on;
            scatter(model.testResult5(:, 1), model.testResult5(:, 2), 5,'MarkerFaceColor',[0, 0.7, 0.7]);
            corrl = corr(model.testResult5(:, 1), model.testResult5(:, 2));
            xlabel('Lateral rectus [%]', 'FontSize', 12);
            ylabel('Medial rectus [%]', 'FontSize', 12);
            title(strcat('Total Muscle Commands (testing)', sprintf('\nCorrelation = %1.2e', corrl)));
            plotpath = sprintf('%s/muscleGraphsScatterTotalTesting', imageSavePath);
            saveas(gcf, plotpath, 'png');

            % Delta
            figure;
            hold on;
            scatter(model.testResult5(:, 3), model.testResult5(:, 4), 5,'MarkerFaceColor',[0, 0.7, 0.7]);
            corrl = corr(model.testResult5(:, 3), model.testResult5(:, 4));
            xlabel('Lateral rectus [%]', 'FontSize', 12);
            ylabel('Medial rectus [%]', 'FontSize', 12);
            title(strcat('\Delta Muscle Commands (testing)', sprintf('\nCorrelation = %1.2e', corrl)));
            plotpath = sprintf('%s/muscleGraphsScatterDeltaTesting', imageSavePath);
            saveas(gcf, plotpath, 'png');
        end

        % save remaining results table
        resultsOverview = {'', '', '', '', '', '', '', '', '', '', '', '', '', '', GetFullPath(imageSavePath)};
        formatSpec = '%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n';
        fprintf(resultsFID, formatSpec, resultsOverview{1 : end});

        % close file
        fclose(resultsFID);
    end
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

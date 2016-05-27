%%% Model testing procedure
%@param model               respective model object to be tested
%@param nStim               # stimuli to be tested
%@pram plotIt               whether plots shall be generated
%@param saveTestResults     whether to save the results (not recommended if model is still trained!)
%@param simulator           simulator handle, provide [] if there is no simulator handle yet
%@param reinitRenderer      1 if renderer was already initialized
%                           0 if renderer wasn't initialized yet
%%%
function testModel2(model, nStim, plotIt, saveTestResults, simulator, reinitRenderer)

    % Vergence error resolution for 2nd testing procedure
    % Needs to be odd number to include vergErr = 0
    test2Resolution = 101;

    % Image processing variables
    textureFile = 'Textures_vanHaterenTest';

    % Prepare Textures
    texture = load(['config/' textureFile]);
    texture = texture.texture;
    nTextures = length(texture);
    if nTextures < nStim
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
        simulator = OpenEyeSim('create');
        if (reinitRenderer == 0)
            simulator.initRenderer();
        else
            % for debugging purposes
            simulator.reinitRenderer();
        end
    end

    % imageSavePath = model.savePath;
    imageSavePath = '.';

    % fixation interval at testing procedure
    testInterval = model.interval;

    command = [0; 0];
    objRange = [model.objDistMin : 0.5 : model.objDistMax];
    testResult = zeros(length(objRange), 7, 66);
    tmpResult1 = zeros(nStim, testInterval + 1);
    tmpResult2 = zeros(nStim, testInterval + 1);
    tmpResult3 = zeros(nStim, testInterval + 1);
    testResult2 = zeros(length(objRange) * 7 * nStim * testInterval, 1 + length(model.scModel));
    testResult3 = zeros(length(objRange) * 7 * nStim, 10);
    % testResult4 = zeros(length(objRange) * 7 * nStim * testInterval, 2);
    % testResult4 = zeros(test2Resolution * length(objRange) * nStim, 3 + length(model.scModel));
    testResult4 = zeros(length(objRange), test2Resolution, 2 + length(model.scModel));

    degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
    % metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

    % muscle function :=  mf(vergence_angle) = muscle force [single muscle]
    resolution = 100001;
    approx = spline(1 : 11, degrees.results_deg(:, 1));

    xValPos = ppval(approx, 1 : 0.0001 : 11)';
    yValPos = linspace(0, 1, resolution)';

    xValNeg = flipud(ppval(approx, 1 : 0.0001 : 11)' * -1);
    yValNeg = linspace(-1, 0, resolution)';

    mfunction = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
    mfunction(:, 1) = mfunction(:, 1) * 2;  % angle for two eyes
    dmf = diff(mfunction(1 : 2, 1));        % delta in angle
    dmf2 = diff(mfunction(1 : 2, 2));       % delta in mf
    indZero = find(mfunction(:, 2) == 0);   % MF == 0_index

    %%% Perfect Response function
    indMaxFix = find(mfunction(:, 1) <= model.vergAngleFixMin + dmf & mfunction(:, 1) >= model.vergAngleFixMin - dmf); % MF(vergAngleFixMin)_index
    indMaxFix = indMaxFix(1);
    indMinFix = find(mfunction(:, 1) <= model.vergAngleFixMax + dmf & mfunction(:, 1) >= model.vergAngleFixMax - dmf); % MF(vergAngleFixMax)_index
    indMinFix = indMinFix(1);

    % perfect_response := [max_fixation_x, max_fixation_y, min_fixation_x, min_fixation_y]
    % x = vergenceError, y = deltaMuscelForce
    perfectResponseMaxFix = [(mfunction(indMaxFix, 1) - flipud(mfunction(indMaxFix : end, 1))), ...
                             (mfunction(indMaxFix, 2) - flipud(mfunction(indMaxFix : end, 2))); ...
                             (mfunction(indMaxFix, 1) - flipud(mfunction(indZero : indMaxFix - 1, 1))), ...
                             (mfunction(indMaxFix, 2) - flipud(mfunction(indZero : indMaxFix - 1, 2)))];

    perfectResponseMinFix = [(mfunction(indMinFix, 1) - flipud(mfunction(indMinFix : end, 1))), ...
                             (mfunction(indMinFix, 2) - flipud(mfunction(indMinFix : end, 2))); ...
                             (mfunction(indMinFix, 1) - flipud(mfunction(indZero : indMinFix - 1, 1))), ...
                             (mfunction(indMinFix, 2) - flipud(mfunction(indZero : indMinFix - 1, 2)))];

    perfectResponse = [perfectResponseMaxFix, perfectResponseMinFix];

    % minimal and maximal angle that can be reached by one-dimensional muscle commands
    angleMin = mfunction(indZero, 1);
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
    % texture:  file path of texture input
    % eyeAngle: angle of single eye (rotation from offspring)
    % objDist:  distance of stimulus
    function refreshImages(texture, eyeAngle, objDist)
        simulator.add_texture(1, texture);
        simulator.set_params(1, eyeAngle, objDist);

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

    tic;
    tr2Ind = 1;
    tr3Ind = 1;
    % tr4Ind = 1;
    % don't repeat testing procedure if nStim == 0, but just plot the results
    if (nStim > 0)
        for odIndex = 1 : size(objRange, 2)
            sprintf('Testing iteration = %d/%d', odIndex, size(objRange, 2))

            % vergence start error
            vergMax = getVergErrMax(objRange(odIndex));
            if vergMax > 3
                vergMax = 3;
            end
            vseRange = [-3, -2, -1, linspace(0, vergMax, 4)];
            angleDes = 2 * atand(model.baseline / (2 * objRange(odIndex)));

            for vseIndex = 1 : size(vseRange, 2)
                tmpResult1(:, 1) = vseRange(vseIndex);

                for stimulusIndex = 1 : nStim
                    currentTexture = texture{stimulusIndex};
                    command(1) = 0;
                    [command(2), angleNew] = getMF(objRange(odIndex), vseRange(vseIndex));

                    for iter = 2 : testInterval + 1
                        % update stimuli
                        refreshImages(currentTexture, angleNew / 2, objRange(odIndex));

                        % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [imagesSavePath '/anaglyph.png']);
                        % generateAnaglyphs(imageSavePath, imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS);

                        % Image patch generation
                        for i = 1 : length(model.scModel)
                            model.preprocessImageFilled(imgGrayLeft, i, 1);
                            model.preprocessImageFilled(imgGrayRight, i, 2);
                            currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
                        end

                        % Generate input feature vector from current images
                        % [feature, ~, ~, errorLarge, errorSmall] = model.generateFR(currentView);
                        [bfFeature, reward, recErrorArray] = model.generateFR(currentView);

                        % Track reconstruction error statistics
                        testResult2(tr2Ind, :) = [angleDes - angleNew, recErrorArray];
                        tr2Ind = tr2Ind + 1;

                        % Absolute command feedback # concatination
                        if (model.rlModel.continuous == 1)
                            feature = [bfFeature; command * model.lambdaMuscleFB]; % single muscle
                        else
                            feature = bfFeature;
                        end

                        % Track Critic's response for single graph
                        % testResult4(tr4Ind, :) = [angleDes - angleNew, model.rlModel.CCritic.v_ji * feature];
                        % tr4Ind = tr4Ind + 1;

                        %%% Calculate metabolic costs
                        % metCost = getMetCost(command) * 2;

                        %%% Action
                        relativeCommand = model.rlModel.act(feature);

                        % add the change in muscle Activities to current ones
                        if (model.rlModel.continuous == 1)
                            % command(1) = 0;
                            % command(2) = command(2) + relativeCommand;  % one muscle
                            command = command + relativeCommand;  % one muscle
                            command = checkCmd(command);                % restrain motor commands to [0,1]
                            angleNew = getAngle(command) * 2;           % resulting angle is used for both eyes
                        else
                            angleNew = angleNew + relativeCommand;
                            if (angleNew > angleMax || angleNew < angleMin)
                                angleNew = model.vergAngleFixMin + (model.vergAngleFixMax - model.vergAngleFixMin) * rand(1, 1);
                            end
                        end

                        % track bad or redundant stimuli
                        % if (iter == 11)
                        %     if (abs(angleDes - angleNew) > 0.5)
                        %         sprintf('VergErr = %.1f\timage = %s\tstimulusIndex = %d\tobjDist = %.1f', (angleDes - angleNew), currentTexture, stimulusIndex, objRange(odIndex))
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
                end

                % final results
                testResult(odIndex, vseIndex, 1 : 11) = mean(tmpResult1);   % vergErr
                testResult(odIndex, vseIndex, 12 : 22) = std(tmpResult1);
                testResult(odIndex, vseIndex, 23 : 33) = mean(tmpResult2);  % deltaMF
                testResult(odIndex, vseIndex, 34 : 44) = std(tmpResult2);
                testResult(odIndex, vseIndex, 45 : 55) = mean(tmpResult3);  % critic's response
                testResult(odIndex, vseIndex, 56 : 66) = std(tmpResult3);
            end
        end

        %% Reconstruction error and Critic's response additional testing procedure
        tmp = zeros(nStim, 2 + length(model.scModel));
        % vergence start error
        vseRange = linspace(-1, 1, test2Resolution);
        angleDes = 2 * atand(model.baseline / (2 * objRange(odIndex)));

        for odIndex = 1 : length(objRange)
            for vseIndex = 1 : length(vseRange)
                if (getVergErrMax(objRange(odIndex)) < vseRange(vseIndex))
                    continue
                end
                command(1) = 0;
                [command(2), angleNew] = getMF(objRange(odIndex), vseRange(vseIndex));
                for stimulusIndex = 1 : nStim
                    % update stimuli
                    currentTexture = texture{stimulusIndex};
                    refreshImages(currentTexture, angleNew / 2, objRange(odIndex));

                    % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [imagesSavePath '/anaglyph.png']);
                    % generateAnaglyphs(imageSavePath, imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS);

                    % Image patch generation
                    for i = 1 : length(model.scModel)
                        model.preprocessImageFilled(imgGrayLeft, i, 1);
                        model.preprocessImageFilled(imgGrayRight, i, 2);
                        currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
                    end

                    % Generate input feature vector from current images
                    % [feature, ~, ~, errorLarge, errorSmall] = model.generateFR(currentView);
                    [bfFeature, reward, recErrorArray] = model.generateFR(currentView);

                    % Absolute command feedback # concatination
                    if (model.rlModel.continuous == 1)
                        % feature = [bfFeature; command(2) * model.lambdaMuscleFB]; % single muscle
                        feature = [bfFeature; command * model.lambdaMuscleFB]; % single muscle
                    else
                        feature = bfFeature;
                    end

                    % Track reconstruction error and Critic's response
                    tmp(stimulusIndex, :) = [model.rlModel.CCritic.v_ji * feature, sum(recErrorArray), recErrorArray];
                end
                testResult4(odIndex, vseIndex, :) = mean(tmp);
            end
        end
        testResult4(testResult4 == 0) = NaN;
        toc

        % save test results
        try
            model.testResult = testResult;
            model.testResult2 = testResult2;
            model.testResult3 = testResult3;
            model.testResult4 = testResult4;
            if (saveTestResults == 1)
                save(strcat(model.savePath, '/model'), 'model');
            end
        catch
            % catch non-existing variables error, occuring in old models
            clone = model.copy();
            delete(model);
            clear model;
            model = clone;
            model.testResult = testResult;
            model.testResult2 = testResult2;
            model.testResult3 = testResult3;
            model.testResult4 = testResult4;
            if (saveTestResults == 1)
                save(strcat(model.savePath, '/model'), 'model');
            end
            delete(clone);
            clear clone;
        end
    end

    if (plotIt == 1)
        % Vergence Error vs. iteration
        rng(0);
        for odIndex = 1 : size(objRange, 2)
            figure;
            hold on;
            grid on;
            grid minor;
            for vseIndex = 1 : size(model.testResult, 2)
                errorbar(0 : testInterval, model.testResult(odIndex, vseIndex, 1 : 11), model.testResult(odIndex, vseIndex, 12 : 22), ...
                         'color', [rand, rand, rand], 'LineWidth', 1.3);
            end
            axis([-1, testInterval + 1, -inf, inf]);
            if (nStim > 0)
                xlabel(sprintf('Iteration step (objDist=%.1fm, #stimuli=%d)', objRange(odIndex), nStim), 'FontSize', 12);
            else
                xlabel(sprintf('Iteration step (objDist=%.1fm)', objRange(odIndex)), 'FontSize', 12);
            end
            ylabel('Vergence Error [deg]', 'FontSize', 12);
            title('Avg Vergence Error over Trial at Testing');
            plotpath = sprintf('%s/AvgVergErrOverTrial_objDist[%.1fm].png', model.savePath, objRange(odIndex));
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
        % plotpath = sprintf('%s/AvgVergErrOverTrial3D', model.savePath);
        % saveas(gcf, plotpath, 'png');

        % deltaMuscleForce vs Vergence Error
        figure;
        hold on;
        grid on;
        % perfect response to vergence error
        hl1 = plot(perfectResponse(:, 1), perfectResponse(:, 2), 'color', [0.5882, 0.9608, 0], ...
                   'DisplayName', 'perfect (fixDist_{max})', 'LineWidth', 1.3);
        hl2 = plot(perfectResponse(:, 3), perfectResponse(:, 4), 'color', [0, 0.5882, 0.9608], ...
                   'DisplayName', 'perfect (fixDist_{min})', 'LineWidth', 1.3);
        lineHandles = [hl1, hl2];

        % actual response
        xmin = 0;
        xmax = 0;
        for odIndex = 1 : size(objRange, 2)
            % delta_mf_t+1(vergAngle_t)
            % hl3 = errorbar(reshape(reshape(model.testResult(odIndex, :, 1 : testInterval), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %                reshape(reshape(model.testResult(odIndex, :, 24 : 33), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %                reshape(reshape(model.testResult(odIndex, :, 35 : 44), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %                'DisplayName', sprintf('%.1fm objDist', objRange(odIndex)), 'Marker', '*', 'MarkerSize', 2.5, ...
            %                'color', [rand, rand, rand], 'LineWidth', 0.7, 'LineStyle', 'none');

            tmpMat = sortrows([reshape(reshape(model.testResult(odIndex, :, 1 : testInterval), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])', ...
                               reshape(reshape(model.testResult(odIndex, :, 24 : 33), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])', ...
                               reshape(reshape(model.testResult(odIndex, :, 35 : 44), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])']);

            [hl3, hp] = boundedline(tmpMat(:, 1), tmpMat(:, 2), tmpMat(:, 3), 'alpha');

            hl3.DisplayName = sprintf('%.1fm objDist', objRange(odIndex));
            hl3.Marker = '*';
            hl3.MarkerSize = 2.5;
            hl3.Color = [rand, rand, rand];
            hp.FaceColor = hl3.Color;
            hl3.LineStyle = 'none';
            % outlinebounds(hl3, hp);
            lineHandles = [lineHandles, hl3];

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
        % l.Title.String = 'objDist [m]';
        % l.Title.FontSize = 12;
        xlabel(sprintf('Vergence Error [deg] (#stimuli=%d)', nStim), 'FontSize', 12);
        ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
        title('\Delta MF(verg_{err}) response at Testing procedure');
        if (~isempty(model.savePath))
            plotpath = sprintf('%s/deltaMFasFktVerErr', model.savePath);
            saveas(gcf, plotpath, 'png');
        end

        % critic's response
        figure;
        hold on;
        grid on;
        grid minor;
        xmin = 0;
        xmax = 0;
        for odIndex = 1 : size(objRange, 2)
            % delta_mf_t+1(vergAngle_t)
            % errorbar(reshape(reshape(model.testResult(odIndex, :, 1 : testInterval), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %          reshape(reshape(model.testResult(odIndex, :, 46 : 55), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %          reshape(reshape(model.testResult(odIndex, :, 57 : 66), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval]), ...
            %          'Marker', '*', 'MarkerSize', 2.5, 'LineWidth', 0.9, 'LineStyle', 'none');

            tmpMat = sortrows([reshape(reshape(model.testResult(odIndex, :, 1 : testInterval), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])', ...
                               reshape(reshape(model.testResult(odIndex, :, 46 : 55), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])', ...
                               reshape(reshape(model.testResult(odIndex, :, 57 : 66), [size(model.testResult, 2), testInterval])', [1, size(model.testResult, 2) * testInterval])']);

            [hl, hp] = boundedline(tmpMat(:, 1), tmpMat(:, 2), tmpMat(:, 3), 'alpha');

            hl.DisplayName = sprintf('%.1fm objDist', objRange(odIndex));
            hl.Marker = '*';
            hl.MarkerSize = 2.5;
            hl.Color = [0, 0.5882, 0.9608];
            hp.FaceColor = hl.Color;
            hl.LineStyle = 'none';

            % for axis adjustment
            tmp = [min(tmpMat(:, 1)), max(tmpMat(:, 1))];
            if (xmin > tmp(1))
                xmin = tmp(1);
            end
            if (xmax < tmp(2))
                xmax = tmp(2);
            end
        end

        % adjust axis to actual response ranges + std deviation
        ymin = -inf;
        ymax = inf;
        axis([xmin * 1.1, xmax * 1.1, ymin, ymax]);
        xlabel(sprintf('Vergence Error [deg] (#stimuli=%d)', nStim), 'FontSize', 12);
        ylabel('Value', 'FontSize', 12);
        title('Critic Value over different disparities');
        if (~isempty(model.savePath))
            plotpath = sprintf('%s/criticValvsVerErr', model.savePath);
            saveas(gcf, plotpath, 'png');
        end

        % critic's response fine resolution
        figure;
        hold on;
        grid on;
        grid minor;

        vseRange = linspace(-1, 1, test2Resolution);
        % errorbar(vseRange, mean(model.testResult4(:, :, 1), 'omitnan'), std(model.testResult4(:, :, 1), 'omitnan'));
        % [hl3, hp] = boundedline(vseRange, mean(model.testResult4(:, :, 1), 'omitnan'), std(model.testResult4(:, :, 1), 'omitnan'), 'alpha');
        [hl, hp] = boundedline(vseRange, mean(model.testResult4(:, :, 1), 'omitnan'), std(model.testResult4(:, :, 1), 'omitnan'));

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
        title('Critic Value over different disparities');
        if (~isempty(model.savePath))
            plotpath = sprintf('%s/criticValvsVerErrFine', model.savePath);
            saveas(gcf, plotpath, 'png');
        end

        %%% Plot the resonstruction error of basis functions over different disparities
        nBins = 1000;
        % calculate mean and std of reconstruction error
        tmpRsp = sortrows(model.testResult2);
        deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nBins;
        % recErrs = nBins x [recErr; total_mean; total_std; scale1_mean; scale1_std; ...]
        recErrs = zeros(nBins, 1 + 2 * (length(model.scModel) + 1));
        tmp = zeros(nBins, 3);

        % total reconstruction error
        for i = 1 : nBins
            tmp(i, 1) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                         & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 1));

            tmp(i, 2) = mean(sum(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                        & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2 : end), 2));

            tmp(i, 3) = std(sum(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                       & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2 : end), 2));
        end
        recErrs(:, 1 : 3) = tmp;

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
        %     plotpath = sprintf('%s/recErrVsVergErr_[%.1fm,%.1fm].png', model.savePath, objRange(1), objRange(end));
        %     saveas(gcf, plotpath, 'png');
        % end

        % reconstruction error fine
        figure;
        hold on;
        grid on;
        grid minor;
        handleArray = zeros(1, 1 + length(model.scModel));
        vseRange = linspace(-1, 1, test2Resolution);

        for i = 2 : size(model.testResult4, 3)
            handleArray(i - 1) = errorbar(vseRange, ...
                                          mean(model.testResult4(:, :, i), 'omitnan'), ...
                                          std(model.testResult4(:, :, i), 'omitnan'), ...
                                          'LineWidth', 0.9);
        end

        captions = cell(1, length(handleArray));
        captions{1} = 'Total Error';
        for i = 2 : length(handleArray)
            captions{i} = sprintf('Scale %d Error', i - 1);
        end
        l = legend(handleArray, captions);

        if(version('-release') == '2015b')
            l.FontSize = 7;
            l.Orientation = 'horizontal';
            l.Location = 'southoutside';
        end
        xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
        ylabel('Resonstruction Error', 'FontSize', 12);
        title(sprintf('Reconstruction Error over different disparities\nobject distances: [%s]', num2str(objRange)));

        if (~ isempty(model.savePath))
            plotpath = sprintf('%s/recErrVsVergErrFine_[%.1fm,%.1fm].png', model.savePath, objRange(1), objRange(end));
            saveas(gcf, plotpath, 'png');
        end

        % Total error
        figure;
        hold on;
        grid on;
        b = boxplot(model.testResult3);

        % remove outliers
        outl = findobj(b,'tag','Outliers');
        set(outl, 'Visible', 'off');

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
        title(sprintf('Total Vergence Error over Trial at Testing\nMean = %5.3f, Median = %5.3f', mean(testResult3(:, 10)), median(testResult3(:, 10))));
        if (~isempty(model.savePath))
            plotpath = sprintf('%s/totalError', model.savePath);
            saveas(gcf, plotpath, 'png');
        end
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

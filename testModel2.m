%%% Model testing procedure
%@param model               respective model object to be tested
%@param nStim               # stimuli to be tested
%@pram plotIt               whether plots shall be generated
%@param saveTestResults     whether to save the results (not recommended if model is still trained!)
%@param reinitRenderer      1 if renderer was already initialized, e.g. when training was conducted before
%                           0 if testModel is called stand-alone
%%%
function testModel2(model, nStim, plotIt, saveTestResults, reinitRenderer)
    % cancel testing procedure
    if (nStim == 0)
        return;
    end

    %%% New renderer
    simulator = OpenEyeSim('create');
    % for debugging or if testing procedure shall be executed
    % without prior training procedure
    if (reinitRenderer == 1)
        simulator.reinitRenderer();
    else
        simulator.initRenderer();
    end

    % imageSavePath = model.savePath;
    imageSavePath = '.';

    % fixation interval at testing procedure
    testInterval = model.interval;

    command = [0; 0];
    objRange = [model.objDistMin : 0.5 : model.objDistMax];
    testResult = zeros(size(objRange, 2), 7, 66);
    tmpResult1 = zeros(nStim, testInterval + 1);
    tmpResult2 = zeros(nStim, testInterval + 1);
    tmpResult3 = zeros(nStim, testInterval + 1);

    % Image processing variables
    textureFile = 'Textures_vanHaterenTest';

    % Prepare Textures
    texture = load(['config/' textureFile]);
    texture = texture.texture;
    nTextures = length(texture);

    degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
    % metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

    % muscle function :=  mf(vergence_angle) = muscle force [single muscle]
    resolution = 100001;
    approx = spline(1:11, degrees.results_deg(:, 1));

    xValPos = ppval(approx, 1:0.0001:11)';
    yValPos = linspace(0, 1, resolution)';

    xValNeg = flipud(ppval(approx, 1:0.0001:11)' * -1);
    yValNeg = linspace(-1, 0, resolution)';

    mfunction = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
    mfunction(:, 1) = mfunction(:, 1) * 2;  % angle for two eyes
    dmf = diff(mfunction(1 : 2, 1));        % delta in angle
    dmf2 = diff(mfunction(1 : 2, 2));       % delta in mf
    indZero = find(mfunction(:, 2) == 0);   % MF == 0_index

    %%% Perfect Response function
    indMaxFix = find(mfunction(:, 1) <= model.vergAngleMin + dmf & mfunction(:, 1) >= model.vergAngleMin - dmf); % MF(vergAngleMin)_index
    indMaxFix = indMaxFix(1);
    indMinFix = find(mfunction(:, 1) <= model.vergAngleMax + dmf & mfunction(:, 1) >= model.vergAngleMax - dmf); % MF(vergAngleMax)_index
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
        cmd = (command * 10) + 1;                               % calculate tabular index
        angle = interp2(degrees.results_deg, cmd(1), cmd(2), 'spline');   % interpolate in tabular
    end

    function angle = getAngle2(command)
        angleIndex = find(mfunction(:, 2) <= command(2) + dmf2 & mfunction(:, 2) >= command(2) - dmf2);
        angle = mfunction(angleIndex, 1);
        angle = angle(ceil(length(angle) / 2));
    end

    %%% Helper function that maps muscle activities to resulting metabolic costs
    % function [tmpMetCost] = getMetCost(command)
    %     cmd = (command * 10) + 1;                               % scale commands to table entries
    %     tmpMetCost = interp2(metCosts.results, cmd(1), cmd(2)); % interpolate in tabular
    % end

    % Generates two new images for both eyes
    function refreshImages(texture, vergAngle, objDist)
        simulator.add_texture(1, texture);
        simulator.set_params(1, vergAngle, objDist);

        result1 = simulator.generate_left();
        result2 = simulator.generate_right();

        k = 1;
        l = 1;
        for i = 1 : 3 : length(result1)
            imgRawLeft(k, l, 1) = result1(i);
            imgRawLeft(k, l, 2) = result1(i + 1);
            imgRawLeft(k, l, 3) = result1(i + 2);

            imgRawRight(k, l, 1) = result2(i);
            imgRawRight(k, l, 2) = result2(i + 1);
            imgRawRight(k, l, 3) = result2(i + 2);

            l = l + 1;
            if (l > 320)
                l = 1;
                k = k + 1;
            end
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

    %%% Helper function for image preprocessing
    %% Patch generation
    function [patches] = preprocessImage(img, imgSize)
        % img = .2989 * img(:,:,1) + .5870 * img(:,:,2) + .1140 * img(:,:,3);

        if (imgSize == 0)
            % small scale
            for i = 1:log2(model.dsRatioS)
                img = impyramid(img, 'reduce');
            end

            % convert to double
            img = double(img);

            % cut fovea in the center
            [h, w, ~] = size(img);
            img = img(fix(h / 2 + 1 - model.foveaS / 2) : fix(h / 2 + model.foveaS / 2), ...
                      fix(w / 2 + 1 - model.foveaS / 2) : fix(w / 2 + model.foveaS / 2));
        else
            % large scale
            for i = 1:log2(model.dsRatioL)
                img = impyramid(img, 'reduce');
            end

            % convert to double
            img = double(img);

            % cut fovea in the center
            [h, w, ~] = size(img);
            img = img(fix(h / 2 + 1 - model.foveaL / 2) : fix(h / 2 + model.foveaL / 2), ...
                      fix(w / 2 + 1 - model.foveaL / 2) : fix(w / 2 + model.foveaL / 2));
        end

        % cut patches and store them as col vectors
        patches = im2col(img, [model.patchSize model.patchSize], 'sliding');            %slide window of 1 px

        % take patches at steps of s (8 px)
        if (imgSize == 0)
            patches = patches(:, model.columnIndS);                               %81 patches
        else
            patches = patches(:, model.columnIndL);                               %81 patches
        end

        % pre-processing steps (0 mean, unit norm)
        patches = patches - repmat(mean(patches), [size(patches, 1) 1]);    %0 mean
        normp = sqrt(sum(patches.^2));                                      %patches norm

        % normalize patches to norm 1
        normp(normp == 0) = eps;                                            %regularizer
        patches = patches ./ repmat(normp, [size(patches, 1) 1]);           %normalized patches
    end

    tic;
    for odIndex = 1 : size(objRange, 2)
        sprintf('Testing iteration = %d/%d', odIndex, size(objRange, 2))

        % vergence start error
        vseRange = [-3, -2, -1, linspace(0, getVergErrMax(objRange(odIndex)), 4)];
        angleDes = 2 * atand(model.baseline / (2 * objRange(odIndex)));

        for vseIndex = 1 : size(vseRange, 2)
            tmpResult1(:, 1) = vseRange(vseIndex);

            for stimulusIndex = 1 : nStim
                currentTexture = texture{stimulusIndex};
                command(1) = 0;
                [command(2), angleNew] = getMF(objRange(odIndex), vseRange(vseIndex));
                refreshImages(currentTexture, angleNew / 2, objRange(odIndex));

                %%% DEBUGGING
                % [status, res] = system(sprintf('./checkEnvironment %s %d %d %s/leftTest.png %s/rightTest.png', ...
                %                                'a.bmp',2, 2,imageSavePath, imageSavePath));

                % refreshImages('b.bmp', 2/2, 2);

                % figure;
                % subplot(1,2,1)
                % imshow(imfuse(imread([imageSavePath '/leftTest.png']), imgRawLeft, 'falsecolor'));
                % subplot(1,2,2)
                % imshow(imfuse(imread([imageSavePath '/rightTest.png']), imgRawRight, 'falsecolor'));

                % figure
                % subplot(2,2,1)
                % imshow(imread([imageSavePath '/leftTest.png']));
                % subplot(2,2,2)
                % imshow(imread([imageSavePath '/rightTest.png']));
                % subplot(2,2,3)
                % imshow(imgRawLeft);
                % subplot(2,2,4)
                % imshow(imgRawRight);
                %%% DEBUGGING END

                for iter = 2 : testInterval + 1
                    % convert images to gray scale
                    imgGrayLeft = .2989 * imgRawLeft(:,:,1) + .5870 * imgRawLeft(:,:,2) + .1140 * imgRawLeft(:,:,3);
                    imgGrayRight = .2989 * imgRawRight(:,:,1) + .5870 * imgRawRight(:,:,2) + .1140 * imgRawRight(:,:,3);

                    % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [imagesSavePath '/anaglyph.png']);
                    % generateAnaglyphs(imageSavePath, imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS);

                    % Image patch generation: left{small scale, large scale}, right{small scale, large scale}
                    [patchesLeftSmall] = preprocessImage(imgGrayLeft, 0);
                    [patchesLeftLarge] = preprocessImage(imgGrayLeft, 1);
                    [patchesRightSmall] = preprocessImage(imgGrayRight, 0);
                    [patchesRightLarge] = preprocessImage(imgGrayRight, 1);

                    % Image patches matrix (input to model)
                    currentView = {[patchesLeftLarge; patchesRightLarge] [patchesLeftSmall; patchesRightSmall]};

                    % Generate input feature vector from current images
                    [feature, ~, ~, errorLarge, errorSmall] = model.generateFR(currentView);

                    % Absolute command feedback # concatination
                    if (model.rlmodel.continuous == 1)
                        feature = [feature; command(2) * model.lambdaMuscleFB]; % single muscle
                        % feature = [feature; command * model.lambdaMuscleFB]; % two muscles
                    end

                    %%% Calculate metabolic costs
                    % metCost = getMetCost(command) * 2;

                    %%% Action
                    relativeCommand = model.rlmodel.act(feature);

                    % add the change in muscle Activities to current ones
                    if (model.rlmodel.continuous == 1)
%                         command = command + relativeCommand;     %two muscels
%                         command(1) = 0;
                        command(2) = command(2) + relativeCommand;  %one muscel
                        command = checkCmd(command);                %restrain motor commands to [0,1]
%                         angleNew = getAngle2(command);        %resulting angle is used for one eye
                        angleNew = getAngle(command) * 2;           %resulting angle is used for both eyes
                    else
                        angleNew = angleNew + relativeCommand;
                        if (angleNew > angleMax || angleNew < angleMin)
                            angleNew = model.vergAngleMin + (model.vergAngleMax - model.vergAngleMin) * rand(1, 1);
                        end
                    end

                    refreshImages(currentTexture, angleNew / 2, objRange(odIndex));

                    % temporary results
                    tmpResult1(stimulusIndex, iter) = angleDes - angleNew;
                    tmpResult2(stimulusIndex, iter) = relativeCommand(1);
                    tmpResult3(stimulusIndex, iter) = model.rlmodel.CCritic.v_ji * feature;
                end
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
    toc

    % save test results
    try
        model.testResult = testResult;
        if (saveTestResults == 1)
            save(strcat(model.savePath, '/model'), 'model');
        end
    catch
        % catch non-existing variables error, occuring at old models
        clone = model.copy();
        delete(model);
        clear model;
        model = clone;
        model.testResult = testResult;
        if (saveTestResults == 1)
            save(strcat(model.savePath, '/model'), 'model');
        end
        delete(clone);
        clear clone;
    end

    if (plotIt == 1)
        % Vergence Error vs. iteration
        for odIndex = 1 : size(objRange, 2)
            figure;
            hold on;
            grid on;
            for vseIndex = 1 : 7
                errorbar(0 : 10, testResult(odIndex, vseIndex, 1 : 11), testResult(odIndex, vseIndex, 12 : 22), ...
                         'color', [rand, rand, rand], 'LineWidth', 1.3);
            end
            axis([-1, testInterval + 1, -inf, inf]);
            xlabel(sprintf('Iteration step (objDist=%.1fm, #stimuli=%d)', objRange(odIndex), nStim), 'FontSize', 12);
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
        % for vseIndex = 1 : 7
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
        for odIndex = 1 : size(objRange, 2)
            % delta_mf_t+1(vergAngle_t)
            % hl3 = errorbar(reshape(reshape(testResult(odIndex, :, 1 : 10), [7, 10])', [1, 7 * 10]), ...
            %                reshape(reshape(testResult(odIndex, :, 24 : 33), [7, 10])', [1, 7 * 10]), ...
            %                reshape(reshape(testResult(odIndex, :, 35 : 44), [7, 10])', [1, 7 * 10]), ...
            %                'DisplayName', sprintf('%.1fm objDist', objRange(odIndex)), 'Marker', '*', 'MarkerSize', 2.5, ...
            %                'color', [rand, rand, rand], 'LineWidth', 0.7, 'LineStyle', 'none');

            tmpMat = sortrows([reshape(reshape(testResult(odIndex, :, 1 : 10), [7, 10])', [1, 7 * 10])', ...
                               reshape(reshape(testResult(odIndex, :, 24 : 33), [7, 10])', [1, 7 * 10])', ...
                               reshape(reshape(testResult(odIndex, :, 35 : 44), [7, 10])', [1, 7 * 10])']);

            [hl3, hp] = boundedline(tmpMat(:, 1), tmpMat(:, 2), tmpMat(:, 3), 'alpha');

            hl3.DisplayName = sprintf('%.1fm objDist', objRange(odIndex));
            hl3.Marker = '*';
            hl3.MarkerSize = 2.5;
            hl3.Color = [rand, rand, rand];
            hp.FaceColor = hl3.Color;
            hl3.LineStyle = 'none';
            % outlinebounds(hl3, hp);
            lineHandles = [lineHandles, hl3];
        end
        l = legend(lineHandles);
        l.Location = 'southeast';
        l.Box = 'off';

        % adjust axis to actual response ranges + std deviation
        xmin = -4;
        xmax = 6;
        ymin = -0.1;
        ymax = 0.1;
        plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.2);
        plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.2);
        axis([xmin, xmax, ymin, ymax]);
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
        for odIndex = 1 : size(objRange, 2)
            % delta_mf_t+1(vergAngle_t)
            % errorbar(reshape(reshape(testResult(odIndex, :, 1 : 10), [7, 10])', [1, 7 * 10]), ...
            %          reshape(reshape(testResult(odIndex, :, 46 : 55), [7, 10])', [1, 7 * 10]), ...
            %          reshape(reshape(testResult(odIndex, :, 57 : 66), [7, 10])', [1, 7 * 10]), ...
            %          'Marker', '*', 'MarkerSize', 2.5, 'LineWidth', 0.9, 'LineStyle', 'none');

            tmpMat = sortrows([reshape(reshape(testResult(odIndex, :, 1 : 10), [7, 10])', [1, 7 * 10])', ...
                               reshape(reshape(testResult(odIndex, :, 46 : 55), [7, 10])', [1, 7 * 10])', ...
                               reshape(reshape(testResult(odIndex, :, 57 : 66), [7, 10])', [1, 7 * 10])']);

            [hl, hp] = boundedline(tmpMat(:, 1), tmpMat(:, 2), tmpMat(:, 3), 'alpha');

            hl.DisplayName = sprintf('%.1fm objDist', objRange(odIndex));
            hl.Marker = '*';
            hl.MarkerSize = 2.5;
            hl.Color = [0, 0.5882, 0.9608];
            hp.FaceColor = hl.Color;
            hl.LineStyle = 'none';
        end

        % adjust axis to actual response ranges + std deviation
        xmin = -4;
        xmax = 6;
        ymin = -inf;
        ymax = inf;
        % plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
        % plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
        axis([xmin, xmax, ymin, ymax]);
        xlabel(sprintf('Vergence Error [deg] (#stimuli=%d)', nStim), 'FontSize', 12);
        ylabel('Value', 'FontSize', 12);
        title('Critic Value over different disparities');
        if (~isempty(model.savePath))
            plotpath = sprintf('%s/criticValvsVerErr', model.savePath);
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

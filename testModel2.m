%%% Model testing procedure
%@param model               respective model object to be tested
%@param nStim               # stimuli to be tested
%@pram plotIt               whether plots shall be generated
%@param saveTestResults     whether to save the results (not recommended if model is still trained!)
%%%
function testModel2(model, nStim, plotIt, saveTestResults)
    % cancel testing procedure
    if (nStim == 0)
        return;
    end

    command = [0, 0];
    objRange = [model.objDistMin : 0.5 : model.objDistMax];
    testResult = zeros(size(objRange, 2), 7, 44);
    tmpResult1 = zeros(nStim, model.interval + 1);
    tmpResult2 = zeros(nStim, model.interval + 1);

    % Image processing variables
    textureFile = 'Textures_vanHaterenTest';
    patchSize = 8;

    dsRatioL = model.scmodel_Large.Dsratio; %downsampling ratio (Large scale) | original 8
    dsRatioS = model.scmodel_Small.Dsratio; %downsampling ratio (Small scale) | original 2

    % fovea = [128 128];
    foveaL = patchSize + patchSize ^ 2 / 2 ^ log2(dsRatioL); %fovea size (Large scale) | 16
    foveaS = patchSize + patchSize ^ 2 / 2 ^ log2(dsRatioS); %fovea size (Small scale) | 40

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

    % Prepare Textures
    texture = load(['config/' textureFile]);
    texture = texture.texture;
    nTextures = length(texture);

    degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
    % metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

    %TODO: define new results variables
    % disZtest = zeros(model.interval, size(objRange, 2), repeat(1));        % desired fixation distance
    % fixZtest = zeros(model.interval, size(objRange, 2), repeat(1));        % actual fixation distance
    % vergErrTest = zeros(model.interval, size(objRange, 2), repeat(1));     % vergence error

    % minimal and maximal angle that can be reached by one-dimensional muscle commands
    angleMin = getAngle([0, 0]) * 2;
    angleMax = getAngle([0, 1]) * 2;

    % muscle function :=  mf(vergence_angle) = muscle force [single muscle]
    resolution = 10001;
    approx = spline(1:11, degrees.results_deg(:, 1));

    xValPos = ppval(approx, 1:0.001:11)';
    yValPos = linspace(0, 1, resolution)';

    xValNeg = flipud(ppval(approx, 1:0.001:11)' * -1);
    yValNeg = linspace(-1, 0, resolution)';

    mfunction = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
    mfunction(:, 1) = mfunction(:, 1) * 2;  % angle for two eyes
    dmf = diff(mfunction(1 : 2, 1));        % delta in angle
    % indZero = find(mfunction(:, 2) == 0); % MF == 0_index

    %%% Helper function that maps {objDist, desiredVergErr} -> {muscleForce, angleInit}
    function [mf, angleInit] = getMF(objDist, desVergErr)
        % correct vergence angle for given object distance
        angleCorrect = 2 * atand(model.baseline / (2 * objDist));
        % desired init angle for given vergence error [deg]
        angleInit = angleCorrect - desVergErr;
        % look up index of angleInit
        indAngleInit = find(mfunction(:, 1) <= angleInit + dmf & mfunction(:, 1) >= angleInit - dmf);
        mf = mfunction(indAngleInit, 2);
        mf = mf(1);
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
    function [angle] = getAngle(command)
        cmd = (command * 10) + 1;                               % calculate tabular index
        angle = interp2(degrees.results_deg, cmd(1), cmd(2));   % interpolate in tabular
    end

    %%% Helper function that maps muscle activities to resulting metabolic costs
    % function [tmpMetCost] = getMetCost(command)
    %     cmd = (command * 10) + 1;                               % scale commands to table entries
    %     tmpMetCost = interp2(metCosts.results, cmd(1), cmd(2)); % interpolate in tabular
    % end

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

                [status, res] = system(sprintf('./checkEnvironment %s %d %d %s/left.png %s/right.png', ...
                                               currentTexture, objRange(odIndex), angleNew, model.savePath, model.savePath));
                % abort execution if error occured
                if (status)
                    sprintf('Error in checkEnvironment:\n%s', res)
                    return;
                end

                for iter = 2 : model.interval + 1
                    % read input images and convert to gray scale
                    imgRawLeft = imread([model.savePath '/left.png']);
                    imgRawRight = imread([model.savePath '/right.png']);
                    imgGrayLeft = .2989 * imgRawLeft(:,:,1) + .5870 * imgRawLeft(:,:,2) + .1140 * imgRawLeft(:,:,3);
                    imgGrayRight = .2989 * imgRawRight(:,:,1) + .5870 * imgRawRight(:,:,2) + .1140 * imgRawRight(:,:,3);

                    % generateAnaglyphs(model, imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS);

                    % Image patch generation: left{small scale, large scale}, right{small scale, large scale}
                    [patchesLeftSmall] = preprocessImage(imgGrayLeft, foveaS, dsRatioS, patchSize, columnIndS);
                    [patchesLeftLarge] = preprocessImage(imgGrayLeft, foveaL, dsRatioL, patchSize, columnIndL);
                    [patchesRightSmall] = preprocessImage(imgGrayRight, foveaS, dsRatioS, patchSize, columnIndS);
                    [patchesRightLarge] = preprocessImage(imgGrayRight, foveaL, dsRatioL, patchSize, columnIndL);

                    % Image patches matrix (input to model)
                    currentView = {[patchesLeftLarge; patchesRightLarge] [patchesLeftSmall; patchesRightSmall]};

                    % Generate input feature vector from current images
                    [feature, ~, ~, errorLarge, errorSmall] = model.generateFR(currentView);

                    % Absolute command feedback # concatination
                    if (model.rlmodel.continuous == 1)
                        feature = [feature; command(2) * model.lambdaMuscleFB];
                    end

                    %%% Calculate metabolic costs
                    % metCost = getMetCost(command) * 2;

                    %%% Action
                    relativeCommand = model.rlmodel.act(feature);

                    % add the change in muscle Activities to current ones
                    if (model.rlmodel.continuous == 1)
                        % command = command + relativeCommand';     %two muscels
                        command(1) = 0;
                        command(2) = command(2) + relativeCommand;  %one muscel
                        command = checkCmd(command);                %restrain motor commands to [0,1]
                        angleNew = getAngle(command) * 2;           %resulting angle is used for both eyes
                    else
                        angleNew = angleNew + relativeCommand;
                        if (angleNew > angleMax || angleNew < angleMin)
                            angleNew = (model.desiredAngleMin + (model.desiredAngleMax - model.desiredAngleMin) * rand(1, 1)) * 2;
                        end
                    end

                    % generate new view (two pictures) with new vergence angle
                    [status, res] = system(sprintf('./checkEnvironment %s %d %d l%s/eft.png %s/right.png', ...
                                                   currentTexture, objRange(odIndex), angleNew, model.savePath, model.savePath));

                    % abort execution if error occured
                    if (status)
                        sprintf('Error in checkEnvironment:\n%s', res)
                        return;
                    end

                    % temporary results
                    tmpResult1(stimulusIndex, iter) = angleDes - angleNew;
                    tmpResult2(stimulusIndex, iter) = relativeCommand;
                end
            end

            % final results
            testResult(odIndex, vseIndex, 1 : 11) = mean(tmpResult1); % vergErr
            testResult(odIndex, vseIndex, 12 : 22) = std(tmpResult1);
            testResult(odIndex, vseIndex, 23 : 33) = mean(tmpResult2); % deltaMF
            testResult(odIndex, vseIndex, 34 : 44) = std(tmpResult2);
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
            xlabel('Iteration step', 'FontSize', 12);
            ylabel('Vergence Error [deg]', 'FontSize', 12);
            title(sprintf('Avg Vergence Error over Trial at Testing (objDist = %.1fm)', objRange(odIndex)));
            plotpath = sprintf('%s/AvgVergErrOverTrial_objDist[%.1fm]', model.savePath, objRange(odIndex));
            saveas(gcf, plotpath, 'png');
        end

        % 3D plot
        % figure;
        % hold on;
        % grid on;
        % [x, y] = meshgrid(0 : model.interval, objRange);
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
        for odIndex = 1 : size(objRange, 2)
            % figure;
            % hold on;
            % grid on;
            errorbar(reshape(reshape(testResult(odIndex, :, 1 : 11), [7, 11])', [1, 7 * 11]), ...
                     reshape(reshape(testResult(odIndex, :, 23 : 33), [7, 11])', [1, 7 * 11]), ...
                     reshape(reshape(testResult(odIndex, :, 34 : 44), [7, 11])', [1, 7 * 11]), ...
                     'DisplayName', num2str(objRange(odIndex)), 'Marker', '*', 'MarkerSize', 3, ...
                     'color', [rand, rand, rand], 'LineWidth', 0.9, 'LineStyle', 'none');
            legend('-DynamicLegend');
            % xlabel('Vergence Error [deg]', 'FontSize', 12);
            % ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
            % title(strcat('\Delta MF(verg_{err}) response at Testing procedure', sprintf(' (objDist = %.1fm)', objRange(odIndex))));
            % if (~isempty(model.savePath))
            %     plotpath = sprintf('%s/deltaMFasFktVerErr_objDist[%.1fm]', model.savePath, objRange(odIndex));
            %     saveas(gcf, plotpath, 'png');
            % end
        end
        xlabel('Vergence Error [deg]', 'FontSize', 12);
        ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
        title('\Delta MF(verg_{err}) response at Testing procedure');
        legend(num2str(objRange));
        if (~isempty(model.savePath))
            plotpath = sprintf('%s/deltaMFasFktVerErr', model.savePath);
            saveas(gcf, plotpath, 'png');
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

%this function generates anaglyphs of the large and small scale fovea and
%one of the two unpreprocessed gray scale images
function generateAnaglyphs(model, leftGray, rightGray, dsRatioL, dsRatioS, foveaL, foveaS)
    anaglyph = imfuse(leftGray, rightGray, 'falsecolor');
    imwrite(anaglyph,  sprintf('%s/anaglyph.png', model.savePath));

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
    imwrite(imresize(anaglyphL, 20), sprintf('%s/anaglyphLargeScale.png', model.savePath));
    largeScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(largeScaleView, 20), sprintf('%s/LargeScaleMontage.png', model.savePath));

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
    imwrite(imresize(anaglyphS, 8), sprintf('%s/anaglyphSmallScale.png', model.savePath));
    smallScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(smallScaleView, 8), sprintf('%s/smallScaleMontage.png', model.savePath));
end
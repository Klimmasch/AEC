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
    
    % imageSavePath = model.savePath;
    imageSavePath = '.'; 

    command = [0, 0];
    objRange = [model.objDistMin : 0.5 : model.objDistMax];
    testResult = zeros(size(objRange, 2), 7, 66);
    tmpResult1 = zeros(nStim, model.interval + 1);
    tmpResult2 = zeros(nStim, model.interval + 1);
    tmpResult3 = zeros(nStim, model.interval + 1);

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

    % minimal and maximal angle that can be reached by one-dimensional muscle commands
    angleMin = mfunction(indZero, 1);
    angleMax = mfunction(end, 1);

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
        angle = interp2(degrees.results_deg, cmd(1), cmd(2));   % interpolate in tabular
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
    
    %%% New renderer
    Simulator = OpenEyeSim('create');
    Simulator.initRenderer();
    % Simulator.reinitRenderer();
    
    imgRawLeft = uint8(zeros(240, 320, 3));
    imgRawRight = uint8(zeros(240, 320, 3));

    function [imLeft, imRight] = refreshImages(texture, vergAngle, objDist)
        Simulator.add_texture(1, texture);
        Simulator.set_params(1, vergAngle, objDist); %2-angle 3-distance

        result = Simulator.generate_left;
        result2 = Simulator.generate_right;

        imLeft=uint8(zeros(240, 320, 3));
        k=1;l=1;
        for i = 1:3:length(result)
            imLeft(k,l,1) = result(i);
            imLeft(k,l,2) = result(i+1);
            imLeft(k,l,3) = result(i+2);

            l=l+1;
            if (l>320)
                l=1;
                k=k+1;
            end
        end
    %     imLeft = COLOR;     %320x240 image

        imRight=uint8(zeros(240, 320, 3));
        k=1;l=1;
        for i = 1:3:length(result2)
            imRight(k,l,1) = result2(i);
            imRight(k,l,2) = result2(i+1);
            imRight(k,l,3) = result2(i+2);

            l=l+1;
            if (l>320)
                l=1;
                k=k+1;
            end
        end
    %     imRight = COLOR2;     %320x240 image
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

%                 [status, res] = system(sprintf('./checkEnvironment %s %d %d %s/leftTest.png %s/rightTest.png', ...
%                                                currentTexture, objRange(odIndex), angleNew, imageSavePath, imageSavePath));
%                 % abort execution if error occured
%                 if (status)
%                     sprintf('Error in checkEnvironment:\n%s', res)
%                     return;
%                 end
                
                for iter = 2 : model.interval + 1
                    % read input images and convert to gray scale
                    [imgRawLeft, imgRawRight] = refreshImages(currentTexture, -angleNew/2, objRange(odIndex));
%                     imgRawLeft = imread([imageSavePath '/leftTest.png']);
%                     imgRawRight = imread([imageSavePath '/rightTest.png']);
                    imgGrayLeft = .2989 * imgRawLeft(:,:,1) + .5870 * imgRawLeft(:,:,2) + .1140 * imgRawLeft(:,:,3);
                    imgGrayRight = .2989 * imgRawRight(:,:,1) + .5870 * imgRawRight(:,:,2) + .1140 * imgRawRight(:,:,3);

                    % imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [imagesSavePath '/anaglyph.png']);
                    % generateAnaglyphs(imageSavePath, imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS);

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
                        angleNew = getAngle2(command);           %resulting angle is used for both eyes
                    else
                        angleNew = angleNew + relativeCommand;
                        if (angleNew > angleMax || angleNew < angleMin)
                            angleNew = model.vergAngleMin + (model.vergAngleMax - model.vergAngleMin) * rand(1, 1);
                        end
                    end

                    % generate new view (two pictures) with new vergence angle
%                     [status, res] = system(sprintf('./checkEnvironment %s %d %d %s/leftTest.png %s/rightTest.png', ...
%                                                    currentTexture, objRange(odIndex), angleNew, imageSavePath, imageSavePath));
% 
%                     % abort execution if error occured
%                     if (status)
%                         sprintf('Error in checkEnvironment:\n%s', res)
%                         return;
%                     end

%                     [imgRawLeft, imgRawRight] = refreshImages(currentTexture, angleNew, objRange(odIndex));
                    
                    % temporary results
                    tmpResult1(stimulusIndex, iter) = angleDes - angleNew;
                    tmpResult2(stimulusIndex, iter) = relativeCommand;
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
            % delta_mf_t+1(vergAngle_t)
            errorbar(reshape(reshape(testResult(odIndex, :, 1 : 10), [7, 10])', [1, 7 * 10]), ...
                     reshape(reshape(testResult(odIndex, :, 24 : 33), [7, 10])', [1, 7 * 10]), ...
                     reshape(reshape(testResult(odIndex, :, 35 : 44), [7, 10])', [1, 7 * 10]), ...
                     'DisplayName', num2str(objRange(odIndex)), 'Marker', '*', 'MarkerSize', 3, ...
                     'color', [rand, rand, rand], 'LineWidth', 0.9, 'LineStyle', 'none');
            % l = legend('-DynamicLegend');
            legend('-DynamicLegend');
        end

        % adjust axis to actual response ranges + std deviation
        xmin = -4;
        xmax = 7;
        ymin = -0.1;
        ymax = 0.1;
        plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
        plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
        axis([xmin, xmax, ymin, ymax]);
        % l.Title.String = 'objDist [m]';
        % l.Title.FontSize = 12;
        xlabel('Vergence Error [deg]', 'FontSize', 12);
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
            errorbar(reshape(reshape(testResult(odIndex, :, 1 : 10), [7, 10])', [1, 7 * 10]), ...
                     reshape(reshape(testResult(odIndex, :, 46 : 55), [7, 10])', [1, 7 * 10]), ...
                     reshape(reshape(testResult(odIndex, :, 57 : 66), [7, 10])', [1, 7 * 10]), ...
                     'LineWidth', 0.9);
        end

        % adjust axis to actual response ranges + std deviation
        % xmin = -4;
        % xmax = 6;
        % ymin = -0.1;
        % ymax = 0.1;
        % plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
        % plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
        % axis([xmin, xmax, ymin, ymax]);
        xlabel('Vergence Error [deg]', 'FontSize', 12);
        ylabel('Value', 'FontSize', 12);
        title('Critic Value over different disparities');
        if (~isempty(model.savePath))
            plotpath = sprintf('%s/criticValvsVerErr', model.savePath);
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

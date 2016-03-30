%%% Model testing procedure
%@param model               respective model object to be tested
%@param randomizationSeed   randomization seed
%@param objRange            stimulus distance range
%@param repeat              # different stimuli for repetition
%%%
function testModel(model, randomizationSeed, objRange, repeat)
    rng(randomizationSeed);
    textureFile = 'Textures_New.mat';

    % Image process variables
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

    % Textures
    texturePath = sprintf('config/%s', textureFile);
    texture = load(texturePath);
    texture = texture.texture;
    nTextures = length(texture);

    degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'
    % metCosts = load('MetabolicCosts.mat');      %loads tabular for metabolic costs as 'results'

    %%% Helper function that maps muscle activities to resulting angle
    function [angle] = getAngle(command)
        cmd = (command * 10) + 1;                               % scale commands to table entries
        angle = interp2(degrees.results_deg, cmd(1), cmd(2));   % interpolate in tabular
    end

    %%% Helper function that maps muscle activities to resulting metabolic costs
    % function [tmpMetCost] = getMetCost(command)
    %     cmd = (command * 10) + 1;                               % scale commands to table entries
    %     tmpMetCost = interp2(metCosts.results, cmd(1), cmd(2)); % interpolate in tabular
    % end

    disZ = zeros(model.interval, size(objRange, 2), repeat);        % desired fixation distance
    fixZ = zeros(model.interval, size(objRange, 2), repeat);        % actual fixation distance
    vergerr = zeros(model.interval, size(objRange, 2), repeat);     % vergence error

    %%% Main execution loop
    for iter1 = 1 : repeat
        sprintf('Testing repetition = %d/%d', iter1, repeat)
        % pick random texture
        currentTexture = texture{(randi(nTextures, 1))};

        for iter2 = 1 : size(objRange, 2)

            % reset muscle activities to random values
            command = [0, 0];
%             command(2) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1, 1);    %only for one muscle
            % command(1) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1, 1);  %two muscles
            % command(2) = command(1);
            % command(2) = 0.1 * rand(1, 1); % random policy

            angleNew = getAngle(command) * 2;
            angleNew = randi(16,1); % for discrete Policy

            [status, res] = system(sprintf('./checkEnvironment %s %s %d %d left.png right.png %d', ...
                                           currentTexture, currentTexture, objRange(iter2), objRange(iter2), angleNew));

            % abort execution if error occured
            if (status)
                sprintf('Error in checkEnvironment:\n%s', res)
                return;
            end

            for iter3 = 1 : model.interval
                % read input images and convert to gray scale
                imgRawLeft = imread('left.png');
                imgRawRight = imread('right.png');
                imgGrayLeft = .2989 * imgRawLeft(:,:,1) + .5870 * imgRawLeft(:,:,2) + .1140 * imgRawLeft(:,:,3);
                imgGrayRight = .2989 * imgRawRight(:,:,1) + .5870 * imgRawRight(:,:,2) + .1140 * imgRawRight(:,:,3);
                
                generateAnaglyphs(imgGrayLeft, imgGrayRight, dsRatioL, dsRatioS, foveaL, foveaS, model.savePath);
                
                % Image patch generation: left{small scale, large scale}, right{small scale, large scale}
                [patchesLeftSmall] = preprocessImage(imgGrayLeft, foveaS, dsRatioS, patchSize, columnIndS);
                [patchesLeftLarge] = preprocessImage(imgGrayLeft, foveaL, dsRatioL, patchSize, columnIndL);
                [patchesRightSmall] = preprocessImage(imgGrayRight, foveaS, dsRatioS, patchSize, columnIndS);
                [patchesRightLarge] = preprocessImage(imgGrayRight, foveaL, dsRatioL, patchSize, columnIndL);

                % Image patches matrix (input to model)
                currentView = {[patchesLeftLarge; patchesRightLarge] [patchesLeftSmall; patchesRightSmall]};

                % Generate input feature vector from current images
                [feature, ~, errorTotal, errorLarge, errorSmall] = model.generateFR(currentView);

                %%% Feedback
                % Absolute command feedback # concatination
%                 feature = [feature; command(2) * model.lambdaMuscleFB];

                %%% Calculate metabolic costs
                % metCost = getMetCost(command) * 2;

                %%% Action
                relativeCommand = model.rlmodel.softmaxAct(feature);

                % add the change in muscle Activities to current ones
                % command = command + relativeCommand';     %two muscels
%                 command(2) = command(2) + relativeCommand;  %one muscel
%                 command = checkCmd(command);                %restrain motor commands to [0,1]
%                 angleNew = getAngle(command) * 2;           %resulting angle is used for both eyes
                angleNew = angleNew + relativeCommand;
                if angleNew > 16 || angleNew < 0.01
                    angleNew = randi(16,1);
                end

                % generate new view (two pictures) with new vergence angle
                [status, res] = system(sprintf('./checkEnvironment %s %s %d %d left.png right.png %d', ...
                                       currentTexture, currentTexture, objRange(iter2), objRange(iter2), angleNew));

                % abort execution if error occured
                if (status)
                    sprintf('Error in checkEnvironment:\n%s', res)
                    return;
                end

                %%% Track results
                % compute desired vergence command, disparity and vergence error
                fixDepth = (model.baseline / 2) / tand(angleNew / 2);           %fixation depth [m]
                angleDes = 2 * atand(model.baseline / (2 * objRange(iter2)));   %desired vergence [deg]
                anglerr = angleDes - angleNew;                                  %vergence error [deg]
                % disparity = 2 * model.focalLength * tand(anglerr / 2);        %current disp [px]

                disZ(iter3, iter2, iter1) = objRange(iter2);
                fixZ(iter3, iter2, iter1) = fixDepth;
                vergerr(iter3, iter2, iter1) = anglerr;
            end
        end
    end

    %%% Vergence error dynamics over 10 iterations
    % mean over repetitions (stimuli) and object distances
    figure;
    hold on;
    grid on;
    % plot(1 : model.interval, mean(mean(vergerr, 3), 2), 'color', [1, 0.549, 0], 'LineWidth', 0.8);
    errorbar([1 : model.interval], mean(mean(vergerr, 3), 2), std(std(vergerr, 0, 3), 0, 2), 'color', [1, 0.549, 0], 'LineWidth', 0.8);

    xlabel('Iteration step', 'FontSize', 12);
    ylabel('Vergence Error [deg]', 'FontSize', 12);
    title('Average Vergence Error over Trial (Testing)');
    plotpath = sprintf('%s/AvgVergErrOverTrial', model.savePath);
    saveas(gcf, plotpath, 'png');

    %%% Vergence error dynamics over whole testing procedure
    figure;
    hold on;
    grid on;
    plot(1 : model.interval * size(objRange, 2) * repeat, reshape(disZ, [numel(disZ), 1]), 'LineWidth', 1.3);
    plot(1 : model.interval * size(objRange, 2) * repeat, reshape(fixZ, [numel(fixZ), 1]), 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.3);
    xlabel('Iteration #', 'FontSize', 12);
    ylabel('Fixation [m]', 'FontSize', 12);
    l = legend('desired', 'actual');
    l.Orientation = 'horizontal';
    l.Location = 'southoutside';
    title('Fixation Distance at Testing Procedure');
    plotpath = sprintf('%s/vergenceAngleTesting', model.savePath);
    saveas(gcf, plotpath, 'png');
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

function generateAnaglyphs(leftGray, rightGray, dsRatioL, dsRatioS, foveaL, foveaS, savePath)
    anaglyph = imfuse(leftGray, rightGray, 'falsecolor');
    imwrite(anaglyph, 'anaglyph.png');

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
    imwrite(imresize(anaglyphL, 20), strcat(savePath, 'anaglyphLargeScale.png'));
    largeScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(largeScaleView, 20), strcat(savePath,'LargeScaleMontage.png'));

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
    imwrite(imresize(anaglyphS, 8), strcat(savePath,'anaglyphSmallScale.png'));
    smallScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(smallScaleView, 8), strcat(savePath,'smallScaleMontage.png'));
end
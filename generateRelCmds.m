%This function goes through all object distances and averages the motor command
% and reconstruction error for different vergence errors.
%TODO what about muscle innervations?
% example usage: averageForVergenceErrors(model, objRange=[0.5, 1, 2], vergRange =
%   -5:0.5:5, repeat = 10)
% takes the model, positiones random textures at the distances of 0.5, 1 and 2
% meters, than goes through the range of vergence "Errors" and generates
% muscle commands and repeats that 10 times
%TODO: analyze whole function
function responseResults = generateRelCmds(model, objRange, vergRange, repeat)
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

    % preparing Textures
    texture = load(['config/' model.textureFile]);
    texture = texture.texture;
    nTextures = length(texture);

    vergErrs = [];
    relCmds = [];
    recErrs = [];
    recErrsSmall = [];
    recErrsLarge = [];

    [~, numDists] = size(objRange);
    % avgCmd = zeros(numDists,1);

    % vergences = vergRange(1):0.1:vergRange(2);
    [~, numVergs] = size(vergRange);
    % avgVergErr = zeros(numVergs,1);

    degrees = load('Degrees.mat');
    resolution = 10001;
    approx = spline(1:11, degrees.results_deg(:, 1));

    xValPos = ppval(approx, 1:0.001:11)';
    yValPos = linspace(0, 1, resolution)';

    xValNeg = flipud(ppval(approx, 1:0.001:11)' * -1);
    yValNeg = linspace(-1, 0, resolution)';

    mf = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
    dmf = diff(mf(1:2, 1)); % delta in angle
    indZero = find(mf(:, 2) == 0); % MF == 0_index

    sprintf('starting to generate vergence commands for different vergence errors ...')
    for rep = 1:repeat
        for objDist = 1:numDists
            angleDes = 2 * atand(model.baseline / (2 * objRange(objDist)));

            for verg = 1:numVergs
                currentTexture = texture{(randi(nTextures, 1))}; %random picture for every iteration
                %generate two new pictures
                [status, res] = system(sprintf('./checkEnvironment %s %d %d left.png right.png', ...
                                               currentTexture, objRange(objDist), angleDes + vergRange(verg)));

                % Abort execution if error occured
                if (status)
                    sprintf('Error in checkEnvironment:\n%s', res)
                    return;
                end

                % Read input images and convert to gray scale
                imgRawLeft = imread('left.png');
                imgRawRight = imread('right.png');
                imgGrayLeft = .2989 * imgRawLeft(:,:,1) + .5870 * imgRawLeft(:,:,2) + .1140 * imgRawLeft(:,:,3);
                imgGrayRight = .2989 * imgRawRight(:,:,1) + .5870 * imgRawRight(:,:,2) + .1140 * imgRawRight(:,:,3);

                % anaglyph = stereoAnaglyph(imgGrayLeft, imgGrayRight);
                % imwrite(anaglyph, 'anaglyph.png');

                % Image patch generation: left{small scale, large scale}, right{small scale, large scale}
                [patchesLeftSmall] = preprocessImage(imgGrayLeft, foveaS, dsRatioS, patchSize, columnIndS);
                [patchesLeftLarge] = preprocessImage(imgGrayLeft, foveaL, dsRatioL, patchSize, columnIndL);
                [patchesRightSmall] = preprocessImage(imgGrayRight, foveaS, dsRatioS, patchSize, columnIndS);
                [patchesRightLarge] = preprocessImage(imgGrayRight, foveaL, dsRatioL, patchSize, columnIndL);

                % Image patches matrix (input to model)
                currentView = {[patchesLeftLarge; patchesRightLarge] [patchesLeftSmall; patchesRightSmall]};

                % Generate input feature vector from current images
                [feature, ~, errorTotal, errorLarge, errorSmall] = model.generateFR(currentView);

                indTemp = find(mf(:, 1) <= angleDes + vergRange(verg) + dmf & mf(:, 1) >= angleDes + vergRange(verg) - dmf);
                if (size(indTemp, 1) < 1)
                    indTemp = indZero;
                end
                feature = [feature; mf(indTemp(1), 2)];
                relCmd = model.rlmodel.softmaxAct(feature);

                %Traking variables
                vergErrs = [vergErrs; vergRange(verg)];
                relCmds = [relCmds; relCmd];
                recErrs = [recErrs; errorTotal];
                recErrsLarge = [recErrsLarge; errorLarge];
                recErrsSmall = [recErrsSmall; errorSmall];
            end
        end
        % sprintf('number of repetitions: %d/%d done', rep, repeat)
    end
    responseResults = struct('relCmds', relCmds, 'vergErrs', vergErrs, 'recErrs', recErrs, 'recErrsLarge', recErrsLarge, 'recErrsSmall', recErrsSmall);
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

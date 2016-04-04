%%% Model testing procedure
%@param model               respective model object to be tested
%@param randomizationSeed   randomization seed
%@param objRange            stimulus distance range
%@param vergRange           vergence range
%@param repeat              [#repetition test, #repetition deltaMF and recErr]
%@param testRandObjRange    whether the object ranges shall be randomized at testing procedure
%%%
function testModel(model, randomizationSeed, objRange, vergRange, repeat, testRandObjRange)
    rng(randomizationSeed);
    textureFile = model.textureFile;

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

    % Prepare Textures
    texture = load(['config/' model.textureFile]);
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

    disZ = zeros(model.interval, size(objRange, 2), repeat(1));        % desired fixation distance
    fixZ = zeros(model.interval, size(objRange, 2), repeat(1));        % actual fixation distance
    vergerr = zeros(model.interval, size(objRange, 2), repeat(1));     % vergence error
    tmpObjRange = 0;

    %%% Average VergErr over Trial loop
    for iter1 = 1 : repeat(1)
        sprintf('Testing repetition = %d/%d', iter1, repeat(1))
        % pick random texture
        currentTexture = texture{(randi(nTextures, 1))};

        for iter2 = 1 : size(objRange, 2)
            % reset muscle activities to random values
            if (model.rlmodel.continuous == 1)
                command = [0, 0];
                command(2) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1, 1); %only for one muscle
                angleNew = getAngle(command) * 2;
            else
                angleNew = randi(16, 1);
            end

            % Object distance = random/determined
            if (testRandObjRange == 1)
                tmpObjRange = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);
            else
                tmpObjRange = objRange(iter2);
            end

            [status, res] = system(sprintf('./checkEnvironment %s %d %d left.png right.png', ...
                                           currentTexture, tmpObjRange, angleNew));
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

                %%% Feedback
                % Absolute command feedback # concatination
                if (model.rlmodel.continuous == 1)
                    feature = [feature; command(2) * model.lambdaMuscleFB];
                end

                %%% Calculate metabolic costs
                % metCost = getMetCost(command) * 2;

                %%% Action
                relativeCommand = model.rlmodel.softmaxAct(feature);

                % add the change in muscle Activities to current ones
                if (model.rlmodel.continuous == 1)
                    % command = command + relativeCommand';     %two muscels
                    command(2) = command(2) + relativeCommand;  %one muscel
                    command = checkCmd(command);                %restrain motor commands to [0,1]
                    angleNew = getAngle(command) * 2;           %resulting angle is used for both eyes
                else
                    angleNew = angleNew + relativeCommand;
                    if (angleNew > 16 || angleNew < 0.01)
                        angleNew = randi(16,1);
                    end
                end

                % generate new view (two pictures) with new vergence angle
                [status, res] = system(sprintf('./checkEnvironment %s %d %d left.png right.png', ...
                                               currentTexture, tmpObjRange, angleNew));

                % abort execution if error occured
                if (status)
                    sprintf('Error in checkEnvironment:\n%s', res)
                    return;
                end

                %%% Track results
                % compute desired vergence command, disparity and vergence error
                fixDepth = (model.baseline / 2) / tand(angleNew / 2);           %fixation depth [m]
                angleDes = 2 * atand(model.baseline / (2 * tmpObjRange));   %desired vergence [deg]
                anglerr = angleDes - angleNew;                                  %vergence error [deg]
                % disparity = 2 * model.focalLength * tand(anglerr / 2);        %current disp [px]

                disZ(iter3, iter2, iter1) = tmpObjRange;
                fixZ(iter3, iter2, iter1) = fixDepth;
                vergerr(iter3, iter2, iter1) = abs(anglerr);
            end
        end
    end

    %%% Vergence error dynamics over 10 iterations
    % mean over repetitions (stimuli) and object distances
    figure;
    hold on;
    grid on;
    errorbar([1 : model.interval], mean(mean(vergerr, 3), 2), std(std(vergerr, 0, 3), 0, 2), 'color', [1, 0.549, 0], 'LineWidth', 0.8);
    axis([0, 11, 0, 3]);
    xlabel('Iteration step', 'FontSize', 12);
    ylabel('Vergence Error [deg]', 'FontSize', 12);
    title('Average Vergence Error over Trial (Testing)');
    plotpath = sprintf('%s/AvgVergErrOverTrial', model.savePath);
    saveas(gcf, plotpath, 'png');

    %%% Vergence error dynamics over whole testing procedure
    figure;
    hold on;
    grid on;
    plot(1 : model.interval * size(objRange, 2) * repeat(1), reshape(disZ, [numel(disZ), 1]), 'LineWidth', 1.3);
    plot(1 : model.interval * size(objRange, 2) * repeat(1), reshape(fixZ, [numel(fixZ), 1]), 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.3);
    xlabel('Iteration #', 'FontSize', 12);
    ylabel('Fixation [m]', 'FontSize', 12);
    l = legend('desired', 'actual');
    l.Orientation = 'horizontal';
    l.Location = 'southoutside';
    title('Fixation Distance at Testing Procedure');
    plotpath = sprintf('%s/vergenceAngleTesting', model.savePath);
    saveas(gcf, plotpath, 'png');

    %%% Relative commands loop
    vergErrs = [];
    relCmds = [];
    recErrs = [];
    recErrsSmall = [];
    recErrsLarge = [];
    criticValue = [];

    resolution = 10001;
    approx = spline(1:11, degrees.results_deg(:, 1));

    xValPos = ppval(approx, 1:0.001:11)';
    yValPos = linspace(0, 1, resolution)';

    xValNeg = flipud(ppval(approx, 1:0.001:11)' * -1);
    yValNeg = linspace(-1, 0, resolution)';

    mf = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
    mf(:, 1) = mf(:, 1) * 2;        % angle for two eyes
    dmf = diff(mf(1:2, 1));         % delta in angle
    indZero = find(mf(:, 2) == 0);  % MF == 0_index

    sprintf('starting to generate vergence commands for different vergence errors ...')
    for iter1 = 1 : repeat(2)
        sprintf('RelCmd repetition = %d/%d', iter1, repeat(2))

        for objDist = 1 : size(objRange, 2)
            angleDes = 2 * atand(model.baseline / (2 * objRange(objDist)));

            for verg = 1 : size(vergRange, 2)
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

                if (model.rlmodel.continuous == 1)
                    feature = [feature; mf(indTemp(1), 2)];
                    value = model.rlmodel.CCritic.v_ji * feature;
                else
                    value = model.rlmodel.Weights{2,1} * feature;
                end

                relCmd = model.rlmodel.softmaxAct(feature);

                %Tracking variables
                relCmds = [relCmds; relCmd];
                vergErrs = [vergErrs; vergRange(verg)];
                recErrs = [recErrs; errorTotal];
                recErrsLarge = [recErrsLarge; errorLarge];
                recErrsSmall = [recErrsSmall; errorSmall];
                criticValue = [criticValue; value];
            end
        end
    end
    responseResults = struct('vergErrs', vergErrs, 'relCmds', relCmds, 'recErrs', recErrs, 'recErrsLarge', recErrsLarge, ...
                             'recErrsSmall', recErrsSmall, 'criticValue', criticValue, 'objRange', objRange, 'vergRange', vergRange);

    deltaMFplotGenDist(model, responseResults);
    recErrPlotGenDist(model, responseResults);
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

%% plot & save delta_MF(Vergence_error)
%@param responseResults model's responses to generated vergErr and objDist
function deltaMFplotGenDist(model, responseResults)
    % desired_angle_min @ 2m = 2 * atand(baseline / (2 * 2)) = 1.6042
    % desired_angle_max @ 0.5m = 2 * atand(baseline / (2 * 0.5)) = 6.4104
    % actual_distance_min = (baseline / 2) / tand(results_deg(11,1)*2 / 2) = 0.0389 [m]
    % actual_distance_max = (baseline / 2) / tand(results_deg(1,1)*2 / 2) = 3.2219 [m]
    % angle_min @ actual_distance_max = results_deg(1,1) * 2 = 0.9958 deg
    % angle_max @ actual_distance_min = results_deg(11,1) * 2 = 71.5164 deg
    % verg_err_min = desired_angle_min - angle_max = 1.6042 - 71.5164 = -69.9122
    % Verg_err_max = desired_angle_max - angle_min = 6.4104 - 0.9958 = 5.4146

    degrees = load('Degrees.mat');
    resolution = 10001;
    approx = spline(1:11, degrees.results_deg(:, 1));

    xValPos = ppval(approx, 1:0.001:11)';
    yValPos = linspace(0, 1, resolution)';

    %TODO: no need for negative values, check it!
    xValNeg = flipud(ppval(approx, 1:0.001:11)' * -1);
    yValNeg = linspace(-1, 0, resolution)';

    % calculate muscle function :=  mf(vergence_angle) = muscle force [single muscle]
    mf = [xValNeg(1 : end - 1), yValNeg(1 : end - 1); xValPos, yValPos];
    dmf = diff(mf(1:2, 1)); % delta in angle
    indZero = find(mf(:, 2) == 0); % MF == 0_index
    %TODO: make a function out of this
    indMaxFix = find(mf(:, 1) <= model.desiredAngleMin + dmf & mf(:, 1) >= model.desiredAngleMin - dmf); % MF(desiredAngleMin)_index
    indMinFix = find(mf(:, 1) <= model.desiredAngleMax + dmf & mf(:, 1) >= model.desiredAngleMax - dmf); % MF(desiredAngleMax)_index

    % perfect_response := [max_fixation_x, max_fixation_y, min_fixation_x, min_fixation_y]
    % x = vergenceError, y = deltaMuscelForce
    perfectResponseMaxFix = [(mf(indMaxFix, 1) - flipud(mf(indMaxFix : end, 1))) * 2, ...
                             (mf(indMaxFix, 2) - flipud(mf(indMaxFix : end, 2))); ...
                             (mf(indMaxFix, 1) - flipud(mf(indZero : indMaxFix - 1, 1))) * 2, ...
                             (mf(indMaxFix, 2) - flipud(mf(indZero : indMaxFix - 1, 2)))];

    perfectResponseMinFix = [(mf(indMinFix, 1) - flipud(mf(indMinFix : end, 1))) * 2, ...
                             (mf(indMinFix, 2) - flipud(mf(indMinFix : end, 2))); ...
                             (mf(indMinFix, 1) - flipud(mf(indZero : indMinFix - 1, 1))) * 2, ...
                             (mf(indMinFix, 2) - flipud(mf(indZero : indMinFix - 1, 2)))];

    perfectResponse = [perfectResponseMaxFix, perfectResponseMinFix];

    %%% generating fixed distances
    % actualResponse = [vergErrs, relCmds];
    actualResponse = [responseResults.vergErrs, responseResults.relCmds];
    nVal = size(responseResults.vergRange, 2); % #bins of statistics

    actualResponse = sortrows(actualResponse);
    deltaVergErr = (abs(actualResponse(1, 1)) + abs(actualResponse(end, 1))) / nVal;
    % actualResponseStat = [index_x = vergence_error angle, mean_muscle_force, std_muscle_force]
    actualResponseStat = zeros(nVal, 3);

    %TODO: whole loop can be simplified here at generated vergErrors, no need to search vergErrors
    for i = 1 : nVal
        %TODO: shift responseResults.vergRange(i) by deltaVergErr/2...maybe not neccessary
        actualResponseStat(i, 1) = responseResults.vergRange(i);
        actualResponseStat(i, 2) = mean(actualResponse(find(actualResponse(:, 1) >= actualResponse(1, 1) + (i - 1) * deltaVergErr ...
                                                            & actualResponse(:, 1) <= actualResponse(1, 1) + i * deltaVergErr), 2));
        actualResponseStat(i, 3) = std(actualResponse(find(actualResponse(:, 1) >= actualResponse(1, 1) + (i - 1) * deltaVergErr ...
                                                           & actualResponse(:, 1) <= actualResponse(1, 1) + i * deltaVergErr), 2));
    end
    actualResponseStat(isnan(actualResponseStat(:, 2)), :) = []; % drop NaN elements

    figure;
    hold on;
    grid on;
    % perfect response to vergence error
    plot(perfectResponse(:, 1), perfectResponse(:, 2), 'color', [0.5882, 0.9608, 0], 'LineWidth', 1.3);
    plot(perfectResponse(:, 3), perfectResponse(:, 4), 'color', [0, 0.5882, 0.9608], 'LineWidth', 1.3);
    % error bars of actual response of the model
    errorbar(actualResponseStat(:, 1), actualResponseStat(:, 2), actualResponseStat(:, 3),'color', [1, 0.5098, 0.1961], 'LineWidth', 0.9);
    % actual response of the model
    plot(actualResponseStat(:, 1), actualResponseStat(:, 2),'color', [1, 0.0784, 0], 'LineWidth', 1.3);
    l = legend('perfect (fixDist_{max})', 'perfect (fixDist_{min})', 'actual');
    l.FontSize = 7;
    l.Orientation = 'horizontal';
    l.Location = 'southoutside';
    % adjust axis to actual response ranges + std deviation
    xmin = min(actualResponseStat(:, 1)) * 1.1;
    xmax = max(actualResponseStat(:, 1)) * 1.1;
    if (xmin >= xmax)
        xmin = -5.5;
        xmax = 5.5;
    end
    ymin = min([-0.1, (min(actualResponseStat(:, 2)) - max(actualResponseStat(:, 3))) * 1.1]);
    ymax = max([max(perfectResponse(:, 4)) * 1.1, (max(actualResponseStat(:, 2)) + max(actualResponseStat(:, 3))) * 1.1]);
    if (ymin >= ymax)
        ymin = -0.1;
        ymax = 0.1;
    end
    plot([xmin, xmax], [0, 0], 'k', 'LineWidth', 0.1);
    plot([0, 0], [ymin, ymax], 'k', 'LineWidth', 0.1);
    axis([xmin, xmax, ymin, ymax]);
    xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
    ylabel('\Delta MF \in [-1, 1]', 'FontSize', 12);
    title('\Delta MF(verg_{err}) response at Testing procedure');
    if (~isempty(model.savePath))
        plotpath = sprintf('%s/deltaMFasFktVerErrGenDist_[%fm,%fm]', model.savePath, responseResults.objRange(1), responseResults.objRange(end));
        saveas(gcf, plotpath, 'png');
    end
end

%% plot & save reconstructionErr(Vergence_error)
%@param responseResults model's responses to generated vergErr and objDist
function recErrPlotGenDist(model, responseResults)
    %plotting the resonstruction error of basis functions over
    %different disparities
    recResponse = [responseResults.vergErrs, responseResults.recErrs];
    recResponseLarge = [responseResults.vergErrs, responseResults.recErrsLarge];
    recResponseSmall = [responseResults.vergErrs, responseResults.recErrsSmall];
    nVal = size(responseResults.vergRange, 2); % #bins of statistics

    %calculate mean and std of reconstruction error
    tmpRsp = sortrows(recResponse);
    deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nVal;
    % tmp = [index_x = vergence_error angle, mean_recError, std_recError]
    tmp = zeros(nVal, 3);

    for i = 1:nVal
        tmp(i, 1) = responseResults.vergRange(i);
        tmp(i, 2) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                     & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
        tmp(i, 3) = std(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                    & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
    end
    recErrs = tmp;
    recErrs(isnan(recErrs(:, 2)), :) = []; % drop NaN elements

    %calculate mean and std of large scale reconstruction error
    tmpRsp = sortrows(recResponseLarge);
    deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nVal;
    % tmp = [index_x = vergence_error angle, mean_recError, std_recError]
    tmp = zeros(nVal, 3);

    for i = 1:nVal
        tmp(i, 1) = responseResults.vergRange(i);
        tmp(i, 2) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                     & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
        tmp(i, 3) = std(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                    & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
    end
    recErrsLarge = tmp;
    recErrsLarge(isnan(recErrsLarge(:, 2)), :) = []; % drop NaN elements

    %calculate mean and std of small scale reconstruction error
    tmpRsp = sortrows(recResponseSmall);
    deltaVergErr = (abs(tmpRsp(1, 1)) + abs(tmpRsp(end, 1))) / nVal;
    % tmp = [index_x = vergence_error angle, mean_recError, std_recError]
    tmp = zeros(nVal, 3);

    for i = 1:nVal
        tmp(i, 1) = responseResults.vergRange(i);
        tmp(i, 2) = mean(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                     & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
        tmp(i, 3) = std(tmpRsp(find(tmpRsp(:, 1) >= tmpRsp(1, 1) + (i - 1) * deltaVergErr ...
                                    & tmpRsp(:, 1) <= tmpRsp(1, 1) + i * deltaVergErr), 2));
    end
    recErrsSmall = tmp;
    recErrsSmall(isnan(recErrsSmall(:, 2)), :) = []; % drop NaN elements

    figure;
    hold on;
    grid on;
    % error bars of reconstruction of the model
    errorbar(recErrs(:, 1), recErrs(:, 2), recErrs(:, 3), 'LineWidth', 0.9); %'color', [1, 0.5098, 0.1961],
    errorbar(recErrsLarge(:, 1), recErrsLarge(:, 2), recErrsLarge(:, 3), 'LineWidth', 0.9);%'color', [1, 0.5098, 0.1961],
    errorbar(recErrsSmall(:, 1), recErrsSmall(:, 2), recErrsSmall(:, 3), 'LineWidth', 0.9);%'color', [1, 0.5098, 0.1961],
    % reconstruction of the model
    % plot(recErrs(:, 1), recErrs(:, 2),'color', [1, 0.0784, 0], 'LineWidth', 1.3);
    l = legend('reconstruction Error', 'small scale recErr', 'large scale recErr');
    l.FontSize = 7;
    l.Orientation = 'horizontal';
    l.Location = 'southoutside';
    xlabel(sprintf('Vergence Error [deg] (bin size = %.2f°)', deltaVergErr), 'FontSize', 12);
    ylabel('resonstruction Error', 'FontSize', 12);
    title(sprintf('Reconstruction Error over different disparities\nobject distances: [%s]', num2str(responseResults.objRange)));

    if (~ isempty(model.savePath))
        plotpath = sprintf('%s/recErrVsVerErrGenDist_[%fm,%fm]', model.savePath, responseResults.objRange(1), responseResults.objRange(end));
        saveas(gcf, plotpath, 'png');
    end
end

%this function generates anaglyphs of the large and small scale fovea and
%one of the two unpreprocessed gray scale images
% TODO: adjust the sizes of the montage view
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
    imwrite(imresize(anaglyphS, 16), sprintf('%s/anaglyphSmallScale.png', model.savePath));
    smallScaleView = imfuse(imgLeftL, imgRightL, 'montage');
    imwrite(imresize(smallScaleView, 8), sprintf('%s/smallScaleMontage.png', model.savePath));
end

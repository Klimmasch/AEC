%%% Model testing procedure that generates movies of all scales in the
%%% model folder. For best results, run it headless!
%@param randomizationSeed   randomization seed
%@param objRange            stimulus distance range
%@param randObjRange        whether the entries in objRange should be random
%@param nStimuli            number of input stimuli
%@param randStimuli         whether the stimuli shall in randomized order
%@param reinitRenderer      choose 0 if it's the first time you start the renderer, other values are usually for debugging purposes
%%%
function generateMovie(objRange, nStimuli, reinitRenderer)

    model = load('/home/aecgroup/aecdata/Results/CriticLR vs ActorLR/16-11-02_2000000iter_1_regul_1e-04_actorRL_[1.00-0.00]_criticRL_[1.00-0.00]/model.mat');
    model = model.model;

    randomizationSeed = 3;
    randObjRange = 1;
    randStimuli = 1;
    testInterval = model.interval*2;
    rng(randomizationSeed);

%    textureFile = 'Textures_vanHaterenTrain';
    % textureFile = 'Textures_vanHaterenTest';
%     textureFile = 'Textures_vanHaterenTestMinimal';
%     textureFile = 'Textures_celine'
    textureFile = 'Textures_mcgillManMadeTest(jpg)'

    imagesSavePath = sprintf('%s/movies', model.savePath);
    timeStamp = datestr(now, 'dd-mmm-yyyy_HH:MM:SS'); % used as a tag for the movies
    mkdir(imagesSavePath);
    frameRate = 3; % frames per second in the resulting video »| 1.5 orig
    markScales = uint8(1); % draws rectangles inside the anaglyphs to indicate foveal regions
    set(0,'DefaulttextInterpreter','none'); % prevents underscores in image files to be interpreted as subscripts



    % Image processing variables
%     patchSize = 8;
%
%     dsRatioL = model.scmodel_Large.Dsratio; %downsampling ratio (Large scale) | original 8
%     dsRatioS = model.scmodel_Small.Dsratio; %downsampling ratio (Small scale) | original 2
%
%     % fovea = [128 128];
%     foveaL = patchSize + patchSize ^ 2 / 2 ^ log2(dsRatioL); %fovea size (Large scale) | 16
%     foveaS = patchSize + patchSize ^ 2 / 2 ^ log2(dsRatioS); %fovea size (Small scale) | 40
%
%     stOvL = patchSize / dsRatioL; %steps of overlap in the ds image | 1
%     stOvS = patchSize / dsRatioS; %steps of overlap in the ds image | 4
%
%     ncL = foveaL - patchSize + 1; %number of patches per column (slide of 1 px) | 9
%     ncS = foveaS - patchSize + 1; %number of patches per column (slide of 1 px) | 33
%
%     % Prepare index matricies for image patches
%     columnIndL = [];
%     for kc = 1:stOvL:ncL
%         tmpInd = (kc - 1) * ncL + 1 : stOvL : kc * ncL;
%         columnIndL = [columnIndL tmpInd];
%     end
%     columnIndS = [];
%     for kc = 1:stOvS:ncS
%         tmpInd = (kc - 1) * ncS + 1 : stOvS : kc * ncS;
%         columnIndS = [columnIndS tmpInd];
%     end

    % Prepare Textures
    texture = load(['config/' textureFile]);
    texture = texture.texture;
    nTextures = length(texture);
    if nStimuli > nTextures
        sprintf('%s only contains %d stimuli, but I will use them all!', textureFile, nTextures)
        nStimuli = nTextures;
    end

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

%     %%% Perfect Response function
%     indMaxFix = find(mfunction(:, 1) <= model.vergAngleMin + dmf & mfunction(:, 1) >= model.vergAngleMin - dmf); % MF(vergAngleMin)_index
%     indMaxFix = indMaxFix(1);
%     indMinFix = find(mfunction(:, 1) <= model.vergAngleMax + dmf & mfunction(:, 1) >= model.vergAngleMax - dmf); % MF(vergAngleMax)_index
%     indMinFix = indMinFix(1);

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

    %%% Helper function for calculating {objDist} -> {maxVergErr}
    function vergErrMax = getVergErrMax(objDist)
        % correct vergence angle for given object distance
        angleCorrect = 2 * atand(model.baseline / (2 * objDist));
        vergErrMax = angleCorrect - angleMin;
    end

    %%% Helper function that maps muscle activities to resulting angle
    function [angle] = getAngle(command)
        cmd = (command * 10) + 1;                               % calculate tabular index
        angle = interp2(degrees.results_deg, cmd(1), cmd(2), 'spline');   % interpolate in tabular
    end


    %%% same function as getAngle but with better performance
%     function angle = getAngle2(command)
%         angleIndex = find(mfunction(:, 2) <= command(2) + dmf2 & mfunction(:, 2) >= command(2) - dmf2);
%         angle = mfunction(angleIndex, 1);
%         angle = angle(ceil(length(angle) / 2));
%     end

    %%% Helper function that maps muscle activities to resulting metabolic costs
    % function [tmpMetCost] = getMetCost(command)
    %     cmd = (command * 10) + 1;                               % scale commands to table entries
    %     tmpMetCost = interp2(metCosts.results, cmd(1), cmd(2)); % interpolate in tabular
    % end

%     disZtest = zeros(model.interval, size(objRange, 2), repeat(1));        % desired fixation distance
%     fixZtest = zeros(model.interval, size(objRange, 2), repeat(1));        % actual fixation distance
%     vergErrTest = zeros(model.interval, size(objRange, 2), repeat(1));     % vergence error
%     tmpObjRange = 0;

    % minimal and maximal angle that can be reached by one-dimensional
    % muscle commands
    angleMin = getAngle([0, 0]) * 2;
    angleMax = getAngle([0, 1]) * 2;

    %%% New renderer
%     simulator = OpenEyeSim('create');
%    simulator = OpenEyeSimV4('create');
%     if reinitRenderer
%         simulator.reinitRenderer();
%     else
%         simulator.initRenderer();
%     end
%     imgRawLeft = uint8(zeros(240, 320, 3));
%     imgRawRight = uint8(zeros(240, 320, 3));
%     imgGrayLeft = uint8(zeros(240, 320, 3));
%     imgGrayRight = uint8(zeros(240, 320, 3));

    % Image patches cell array (input to model)
%     currentView = cell(1, length(model.scModel));

%     function refreshImages(texture, vergAngle, objDist)
% 
%         simulator.add_texture(1, texture);
%         simulator.set_params(1, vergAngle, objDist, 0, 3);
% 
%         result1 = simulator.generate_left();
%         result2 = simulator.generate_right();
% 
%         imgRawLeft = permute(reshape(result1, ...
%                                      [size(imgRawLeft, 3), ...
%                                       size(imgRawLeft, 2), ...
%                                       size(imgRawLeft, 1)]), ...
%                                      [3, 2, 1]);
% 
%         imgRawRight = permute(reshape(result2, ...
%                                       [size(imgRawRight, 3), ...
%                                        size(imgRawRight, 2), ...
%                                        size(imgRawRight, 1)]), ...
%                                       [3, 2, 1]);
% 
%         % convert images to gray scale
%         imgGrayLeft = 0.2989 * imgRawLeft(:, :, 1) + 0.5870 * imgRawLeft(:, :, 2) + 0.1140 * imgRawLeft(:, :, 3);
%         imgGrayRight = 0.2989 * imgRawRight(:, :, 1) + 0.5870 * imgRawRight(:, :, 2) + 0.1140 * imgRawRight(:, :, 3);
%     end

    t = 1;
    %%% Average VergErr over Trial loop
    for stimulus = 1 : nStimuli
        % pick texture
        if randStimuli
            currentTexture = texture{(randi(nTextures, 1))}
        else
            currentTexture = texture{stimulus}
        end
        sprintf('Iteration %d of %d', stimulus, nStimuli)

        for objDist = 1 : size(objRange, 2)

            vergErrMax = getVergErrMax(objRange(objDist));
            if vergErrMax > 2
                vergErrMax = 2;
            end
%             vseRange = [-3, -2, -1, linspace(0, vergErrMax, 4)];
            vseRange = [-2, 0, vergErrMax];
%             vseRange = [-3, 0, 3];
            if stimulus == nStimuli
                vseRange = [vseRange, vseRange(end)]; % this tries to tackle a bug with removing the last few images from the video
            end

            for vseIndex = 1 : size(vseRange, 2)

%                 if (model.rlmodel.continuous == 1)
%                 command = [0; 0];
%                 [command(2), angleNew] = getMF(objRange(objDist), vseRange(vseIndex));
                [command, angleNew] = getMF2(objRange(objDist), vseRange(vseIndex));
%                     command(2) = model.muscleInitMin + (model.muscleInitMax - model.muscleInitMin) * rand(1, 1); %only for one muscle
%                     angleNew = getAngle(command) * 2;
%                else
%                   angleNew = (model.desiredAngleMin + (model.desiredAngleMax - model.desiredAngleMin) * rand(1,1)) * 2; % same init range as above
%                 end

                % Object distance = random/determined
                if randObjRange
                    tmpObjRange = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1, 1);
                else
                    tmpObjRange = objRange(objDist);
                end
                %desired vergence [deg]
                angleDes = 2 * atand(model.baseline / (2 * tmpObjRange));

%                 savePathLeft = sprintf('%s/left%.3d.png', imagesSavePath,t);
%                 savePathRight = sprintf('%s/right%.3d.png', imagesSavePath,t);
%                 savePathLeft = sprintf('%s/leftMovie.png', imagesSavePath);
%                 savePathRight = sprintf('%s/rightMovie.png', imagesSavePath);
%                 [status, res] = system(sprintf('./checkEnvironment %s %d %d %s %s', ...
%                                                currentTexture, tmpObjRange, angleNew, savePathLeft, savePathRight));
%                 % abort execution if error occured
%                 if (status)
%                     sprintf('Error in checkEnvironment:\n%s', res)
%                     return;
%                 end
%
%                 [imgRawLeft, imgRawRight] = refreshImages(currentTexture, angleNew, tmpObjRange);

                for iteration = 1 : testInterval
                    % read input images and convert to gray scale
                    model.refreshImages(currentTexture, angleNew / 2, tmpObjRange);

                    for i = 1 : length(model.scModel)
                        model.preprocessImage(imgGrayLeft, i, 1);
                        model.preprocessImage(imgGrayRight, i, 2);
                        currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
                    end

                    [bfFeature, reward, recErrorArray] = model.generateFR(currentView);

                    feature = [bfFeature; command * model.lambdaMuscleFB];

                    %%% Calculate metabolic costs
                    % metCost = getMetCost(command) * 2;
                    
                    relativeCommand = model.rlModel.act(feature);

                    %% now generate anaglyphs before the movement is executed
    %                 imwrite(imfuse(imgGrayLeft, imgGrayRight, 'falsecolor'), [imagesSavePath '/anaglyph.png']);
                    infos = {iteration, tmpObjRange, vseRange(vseIndex), angleDes, angleNew, command', relativeCommand', reward, recErrorArray};
                    generateAnaglyphs(t, markScales, infos)


                    % add the change in muscle Activities to current ones
                    command = command + relativeCommand;     %two muscles
                    command = checkCmd(command);                %restrain motor commands to [0,1]



                    angleNew = getAngle(command) * 2;


                    t = t + 1;
                end
            end
        end
    end



[status, sysout] = system(sprintf('avconv -r %d -i %s/anaglyphs%%03d.png %s/anaglyphs_%s.mp4', frameRate, imagesSavePath, imagesSavePath, timeStamp));
if status
    sprintf('Error in generating movie anaglyphs.mp4:\n%s',sysout)
end
% after the movies are produced, delete all images
system(sprintf('rm %s/anaglyph*.png', imagesSavePath));

%%% Saturation function that keeps motor commands in [0, 1]
%   corresponding to the muscelActivity/metabolicCost tables
function [cmd] = checkCmd(cmd)
    i0 = cmd < 0;
    cmd(i0) = 0;
    i1 = cmd > 1;
    cmd(i1) = 1;
end


% function generateAnaglyphs(leftGray, rightGray, dsRatioL, dsRatioS, foveaL, foveaS, savePath, identifier, markScales, infos)
function generateAnaglyphs(identifier, markScales, infos)

    numberScales = length(model.scModel);

    %defining colors in the image: (from larges to smallest scale)
    scalingColors = {'blue', 'red', 'green'};

%     cLarge = 'blue';
%     cSmall = 'red';
    cText = 'yellow';

    scaleImages = cell(numberScales);

    for scale = 1:numberScales
        %Downsampling Large
        imgLeft = imgGrayLeft(:);
        imgLeft = reshape(imgLeft, size(imgGrayLeft));
        imgRight = imgGrayRight(:);
        imgRight = reshape(imgRight, size(imgGrayRight));

        for ds = 1:log2(model.dsRatio(scale))
            imgLeft = impyramid(imgLeft, 'reduce');
            imgRight = impyramid(imgRight, 'reduce');
        end

        % cut fovea in the center
        [h, w, ~] = size(imgLeft);
        cutOutVertical = [fix(h / 2 + 1 - model.pxFieldOfView(scale) / 2), fix(h / 2 + model.pxFieldOfView(scale) / 2)];
        cutOutHorizontal = [fix(w / 2 + 1 - model.pxFieldOfView(scale) / 2), fix(w / 2 + model.pxFieldOfView(scale) / 2)];


        imgLeft = imgLeft(cutOutVertical(1) : cutOutVertical(2) , cutOutHorizontal(1) : cutOutHorizontal(2));
        imgRight = imgRight(cutOutVertical(1) : cutOutVertical(2) , cutOutHorizontal(1) : cutOutHorizontal(2));

        %create an anaglyph of the two pictures, scale it up and save it
        anaglyph = imfuse(imgLeft, imgRight, 'falsecolor');

        if markScales
            anaglyph = insertShape(anaglyph, 'rectangle', [model.pxFieldOfView(scale) + 1 - model.patchSize, 1, model.patchSize, model.patchSize], 'color', scalingColors(scale));
        end

        scaleImages{scale} = anaglyph;
    end


%     imwrite(anaglyph, [savePath '/anaglyph.png']);
    anaglyph = imfuse(imgGrayLeft, imgGrayRight, 'falsecolor');
    [h, w, ~] = size(anaglyph);

    if markScales
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
%     anaglyph = imfuse(leftGray, rightGray, 'falsecolor');

%     imwrite(anaglyph, sprintf('%s/anaglyph%.3d.png', savePath, identifier));


%     imwrite(imresize(anaglyphL, 20), [savePath '/anaglyphLargeScale.png']);
%     imwrite(anaglyphL, sprintf('%s/anaglyphLargeScale%.3d.png', savePath, identifier));
%     largeScaleView = imfuse(imgLeftL, imgRightL, 'montage');
%     imwrite(imresize(largeScaleView, 20), sprintf('%s/LargeScaleMontage%.3d.png', savePath, identifier));

%     imwrite(imresize(anaglyphS, 8), [savePath '/anaglyphSmallScale.png']);

%     imwrite(anaglyphS, sprintf('%s/anaglyphSmallScale%.3d.png', savePath, identifier));
%     smallScaleView = imfuse(imgLeftL, imgRightL, 'montage');
%     imwrite(imresize(smallScaleView, 8), sprintf('%s/smallScaleMontage%.3d.png', savePath, identifier));

%     xRange = [0, 320]; yRange = [0, 240];
    %% todo: insert reconstruction error for all and for each scale
    %% infos = {iteration, tmpObjRange, vseRange(vseIndex), angleDes, angleNew, command', relativeCommand', reward, recErrorArray};
    angleFix = model.baseline / (2 * tand(angleNew / 2));
    xPos = [10, 220, 220, 10, 10, 260, 10, 10, 240]; %display in headful modus: [10, 200, 10, 10, 260, 10, 10]
    yPos = [10, 10, 20, 230, 220, 230, 190, 200, 220];
    imName = strsplit(currentTexture, '/');
    imName = imName{end};
    insert = {sprintf('Image: \t%s', imName), ...
                sprintf('Object distance: \t%.2f', infos{2}), ...
                sprintf('Fixation depth:   \t%.2f', angleFix), ...
                sprintf('Vergence Error:          \t%.3f', infos{4} - infos{5}), ...
                sprintf('Start Vergence Error: \t%.3f', infos{3}), ...
                sprintf('Iteration: \t%d', infos{1}), ...
                sprintf('Muscle Activation:    \t%.5f,  \t%.5f', infos{6}(1), infos{6}(2)), ...
                sprintf('Relative Command: \t%.5f,  \t%.5f', infos{7}(1), infos{7}(2)), ...
                sprintf('Reward: \t%.3f', infos{8})};

    f = figure('Position', [0, 0, 960, 640]); % create a figure with a specific size
%     subplot(2,3,[1,2,4,5]);
%     title(sprintf('Object distance: %d,\titeration: %d,\tvergence error: %d', infos(1), infos(2), infos(3)));
    subplot('Position', [0.05, 0.12, 0.6, 0.8]);
    imshow(anaglyph);
    text(xPos, yPos, insert, 'color', cText); % these numbers resemble white, but white will appear black after saving ;P
%     text(10, 20, num2str(size(anaglyph)), 'color', 'yellow'); %[1-eps, 1-eps, 1-eps]
%     title('Anaglyph');
%     subplot(2,3,3);
    % positionVector = [left, bottom, width, height], all between 0 and 1
    sizeScaleImg = 0.6/numberScales;
    for sp = 1:numberScales
        if sp == 1
            subplot('Position', [0.7, 0.9 - sizeScaleImg, sizeScaleImg, sizeScaleImg])
        else
            subplot('Position', [0.7, 0.9 - ((sizeScaleImg + 0.05) * sp), sizeScaleImg, sizeScaleImg])
        end
        imshow(scaleImages{sp});
        descr = {sprintf('Scale %d', sp), sprintf('Rec. Error: %.3f', infos{9}(sp))};
        % todo: the fist two values in text have to be adapted, especially
        % for more than 2 scales, also the size of the letters,
        text(0.03, 0.1, descr, 'color', cText, 'Units', 'normalized')
    end
%     subplot('Position', [0.7, 0.6, 0.25, 0.25])
%     imshow(anaglyphS);
%     title('fine scale');
%     text(3, 37, sprintf('Fine Scale'), 'color', cText); % size = 40 40
%     subplot(2,3,6);
%     subplot('Position', [0.7, 0.2, 0.25, 0.25])
%     imshow(anaglyphL);
%     title('coarse scale');
%     text(1.5, 15, sprintf('Coarse Scale'), 'color', cText); % size = 16 16
    saveas(f, sprintf('%s/movies/anaglyphs%03d.png', model.savePath, identifier), 'png');
    close 1;
end

end

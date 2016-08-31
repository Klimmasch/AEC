%%% This method creates a trajectory from the given paramters and plots it
%%% on the plane of object depth.
%%% Note that this script is intended solely for continuous models.
% @param model                  the model to be tested
% @param objDist                the object distance
% @param startVergErr           the vergence error to muscles start with
% @param initMethod             either 'simple' or 'random'
% @param numIters               number of iterations that are executed
% @param stimuli                an array of indizes from the texture files
% @param simulator              either a simulator object or [] for a new one
% @param titleStr               string identifier that is used for the title and the saved image
% @param savePlot               true or false if the resulting plots should be saved
%%%
%%TODO: 
% enable multiple fixation dists in one plot with same init values
% idea: instead of contourf, just plot single lines that correspond to
% spec. obj. dists
function plotTrajectory(model, objDist, startVergErr, initMethod, numIters, stimuli, simulator, titleStr, savePlot)

    rng(1);
    angleDes = atand(model.baseline / (2 * objDist)); % only for one eye this time
    cmdInit = [[0.003; 0.012], [0.003; 0.004], [0.01; 0.004]]; % hand-picked inits for muscles, used in initMethod 'random'

    plotAnaglyphs = true;
%     plotAnaglyphs = false;
    
    %% everything concerning rendering:
%     texture = load(['config/' 'Textures_mcGillTest']);
    texture = load(['config/' 'Textures_vanHaterenTest']);
    texture = texture.texture;
    nTextures = length(texture);
    nStimuli = length(stimuli);
    
    if length(stimuli) > nTextures
        sprintf('Texture file does not contain that many images, but I will continue anyways.')
        stimuli = stimuli(1:length(texture));
    end
    
    if (isempty(simulator)) 
        sprintf('Please execute prepareSimulator.m to create a new simulator!')
        return;
%         simulator = OpenEyeSimV4('create');
%         simulator = OpenEyeSim('create');
%         if (reinit == 0)
%             simulator.initRenderer();
%         else
%             % for debugging purposes
%             simulator.reinitRenderer();
%         end
    end
    
    imgRawLeft = uint8(zeros(240, 320, 3));
    imgRawRight = uint8(zeros(240, 320, 3));
    imgGrayLeft = uint8(zeros(240, 320, 3));
    imgGrayRight = uint8(zeros(240, 320, 3));
    currentView = cell(1, length(model.scModel));
    
    function refreshImages(texture, eyeAngle, objDist, scalingFactor)
%         simulator.add_texture(1, texture);
        simulator.set_params(texture, eyeAngle, objDist, 0, scalingFactor);

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
    
    %% everything you need for muscle commands:
    degrees = load('Degrees.mat');
    
    function angle = getAngle(command)
        cmd = (command * 10) + 1;                                       % scale commands to table entries
        angle = interp2(degrees.results_deg, cmd(1), cmd(2), 'spline'); % interpolate in table by cubic splines
    end

    function [cmd] = checkCmd(cmd)
        i0 = cmd < 0;
        cmd(i0) = 0;
        i1 = cmd > 1;
        cmd(i1) = 1;
    end
    
    resolution = 100001;
    approx = spline(1 : 11, degrees.results_deg(:, 1));

    xValPos = ppval(approx, 1 : 0.0001 : 11)';
    yValPos = linspace(0, 1, resolution)';
    
    mfunction = [xValPos, yValPos];
    mfunction(:, 1) = mfunction(:, 1) * 2;  % angle for two eyes
    dmf = abs(diff(mfunction(1 : 2, 1)));   % delta in angle
    dmf2 = diff(mfunction(1 : 2, 2));       % delta in mf
    indZero = find(mfunction(:, 2) == 0);   % MF == 0_index
    
    approx = spline(1 : 11, degrees.results_deg(1, :));
    xValPos = ppval(approx, 1 : 0.0001 : 11)';
    yValPos = linspace(0, 1, resolution)';
    
    mfunction2 = [xValPos, yValPos];
    mfunction2(:, 1) = mfunction2(:, 1) * 2;    % angle for two eyes
    dmf3 = abs(diff(mfunction2(1 : 2, 1)));     % delta in angle
    dmf4 = diff(mfunction2(1 : 2, 2));          % delta in mf
    indZero = find(mfunction2(:, 2) == 0);      % MF == 0_index
    
    angleMin = min(mfunction2(mfunction2(:, 1) > 0));
    angleMax = mfunction(end, 1);
    
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
    
    %% main loop:
    trajectory = zeros(length(stimuli), numIters, 2);
    figure;
%     title(sprintf('Trajectory of Oject Fixation at %1.1fm\n%s', objDist, titleStr));
    figIter = 1;
    for stim = 1:nStimuli
%         currentTexture = texture{stimuli(stim)};
        currentTexture = stimuli(stim);
        
        if strcmp(initMethod, 'simple')
            [command, angleNew] = getMF2(objDist, startVergErr);
            angleNew = angleNew / 2;
        elseif strcmp(initMethod, 'random')
%             initFixation = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1,1)
%             command = [0; 0];
%             command(1) = model.muscleInitMin(1) + (model.muscleInitMax(1) - model.muscleInitMin(1)) * rand(1, 1);
%             command(2) = model.muscleInitMin(2) + (model.muscleInitMax(2) - model.muscleInitMin(2)) * rand(1, 1); 
            command = cmdInit(:, stim);
            angleNew = getAngle(command);
        end
            
        for iter = 1:numIters
            trajectory(stim, iter, :) = command;
            refreshImages(currentTexture, angleNew, objDist, 3)
            
            % show anaglyphs for quit performance check
            if (plotAnaglyphs) && ((iter == 1) || (iter == numIters))
                subplot(nStimuli, 2, figIter);
                imshow(stereoAnaglyph(imgGrayLeft, imgGrayRight))
                title(sprintf('fix. depth = %1.1fm (%.3f°)\nverg. error = %.3f', (model.baseline / (2 * tand(angleNew))), angleNew, angleDes - angleNew)) 
                figIter = figIter + 1;
%             elseif (plotAnaglyphs) && (iter == numIters)
%                 subplot(nStimuli, 2, figIter);
%                 imshow(stereoAnaglyph(imgGrayLeft, imgGrayRight))
%                 title(sprintf('angle = %1.1fm (%.3f°),\nverg. error = %.3f', (model.baseline / (2 * tand(angleNew))), angleNew, angleDes - angleNew))   
%                 figIter = figIter + 1;
            end
            
            for i = 1 : length(model.scModel)
                model.preprocessImage(imgGrayLeft, i, 1);
                model.preprocessImage(imgGrayRight, i, 2);
                currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
            end

            [bfFeature, ~, recErrorArray] = model.generateFR(currentView);  % encode image patches
            feature = [bfFeature; command * model.lambdaMuscleFB];          % append muscle activities to feature vector
            relativeCommand = model.rlModel.act(feature);                   % generate change in muscle activity
            command = checkCmd(command + relativeCommand);                  % calculate new muscle activities
            angleNew = getAngle(command);                               % transform into angle
            
            
        end
    end
    
    %% plotting results
%     figure('size', [500, 500]);
    h = figure();
    hold on;
    title(sprintf('Oject Fixation at %1.1fm (%.3f°)\n%s', objDist, angleDes, titleStr));
    resolutionFactor = 6;
    tableSize = size(degrees.results_deg)-1;
    
    % set maximum x and y values of the tabular that should be plotted
    xlim = 0.1;             % unfortunately, these are the smallest values since the tabular only has 11 entries
    ylim = 0.1;
    plotRange = [xlim, ylim];
    
    % only this part of the tabular in increased in resolution
    rangeIndizes = ceil(plotRange .* tableSize)+1;
    degreesX = interp2(degrees.results_deg(1:rangeIndizes(1), 1:rangeIndizes(2)), resolutionFactor);
    degSize = size(degreesX);                                       % corresp. to ((rangeIndizes-1)*2^resolutionFactor)+1
    
    % now we cut out an even smaller part of this tabular.
    factor = 1      % 1/factor * axisLim is the part that will be displayed
    maxIndexX = ceil(degSize(1)/factor);
    maxIndexY = ceil(degSize(2)/factor);
    
    scaleSize = ((tableSize(1)-1)*2^resolutionFactor)+1;            % corresp. to the size of the whole tabular with increased resolution
%     fun = surf(degreesX);                                           % TODO: change axis labels 
%     fun = contourf(linspace(0, xlim / factor, maxIndexX), linspace(0, ylim / factor, maxIndexY), degreesX(1:maxIndexX, 1:maxIndexY), 'LineWidth', 0.001);  % temporary solution to enable plotRange = [0.025, 0.025]
    fun = contourf(degreesX(1:maxIndexX, 1:maxIndexY), 20, 'LineWidth', 0.001);  % temporary solution to enable plotRange = [0.025, 0.025]
    colBar = colorbar();
%     colorbar('Ticks',[-5,-2,1,4,7],...
%          'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'}) %% todo:
%          replace vergence degrees with fixation distance
    colBar.Label.String = 'vergence degree';
    
    %     alpha(fun, 0.5)                                                 % add transparency to degree plot
    
    scalingFactor = degSize ./ plotRange;
    for stim = 1:length(stimuli)
        plot(trajectory(stim,1,1) * scalingFactor + 1, trajectory(stim,1,2) * scalingFactor + 1, 'g.', 'MarkerSize', 50); % first value gets a bigger dot
        plot(trajectory(stim,:,1)' * scalingFactor + 1, trajectory(stim,:,2)' * scalingFactor + 1, '.-', 'LineWidth', 2, 'MarkerSize', 20);
        plot(trajectory(stim,end,1) * scalingFactor + 1, trajectory(stim,end,2) * scalingFactor + 1, 'r.', 'MarkerSize', 50); % , 'MarkerSize', 100
    end
    
    fun = gca();    
    fun.XAxis.TickLabels = linspace(0, xlim/factor, maxIndexX); %TODO: get auf jeden Fall noch schoener --> weniger ticks!
    fun.YAxis.TickLabels = linspace(0, ylim/factor, maxIndexY);

    fun.XTickLabelRotation = -20;
    
    xlabel('lateral rectus activation [%]');
    ylabel('medial rectus activation [%]');
%     fun.XTick = 1 : length(linspace(0, xlim/factor, maxIndexX)); %TODO: get auf jeden Fall noch schoener --> weniger ticks!
%     fun.YTick = 1 : length(linspace(0, ylim/factor, maxIndexY));
    
%     fun.XTickLabel = sprintf('%.4f ', linspace(0, xlim/factor, maxIndexX)); %TODO: get auf jeden Fall noch schoener --> weniger ticks!
%     fun.YTickLabel = sprintf('%.4f ', linspace(0, ylim/factor, maxIndexY));
    
    
%     xt = get(fun,'xtick');
% %     xt_label = sprintf('%.2e|',xt);
%     xt_label = sprintf('%.4f',xt);
%     set(fun,'xticklabel',xt_label);
    
%     yt = get(fun,'ytick');
% %     yt_label = sprintf('%.2e|',xt);
%     yt_label = sprintf('%.4f|',yt);
%     set(fun,'yticklabel',yt_label);
    if savePlot
        timestamp = datestr(now, 'dd-mm-yyyy_HH:MM:SS_');
        savePath = strcat(model.savePath, '/', timestamp, titelStr);
        saveas(h, savePath, 'png'); % could be saved as 'fig' as well ...
    end
    
end
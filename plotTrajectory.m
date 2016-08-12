%%% This method creates a trajectory from the given paramters and plots it
%%% on the plane of object depth.
%%% Note that this script is intended solely for continuous models.

function plotTrajectory(model, objDist, startVergErr, initMethod, numIters, stimuli, simulator, reinit, savePlot)

    %% everything concerning rendering:
%     texture = load(['config/' 'Textures_mcGillTest']);
    texture = load(['config/' 'Textures_vanHaterenTest']);
    texture = texture.texture;
    if length(stimuli) > length(texture)
        sprintf('Texture file does not contain that many images, but I will continue anyways.')
        stimuli = stimuli(1:length(texture));
    end
    
    if (isempty(simulator))
%         simulator = OpenEyeSimV4('create');
        simulator = OpenEyeSim('create');
        if (reinit == 0)
            simulator.initRenderer();
        else
            % for debugging purposes
            simulator.reinitRenderer();
        end
    end
    
    imgRawLeft = uint8(zeros(240, 320, 3));
    imgRawRight = uint8(zeros(240, 320, 3));
    imgGrayLeft = uint8(zeros(240, 320, 3));
    imgGrayRight = uint8(zeros(240, 320, 3));
    currentView = cell(1, length(model.scModel));
    
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
    for stim = 1:length(stimuli)
        currentTexture = texture{stimuli(stim)};
        
        if initMethod == 'simple'
            [command, angleNew] = getMF2(objDist, startVergErr);
        elseif initMethod == 'random'
%             initFixation = model.objDistMin + (model.objDistMax - model.objDistMin) * rand(1,1)
            command = [0; 0];
            command(1) = model.muscleInitMin(1) + (model.muscleInitMax(1) - model.muscleInitMin(1)) * rand(1, 1);
            command(2) = model.muscleInitMin(2) + (model.muscleInitMax(2) - model.muscleInitMin(2)) * rand(1, 1); 
            angleNew = getAngle(command) * 2;
        end
            
        for iter = 1:numIters
            trajectory(stim, iter, :) = command;
            refreshImages(currentTexture, angleNew / 2, objDist, 3)

            for i = 1 : length(model.scModel)
                model.preprocessImage(imgGrayLeft, i, 1);
                model.preprocessImage(imgGrayRight, i, 2);
                currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
            end

            [bfFeature, ~, recErrorArray] = model.generateFR(currentView);  % encode image patches
            feature = [bfFeature; command * model.lambdaMuscleFB];          % append muscle activities to feature vector
            relativeCommand = model.rlModel.act(feature);                   % generate change in muscle activity
            command = checkCmd(command + relativeCommand);                  % calculate new muscle activities
            angleNew = getAngle(command) * 2;                               % transform into angle
        end
    end
    
    %% plotting results
    figure;
    hold on;
    resolutionFactor = 6;
    tableSize = size(degrees.results_deg)-1;
    
    plotRange = [0.1, 0.1];                                         % maximum x and y values of the tabular that should be plotted
    rangeIndizes = ceil(plotRange .* tableSize)+1;
    
    degreesX = interp2(degrees.results_deg(1:rangeIndizes(1), 1:rangeIndizes(2)), resolutionFactor);
    degSize = size(degreesX);                                       % corresp. to ((rangeIndizes-1)*2^resolutionFactor)+1
    scaleSize = ((tableSize(1)-1)*2^resolutionFactor)+1;            % corresp. to the size of the whole tabular with increased resolution
%     fun = surf(degreesX);                                           % TODO: change axis labels
    fun = contourf(degreesX(1:ceil(degSize(1)/8), 1:ceil(degSize/8)), 'LineWidth', 0.001);  % temporary solution to enable plotRange = [0.025, 0.025]
%     alpha(fun, 0.5)                                                 % add transparency to degree plot
    
    scalingFactor = degSize ./ plotRange;
    for stim = 1:length(stimuli)
        plot(trajectory(stim,1,1) * scalingFactor + 1, trajectory(stim,1,2) * scalingFactor + 1, 'g.', 'MarkerSize', 50); % first value gets a bigger dot
        plot(trajectory(stim,:,1)' * scalingFactor + 1, trajectory(stim,:,2)' * scalingFactor + 1, '.-', 'LineWidth', 2, 'MarkerSize', 20);
        plot(trajectory(stim,end,1) * scalingFactor + 1, trajectory(stim,end,2) * scalingFactor + 1, 'r.', 'MarkerSize', 50); % , 'MarkerSize', 100
    end
    
    
%     plot(trajectory);
end
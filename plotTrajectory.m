%%% This method creates a trajectory from the given paramters and plots it
%%% on the plane of object depth.
%%% Note that this script is intended solely for continuous models.
% @param model                  the model to be tested
% @param objDist                the object distance
% @param startVergErr           the vergence error to muscles start with
% @param initMethod             either 'simple' or 'random'
% @param numIters               number of iterations that are executed
% @param stimuliIndices         an array of indizes from the texture files
% @param simulator              either a simulator object or [] for a new one
% @param titleStr               string identifier that is used for the title and the saved image
% @param savePlot               true or false if the resulting plots should be saved
%%%
%%TODO:
% enable multiple fixation dists in one plot with same init values
function plotTrajectory(model, objDist, startVergErr, initMethod, numIters, stimuliIndices, simulator, titleStr, savePlot)

    % simulator check
    if (isempty(simulator))
        sprintf('An initialized simulator is necessary to continue.\nPlease execute simulator = prepareSimulator();')
        return;
    end

    %%% Saturation function that keeps motor commands in [0, 1]
    %   corresponding to the muscelActivity/metabolicCost tables
    function [cmd] = checkCmd(cmd)
        i0 = cmd < 0;
        cmd(i0) = 0;
        i1 = cmd > 1;
        cmd(i1) = 1;
    end

    % preperation
    rng(1);

    if strcmp(initMethod, 'advanced')
        initMethod = uint8(0);

    elseif strcmp(initMethod, 'fixed')
        initMethod = uint8(1);
        cmdInit = [[0.003; 0.012], [0.003; 0.004], [0.01; 0.004]]; % hand-picked inits for muscles, used in initMethod 'random'

    elseif strcmp(initMethod, 'simple')
        initMethod = uint8(2);

    else
        sprintf('Muscle initialization method %s not supported.', initMethod)
        return;
    end

    plotAnaglyphs = true;
    nStimuli = length(stimuliIndices);
    trajectory = zeros(nStimuli, numIters + 1, 2);

    %% main loop:
    angleDes = 2 * atand(model.baseline / (2 * objDist));
    figure;
    figIter = 1;
    for stimIter = 1 : nStimuli
        currentTexture = stimuliIndices(stimIter);

        % muscle init
        if (initMethod == 0)
            try
                [command, angleNew] = model.getMFedood(objDist, startVergErr);
            catch
                % catch non-existing variables error, occuring in non-up-to-date models
                try
                    clone = model.copy();
                    delete(model);
                    clear model;
                    model = clone;
                    [command, angleNew] = model.getMFedood(objDist, startVergErr);
                    delete(clone);
                    clear clone;
                catch
                    % catch when new model property isn't present in Model class yet
                    sprintf('Error: One or more new model properties (variables) are not present in Model.m class yet!')
                    return;
                end
            end
        elseif (initMethod == 1)
            command = cmdInit(:, stimIter);
            angleNew = model.getAngle(command);
        elseif (initMethod == 2)
            [command, angleNew] = model.getMF2(objDist, startVergErr);
        end
        trajectory(stimIter, 1, :) = command;

        for iter = 1 : numIters
            model.refreshImagesNew(simulator, currentTexture, angleNew / 2, objDist, 3);

            % show anaglyphs for quit performance check
            if (plotAnaglyphs && ((iter == 1) || (iter == numIters)))
                subplot(nStimuli, 2, figIter);
                imshow(stereoAnaglyph(model.imgGrayLeft, model.imgGrayRight))
                if (iter == 1)
                    title(sprintf('fix. depth = %1.1fm (%.3f°)\nverg. error = %.3f', (model.baseline / 2) / tand(angleNew / 2), angleNew, angleDes - angleNew));
                end
                figIter = figIter + 1;
            end

            for i = 1 : length(model.scModel)
                model.preprocessImage(i, 1);
                model.preprocessImage(i, 2);
                currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
            end

            [bfFeature, ~, ~] = model.generateFR(currentView);              % encode image patches
            feature = [bfFeature; command * model.lambdaMuscleFB];          % append muscle activities to feature vector
            relativeCommand = model.rlModel.act(feature);                   % generate change in muscle activity
            command = checkCmd(command + relativeCommand);                  % calculate new muscle activities
            angleNew = model.getAngle(command) * 2;                         % transform into angle

            trajectory(stimIter, iter + 1, :) = command;

            if (iter == numIters)
                title(sprintf('fix. depth = %1.1fm (%.3f°)\nverg. error = %.3f', (model.baseline / 2) / tand(angleNew / 2), angleNew, angleDes - angleNew));
            end
        end
    end

    %% Plotting results
    h = figure();
    hold on;
    title(sprintf('Oject Fixation Trajectories at %1.1fm (%.3f°)\n%s', objDist, angleDes, titleStr));

    % pcHandle = pcolor(model.degreesIncRes); % use vergence degree as color dimension (background)
    pcHandle = pcolor(model.metCostsIncRes);  % use metabolic costs as color dimension (background)
    % shading interp;
    set(pcHandle, 'EdgeColor', 'none');

    cb = colorbar();
    % cb.Label.String = 'vergence degree'; % use vergence degree as color dimension (background)
    cb.Label.String = 'metabolic costs';   % use metabolic costs as color dimension (background)

    ax = gca;
    ax.XTick = linspace(1, size(model.degreesIncRes, 2), 8);
    ax.YTick = linspace(1, size(model.degreesIncRes, 1), 8);

    ax.XTickLabel = strsplit(num2str(linspace(1, size(model.degreesIncRes, 2), 8) * model.scaleFacLR, '%4.2f '));
    ax.YTickLabel = strsplit(num2str(linspace(1, size(model.degreesIncRes, 1), 8) * model.scaleFacMR, '%4.2f '));

    ax.XTickLabelRotation = 45;
    ax.YTickLabelRotation = 45;

    axis([1, size(model.degreesIncRes, 2), 1, size(model.degreesIncRes, 1)]);

    % draw a line of points into the plane that represent the desired vergence
    [lateralDes, medialDes] = model.getAnglePoints(objDist, 0);
    plot(lateralDes ./ model.scaleFacLR, medialDes ./ model.scaleFacMR, 'g', 'LineWidth', 1.8);

    % add corresponding distance value to desired vergence graph
    text(lateralDes(end - ceil(length(lateralDes) / 10)) / model.scaleFacLR, ...
         medialDes(end - ceil(length(medialDes) / 10)) / model.scaleFacMR, ...
         sprintf('%3.1fm', (model.baseline / 2) / tand(angleDes / 2)));

    % draw trajectories
    for stim = 1 : length(stimuliIndices)
        plot(trajectory(stim, 1, 1) ./ model.scaleFacLR, trajectory(stim, 1, 2)./ model.scaleFacMR, 'r.', 'MarkerSize', 40);
        plot(trajectory(stim, :, 1)' ./ model.scaleFacLR, trajectory(stim, :, 2)'./ model.scaleFacMR, '.-', 'LineWidth', 2, 'MarkerSize', 20);
        plot(trajectory(stim, end, 1) ./ model.scaleFacLR, trajectory(stim, end, 2)./ model.scaleFacMR, 'g.', 'MarkerSize', 40);
    end

    xlabel('lateral rectus activation [%]');
    ylabel('medial rectus activation [%]');

    if savePlot
        timestamp = datestr(now, 'dd-mm-yyyy_HH:MM:SS_');
        savePath = strcat(model.savePath, '/', timestamp, titelStr);
        saveas(h, savePath, 'png');
    end
end

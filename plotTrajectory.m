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

    %%% Saturation function that keeps motor commands in [0, 1]
    %   corresponding to the muscelActivity/metabolicCost tables
    function [cmd] = checkCmd(cmd)
        i0 = cmd < 0;
        cmd(i0) = 0;
        i1 = cmd > 1;
        cmd(i1) = 1;
    end

    rng(1);
    angleDes = 2 * atand(model.baseline / (2 * objDist)); % only for one eye this time
    cmdInit = [[0.003; 0.012], [0.003; 0.004], [0.01; 0.004]]; % hand-picked inits for muscles, used in initMethod 'random'

    plotAnaglyphs = true;
%     plotAnaglyphs = false;

    % simulator = prepareSimulator();
    if (isempty(simulator))
        sprintf('Please execute simulator = prepareSimulator();')
        return;
    end

    nStimuli = length(stimuli);

    %% main loop:
    trajectory = zeros(length(stimuli), numIters + 1, 2);
    figure;
    figIter = 1;
    for stim = 1 : nStimuli
        currentTexture = stimuli(stim);

        if strcmp(initMethod, 'simple')
            [command, angleNew] = model.getMF2(objDist, startVergErr);
        elseif strcmp(initMethod, 'random')
            command = cmdInit(:, stim);
            angleNew = model.getAngle(command);
        elseif strcmp(initMethod, 'advanced')
            try
                [command, angleNew] = model.getMFedood(objDist, startVergErr, false);
            catch
                % catch non-existing variables error, occuring in non-up-to-date models
                try
                    clone = model.copy();
                    delete(model);
                    clear model;
                    model = clone;
                    [command, angleNew] = model.getMFedood(objDist, startVergErr, false);
                    delete(clone);
                    clear clone;
                catch
                    % catch when new model property isn't present in Model class yet
                    sprintf('Error: One or more new model properties (variables) are not present in Model.m class yet!')
                    return;
                end
            end
        end
        trajectory(stim, 1, :) = command;

        model.refreshImagesNew(simulator, currentTexture, angleNew / 2, objDist, 3);

        subplot(nStimuli, 2, figIter);
        imshow(stereoAnaglyph(model.imgGrayLeft, model.imgGrayRight))
        title(sprintf('fix. depth = %1.1fm (%.3f°)\nverg. error = %.3f', (model.baseline / 2) / (tand(angleNew / 2)), angleNew, angleDes - angleNew));
        figIter = figIter + 1;

        for iter = 1 : numIters
            model.refreshImagesNew(simulator, currentTexture, angleNew / 2, objDist, 3);

            % show anaglyphs for quit performance check
            % if (plotAnaglyphs) && ((iter == 1) || (iter == numIters))
            %     subplot(nStimuli, 2, figIter);
            %     imshow(stereoAnaglyph(model.imgGrayLeft, model.imgGrayRight))
            %     title(sprintf('fix. depth = %1.1fm (%.3f°)\nverg. error = %.3f', (model.baseline / 2) / (tand(angleNew / 2)), angleNew, angleDes - angleNew));
            %     figIter = figIter + 1;
            % end

            for i = 1 : length(model.scModel)
                model.preprocessImage(i, 1);
                model.preprocessImage(i, 2);
                currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
            end

            [bfFeature, ~, ~] = model.generateFR(currentView);              % encode image patches
            feature = [bfFeature; command * model.lambdaMuscleFB];          % append muscle activities to feature vector
            relativeCommand = model.rlModel.act(feature);                   % generate change in muscle activity
            command = checkCmd(command + relativeCommand);                  % calculate new muscle activities
            angleNew = model.getAngle(command) * 2;                                   % transform into angle

            trajectory(stim, iter + 1, :) = command;
        end

        subplot(nStimuli, 2, figIter);
        imshow(stereoAnaglyph(model.imgGrayLeft, model.imgGrayRight))
        title(sprintf('fix. depth = %1.1fm (%.3f°)\nverg. error = %.3f', (model.baseline / 2) / tand(angleNew / 2), angleNew, angleDes - angleNew));
        figIter = figIter + 1;
    end

    %% plotting results
    resolutionFactor = 6;
    tableSize = size(model.degrees.results_deg)-1;

    % set maximum x and y values of the tabular that should be plotted
    xlim = 0.1;             % unfortunately, these are the smallest values since the tabular only has 11 entries
    ylim = 0.1;
    plotRange = [xlim, ylim];

    % only this part of the tabular in increased in resolution
    rangeIndizes = ceil(plotRange .* tableSize)+1;
    degreesX = interp2(model.degrees.results_deg(1:rangeIndizes(1), 1:rangeIndizes(2)), resolutionFactor);
    degSize = size(degreesX);                                       % corresp. to ((rangeIndizes-1)*2^resolutionFactor)+1

    % now we cut out an even smaller part of this tabular.
    factor = 1;      % 1/factor * axisLim is the part that will be displayed
    maxIndexX = ceil(degSize(1)/factor);
    maxIndexY = ceil(degSize(2)/factor);

    scaleSize = ((tableSize(1)-1)*2^resolutionFactor)+1;            % corresp. to the size of the whole tabular with increased resolution
%     fun = surf(degreesX);                                           % TODO: change axis labels
%     fun = contourf(linspace(0, xlim / factor, maxIndexX), linspace(0, ylim / factor, maxIndexY), degreesX(1:maxIndexX, 1:maxIndexY), 'LineWidth', 0.001);  % temporary solution to enable plotRange = [0.025, 0.025]
%     fun = contourf(degreesX(1:maxIndexX, 1:maxIndexY), 20, 'LineWidth', 0.001);  % temporary solution to enable plotRange = [0.025, 0.025]
%     colBar = colorbar();
% %     colorbar('Ticks',[-5,-2,1,4,7],...
% %          'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'}) %% todo:
% %          replace vergence degrees with fixation distance
%     colBar.Label.String = 'vergence degree';

    h = figure();
    hold on;
    title(sprintf('Oject Fixation at %1.1fm (%.3f°)\n%s', objDist, angleDes, titleStr));

    pcHandle = pcolor(1 : maxIndexX, 1 : maxIndexY, degreesX);
    % axis([xb(1), xb(end), yb(1), yb(end)]);
    % shading interp;
    set(pcHandle, 'EdgeColor', 'none');

    % colormap(createCM());
    cb = colorbar();
    cb.Label.String = 'vergence degree';

    scalingFactor = degSize ./ plotRange;

    % draw a line of points into the plane that represent the desired vergence
    [lateralDes, medialDes] = model.getAnglePoints(objDist, 0);
    plot(lateralDes * scalingFactor + 1, medialDes * scalingFactor + 1, 'g', 'LineWidth', 1.8); %[0, 0.7255, 0.1765]

    % draw trajectories
    for stim = 1:length(stimuli)
        plot(trajectory(stim,1,1) * scalingFactor + 1, trajectory(stim,1,2) * scalingFactor + 1, 'r.', 'MarkerSize', 40); % first value gets a bigger dot
        plot(trajectory(stim,:,1)' * scalingFactor + 1, trajectory(stim,:,2)' * scalingFactor + 1, '.-', 'LineWidth', 2, 'MarkerSize', 20);
        plot(trajectory(stim,end,1) * scalingFactor + 1, trajectory(stim,end,2) * scalingFactor + 1, 'g.', 'MarkerSize', 40); % , 'MarkerSize', 100
    end

    fun = gca();
    fun.XAxis.TickLabels = linspace(0, xlim/factor, maxIndexX); %TODO: get auf jeden Fall noch schoener --> weniger ticks!
    fun.YAxis.TickLabels = linspace(0, ylim/factor, maxIndexY);

    fun.XTickLabelRotation = -20;

    xlabel('lateral rectus activation [%]');
    ylabel('medial rectus activation [%]');

    if savePlot
        timestamp = datestr(now, 'dd-mm-yyyy_HH:MM:SS_');
        savePath = strcat(model.savePath, '/', timestamp, titelStr);
        saveas(h, savePath, 'png'); % could be saved as 'fig' as well ...
    end

end

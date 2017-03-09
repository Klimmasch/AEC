% Generates results plots for ICDL conference publication 2017
function generateICDLPlots(modelWoMetCostsFullPath, modelWMetCostsFullPath) %, testAtiter)

    try
        modelHandle = [load(strcat(modelWoMetCostsFullPath, '/model.mat')), ...
                       load(strcat(modelWMetCostsFullPath, '/model.mat'))];
    catch
        error('Model(s) could not be loaded.');
    end

    savePath = strcat('/home/aecgroup/aecdata/ICDLPlots', datestr(now, '/dd-mm-yy_HH:MM:SS'));
    mkdir(savePath);

    % in respect of old testHist(end) == 0 bug
    adjust = 0;
    for i = 1 : length(modelHandle)
        if (modelHandle(i).model.testHist(end, 1) == 0)
            adjust = 1;
        end
    end

    %%% Figure A
    % RMSE vergence error [deg] & delta MC opt [%] @ testing vs. traintime
    figA = figure();
    hold on;
    grid on;

    lineHandles = [0, 0];       % [w/ metCosts, w/o metCosts]
    lineStyles = [':', '-'];    % [w/ metCosts, w/o metCosts]
    lineWidths = [1.3, 1.3];    % [vergErr, metCosts]

    markerStyles = ['x', 'x'];  % [w/ metCosts, w/o metCosts]
    markerSizes = [5, 5];       % [vergErr, metCosts]
    colors = ['b', 'r'];        % [vergErr, metCosts]

    for i = 1 : length(modelHandle)
        % sort fields in ascending order
        % [modelHandle(i).model.testAt, sortIndex] = sort(modelHandle(i).model.testAt);
        %  modelHandle(i).model.testHist = modelHandle(i).model.testHist(sortIndex, :);

        %%% RMSE vergence error [deg] -> 1st y-axis
        % color trick -> black entries in legend
        if (i == 1)
            axTmp = plot(modelHandle(i).model.testAt(1 : end - adjust), modelHandle(i).model.testHist(1 : end - adjust, 1), ...
                         'LineStyle', lineStyles(i), 'Marker', markerStyles(i), 'MarkerSize', markerSizes(1), 'Color', 'k', 'LineWidth', lineWidths(1));

            plot(modelHandle(i).model.testAt(1 : end - adjust), modelHandle(i).model.testHist(1 : end - adjust, 1), ...
                 'LineStyle', lineStyles(i), 'Marker', markerStyles(i), 'MarkerSize', markerSizes(1), 'Color', colors(1), 'LineWidth', lineWidths(1));
        else
            axTmp = plot(ax1, modelHandle(i).model.testAt(1 : end - adjust), modelHandle(i).model.testHist(1 : end - adjust, 1), ...
                         'LineStyle', lineStyles(i), 'Marker', markerStyles(i), 'MarkerSize', markerSizes(1), 'Color', 'k', 'LineWidth', lineWidths(1));

            plot(ax1, modelHandle(i).model.testAt(1 : end - adjust), modelHandle(i).model.testHist(1 : end - adjust, 1), ...
                 'LineStyle', lineStyles(i), 'Marker', markerStyles(i), 'MarkerSize', markerSizes(1), 'Color', colors(1), 'LineWidth', lineWidths(1));
        end

        lineHandles(i) = axTmp;

        if (i == 1)
            ax1 = gca; % current axes
            ax1.YColor = colors(1);
            ax1.XAxis.Label.String = 'Traintime';
            ax1.XAxis.Label.FontSize = 12;
            ax1.Title.String = 'Test Performance & Metabolic Costs vs. Traintime';
            ax1.YAxis.Label.String = 'RMSE(verg_{err}) [deg]';
            ax1.YAxis.Label.FontSize = 12;

            ax2 = axes('Position', ax1.Position, ...
                       'YAxisLocation', 'right', ...
                       'Color', 'none');
        end

        %%% mean, std deltaMetCost -> 2nd y-axis
        % TODO scale to %
        [hl, hp] = boundedline(modelHandle(i).model.testAt(1 : end - adjust), modelHandle(i).model.testHist(1 : end - adjust, 5), modelHandle(i).model.testHist(1 : end - adjust, 6), 'alpha');

        hl.Parent = ax2;
        hp.Parent = ax2;

        hl.Marker = markerStyles(i);
        hl.MarkerSize = markerSizes(2);

        hl.Color = colors(2);
        hp.FaceColor = colors(2);
        hl.LineWidth = lineWidths(2);
        hl.LineStyle = lineStyles(i);

        if (i == 1)
            % |\DeltaMC_{opt}| = |MC_{actual} - MC_{optimal}| / |MC_{start} - MC_{optimal}|
            ax2.YAxis.Label.String = '|\DeltaMC_{opt}| [%]';
            ax2.YAxis.Label.FontSize = 12;
            ax2.YColor = colors(2);
        end
    end

    % ax1.YAxis.Limits = [0, inf];
    ax2.YAxis.Limits = [0, inf];

    l = legend(lineHandles);
    l.Box = 'off';
    l.String{1} = 'w/o met. costs';
    l.String{2} = 'w/ met. costs';

    plotpath = sprintf('%s/testPerformanceDMCVsTraintime', savePath);
    saveas(figA, plotpath, 'png');
    close(figA);

    %%% Figure B
    % vergence error [deg] & MC opt approach [%] @ testing vs. iteration step

    %%% Figure C
    % muscleplain trajectories
end

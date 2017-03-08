% Generates results plots for ICDL conference publication 2017
function generateICDLPlots(modelWoMetCostsFullPath, modelWMetCostsFullPath) %, testAtiter)

    try
        modelHandle = [load(strcat(modelWoMetCostsFullPath, '/model.mat')), ...
                       load(strcat(modelWMetCostsFullPath, '/model.mat'))];
    catch
        warning('Model(s) could not be loaded.');
    end

    savePath = strcat('/home/aecgroup/aecdata/ICDLPlots', datestr(now, '/dd-mm-yy_HH:MM:SS'));
    mkdir(savePath);

    %%% Figure A
    % RMSE vergence error [deg] & delta MC opt [%] @ testing vs. traintime
    figA = figure();
    hold on;
    grid on;
    lineStyles = [':', '-'];
    markerStyles = ['x', 'x'];
    markerSizes = [5, 5];
    colors = ['b', 'r'];
    lineHandles = [0, 0];

    for i = 1 : length(modelHandle)
        % sort fields in ascending order
%         [modelHandle(i).model.testAt, sortIndex] = sort(modelHandle(i).model.testAt);
%          modelHandle(i).model.testHist = modelHandle(i).model.testHist(sortIndex, :);

        %%% RMSE vergence error [deg] -> 1st y-axis
        % color trick -> black entries in legend
        axTmp = plot(modelHandle(i).model.testAt, modelHandle(i).model.testHist(:, 1), ...
                     'LineStyle', lineStyles(i), 'Marker', markerStyles(i), 'MarkerSize', markerSizes(1), 'Color', 'k', 'LineWidth', 1.3);
        lineHandles(i) = axTmp;

        plot(modelHandle(i).model.testAt, modelHandle(i).model.testHist(:, 1), ...
             'LineStyle', lineStyles(i), 'Marker', markerStyles(i), 'MarkerSize', markerSizes(1), 'Color', colors(1), 'LineWidth', 1.3);

        if (i == 1)
            ax1 = gca; % current axes
            ax1.YColor = colors(1);
            ax1Pos = ax1.Position; % position of first axes

            ax1.XAxis.Label.String = 'Traintime';
            ax1.XAxis.Label.FontSize = 12;
            ax1.Title.String = 'Test Performance & Metabolic Costs vs. Traintime';
            ax1.YAxis.Label.String = 'RMSE(verg_{err}) [deg]';
            ax1.YAxis.Label.FontSize = 12;

            ax2 = axes('Position', ax1Pos, ...
                       'YAxisLocation', 'right', ...
                       'Color', 'none');
        end

        %%% mean, std deltaMetCost -> 2nd y-axis
        % TODO scale to %
        [hl, hp] = boundedline(modelHandle(i).model.testAt, modelHandle(i).model.testHist(:, 5), modelHandle(i).model.testHist(:, 6), 'alpha');

        hl.Parent = ax2;
        hp.Parent = ax2;

        hl.Marker = 'x';
        hl.MarkerSize = markerSizes(2);

        hl.Color = colors(2);
        hp.FaceColor = colors(2);
        hl.LineWidth = 1.6;
        hl.LineStyle = lineStyles(i);

        if (i == 1)
            ax2.YColor = colors(2);
            % |\DeltaMC_{opt}| = |MC_{actual} - MC_{optimal}| / |MC_{start} - MC_{optimal}|
            ax2.YAxis.Label.String = '|\DeltaMC_{opt}| [%]';
            ax2.YAxis.Label.FontSize = 12;
        end
    end

    % ax1.YAxis.Limits = [0, inf];
    ax2.YAxis.Limits = [0, inf];

    l = legend(lineHandles);
    l.Box = 'off';
    l.String{1} = 'w/o met. costs';
    l.String{2} = 'w/ met. costs';

    % ax1.YAxis.Label.String = 'RMSE(verg_{err}) [deg]';
    % ax1.YAxis.Label.FontSize = 12;

    % % |\DeltaMC_{opt}| = |MC_{actual} - MC_{optimal}| / |MC_{start} - MC_{optimal}|
    % ax2.YAxis.Label.String = '|\DeltaMC_{opt}|';
    % ax2.YAxis.Label.FontSize = 12;

    % ax1.XAxis.Label.String = 'Traintime';
    % ax1.XAxis.Label.FontSize = 12;
    % ax1.Title.String = 'Test Performance & Metabolic Costs vs. Traintime';

    plotpath = sprintf('%s/testPerformanceDMCVsTraintime', savePath);
    saveas(figA, plotpath, 'png');
    close(figA);





    % mean, std vergErr
    % subplot(2, 2, 2);
    % hold on;
    % grid on;
    % [hl, hp] = boundedline(modelHandle.testAt, modelHandle.testHist(:, 2), modelHandle.testHist(:, 3), 'alpha');

    % hl.Marker = 'x';
    % hl.MarkerSize = 5;

    % hl.Color = [rand, rand, rand];
    % hp.FaceColor = hl.Color;
    % hl.LineWidth = 1.6;

    % xlabel('Traintime', 'FontSize', 12);
    % ylabel('|verg_{err}| [deg]', 'FontSize', 12);

    % % RMSE deltaMetCost
    % subplot(2, 2, 3);
    % hold on;
    % grid on;
    % plot(modelHandle.testAt, modelHandle.testHist(:, 4), 'x-', 'LineWidth', 1.3);

    % xlabel('Traintime', 'FontSize', 12);
    % ylabel('RMSE(|\Deltamc|)', 'FontSize', 12);

    % % mean, std deltaMetCost
    % subplot(2, 2, 4);
    % hold on;
    % grid on;
    % [hl, hp] = boundedline(modelHandle.testAt, modelHandle.testHist(:, 5), modelHandle.testHist(:, 6), 'alpha');

    % hl.Marker = 'x';
    % hl.MarkerSize = 5;

    % hl.Color = [rand, rand, rand];
    % hp.FaceColor = hl.Color;
    % hl.LineWidth = 1.6;

    % xlabel('Traintime', 'FontSize', 12);
    % ylabel('|\Deltamc| = |mc_{actual} - mc_{desired}|', 'FontSize', 12);

    % % Subplot overall title
    % suptitle('Test Performance vs. Traintime');

    % plotpath = sprintf('%s/testPerformanceVsTraintime', modelHandle.savePath);
    % saveas(gcf, plotpath, 'png');

    %%% Figure B
    % vergence error [deg] & MC opt approach [%] @ testing vs. iteration step

    %%% Figure C
    % muscleplain trajectories
end

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

    % member models? III member!
    fileID = fopen(strcat(savePath, '/README.txt'), 'wt' );
    fprintf(fileID, 'model w/o MetCosts: %s\n', modelWoMetCostsFullPath);
    fprintf(fileID, 'model w/  MetCosts: %s\n', modelWMetCostsFullPath);
    fclose(fileID);

    % in respect of old testHist(end) == 0 bug
    adjust = 0;
    for i = 1 : length(modelHandle)
        if (modelHandle(i).model.testHist(end, 1) == 0)
            adjust = 1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % close(figA);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Figure B
    % vergence error [deg] & MC opt approach [%] @ testing vs. iteration step
    % Total error
    % totelErrorHandle = figure();
    % hold on;
    % grid on;
    % b = boxplot(model.testResult3);

    % % remove outliers
    % outl = findobj(b,'tag','Outliers');
    % set(outl, 'Visible', 'off');

    % % rescale axis to whiskers + offset
    % upWi = findobj(b, 'tag', 'Upper Whisker');
    % lowWi = findobj(b, 'tag', 'Lower Whisker');
    % axis([0, model.testInterval + 1, ...
    %       min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
    %       max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    % if (nStim > 0)
    %     xlabel(sprintf('Iteration step (#stimuli=%d)', nStim), 'FontSize', 12);
    % else
    %     xlabel('Iteration step', 'FontSize', 12);
    % end
    % ylabel('Vergence Error [deg]', 'FontSize', 12);
    % title(sprintf('Total Vergence Error over Trial at Testing\nMean = %.4f°, Median = %.4f°,\n4*IQR = %.4f, RMSE = %.4f° at %dth step', ...
    %               mean(model.testResult3(:, model.testInterval)), median(model.testResult3(:, model.testInterval)), ...
    %               iqr(model.testResult3(:, model.testInterval)) * 4, sqrt(mean(model.testResult3(:, model.testInterval) .^ 2)), model.testInterval));

    %%%%%%%%%%%%%%%%%%%%
    figB = figure();
    hold on;
    grid on;

    steps = 3;              % show just first steps iterations & last iteration
    colors = ['b', 'r'];    % [vergErr, metCosts]

    % tmpMatrix = [vergErr_woMetCosts, vergErr_wMetCosts, metCostsApproach_woMetCosts, metCostsApproach_wMetCosts]
    tmpMatrix = horzcat(modelHandle(1).model.testResult3, modelHandle(2).model.testResult3, ...
                        modelHandle(1).model.metCostsApproach, modelHandle(2).model.metCostsApproach);

    % sort by iteration step
    idx = [];
    for (i = 1 : 20)
        idx(end + 1 : end + 4) = i : 20 : 4 * 20;
    end
    tmpMatrix = tmpMatrix(:, idx);
    tmpMatrix = [tmpMatrix(:, 1 : 4 * steps), tmpMatrix(:, end - 3 : end)];

    boxHandl = boxplot(tmpMatrix);

    % remove outliers
    outl = findobj(boxHandl, 'tag', 'Outliers');
    set(outl, 'Visible', 'off');

    % rescale axis to whiskers + offset
    upWi = findobj(boxHandl, 'tag', 'Upper Whisker');
    lowWi = findobj(boxHandl, 'tag', 'Lower Whisker');
    axis([-inf, inf, ... %0, 4 * steps + 1, ...
          min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
          max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    xlabel('Iteration step', 'FontSize', 12);
    ylabel('Vergence Error [deg]', 'FontSize', 12);
    title(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));

    plotpath = sprintf('%s/totalVergErrMetCostsApproachVsTraintimeALL', savePath);
    saveas(figB, plotpath, 'png');
    % close(figB);

    figB2 = figure();
    hold on;

    steps = 3;              % show just first steps iterations & last iteration
    colors = ['b', 'r'];    % [vergErr, metCosts]

    % tmpMatrix = [vergErr_woMetCosts, vergErr_wMetCosts, metCostsApproach_woMetCosts, metCostsApproach_wMetCosts]
    tmpMatrixVergErr = horzcat(modelHandle(1).model.testResult3, modelHandle(2).model.testResult3);
    tmpMatrixMetApp = horzcat(modelHandle(1).model.metCostsApproach, modelHandle(2).model.metCostsApproach);

    % sort by iteration step
    idx = [];
    for (i = 1 : 20)
        idx(end + 1 : end + 2) = i : 20 : 2 * 20;
    end
    tmpMatrixVergErr = tmpMatrixVergErr(:, idx);
    tmpMatrixVergErr = [tmpMatrixVergErr(:, 1 : 2 * steps), tmpMatrixVergErr(:, end - 3 : end)];
    tmpMatrixMetApp = tmpMatrixMetApp(:, idx);
    tmpMatrixMetApp = [tmpMatrixMetApp(:, 1 : 2 * steps), tmpMatrixMetApp(:, end - 3 : end)];

    sub1 = subplot(2, 1, 1);
    boxHandl = boxplot(tmpMatrixVergErr);
    grid minor;

    % remove outliers
    outl = findobj(boxHandl, 'tag', 'Outliers');
    set(outl, 'Visible', 'off');

    % rescale axis to whiskers + offset
    upWi = findobj(boxHandl, 'tag', 'Upper Whisker');
    lowWi = findobj(boxHandl, 'tag', 'Lower Whisker');
    axis([-inf, inf, ... %0, 4 * steps + 1, ...
          min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
          max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    xlabel('Iteration step', 'FontSize', 12);
    ylabel(sprintf('Vergence\nError [deg]'), 'FontSize', 12);

    sub2 = subplot(2, 1, 2);
    boxHandl2 = boxplot(tmpMatrixMetApp);
    grid minor;

    % remove outliers
    outl = findobj(boxHandl2, 'tag', 'Outliers');
    set(outl, 'Visible', 'off');

    % rescale axis to whiskers + offset
    upWi = findobj(boxHandl2, 'tag', 'Upper Whisker');
    lowWi = findobj(boxHandl2, 'tag', 'Lower Whisker');
    axis([-inf, inf, ... %0, 4 * steps + 1, ...
          min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
          max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    xlabel('Iteration step', 'FontSize', 12);
    ylabel(sprintf('Opt. Metabolic\nCosts Approach [%%]'), 'FontSize', 12);

    suptitle(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));

    plotpath = sprintf('%s/totalVergErrMetCostsApproachVsTraintimeSub2', savePath);
    saveas(figB2, plotpath, 'png');

    %%% Figure C
    % muscleplain trajectories
end

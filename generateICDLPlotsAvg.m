% Generates results plots for ICDL conference publication 2017
function generateICDLPlotsAvg(simulator, modelAt)

    % load given models
    try
        modelHandle = [load(strcat('/home/aecgroup/aecdata/Results/BmsfX5/17-03-19_1000000iter_3_gamma_0.1_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/BmsfX5/17-03-18_1000000iter_3_gamma_0.1_metCost_[0.075]', '/model.mat')); ...
                       load(strcat('/home/aecgroup/aecdata/Results/BmsfX5/17-03-19_1000000iter_4_gamma_0.1_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/BmsfX5/17-03-18_1000000iter_4_gamma_0.1_metCost_[0.075]', '/model.mat')); ...
                       load(strcat('/home/aecgroup/aecdata/Results/BmsfX5/17-03-19_1000000iter_5_gamma_0.1_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/BmsfX5/17-03-18_1000000iter_5_gamma_0.1_metCost_[0.075]', '/model.mat')); ...
                       load(strcat('/home/aecgroup/aecdata/Results/BmsfX5/17-03-19_1000000iter_6_gamma_0.1_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/BmsfX5/17-03-18_1000000iter_6_gamma_0.1_metCost_[0.075]', '/model.mat')); ...
                       load(strcat('/home/aecgroup/aecdata/Results/BmsfX5/17-03-19_1000000iter_7_gamma_0.1_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/BmsfX5/17-03-18_1000000iter_7_gamma_0.1_metCost_[0.075]', '/model.mat'))];
    catch
        error('Model(s) could not be loaded.');
    end

    % figure directory
    savePath = strcat('/home/aecgroup/aecdata/ICDLPlots', datestr(now, '/dd-mm-yy_HH:MM:SS'));
    mkdir(savePath);

    % in respect of old testHist(end) == 0 bug
    adjust = [0, 0];
    for i = 1 : size(modelHandle, 2)
        j = 1;
        for k = flip(2 : size(modelHandle(1, i).model.testHist, 1))
            if (modelHandle(1, i).model.testHist(k, 1) == 0)
                adjust(1, i) = j;
                j = j + 1;
            end
        end
    end
    if (any(adjust))
        warning('testHist contains zero entries. %d of w/o and %d of w/ model.testHist entries will be discarded.', adjust(1), adjust(2));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Figure A
    % RMSE vergence error [deg] & delta MC opt [%] @ testing vs. traintime
    % lineHandles = [];                               % [w/ metCosts, w/o metCosts]
    lineStyles = [':', '-'];                        % [w/ metCosts, w/o metCosts]
    lineWidths = [1.3, 1.3];                        % [vergErr, metCosts]

    markerStyles = ['x', 'x'];                      % [w/ metCosts, w/o metCosts]
    markerSizes = [5, 5];                           % [vergErr, metCosts]
    colors = {'b', 'r', 'b', [0, 1, 128/255]};      % median[vergErr, metCosts] IQR_patches[w/ metCosts, w/o metCosts]

    nTicks = 10;

    objRange = [0.5, 1 : 6];
    tmpModelHandle = cell(2, length(modelAt));

    % dataMatrix = {w/o[vergErrMatrix, metCostsMatrix], w/[vergErrMatrix, metCostsMatrix]}
    % vergErrMatrix = metCostsMatrix = objDist * VSE/0° * nStim * testInterval x testAt
    dataMatrix = {{zeros(1680 * size(modelHandle, 1), length(modelAt)), []}, {zeros(1680 * size(modelHandle, 1), length(modelAt)), []}};
    dataMatrixEnd = {{zeros(1680 * size(modelHandle, 1), 20), []}, {zeros(1680 * size(modelHandle, 1), 20), []}}; % final values @ modelAt = 1mio & iter = [1, 20]

    at0Matrix = {{zeros(1680, 1), []}, {zeros(1680, 1), []}};
    iqrLine = zeros(4, length(modelAt) + 1); % +modelAt0

    % extract all relevant data from all sub-experiments
    for i = 1 : size(modelHandle, 2)
        for seedIter = 1 : size(modelHandle, 1)
            for trainedUntil = 1 : length(modelAt)
                try
                    subFolder = sprintf('modelAt%d', modelAt(trainedUntil));
                    tmpModelHandle{i, trainedUntil} = load(sprintf('%s/%s/model.mat', modelHandle(seedIter, i).model.savePath, subFolder));
                catch
                   % catch case when (sub-)experiment started, but has no test results yet
                   error('%s/%s/model.mat\ncould not be loaded.', modelHandle(seedIter, i).model.savePath, subFolder);
                end

                % fill data matrix & exclude VSA = 0° trials
                nStim = size(tmpModelHandle{i, trainedUntil}.model.testResult3, 1) / (length(objRange) * 7);
                currIdx = 1 + 1680 * (seedIter - 1);
                startInd = 1;
                startInd0 = nStim * 3 + 1;
                for j = 1 : length(objRange)
                    colSize1 = length(tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : startInd0 - 1, tmpModelHandle{i, trainedUntil}.model.testInterval));
                    dataMatrix{i}{1}(currIdx : currIdx + colSize1 - 1, trainedUntil) = tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : startInd0 - 1, ...
                                                                                                      tmpModelHandle{i, trainedUntil}.model.testInterval);

                    % end matrix creation
                    if (trainedUntil == length(modelAt))
                        for k = 1 : 20
                            dataMatrixEnd{i}{1}(currIdx : currIdx + colSize1 - 1, k) = tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : startInd0 - 1, k);
                        end
                    end

                    currIdx = currIdx + colSize1;
                    endInd0 = startInd0 + nStim - 1;
                    startInd = endInd0 + 1;
                    startInd0 = endInd0 + nStim * 6 + 1;
                end
                % concatinate remainder
                colSize2 = length(tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : end, tmpModelHandle{i, trainedUntil}.model.testInterval));
                dataMatrix{i}{1}(currIdx : currIdx + colSize2 - 1, trainedUntil) = tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : end, ...
                                                                                                  tmpModelHandle{i, trainedUntil}.model.testInterval);

                % delta metCosts
                dataMatrix{i}{2}(1 + 1960 * (seedIter - 1) : 1960 * seedIter, trainedUntil) = tmpModelHandle{i, trainedUntil}.model.testResult7(:, tmpModelHandle{i, trainedUntil}.model.testInterval);

                % end matrix creation
                if (trainedUntil == length(modelAt))
                    for k = 1 : 20
                        dataMatrixEnd{i}{1}(currIdx : currIdx + colSize2 - 1, k) = tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : end, k);
                    end
                    dataMatrixEnd{i}{2}(1 + 1960 * (seedIter - 1) : 1960 * seedIter, :) = tmpModelHandle{i, trainedUntil}.model.testResult7;
                end
            end
        end

        % add modelAt0 entries
        if (i == 1)
            hm = load('/home/aecgroup/aecdata/Results/17-03-08_300000iter_1_newStandard_0,3Mio/modelAt0/model.mat');
        end

        % fill data matrix & exclude VSA = 0° trials
        nStim = size(hm.model.testResult3, 1) / (length(objRange) * 7);
        currIdx = 1;
        startInd = 1;
        startInd0 = nStim * 3 + 1;
        for j = 1 : length(objRange)
            colSize1 = length(hm.model.testResult3(startInd : startInd0 - 1, hm.model.testInterval));
            at0Matrix{i}{1}(currIdx : currIdx + colSize1 - 1) = hm.model.testResult3(startInd : startInd0 - 1, hm.model.testInterval);
            currIdx = currIdx + colSize1;
            endInd0 = startInd0 + nStim - 1;
            startInd = endInd0 + 1;
            startInd0 = endInd0 + nStim * 6 + 1;
        end
        % concatinate remainder
        colSize2 = length(hm.model.testResult3(startInd : end, hm.model.testInterval));
        at0Matrix{i}{1}(currIdx : currIdx + colSize2 - 1) = hm.model.testResult3(startInd : end, hm.model.testInterval);

        % delta metCosts
        at0Matrix{i}{2} = hm.model.testResult7(:, hm.model.testInterval);

        % concatinate both matricies
        dataMatrix{i}{1} = horzcat(repmat(at0Matrix{i}{1}, [size(modelHandle, 1), 1]), dataMatrix{i}{1});
        dataMatrix{i}{2} = horzcat(repmat(at0Matrix{i}{2}, [size(modelHandle, 1), 1]), dataMatrix{i}{2});

        % extract IQR edge coordinates
        tmpFig = figure();
        boxHandle = boxplot(dataMatrix{i}{1});
        upWi = findobj(boxHandle, 'tag', 'Upper Whisker');
        lowWi = findobj(boxHandle, 'tag', 'Lower Whisker');
        iqrLine(i * 2 - 1 : i * 2, :) = [arrayfun(@(x) x.YData(1), upWi)'; arrayfun(@(x) x.YData(2), lowWi)'];
        close(tmpFig);
    end

    figA = figure();
    hold on;
    modelAt = horzcat(0, modelAt);

    for i = 1 : size(modelHandle, 2)
        if (i == 1)
            ax1 = gca; % current axes
            ax1.YColor = colors{1};
            ax1.XAxis.Label.String = 'Traintime';
            ax1.XAxis.Label.FontSize = 12;
            ax1.Title.String = 'Test Performance & Metabolic Costs vs. Traintime';
            % ax1.YAxis.Label.String = 'RMSE(verg_{err}) [deg]';
            % ax1.YAxis.Label.String = sprintf('Vergence Error\nw.r.t. Opt. median [%%]');
            ax1.YAxis.Label.String = 'median(verg_{err}) [deg]';
            ax1.YAxis.Label.FontSize = 12;

            ax2 = axes('Position', ax1.Position, ...
                       'YAxisLocation', 'right', ...
                       'Color', 'none');
        end

        % fill area defined by upper & lower IQR bounds
        if (i == 1)
            hp1 = patch([modelAt, flip(modelAt)], [iqrLine(i * 2 - 1, :), flip(iqrLine(i * 2, :))], ...
                        colors{3 + i - 1}, 'LineStyle', 'none', 'FaceAlpha', 0.2);
            hp1.Parent = ax1;

            % color trick -> black entries in legend
            hl1 = plot(ax1, modelAt, median(dataMatrix{i}{1}), ...
                         'LineStyle', lineStyles(i), 'Marker', markerStyles(i), 'MarkerSize', markerSizes(1), 'Color', 'k', 'LineWidth', lineWidths(1));
        else
            hp2 = patch([modelAt, flip(modelAt)], [iqrLine(i * 2 - 1, :), flip(iqrLine(i * 2, :))], ...
                        colors{3 + i - 1}, 'LineStyle', 'none', 'FaceAlpha', 0.3);
            hp2.Parent = ax1;

            % color trick -> black entries in legend
            hl2 = plot(ax1, modelAt, median(dataMatrix{i}{1}), ...
                         'LineStyle', lineStyles(i), 'Marker', markerStyles(i), 'MarkerSize', markerSizes(1), 'Color', 'k', 'LineWidth', lineWidths(1));
        end
        hold on;

        hl3 = plot(ax1, modelAt, median(dataMatrix{i}{1}), ...
                      'LineStyle', lineStyles(i), 'Marker', markerStyles(i), 'MarkerSize', markerSizes(1), 'Color', colors{1}, 'LineWidth', lineWidths(1));
        hold on;

        hl4 = plot(ax2, modelAt, median(dataMatrix{i}{2}), ...
                  'LineStyle', lineStyles(i), 'Marker', markerStyles(i), 'MarkerSize', markerSizes(2), 'Color', colors{2}, 'LineWidth', lineWidths(2));
        hold on;

        if (i == 1)
            % |\DeltaMC_{opt}| = |MC_{actual} - MC_{optimal}| / |MC_{start} - MC_{optimal}|
            ax2.YAxis.Label.String = 'median(\Deltamet. costs) [J]';
            ax2.YAxis.Label.FontSize = 12;
            ax2.YColor = colors{2};
            % ax2.YAxis.Label.Rotation = -90;
            ax2.YAxisLocation = 'right';
        end
    end

    grid(ax1, 'on');
    % ax1.YMinorGrid = 'on';

    % ax1.YAxis.Limits = [-0.05, 0.1];
    % ax2.YAxis.Limits = [-0.2, 1.1];

    % set #nTicks ticks for y-axis
    % set(ax1, 'YTick', round(linspace(ax1.YAxis.Limits(1), ax1.YAxis.Limits(2), nTicks - 2), 2));
    % set(ax2, 'YTick', round(linspace(ax2.YAxis.Limits(1), ax2.YAxis.Limits(2), nTicks), 2));

    % gKey = {'w/o met. costs', 'w/  met. costs'};
    % l = gridLegend(lineHandles, 1, gKey, 'Location', 'southwest');
    %     %, 'Orientation', 'Horizontal', 'Location', 'southoutside', 'Fontsize', 8);

    lineHandles = [hl1, hp1, hl2, hp2];
    gKey = {'w/o met. costs', 'IQR', ...
            'w/  met. costs', 'IQR'};

    % l = gridLegend(lineHandles, 2, gKey, 'Location', 'southwest');
        %, 'Orientation', 'Horizontal', 'Location', 'southoutside', 'Fontsize', 8);
    l = legend(lineHandles, gKey, 'Orientation', 'horizontal', 'Location', 'south');
    % l.Box = 'off';

    % ax2.YAxis.Label.Position(1) = ax2.YAxis.Label.Position(1) * 1.5;

    plotpath = sprintf('%s/FigA_VergErrMetCostsVsTraintime', savePath);
    saveas(figA, plotpath, 'png');
    close(figA);

    % member models? III member!
    % TODO: add into file: IQR & median of vergerr and metcosts @ 20th iteration
    fileID = fopen(strcat(savePath, '/README.txt'), 'at' );
    fprintf(fileID, 'model w/o MetCosts: %s\n', modelHandle(1, 1).model.savePath);
    fprintf(fileID, 'model w/  MetCosts: %s\n\n', modelHandle(1, 2).model.savePath);

    fprintf(fileID, '========================================================================================================\n');
    fprintf(fileID, 'model w/o MetCosts\n');
    fprintf(fileID, 'modelAt: %s\n\n', int2str(modelAt));

    fprintf(fileID, 'figA vergErr median:\t%f %f %f %f %f %f\n', median(dataMatrix{1}{1}));
    fprintf(fileID, 'figA vergErr iqr:\t%f %f %f %f %f %f\n', iqr(dataMatrix{1}{1}));
    fprintf(fileID, 'figA vergErr mean:\t%f %f %f %f %f %f\n', mean(dataMatrix{1}{1}));
    fprintf(fileID, 'figA vergErr std:\t%f %f %f %f %f %f\n\n', std(dataMatrix{1}{1}));

    fprintf(fileID, 'figA metCosts median:\t%f %f %f %f %f %f\n', median(dataMatrix{1}{2}));
    fprintf(fileID, 'figA metCosts iqr:\t%f %f %f %f %f %f\n', iqr(dataMatrix{1}{2}));
    fprintf(fileID, 'figA metCosts mean:\t%f %f %f %f %f %f\n', mean(dataMatrix{1}{2}));
    fprintf(fileID, 'figA metCosts std:\t%f %f %f %f %f %f\n\n', std(dataMatrix{1}{2}));

    fprintf(fileID, 'model w/ MetCosts\n');
    fprintf(fileID, 'modelAt: %s\n\n', int2str(modelAt));

    fprintf(fileID, 'figA vergErr median:\t%f %f %f %f %f %f\n', median(dataMatrix{2}{1}));
    fprintf(fileID, 'figA vergErr iqr:\t%f %f %f %f %f %f\n', iqr(dataMatrix{2}{1}));
    fprintf(fileID, 'figA vergErr mean:\t%f %f %f %f %f %f\n', mean(dataMatrix{2}{1}));
    fprintf(fileID, 'figA vergErr std:\t%f %f %f %f %f %f\n\n', std(dataMatrix{2}{1}));

    fprintf(fileID, 'figA metCosts median:\t%f %f %f %f %f %f\n', median(dataMatrix{2}{2}));
    fprintf(fileID, 'figA metCosts iqr:\t%f %f %f %f %f %f\n', iqr(dataMatrix{2}{2}));
    fprintf(fileID, 'figA metCosts mean:\t%f %f %f %f %f %f\n', mean(dataMatrix{2}{2}));
    fprintf(fileID, 'figA metCosts std:\t%f %f %f %f %f %f\n\n', std(dataMatrix{2}{2}));

    fclose(fileID);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% Figure B
    % %% all the data in a single plot

    % % figB0 = figure();
    % % hold on;
    % % grid on;

    % % steps = 3;              % show just first steps iterations & last iteration
    % % % colors = ['b', 'r'];    % [vergErr, metCosts]

    % % % tmpMatrix = [vergErr_woMetCosts, vergErr_wMetCosts, metCostsApproach_woMetCosts, metCostsApproach_wMetCosts]
    % % tmpMatrix = horzcat(modelHandle(1).model.testResult3, modelHandle(2).model.testResult3, ...
    % %                     modelHandle(1).model.metCostsApproach, modelHandle(2).model.metCostsApproach);

    % % % sort by iteration step
    % % idx = [];
    % % for (i = 1 : 20)
    % %     idx(end + 1 : end + 4) = i : 20 : 4 * 20;
    % % end
    % % tmpMatrix = tmpMatrix(:, idx);
    % % tmpMatrix = [tmpMatrix(:, 1 : 4 * steps), tmpMatrix(:, end - 3 : end)];

    % % boxHandl = boxplot(tmpMatrix);

    % % % remove outliers
    % % outl = findobj(boxHandl, 'tag', 'Outliers');
    % % set(outl, 'Visible', 'off');

    % % % rescale axis to whiskers + offset
    % % upWi = findobj(boxHandl, 'tag', 'Upper Whisker');
    % % lowWi = findobj(boxHandl, 'tag', 'Lower Whisker');
    % % axis([-inf, inf, ... %0, 4 * steps + 1, ...
    % %       min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
    % %       max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    % % xlabel('Iteration step', 'FontSize', 12);
    % % ylabel('Vergence Error [deg]', 'FontSize', 12);
    % % title(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));

    % % plotpath = sprintf('%s/FigB0_VergErrMetCostsApproachVsTestIterALL', savePath);
    % % saveas(figB0, plotpath, 'png');
    % % close(figB0);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %% data separated into 2 subplots
    % figB1 = figure();
    % hold on;

    % steps = 3;              % show just first steps iterations & last iteration
    % colors = {[0, 100/255, 200/255], [0, 95/255, 0]};    % [w/o metcosts, w/ metcosts] for boxes
    % captions = cell(1, 2);
    % captions{1} = 'w/o met. costs';
    % captions{2} = 'w/ met. costs';

    % % tmpMatrix = [vergErr_woMetCosts, vergErr_wMetCosts, metCostsApproach_woMetCosts, metCostsApproach_wMetCosts]
    % tmpMatrixVergErr = vertcat(horzcat(modelHandle(1, 1).model.vergenceAngleApproach, modelHandle(1, 2).model.vergenceAngleApproach), ...
    %                            horzcat(modelHandle(2, 1).model.vergenceAngleApproach, modelHandle(2, 2).model.vergenceAngleApproach), ...
    %                            horzcat(modelHandle(3, 1).model.vergenceAngleApproach, modelHandle(3, 2).model.vergenceAngleApproach), ...
    %                            horzcat(modelHandle(4, 1).model.vergenceAngleApproach, modelHandle(4, 2).model.vergenceAngleApproach), ...
    %                            horzcat(modelHandle(5, 1).model.vergenceAngleApproach, modelHandle(5, 2).model.vergenceAngleApproach));

    % tmpMatrixMetApp = vertcat(horzcat(modelHandle(1, 1).model.metCostsApproach, modelHandle(1, 2).model.metCostsApproach), ...
    %                           horzcat(modelHandle(2, 1).model.metCostsApproach, modelHandle(2, 2).model.metCostsApproach), ...
    %                           horzcat(modelHandle(3, 1).model.metCostsApproach, modelHandle(3, 2).model.metCostsApproach), ...
    %                           horzcat(modelHandle(4, 1).model.metCostsApproach, modelHandle(4, 2).model.metCostsApproach), ...
    %                           horzcat(modelHandle(5, 1).model.metCostsApproach, modelHandle(5, 2).model.metCostsApproach));

    % % sort by iteration step
    % idx = [];
    % for (i = 1 : 20)
    %     idx(end + 1 : end + 2) = i : 20 : 2 * 20;
    % end
    % tmpMatrixVergErr = tmpMatrixVergErr(:, idx);
    % tmpMatrixVergErr = [tmpMatrixVergErr(:, 1 : 2 * steps), tmpMatrixVergErr(:, end - 1 : end)];
    % tmpMatrixMetApp = tmpMatrixMetApp(:, idx);
    % tmpMatrixMetApp = [tmpMatrixMetApp(:, 1 : 2 * steps), tmpMatrixMetApp(:, end - 1 : end)];

    % sub1 = subplot(2, 1, 1);
    % pos = [1 1.2 1.5 1.7 2 2.2 2.5 2.7];
    % boxHandl = boxplot(tmpMatrixVergErr, 'labels', {'1','','2','','3','','20',''}, 'positions', pos);
    % tmpHandle = findobj(boxHandl, 'type', 'text');
    % set(tmpHandle, 'Interpreter', 'tex');
    % grid minor;

    % subBoxHandl = findobj(gca,'Tag','Box');
    % % subBoxHandl = findobj(boxHandl,'Tag','Box');

    % boxesArray = findobj(boxHandl);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray(idx2(j)).Color = colors{2};
    %         end
    %     end
    % end

    % % remove outliers
    % outl = findobj(boxHandl, 'tag', 'Outliers');
    % set(outl, 'Visible', 'off');

    % % rescale axis to whiskers + offset
    % upWi = findobj(boxHandl, 'tag', 'Upper Whisker');
    % lowWi = findobj(boxHandl, 'tag', 'Lower Whisker');
    % axis([0.9, 2.8, ...
    %       min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
    %       max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    % % xlabel('Iteration step', 'FontSize', 12);
    % % ylabel(sprintf('Vergence\nError [deg]'), 'FontSize', 12);
    % ylabel(sprintf('Vergence Error\nReduction [%%]'), 'FontSize', 12);

    % %% put ylabel right and rotate text
    % % set(sub1, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('Metabolic Costs\nReduction [%%]'), 'rot', -90, 'FontSize', 12);
    % % p = get(lh, 'position');
    % % set(sub1,'yaxislocation','left');
    % % set(lh,'position', p);

    % sub2 = subplot(2, 1, 2);
    % boxHandl2 = boxplot(tmpMatrixMetApp, 'labels', {'1','','2','','3','','20',''}, 'positions', pos);
    % grid minor;

    % boxesArray = findobj(boxHandl2);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray(idx2(j)).Color = colors{2};
    %         end
    %     end
    % end

    % % remove outliers
    % outl = findobj(boxHandl2, 'tag', 'Outliers');
    % set(outl, 'Visible', 'off');

    % % rescale axis to whiskers + offset
    % upWi = findobj(boxHandl2, 'tag', 'Upper Whisker');
    % lowWi = findobj(boxHandl2, 'tag', 'Lower Whisker');
    % axis([0.9, 2.8, ...
    %       min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
    %       max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    % xlabel('Iteration step', 'FontSize', 12);
    % % ylabel(sprintf('Opt. Metabolic\nCosts Approach [%%]'), 'FontSize', 12);
    % ylabel(sprintf('Metabolic Costs\nReduction [%%]'), 'FontSize', 12);

    % %% put ylabel right and rotate text
    % % set(sub2, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('Metabolic Costs\nReduction [%%]'), 'rot', -90, 'FontSize', 12);
    % % p = get(lh, 'position');
    % % set(sub2, 'yaxislocation', 'left');
    % % set(lh, 'position', p);

    % % suptitle(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));
    % suptitle(sprintf('Reduction of Vergence Error & Metabolic Costs\nvs. Iteration at Testing'));

    % l = legend(subBoxHandl([2, 1]), captions);
    % l.FontSize = 7;
    % l.Orientation = 'horizontal';
    % l.Location = 'southoutside';

    % %% repositioning subfigures
    % sub1.Position(3 : 4) = sub2.Position(3 : 4);
    % % sub1.Position(2) = sub1.Position(2) * 0.95;
    % % sub2.Position(2) = sub2.Position(2) * 0.9;
    % l.Position(2) = l.Position(2) * 1.075;

    % plotpath = sprintf('%s/FigB1_VergErrMetCostsApproachVsTestIter', savePath);
    % saveas(figB1, plotpath, 'png');
    % close(figB1);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % alternative data

    % figB2 = figure();
    % hold on;

    % steps = 3;              % show just first steps iterations & last iteration
    % colors = {[0, 100/255, 200/255], [0, 95/255, 0]};    % [w/o metcosts, w/ metcosts] for boxes
    % captions = cell(1, 2);
    % captions{1} = 'w/o met. costs';
    % captions{2} = 'w/ met. costs';

    % % % tmpMatrix = [vergErr_woMetCosts, vergErr_wMetCosts, metCostsApproach_woMetCosts, metCostsApproach_wMetCosts]
    % % % tmpMatrixVergErr = horzcat(modelHandle(1).model.testResult3, modelHandle(2).model.testResult3);
    % tmpMatrixVergErr = horzcat(dataMatrixEnd{1}{1}, dataMatrixEnd{2}{1});
    % tmpMatrixMetApp = horzcat(dataMatrixEnd{1}{2}, dataMatrixEnd{2}{2});

    % % % sort by iteration step
    % idx = [];
    % for (i = 1 : 20)
    %     idx(end + 1 : end + 2) = [i, i + 20];
    % end
    % tmpMatrixVergErr = tmpMatrixVergErr(:, idx);
    % tmpMatrixVergErr = [tmpMatrixVergErr(:, 1 : 2 * steps), tmpMatrixVergErr(:, end - 1 : end)];
    % tmpMatrixMetApp = tmpMatrixMetApp(:, idx);
    % tmpMatrixMetApp = [tmpMatrixMetApp(:, 1 : 2 * steps), tmpMatrixMetApp(:, end - 1 : end)];

    % sub1 = subplot(2, 1, 1);
    % % pos = [1 1.33 2 2.33 3 3.33 4 4.33];
    % pos = [1 1.2 1.5 1.7 2 2.2 2.5 2.7];
    % boxHandl = boxplot(tmpMatrixVergErr, 'labels', {'1','','2','','3','','20',''}, 'positions', pos);
    % tmpHandle = findobj(boxHandl, 'type', 'text');
    % set(tmpHandle, 'Interpreter', 'tex');
    % grid minor;

    % subBoxHandl = findobj(gca,'Tag','Box');
    % % subBoxHandl = findobj(boxHandl,'Tag','Box');

    % boxesArray = findobj(boxHandl);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray(idx2(j)).Color = colors{2};
    %         end
    %     end
    % end

    % % remove outliers
    % outl = findobj(boxHandl, 'tag', 'Outliers');
    % set(outl, 'Visible', 'off');

    % % rescale axis to whiskers + offset
    % upWi = findobj(boxHandl, 'tag', 'Upper Whisker');
    % lowWi = findobj(boxHandl, 'tag', 'Lower Whisker');
    % axis([0.9, 2.8, ...
    %       min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
    %       max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    % ylabel('verg_{err} [deg]', 'FontSize', 12);

    % %% put ylabel right and rotate text
    % % set(sub1, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('Metabolic Costs\nReduction [%%]'), 'rot', -90, 'FontSize', 12);
    % % p = get(lh, 'position');
    % % set(sub1,'yaxislocation','left');
    % % set(lh,'position', p);

    % sub2 = subplot(2, 1, 2);
    % boxHandl2 = boxplot(tmpMatrixMetApp, 'labels', {'1','','2','','3','','20',''}, 'positions', pos);
    % grid minor;

    % boxesArray = findobj(boxHandl2);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray(idx2(j)).Color = colors{2};
    %         end
    %     end
    % end

    % % remove outliers
    % outl = findobj(boxHandl2, 'tag', 'Outliers');
    % set(outl, 'Visible', 'off');

    % % rescale axis to whiskers + offset
    % upWi = findobj(boxHandl2, 'tag', 'Upper Whisker');
    % lowWi = findobj(boxHandl2, 'tag', 'Lower Whisker');
    % axis([0.9, 2.8, ...
    %       min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
    %       max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    % xlabel('Iteration step', 'FontSize', 12);
    % ylabel('\Deltamet. costs [J]', 'FontSize', 12);

    % %% put ylabel right and rotate text
    % % set(sub2, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('Metabolic Costs\nReduction [%%]'), 'rot', -90, 'FontSize', 12);
    % % p = get(lh, 'position');
    % % set(sub2, 'yaxislocation', 'left');
    % % set(lh, 'position', p);

    % % suptitle(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));
    % suptitle(sprintf('Reduction of Vergence Error & Metabolic Costs\nvs. Iteration at Testing'));

    % % l = legend(subBoxHandl([2, 1]), captions);
    % % l.FontSize = 7;
    % % l.Orientation = 'horizontal';
    % % l.Location = 'southoutside';

    % [l, objh, ~, ~] = legend(subBoxHandl([2, 1]), captions, 'Orientation', 'horizontal', 'Location', 'southoutside');
    % set(objh, 'linewidth', 2);
    % l.Position(2) = 0.465;

    % %% repositioning subfigures
    % sub1.Position(3 : 4) = sub2.Position(3 : 4);
    % sub1.Position(2) = 0.6;

    % plotpath = sprintf('%s/FigB2_VergErrMetCostsVsTestIter', savePath);
    % saveas(figB2, plotpath, 'png');
    % close(figB2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % combines plot

    figB3 = figure();
    hold on;

    steps = 3;                                          % show just first steps iterations & last iteration
    colors = {[0, 100/255, 200/255], [0, 95/255, 0]};   % [w/o metcosts, w/ metcosts] for boxes
    captions = cell(1, 2);
    captions{1} = 'w/o met. costs';
    captions{2} = 'w/ met. costs';

    % tmpMatrix = [vergErr_woMetCosts, vergErr_wMetCosts, metCostsApproach_woMetCosts, metCostsApproach_wMetCosts]
    % tmpMatrixVergErr = horzcat(modelHandle(1).model.testResult3, modelHandle(2).model.testResult3);
    tmpMatrixVergErr = horzcat(dataMatrixEnd{1}{1}, dataMatrixEnd{2}{1});

    % sort by iteration step
    idx = [];
    for (i = 1 : 20)
        idx(end + 1 : end + 2) = [i, i + 20];
    end
    tmpMatrixVergErr = tmpMatrixVergErr(:, idx);
    tmpMatrixVergErr = [tmpMatrixVergErr(:, 1 : 2 * steps), tmpMatrixVergErr(:, end - 1 : end)];

    tmpMatrixMetApp = vertcat(horzcat(modelHandle(1, 1).model.metCostsApproach, modelHandle(1, 2).model.metCostsApproach), ...
                              horzcat(modelHandle(2, 1).model.metCostsApproach, modelHandle(2, 2).model.metCostsApproach), ...
                              horzcat(modelHandle(3, 1).model.metCostsApproach, modelHandle(3, 2).model.metCostsApproach), ...
                              horzcat(modelHandle(4, 1).model.metCostsApproach, modelHandle(4, 2).model.metCostsApproach), ...
                              horzcat(modelHandle(5, 1).model.metCostsApproach, modelHandle(5, 2).model.metCostsApproach));

    % sort by iteration step
    idx = [];
    for (i = 1 : 20)
        idx(end + 1 : end + 2) = i : 20 : 2 * 20;
    end
    tmpMatrixMetApp = tmpMatrixMetApp(:, idx);
    tmpMatrixMetApp = [tmpMatrixMetApp(:, 1 : 2 * steps), tmpMatrixMetApp(:, end - 1 : end)];

    sub1 = subplot(2, 1, 1);
    pos = [1 1.2 1.5 1.7 2 2.2 2.5 2.7];
    boxHandl = boxplot(tmpMatrixVergErr, 'labels', {'1','','2','','3','','20',''}, 'positions', pos);
    tmpHandle = findobj(boxHandl, 'type', 'text');
    set(tmpHandle, 'Interpreter', 'tex');
    grid minor;

    subBoxHandl = findobj(gca,'Tag','Box');
    % subBoxHandl = findobj(boxHandl,'Tag','Box');

    boxesArray = findobj(boxHandl);
    for i = 1 : size(tmpMatrixVergErr, 2)
        idx2 = (1 : 7) + (i - 1) * 7;
        idx2(6 : 7) = [];
        if (mod(i, 2) == 1)
            for j = 1 : length(idx2)
                boxesArray(idx2(j)).Color = colors{1};
            end
        else
            for j = 1 : length(idx2)
                boxesArray(idx2(j)).Color = colors{2};
            end
        end
    end

    % remove outliers
    outl = findobj(boxHandl, 'tag', 'Outliers');
    set(outl, 'Visible', 'off');

    % rescale axis to whiskers + offset
    upWi = findobj(boxHandl, 'tag', 'Upper Whisker');
    lowWi = findobj(boxHandl, 'tag', 'Lower Whisker');
    axis([0.9, 2.8, ...
          min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
          max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    ylabel('verg_{err} [deg]', 'FontSize', 12);

    %% put ylabel right and rotate text
    % set(sub1, 'yaxislocation', 'right');
    % lh = ylabel(sprintf('Metabolic Costs\nReduction [%%]'), 'rot', -90, 'FontSize', 12);
    % p = get(lh, 'position');
    % set(sub1,'yaxislocation','left');
    % set(lh,'position', p);

    sub2 = subplot(2, 1, 2);
    boxHandl2 = boxplot(tmpMatrixMetApp, 'labels', {'1','','2','','3','','20',''}, 'positions', pos);
    grid minor;

    boxesArray = findobj(boxHandl2);
    for i = 1 : size(tmpMatrixVergErr, 2)
        idx2 = (1 : 7) + (i - 1) * 7;
        idx2(6 : 7) = [];
        if (mod(i, 2) == 1)
            for j = 1 : length(idx2)
                boxesArray(idx2(j)).Color = colors{1};
            end
        else
            for j = 1 : length(idx2)
                boxesArray(idx2(j)).Color = colors{2};
            end
        end
    end

    % remove outliers
    outl = findobj(boxHandl2, 'tag', 'Outliers');
    set(outl, 'Visible', 'off');

    % rescale axis to whiskers + offset
    upWi = findobj(boxHandl2, 'tag', 'Upper Whisker');
    lowWi = findobj(boxHandl2, 'tag', 'Lower Whisker');
    axis([0.9, 2.8, ...
          min(arrayfun(@(x) x.YData(1), lowWi)) + min(arrayfun(@(x) x.YData(1), lowWi)) * 0.1, ...
          max(arrayfun(@(x) x.YData(2), upWi)) * 1.1]);

    xlabel('Iteration step', 'FontSize', 12);
    ylabel(sprintf('Metabolic Costs\nReduction [%%]'), 'FontSize', 12);

    %% put ylabel right and rotate text
    % set(sub2, 'yaxislocation', 'right');
    % lh = ylabel(sprintf('Metabolic Costs\nReduction [%%]'), 'rot', -90, 'FontSize', 12);
    % p = get(lh, 'position');
    % set(sub2, 'yaxislocation', 'left');
    % set(lh, 'position', p);

    % suptitle(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));
    suptitle(sprintf('Reduction of Vergence Error & Metabolic Costs\nvs. Iteration at Testing'));

    [l, objh, ~, ~] = legend(subBoxHandl([2, 1]), captions, 'Orientation', 'horizontal', 'Location', 'southoutside');
    set(objh, 'linewidth', 2);

    %% repositioning subfigures
    sub1.Position(3 : 4) = sub2.Position(3 : 4);
    sub1.Position(2) = 0.6;
    l.Position(2) = 0.465;

    plotpath = sprintf('%s/FigB3_VergErrMetCostsVsTestIter', savePath);
    saveas(figB3, plotpath, 'png');
    close(figB3);

    fileID = fopen(strcat(savePath, '/README.txt'), 'at' );
    fprintf(fileID, '========================================================================================================\n');
    fprintf(fileID, 'testIter: %s %s\n\n', int2str([1, 1, 2, 2, 3, 3, 20, 20]), '= [w/o MetCosts, w/ MetCosts, w/o MetCosts, w/ MetCosts, ...]');

    fprintf(fileID, 'figB3 vergErr median:\t%f %f %f %f %f %f %f %f\n', median(tmpMatrixVergErr));
    fprintf(fileID, 'figB3 vergErr iqr:\t%f %f %f %f %f %f %f %f\n', iqr(tmpMatrixVergErr));
    fprintf(fileID, 'figB3 vergErr mean:\t%f %f %f %f %f %f %f %f\n', mean(tmpMatrixVergErr));
    fprintf(fileID, 'figB3 vergErr std:\t%f %f %f %f %f %f %f %f\n\n', std(tmpMatrixVergErr));

    fprintf(fileID, 'figB3 metCosts median:\t%f %f %f %f %f %f %f %f\n', median(tmpMatrixMetApp));
    fprintf(fileID, 'figB3 metCosts iqr:\t%f %f %f %f %f %f %f %f\n', iqr(tmpMatrixMetApp));
    fprintf(fileID, 'figB3 metCosts mean:\t%f %f %f %f %f %f %f %f\n', mean(tmpMatrixMetApp));
    fprintf(fileID, 'figB3 metCosts std:\t%f %f %f %f %f %f %f %f\n', std(tmpMatrixMetApp));

    fclose(fileID);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Figure C
    % muscleplain trajectories

    % simulator check
    if (isempty(simulator))
        simulator = prepareSimulator([]);
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
    rng(667);
    objDist = [0.5, 6];
    startVergErr = [-2, 2];
    numIters = 20;
    stimuliIndices = [1];

    % initMethod elem. {0, 1}; 0 = getMFedoodD, 1 = fixed
    initMethod = uint8(0);

    % hand-picked inits for medial rectus for initMethod = 0
    % mrVal := {'w/o met. costs', 'w/  met. costs'} x objDist * startVergErr
    mrVal = [[0.055, 0.085, 0.055, 0.085]; [0.03, 0.055, 0.03, 0.05]];
    % mrVal = [[0.04, 0.05, 0.04, 0.05]; [0.03, 0.04, 0.03, 0.04]];

    if ((initMethod == 0) && (size(mrVal, 1) ~= size(modelHandle, 2)))
        error('It must hold size(mrVal, 1) = %d = size(modelHandle, 2) = %d', ...
              size(mrVal, 1), size(modelHandle, 2));
    elseif ((initMethod == 0) && (size(mrVal, 2) ~= length(objDist) * length(startVergErr)))
        error('It must hold size(mrVal, 2) = %d = length(objDist) * length(startVergErr) = %d', ...
              size(mrVal, 2), length(objDist) + length(startVergErr));
    end

    % hand-picked inits for muscles for initMethod = 1
    % cmdInit = [[0.03; 0.16], [0.05; 0.12], [0.07; 0.12], [0.04; 0.08], [0.06; 0.06], [0.08; 0.06]];
    cmdInit = [[0.075; 0.12], [0.075; 0.12], [0.075; 0.12], [0.06; 0.1], [0.06; 0.1], [0.06; 0.1]];

    if ((initMethod == 1) && (size(cmdInit, 2) ~= length(objDist) + length(stimuliIndices) + length(startVergErr)))
        error('It must hold size(cmdInit, 2) = %d = length(objDist) + length(stimuliIndices) + length(startVergErr) = %d', ...
              size(cmdInit, 2), length(objDist) + length(stimuliIndices) + length(startVergErr));
    end

    nStimuli = length(stimuliIndices);
    trajectory = zeros(length(objDist), length(startVergErr), nStimuli, numIters + 1, 2);
    trajectory2 = zeros(length(objDist), length(startVergErr), nStimuli, numIters + 1, 2);

    %% Data generation
    % vergErrMax = 2;
    % angleMin = (atand(modelHandle(1, 1).model.baseline / (2 * modelHandle(1, 1).model.objDistMax)) * 2) - vergErrMax; %angle for both eyes
    % angleMax = (atand(modelHandle(1, 1).model.baseline / (2 * modelHandle(1, 1).model.objDistMin)) * 2) + vergErrMax;
    angleMinT = 0;
    angleMaxT = (atand(modelHandle(1, 1).model.baseline / (2 * modelHandle(1, 1).model.objDistMax)) * 2) + 6;
    iter1 = 1;
    for odIndex = 1 : length(objDist)
        angleDes = 2 * atand(modelHandle(1, 1).model.baseline / (2 * objDist(odIndex)));

        for stimIter = 1 : nStimuli
            currentTexture = stimuliIndices(stimIter);

            for vergErrIndex = 1 : length(startVergErr)
                % muscle init
                if (initMethod == 0)
                    % catch negative/divergent vergence angle
                    % [command, angleNew] = modelHandle(1, 1).model.getMFedood(objDist(odIndex), min(startVergErr(vergErrIndex), modelHandle(1, 1).model.getVergErrMax(objDist(odIndex))));
                    [command, angleNew] = modelHandle(1, 1).model.getMFedoodD(objDist(odIndex), ...
                                                                           min(startVergErr(vergErrIndex), modelHandle(1, 1).model.getVergErrMax(objDist(odIndex))), ...
                                                                           mrVal(1, odIndex * 2 - 1 + vergErrIndex - 1));
                elseif (initMethod == 1)
                    command = cmdInit(:, iter1);
                    angleNew = modelHandle(1, 1).model.getAngle(command) * 2;
                    iter1 = iter1 + 1;
                end
                trajectory(odIndex, vergErrIndex, stimIter, 1, :) = command;

                for iter = 1 : numIters
                    modelHandle(1, 1).model.refreshImagesNew(simulator, currentTexture, angleNew / 2, objDist(odIndex), 3);

                    for i = 1 : length(modelHandle(1, 1).model.scModel)
                        modelHandle(1, 1).model.preprocessImage(i, 1);
                        modelHandle(1, 1).model.preprocessImage(i, 2);
                        currentView{i} = vertcat(modelHandle(1, 1).model.patchesLeft{i}, modelHandle(1, 1).model.patchesRight{i});
                    end

                    [bfFeature, ~, ~] = modelHandle(1, 1).model.generateFR(currentView); % encode image patches

                    if (modelHandle(1, 1).model.normFeatVect == 0)
                        %% Standard feature vector compilation:
                        % append muscle activities to basis function vector
                        feature = [bfFeature; command * modelHandle(1, 1).model.lambdaMuscleFB];
                    else
                        %% Normalized feature vector:
                        % z-transform raw feature vector (no muscle feedback scaling)
                        feature = [bfFeature; command];
                        for i = 1 : length(feature)
                            feature(i) = modelHandle(1, 1).model.onlineNormalize(modelHandle(1, 1).model.trainedUntil, feature(i), i, 0);
                        end
                        feature = [feature(1 : end - 2); feature(end - 1 : end) * modelHandle(1, 1).model.lambdaMuscleFB];
                    end

                    %% concatinate bias entry
                    if (modelHandle(1, 1).model.rlModel.bias > 0)
                        feature = [feature; modelHandle(1, 1).model.rlModel.bias];
                    end

                    relativeCommand = modelHandle(1, 1).model.rlModel.act(feature);    % generate change in muscle activity
                    command = checkCmd(command + relativeCommand);                  % calculate new muscle activities
                    angleNew = modelHandle(1, 1).model.getAngle(command) * 2;          % transform into angle

                    trajectory(odIndex, vergErrIndex, stimIter, iter + 1, :) = command;
                end
            end
        end
    end

    for odIndex = 1 : length(objDist)
        angleDes = 2 * atand(modelHandle(1, 2).model.baseline / (2 * objDist(odIndex)));

        for stimIter = 1 : nStimuli
            currentTexture = stimuliIndices(stimIter);

            for vergErrIndex = 1 : length(startVergErr)
                % muscle init
                if (initMethod == 0)
                    % catch negative/divergent vergence angle
                    % [command, angleNew] = modelHandle(1, 2).model.getMFedood(objDist(odIndex), min(startVergErr(vergErrIndex), modelHandle(1, 2).model.getVergErrMax(objDist(odIndex))));
                    [command, angleNew] = modelHandle(1, 2).model.getMFedoodD(objDist(odIndex), ...
                                                                           min(startVergErr(vergErrIndex), modelHandle(1, 2).model.getVergErrMax(objDist(odIndex))), ...
                                                                           mrVal(2, odIndex * 2 - 1 + vergErrIndex - 1));
                elseif (initMethod == 1)
                    command = cmdInit(:, iter1);
                    angleNew = modelHandle(1, 2).model.getAngle(command) * 2;
                    iter1 = iter1 + 1;
                end
                trajectory2(odIndex, vergErrIndex, stimIter, 1, :) = command;

                for iter = 1 : numIters
                    modelHandle(1, 2).model.refreshImagesNew(simulator, currentTexture, angleNew / 2, objDist(odIndex), 3);

                    for i = 1 : length(modelHandle(1, 2).model.scModel)
                        modelHandle(1, 2).model.preprocessImage(i, 1);
                        modelHandle(1, 2).model.preprocessImage(i, 2);
                        currentView{i} = vertcat(modelHandle(1, 2).model.patchesLeft{i}, modelHandle(1, 2).model.patchesRight{i});
                    end

                    [bfFeature, ~, ~] = modelHandle(1, 2).model.generateFR(currentView); % encode image patches

                    if (modelHandle(1, 2).model.normFeatVect == 0)
                        %% Standard feature vector compilation:
                        % append muscle activities to basis function vector
                        feature = [bfFeature; command * modelHandle(1, 2).model.lambdaMuscleFB];
                    else
                        %% Normalized feature vector:
                        % z-transform raw feature vector (no muscle feedback scaling)
                        feature = [bfFeature; command];
                        for i = 1 : length(feature)
                            feature(i) = modelHandle(1, 2).model.onlineNormalize(modelHandle(1, 2).model.trainedUntil, feature(i), i, 0);
                        end
                        feature = [feature(1 : end - 2); feature(end - 1 : end) * modelHandle(1, 2).model.lambdaMuscleFB];
                    end

                    %% bias analysis
                    if (modelHandle(1, 2).model.rlModel.bias > 0)
                        feature = [feature; modelHandle(1, 2).model.rlModel.bias];
                    end

                    relativeCommand = modelHandle(1, 2).model.rlModel.act(feature);    % generate change in muscle activity
                    command = checkCmd(command + relativeCommand);                  % calculate new muscle activities
                    angleNew = modelHandle(1, 2).model.getAngle(command) * 2;          % transform into angle

                    trajectory2(odIndex, vergErrIndex, stimIter, iter + 1, :) = command;
                end
            end
        end
    end

    %% Plot results
    % colors{1} = 'w/o met. costs'
    % colors{2} = 'w/  met. costs'
    % colors{i} = {traj_line(1:10)_line, traj_line(1:10)_MarkerEdgeColor, traj_line(1:10)_MarkerFaceColor,
    %              traj_line(11:20),  traj_line(1:10)_MarkerEdgeColor, traj_line(1:10)_MarkerFaceColor
    %              end_fixation_MarkerEdgeColor, end_fixation_MarkerFaceColor}
    %
    % orange = [1, 94 / 255, 41 / 255],
    % lightOrange = [1, 201 / 255, 41 / 255]
    colors = {{colors{1}, colors{1}, colors{1}, [1, 201 / 255, 41 / 255], 'k', 'y'}, ...
              {colors{2}, colors{2}, colors{2}, [1, 94 / 255, 41 / 255], 'k', 'm'}};
              % {[0, 95/255, 0], [0, 95/255, 0], [0, 95/255, 0], [1, 94 / 255, 41 / 255], 'k', 'm'}};

    lineStyles = {'-', '-'};
    % markers = {'o', 'o'}; % {traj_line(1:10), end_fixation}
    markers = {'*', 'o'};

    figC = figure('OuterPosition', [100, 100, 600, 600]);
    hold on;

    title('Object Fixation Trajectories');
    xlabel('lateral rectus activation [%]');
    ylabel('medial rectus activation [%]');

    % pcHandle = pcolor(modelHandle(1, 1).model.degreesIncRes); % use vergence degree as color dimension (background)
    pcHandle = pcolor(modelHandle(1, 1).model.metCostsIncRes);  % use metabolic costs as color dimension (background)
    % shading interp;
    set(pcHandle, 'EdgeColor', 'none');

    colormap(createCM(7));
    cb = colorbar();
    % cb.Label.String = 'vergence degree'; % use vergence degree as color dimension (background)
    cb.Label.String = 'metabolic costs [J]';   % use metabolic costs as color dimension (background)

    ax = gca;
    set(ax, 'Layer','top'); % bring axis to the front

    ax.XTick = linspace(1, size(modelHandle(1, 1).model.degreesIncRes, 2), 11);
    ax.YTick = linspace(1, size(modelHandle(1, 1).model.degreesIncRes, 1), 11);

    ax.XTickLabel = strsplit(num2str(linspace(0, 10, 11)));
    ax.YTickLabel = strsplit(num2str(linspace(0, 20, 11)));

    axis([1, size(modelHandle(1, 1).model.degreesIncRes, 2), 1, size(modelHandle(1, 1).model.degreesIncRes, 1)]);
    axPos = ax.Position;

    ht = cell(1, length(objDist)); % objDist label handles

    % draw objects + offsets
    for odIndex = 1 : length(objDist)
        % draw +1 pixel offset in respect to desired vergence distance
        [lateralDes, medialDes] = modelHandle(1, 1).model.getAnglePoints(objDist(odIndex), 0.22);
        plot(lateralDes ./ modelHandle(1, 1).model.scaleFacLR, medialDes ./ modelHandle(1, 1).model.scaleFacMR, ...
             'color', 'k', 'LineStyle', ':', 'LineWidth', 1);%[0, 0.5882, 0]
             %'color', [0, 0.5882, 0], 'LineStyle', ':', 'LineWidth', 1);%[0, 0.5882, 0]

        % draw -1 pixel offset in respect to desired vergence distance
        [lateralDes, medialDes] = modelHandle(1, 1).model.getAnglePoints(objDist(odIndex), -0.22);
        plot(lateralDes ./ modelHandle(1, 1).model.scaleFacLR, medialDes ./ modelHandle(1, 1).model.scaleFacMR, ...
             'color', 'k', 'LineStyle', ':', 'LineWidth', 1);%[0, 0.5882, 0]
             %'color', [0, 0.5882, 0], 'LineStyle', ':', 'LineWidth', 1);%[0, 0.5882, 0]

        % draw a line of points into the plane that represent the desired vergence
        [lateralDes, medialDes] = modelHandle(1, 1).model.getAnglePoints(objDist(odIndex), 0);
        plot(lateralDes ./ modelHandle(1, 1).model.scaleFacLR, medialDes ./ modelHandle(1, 1).model.scaleFacMR, ...
                            'color', 'k', 'LineWidth', 0.75);
    % 'color', [0.6510, 1.0000, 0.6588], 'LineWidth', 1.8);


        % add corresponding distance value to desired vergence graph
        ht{odIndex} = text(lateralDes(end - ceil(length(lateralDes) / 10)) / modelHandle(1, 1).model.scaleFacLR, ...
                           medialDes(end - ceil(length(medialDes) / 10)) / modelHandle(1, 1).model.scaleFacMR, ...
                           sprintf('%3.1f m', objDist(odIndex)));
    end

    % draw trajectories
    for odIndex = 1 : length(objDist)
        for stim = 1 : length(stimuliIndices)
            for vergErrIndex = 1 : length(startVergErr)
                % first plot whole trajectory without metabolic costs
                hl2 = plot(reshape(trajectory(odIndex, vergErrIndex, stim, :, 1), [numIters + 1, 1]) ./ modelHandle(1, 1).model.scaleFacLR + 1, ...
                           reshape(trajectory(odIndex, vergErrIndex, stim, :, 2), [numIters + 1, 1]) ./ modelHandle(1, 1).model.scaleFacMR + 1, ...
                           'Color', colors{1}{1}, 'LineStyle', lineStyles{1}, 'LineWidth', 2.1, ...%1
                           'Marker', markers{1}, 'MarkerEdgeColor', colors{1}{2}, 'MarkerFaceColor',  colors{1}{3}, 'MarkerSize', 3);%4

                % % plot iter 1-interval in differen color if numIters >= model.interval
                % if (numIters >= modelHandle(1, 1).model.interval)
                %     hl1 = plot(reshape(trajectory(odIndex, vergErrIndex, stim, 1 : (modelHandle(1, 1).model.interval), 1), [modelHandle(1, 1).model.interval, 1]) ./ modelHandle(1, 1).model.scaleFacLR + 1, ...
                %                reshape(trajectory(odIndex, vergErrIndex, stim, 1 : (modelHandle(1, 1).model.interval), 2), [modelHandle(1, 1).model.interval, 1]) ./ modelHandle(1, 1).model.scaleFacMR + 1, ...
                %                'Color', colors{1}{1}, 'LineStyle', lineStyles{1}, 'LineWidth', 1, 'Marker', markers{1}, 'MarkerEdgeColor', colors{1}{2}, 'MarkerFaceColor',  colors{1}{3}, 'MarkerSize', 4);
                % else
                %     hl1 = [];
                % end

                % % plot init point
                % plot(trajectory(odIndex, vergErrIndex, stim, 1, 1) / modelHandle(1, 1).model.scaleFacLR + 1, ...
                %      trajectory(odIndex, vergErrIndex, stim, 1, 2) / modelHandle(1, 1).model.scaleFacMR + 1, ...
                %     'MarkerEdgeColor','k', 'MarkerFaceColor', 'r', 'MarkerSize', 4);

                % plot destination point
                hl3 = plot(trajectory(odIndex, vergErrIndex, stim, end, 1) / modelHandle(1, 1).model.scaleFacLR + 1, ...
                           trajectory(odIndex, vergErrIndex, stim, end, 2) / modelHandle(1, 1).model.scaleFacMR + 1, ...
                           'LineStyle', 'none', 'Marker', markers{2}, 'MarkerEdgeColor', colors{1}{5}, 'MarkerFaceColor',  colors{1}{6}, 'MarkerSize', 4);%'LineWidth', 1, 'Marker', markers{2}, 'MarkerEdgeColor', colors{1}{5}, 'MarkerFaceColor',  colors{1}{6}, 'MarkerSize', 4);
            end
        end
    end

    for odIndex = 1 : length(objDist)
        for stim = 1 : length(stimuliIndices)
            for vergErrIndex = 1 : length(startVergErr)
                % first plot whole trajectory with metabolic costs
                hl5 = plot(reshape(trajectory2(odIndex, vergErrIndex, stim, :, 1), [numIters + 1, 1]) ./ modelHandle(1, 2).model.scaleFacLR + 1, ...
                           reshape(trajectory2(odIndex, vergErrIndex, stim, :, 2), [numIters + 1, 1]) ./ modelHandle(1, 2).model.scaleFacMR + 1, ...
                           'Color', colors{2}{1}, 'LineStyle', lineStyles{2}, 'LineWidth', 2.1, ...%1
                           'Marker', markers{1}, 'MarkerEdgeColor', colors{2}{2}, 'MarkerFaceColor',  colors{2}{3}, 'MarkerSize', 3);%4

                % % plot iter 1-interval in differen color if numIters >= model.interval
                % if (numIters >= modelHandle(1, 2).model.interval)
                %     hl4 = plot(reshape(trajectory2(odIndex, vergErrIndex, stim, 1 : (modelHandle(1, 2).model.interval), 1), [modelHandle(1, 2).model.interval, 1]) ./ modelHandle(1, 2).model.scaleFacLR + 1, ...
                %                reshape(trajectory2(odIndex, vergErrIndex, stim, 1 : (modelHandle(1, 2).model.interval), 2), [modelHandle(1, 2).model.interval, 1]) ./ modelHandle(1, 2).model.scaleFacMR + 1, ...
                %                'Color', colors{2}{1}, 'LineStyle', lineStyles{2}, 'LineWidth', 1, 'Marker', markers{1}, 'MarkerEdgeColor', colors{2}{2}, 'MarkerFaceColor',  colors{2}{3}, 'MarkerSize', 4);
                % else
                %     hl4 = [];
                % end

                % % plot init point
                % plot(trajectory2(odIndex, vergErrIndex, stim, 1, 1) / modelHandle(1, 2).model.scaleFacLR + 1, ...
                %      trajectory2(odIndex, vergErrIndex, stim, 1, 2) / modelHandle(1, 2).model.scaleFacMR + 1, ...
                %     'MarkerEdgeColor','k', 'MarkerFaceColor', 'r', 'MarkerSize', 4);

                % plot destination point
                hl6 = plot(trajectory2(odIndex, vergErrIndex, stim, end, 1) / modelHandle(1, 2).model.scaleFacLR + 1, ...
                           trajectory2(odIndex, vergErrIndex, stim, end, 2) / modelHandle(1, 2).model.scaleFacMR + 1, ...
                           'LineStyle', 'none', 'Marker', markers{2}, 'MarkerEdgeColor', colors{2}{5}, 'MarkerFaceColor',  colors{2}{6}, 'MarkerSize', 4);%'LineWidth', 1, 'Marker', markers{2}, 'MarkerEdgeColor', colors{2}{5}, 'MarkerFaceColor',  colors{2}{6}, 'MarkerSize', 4);
            end
        end
    end

    % gKey = {sprintf('1..%dth  iteration', modelHandle(1, 1).model.interval), ...
    %         sprintf('%d..%dth iteration', modelHandle(1, 1).model.interval, numIters), ...
    %         'end fixation w/o met. costs', ...
    %         sprintf('1..%dth  iteration', modelHandle(1, 1).model.interval), ...
    %         sprintf('%d..%dth iteration', modelHandle(1, 1).model.interval, numIters), ...
    %         'end fixation w/  met. costs'};

    gKey = {sprintf('0th..%dth iteration', numIters), ...
            'end fixation  \bfw/o met. costs', ...
            sprintf('0th..%dth iteration', numIters), ...
            'end fixation  \bfw/   met. costs'};

    % hDummy1 = plot(NaN, NaN, 'LineStyle', 'none');
    % hDummy2 = plot(NaN, NaN, 'LineStyle', 'none');

    % with dummy entries for last 2 columns
    % gKey = {sprintf('0th..%dth iteration', numIters), ...
    %         'end fixation', ...
    %         '\bfw/o met. costs', ...
    %         sprintf('0th..%dth iteration', numIters), ...
    %         'end fixation', ...
    %         '\bfw/   met. costs'};

    l = gridLegend([hl2, hl3, hl5, hl6], 2, gKey, 'Orientation', 'Horizontal', 'Location', 'southoutside', 'Fontsize', 8);
    % l = gridLegend([hl2, hl3, hDummy1, hl5, hl6, hDummy2], 3, gKey, 'Orientation', 'Horizontal', 'Location', 'southoutside', 'Fontsize', 8);
    % l.Box = 'off';

    %% repositioning subfigures
    % ax.Position = axPos;
    % ax.PlotBoxAspectRatioMode = 'manual';
    % ax.DataAspectRatioMode = 'manual';
    % ax.ActivePositionProperty = 'position';
    % ax.OuterPosition(4) = ax.OuterPosition(4) * 1.1;
    ax.DataAspectRatioMode = 'manual';

    tmp = ht{1}.Position(2) - ht{1}.Position(2) * 0.96;
    ht{1}.Position(2) = ht{1}.Position(2) * 0.96;
    ht{2}.Position(2) = ht{2}.Position(2) - tmp;

    % manual positioning
    ax.Position = [0.085, 0.2, 0.75, 0.75];
    l.Position(1) = 0.15;
    l.Position(2) = 0.04;

    % set(figC,'PaperPositionMode','auto'); % keep aspect ratio
    plotpath = sprintf('%s/FigC_muscleActivityTrajectories', savePath);
    saveas(figC, plotpath, 'png');
    close(figC);
end

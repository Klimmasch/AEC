% Generates results plots for ICDL conference publication 2017
function generateICDLPlotsAvg(simulator)

    % dataVergErr := {0 = median(vergerr), 1 = median(abs(vergerr))}
    dataVergErr = 1;

    %% TODO: implement flexible modelAt
    modelAt1 = [250000 : 250000 : 1000000];
    modelAt2 = [200000 : 200000 : 1000000];
    % modelAt1 = [50000 : 50000 : 200000, 300000 : 100000 : 1000000];
    % modelAt2 = modelAt1;

    % check modelAt
    checked = false;
%     while (~checked)
%         if (exist(sprintf('%s/%s', modelWoMetCostsFullPath, sprintf('modelAt%d', modelAt1(1))), 'dir') == 7)
%             if (exist(sprintf('%s/%s', modelWMetCostsFullPath, sprintf('modelAt%d', modelAt2(1))), 'dir') == 7)
%                 checked = true;
%             elseif (modelAt2(1) == 250000)
%                 error('modelWMetCosts could not be checked.');
%             else
%                 modelAt2 = [250000 : 250000 : 1000000];
%             end
%         elseif (modelAt1(1) == 250000)
%             error('modelWoMetCosts could not be checked.');
%         else
%             modelAt1 = [250000 : 250000 : 1000000];
%         end
%     end
    testPoints = {modelAt1, modelAt2};

    % load given models
    try
        % modelHandle = [load(strcat('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_8_gamma_0.3_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_8_gamma_0.3_metCost_[0.035]', '/model.mat')); ...
        %                load(strcat('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_9_gamma_0.3_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_9_gamma_0.3_metCost_[0.035]', '/model.mat')); ...
        %                load(strcat('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_10_gamma_0.3_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_10_gamma_0.3_metCost_[0.035]', '/model.mat')); ...
        %                load(strcat('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_11_gamma_0.3_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_11_gamma_0.3_metCost_[0.035]', '/model.mat')); ...
        %                load(strcat('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_12_gamma_0.3_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_12_gamma_0.3_metCost_[0.035]', '/model.mat'))];
        modelHandle = [load(strcat('/home/aecgroup/aecdata/Results/GammaVsMetCosts_FinerGrainBias002/17-03-09_1000000iter_1_gamma_0.3_metCost_[0.000]', '/model.mat')), load(strcat('/home/aecgroup/aecdata/Results/GammaVsMetCosts_FinerGrainBias002/17-03-17_1000000iter_3_gamma_0.3_metCost_[0.035]', '/model.mat'))];
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

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% Figure A
    % % RMSE vergence error [deg] & delta MC opt [%] @ testing vs. traintime
    % objHandels = {};
    % lineStyles = {'-', '--'};                       % [left y-axis, right y-axis]
    % lineWidths = [2, 2];                            % [vergErr, metCosts]
    fontSizes = [12];                               % global font size for axis labels, axis ticks and legends

    % markerStyles = ['x', 'o'];                      % [left y-axis, right y-axis]
    % markerSizes = [5, 5];                           % [vergErr, metCosts]

    % % median vergErr    [w/o metcosts, w/ metcosts]
    % % median metCosts   [w/o metcosts, w/ metcosts]
    % % IQR_patches       [w/o metCosts, w/ metCosts]
    % colors = {[0, 100/255, 200/255], [0, 95/255, 0], ...
    %           [0, 100/255, 200/255], [0, 95/255, 0], ...
    %           [0, 100/255, 200/255], [0, 1, 128/255]};

    % alphas = {0.2, 0.3};

    % nTicks = 10;

    % objRange = [0.5, 1 : 6];
    % tmpModelHandle = cell(2, max(length(modelAt1), length(modelAt2)));

    % % dataMatrix = {w/o[vergErrMatrix, metCostsMatrix], w/[vergErrMatrix, metCostsMatrix]}
    % % vergErrMatrix = metCostsMatrix = objDist * VSE/0° * nStim * testInterval x testAt
    % dataMatrix = {{zeros(1680 * size(modelHandle, 1), length(modelAt1)), []}, {zeros(1680 * size(modelHandle, 1), length(modelAt2)), []}};
    % dataMatrixEnd = {{zeros(1680 * size(modelHandle, 1), 20), []}, {zeros(1680 * size(modelHandle, 1), 20), []}}; % final values @ modelAt = 1mio & iter = [1, 20]

    % at0Matrix = {{zeros(1680, 1), []}, {zeros(1680, 1), []}};
    % iqrLine = zeros(4, max(length(modelAt1), length(modelAt2)) + 1); % +modelAt0

    % % extract all relevant data from all sub-experiments
    % for i = 1 : size(modelHandle, 2)
    %     for seedIter = 1 : size(modelHandle, 1)
    %         for trainedUntil = 1 : length(testPoints{i})
    %             try
    %                 subFolder = sprintf('modelAt%d', testPoints{i}(trainedUntil));
    %                 tmpModelHandle{i, trainedUntil} = load(sprintf('%s/%s/model.mat', modelHandle(seedIter, i).model.savePath, subFolder));
    %             catch
    %                 % catch case when (sub-)experiment started, but has no test results yet
    %                 error('%s/%s/model.mat\ncould not be loaded.', modelHandle(seedIter, i).model.savePath, subFolder);
    %             end

    %             % fill data matrix & exclude VSA = 0° trials
    %             nStim = size(tmpModelHandle{i, trainedUntil}.model.testResult3, 1) / (length(objRange) * 7);
    %             currIdx = 1 + 1680 * (seedIter - 1);
    %             startInd = 1;
    %             startInd0 = nStim * 3 + 1;
    %             for j = 1 : length(objRange)
    %                 colSize1 = length(tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : startInd0 - 1, tmpModelHandle{i, trainedUntil}.model.testInterval));
    %                 dataMatrix{i}{1}(currIdx : currIdx + colSize1 - 1, trainedUntil) = tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : startInd0 - 1, ...
    %                                                                                                   tmpModelHandle{i, trainedUntil}.model.testInterval);

    %                 % end matrix creation
    %                 if (trainedUntil == length(testPoints{i}))
    %                     for k = 1 : 20
    %                         dataMatrixEnd{i}{1}(currIdx : currIdx + colSize1 - 1, k) = tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : startInd0 - 1, k);
    %                     end
    %                 end

    %                 currIdx = currIdx + colSize1;
    %                 endInd0 = startInd0 + nStim - 1;
    %                 startInd = endInd0 + 1;
    %                 startInd0 = endInd0 + nStim * 6 + 1;
    %             end
    %             % concatinate remainder
    %             colSize2 = length(tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : end, tmpModelHandle{i, trainedUntil}.model.testInterval));
    %             dataMatrix{i}{1}(currIdx : currIdx + colSize2 - 1, trainedUntil) = tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : end, ...
    %                                                                                               tmpModelHandle{i, trainedUntil}.model.testInterval);

    %             % delta metCosts
    %             dataMatrix{i}{2}(1 + 1960 * (seedIter - 1) : 1960 * seedIter, trainedUntil) = tmpModelHandle{i, trainedUntil}.model.testResult7(:, tmpModelHandle{i, trainedUntil}.model.testInterval);

    %             % end matrix creation
    %             if (trainedUntil == length(testPoints{i}))
    %                 for k = 1 : 20
    %                     dataMatrixEnd{i}{1}(currIdx : currIdx + colSize2 - 1, k) = tmpModelHandle{i, trainedUntil}.model.testResult3(startInd : end, k);
    %                 end
    %                 dataMatrixEnd{i}{2}(1 + 1960 * (seedIter - 1) : 1960 * seedIter, :) = tmpModelHandle{i, trainedUntil}.model.testResult7;
    %             end
    %         end
    %     end

    %     % add modelAt0 entries
    %     if (i == 1) && (length(testPoints{i}) == 4)
    %         hm = load('/home/aecgroup/aecdata/Results/17-03-08_300000iter_1_newStandard_0,3Mio/modelAt0/model.mat');
    %     elseif (length(testPoints{i}) == 12)
    %         if (i == 1)
    %             hm = load('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_8_gamma_0.3_metCost_[0.000]/modelAt0/model.mat');
    %         else
    %             hm = load('/home/aecgroup/aecdata/Results/BiasedBmsfX5NewTest/17-03-23_1000000iter_8_gamma_0.3_metCost_[0.035]/modelAt0/model.mat');
    %         end
    %     end

    %     % fill data matrix & exclude VSA = 0° trials
    %     nStim = size(hm.model.testResult3, 1) / (length(objRange) * 7);
    %     currIdx = 1;
    %     startInd = 1;
    %     startInd0 = nStim * 3 + 1;
    %     for j = 1 : length(objRange)
    %         colSize1 = length(hm.model.testResult3(startInd : startInd0 - 1, hm.model.testInterval));
    %         at0Matrix{i}{1}(currIdx : currIdx + colSize1 - 1) = hm.model.testResult3(startInd : startInd0 - 1, hm.model.testInterval);
    %         currIdx = currIdx + colSize1;
    %         endInd0 = startInd0 + nStim - 1;
    %         startInd = endInd0 + 1;
    %         startInd0 = endInd0 + nStim * 6 + 1;
    %     end
    %     % concatinate remainder
    %     colSize2 = length(hm.model.testResult3(startInd : end, hm.model.testInterval));
    %     at0Matrix{i}{1}(currIdx : currIdx + colSize2 - 1) = hm.model.testResult3(startInd : end, hm.model.testInterval);

    %     % delta metCosts
    %     at0Matrix{i}{2} = hm.model.testResult7(:, hm.model.testInterval);

    %     % concatinate both matricies
    %     dataMatrix{i}{1} = horzcat(repmat(at0Matrix{i}{1}, [size(modelHandle, 1), 1]), dataMatrix{i}{1});
    %     dataMatrix{i}{2} = horzcat(repmat(at0Matrix{i}{2}, [size(modelHandle, 1), 1]), dataMatrix{i}{2});

    %     if (dataVergErr == 1)
    %         dataMatrix{i}{1} = abs(dataMatrix{i}{1});
    %     end

    %     % extract IQR edge coordinates
    %     tmpFig = figure();
    %     boxHandle = boxplot(dataMatrix{i}{1});
    %     upWi = findobj(boxHandle, 'tag', 'Upper Whisker');
    %     lowWi = findobj(boxHandle, 'tag', 'Lower Whisker');
    %     iqrLine(i * 2 - 1 : i * 2, :) = [arrayfun(@(x) x.YData(1), upWi)'; arrayfun(@(x) x.YData(2), lowWi)'];
    %     close(tmpFig);
    % end

    % % plot
    % figA = figure();
    % hold on;
    % modelAt1 = horzcat(0, modelAt1);
    % modelAt2 = horzcat(0, modelAt2);
    % testPoints = {modelAt1, modelAt2};

    % for i = 1 : size(modelHandle, 2)
    %     if (i == 1)
    %         ax1 = gca; % current axes
    %         % ax1.YColor = colors{1};
    %         ax1.XAxis.Label.String = 'Traintime';
    %         ax1.XAxis.Label.FontSize = fontSizes(1);
    %         % ax1.Title.String = 'Test Performance & Metabolic Costs vs. Traintime';
    %         if (dataVergErr == 0)
    %             % ax1.YAxis.Label.String = 'median(\Delta\gamma) [deg]';
    %             ax1.YAxis.Label.String = '|\Delta\gamma| [deg]';
    %         elseif (dataVergErr == 1)
    %             % ax1.YAxis.Label.String = 'median(|\Delta\gamma|) [deg]';
    %             ax1.YAxis.Label.String = '|\Delta\gamma| [deg]';
    %         else
    %             error('dataVergErr = %d is not supported.', dataVergErr);
    %         end
    %         ax1.YAxis.Label.FontSize = fontSizes(1);

    %         ax2 = axes('Position', ax1.Position, ...
    %                    'YAxisLocation', 'right', ...
    %                    'Color', 'none');
    %     end

    %     % fill area defined by upper & lower IQR bounds
    %     objHandels{end + 1} = patch([testPoints{i}, flip(testPoints{i})], ...
    %                                 [iqrLine(i * 2 - 1, 1 : length(testPoints{i})), flip(iqrLine(i * 2, 1 : length(testPoints{i})))], ...
    %                                 colors{i + 4}, 'LineStyle', 'none', 'FaceAlpha', alphas{i});
    %     objHandels{end}.Parent = ax1;
    %     hold on;

    %     % median(|verg_{err}|) [deg]
    %     objHandels{end + 1} = plot(ax1, testPoints{i}, median(dataMatrix{i}{1}), ...
    %                                'LineStyle', lineStyles{1}, ...
    %                                'Marker', markerStyles(1), 'MarkerSize', markerSizes(1), ...
    %                                'Color', colors{1 + i - 1}, 'LineWidth', lineWidths(1));
    %     hold on;

    %     % median(\Delta C^{met}) [W]
    %     objHandels{end + 1} = plot(ax2, testPoints{i}, median(dataMatrix{i}{2}), ...
    %                                'LineStyle', lineStyles{2}, ...
    %                                'Marker', markerStyles(2), 'MarkerSize', markerSizes(2), ...
    %                                'MarkerEdgeColor', colors{i}, 'MarkerFaceColor', colors{i}, ...
    %                                'Color', colors{3 + i - 1}, 'LineWidth', lineWidths(2));
    %     hold on;

    %     if (i == 1)
    %         % |\DeltaMC_{opt}| = |MC_{actual} - MC_{optimal}| / |MC_{start} - MC_{optimal}|
    %         ax2.YAxis.Label.String = '\DeltaC [W]';
    %         ax2.YAxis.Label.FontSize = fontSizes(1);
    %         % ax2.YColor = colors{3};
    %         % ax2.YAxis.Label.Rotation = -90;
    %         ax2.YAxisLocation = 'right';
    %     end
    % end

    % % grid(ax1, 'on');
    % % ax1.YMinorGrid = 'on';

    % ax1.YAxis.Limits = [-0.1, 2.5];
    % ax2.YAxis.Limits = [-0.1, 1.2];

    % % set #nTicks ticks for y-axis
    % % set(ax1, 'YTick', round(linspace(ax1.YAxis.Limits(1), ax1.YAxis.Limits(2), nTicks - 2), 2));
    % % set(ax2, 'YTick', round(linspace(ax2.YAxis.Limits(1), ax2.YAxis.Limits(2), nTicks), 2));

    % %% align both y-axis to zero
    % ratio = ax2.YAxis.Limits(1) / ax2.YAxis.Limits(2);
    % if ax1.YAxis.Limits(2) * ratio < ax1.YAxis.Limits(1)
    %     ax1.YAxis.Limits = [ax1.YAxis.Limits(2) * ratio, ax1.YAxis.Limits(2)];
    % else
    %     ax1.YAxis.Limits = [ax1.YAxis.Limits(1), ax1.YAxis.Limits(1) / ratio];
    % end

    % % realign plot order
    % uistack(objHandels{1},'bottom');
    % uistack(objHandels{4},'bottom');
    % uistack(objHandels{2},'top');
    % uistack(objHandels{5},'top');

    % % lineHandles = [objHandels{2}, objHandels{1}, objHandels{3}, objHandels{5}, objHandels{4}, objHandels{6}];
    % % gKey = {'M_{ } median(|\Delta\gamma|)', 'M_{ } IQR', 'M_{ } median(\DeltaC)', ...
    % %         'M_{C} median(|\Delta\gamma|)', 'M_{C} IQR', 'M_{C} median(\DeltaC)'};
    % lineHandles = [objHandels{2}, objHandels{3}, objHandels{5}, objHandels{6}];
    % gKey = {'M_{ } |\Delta\gamma|', 'M_{ } \DeltaC', ...
    %         'M_{C} |\Delta\gamma|', 'M_{C} \DeltaC'};


    % l = legend(lineHandles, gKey, 'Location', 'east');
    % l.FontSize = fontSizes(1);
    % l.Position(2) = 0.4;
    % l.Box = 'off';

    % % ax2.YAxis.Label.Position(1) = ax2.YAxis.Label.Position(1) * 1.5;

    % plotpath = sprintf('%s/FigA_VergErrMetCostsVsTraintime', savePath);
    % saveas(figA, plotpath, 'png');
    % close(figA);

    % % member models? III member!
    % fileID = fopen(strcat(savePath, '/README.txt'), 'at' );
    % fprintf(fileID, 'model w/o MetCosts: %s\n', modelHandle(1, 1).model.savePath);
    % fprintf(fileID, 'model w/  MetCosts: %s\n\n', modelHandle(1, 2).model.savePath);

    % fprintf(fileID, '========================================================================================================\n');
    % fprintf(fileID, 'model w/o MetCosts\n');
    % fprintf(fileID, 'modelAt: %s\n\n', int2str(modelAt1));

    % if (length(modelAt1) == 6)
    %     fprintf(fileID, 'figA vergErr median:\t%f %f %f %f %f %f\n', median(dataMatrix{1}{1}));
    %     fprintf(fileID, 'figA vergErr iqr:\t%f %f %f %f %f %f\n', iqr(dataMatrix{1}{1}));
    %     fprintf(fileID, 'figA vergErr mean:\t%f %f %f %f %f %f\n', mean(dataMatrix{1}{1}));
    %     fprintf(fileID, 'figA vergErr std:\t%f %f %f %f %f %f\n\n', std(dataMatrix{1}{1}));

    %     fprintf(fileID, 'figA metCosts median:\t%f %f %f %f %f %f\n', median(dataMatrix{1}{2}));
    %     fprintf(fileID, 'figA metCosts iqr:\t%f %f %f %f %f %f\n', iqr(dataMatrix{1}{2}));
    %     fprintf(fileID, 'figA metCosts mean:\t%f %f %f %f %f %f\n', mean(dataMatrix{1}{2}));
    %     fprintf(fileID, 'figA metCosts std:\t%f %f %f %f %f %f\n\n', std(dataMatrix{1}{2}));
    % else
    %     fprintf(fileID, 'figA vergErr median:\t%f %f %f %f %f\n', median(dataMatrix{1}{1}));
    %     fprintf(fileID, 'figA vergErr iqr:\t%f %f %f %f %f\n', iqr(dataMatrix{1}{1}));
    %     fprintf(fileID, 'figA vergErr mean:\t%f %f %f %f %f\n', mean(dataMatrix{1}{1}));
    %     fprintf(fileID, 'figA vergErr std:\t%f %f %f %f %f\n\n', std(dataMatrix{1}{1}));

    %     fprintf(fileID, 'figA metCosts median:\t%f %f %f %f %f\n', median(dataMatrix{1}{2}));
    %     fprintf(fileID, 'figA metCosts iqr:\t%f %f %f %f %f\n', iqr(dataMatrix{1}{2}));
    %     fprintf(fileID, 'figA metCosts mean:\t%f %f %f %f %f\n', mean(dataMatrix{1}{2}));
    %     fprintf(fileID, 'figA metCosts std:\t%f %f %f %f %f\n\n', std(dataMatrix{1}{2}));
    % end

    % fprintf(fileID, 'model w/ MetCosts\n');
    % fprintf(fileID, 'modelAt: %s\n\n', int2str(modelAt2));

    % if (length(modelAt2) == 6)
    %     fprintf(fileID, 'figA vergErr median:\t%f %f %f %f %f %f\n', median(dataMatrix{2}{1}));
    %     fprintf(fileID, 'figA vergErr iqr:\t%f %f %f %f %f %f\n', iqr(dataMatrix{2}{1}));
    %     fprintf(fileID, 'figA vergErr mean:\t%f %f %f %f %f %f\n', mean(dataMatrix{2}{1}));
    %     fprintf(fileID, 'figA vergErr std:\t%f %f %f %f %f %f\n\n', std(dataMatrix{2}{1}));

    %     fprintf(fileID, 'figA metCosts median:\t%f %f %f %f %f %f\n', median(dataMatrix{2}{2}));
    %     fprintf(fileID, 'figA metCosts iqr:\t%f %f %f %f %f %f\n', iqr(dataMatrix{2}{2}));
    %     fprintf(fileID, 'figA metCosts mean:\t%f %f %f %f %f %f\n', mean(dataMatrix{2}{2}));
    %     fprintf(fileID, 'figA metCosts std:\t%f %f %f %f %f %f\n\n', std(dataMatrix{2}{2}));
    % else
    %     fprintf(fileID, 'figA vergErr median:\t%f %f %f %f %f\n', median(dataMatrix{2}{1}));
    %     fprintf(fileID, 'figA vergErr iqr:\t%f %f %f %f %f\n', iqr(dataMatrix{2}{1}));
    %     fprintf(fileID, 'figA vergErr mean:\t%f %f %f %f %f\n', mean(dataMatrix{2}{1}));
    %     fprintf(fileID, 'figA vergErr std:\t%f %f %f %f %f\n\n', std(dataMatrix{2}{1}));

    %     fprintf(fileID, 'figA metCosts median:\t%f %f %f %f %f\n', median(dataMatrix{2}{2}));
    %     fprintf(fileID, 'figA metCosts iqr:\t%f %f %f %f %f\n', iqr(dataMatrix{2}{2}));
    %     fprintf(fileID, 'figA metCosts mean:\t%f %f %f %f %f\n', mean(dataMatrix{2}{2}));
    %     fprintf(fileID, 'figA metCosts std:\t%f %f %f %f %f\n\n', std(dataMatrix{2}{2}));
    % end

    % fclose(fileID);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% Figure B
    % %% data separated into 2 subplots
    % figB1 = figure();
    % hold on;

    % steps = 3;              % show just first steps iterations & last iteration
    % colors = colors(1 : 2);
    % % colors{1} = colors{1} .* 0.9;
    % % colors{2} = colors{2} .* 0.9;
    % captions = cell(1, 2);
    % captions{1} = 'M_{ }';
    % captions{2} = 'M_{C}';

    % % tmpMatrix = [vergErr_woMetCosts, vergErr_wMetCosts, metCostsApproach_woMetCosts, metCostsApproach_wMetCosts]
    % tmpMatrixVergErr = vertcat(horzcat(modelHandle(1, 1).model.vergenceAngleApproach, modelHandle(1, 2).model.vergenceAngleApproach), ...
    %                           horzcat(modelHandle(2, 1).model.vergenceAngleApproach, modelHandle(2, 2).model.vergenceAngleApproach), ...
    %                           horzcat(modelHandle(3, 1).model.vergenceAngleApproach, modelHandle(3, 2).model.vergenceAngleApproach), ...
    %                           horzcat(modelHandle(4, 1).model.vergenceAngleApproach, modelHandle(4, 2).model.vergenceAngleApproach), ...
    %                           horzcat(modelHandle(5, 1).model.vergenceAngleApproach, modelHandle(5, 2).model.vergenceAngleApproach));

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
    % boxHandl = boxplot(tmpMatrixVergErr, 'positions', pos);
    % tmpHandle = findobj(boxHandl, 'type', 'text');
    % set(tmpHandle, 'Interpreter', 'tex');
    % grid minor;

    % subBoxHandl = findobj(gca,'Tag','Box');
    % % subBoxHandl = findobj(boxHandl,'Tag','Box');

    % boxesArray1 = findobj(boxHandl);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray1(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray1(idx2(j)).Color = colors{2};
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

    % % xlabel('Iteration step', 'FontSize', fontSizes(1));
    % % ylabel(sprintf('Vergence\nError [deg]'), 'FontSize', fontSizes(1));
    % ylabel(sprintf('Vergence Error\nReduction [%%]'), 'FontSize', fontSizes(1));

    % %% put ylabel right and rotate text
    % % set(sub1, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('C^{met}_{reduction} [%%]'), 'rot', -90, 'FontSize', fontSizes(1));
    % % p = get(lh, 'position');
    % % set(sub1,'yaxislocation','left');
    % % set(lh,'position', p);

    % sub2 = subplot(2, 1, 2);
    % boxHandl2 = boxplot(tmpMatrixMetApp, 'positions', pos);
    % grid minor;

    % boxesArray2 = findobj(boxHandl2);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray2(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray2(idx2(j)).Color = colors{2};
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

    % % manually adjust XTicks
    % sub1.XTick = [1.1, 1.6, 2.1, 2.6];
    % sub1.XTickLabel = {'1','2','3', '20'};
    % sub2.XTick = [1.1, 1.6, 2.1, 2.6];
    % sub2.XTickLabel = {'1','2','3', '20'};

    % xlabel('Iteration step', 'FontSize', fontSizes(1));
    % % ylabel(sprintf('Opt. Metabolic\nCosts Approach [%%]'), 'FontSize', fontSizes(1));
    % ylabel(sprintf('c [%%]'), 'FontSize', fontSizes(1));

    % % manual Figure adjustments
    % set(boxesArray1,'LineWidth', 1);
    % set(boxesArray2,'LineWidth', 1);

    % %% put ylabel right and rotate text
    % % set(sub2, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('C^{met}_{reduction} [%%]'), 'rot', -90, 'FontSize', fontSizes(1));
    % % p = get(lh, 'position');
    % % set(sub2, 'yaxislocation', 'left');
    % % set(lh, 'position', p);

    % % suptitle(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));
    % % suptitle(sprintf('Reduction of Vergence Error & Metabolic Costs\nvs. Iteration at Testing'));

    % l = legend(subBoxHandl([2, 1]), captions);
    % l.FontSize = fontSizes(1);
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
    % colors = colors(1 : 2);
    % % colors{1} = colors{1} .* 0.9;
    % % colors{2} = colors{2} .* 0.9;
    % captions = cell(1, 2);
    % captions{1} = 'M_{ }';
    % captions{2} = 'M_{C}';

    % % % tmpMatrix = [vergErr_woMetCosts, vergErr_wMetCosts, metCostsApproach_woMetCosts, metCostsApproach_wMetCosts]
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
    % pos = [1 1.2 1.5 1.7 2 2.2 2.5 2.7];
    % boxHandl = boxplot(tmpMatrixVergErr, 'positions', pos);
    % tmpHandle = findobj(boxHandl, 'type', 'text');
    % set(tmpHandle, 'Interpreter', 'tex');
    % grid minor;

    % subBoxHandl = findobj(gca,'Tag','Box');
    % % subBoxHandl = findobj(boxHandl,'Tag','Box');

    % boxesArray1 = findobj(boxHandl);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray1(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray1(idx2(j)).Color = colors{2};
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

    % ylabel('\Delta\gamma [deg]', 'FontSize', fontSizes(1));

    % %% put ylabel right and rotate text
    % % set(sub1, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('C^{met}_{reduction} [%%]'), 'rot', -90, 'FontSize', fontSizes(1));
    % % p = get(lh, 'position');
    % % set(sub1,'yaxislocation','left');
    % % set(lh,'position', p);

    % sub2 = subplot(2, 1, 2);
    % boxHandl2 = boxplot(tmpMatrixMetApp, 'positions', pos);
    % grid minor;

    % boxesArray2 = findobj(boxHandl2);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray2(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray2(idx2(j)).Color = colors{2};
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

    % % manually adjust XTicks
    % sub1.XTick = [1.1, 1.6, 2.1, 2.6];
    % sub1.XTickLabel = {'1','2','3', '20'};
    % sub2.XTick = [1.1, 1.6, 2.1, 2.6];
    % sub2.XTickLabel = {'1','2','3', '20'};

    % xlabel('Iteration step', 'FontSize', fontSizes(1));
    % ylabel('\DeltaC [W]', 'FontSize', fontSizes(1));

    % % manual Figure adjustments
    % set(boxesArray1,'LineWidth', 1);
    % set(boxesArray2,'LineWidth', 1);

    % %% put ylabel right and rotate text
    % % set(sub2, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('C^{met}_{reduction} [%%]'), 'rot', -90, 'FontSize', fontSizes(1));
    % % p = get(lh, 'position');
    % % set(sub2, 'yaxislocation', 'left');
    % % set(lh, 'position', p);

    % % suptitle(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));
    % % suptitle(sprintf('Reduction of Vergence Error & Metabolic Costs\nvs. Iteration at Testing'));

    % % l = legend(subBoxHandl([2, 1]), captions);
    % % l.FontSize = fontSizes(1);
    % % l.Orientation = 'horizontal';
    % % l.Location = 'southoutside';

    % [l, objh, ~, ~] = legend(subBoxHandl([2, 1]), captions, 'Orientation', 'horizontal', 'Location', 'southoutside');
    % set(objh, 'linewidth', 2);
    % objh(3).Color = colors{1} / 0.9;
    % objh(5).Color = colors{2} / 0.9;

    % %% repositioning subfigures
    % sub1.Position(3 : 4) = sub2.Position(3 : 4);
    % sub1.Position(2) = 0.65;
    % l.Position(2) = 0.5;

    % plotpath = sprintf('%s/FigB2_VergErrMetCostsVsTestIter', savePath);
    % saveas(figB2, plotpath, 'png');
    % close(figB2);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % combines plot

    % figB3 = figure();
    % hold on;

    % steps = 3;              % show just first steps iterations & last iteration
    % colors = colors(1 : 2);
    % % colors{1} = colors{1} .* 0.9;
    % % colors{2} = colors{2} .* 0.9;
    % captions = {'M_{ }', 'M_{C}'};

    % % tmpMatrix = [vergErr_woMetCosts, vergErr_wMetCosts, metCostsApproach_woMetCosts, metCostsApproach_wMetCosts]
    % tmpMatrixVergErr = horzcat(dataMatrixEnd{1}{1}, dataMatrixEnd{2}{1});

    % % sort by iteration step
    % idx = [];
    % for (i = 1 : 20)
    %     idx(end + 1 : end + 2) = [i, i + 20];
    % end
    % tmpMatrixVergErr = tmpMatrixVergErr(:, idx);
    % tmpMatrixVergErr = [tmpMatrixVergErr(:, 1 : 2 * steps), tmpMatrixVergErr(:, end - 1 : end)];

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
    % tmpMatrixMetApp = tmpMatrixMetApp(:, idx);
    % tmpMatrixMetApp = [tmpMatrixMetApp(:, 1 : 2 * steps), tmpMatrixMetApp(:, end - 1 : end)];

    % sub1 = subplot(2, 1, 1);
    % pos = [1 1.2 1.5 1.7 2 2.2 2.5 2.7];
    % boxHandl = boxplot(tmpMatrixVergErr, 'positions', pos);
    % tmpHandle = findobj(boxHandl, 'type', 'text');
    % set(tmpHandle, 'Interpreter', 'tex');
    % grid minor;

    % subBoxHandl = findobj(gca,'Tag','Box');
    % % subBoxHandl = findobj(boxHandl,'Tag','Box');

    % boxesArray1 = findobj(boxHandl);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray1(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray1(idx2(j)).Color = colors{2};
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

    % ylabel('\Delta\gamma [deg]', 'FontSize', fontSizes(1));

    % %% put ylabel right and rotate text
    % % set(sub1, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('C^{met}_{reduction} [%%]'), 'rot', -90, 'FontSize', fontSizes(1));
    % % p = get(lh, 'position');
    % % set(sub1,'yaxislocation','left');
    % % set(lh,'position', p);

    % sub2 = subplot(2, 1, 2);
    % boxHandl2 = boxplot(tmpMatrixMetApp, 'positions', pos);
    % grid minor;

    % boxesArray2 = findobj(boxHandl2);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray2(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray2(idx2(j)).Color = colors{2};
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

    % % manually adjust XTicks
    % sub1.XTick = [1.1, 1.6, 2.1, 2.6];
    % sub1.XTickLabel = {'1','2','3', '20'};
    % sub2.XTick = [1.1, 1.6, 2.1, 2.6];
    % sub2.XTickLabel = {'1','2','3', '20'};

    % xlabel('Iteration step', 'FontSize', fontSizes(1));
    % ylabel(sprintf('c [%%]'), 'FontSize', fontSizes(1));

    % % manual Figure adjustments
    % set(boxesArray1,'LineWidth', 1);
    % set(boxesArray2,'LineWidth', 1);

    % %% put ylabel right and rotate text
    % % set(sub2, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('C^{met}_{reduction} [%%]'), 'rot', -90, 'FontSize', fontSizes(1));
    % % p = get(lh, 'position');
    % % set(sub2, 'yaxislocation', 'left');
    % % set(lh, 'position', p);

    % % suptitle(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));
    % % suptitle(sprintf('Reduction of Vergence Error & Metabolic Costs\nvs. Iteration at Testing'));

    % [l, objh, ~, ~] = legend(subBoxHandl([2, 1]), captions, 'Orientation', 'horizontal', 'Location', 'southoutside');
    % set(objh, 'linewidth', 2);
    % objh(3).Color = colors{1} / 0.9;
    % objh(5).Color = colors{2} / 0.9;

    % %% repositioning subfigures
    % sub1.Position(3 : 4) = sub2.Position(3 : 4);
    % sub1.Position(2) = 0.65;
    % l.Position(2) = 0.5;

    % plotpath = sprintf('%s/FigB3_VergErrMetCostsVsTestIter', savePath);
    % saveas(figB3, plotpath, 'png');
    % close(figB3);

    % fileID = fopen(strcat(savePath, '/README.txt'), 'at' );
    % fprintf(fileID, '========================================================================================================\n');
    % fprintf(fileID, 'testIter: %s %s\n\n', int2str([1, 1, 2, 2, 3, 3, 20, 20]), '= [w/o MetCosts, w/ MetCosts, w/o MetCosts, w/ MetCosts, ...]');

    % fprintf(fileID, 'figB3 vergErr median:\t%f %f %f %f %f %f %f %f\n', median(tmpMatrixVergErr));
    % fprintf(fileID, 'figB3 vergErr iqr:\t%f %f %f %f %f %f %f %f\n', iqr(tmpMatrixVergErr));
    % fprintf(fileID, 'figB3 vergErr mean:\t%f %f %f %f %f %f %f %f\n', mean(tmpMatrixVergErr));
    % fprintf(fileID, 'figB3 vergErr std:\t%f %f %f %f %f %f %f %f\n\n', std(tmpMatrixVergErr));

    % fprintf(fileID, 'figB3 metCosts median:\t%f %f %f %f %f %f %f %f\n', median(tmpMatrixMetApp));
    % fprintf(fileID, 'figB3 metCosts iqr:\t%f %f %f %f %f %f %f %f\n', iqr(tmpMatrixMetApp));
    % fprintf(fileID, 'figB3 metCosts mean:\t%f %f %f %f %f %f %f %f\n', mean(tmpMatrixMetApp));
    % fprintf(fileID, 'figB3 metCosts std:\t%f %f %f %f %f %f %f %f\n\n', std(tmpMatrixMetApp));

    % fclose(fileID);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % combines plot

    % figB4 = figure();
    % hold on;

    % steps = 3;                                          % show just first steps iterations & last iteration
    % colors = colors(1 : 2);
    % % colors{1} = colors{1} .* 0.9;
    % % colors{2} = colors{2} .* 0.9;
    % captions = {'M_{ }', 'M_{C}'};

    % % tmpMatrix = [vergErr_woMetCosts, vergErr_wMetCosts, metCostsApproach_woMetCosts, metCostsApproach_wMetCosts]
    % tmpMatrixVergErr = horzcat(dataMatrixEnd{1}{1}, dataMatrixEnd{2}{1});
    % tmpMatrixVergErr = abs(tmpMatrixVergErr);

    % iqrLine2 = zeros(4, 4);

    % % sort by iteration step
    % idx = [];
    % for (i = 1 : 20)
    %     idx(end + 1 : end + 2) = [i, i + 20];
    % end
    % tmpMatrixVergErr = tmpMatrixVergErr(:, idx);
    % tmpMatrixVergErr = [tmpMatrixVergErr(:, 1 : 2 * steps), tmpMatrixVergErr(:, end - 1 : end)];

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
    % tmpMatrixMetApp = tmpMatrixMetApp(:, idx);
    % tmpMatrixMetApp = [tmpMatrixMetApp(:, 1 : 2 * steps), tmpMatrixMetApp(:, end - 1 : end)];

    % sub1 = subplot(2, 1, 1);
    % pos = [1 1.2 1.5 1.7 2 2.2 2.5 2.7];
    % boxHandl = boxplot(tmpMatrixVergErr, 'positions', pos);
    % tmpHandle = findobj(boxHandl, 'type', 'text');
    % set(tmpHandle, 'Interpreter', 'tex');
    % % grid minor;

    % subBoxHandl = findobj(gca,'Tag','Box');
    % % subBoxHandl = findobj(boxHandl,'Tag','Box');

    % boxesArray1 = findobj(boxHandl);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray1(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray1(idx2(j)).Color = colors{2};
    %         end
    %     end
    % end

    % % remove outliers
    % outl = findobj(boxHandl, 'tag', 'Outliers');
    % set(outl, 'Visible', 'off');

    % % rescale axis to whiskers + offset
    % upWi = findobj(boxHandl, 'tag', 'Upper Whisker');
    % lowWi = findobj(boxHandl, 'tag', 'Lower Whisker');

    % % tmpMinY = round(max(arrayfun(@(x) x.YData(2), upWi)) * -0.11, 1);
    % % tmpMaxY = round(max(arrayfun(@(x) x.YData(2), upWi)) * 1.1, 1);

    % % axis([0.9, 2.8, tmpMinY, tmpMaxY]);
    % % set(gca, 'Ytick', linspace(tmpMinY, tmpMaxY, 6));
    % tmpMinY = 0;
    % tmpMaxY = 2.4;
    % axis([0.9, 2.8, tmpMinY, tmpMaxY]);
    % set(gca, 'Ytick', linspace(tmpMinY, tmpMaxY, 7));

    % ylabel('|\Delta\gamma| [deg]', 'FontSize', fontSizes(1));

    % % store values for FigB5
    % tmp = [arrayfun(@(x) x.YData(1), upWi)'; arrayfun(@(x) x.YData(2), lowWi)'];
    % iqrLine2(1, :) = tmp(1, 1 : 2 : end);
    % iqrLine2(2, :) = tmp(2, 1 : 2 : end);
    % iqrLine2(3, :) = tmp(1, 2 : 2 : end);
    % iqrLine2(4, :) = tmp(2, 2 : 2 : end);

    % %% put ylabel right and rotate text
    % % set(sub1, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('C^{met}_{reduction} [%%]'), 'rot', -90, 'FontSize', fontSizes(1));
    % % p = get(lh, 'position');
    % % set(sub1,'yaxislocation','left');
    % % set(lh,'position', p);

    % sub2 = subplot(2, 1, 2);
    % boxHandl2 = boxplot(tmpMatrixMetApp, 'positions', pos);
    % % grid minor;

    % boxesArray2 = findobj(boxHandl2);
    % for i = 1 : size(tmpMatrixVergErr, 2)
    %     idx2 = (1 : 7) + (i - 1) * 7;
    %     idx2(6 : 7) = [];
    %     if (mod(i, 2) == 1)
    %         for j = 1 : length(idx2)
    %             boxesArray2(idx2(j)).Color = colors{1};
    %         end
    %     else
    %         for j = 1 : length(idx2)
    %             boxesArray2(idx2(j)).Color = colors{2};
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

    % % manually adjust XTicks
    % sub1.XTick = [1.1, 1.6, 2.1, 2.6];
    % sub1.XTickLabel = {'1','2','3', '20'};
    % sub2.XTick = [1.1, 1.6, 2.1, 2.6];
    % sub2.XTickLabel = {'1','2','3', '20'};

    % xlabel('Iteration step', 'FontSize', fontSizes(1));
    % ylabel(sprintf('c [%%]'), 'FontSize', fontSizes(1));

    % % add 1-px patch as background
    % hp1 = patch([0, 3, 3, 0], [0.22, 0.22, 0, 0], ...
    %             [150 / 255, 150 / 255, 150 / 255], 'LineStyle', 'none', 'LineWidth', 2, 'FaceAlpha', 0.2);

    % hp1.Parent = sub1;
    % uistack(hp1,'bottom');

    % % manual Figure adjustments
    % set(boxesArray1,'LineWidth', 1);
    % set(boxesArray2,'LineWidth', 1);

    % %% put ylabel right and rotate text
    % % set(sub2, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('C^{met}_{reduction} [%%]'), 'rot', -90, 'FontSize', fontSizes(1));
    % % p = get(lh, 'position');
    % % set(sub2, 'yaxislocation', 'left');
    % % set(lh, 'position', p);

    % % suptitle(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));
    % % suptitle(sprintf('Reduction of Vergence Error & Metabolic Costs\nvs. Iteration at Testing'));

    % [l, objh, ~, ~] = legend(subBoxHandl([2, 1]), captions, 'Orientation', 'horizontal', 'Location', 'northeast');
    % set(objh, 'linewidth', 2);
    % objh(3).Color = colors{1} / 0.9;
    % objh(5).Color = colors{2} / 0.9;

    % %% repositioning subfigures
    % % sub1.Position(3 : 4) = sub2.Position(3 : 4);
    % % sub1.Position(2) = 0.65;
    % % l.Position(2) = 0.5;
    % l.Box = 'off';

    % plotpath = sprintf('%s/FigB4_VergErrMetCostsVsTestIter', savePath);
    % saveas(figB4, plotpath, 'png');
    % close(figB4);

    % fileID = fopen(strcat(savePath, '/README.txt'), 'at' );
    % fprintf(fileID, '========================================================================================================\n');
    % fprintf(fileID, 'testIter: %s %s\n\n', int2str([1, 1, 2, 2, 3, 3, 20, 20]), '= [w/o MetCosts, w/ MetCosts, w/o MetCosts, w/ MetCosts, ...]');

    % fprintf(fileID, 'figB4 vergErr median:\t%f %f %f %f %f %f %f %f\n', median(tmpMatrixVergErr));
    % fprintf(fileID, 'figB4 vergErr iqr:\t%f %f %f %f %f %f %f %f\n', iqr(tmpMatrixVergErr));
    % fprintf(fileID, 'figB4 vergErr mean:\t%f %f %f %f %f %f %f %f\n', mean(tmpMatrixVergErr));
    % fprintf(fileID, 'figB4 vergErr std:\t%f %f %f %f %f %f %f %f\n\n', std(tmpMatrixVergErr));

    % fprintf(fileID, 'figB4 metCosts median:\t%f %f %f %f %f %f %f %f\n', median(tmpMatrixMetApp));
    % fprintf(fileID, 'figB4 metCosts iqr:\t%f %f %f %f %f %f %f %f\n', iqr(tmpMatrixMetApp));
    % fprintf(fileID, 'figB4 metCosts mean:\t%f %f %f %f %f %f %f %f\n', mean(tmpMatrixMetApp));
    % fprintf(fileID, 'figB4 metCosts std:\t%f %f %f %f %f %f %f %f\n\n', std(tmpMatrixMetApp));

    % fclose(fileID);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figB5 = figure();
    % hold on;

    % steps = 3;                        % show just first steps iterations & last iteration
    % colors = colors(1 : 2);
    % % colors{1} = colors{1} .* 0.9;
    % % colors{2} = colors{2} .* 0.9;
    % captions = {'M_{ }', 'M_{C}'};
    % alphas = {0.2, 0.3};
    % xVal = 1 : 4;

    % objHandels = {};

    % sub1 = subplot(2, 1, 1);
    % pos = [1 1.2 1.5 1.7 2 2.2 2.5 2.7];
    % grid minor;

    % % fill area defined by upper & lower IQR bounds
    % objHandels{end + 1} = patch([xVal, flip(xVal)], [iqrLine2(1, xVal), flip(iqrLine2(2, xVal))], ...
    %                             colors{1}, 'LineStyle', 'none', 'LineWidth', 2, 'FaceAlpha', alphas{1});
    % hold on;

    % % median(|verg_{err}|) [deg]
    % objHandels{end + 1} = plot(sub1, xVal, median(tmpMatrixVergErr(:, 1 : 2 : end)), ...
    %                            'Color', colors{1});
    % hold on;

    % objHandels{end + 1} = patch([1 : 4, flip(1 : 4)], [iqrLine2(3, 1 : 4), flip(iqrLine2(4, 1 : 4))], ...
    %                             colors{2}, 'LineStyle', 'none', 'LineWidth', 2, 'FaceAlpha', alphas{2});
    % hold on;

    % objHandels{end + 1} = plot(sub1, xVal, median(tmpMatrixVergErr(:, 2 : 2 : end)), ...
    %                            'Color', colors{2});
    % hold on;
    % ylabel('|\Delta\gamma| [deg]', 'FontSize', fontSizes(1));

    % %% put ylabel right and rotate text
    % % set(sub1, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('C^{met}_{reduction} [%%]'), 'rot', -90, 'FontSize', fontSizes(1));
    % % p = get(lh, 'position');
    % % set(sub1,'yaxislocation','left');
    % % set(lh,'position', p);

    % sub2 = subplot(2, 1, 2);
    % boxHandl2 = boxplot(tmpMatrixMetApp, 'positions', pos);
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

    % % manually adjust XTicks
    % % sub1.XTick = [1.1, 1.6, 2.1, 2.6];
    % sub1.XTickLabel = {'1','','2','','3','', '20'};
    % sub2.XTick = [1.1, 1.6, 2.1, 2.6];
    % sub2.XTickLabel = {'1','2','3', '20'};

    % xlabel('Iteration step', 'FontSize', fontSizes(1));
    % ylabel(sprintf('c [%%]'), 'FontSize', fontSizes(1));

    % %% put ylabel right and rotate text
    % % set(sub2, 'yaxislocation', 'right');
    % % lh = ylabel(sprintf('C^{met}_{reduction} [%%]'), 'rot', -90, 'FontSize', fontSizes(1));
    % % p = get(lh, 'position');
    % % set(sub2, 'yaxislocation', 'left');
    % % set(lh, 'position', p);

    % % suptitle(sprintf('Total Vergence Error & Metabolic Costs Approach\nvs. Trial at Testing'));
    % % suptitle(sprintf('Reduction of Vergence Error & Metabolic Costs\nvs. Iteration at Testing'));

    % [l, objh, ~, ~] = legend([objHandels{2}, objHandels{4}], captions, 'Orientation', 'horizontal', 'Location', 'southoutside');
    % set(objh, 'linewidth', 2);
    % objh(3).Color = colors{1} / 0.9;
    % objh(5).Color = colors{2} / 0.9;

    % %% repositioning subfigures
    % sub1.Position(3 : 4) = sub2.Position(3 : 4);
    % sub1.Position(2) = 0.65;
    % l.Position(2) = 0.5;

    % plotpath = sprintf('%s/FigB5_VergErrMetCostsVsTestIter', savePath);
    % saveas(figB5, plotpath, 'png');
    % close(figB5);

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
    startVergErr = [-2, 1];
    startVergErr2 = [-2, 1];
    numIters = 20;
    stimuliIndices = [1];

    % initMethod elem. {0, 1}; 0 = getMFedoodD, 1 = fixed
    initMethod = uint8(0);

    % hand-picked inits for medial rectus for initMethod = 0 and 1
    % mrVal := {'w/o met. costs', 'w/  met. costs'} x objDist * startVergErr
    % mrVal = [[0.055, 0.085, 0.055, 0.085]; [0.03, 0.055, 0.03, 0.05]];
    mrVal = [[0.055, 0.08, 0.055, 0.08]; [0.03, 0.045, 0.03, 0.045]];

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
                    % [command, angleNew] = modelHandle(1, 1).model.getMFedoodD(objDist(odIndex), ...
                    %                                                        min(startVergErr(vergErrIndex), modelHandle(1, 1).model.getVergErrMax(objDist(odIndex))), ...
                    %                                                        mrVal(1, odIndex * 2 - 1 + vergErrIndex - 1));
                    [command, angleNew] = modelHandle(1, 1).model.getMFedoodD(objDist(odIndex), ...
                                                                           startVergErr2(vergErrIndex), ...
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
                    % [command, angleNew] = modelHandle(1, 2).model.getMFedoodD(objDist(odIndex), ...
                    %                                                        min(startVergErr(vergErrIndex), modelHandle(1, 2).model.getVergErrMax(objDist(odIndex))), ...
                    %                                                        mrVal(2, odIndex * 2 - 1 + vergErrIndex - 1));
                    [command, angleNew] = modelHandle(1, 2).model.getMFedoodD(objDist(odIndex), ...
                                                                           startVergErr(vergErrIndex), ...
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
    % colors = {{colors{1}, colors{1}, colors{1}, [1, 201 / 255, 41 / 255], 'k', 'y'}, ...
    %           {colors{2}, colors{2}, colors{2}, [1, 94 / 255, 41 / 255], 'k', 'm'}};

    % white end points
    % colors = {{[70 / 255, 160 / 255, 255 / 255], [70 / 255, 160 / 255, 255 / 255], [70 / 255, 160 / 255, 255 / 255], [1, 201 / 255, 41 / 255], 'k', [1, 1, 1]}, ...
    %           {[0, 200/255, 0], [0, 200/255, 0], [0, 200/255, 0], [1, 94 / 255, 41 / 255], 'k', [1, 1, 1]}};

    % dark end points
    % colors = {{[70 / 255, 160 / 255, 255 / 255], [70 / 255, 160 / 255, 255 / 255], [70 / 255, 160 / 255, 255 / 255], [1, 201 / 255, 41 / 255], 'k', [0, 100/255, 200/255]}, ...
    %           {[0, 200/255, 0], [0, 200/255, 0], [0, 200/255, 0], [1, 94 / 255, 41 / 255], 'k', [0, 95/255, 0]}};

    % colors = {{[70 / 255, 160 / 255, 255 / 255], [70 / 255, 160 / 255, 255 / 255], [70 / 255, 160 / 255, 255 / 255], [1, 201 / 255, 41 / 255], 'k', [0, 100/255, 200/255]}, ...
    %           {[0, 173/255, 87 / 255], [0, 173/255, 87 / 255], [0, 173/255, 87 / 255], [1, 94 / 255, 41 / 255], 'k', [0, 95/255, 0]}};

    % colors = {{[46 / 255, 147 / 255, 255 / 255], [46 / 255, 147 / 255, 255 / 255], [46 / 255, 147 / 255, 255 / 255], [1, 201 / 255, 41 / 255], 'k', [0, 110/255, 200/255]}, ...
    %           {[0, 173/255, 87 / 255], [0, 173/255, 87 / 255], [0, 173/255, 87 / 255], [1, 94 / 255, 41 / 255], 'k', [0, 110/255, 0]}};

    colors = {{[46 / 255, 147 / 255, 255 / 255], [46 / 255, 147 / 255, 255 / 255], [46 / 255, 147 / 255, 255 / 255], [1, 201 / 255, 41 / 255], 'k', 'r'}, ...
              {[0, 173/255, 87 / 255], [0, 173/255, 87 / 255], [0, 173/255, 87 / 255], [1, 94 / 255, 41 / 255], 'k', 'r'}};

    lineStyles = {'-', '-'};
    % markers = {'o', 'o'}; % {traj_line(1:10), end_fixation}
    markers = {'*', 'o'};
    % w/o met. costs: [line, end_point]
    % w/  met. costs: [line, end_point]
    markersizes = {3, 5, 3, 5};

    figC = figure('OuterPosition', [100, 100, 600, 600]);
    hold on;

    % title('Object Fixation Trajectories');
    xlabel('lateral rectus activation [%]');
    ylabel('medial rectus activation [%]');

    % pcHandle = pcolor(modelHandle(1, 1).model.degreesIncRes .* 2);    % use vergence degree as color dimension (background)
    pcHandle = pcolor(modelHandle(1, 1).model.metCostsIncRes .* 2);     % use metabolic costs as color dimension (background)
    % shading interp;
    set(pcHandle, 'EdgeColor', 'none');

    colormap(createCM(7));
    cb = colorbar();
    % cb.Label.String = 'vergence degree';  % use vergence degree as color dimension (background)
    cb.Label.String = 'C [W]';              % use metabolic costs as color dimension (background)
    cb.Label.FontSize = fontSizes(1);

    ax = gca;
    set(ax, 'Layer','top'); % bring axis to the front

    ax.XTick = linspace(1, size(modelHandle(1, 1).model.degreesIncRes, 2), 11);
    ax.YTick = linspace(1, size(modelHandle(1, 1).model.degreesIncRes, 1), 11);

    ax.XTickLabel = strsplit(num2str(linspace(0, 10, 11)));
    ax.YTickLabel = strsplit(num2str(linspace(0, 20, 11)));

    axis([1, size(modelHandle(1, 1).model.degreesIncRes, 2), 1, size(modelHandle(1, 1).model.degreesIncRes, 1)]);
    axPos = ax.Position;

    ht = cell(1, length(objDist)); % objDist label handles

    % draw lines of equal metCosts
    hlconstMetCosts = [];
    constMetCosts = [0.65, 0.95, 1.15];
    metCostsDiff = abs(max(max(diff(modelHandle(1, 2).model.metCostsIncRes))));
    modelHandle(1, 2).model.metCostsIncRes = interp2(modelHandle(1, 2).model.metCosts.results(1 : 3, 1 : 2), 10, 'spline');
    for k = 1 : length(constMetCosts)
        [yi, xi] = find(modelHandle(1, 2).model.metCostsIncRes >= constMetCosts(k) - metCostsDiff & modelHandle(1, 2).model.metCostsIncRes <= constMetCosts(k) + metCostsDiff);
        hlconstMetCosts(k) = plot(xi, yi, 'w');
    end

    % draw objects + offsets
    for odIndex = 1 : length(objDist)
        % draw +1 pixel offset in respect to desired vergence distance
        [lateralDes, medialDes] = modelHandle(1, 1).model.getAnglePoints(objDist(odIndex), 0.22);
        plot(lateralDes ./ modelHandle(1, 1).model.scaleFacLR, medialDes ./ modelHandle(1, 1).model.scaleFacMR, ...
             'color', 'k', 'LineStyle', '--', 'LineWidth', 0.75);

        % draw -1 pixel offset in respect to desired vergence distance
        [lateralDes, medialDes] = modelHandle(1, 1).model.getAnglePoints(objDist(odIndex), -0.22);
        plot(lateralDes ./ modelHandle(1, 1).model.scaleFacLR, medialDes ./ modelHandle(1, 1).model.scaleFacMR, ...
             'color', 'k', 'LineStyle', '--', 'LineWidth', 0.75);

        % draw a line of points into the plane that represent the desired vergence
        [lateralDes, medialDes] = modelHandle(1, 1).model.getAnglePoints(objDist(odIndex), 0);
        plot(lateralDes ./ modelHandle(1, 1).model.scaleFacLR, medialDes ./ modelHandle(1, 1).model.scaleFacMR, ...
                            'color', 'k', 'LineWidth', 1);

        % add corresponding distance value to desired vergence graph
        % ht{odIndex} = text(lateralDes(end - ceil(length(lateralDes) / 8.9)) / modelHandle(1, 1).model.scaleFacLR, ...
        %                    medialDes(end - ceil(length(medialDes) / 8.9)) / modelHandle(1, 1).model.scaleFacMR, ...
        %                    sprintf('%3.1f m', objDist(odIndex)), ...
        %                    'Fontsize', fontSizes(1));
    end

    % draw trajectories
    for odIndex = 1 : length(objDist)
        for stim = 1 : length(stimuliIndices)
            for vergErrIndex = 1 : length(startVergErr)
                % plot modification hack
                if (odIndex == 1 && vergErrIndex == 2) || (odIndex == 2 && vergErrIndex == 1)
                    continue;
                end

                % first plot whole trajectory without metabolic costs
                hl2 = plot(reshape(trajectory(odIndex, vergErrIndex, stim, :, 1), [numIters + 1, 1]) ./ modelHandle(1, 1).model.scaleFacLR + 1, ...
                           reshape(trajectory(odIndex, vergErrIndex, stim, :, 2), [numIters + 1, 1]) ./ modelHandle(1, 1).model.scaleFacMR + 1, ...
                           'Color', colors{1}{1}, 'LineStyle', lineStyles{1}, 'LineWidth', 2.1, ...%1
                           'Marker', markers{1}, 'MarkerEdgeColor', colors{1}{2}, 'MarkerFaceColor',  colors{1}{3}, 'MarkerSize', markersizes{1});%4

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
                           'LineStyle', 'none', 'Marker', markers{2}, 'MarkerEdgeColor', colors{1}{5}, 'MarkerFaceColor',  colors{1}{6}, 'MarkerSize', markersizes{2});%'LineWidth', 1, 'Marker', markers{2}, 'MarkerEdgeColor', colors{1}{5}, 'MarkerFaceColor',  colors{1}{6}, 'MarkerSize', 4);
            end
        end
    end

    for odIndex = 1 : length(objDist)
        for stim = 1 : length(stimuliIndices)
            for vergErrIndex = 1 : length(startVergErr)
                % plot modification hack
                if (odIndex == 1 && vergErrIndex == 2) || (odIndex == 2 && vergErrIndex == 1)
                    continue;
                end

                % first plot whole trajectory with metabolic costs
                hl5 = plot(reshape(trajectory2(odIndex, vergErrIndex, stim, :, 1), [numIters + 1, 1]) ./ modelHandle(1, 2).model.scaleFacLR + 1, ...
                           reshape(trajectory2(odIndex, vergErrIndex, stim, :, 2), [numIters + 1, 1]) ./ modelHandle(1, 2).model.scaleFacMR + 1, ...
                           'Color', colors{2}{1}, 'LineStyle', lineStyles{2}, 'LineWidth', 2.1, ...%1
                           'Marker', markers{1}, 'MarkerEdgeColor', colors{2}{2}, 'MarkerFaceColor',  colors{2}{3}, 'MarkerSize', markersizes{3});%4

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
                           'LineStyle', 'none', 'Marker', markers{2}, 'MarkerEdgeColor', colors{2}{5}, 'MarkerFaceColor',  colors{2}{6}, 'MarkerSize', markersizes{4});%'LineWidth', 1, 'Marker', markers{2}, 'MarkerEdgeColor', colors{2}{5}, 'MarkerFaceColor',  colors{2}{6}, 'MarkerSize', 4);
            end
        end
    end

    % realign plot order
    uistack(hl3,'top');
    uistack(hl6,'top');

    % gKey = {sprintf('1..%dth  iteration', modelHandle(1).model.interval), ...
    %         sprintf('%d..%dth iteration', modelHandle(1).model.interval, numIters), ...
    %         'end fixation w/o met. costs', ...
    %         sprintf('1..%dth  iteration', modelHandle(1).model.interval), ...
    %         sprintf('%d..%dth iteration', modelHandle(1).model.interval, numIters), ...
    %         'end fixation w/  met. costs'};

    % gKey = {sprintf('0th..%dth iteration', numIters), ...
    %         'end fixation  \bfw/o met. costs', ...
    %         sprintf('0th..%dth iteration', numIters), ...
    %         'end fixation  \bfw/   met. costs'};

    % gKey = {strcat(sprintf('0th..%dth iteration', numIters), ' \bfw/o met. costs'),
    %         strcat(sprintf('0th..%dth iteration', numIters), ' \bfw/   met. costs')};
    % gKey = {strcat(sprintf('0th..%dth iteration', numIters), ' M'),
    %         strcat(sprintf('0th..%dth iteration', numIters), ' M_{C}')};
    gKey = {' M_{ }', 'M_{C}'};

    % hDummy1 = plot(NaN, NaN, 'LineStyle', 'none');
    % hDummy2 = plot(NaN, NaN, 'LineStyle', 'none');

    % with dummy entries for last 2 columns
    % gKey = {sprintf('0th..%dth iteration', numIters), ...
    %         'end fixation', ...
    %         '\bfw/o met. costs', ...
    %         sprintf('0th..%dth iteration', numIters), ...
    %         'end fixation', ...
    %         '\bfw/   met. costs'};

    l = gridLegend([hl2, hl5], 1, gKey, 'Location', 'northwest', 'Fontsize', fontSizes(1));
    % l = gridLegend([hl2, hl3, hl5, hl6], 2, gKey, 'Orientation', 'Horizontal', 'Location', 'southoutside', 'Fontsize', fontSizes(1));
    % l = gridLegend([hl2, hl3, hDummy1, hl5, hl6, hDummy2], 3, gKey, 'Orientation', 'Horizontal', 'Location', 'southoutside', 'Fontsize', fontSizes(1));
    % l.Box = 'off';

    %% repositioning subfigures
    % ax.Position = axPos;
    % ax.PlotBoxAspectRatioMode = 'manual';
    % ax.DataAspectRatioMode = 'manual';
    % ax.ActivePositionProperty = 'position';
    % ax.OuterPosition(4) = ax.OuterPosition(4) * 1.1;
    ax.DataAspectRatioMode = 'manual';

%     tmp = ht{1}.Position(2) - ht{1}.Position(2) * 0.96;
%     ht{1}.Position(2) = ht{1}.Position(2) * 0.96;
%     ht{2}.Position(2) = ht{2}.Position(2) - tmp;

    % manual positioning
    ax.Position = [0.085, 0.2, 0.75, 0.75];
    % l.Position(1) = 0.15;
    % l.Position(2) = 0.04;

    % set(figC,'PaperPositionMode','auto'); % keep aspect ratio
    plotpath = sprintf('%s/FigC_muscleActivityTrajectories', savePath);
    saveas(figC, plotpath, 'png');
    close(figC);

    % store start- & end-points of trajectories
    fileID = fopen(strcat(savePath, '/README.txt'), 'at' );
    fprintf(fileID, '========================================================================================================\n');
    fprintf(fileID, '[mr, lr activation] start, metCost start, [mr, lr activation] end, metCost end\n\n');

    for odIndex = 1 : length(objDist)
        for stim = 1 : length(stimuliIndices)
            for vergErrIndex = 1 : length(startVergErr)
                fprintf(fileID, 'figC model w/o MetCosts:\t[%f %f] %f, [%f %f] %f\n', ...
                    trajectory(odIndex, vergErrIndex, stim, 1, 2), ...
                    trajectory(odIndex, vergErrIndex, stim, 1, 1), ...
                    modelHandle(1, 1).model.getMetCost(trajectory(odIndex, vergErrIndex, stim, 1, :)) * 2, ...
                    trajectory(odIndex, vergErrIndex, stim, end, 2), ...
                    trajectory(odIndex, vergErrIndex, stim, end, 1), ...
                    modelHandle(1, 1).model.getMetCost(trajectory(odIndex, vergErrIndex, stim, end, :)) * 2);
            end
        end
    end

    fprintf(fileID, '\n');

    for odIndex = 1 : length(objDist)
        for stim = 1 : length(stimuliIndices)
            for vergErrIndex = 1 : length(startVergErr)
                fprintf(fileID, 'figC model w/  MetCosts:\t[%f %f] %f, [%f %f] %f\n', ...
                    trajectory2(odIndex, vergErrIndex, stim, 1, 2), ...
                    trajectory2(odIndex, vergErrIndex, stim, 1, 1), ...
                    modelHandle(1, 2).model.getMetCost(trajectory2(odIndex, vergErrIndex, stim, 1, :)) * 2, ...
                    trajectory2(odIndex, vergErrIndex, stim, end, 2), ...
                    trajectory2(odIndex, vergErrIndex, stim, end, 1), ...
                    modelHandle(1, 2).model.getMetCost(trajectory2(odIndex, vergErrIndex, stim, end, :)) * 2);
            end
        end
    end

    fprintf(fileID, '\n\n');
    fclose(fileID);

    function updateMetcostApproach(model)
        %% Desired (updated) metabolic costs approach [%] vs. iteration
        objRange = [0.5, 1 : 6];

        metCostsApproach = zeros(size(model.testResult5, 1) / model.testInterval, model.testInterval);
        metCost = zeros(1, model.testInterval);

        % calculate all desired vergence angles and metabolic costs
        angleDes = objRange';
        metCostDesired = objRange';
        for odIndex = 1 : length(objRange)
            [cmdDesired, ~] = model.getMF2(objRange(odIndex), 0);
            metCostDesired(odIndex) = model.getMetCost(cmdDesired) * 2;
        end
        metCostDesired = repelem(metCostDesired, 7 * nStim, 1);

        for trial = 1 : size(metCostsApproach, 1)

            % starting point in muscle space
            cmdStart = [model.testResult5(trial * model.testInterval - model.testInterval + 1, 1) - model.testResult5(trial * model.testInterval - model.testInterval + 1, 3); ...
                        model.testResult5(trial * model.testInterval - model.testInterval + 1, 2) - model.testResult5(trial * model.testInterval - model.testInterval + 1, 4)];

            metCostStart = model.getMetCost(cmdStart) * 2;

            % deltaMetabolicCostsDesired = metCostDesired - metCostStart
            deltaMetCostDesired = metCostDesired(trial) - metCostStart;

            for iter = 1 : model.testInterval
                % metabolicCosts(iteration)
                metCost(iter) = model.getMetCost([model.testResult5(trial * model.testInterval - model.testInterval + iter, 1); ...
                                                  model.testResult5(trial * model.testInterval - model.testInterval + iter, 2)]) * 2;
            end

            % metCostsApproach(iteration) = (metCost(iteration) - metCostStart) / deltaDesired
            % metCostsApproach(trial, :) = (metCost - metCostStart) / deltaMetCostDesired;
            metCostsApproach(trial, :) = (metCost - metCostDesired(trial)) / metCostDesired(trial);
        end
        % scale to [%]
        % metCostsApproach = metCostsApproach .* 100;
        model.metCostsApproach = metCostsApproach;
    end
end

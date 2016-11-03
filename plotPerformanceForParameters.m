%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script reads in a number of models
%% and plots their performance according
%% to a specified set of two paramters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotPerformanceForParameters(modelAt)
    % folder with all subfolders containing the experiments
    % parentFolder = '/home/aecgroup/aecdata/Results/';
    % parentFolder = '/home/aecgroup/aecdata/Results/Regularizer vs Actor Learning Rate';
    % parentFolder = '/home/aecgroup/aecdata/Results/Discount Factor vs Interval';
    parentFolder = '/home/aecgroup/aecdata/Results/CriticLR vs ActorLR';

    % a string (or part of it) all relevant folders share
    commonName = 'CriticLR_';
    % commonName = '0';

    files = dir(sprintf('%s/*%s*', parentFolder, commonName));
    if isempty(modelAt)
        modelAt = 2000000;
    end
    nFiles = length(files);

    %% here, specify the parameter ranges that should be used
    %  these may simply be copied from parOES.m and putting ';'' after every set of params

    %% regularizer vs actor leaning range

    % labelVar1 = 'Actor weight regularizer';
    % labelVar2 = 'Actor LR [start, end]';

    % var1 = [1e-2; 1e-3; 5e-4; 1e-4; 5e-5; 1e-5; 1e-6]; % regularizer
    % var2 = [[1, 0]; [0.5, 0]; [0.5, 0.5]]; % actorLearningRange

    % numberFormatVar1 = '%1.0e';
    % numberFormatVar2 = '[%1.2f - %1.2f]';

    % plotSavePath = strcat(parentFolder, '/hiddenLayerRegulActorLRComparison');


    %% discount factor vs interval

    % var1 = [0.1; 0.3; 0.9];
    % var2 = [10; 50; 100];

    % labelVar1 = 'discount factor';
    % labelVar2 = 'interval';

    % numberFormatVar1 = '%1.1f';
    % numberFormatVar2 = '%d';

    % plotSavePath = strcat(parentFolder, '/DiscountVsInterval');

    %% actor LR vs critic LR

    labelVar1 = 'actor learning range';
    labelVar2 = 'critic learning range';

    var1 = [[1, 1]; [1, 0]; [0.75, 0.75]; [0.75, 0]; [0.5, 0.5]; [0.5, 0]; [0.25, 0.25]; [0.25, 0]];
    var2 = [[1, 1]; [1, 0]; [0.75, 0.75]; [0.75, 0]; [0.5, 0.5]; [0.5, 0]; [0.25, 0.25]; [0.25, 0]];
    var1 = flipud(var1); % setting the offspring to lower left had corner

    numberFormatVar1 = '[%1.2f - %1.2f]';
    numberFormatVar2 = '[%1.2f - %1.2f]';

    plotSavePath = strcat(parentFolder, '/CriticLRActorLR');

    %%%%%%% after this, you just have to update the extraction of values in the for loop %%%%%%%%%%%

    length1 = length(var1);
    length2 = length(var2);

    results = zeros(length1, length2, 5); % results are the three measurements rmse, median, and iqr, critValDelta, critValNiveau
    subFolder = sprintf('modelAt%d', modelAt);

    % iter = 1;
    for f = 1 : nFiles
        try
            % model = load(files{f});
            model = load(sprintf('%s/%s/%s/model.mat', parentFolder, files(f).name, subFolder));
            model = model.model;
            testInterval = model.interval * 2;
            % testInterval = 20;

            % hack for some older simulations:
            % if length(model.rlModel.actorLearningRange) == 1
            %     value = model.rlModel.actorLearningRange;
            %     model.rlModel.actorLearningRange = [value, value];
            % end

            % sprintf(model.savePath)
            % sprintf('%d', model.rlModel.CActor.regularizer)

            %% finding indizes: also needs to be updated everytime
            % ind = find(var1 == model.rlModel.CActor.varianceRange(1));
            % jnd = find(var2 == model.rlModel.CActor.varianceRange(2));

            % ind = find(ismember(var2, model.rlModel.actorLearningRange, 'rows'));
            % jnd = find(ismember(var1, model.rlModel.CActor.regularizer, 'rows'));

            % note: different order than expected
            ind = find(ismember(var2, model.rlModel.criticLearningRange, 'rows'));
            jnd = find(ismember(var1, model.rlModel.actorLearningRange, 'rows'));

            % ind = find(ismember(var2, model.interval, 'rows'));
            % jnd = find(ismember(var1, model.rlModel.CCritic.gamma, 'rows'));

            results(ind, jnd, 1) = sqrt(mean(model.testResult3(:, testInterval) .^ 2));
            results(ind, jnd, 2) = iqr(model.testResult3(:, testInterval)) * 4;
            results(ind, jnd, 3) = median(model.testResult3(:, testInterval));

            %%% Critic value function steepness
            % critValDelta = mean(critic_value(vergErr = 0) - (critic_value(vergErr = -0.5) + critic_value(vergErr = 0.5)) / 2)
            % critValNiveau = mean(critic_value(vergErr = 0))
            % average over all stimuli and each objDist
            test2Resolution = 101;
            vseRange = linspace(-1, 1, test2Resolution);

            % critValDelta
            results(ind, jnd, 4) = mean(mean(abs(abs(model.testResult4(:, vseRange == 0, 1 : 2 + length(model.scModel) : end)) ...
                                               - abs((model.testResult4(:, vseRange == -0.5, 1 : 2 + length(model.scModel) : end) ...
                                                     + model.testResult4(:, vseRange == 0.5, 1 : 2 + length(model.scModel) : end)) / 2)), 3));
            % critValNiveau
            results(ind, jnd, 5) = mean(mean(abs(model.testResult4(:, vseRange == 0, 1 : 2 + length(model.scModel) : end)), 3));
        catch
           % catch case when (sub-)experiment started, but has no test results yet
           continue;
        end
    end

    %% plotting section
    % var1descr = [];
    % var2descr = '';

    % [~, nEntries] = size(var1);
    % for ind1 = 1 : length1
    %   var1descr(ind1) = '[';
    %   for ind2 = 1 : nEntries
    %       var1descr(ind1) = strcat(var1descr(ind), var1(ind1, ind2))
    %       if ind2 < nEntries
    %           var1descr(ind1) = strcat(var1descr(ind1), ',');
    %       end
    %   end
    %   var1descr(ind1) = strcat(var1descr(ind1), ']');
    % end

    % mark unfinished experiments with white color code
    colordata = createCM(3);
    if (~isempty(results(results == 0)))
        colordata(1, :) = [1, 1, 1];
    end

    figure;
    % suptitle(sprintf('Parameter Comparison at %d iterations', modelAt));

    % plot the RMSE
    subplot(3, 1, 1);
    colormap(colordata);
    imagesc(results(:, :, 1));

    txt = results(:, :, 1);
    txt(txt == 0) = Inf;
    txt = num2str(txt(:),'%0.2f');
    txt = strtrim(cellstr(txt));
    [x, y] = meshgrid(1 : max(length1, length2));
    hStrings = text(x(:), y(:), txt(:), 'HorizontalAlignment', 'center');

    % one dim. case
    % set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(flip(var2(:)), numberFormatVar1), ...
    %          'YTick', 1:length2, 'YTickLabel', num2str(var1(:), numberFormatVar2), ...
    %          'TickLength', [0, 0]);
    % two dim case
    set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(var1(:, :), numberFormatVar1), ...
             'YTick', 1:length2, 'YTickLabel', num2str(var2(:, :), numberFormatVar2), ...
             'TickLength', [0, 0], 'FontSize', 7);
    % set(gca, 'XTick', 1:length1, 'XTickLabel', var1, ...
    %          'YTick', 1:length2, 'YTickLabel', var2, ...
    %          'TickLength', [0, 0]);

    title('RMSE');
    xlabel(sprintf(labelVar1));
    ylabel(sprintf(labelVar2));
    colorbar();
    % saveas(gca, sprintf('%s_rmse.png', plotSavePath));

    % plot the interquartile range
    subplot(3, 1, 2);
    colormap(colordata);
    imagesc(results(:, :, 2));

    txt = results(:, :, 2);
    txt(txt == 0) = Inf;
    txt = num2str(txt(:),'%0.2f');
    % txt(find(txt == '0.00')) = ' '; % this may cause errors.
    txt = strtrim(cellstr(txt));
    [x, y] = meshgrid(1 : max(length1, length2));
    hStrings = text(x(:), y(:), txt(:), 'HorizontalAlignment', 'center');

    % one dim. case
    % set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(flip(var2(:)), numberFormatVar1), ...
    %          'YTick', 1:length2, 'YTickLabel', num2str(var1(:), numberFormatVar2), ...
    %          'TickLength', [0, 0]);
    % two dim case
    set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(var1(:, :), numberFormatVar1), ...
             'YTick', 1:length2, 'YTickLabel', num2str(var2(:, :), numberFormatVar2), ...
             'TickLength', [0, 0], 'FontSize', 7);
    % set(gca, 'XTick', 1:length1, 'XTickLabel', var1, ...
    %          'YTick', 1:length2, 'YTickLabel', var2, ...
    %          'TickLength', [0, 0]);

    title('IQR*4');
    xlabel(sprintf(labelVar1));
    ylabel(sprintf(labelVar2));
    colorbar();
    % saveas(gca, sprintf('%s_iqr.png', plotSavePath));

    % plot the median
    subplot(3, 1, 3);
    colormap(colordata);
    imagesc(abs(results(:, :, 3))); % plot absolute values

    txt = results(:, :, 3);
    txt(txt == 0) = Inf;
    txt = num2str(txt(:),'%0.2f');
    % txt(find(txt == '0.00')) = ' '; % this may cause errors.
    txt = strtrim(cellstr(txt));
    [x, y] = meshgrid(1 : max(length1, length2));
    hStrings = text(x(:), y(:), txt(:), 'HorizontalAlignment', 'center');

    % one dim. case
    % set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(flip(var2(:)), numberFormatVar1), ...
    %          'YTick', 1:length2, 'YTickLabel', num2str(var1(:), numberFormatVar2), ...
    %          'TickLength', [0, 0]);
    % two dim case
    set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(var1(:, :), numberFormatVar1), ...
             'YTick', 1:length2, 'YTickLabel', num2str(var2(:, :), numberFormatVar2), ...
             'TickLength', [0, 0], 'FontSize', 7);
    % set(gca, 'XTick', 1:length1, 'XTickLabel', var1, ...
    %          'YTick', 1:length2, 'YTickLabel', var2, ...
    %          'TickLength', [0, 0]);

    title('Median');
    xlabel(sprintf(labelVar1));
    ylabel(sprintf(labelVar2));
    colorbar();
    % saveas(gca, sprintf('%s_median.png', plotSavePath));

    saveas(gca, sprintf('%s_at%diter.png', plotSavePath, modelAt));

    %%% Critic's value function
    colordata = createCM(3);
    colordata = flipud(colordata);
    if (~isempty(results(results == 0)))
        colordata(1, :) = [1, 1, 1];
    end

    figure;
    % suptitle(sprintf('Parameter Comparison at %d iterations', modelAt));

    % plot the critValDelta
    subplot(2, 1, 1);
    colormap(colordata);
    imagesc(results(:, :, 4));

    txt = results(:, :, 4);
    txt(txt == 0) = Inf;
    txt = num2str(txt(:),'%0.3f');
    txt = strtrim(cellstr(txt));
    [x, y] = meshgrid(1 : max(length1, length2));
    hStrings = text(x(:), y(:), txt(:), 'HorizontalAlignment', 'center');

    % one dim. case
    % set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(flip(var2(:)), numberFormatVar1), ...
    %          'YTick', 1:length2, 'YTickLabel', num2str(var1(:), numberFormatVar2), ...
    %          'TickLength', [0, 0]);
    % two dim case
    set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(var1(:, :), numberFormatVar1), ...
             'YTick', 1:length2, 'YTickLabel', num2str(var2(:, :), numberFormatVar2), ...
             'TickLength', [0, 0], 'FontSize', 7);
    % set(gca, 'XTick', 1:length1, 'XTickLabel', var1, ...
    %          'YTick', 1:length2, 'YTickLabel', var2, ...
    %          'TickLength', [0, 0]);

    title(strcat('\Deltacritic_{val}', sprintf(' = |mean(critic_{val}(verg_{Err} = 0)\n - (critic_{val}(verg_{Err} = -0.5) + critic_{val}(verg_{Err} = 0.5)) / 2)|')));
    xlabel(sprintf(labelVar1));
    ylabel(sprintf(labelVar2));
    colorbar();

    % plot the critValNiveau
    subplot(2, 1, 2);
    colormap(colordata);
    imagesc(results(:, :, 5));

    txt = results(:, :, 5);
    txt(txt == 0) = Inf;
    txt = num2str(txt(:),'%0.3f');
    % txt(find(txt == '0.00')) = ' '; % this may cause errors.
    txt = strtrim(cellstr(txt));
    [x, y] = meshgrid(1 : max(length1, length2));
    hStrings = text(x(:), y(:), txt(:), 'HorizontalAlignment', 'center');

    % one dim. case
    % set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(flip(var2(:)), numberFormatVar1), ...
    %          'YTick', 1:length2, 'YTickLabel', num2str(var1(:), numberFormatVar2), ...
    %          'TickLength', [0, 0]);
    % two dim case
    set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(var1(:, :), numberFormatVar1), ...
             'YTick', 1:length2, 'YTickLabel', num2str(var2(:, :), numberFormatVar2), ...
             'TickLength', [0, 0], 'FontSize', 7);
    % set(gca, 'XTick', 1:length1, 'XTickLabel', var1, ...
    %          'YTick', 1:length2, 'YTickLabel', var2, ...
    %          'TickLength', [0, 0]);

    title(sprintf('critic_{val} Niveau = |mean(critic_{val}(verg_{Err} = 0))|'));
    xlabel(sprintf(labelVar1));
    ylabel(sprintf(labelVar2));
    colorbar();

    saveas(gca, sprintf('%s_CriticVal_at%diter.png', plotSavePath, modelAt));
end

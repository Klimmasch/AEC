%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script reads in a number of models
%% and plots their performance according
%% to a specified set of two paramters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotPerformanceForParameters(modelAt)
    % ---------------------
    % Experiment definition
    % ---------------------

    %%% Definition Structure
    % 1)    folder with all subfolders containing the experiments
    % 2,3)  x (=labelVar1), y (=labelVar2) axis labels
    % 4,5)  x (=var1), y (=var2) tick values
    % 6,7)  x (=numberFormatVar1), y (=numberFormatVar2) tick values number formating
    % 8)    figure name

    %%% actor LR vs critic LR
    parentFolder = '/home/aecgroup/aecdata/Results/CriticLR vs ActorLR';

    labelVar1 = 'actor learning range';
    labelVar2 = 'critic learning range';

    var1 = [[1, 1]; [1, 0]; [0.75, 0.75]; [0.75, 0]; [0.5, 0.5]; [0.5, 0]; [0.25, 0.25]; [0.25, 0]];
    var2 = [[1, 1]; [1, 0]; [0.75, 0.75]; [0.75, 0]; [0.5, 0.5]; [0.5, 0]; [0.25, 0.25]; [0.25, 0]];
    var1 = flipud(var1); % setting the offspring to lower left hand corner

    numberFormatVar1 = '[%1.2f - %1.2f]';
    numberFormatVar2 = '[%1.2f - %1.2f]';

    figName = 'CriticLRActorLR';


    %%% regularizer vs actor leaning range
    % parentFolder = '/home/aecgroup/aecdata/Results/Regularizer vs Actor Learning Rate';

    labelVar1 = 'Actor weight regularizer';
    labelVar2 = 'Actor LR [start, end]';

    var1 = [1e-2; 1e-3; 5e-4; 1e-4; 5e-5; 1e-5; 1e-6]; % regularizer
    var2 = [[1, 0]; [0.5, 0]; [0.5, 0.5]]; % actorLearningRange

    numberFormatVar1 = '%1.0e';
    numberFormatVar2 = '[%1.2f - %1.2f]';


    % figName = 'hiddenLayerRegulActorLRComparison';

    %%% discount factor vs interval
    % parentFolder = '/home/aecgroup/aecdata/Results/Gamma_vs_Interval_fewerResources';

    % labelVar1 = 'discount factor';
    % labelVar2 = 'interval';

    % var1 = [0.1; 0.3; 0.9];
    % var2 = [10; 50; 100];

    % numberFormatVar1 = '%1.1f';
    % numberFormatVar2 = '%d';

    % figName = 'DiscountVsInterval';

    %%% metabolic costs
    % parentFolder = '/home/klimmasch/projects/results/exploringMetCost';


    % labelVar1 = 'metabolic Costs';
    % labelVar2 = '';

    % var1 = [[0.0324, 0.0324]; [0.0162, 0.0162]; [0.0121, 0.0121]; [0.0081, 0.0081]; [0.0040, 0.0040]; [0,0]];
    % var2 = [1e-4]; % only chooses values with the right regularizer

    % numberFormatVar1 = '[%1.4f,%1.4f]';
    % numberFormatVar2 = '%';

    % figName = 'metCostComp';

    % ------------------
    % Results extraction
    % ------------------

    subExperiments = dir(parentFolder);             % get a list of all files and folders in this folder
    dirFlags = [subExperiments.isdir];              % get a logical vector that tells which is a directory
    subExperiments = subExperiments(dirFlags);      % extract only those that are directories
    subExperiments = subExperiments(3 : end);       % exclude '.' and '..' directories

    if (~exist('modelAt'))
        modelAt = 2000000;
    end
    subFolder = sprintf('modelAt%d', modelAt);

    plotSavePath = strcat(parentFolder, '/', figName);

    length1 = length(var1);
    length2 = length(var2);
    criticResolution = 101;
    vseRange = linspace(-1, 1, criticResolution);


    % results := |var1| x |var2| x
    %            {rmse(vergErr), median(vergErr), and iqr(vergErr),
    %             rmse(deltaMC), median(deltaMC), and iqr(deltaMC),
    %             critValDelta, critValNiveau}
    results = zeros(length1, length2, 8);

    for f = 1 : length(subExperiments)
        try
            model = load(sprintf('%s/%s/%s/model.mat', parentFolder, subExperiments(f).name, subFolder));
        catch
           % catch case when (sub-)experiment started, but has no test results yet
           continue;
        end
        model = model.model;

        testInterval = model.interval * 2;
        % testInterval = 200;

        % hack for some older simulations:
        % if length(model.rlModel.actorLearningRange) == 1
        %     value = model.rlModel.actorLearningRange;
        %     model.rlModel.actorLearningRange = [value, value];
        % end

        %%% Extract indicies
        %%% actor exploration variance range
        % ind = find(var1 == model.rlModel.CActor.varianceRange(1));
        % jnd = find(var2 == model.rlModel.CActor.varianceRange(2));

        %%% regularizer vs actor leaning range
        % ind = find(ismember(var2, model.rlModel.actorLearningRange, 'rows'));
        % jnd = find(ismember(var1, model.rlModel.CActor.regularizer, 'rows'));

        %%% actor LR vs critic LR
        ind = find(ismember(var2, model.rlModel.criticLearningRange, 'rows'));
        jnd = find(ismember(var1, model.rlModel.actorLearningRange, 'rows'));

        %%% discount factor vs interval
        % ind = find(ismember(var2, model.interval, 'rows'));
        % jnd = find(ismember(var1, model.rlModel.CCritic.gamma, 'rows'));

        %%% regularizer vs metabolic costs
        % ind = find(ismember(var2, model.rlModel.CActor.regularizer, 'rows'));
        % jnd = find(ismember(var1, model.metCostRange, 'rows'));

        %%% Extract results
        %%% Vergence Error
        results(ind, jnd, 1) = sqrt(mean(model.testResult3(:, testInterval) .^ 2));
        results(ind, jnd, 2) = iqr(model.testResult3(:, testInterval)) * 4;
        results(ind, jnd, 3) = median(model.testResult3(:, testInterval));

        %%% delta Metabolic costs
        results(ind, jnd, 4) = sqrt(mean(model.testResult7(:, testInterval) .^ 2));
        results(ind, jnd, 5) = iqr(model.testResult7(:, testInterval)) * 4;
        results(ind, jnd, 6) = median(model.testResult7(:, testInterval));

        %%% Critic value function steepness
        % critValDelta = mean(critic_value(vergErr = 0) - (critic_value(vergErr = -0.5) + critic_value(vergErr = 0.5)) / 2)
        % critValNiveau = mean(critic_value(vergErr = 0))
        % average over all stimuli and each objDist

        % critValDelta
        results(ind, jnd, 7) = mean(mean(abs(abs(model.testResult4(:, vseRange == 0, 1 : 2 + length(model.scModel) : end)) ...
                                           - abs((model.testResult4(:, vseRange == -0.5, 1 : 2 + length(model.scModel) : end) ...
                                                 + model.testResult4(:, vseRange == 0.5, 1 : 2 + length(model.scModel) : end)) / 2)), 3));
        % critValNiveau
        results(ind, jnd, 8) = mean(mean(abs(model.testResult4(:, vseRange == 0, 1 : 2 + length(model.scModel) : end)), 3));
    end

    % ----------------
    % Results plotting
    % ----------------

    % mark unfinished experiments with white color code
    colordata = createCM(3);
    if (~isempty(results(results == 0)))
        colordata(1, :) = [1, 1, 1];
    end

    [x, y] = meshgrid(1 : length1, 1 : length2);
    titleStrings = {'RMSE', 'IQR*4', 'Median'};
    resultsIter = {[1 : 3]; [4 : 6]; [7, 8]};

    for figIter = 1 : 3
        if (figIter == 3)
            colordata = createCM(3);
            colordata = flipud(colordata);
            if (~isempty(results(results == 0)))
                colordata(1, :) = [1, 1, 1];
            end
        end

        figure;
        for subfigIter = 1 : length(resultsIter{figIter})
            subplot(length(resultsIter{figIter}), 1, subfigIter);
            colormap(colordata);
            if (figIter == 1 && subfigIter == 3)
                % apply color to abs values in median plot
                imagesc(abs(results(:, :, resultsIter{figIter}(subfigIter))));
            else
                imagesc(results(:, :, resultsIter{figIter}(subfigIter)));
            end

            txt = results(:, :, resultsIter{figIter}(subfigIter));
            txt(txt == 0) = Inf;
            if (figIter < 3)
                txt = num2str(txt(:), '%0.2f');
            else
                txt = num2str(txt(:), '%0.3f');
            end
            txt = strtrim(cellstr(txt));
            text(x(:), y(:), txt(:), 'HorizontalAlignment', 'center');

            if (size(var1, 2) == 2)
                % two dim case
                set(gca, 'XTick', 1 : length1, 'XTickLabel', num2str(var1(:, :), numberFormatVar1), ...
                         'YTick', 1 : length2, 'YTickLabel', num2str(var2(:, :), numberFormatVar2), ...
                         'TickLength', [0, 0], 'FontSize', 7);
            elseif (size(var1, 2) == 1)
                % one dim. case
                set(gca, 'XTick', 1 : length1, 'XTickLabel', var1, ...
                         'YTick', 1 : length2, 'YTickLabel', var2, ...
                         'TickLength', [0, 0]);
            else
                error('Parameter range var1 has unsupported size of %d x %d.', size(var1));
            end

            if (figIter < 3)
                title(titleStrings{subfigIter});
            elseif (subfigIter == 1)
                title(strcat('\Deltacritic_{val}', ...
                             sprintf(' = |mean(critic_{val}(verg_{Err} = 0)\n - (critic_{val}(verg_{Err} = -0.5)'), ...
                             ' + critic_{val}(verg_{Err} = 0.5)) / 2)|'));
            else
                title(sprintf('critic_{val} Niveau = |mean(critic_{val}(verg_{Err} = 0))|'));
            end

            xlabel(sprintf(labelVar1));
            ylabel(sprintf(labelVar2));
            colorbar();
        end


        if (figIter == 1)
            % suptitle(sprintf('Parameter Comparison at %d iterations', modelAt));
            saveas(gca, sprintf('%s_at%diter.png', plotSavePath, modelAt));
        elseif (figIter == 2)
            suptitle(sprintf('Metabolic Costs at %d iterations', modelAt));
            saveas(gca, sprintf('%s_metCost_at%diter.png', plotSavePath, modelAt));
        elseif (figIter == 3)
            saveas(gca, sprintf('%s_CriticVal_at%diter.png', plotSavePath, modelAt));
        end
    end
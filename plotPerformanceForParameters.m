%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script reads in a number of models
%% and plots their performance according
%% to a specified set of two paramters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotPerformanceForParameters()
    files = { ...
        '/home/aecgroup/aecdata/Results/model_11-Oct-2016_15:35:13_2000000_1_cluster_varDec1e4-1e4/modelAt2000000/model.mat', ...
        '/home/aecgroup/aecdata/Results/model_11-Oct-2016_15:40:58_2000000_1_cluster_varDec1e4-1e5/modelAt2000000/model.mat', ...
        '/home/aecgroup/aecdata/Results/model_11-Oct-2016_15:41:48_2000000_1_cluster_varDec1e4-1e6/modelAt2000000/model.mat', ...
        '/home/aecgroup/aecdata/Results/model_11-Oct-2016_15:48:21_2000000_1_cluster_varDec1e5-1e5/modelAt2000000/model.mat', ...
        '/home/aecgroup/aecdata/Results/model_11-Oct-2016_15:49:55_2000000_1_cluster_varDec1e5-1e6/modelAt2000000/model.mat', ...
        '/home/aecgroup/aecdata/Results/model_11-Oct-2016_15:51:04_2000000_1_cluster_varDec1e6-1e6/modelAt2000000/model.mat', ...
        % '/home/aecgroup/aecdata/Results/16-10-15_2000000iter_1_regul_1e-2_actorLR_1To0/modelAt2000000/model.mat', ...
        % '/home/aecgroup/aecdata/Results/16-10-15_2000000iter_1_regul_1e-2_actorLR_0.5To0/modelAt2000000/model.mat', ...
        % '/home/aecgroup/aecdata/Results/16-10-15_2000000iter_1_regul_1e-2_actorLR_0.5/modelAt2000000/model.mat', ...
        % '/home/aecgroup/aecdata/Results/16-10-15_2000000iter_1_regul_1e-3_actorLR_1To0/modelAt2000000/model.mat', ...
        % '/home/aecgroup/aecdata/Results/16-10-15_2000000iter_1_regul_1e-3_actorLR_0.5To0/modelAt2000000/model.mat', ...
        % '/home/aecgroup/aecdata/Results/16-10-15_2000000iter_1_regul_1e-3_actorLR_0.5/modelAt2000000/model.mat', ...
        % '/home/aecgroup/aecdata/Results/16-10-15_2000000iter_1_regul_1e-4_actorLR_1To0/modelAt2000000/model.mat', ...
        % '/home/aecgroup/aecdata/Results/16-10-15_2000000iter_1_regul_1e-4_actorLR_0.5To0/modelAt2000000/model.mat', ...
        % '/home/aecgroup/aecdata/Results/16-10-15_2000000iter_1_regul_1e-4_actorLR_0.5/modelAt2000000/model.mat', ...
    };

nFiles = length(files);

%% here, specify the parameter ranges that should be used
varUpper = [1e-4, 1e-5, 1e-6];
varLower = [1e-4, 1e-5, 1e-6];

% varUpper = [1e-4, 1e-3, 1e-2];
% varLower = [[1, 0], [0.5, 0], [0.5]];

length1 = length(varUpper);
length2 = length(varLower);

results = zeros(length1, length2, 3); % results are the three measurements rmse, median, and iqr

for f = 1 : nFiles
    model = load(files{f});
    model = model.model;

    % find indizes
    ind = find(varUpper == model.rlModel.CActor.varianceRange(1));
    jnd = find(varLower == model.rlModel.CActor.varianceRange(2));

    testInterval = model.interval * 2;

    results(ind, jnd, 1) = sqrt(mean(model.testResult3(:, testInterval) .^ 2));
    results(ind, jnd, 2) = iqr(model.testResult3(:, testInterval));
    results(ind, jnd, 3) = median(model.testResult3(:, testInterval));

end

results(: ,:, 1) = results(:, :, 1)';
results(: ,:, 2) = results(:, :, 2)';
results(: ,:, 3) = results(:, :, 3)';

%% plotting section
plotSavePath = './varianceRangeComparison';

numberFormatVar1 = '%1.0e';
numberFormatVar2 = '%1.0e';

labelVar1 = 'variance\nstart value';
labelVar2 = 'variance end value';

% results(find(results == 0)) = inf; % remove zeros?

figure;
title('Parameter Comparison');

% plot the RMSE
subplot(3, 1, 1);
imagesc(results(:, :, 1));
% colormap(jet); %(flipud(autumn)), seems to be the global map for the whole plot

colordata = createCM(3);
colordata(1, :) = [1, 1, 1]; % plot NaNs white
colormap(colordata);

txt = results(:, :, 1);
txt(find(txt == 0)) = inf;
txt = num2str(txt(:),'%0.2f');
% txt(find(txt == '0.00')) = ' '; % this may cause errors.
txt = strtrim(cellstr(txt));
[x, y] = meshgrid(1 : max(length1, length2));
hStrings = text(x(:), y(:), txt(:), 'HorizontalAlignment', 'center');

set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(flipud(varUpper(:)), numberFormatVar1), ...
         'YTick', 1:length2, 'YTickLabel', num2str(varLower(:), numberFormatVar2), ...
         'TickLength', [0, 0]);

title('RMSE');
xlabel(sprintf(labelVar2));
ylabel(sprintf(labelVar1));
colorbar();

% saveas(gca, sprintf('%s_rmse.png', plotSavePath));

% plot the interquartile range
subplot(3, 1, 2);
imagesc(results(:, :, 2));
% colormap(jet); %(flipud(autumn))

txt = results(:, :, 2);
txt(find(txt == 0)) = inf;
txt = num2str(txt(:),'%0.2f');
% txt(find(txt == '0.00')) = ' '; % this may cause errors.
txt = strtrim(cellstr(txt));
[x, y] = meshgrid(1 : max(length1, length2));
hStrings = text(x(:), y(:), txt(:), 'HorizontalAlignment', 'center');

set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(flipud(varUpper(:)), numberFormatVar1), ...
         'YTick', 1:length2, 'YTickLabel', num2str(varLower(:), numberFormatVar2), ...
         'TickLength', [0, 0]);

title('IQR*4');
xlabel(sprintf(labelVar2));
ylabel(sprintf(labelVar1));
colorbar();

% saveas(gca, sprintf('%s_iqr.png', plotSavePath));

% plot the median
subplot(3, 1, 3);
imagesc(abs(results(:, :, 3))); % plot absolute values
% colormap(hsv); %(flipud(autumn))

txt = results(:, :, 3);
txt(find(txt == 0)) = inf;
txt = num2str(txt(:),'%0.2f');
% txt(find(txt == '0.00')) = ' '; % this may cause errors.
txt = strtrim(cellstr(txt));
[x, y] = meshgrid(1 : max(length1, length2));
hStrings = text(x(:), y(:), txt(:), 'HorizontalAlignment', 'center');

set(gca, 'XTick', 1:length1, 'XTickLabel', num2str(flipud(varUpper(:)), numberFormatVar1), ...
         'YTick', 1:length2, 'YTickLabel', num2str(varLower(:), numberFormatVar2), ...
         'TickLength', [0, 0]);

title('Median');
xlabel(sprintf(labelVar2));
ylabel(sprintf(labelVar1));
colorbar();

% single colorbar
% hp4 = get(subplot(3, 1, 3),'Position');
% colorbar('Position', [hp4(1) + hp4(3) + 0.01  hp4(2)  0.1  hp4(2) + hp4(3) * 2.1])

% saveas(gca, sprintf('%s_median.png', plotSavePath));

saveas(gca, sprintf('%s.png', plotSavePath));


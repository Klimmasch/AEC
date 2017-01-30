% parent = '/home/klimmasch/projects/results/Discount Factor vs Interval/'
% files = dir('/home/klimmasch/projects/results/Discount Factor vs Interval/*00*');
% folder = 'CriticLR vs ActorLR';
% folder = 'Discount Factor vs Interval';
% folder = 'Regularizer vs Actor Learning Rate';
% folder = 'Regularizer vs ActorLR';
% folder = 'Regularizer_vs_ActorLR';
% folder = 'exploringMetCost'
% folder = 'Gamma_vs_Interval_fewerResources/'
% folder = 'lambdaMuscleFB_vs_desiredStdZT';
folder = 'lambdaMuscleFB_vs_desiredStdZT_seed2';

parent = strcat('/home/aecgroup/aecdata/Results/', folder);
% parent = strcat('/home/klimmasch/projects/results/', folder);
% files = dir(sprintf('%s/*00iter_1_*', parent));
files = dir(sprintf('%s/*_2_*', parent));
simulator = prepareSimulator([]);

for f = 1:length(files)
    savePath = sprintf('%s/%s', parent, files(f).name)
    % model = load(strcat(savePath, '/modelAt2000000/model.mat'));

    % model = model.model;
    % model.savePath = savePath;
    % if length(model.rlModel.actorLearningRange) == 1
    %     val = model.rlModel.actorLearningRange(1);
    %     model.rlModel.actorLearningRange = [val, val];
    % end
    % display(model.savePath)
%     savePathNew = sprintf('%s/16-10-19_2000000iter_1_cluster_regul_%1.0e_actorLR_[%1.2f-%1.2f]', parent, model.rlModel.CActor.regularizer, model.rlModel.actorLearningRange(1), model.rlModel.actorLearningRange(2))
%     model.savePath = savePathNew;
    % save(strcat(model.savePath, '/model'), 'model');
    nStimTest = 40;
    t = [500000, 1000000, 1500000, 2000000];
    for i = 1 : length(t)
        model = load(strcat(savePath, sprintf('/modelAt%d/model.mat',t(i))));
        model = model.model;
        testModelContinuous(model, nStimTest, 1, 1, 0, simulator, 0, sprintf('modelAt%d', t(i)));
        close all;
    end
end

% fixing missing testAt500000 entry in testPerformanceVsTrainTime plot
% clusterRuns = {'CriticLR vs ActorLR', 'Discount Factor vs Interval', 'Regularizer vs Actor Learning Rate', 'Regularizer vs ActorLR', 'Regularizer_vs_ActorLR'};
% clusterRuns = {'exploringMetCost'};
%
% for k = 1 : length(clusterRuns)
%     sprintf('Working on %s', clusterRuns{k})
%     parent = strcat('/home/aecgroup/aecdata/Results/', clusterRuns{k});
%     subFolders = dir(sprintf('%s/*_1_*', parent));
%
%     for l = 1 : length(subFolders)
%         % load final model instance
%         try
%             load(strcat(parent, '/', subFolders(l).name, '/model.mat'));
%         catch
%             warning('%s does not exist (yet), continue...', strcat(parent, '/', subFolders(l).name, '/model.mat'));
%             continue;
%         end
%
%         % fix index
%         if (model.testHist(2, 1) == 0)
%             model.testHist = [model.testHist(1, :); model.testHist(3 : 5, :); zeros(1, 6)];
%         end
%
%         % load testAt2000000 model
%         try
%             tmpModel = load(strcat(parent, '/', subFolders(l).name, '/modelAt2000000/model.mat'));
%         catch
%             warning('%s does not exist (yet), continue...', strcat(parent, '/', subFolders(l).name, '/modelAt2000000/model.mat'));
%             continue;
%         end
%
%         % include missing testAt2000000 entry
%         model.testHist(5, :) = [sqrt(mean(tmpModel.model.testResult3(:, tmpModel.model.interval * 2) .^ 2)), ...
%                                 mean(abs(tmpModel.model.testResult3(:, tmpModel.model.interval * 2) .^ 2)), ...
%                                 std(abs(tmpModel.model.testResult3(:, tmpModel.model.interval * 2) .^ 2)), ...
%                                 sqrt(mean(tmpModel.model.testResult7(:, tmpModel.model.interval * 2) .^ 2)), ...
%                                 mean(abs(tmpModel.model.testResult7(:, tmpModel.model.interval * 2) .^ 2)), ...
%                                 std(abs(tmpModel.model.testResult7(:, tmpModel.model.interval * 2) .^ 2))];
%
%         % update model's savePath
%         model.savePath = strcat(parent, '/', subFolders(l).name);
%
%         % save updated model
%         save(strcat(model.savePath, '/model'), 'model');
%
%         % replot
%         model.allPlotSave(7);
%         close all;
%     end
% end

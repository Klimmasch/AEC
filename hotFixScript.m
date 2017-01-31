% parent = '/home/klimmasch/projects/results/Discount Factor vs Interval/'
% files = dir('/home/klimmasch/projects/results/Discount Factor vs Interval/*00*');
% folder = 'CriticLR vs ActorLR';
% folder = 'Discount Factor vs Interval';
% folder = 'Regularizer vs Actor Learning Rate';
% folder = 'Regularizer vs ActorLR';
% folder = 'Regularizer_vs_ActorLR';
% folder = 'exploringMetCost'
% folder = 'Gamma_vs_Interval_fewerResources'
% folder = 'varDec_new'
% folder = 'steplength_actorVsWeightInit'
% folder = 'steplength_actorVsVariance_reg1e-5';
% folder = 'steplength_actorVsRegul';
% folder = 'Gamma_vs_Interval_fewerResources';

folders = {'steplength_actorVsRegul', 'steplength_actorVsRegul_1mio', 'steplength_actorVsRegul_reg1e-5', 'steplength_actorVsWeightInit', 'increase_step_width'};

for k = 1 : length(folders)
    
    folder = folders{k};
    
    parent = strcat('/home/aecgroup/aecdata/Results/', folder);
    % parent = strcat('/home/klimmasch/projects/results/', folder);
    files = dir(sprintf('%s/*00iter_1_*', parent));
    % simulator = prepareSimulator([]);

    modelAt = [500000, 1000000, 1500000, 2000000];

    for t = 1 : length(modelAt) 
        for f = 1 : length(files)
            savePath = sprintf('%s/%s', parent, files(f).name)
            filePath = strcat(savePath, sprintf('/modelAt%d/', modelAt(t)))
        %     filePath = strcat(savePath, '/')
            try
                model = load(strcat(filePath, 'model.mat'));
            catch
                sprintf(filePath, ' \n contains no valid file')
                continue
            end
            % model = load(strcat(savePath, '/model.mat'));

            model = model.model;
            model.savePath = savePath;
            % if length(model.rlModel.actorLearningRange) == 1
            %     val = model.rlModel.actorLearningRange(1);
            %     model.rlModel.actorLearningRange = [val, val];
            % end
            % display(model.savePath)
            % save(strcat(model.savePath, '/model'), 'model');

            if ~exist(strcat(filePath, 'metCostsApproach.png'))
                nStimTest = 0;
                try
                    testModelContinuous(model, nStimTest, 1, 1, 1, simulator, 0, sprintf('modelAt%d', modelAt(t)));
                catch 
                    testModelContinuous(model, 40, 1, 1, 1, simulator, 0, sprintf('modelAt%d', modelAt(t)));
    %                 sprintf('%s\n seem to have no training data.', model.savePath)
                end
        %         model.allPlotSave([1:7]);
            else
                sprintf('skipping testing')
            end


        %     model.allPlotSave([7]);
        %     display(model.trainedUntil)
        %     display(size(model.testResult7))
        %     model.allPlotSave([1:3,5:7]); % level == 4 seems to cause problems in some interval>10 cases
            close all;
        end
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

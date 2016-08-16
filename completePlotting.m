% This script takes some files and generated the plots and makes some testing on them
function completePlotting()
    files = {...
        '/home/klimmasch/projects/results/model_12-Aug-2016_20:05:06_3000000_nonhomeo_1_bmsf_newInit_steeperReward/model.mat', ...
        '/home/klimmasch/projects/results/model_12-Aug-2016_19:42:43_3000000_nonhomeo_1_bmsf_newInit_20iters/model.mat', ...
        '/home/klimmasch/projects/results/model_12-Aug-2016_19:46:26_3000000_nonhomeo_1_bmsf_newInit_metCost50pc/model.mat', ...
        '/home/klimmasch/projects/results/model_12-Aug-2016_19:47:48_3000000_nonhomeo_1_bmsf_newInit_metCost25pc/model.mat',...
        '/home/klimmasch/projects/results/model_12-Aug-2016_20:10:36_3000000_nonhomeo_1_bmsf_newInit_newRenderen_newImages/model.mat',...
        '/home/klimmasch/projects/results/model_12-Aug-2016_19:39:30_3000000_nonhomeo_1_bmsf_newInit_lMF1/model.mat', ...
        '/home/klimmasch/projects/results/model_12-Aug-2016_19:36:14_3000000_nonhomeo_1_bmsf_newInit_400BF/model.mat', ...
        '/home/klimmasch/projects/results/model_13-Aug-2016_03:55:38_3000000_nonhomeo_1_bmsf_newInit_newRender_InitNoReset/model.mat', ...
        };

    simulator = OpenEyeSim('create');
%     simulator = OpenEyeSimV4('create');
%     simulator.reinitRenderer();
    simulator.initRenderer();
    
    stimulusRange = 50;
    vergRange = 7;
    objDists = [0.5:0.5:2];
    objRange = length(objDists);
    zeroVergInd = ceil(vergRange/2);
    takeLastValues = 5;
    imageFlag = sprintf('last%dvalues_newScale', takeLastValues);

    meanOffset = zeros(objRange, length(files), stimulusRange * takeLastValues);
    
    for i = 1:length(files)
        model = load(files{i});
        model = model.model;
%         model.allPlotSave;
%         testModel(model, 23, [0.5, 1, 1.5, 2], [-2 : 0.2 : 2], [20, 20], 1, 0, 1, 1);
%         testModel2(model, 50, 1, 1, simulator, 0); % use more stimuli than in the textures file, just in case
        if i == 1
            testModelContinuous(model, stimulusRange, 1, 1, simulator, 0, sprintf('modelAt%d,oldRenderer', model.trainedUntil));
        else
            testModelContinuous(model, stimulusRange, 1, 1, simulator, 1, sprintf('modelAt%d,oldRenderer', model.trainedUntil));
        end
        sprintf('###### plotting and testing completed in %s #######', files{i});
        close all;
        
%         nSamples = size(model.testResult3, 1);
%         if nSamples == 0
%             simulator = OpenEyeSim('create');
%             simulator.initRenderer();
%             simulator.reinitRenderer();
%             testModel2(model, 50, 0, 1, simulator, 0);
%         end
%         nSamples = size(model.testResult3, 1)
%         if nSamples == 0
%             testModel2(model, 50, 0, 1, simulator, 0);
%         end
        % extract all values that are generated with zero vergErr init
%         tmpResult = []; % these are the values of the 33 input stimuli that started with 0 vergErr at the end of 10 iterations
%         for j = 1:nSamples


%             if mod(j, vergRange) == zeroVergInd % if we started at zero vergence Error
%                 % append the last value(s)
%                 for v = 1:takeLastValues
%                     tmpResult = [tmpResult, model.testResult3(j, end - v + 1)];
%                 end
%             end
% 
%             if length(tmpResult) == stimulusRange * takeLastValues
%                 objInd = ceil((j * objRange) / nSamples);
%                 meanOffset(objInd, i, :) = tmpResult;
%                 tmpResult = [];
%             end
%         end
    end
    %% todo: remove outliers from the plot and scale vergence error axis to the same range --> this is now part of the testmodelContinuous
%     for k = 1:objRange
%         figure; hold on; grid on; grid minor;
%         [~, tmpSize1, tmpSize2] = size(meanOffset(k, :, :));
%         boxplot(reshape(meanOffset(k, :, :), [tmpSize1, tmpSize2])', 'labels', [0 1 5 10 30 50]);
%         title(sprintf('Offset from zero vergence error at %.1f m', objDists(k)));
%         xlabel('lambda_{MF} in %');
%         ylabel('Vergence Error in deg');
%         ylim([-1, 2]);
%         saveas(gcf, sprintf('../results/biasesAt%.2f_%s.png', objDists(k), imageFlag), 'png');
% 
%     end
end

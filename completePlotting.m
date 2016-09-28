% This script takes some files and generated the plots and makes some testing on them
function completePlotting()
    files = {...
        '/home/klimmasch/projects/results/model_17-Sep-2016_09:55:47_5000000_1_bmsf_advancedInit_metCost25/modelAt5000000/model.mat', ...
        '/home/klimmasch/projects/results/model_20-Sep-2016_17:50:49_5000000_1_bmsf_advancedInit_metCost0_BFActivation5-5/modelAt5000000/model.mat', ...
        '/home/klimmasch/projects/results/model_17-Sep-2016_09:32:14_5000000_1_bmsf_advancedInit_metCost0/modelAt5000000/model.mat', ...
        '/home/klimmasch/projects/results/model_13-Sep-2016_17:43:39_6000000_1_bmsf_advancedInit_metCost1/modelAt5500000/model.mat',...
        '/home/klimmasch/projects/results/model_13-Sep-2016_17:43:39_6000000_1_bmsf_advancedInit_metCost1/modelAt5000000/model.mat',...'
        };

    % simulator = OpenEyeSim('create');
    % simulator = OpenEyeSimV4('create');
    % simulator.reinitRenderer();
    % simulator.initRenderer();
    simulator = prepareSimulator();

    stimulusRange = 40;
    vergRange = 7;
    objDists = [0.5 1 : 6];
    objRange = length(objDists);
    zeroVergInd = ceil(vergRange/2);
    takeLastValues = 5;
    imageFlag = sprintf('last%dvalues', takeLastValues);

    meanOffset = zeros(objRange, length(files), stimulusRange * takeLastValues);
    
    for i = 1:length(files)
        model = load(files{i});
        model = model.model;
        model.allPlotSave(6);
        testModelContinuous(model, stimulusRange, 1, 1, simulator, 0, sprintf('modelAt%d', model.trainedUntil));

        sprintf('###### plotting and testing completed in %s #######', files{i});
        close all;
    end
end

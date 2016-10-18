% This script takes some files and generated the plots and makes some testing on them
function completePlotting(simulator)
    files = {...
        '/home/klimmasch/projects/results/model_20-Sep-2016_17:51:44_5000000_1_bmsf_advancedInit_metCost0_criticLR075/modelAt1000000/model.mat', ...
        '/home/klimmasch/projects/results/model_20-Sep-2016_17:51:44_5000000_1_bmsf_advancedInit_metCost0_criticLR075/modelAt2000000/model.mat', ...
        '/home/klimmasch/projects/results/model_20-Sep-2016_17:51:44_5000000_1_bmsf_advancedInit_metCost0_criticLR075/modelAt3000000/model.mat', ...
        '/home/klimmasch/projects/results/model_20-Sep-2016_17:51:44_5000000_1_bmsf_advancedInit_metCost0_criticLR075/modelAt4000000/model.mat', ...
        '/home/klimmasch/projects/results/model_20-Sep-2016_17:51:44_5000000_1_bmsf_advancedInit_metCost0_criticLR075/modelAt5000000/model.mat', ...
        };

    % simulator = OpenEyeSim('create');
    % simulator = OpenEyeSimV4('create');
    % simulator.reinitRenderer();
    % simulator.initRenderer();
    if isempty(simulator)
        simulator = prepareSimulator();
    end

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
%         model.allPlotSave(6);
        testModelContinuous(model, stimulusRange, 1, 1, 0, simulator, 0, sprintf('modelAt%d', model.trainedUntil));

        sprintf('###### plotting and testing completed in %s #######', files{i});
        close all;
    end
end

% This script takes some files and generated the plots and makes some testing on them
function completePlotting()
    files = {...
        '/model.mat', ...
        '/model.mat', ...
        '/model.mat', ...
        '/model.mat',...
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
        % model.allPlotSave;
        testModelContinuous(model, stimulusRange, 1, 1, simulator, 0, sprintf('modelAt%d,oldRenderer', model.trainedUntil));

        sprintf('###### plotting and testing completed in %s #######', files{i});
        close all;
    end
end

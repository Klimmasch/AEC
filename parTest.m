

function parTest(folders, nWorkers, runParallel)

%     folders = {'GammaVsMetCosts_0,5mio'}; % cluster run tier
%     folders = {'CritLRVsMetCosts_1mio'};

    nStim = 40;
    testpaths = {};
    iter = 1;
    
    for k = 1 : length(folders) % single run tier
        folder = folders{k};
        parent = strcat('/home/aecgroup/aecdata/Results/', folder);
        subfolders = dir(sprintf('%s/*iter*', parent));
        
        for s = 1 : length(subfolders) % test folder tier
            subfolder = subfolders(s)

%             if subfolder.isdir
%                 testpaths{iter} = sprintf('%s/%s/model.mat', parent, subfolder.name);
%                 iter = iter + 1;
%             end
            testfolders = dir(sprintf('%s/%s/*modelAt*',parent, subfolder.name));
            
            for t = 1 : length(testfolders)
                testfolder = testfolders(t);
                testpaths{iter} = sprintf('%s/%s/%s/model.mat', parent, subfolder.name, testfolder.name);
                iter = iter + 1;
            end
        end
    end

%     testpaths% = fliplr(testpaths); % change order to run on different nodes
    
    if runParallel
        myCluster = parcluster('local');
        myCluster.NumWorkers = nWorkers;
        % saveProfile(myCluster);

        if (isempty(gcp('nocreate')))
            mPool = parpool(nWorkers);
        end
   
        parfor tp = 1 : length(testpaths)
            savePath = testpaths{tp}
            try
                model = load(strcat(savePath));
            catch
                sprintf('no valid file found')
                continue
            end
            model = model.model;
            simulator = [];  % lets see if it works with the parallel pool like this
            % if ~exist(sprintf('%s/modelAt%d/muscleActivityTrajectory.png', model.savePath, model.trainedUntil), 'file') % if the last image from the test procedure does not exists ...
            % if ~any(model.testHist(find(model.testAt == model.trainedUntil))) % test if according field in testHist is empty
            if isempty(model.testResult7)    
                sprintf('could not find test results in %s\n starting test procedure.', savePath)
                testModelContinuous(model, nStim, 1, 1, 1, simulator, 0, sprintf('modelAt%d', model.trainedUntil), [1, 3 : 6]);
%                 testModelContinuous(model, 0, 1, 1, 1, simulator, 0, sprintf('modelAt%d', model.trainedUntil), [1, 3 : 4, 6]);
%                 testModelContinuous(model, 0, 0, 1, 1, simulator, 0, sprintf('modelAt%d', model.trainedUntil), [1, 3, 4, 6]); % only save testHist
            else
                sprintf('skip testing for\n%s', savePath)
            end
        end
        
    else
        textureFile = 'Textures_mcgillManMade40.mat';
        simulator = OpenEyeSimV5('create'); % latest renderer

        simulator.initRenderer();

        texture = load(sprintf('config/%s', textureFile));
        texture = texture.texture;
        for i = 1 : nStim
            simulator.add_texture(i, texture{i});
        end
        sprintf('%d textures added to the testing simulator', nStim)
        
        for tp = 1 : length(testpaths)
            savePath = testpaths{tp}
            try
                model = load(strcat(savePath));
            catch 
                sprintf('no valid file found')
                continue
            end
            model = model.model;
            
            % nur fuer alex:
            if model.trainedUntil ~= 500000
                continue
            end
            
            if model.trainedUntil == model.trainTime
                model.allPlotSave(1:7);
                close all;
                sprintf('finisched plotting')
            else
                sprintf('skip testing')
            end
            % sprintf('trainedUntil: %d', model.trainedUntil)
            % if ~exist(sprintf('%s/modelAt%d/muscleActivityTrajectory.png', model.savePath, model.trainedUntil), 'file') % if the last image from the test procedure does not exists ...
            % if ~any(model.testHist(find(model.testAt == model.trainedUntil))) % test if according field in testHist is empty
%             if isempty(model.testResult7)
%                 sprintf('could not find test results in %s\n starting test procedure.', savePath)
%                 testModelContinuous(model, nStim, 1, 1, 2, simulator, 0, sprintf('modelAt%d', model.trainedUntil), [1, 3 : 6]);
%                 close all;
%             else
%                 sprintf('skip testing for\n%s', savePath)
%             end
        end
    end
end



function parTest(nWorkers)

    folders = {'GammaVsMetCosts_0,5mio'}; % cluster run tier
    testpaths = {};
    iter = 1;
    
    for k = 1 : length(folders) % single run tier
        folder = folders{k};
        parent = strcat('/home/aecgroup/aecdata/Results/', folder);
        subfolders = dir(sprintf('%s/*iter*', parent));
        
        for s = 1 : length(subfolders) % test folder tier
            subfolder = subfolders(s);
            testfolders = dir(sprintf('%s/%s/*modelAt*',parent, subfolder.name));
            
            for t = 1 : length(testfolders)
                testfolder = testfolders(t);
                testpaths{iter} = sprintf('%s/%s/%s/model.mat', parent, subfolder.name, testfolder.name);
                iter = iter + 1;
            end
        end
    end
    
    myCluster = parcluster('local');
    myCluster.NumWorkers = nWorkers;
    % saveProfile(myCluster);

    if (isempty(gcp('nocreate')))
        mPool = parpool(nWorkers);
    end
                
    parfor tp = 1 : length(testpaths)
        savePath = testpaths{tp}
        model = load(strcat(savePath));
        model = model.model;

        simulator = []; % lets see if it works with the parallel pool like this
        testModelContinuous(model, 40, 1, 1, 0, simulator, 0, sprintf('modelAt%d', model.trainedUntil), [1 : 6]);
    end
        
    end
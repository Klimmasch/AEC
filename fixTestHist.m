% this script fixes a bug in the saving of testing history

function fixTestHist(folders)

    if isempty(folders)
%     folders = {'GammaVsMetCosts_0,5mio'}; % cluster run tier
%     folders = {'CritLRVsMetCosts_1mio'};
    end

    testpaths = {};
    iter = 1;
    
    for k = 1 : length(folders) % single run tier
        folder = folders{k};
        parent = strcat('/home/aecgroup/aecdata/Results/', folder);
        subfolders = dir(sprintf('%s/*iter*', parent));
        
        for s = 1 : length(subfolders) % test folder tier
            subfolder = subfolders(s);

            if subfolder.isdir
                testpaths{iter} = sprintf('%s/%s/model.mat', parent, subfolder.name);
                iter = iter + 1;
            end
%             testfolders = dir(sprintf('%s/%s/*modelAt*',parent, subfolder.name));
            
%             for t = 1 : length(testfolders)
%                 testfolder = testfolders(t);
%                 testpaths{iter} = sprintf('%s/%s/%s/model.mat', parent, subfolder.name, testfolder.name);
%                 iter = iter + 1;
%             end
        end
    end

%     testpaths% = fliplr(testpaths); % change order to run on different nodes
    

%         textureFile = 'Textures_mcgillManMade40.mat';
%         simulator = OpenEyeSimV5('create'); % latest renderer
% 
%         simulator.initRenderer();
% 
%         texture = load(sprintf('config/%s', textureFile));
%         texture = texture.texture;
%         for i = 1 : nStim
%             simulator.add_texture(i, texture{i});
%         end
%         sprintf('%d textures added to the testing simulator', nStim)
        
        for tp = 1 : length(testpaths)
            savePath = testpaths{tp}
            try
                model = load(strcat(savePath));
            catch 
                sprintf('no valid file found')
                continue
            end
            model = model.model;
            
            % case testhist contains a 0
            if size(model.testHist, 1) == length(model.testAt)
                % new standard case
                if model.testAt(1) ~= 0
                    obj.testAt = horzcat(0, obj.testAt);
                end
                    
                if model.normFeatVect == 0
                    model.testHist(1, :) = [1.1593, 1.3440, 1.5243, 1.0736, 1.1527, 0.9517];
                    sprintf('no normalized feature vector: warning: these are the old values in model.testhist(1, :)')
                else
%                     model.testHist(1, :) = [1.1557, 1.3355, 1.5078, 1.0812, 1.1689, 0.9552]; % for the rmse
                    model.testHist(1, :) = [0.0028, 0.0000, 0.0085, 0.0007, 0.0018, 0.0010];
                end

                for testInd = 1 : length(model.testAt)
                    if testInd > 1
                        tmpFilePath = sprintf('%s/modelAt%d/model.mat', model.savePath, model.testAt(testInd));
                        try
                            tempModel = load(tmpFilePath);
                            tempModel = tempModel.model;
                        catch
                            sprintf('%s\ndoes not contain a valid file', tmpFilePath)
                            continue
                        end

                        model.testHist(testInd, :) = ...
                                  [mean(tempModel.vergenceAngleApproach(:, tempModel.testInterval)), ...
                                   median(tempModel.vergenceAngleApproach(:, tempModel.testInterval)), ...
                                   iqr(tempModel.vergenceAngleApproach(:, tempModel.testInterval)), ...
                                   mean(tempModel.metCostsApproach(:, tempModel.testInterval)), ...
                                   median(tempModel.metCostsApproach(:, tempModel.testInterval)), ...
                                   iqr(tempModel.metCostsApproach(:, tempModel.testInterval))];
                    end
                end 
                
                model.allPlotSave(7);
                save(strcat(model.savePath, '/model'), 'model');
                
            else
                sprintf('length(model.testHist, 1) ~= length(model.testAt) !') 
            end
          
%             if model.trainedUntil == model.trainTime
%                 model.allPlotSave(1:7);
%                 close all;
%                 sprintf('finisched plotting')
%             else
%                 sprintf('skip testing')
%             end
%             testModelContinuous(model, 0, 1, 0, 1, [], 0, sprintf('modelAt%d', model.trainedUntil), [1, 3, 4, 6]);
            % sprintf('trainedUntil: %d', model.trainedUntil)
            % if ~exist(sprintf('%s/modelAt%d/muscleActivityTrajectory.png', model.savePath, model.trainedUntil), 'file') % if the last image from the test procedure does not exists ...
%             if ~any(model.testHist(find(model.testAt == model.trainedUntil))) % test if according field in testHist is empty
%             if isempty(model.testResult7)
%                 sprintf('could not find test results in %s\n starting test procedure.', savePath)
%                 testModelContinuous(model, nStim, 1, 1, 2, simulator, 0, sprintf('modelAt%d', model.trainedUntil), [1, 3 : 6]);
%                 close all;
%             else
%                 sprintf('skip testing for\n%s', savePath)
%             end
            close all
        end
end

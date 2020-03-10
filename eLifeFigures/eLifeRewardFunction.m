% Plotting the reward function for the different reared models

function eLifeRewardFunction(scale, generateNewData, simulator, saveTag)
    filePath = '/home/aecgroup/aecdata/Results/'; %'/SAB2018/' or 'eLifePaper'
%     filePath = '/home/mrklim/fiasAECGroup/aecgroup/aecdata/Results/'
    savePath = '/home/aecgroup/aecdata/Results/eLifePaper/plots/Rewards/'
%     savePath = '/home/mrklim/fiasAECGroup/aecdata/Results/eLifePaper/plots/Rewards/'

    plottingScheme = 'se'
%     plottingScheme = 'simpleStd'
    
    objDists = [0.5, 3, 6]%[0.5, 1 : 6]
    dispRange = linspace(-2,2,11);
    nStim = 10%40
    colrs = {[27,158,119]./256, [ 217,95,2]./256}; %green and orange for coarse and fine scale, color-blind safe
    
    models = {
        {   % normal case
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1'},
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2'},
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3'},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4'},
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5'},
        },{ % vertical only (filter size 33 px)
            {'eLifePaper/explFilterSizes/19-05-23_500000iter_1_fsize6std_filtBoth_30_prob_1_seed1'},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_30_prob_1_seed2'},
            {'eLifePaper/explFilterSizes/18-10-01_500000iter_3_fsize6std_filtBoth_30_prob_1_seed3'},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_30_prob_1_seed4'},
            {'eLifePaper/explFilterSizes/19-05-24_500000iter_5_fsize6std_filtBoth_30_prob_1_seed5'},
            % vertical only (filter size 100 px, filter name: 31)
%             {'eLifePaper/explFilterSizes/'},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_31_prob_1_seed2'},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_31_prob_1_seed3'},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_31_prob_1_seed4'},
%             {'eLifePaper/explFilterSizes/'},
        },{ % horizontal only (33 px)
            {'eLifePaper/explFilterSizes/19-05-23_500000iter_1_fsize6std_filtBoth_34_prob_1_seed1'},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_34_prob_1_seed2'},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_3_fsize6std_filtBoth_34_prob_1_seed3'},
            {'eLifePaper/explFilterSizes/18-10-01_500000iter_4_fsize6std_filtBoth_34_prob_1_seed4'},
            {'eLifePaper/explFilterSizes/19-05-24_500000iter_5_fsize6std_filtBoth_34_prob_1_seed5'},
            % horizontal only (100 px, filter name: 35)
%             {'eLifePaper/explFilterSizes/'},
%             {'eLifePaper/explFilterSizes/18-09-29_500000iter_2_fsize6std_filtBoth_35_prob_1_seed2'},
%             {'eLifePaper/explFilterSizes/18-09-29_500000iter_3_fsize6std_filtBoth_35_prob_1_seed3'},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_4_fsize6std_filtBoth_35_prob_1_seed4'},
%             {'eLifePaper/explFilterSizes/'},
        },{ % orthogonal
            {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_46_prob_1_seed1'},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_2_fsize6std_filtBoth_46_prob_1_seed2'},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_3_fsize6std_filtBoth_46_prob_1_seed3'},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_4_fsize6std_filtBoth_46_prob_1_seed4'},
            {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_46_prob_1_seed5'},
        },{ 
            % monocular blur
            {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_45_prob_1_seed1'},
            {'eLifePaper/explFilterSizes/18-10-18_500000iter_2_fsize6std_filtBoth_45_prob_1_seed2'},
            {'eLifePaper/explFilterSizes/18-10-19_500000iter_3_fsize6std_filtBoth_45_prob_1_seed3'},
            {'eLifePaper/explFilterSizes/18-10-18_500000iter_4_fsize6std_filtBoth_45_prob_1_seed4'},
            {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_45_prob_1_seed5'},
            % },{ % complete monocular deprivation ?
            %     {'eLifePaper/explFilterSizes/', scale, 1},
            %     {'eLifePaper/explFilterSizes/', scale, 1},
            %     {'eLifePaper/explFilterSizes/', scale, 1},
        },{  % strabismic case (5 degree)
% %             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_5_seed2'},
% %             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_5_seed3'},
% %             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_5_seed4'},
%             % 10 degree strabism angle, vergErr and objDist fixed at 6m
%             {'eLifePaper/strabism/19-05-29_500000iter_1_fixAllAt6m_filtB_29_strabAngle_10_seed1'},
%             {'eLifePaper/strabism/19-02-18_500000iter_2_fixAllAt6m_filtB_29_strabAngle_10_seed2'},
%             {'eLifePaper/strabism/19-02-18_500000iter_3_fixAllAt6m_filtB_29_strabAngle_10_seed3'},
%             {'eLifePaper/strabism/19-02-18_500000iter_4_fixAllAt6m_filtB_29_strabAngle_10_seed4'},
%             {'eLifePaper/strabism/19-05-29_500000iter_5_fixAllAt6m_filtB_29_strabAngle_10_seed5'},
% %             % non-fixed version
%             {'eLifePaper/strabism/'},
            {'eLifePaper/strabism/19-06-26_500000iter_1_learingActivefiltB_29_strabAngle_10_seed1'},
            {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2'},
            {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3'},
            {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4'},
            {'eLifePaper/strabism/19-06-26_500000iter_5_learingActivefiltB_29_strabAngle_10_seed5'},
%             {'eLifePaper/strabism/'},
%             % 10 degree strabism angle + fix obj.Dist and vergAngle
% %             {'eLifePaper/strabism/19-02-18_500000iter_2_fixAllAt6m_filtB_29_strabAngle_10_seed2'},
% %             {'eLifePaper/strabism/19-02-18_500000iter_3_fixAllAt6m_filtB_29_strabAngle_10_seed3'},
% %             {'eLifePaper/strabism/19-02-18_500000iter_4_fixAllAt6m_filtB_29_strabAngle_10_seed4'},
        }
    };
    names = {'normal', 'vertical', 'horizontal', 'orthogonal', 'monocular', 'strabismic'};
    
    if generateNewData
        
        if isempty(simulator)
            simulator = prepareSimulator({'Textures_mcgillManMade40.mat'});
        end
        
        allRecs = {};   % for single models per condition
        allCrits = {};
        
        recMeans = {};  % for averaging over multiple models per condition
        recSEs = {};
        
        for m = 1:length(models)
            
            recErrs = zeros(length(models{m}), length(objDists), length(dispRange), nStim, 2);
            critResp = zeros(length(models{m}), length(objDists), length(dispRange), nStim);
            
            for n = 1:length(models{m})

                model = load(strcat(filePath, string(models{m}{n}), '/model.mat'));
                model = model.model;
                nScales = length(model.scModel);

                for od = 1:length(objDists)
                    for disp = 1:length(dispRange)
                        for stim = 1:nStim
                            angleDes = 2 * atand(model.baseline / (2 * objDists(od)));   % desired vergence [deg]
    %                         angleNew = angleDes + dispRange(disp);
                            [command, angleNew] = model.getMF2(objDists(od), dispRange(disp));
    %                         [command, angleNew] = model.getMFedood(objDists(od), dispRange(disp));

                            if isempty(model.objSize)
                                model.objSize = 3; %default value
                            end
                            
                            %% for applying corrections:
%                             model.strabAngle = 0;
%                             model.filterLeft = [];
%                             model.filterRight = [];
                            
                            
                            model.refreshImagesNew(simulator, stim, angleNew / 2, objDists(od), model.objSize, [0, 0, 0]);

                            % change left and right images to simulate altered rearing conditions
    %                         randForLeftFilt = rand(1,1);
    %                         randForRightFilt = rand(1,1);

                            if ~isempty(model.filterLeft)
    %                             if randForLeftFilt < model.filterLeftProb
                                    model.imgGrayLeft = conv2(model.filterLeft{1}, model.filterLeft{2}, model.imgGrayLeft, 'same');
    %                             end
                            end
                            if ~isempty(model.filterRight)
    %                             if randForRightFilt < model.filterRightProb
                                    model.imgGrayRight = conv2(model.filterRight{1}, model.filterRight{2}, model.imgGrayRight, 'same');
                %                     model.imgGrayRight = ones(size(model.imgGrayRight)).*mean(mean(model.imgGrayRight));
    %                             end
                            end


                            % imshow(imfuse(model.imgGrayLeft, model.imgGrayRight))

                            % Image patch generation
                            for i = 1 : length(model.scModel)
                                model.preprocessImage(i, 1);
                                model.preprocessImage(i, 2);
                                currentView{i} = vertcat(model.patchesLeft{i}, model.patchesRight{i});
                            end

                            % Generate basis function feature vector from current images
                            [bfFeature, reward, recErrorArray] = model.generateFR(currentView);


                            for s = 1:length(model.scModel)
                                recErrs(n, od, disp, stim, s) = recErrorArray(s); % from coarse to fine
                            end

                            if (model.normFeatVect == 0)
                                feature = [bfFeature; command * model.lambdaMuscleFB];
                            else
                                %% Normalized feature vector
                                % z-transform raw feature vector (no muscle feedback scaling)
                                feature = [bfFeature; command];
                                for i = 1 : length(feature)
                                    feature(i) = model.onlineNormalize(model.trainedUntil, feature(i), i, 1);
                                end
                                feature = [feature(1 : end - 2); feature(end - 1 : end) * model.lambdaMuscleFB];
                            end

                            % bias analysis
                            if (model.rlModel.bias > 0)
                                feature = [feature; model.rlModel.bias];
                            end

    %                         %%% Calculate metabolic costs
    %                         metCost = model.getMetCost(command) * 2;
    % 
    %                         %%% Calculate reward function
    %                         % Standard reward
    %                         rewardFunction = model.lambdaRec * reward - model.lambdaMet * metCost;

    %                         for s = 1:length(model.scModel)
    %                             model.scModel{s}.CCritic.calculate(feature)
    %                             critResp(objDists(od), dispRange(disp), stim, s) = model.scModel{s}.CCritic.value;
    %                         end
                            model.rlModel.CCritic.calculate(feature);
                            critResp(n, od, disp, stim) = model.rlModel.CCritic.value;
                        end
                    end
                end
                allRecs{m} = recErrs;
                allCrits{m} = critResp;
            end
            
            means = zeros(length(dispRange),1);
            stds = zeros(length(dispRange),1);
            for d = 1:length(dispRange)
                dat = reshape(recErrs(:,:,d,:,scale),length(objDists)*nStim*length(models{m}), 1);
                means(d) = mean(dat);
%                 stds(d) = std(dat);
                stds(d) = std(dat) / sqrt(length(dat));
            end
                
            recMeans{m} = means;
%             recSEs{m} = std(cell2mat(conditionData)')/sqrt(length(models{m}));
            recSEs{m} = stds;
        end
        
        if ~isempty(saveTag)
            save(strcat(savePath,'rewardData_',saveTag,'.mat'), 'allRecs', 'allCrits', 'recMeans', 'recSEs');
        end
    else
        %% load saved data
        load(strcat(savePath,'rewardData_',saveTag))
    end
    
    %% Plotting section: reconstruction error for coarse and fine scale
    colors = [[77,175,74]; [55,126,184]; [152,78,163]; [255,127,0]; [228,26,28]; [166,86,40]; [118,42,131]]./256;
%     scale = 1 % scale that will be plotted
    
%     figure('position', [300 350 1000 625]); % single subplots
    figure('position', [300 350 780 500]);
    plotHandles = {};
    for m = 1:length(models)
%         subplot(2,3,m);
        recs = allRecs{m};
        recsC = zeros(length(dispRange), 1);
        stdC = zeros(length(dispRange), 1);
        recsF = zeros(length(dispRange), 1);
        stdF = zeros(length(dispRange), 1);
        for d = 1:length(dispRange)
            coarseData = reshape(recs(:,:,d,:,1),1,[]);
            recsC(d) = mean(coarseData);
            stdC(d) = std(coarseData);
            stdC(d) = std(coarseData);
%             stdC(d) = std(coarseData)/sqrt(length(coarseData)); % se looks very nice!
            
            fineData = reshape(recs(:,:,d,:,2),1,[]);
            recsF(d) = mean(fineData);
            stdF(d) = std(fineData);
%             stdF(d) = std(fineData)/sqrt(length(fineData));
        end
        hold on;
% %         plot(dispRange, recsC, 'Color', colrs{1});
%         errorbar(dispRange, recsC, stdC,'Color', colrs{1}) % replace with boundedline?
% %         plot(dispRange, recsF, 'Color', colrs{2});
%         errorbar(dispRange, recsF, stdF, 'Color', colrs{2})
        
        % plotting mean and std over one simulation per condition
        if strcmp(plottingScheme, 'simpleStd')
            % plotting mean and std over one simulation per condition
            if scale == 1
                [hl, hp] = boundedline(dispRange, recsC, stdC, 'alpha');
            elseif scale == 2
                [hl, hp] = boundedline(dispRange, recsF, stdF, 'alpha');
            end
            hl.LineWidth = 1;
        elseif strcmp(plottingScheme, 'se')
            % plotting mean and se over multiple simulations
            [hl, hp] = boundedline(dispRange, recMeans{m}, recSEs{m}, 'alpha');
            hl.LineWidth = 2;
        end
        hl.Color = colors(m, :);%./255;
        hp.FaceColor = hl.Color;
        hp.FaceAlpha = 0.3;
        plotHandles{m} = hl;
        
%         title(names{m});
%         ylim([0, 0.3]);
%         if m == 3
%             legend('coarse scale', 'fine scales')
%         end
%         if m == 5
%             xlabel('disparity [deg]');
%         end
%         if (m == 1) || (m == 4)
%             ylabel('reconstruction error [au]')
%         end
    end
    
    xlabel('disparity [deg]');
    ylabel('reconstruction error [au]')
    
    legend({'1', names{1}, '2', names{2}, '3', names{3}, '4', names{4}, '5', names{5}, '6', names{6}}, 'Location', 'east');
%     legend({'1', names{1}, '2', names{2}, '3', names{3}, '4', names{4}}, 'Location', 'east');
    [~,hObj] = legend([plotHandles{:}]);
    hL=findobj(hObj, 'type', 'line');
    set(hL, 'linewidth', 3);
    
    if strcmp(plottingScheme, 'simpleStd')
        set(gca, 'YScale', 'log');
        ylim([0, 0.3]);
    elseif strcmp(plottingScheme, 'se')
        set(gca, 'YScale', 'log');
        ylim([0, 0.1]);
    end
    
    if ~isempty(saveTag)
        saveas(gcf, sprintf('%sRecErrs_%s_scale%d.png', savePath, saveTag, scale));
%         ylim([0,0.008])
%         xlim([-0.5,0.5])
%         saveas(gcf, sprintf('%sRecErrs_%s_scale%d_zoomedIn.png', savePath, saveTag, scale));
    end
    
%     %% Plotting section: citics response
%     figure('position', [300 350 1000 625]); % maybe update size
%     for m = 1:length(models)
% %         subplot(2,3,m);
%         crits = allCrits{m};
%         critsC = zeros(length(dispRange), 1);
%         stdC = zeros(length(dispRange), 1);
%         for d = 1:length(dispRange)
%             critsC(d) = mean(reshape(crits(:,d,:,1),1,[]));
%             stdC(d) = std(reshape(crits(:,d,:,1),1,[]));
%         end
%         hold on;
% %         plot(dispRange, recsC, 'Color', colrs{1});
% %         errorbar(dispRange, critsC, stdC, 'Color', colrs{1}) % replace with boundedline?
% %         plot(dispRange, recsF, 'Color', colrs{2});
% %         errorbar(dispRange, critsF, stdF, 'Color', colrs{2})
% 
%         [hl, hp] = boundedline(dispRange, critsC, stdC, 'alpha');
%         hl.Color = colors(m, :);%./255;
%         hp.FaceColor = hl.Color;
%         plotHandles{m} = hl;
% 
% %         title(names{m});
% %         ylim([-0.5, 0]);
% %         if m == 3
% %             legend('coarse scale', 'fine scales')
% %         end
% %         if m == 5
% %             xlabel('disparity [deg]');
% %         end
% %         if (m == 1) || (m == 4)
% %             ylabel("critic's response [au]")
% %         end
%     end
%     xlabel('disparity [deg]');
%     ylabel("critic's response [au]")
%     ylim([-0.5, 0]);
%     
%     legend({'1', names{1}, '2', names{2}, '3', names{3}, '4', names{4}, '5', names{5}, '6', names{6}}, 'Location', 'east');
% %     legend({'1', names{1}, '2', names{2}, '3', names{3}, '4', names{4}}, 'Location', 'east');
%     [~,hObj] = legend([plotHandles{:}]);
%     hL=findobj(hObj, 'type', 'line');
%     set(hL, 'linewidth', 3);
%     
%     if ~isempty(saveTag)
%         saveas(gcf, sprintf('%sCritResp_%s.png', savePath, saveTag));
%     end

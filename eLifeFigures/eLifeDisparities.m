%% This script gathers the data from all accoding simulations and generates
%% an overview of the edge representation

%% depicted are averaged results from normal, vertical, horizontal, orthogonal,
%% strabismic and monocular rearing
function eLifeDisparities(scale, saveTag)
    filePath = '/home/aecgroup/aecdata/Results/'; %'/SAB2018/' or 'eLifePaper'
    savePath = '/home/aecgroup/aecdata/Results/eLifePaper/plots/Disparity/'

%     scale = 1; % now try to combine both scales in one plot
    eye = 1; % take binocular fits
    
    models = {
%         {   % normal case
%             {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1', scale, eye},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2', scale, eye},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4', scale, eye},
%             {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5', scale, eye},
% % %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2', 2, eye},
% % %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3', 2, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4', 2, eye},
% % 
%         },{ % vertical only (filter size 33 px)
%             {'eLifePaper/explFilterSizes/19-05-23_500000iter_1_fsize6std_filtBoth_30_prob_1_seed1', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_30_prob_1_seed2', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-01_500000iter_3_fsize6std_filtBoth_30_prob_1_seed3', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_30_prob_1_seed4', scale, eye},
%             {'eLifePaper/explFilterSizes/19-05-24_500000iter_5_fsize6std_filtBoth_30_prob_1_seed5', scale, eye},
% %             {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_30_prob_1_seed2', 2, eye},
% %             {'eLifePaper/explFilterSizes/18-10-01_500000iter_3_fsize6std_filtBoth_30_prob_1_seed3', 2, eye},
% %             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_30_prob_1_seed4', 2, eye},
%         },{ % horizontal only (33 px)
%             {'eLifePaper/explFilterSizes/19-05-23_500000iter_1_fsize6std_filtBoth_34_prob_1_seed1', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_34_prob_1_seed2', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_3_fsize6std_filtBoth_34_prob_1_seed3', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-01_500000iter_4_fsize6std_filtBoth_34_prob_1_seed4', scale, eye},
%             {'eLifePaper/explFilterSizes/19-05-24_500000iter_5_fsize6std_filtBoth_34_prob_1_seed5', scale, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_34_prob_1_seed2', 2, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-02_500000iter_3_fsize6std_filtBoth_34_prob_1_seed3', 2, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-01_500000iter_4_fsize6std_filtBoth_34_prob_1_seed4', 2, eye},
% % 
%         },{ 
% %             % orthogonal left eye
%             {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_46_prob_1_seed1', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-22_500000iter_2_fsize6std_filtBoth_46_prob_1_seed2', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-22_500000iter_3_fsize6std_filtBoth_46_prob_1_seed3', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-22_500000iter_4_fsize6std_filtBoth_46_prob_1_seed4', scale, eye},
%             {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_46_prob_1_seed5', scale, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-22_500000iter_2_fsize6std_filtBoth_46_prob_1_seed2', 2, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-22_500000iter_3_fsize6std_filtBoth_46_prob_1_seed3', 2, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-22_500000iter_4_fsize6std_filtBoth_46_prob_1_seed4', 2, eye},
% % %         },{ % orthogonal right eye
% % %             {'eLifePaper/explFilterSizes/18-10-22_500000iter_2_fsize6std_filtBoth_46_prob_1_seed2', scale, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-22_500000iter_3_fsize6std_filtBoth_46_prob_1_seed3', scale, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-22_500000iter_4_fsize6std_filtBoth_46_prob_1_seed4', scale, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-22_500000iter_2_fsize6std_filtBoth_46_prob_1_seed2', 2, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-22_500000iter_3_fsize6std_filtBoth_46_prob_1_seed3', 2, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-22_500000iter_4_fsize6std_filtBoth_46_prob_1_seed4', 2, eye},
%         },{ 
% % %             % monocular blur --> left eye only
% %             {'eLifePaper/explFilterSizes/18-10-18_500000iter_2_fsize6std_filtBoth_45_prob_1_seed2', scale, 2},
% %             {'eLifePaper/explFilterSizes/18-10-19_500000iter_3_fsize6std_filtBoth_45_prob_1_seed3', scale, 2},
% %             {'eLifePaper/explFilterSizes/18-10-18_500000iter_4_fsize6std_filtBoth_45_prob_1_seed4', scale, 2},
% %             {'eLifePaper/explFilterSizes/18-10-18_500000iter_2_fsize6std_filtBoth_45_prob_1_seed2', 2, 2},
% %             {'eLifePaper/explFilterSizes/18-10-19_500000iter_3_fsize6std_filtBoth_45_prob_1_seed3', 2, 2},
% %             {'eLifePaper/explFilterSizes/18-10-18_500000iter_4_fsize6std_filtBoth_45_prob_1_seed4', 2, 2},
% %             % monocular blur --> both eyes
%             {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_45_prob_1_seed1', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-18_500000iter_2_fsize6std_filtBoth_45_prob_1_seed2', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-19_500000iter_3_fsize6std_filtBoth_45_prob_1_seed3', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-18_500000iter_4_fsize6std_filtBoth_45_prob_1_seed4', scale, eye},
%             {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_45_prob_1_seed5', scale, eye},
% %             % },{ % complete monocular deprivation ?
% %             %     {'eLifePaper/explFilterSizes/', scale, 1},
% %             %     {'eLifePaper/explFilterSizes/', scale, 1},
% %             %     {'eLifePaper/explFilterSizes/', scale, 1},
%          },{ 
% %             % 10 degree strabism angle + const. object distance + const
% %             % vergence angle
% % %             {'eLifePaper/strabism/19-05-29_500000iter_1_fixAllAt6m_filtB_29_strabAngle_10_seed1', scale, eye},
% % %             {'eLifePaper/strabism/19-02-18_500000iter_2_fixAllAt6m_filtB_29_strabAngle_10_seed2', scale, eye},
% % %             {'eLifePaper/strabism/19-02-18_500000iter_3_fixAllAt6m_filtB_29_strabAngle_10_seed3', scale, eye},
% % %             {'eLifePaper/strabism/19-02-18_500000iter_4_fixAllAt6m_filtB_29_strabAngle_10_seed4', scale, eye},
% % %             {'eLifePaper/strabism/19-05-29_500000iter_5_fixAllAt6m_filtB_29_strabAngle_10_seed5', scale, eye},
% % %             {'eLifePaper/strabism/19-02-18_500000iter_2_fixAllAt6m_filtB_29_strabAngle_10_seed2', 2, eye},
% % %             {'eLifePaper/strabism/19-02-18_500000iter_3_fixAllAt6m_filtB_29_strabAngle_10_seed3', 2, eye},
% % %             {'eLifePaper/strabism/19-02-18_500000iter_4_fixAllAt6m_filtB_29_strabAngle_10_seed4', 2, eye},
% %             % 10 deg strabism with regular model and all object distances
%             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2', scale, eye},
%             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3', scale, eye},
%             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4', scale, eye},
%             {'eLifePaper/strabism/19-06-26_500000iter_1_learingActivefiltB_29_strabAngle_10_seed1', scale, eye},
%             {'eLifePaper/strabism/19-06-26_500000iter_5_learingActivefiltB_29_strabAngle_10_seed5', scale, eye},
% %             },
            {
%             % strabismic case (5 degree)
% %             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_5_seed2', scale, eye},
% %             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_5_seed3', scale, eye},
% %             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_5_seed4', scale, eye},
% %             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_5_seed2', 2, eye},
% %             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_5_seed3', 2, eye},
% %             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_5_seed4', 2, eye},
%             {'eLifePaper/strabism/19-02-11_500000iter_2_fixVergAngleAt6m_filtB_29_strabAngle_5_seed2', scale, eye},
%             {'eLifePaper/strabism/19-02-11_500000iter_3_fixVergAngleAt6m_filtB_29_strabAngle_5_seed3', scale, eye},
%             {'eLifePaper/strabism/19-02-11_500000iter_4_fixVergAngleAt6m_filtB_29_strabAngle_5_seed4', scale, eye},
%             % strabismic case (3 degree)
            % closer range and laplacian disparity distribution in the input
            {'eLifePaper/inducedStrabism/20-09-22_500000iter_2_inducedStrab_3deg_lapSig02_od05-1m', scale, eye},
%             {'eLifePaper/strabism/', scale, eye},
%             {'eLifePaper/strabism/', scale, eye},

%             {'eLifePaper/strabism/19-05-24_500000iter_1_fixAllAt6m_filtB_29_strabAngle_3_seed1', scale, eye}, % !!!
%             {'eLifePaper/strabism/19-02-18_500000iter_2_AllfixAt6m_filtB_29_strabAngle_3_seed2', scale, eye},
%             {'eLifePaper/strabism/19-02-18_500000iter_3_AllfixAt6m_filtB_29_strabAngle_3_seed3', scale, eye},
%             {'eLifePaper/strabism/19-02-18_500000iter_4_AllfixAt6m_filtB_29_strabAngle_3_seed4', scale, eye},
%             {'eLifePaper/strabism/19-05-24_500000iter_5_fixAllAt6m_filtB_29_strabAngle_3_seed5', scale, eye},
% %             {'eLifePaper/strabism/19-02-18_500000iter_2_AllfixAt6m_filtB_29_strabAngle_3_seed2', 2, eye},
% %             {'eLifePaper/strabism/19-02-18_500000iter_3_AllfixAt6m_filtB_29_strabAngle_3_seed3', 2, eye},
% %             {'eLifePaper/strabism/19-02-18_500000iter_4_AllfixAt6m_filtB_29_strabAngle_3_seed4', 2, eye},
%             % 10 degree strabism angle
% %             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2', scale, eye},
% %             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3', scale, eye},
% %             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4', scale, eye},
% %             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2', 2, eye},
% %             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3', 2, eye},
% %             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4', 2, eye},
        }

%        {   % normal case
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2', scale, eye},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3', scale, eye},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4', scale, eye},
% %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2', 2, eye},
% %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3', 2, eye},
% %             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4', 2, eye},
%         },{ % no RL
%             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_2_noLearninginitMethod_2_lapSigma_0_seed2', scale, eye},
%             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_3_noLearninginitMethod_2_lapSigma_0_seed3', scale, eye},
%             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_4_noLearninginitMethod_2_lapSigma_0_seed4', scale, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_2_noLearninginitMethod_2_lapSigma_0_seed2', 2, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_3_noLearninginitMethod_2_lapSigma_0_seed3', 2, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_4_noLearninginitMethod_2_lapSigma_0_seed4', 2, eye},
%         },{ % no RL & no Disparities
%             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_2_noLearning_initMethod_4_lapSigma_0_seed2', scale, eye},
%             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_3_noLearning_initMethod_4_lapSigma_0_seed3', scale, eye},
%             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_4_noLearning_initMethod_4_lapSigma_0_seed4', scale, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_2_noLearning_initMethod_4_lapSigma_0_seed2', 2, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_3_noLearning_initMethod_4_lapSigma_0_seed3', 2, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_4_noLearning_initMethod_4_lapSigma_0_seed4', 2, eye},
%         }
    }
%     names = {'normal', 'monocular', 'vertical', 'horizontal', 'strabismic', 'orthogonal'};
%     names = {{'Normal',''}, {'Vertical',''}, {'Horizontal',''}, {'Orthogonal',''}, {'Monocular',''}, {'Strabismic', '(10 deg)'}, {'Strabismic', '(3 deg)'}};
    names = {'normal', 'vertical', 'horizontal', 'orthogonal', 'monocular', 'strabismic', 'strabismic'};
%     names = {{'normal'}, 'vertical', 'horizontal', 'orthogonal',
%     'monocular', {'strabismic (10 deg)'}, {'strabismic (3 deg)'}}; % one
%     line
%     names = {'normal', 'vertical', 'horizontal', 'orthogonal', 'monocular', 'strabismic (10 deg)'};
%     names = {'normal', 'no RL', 'no RL & no disp.'};
%     names = {'','','','','','',''}
%     names = {['strabismic (3' char(176) ')']};
    names = {['strabismic (3 deg)']};
%    names = {['strabismic (5 deg)']};

    
    colors = [[77,175,74]; [55,126,184]; [152,78,163]; [255,127,0]; [228,26,28]; [166,86,40]; [118,42,131]]./256;
    colors = [[166,86,40]]./256; % just for strab 3 case
    
    %% calculate mean histograms of orientation representation
    threshold = 0.2;
    nBins = 11;
    useFreq = 1
    if scale == 1
        maxDisp = (50/240) * 8 * 4;
    elseif scale == 2
        maxDisp = (50/240) * 8;
    end
    
%     bins_h = -7.5:15:172.5;
    nBins = 12;
    bins_h = linspace(-maxDisp, maxDisp, nBins);
    
    histMeans = zeros(size(models,2), size(bins_h,2)-1);
    histStds = zeros(size(models,2), size(bins_h,2)-1);
    nFits = zeros(size(models,2), 1);
    
%     disparityVals = {};

    for i = 1:length(models)
        tmp_N = 1:size(models{i}); % relative nums for each model
        allVals = zeros(size(bins_h,2)-1, size(tmp_N,2));
        fits = zeros(size(tmp_N,2), 1);
        
%         dispVals = [];
        
        for j = 1:length(tmp_N)
            if ~useFreq
                para = load(strcat(filePath, models{i}{j}{1}, '/modelAt500000/scModel', num2str(models{i}{j}{2}), 'eyes', num2str(models{i}{j}{3}), '.mat'));
            else
                para = load(strcat(filePath, models{i}{j}{1}, '/modelAt500000/scModel', num2str(models{i}{j}{2}), 'eyes', num2str(models{i}{j}{3}), '_freq.mat'));
            end
            Resnorm_Set = para.Error;
            pB = para.Fitted_Gabor_Filter_Parameters;
            
            idx = find(Resnorm_Set < threshold);
            
            if ~useFreq
                dispPref = (pB(idx, 3).*(pB(idx, 6)-pB(idx, 5)))./(2*pi*cos(pB(idx, 2)));
            else
%                 dispPref = ((pB(idx, 6)-pB(idx, 5)))./(2*pi*cos(pB(idx, 2)))./pB(idx, 3);
                dispPref = ((pB(idx, 6)-pB(idx, 5)))./(2*pi*cos(pB(idx, 2)).*pB(idx, 3));
            end
            
%             filteredVals = [];
%             for d = 1:length(dispPref)
%                 if ~(dispPref(d) > maxDisp) & ~(dispPref(d) < -maxDisp)
%                     filteredVals = [filteredVals; dispPref(d)];
%                 end
%             end
                        
%             dispVals = [dispVals; filteredVals];
            N = histcounts(dispPref,bins_h);
            
%             if (j == 2) % s4 num of used fits
%                 bf_nums(i) = length(idx);
%             end

%             bins_h = -7.5:15:172.5;

%             N = histcounts(mod(pB(idx,2)*180/pi+7.5, 180)-7.5,bins_h);
            fits(j) = sum(N);
%             tmp_N(j) = N(orientation)/sum(N);
            allVals(:, j) = N./fits(j);
        end

%         if (i == 1)
%             ttest_sample1 = tmp_N;
%         end
% 
%         if (i == size(models, 1))
%             ttest_sample2 = tmp_N;
%         end
% 
%         Ns(i) = mean(tmp_N);
%         Ns_err(i) = std(tmp_N);
        histMeans(i, :) = mean(allVals,2);
        histStds(i, :) = std(allVals,0,2);
        nFits(i) = mean(fits);
%         disparityVals{i} = dispVals;
    end
    
%     figure('Position', [100 100 1200 800]);
%     figure('Position', [100 100 1200 400]);
%     figure('Position', [100 100 1500 250]); % cuts of some of the labels??
%     figure('Position', [100 200 1500 350]);
    
%     figure('Position',[100 100 1166 534]); % for two columns
    figure('position', [300 350 1000 625]); % for consistent layout
%     hold on;
    for i = 1:length(models)
%     for i = 1:length(models)-1
%         if (i < 5)
%             subplot(2, 4, i);
%         else
%             subplot(2, 4, i+1);
%         end
%         subplot(2,length(models), i)
%         subplot(1,length(models), i)


        subplot(2,3,i)
%         subplot(1,1,i)
        
        hold on;
        title(names{i});
        bar(histMeans(i,:), 'FaceColor', colors(i,:));
        errorbar(histMeans(i,:), histStds(i,:),'k', 'CapSize', 4, 'linestyle', 'none')
%         bar(histMeans(i,:), 'FaceColor', colors(i,:), 'EdgeColor', colors(i,:));
%         errorbar(histMeans(i,:), histStds(i,:), 'color', colors(i,:), 'CapSize', 4, 'linestyle', 'none')
%         xticks([1, 2, 3, 4])%, 5, 6, 7])
%         xticklabels({[sprintf('Hunt (v)')], sprintf('McG (v)'), sprintf('Hunt (h)'), sprintf('McG (h)')});
%         xlim([0.5, 4.5])
%         ylim([0 50])
%         if (i == 5) %|| (i == 4)
            xlabel('pref. disparity [deg]')
%         end
%         if rem(i, 4) == 1
        if (i == 1) || (i == 4)
%             ylabel('number of bases [%]')
            ylabel('number of RFs [%]')
        end
        % xticks([1, nBins/2, nBins-1])
        % xticklabels({'0', '90'});
%         xticks(linspace(1, nBins-1, 5))
%         prec = 2;
        % xticklabels({num2str(-maxDisp, prec), num2str(-maxDisp/2, prec), 0, num2str(maxDisp/2, prec),num2str(-maxDisp, prec)});
        if scale == 1
            xticks(linspace(1, nBins-1, 3))
            xticklabels({'-6.6', '0', '6.6'})
%             xticks(linspace(1, nBins-1, 5))
%             xticklabels({'-6.6', '-3.3', '0', '3.3', '6.6'})
        elseif scale == 2
%             xticks(linspace(1, nBins-1, 3))
%             xticklabels({'-1.6', '0', '1.6'})
            xticks(linspace(1, nBins-1, 5))
            xticklabels({'-1.6', '-0.8', '0', '0.8', '1.6'})
        end
        
        if scale == 1
            % for wave fits
%             ylim([0, 0.9]);
%             yticks([0, 0.3, 0.6, 0.9]);
%             yticklabels({0, 30, 60, 90});
%             t = text(7, 0.85, strcat(['$\bar{N}$=',num2str(nFits(i), '%.0f')]));
% %             t = text(7, 0.5, strcat(['$\bar{N}$=',num2str(nFits(i), '%.0f')]));
%             set(t,'interpreter','Latex','FontSize',12);
            % for frequencies
            ylim([0, 0.8]);
            yticks([0, 0.2, 0.4, 0.6, 0.8]);
            yticklabels({0, 20, 40, 60, 80});
%             t = text(7, 0.725, strcat(['$\bar{N}\mathsf{=',num2str(nFits(i), '%.0f}$')]));
            t = text(8, 0.725, strcat(['$\bar{N}\textsf{=',num2str(nFits(i), '%.0f}$')]));
%             t = text(7, 0.5, strcat(['$\bar{N}$=',num2str(nFits(i), '%.0f')]));
            set(t,'interpreter','Latex','FontSize',11);
        elseif scale == 2
            % for wave fits
%             ylim([0, 0.6]);
%             yticks([0, 0.2, 0.4, 0.6]);
%             yticklabels({0, 20, 40, 60});
%             t = text(7, 0.55, strcat(['$\bar{N}$=',num2str(nFits(i), '%.0f')]));
%             set(t,'interpreter','Latex','FontSize',12);
            % for frequency fits
            ylim([0, 0.6])
            yticks([0, 0.2, 0.4, 0.6]);
            yticklabels({0, 20, 40, 60});
            t = text(8, 0.55, strcat(['$\bar{N}\mathsf{=',num2str(nFits(i), '%.0f}$')]));
            set(t,'interpreter','Latex','FontSize',11);
        end
%         ylim([0, inf]);
%         grid on;
        % text(9, 0.55, sprintf('|N|=%.0f',nFits(i)))
        % t = text(9, 0.55, sprintf('$\bar{N}$=%.0f',nFits(i)))
    end
    
    sprintf('Done')
    if ~isempty(saveTag)
        saveas(gcf, sprintf('%sdisparities_scale%d_%s.png', savePath, scale, saveTag));
    end
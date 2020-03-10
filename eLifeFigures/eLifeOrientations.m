%% This script gathers the data from all accoding simulations and generates
%% an overview of the edge representation

%% depicted are averaged results from normal, vertical, horizontal, orthogonal,
%% strabismic and monocular rearing
function eLifeOrientations(scale, saveTag)
    filePath = '/home/aecgroup/aecdata/Results/'; %'/SAB2018/' or 'eLifePaper'
    savePath = '/home/aecgroup/aecdata/Results/eLifePaper/plots/Orientations/'

    scale = 1; % now try to combine both scales in one plot
    eye = 2; % 1: binocular fits, 2: monocular left, 3: monocular right
    displayAvgNumbers = 1;
    
    models = {
        {   % normal case
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1', scale, eye},
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2', scale, eye},
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3', scale, eye},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4', scale, eye},
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5', scale, eye},
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1', 2, eye},
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2', 2, eye},
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3', 2, eye},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4', 2, eye},
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5', 2, eye},
            

%         },{ % normal case, natural image database from Hunt et al.
%             
%             {'eLifePaper/HHDatabase/19-05-30_500000iter_1_HHDB_filtBoth_29_prob_1_seed1', scale, eye},
%             {'eLifePaper/HHDatabase/19-05-30_500000iter_2_HHDB_filtBoth_29_prob_1_seed2', scale, eye},
%             {'eLifePaper/HHDatabase/19-05-30_500000iter_3_HHDB_filtBoth_29_prob_1_seed3', scale, eye},
%             {'eLifePaper/HHDatabase/19-05-30_500000iter_4_HHDB_filtBoth_29_prob_1_seed4', scale, eye},
%             {'eLifePaper/HHDatabase/19-05-30_500000iter_5_HHDB_filtBoth_29_prob_1_seed5', scale, eye},
%             {'eLifePaper/HHDatabase/19-05-30_500000iter_1_HHDB_filtBoth_29_prob_1_seed1', 2, eye},
%             {'eLifePaper/HHDatabase/19-05-30_500000iter_2_HHDB_filtBoth_29_prob_1_seed2', 2, eye},
%             {'eLifePaper/HHDatabase/19-05-30_500000iter_3_HHDB_filtBoth_29_prob_1_seed3', 2, eye},
%             {'eLifePaper/HHDatabase/19-05-30_500000iter_4_HHDB_filtBoth_29_prob_1_seed4', 2, eye},
%             {'eLifePaper/HHDatabase/19-05-30_500000iter_5_HHDB_filtBoth_29_prob_1_seed5', 2, eye},
% 
%          }
        },{ % vertical only (filter size 33 px, filter name: 30)
            {'eLifePaper/explFilterSizes/19-05-23_500000iter_1_fsize6std_filtBoth_30_prob_1_seed1', scale, eye},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_30_prob_1_seed2', scale, eye},
            {'eLifePaper/explFilterSizes/18-10-01_500000iter_3_fsize6std_filtBoth_30_prob_1_seed3', scale, eye},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_30_prob_1_seed4', scale, eye},
            {'eLifePaper/explFilterSizes/19-05-24_500000iter_5_fsize6std_filtBoth_30_prob_1_seed5', scale, eye},
            {'eLifePaper/explFilterSizes/19-05-23_500000iter_1_fsize6std_filtBoth_30_prob_1_seed1', 2, eye},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_30_prob_1_seed2', 2, eye},
            {'eLifePaper/explFilterSizes/18-10-01_500000iter_3_fsize6std_filtBoth_30_prob_1_seed3', 2, eye},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_30_prob_1_seed4', 2, eye},
            {'eLifePaper/explFilterSizes/19-05-24_500000iter_5_fsize6std_filtBoth_30_prob_1_seed5', 2, eye},
% %             % vertical only (filter size 100 px, filter name: 31)
% %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_31_prob_1_seed2', scale, eye},
% %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_31_prob_1_seed3', scale, eye},
% %             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_31_prob_1_seed4', scale, eye},
% %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_31_prob_1_seed2', 2, eye},
% %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_31_prob_1_seed3', 2, eye},
% %             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_31_prob_1_seed4', 2, eye},
        },{ % horizontal only (33 px, filter name: 34)
            {'eLifePaper/explFilterSizes/19-05-23_500000iter_1_fsize6std_filtBoth_34_prob_1_seed1', scale, eye},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_34_prob_1_seed2', scale, eye},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_3_fsize6std_filtBoth_34_prob_1_seed3', scale, eye},
            {'eLifePaper/explFilterSizes/18-10-01_500000iter_4_fsize6std_filtBoth_34_prob_1_seed4', scale, eye},
            {'eLifePaper/explFilterSizes/19-05-24_500000iter_5_fsize6std_filtBoth_34_prob_1_seed5', scale, eye},
            {'eLifePaper/explFilterSizes/19-05-23_500000iter_1_fsize6std_filtBoth_34_prob_1_seed1', 2, eye},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_34_prob_1_seed2', 2, eye},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_3_fsize6std_filtBoth_34_prob_1_seed3', 2, eye},
            {'eLifePaper/explFilterSizes/18-10-01_500000iter_4_fsize6std_filtBoth_34_prob_1_seed4', 2, eye},
            {'eLifePaper/explFilterSizes/19-05-24_500000iter_5_fsize6std_filtBoth_34_prob_1_seed5', 2, eye},
% %             % horizontal only (100 px, filter name: 35)
% %             {'eLifePaper/explFilterSizes/18-09-29_500000iter_2_fsize6std_filtBoth_35_prob_1_seed2', scale, eye},
% %             {'eLifePaper/explFilterSizes/18-09-29_500000iter_3_fsize6std_filtBoth_35_prob_1_seed3', scale, eye},
% %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_4_fsize6std_filtBoth_35_prob_1_seed4', scale, eye},
% %             {'eLifePaper/explFilterSizes/18-09-29_500000iter_2_fsize6std_filtBoth_35_prob_1_seed2', 2, eye},
% %             {'eLifePaper/explFilterSizes/18-09-29_500000iter_3_fsize6std_filtBoth_35_prob_1_seed3', 2, eye},
% %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_4_fsize6std_filtBoth_35_prob_1_seed4', 2, eye},
% 
        },{ 
            % orthogonal left eye (33 px, filter name: 46)
            {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_46_prob_1_seed1', scale, 2},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_2_fsize6std_filtBoth_46_prob_1_seed2', scale, 2},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_3_fsize6std_filtBoth_46_prob_1_seed3', scale, 2},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_4_fsize6std_filtBoth_46_prob_1_seed4', scale, 2},
            {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_46_prob_1_seed5', scale, 2},
            {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_46_prob_1_seed1', 2, 2},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_2_fsize6std_filtBoth_46_prob_1_seed2', 2, 2},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_3_fsize6std_filtBoth_46_prob_1_seed3', 2, 2},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_4_fsize6std_filtBoth_46_prob_1_seed4', 2, 2},
            {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_46_prob_1_seed5', 2, 2},
%             
%         },{ % orthogonal right eye
%             {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_46_prob_1_seed1', scale, 3},
%             {'eLifePaper/explFilterSizes/18-10-22_500000iter_2_fsize6std_filtBoth_46_prob_1_seed2', scale, 3},
%             {'eLifePaper/explFilterSizes/18-10-22_500000iter_3_fsize6std_filtBoth_46_prob_1_seed3', scale, 3},
%             {'eLifePaper/explFilterSizes/18-10-22_500000iter_4_fsize6std_filtBoth_46_prob_1_seed4', scale, 3},
%             {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_46_prob_1_seed5', scale, 3},
%             
%             {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_46_prob_1_seed1', 2, 3},
%             {'eLifePaper/explFilterSizes/18-10-22_500000iter_2_fsize6std_filtBoth_46_prob_1_seed2', 2, 3},
%             {'eLifePaper/explFilterSizes/18-10-22_500000iter_3_fsize6std_filtBoth_46_prob_1_seed3', 2, 3},
%             {'eLifePaper/explFilterSizes/18-10-22_500000iter_4_fsize6std_filtBoth_46_prob_1_seed4', 2, 3},
%             {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_46_prob_1_seed5', 2, 3},
           },{ 
% %             % monocular blur --> left eye only
% %             {'eLifePaper/explFilterSizes/18-10-18_500000iter_2_fsize6std_filtBoth_45_prob_1_seed2', scale, 2},
% %             {'eLifePaper/explFilterSizes/18-10-19_500000iter_3_fsize6std_filtBoth_45_prob_1_seed3', scale, 2},
% %             {'eLifePaper/explFilterSizes/18-10-18_500000iter_4_fsize6std_filtBoth_45_prob_1_seed4', scale, 2},
% %             {'eLifePaper/explFilterSizes/18-10-18_500000iter_2_fsize6std_filtBoth_45_prob_1_seed2', 2, 2},
% %             {'eLifePaper/explFilterSizes/18-10-19_500000iter_3_fsize6std_filtBoth_45_prob_1_seed3', 2, 2},
% %             {'eLifePaper/explFilterSizes/18-10-18_500000iter_4_fsize6std_filtBoth_45_prob_1_seed4', 2, 2},
%             % monocular blur (33 px, filter name: 45)
            {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_45_prob_1_seed1', scale, eye},
            {'eLifePaper/explFilterSizes/18-10-18_500000iter_2_fsize6std_filtBoth_45_prob_1_seed2', scale, eye},
            {'eLifePaper/explFilterSizes/18-10-19_500000iter_3_fsize6std_filtBoth_45_prob_1_seed3', scale, eye},
            {'eLifePaper/explFilterSizes/18-10-18_500000iter_4_fsize6std_filtBoth_45_prob_1_seed4', scale, eye},
            {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_45_prob_1_seed5', scale, eye},
            {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_45_prob_1_seed1', 2, eye},
            {'eLifePaper/explFilterSizes/18-10-18_500000iter_2_fsize6std_filtBoth_45_prob_1_seed2', 2, eye},
            {'eLifePaper/explFilterSizes/18-10-19_500000iter_3_fsize6std_filtBoth_45_prob_1_seed3', 2, eye},
            {'eLifePaper/explFilterSizes/18-10-18_500000iter_4_fsize6std_filtBoth_45_prob_1_seed4', 2, eye},
            {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_45_prob_1_seed5', 2, eye},
%             
%             % },{ % complete monocular deprivation ?
%             %     {'eLifePaper/explFilterSizes/', scale, 1},
%             %     {'eLifePaper/explFilterSizes/', scale, 1},
%             %     {'eLifePaper/explFilterSizes/', scale, 1},
% 
        },{ 
%         
%             % strabismic case (5 degree)
% %             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_5_seed2', scale, 1},
% %             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_5_seed3', scale, 1},
% %             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_5_seed4', scale, 1},
% %             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_5_seed2', 2, 1},
% %             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_5_seed3', 2, 1},
% %             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_5_seed4', 2, 1},
%             % 10 degree strabism angle
            {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2', scale, eye},
            {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3', scale, eye},
            {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4', scale, eye},
            {'eLifePaper/strabism/19-06-26_500000iter_1_learingActivefiltB_29_strabAngle_10_seed1', scale, eye},
            {'eLifePaper/strabism/19-06-26_500000iter_5_learingActivefiltB_29_strabAngle_10_seed5', scale, eye},
            {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2', 2, eye},
            {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3', 2, eye},
            {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4', 2, eye},
            {'eLifePaper/strabism/19-06-26_500000iter_1_learingActivefiltB_29_strabAngle_10_seed1', 2, eye},
            {'eLifePaper/strabism/19-06-26_500000iter_5_learingActivefiltB_29_strabAngle_10_seed5', 2, eye},
%             % 10 degree strabism angle, vergence fixed at 6 m
%             {'eLifePaper/strabism/19-05-29_500000iter_1_fixAllAt6m_filtB_29_strabAngle_10_seed1', scale, eye},
%             {'eLifePaper/strabism/19-02-18_500000iter_2_fixAllAt6m_filtB_29_strabAngle_10_seed2', scale, eye},
%             {'eLifePaper/strabism/19-02-18_500000iter_3_fixAllAt6m_filtB_29_strabAngle_10_seed3', scale, eye},
%             {'eLifePaper/strabism/19-02-18_500000iter_4_fixAllAt6m_filtB_29_strabAngle_10_seed4', scale, eye},
%             {'eLifePaper/strabism/19-05-29_500000iter_5_fixAllAt6m_filtB_29_strabAngle_10_seed5', scale, eye},
%             {'eLifePaper/strabism/19-05-29_500000iter_1_fixAllAt6m_filtB_29_strabAngle_10_seed1', 2, eye},
%             {'eLifePaper/strabism/19-02-18_500000iter_2_fixAllAt6m_filtB_29_strabAngle_10_seed2', 2, eye},
%             {'eLifePaper/strabism/19-02-18_500000iter_3_fixAllAt6m_filtB_29_strabAngle_10_seed3', 2, eye},
%             {'eLifePaper/strabism/19-02-18_500000iter_4_fixAllAt6m_filtB_29_strabAngle_10_seed4', 2, eye},
%             {'eLifePaper/strabism/19-05-29_500000iter_5_fixAllAt6m_filtB_29_strabAngle_10_seed5', 2, eye},
        }
%     
% %            {   % normal case
% %             {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1', scale, eye},           
% %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2', scale, eye},
% %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3', scale, eye},
% %             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4', scale, eye},
% %             {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5', scale, eye},
% % %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2', 2, eye},
% % %             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3', 2, eye},
% % %             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4', 2, eye},
% %         },{ % no RL
% %             {'eLifePaper/vergenceInfluence/19-05-15_500000iter_1_noLearning_initMethod_2_lapSigma_0_seed1', scale, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_2_noLearninginitMethod_2_lapSigma_0_seed2', scale, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_3_noLearninginitMethod_2_lapSigma_0_seed3', scale, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_4_noLearninginitMethod_2_lapSigma_0_seed4', scale, eye},
% %             {'eLifePaper/vergenceInfluence/19-05-15_500000iter_5_noLearning_initMethod_2_lapSigma_0_seed5', scale, eye},
% % %             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_2_noLearninginitMethod_2_lapSigma_0_seed2', 2, eye},
% % %             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_3_noLearninginitMethod_2_lapSigma_0_seed3', 2, eye},
% % %             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_4_noLearninginitMethod_2_lapSigma_0_seed4', 2, eye},
% %         },{ % no RL & no Disparities
% %             {'eLifePaper/vergenceInfluence/19-05-16_500000iter_1_noLearning_initMethod_4_lapSigma_0_seed1', scale, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_2_noLearning_initMethod_4_lapSigma_0_seed2', scale, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_3_noLearning_initMethod_4_lapSigma_0_seed3', scale, eye},
% %             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_4_noLearning_initMethod_4_lapSigma_0_seed4', scale, eye},
% %             {'eLifePaper/vergenceInfluence/19-05-16_500000iter_5_noLearning_initMethod_4_lapSigma_0_seed5', scale, eye},
% % %             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_2_noLearning_initMethod_4_lapSigma_0_seed2', 2, eye},
% % %             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_3_noLearning_initMethod_4_lapSigma_0_seed3', 2, eye},
% % %             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_4_noLearning_initMethod_4_lapSigma_0_seed4', 2, eye},
% %         }

    }
%     names = {'normal', 'monocular', 'vertical', 'horizontal', 'strabismic', 'orthogonal (left)', 'orthogonal (right)'};
%     names = {'Normal', 'Monocular', 'Vertical', 'Horizontal', 'Strabismic', 'Orthogonal (left)', 'Orthogonal (right)'};
%     names = {'Normal', 'Vertical', 'Horizontal', 'Monocular', 'Orthogonal (left eye)', 'Orthogonal (right eye)','Strabismic'};
%     names = {{'Normal',''}, {'Vertical',''}, {'Horizontal',''}, {'Monocular',''}, {'Orthogonal','(Left Eye)'}, {'Orthogonal','(Right Eye)'},  {'Strabismic', ''}};
    names = {'normal', 'vertical', 'horizontal', 'orthogonal', 'monocular',  'strabismic'};
%     names = {'McGill image database (man-made)', 'Natural images as used in Hoyer & Hyvärinen (2000)'};
%     names = {{'Man-made structures from', 'McGill image database'}, {'Natural images', 'as used in Hoyer & Hyvärinen (2000)'}};
    colors = [[77,175,74]; [55,126,184]; [152,78,163]; [255,127,0]; [228,26,28]; [166,86,40]; [118,42,131]]./256;
%     colors = [[77,175,74]; [77,175,74]]./256;
    
    %% calculate mean histograms of orientation representation
    threshold = 0.2;
    bins_h = -7.5:15:172.5;
%     bins_h = -15:15:180-15;
    histMeans = zeros(size(models,2), size(bins_h,2)-1);
    histStds = zeros(size(models,2), size(bins_h,2)-1);
    nFits = zeros(size(models,2), 1);
%     orientationVals = {}; 
    
    for i = 1:length(models)
        tmp_N = 1:size(models{i}); % relative nums for each model
        allVals = zeros(size(bins_h,2)-1, size(tmp_N,2));
        fits = zeros(size(tmp_N,2), 1);
%         rawV = [];
        
        for j = 1:length(tmp_N)
%             para = load(strcat(filePath, models{i}{j}{1}, '/modelAt500000/scModel', num2str(models{i}{j}{2}), "eyes", num2str(models{i}{j}{3}), '.mat'));
            para = load(strcat(filePath, models{i}{j}{1}, '/modelAt500000/scModel', num2str(models{i}{j}{2}), "eyes", num2str(models{i}{j}{3}), '_freq.mat'));
            Resnorm_Set = para.Error;
            Parameter_Set = para.Fitted_Gabor_Filter_Parameters;
            idx = find(Resnorm_Set < threshold);

%             if (j == 2) % s4 num of used fits
%                 bf_nums(i) = length(idx);
%             end

%             bins_h = -7.5:15:172.5;
%             rawV = [rawV; mod(Parameter_Set(idx,2)*180/pi+7.5, 180)-7.5];
            N = histcounts(mod(Parameter_Set(idx,2)*180/pi+7.5, 180)-7.5,bins_h);
%             N = histcounts(mod(Parameter_Set(idx,2)*180/pi, 180),bins_h);
            
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
%         orientationVals{i} = rawV;
    end
    
%     figure('Position', [100 100 1200 800]);
%     figure('Position', [100 100 1200 300]);
%     figure('Position',[100 100 1166 534]); % for separate ortho case
    
%     h = figure('position', [300 400 1400 800]); % for consistend layout
    figure('position', [300 350 1000 625]); 
%     set(gca, 'TickLabelInterpreter', 'latex')
    
%     hold on;
    for j = 1:length(models)%-1
%         if (i <= 5)
%         if (i <= 4)
%             subplot(2, 4, i);
%         else
%             subplot(2, 4, i+1);
%         end
%         subplot(1, length(models), i);
        
        subplot(2,3,j)
%         subplot(1,2,j)
        
%         % hotfix to remove right eye ortho case
%         if j < 5
            i = j;
%         else
%             i = j+1;
%         end
        
        hold on;
        title(names{i});
        bar(histMeans(i,:), 'FaceColor', colors(j,:));
        errorbar(histMeans(i,:), histStds(i,:), 'k', 'linestyle', 'none', 'CapSize', 4)
%         bar(histMeans(i,:), 'FaceColor', colors(j,:), 'EdgeColor', colors(j,:));
%         errorbar(histMeans(i,:), histStds(i,:), 'color', colors(j,:), 'linestyle', 'none', 'CapSize', 4)
%         xticks([1, 2, 3, 4])%, 5, 6, 7])
%         xticklabels({[sprintf('Hunt (v)')], sprintf('McG (v)'), sprintf('Hunt (h)'), sprintf('McG (h)')});
%         xlim([0.5, 4.5])
%         ylim([0 50])
%         if i == 1
        if j == 5
            xlabel('orientation [deg]')
        end
%         if rem(i, 4) == 1
%         if (i == 1) || (i == 2) || (i == 5) % only for presentations
        if (j == 1) || (j == 4)
%             ylabel('Number of Bases [%]')
            ylabel('number of bases [%]')
        end
        xticks([1, 7])
        xticklabels({'0', '90'});
        ylim([0, 0.51]);
%         yticks([0, 0.2, 0.4, 0.6]);
%         yticklabels({0, 20, 40, 60});
%         ylim([0, 0.4]);
%         yticks([0, 0.25, 0.5]);
%         yticklabels({0, 25, 50});
        yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5]);
        yticklabels({0, 10, 20, 30, 40, 50});
%         grid minor
        
%         grid on
%         text(9, 0.55, sprintf('|N|=%.0f',nFits(i)))
%         t = text(8, 0.55, strcat(['$\bar{N}$=',num2str(nFits(i), '%.0f')]));

        if displayAvgNumbers
%             t = text(9, 0.55, strcat(['$\bar{N}$=',num2str(nFits(i), '%.0f')]));
%             t = text(9, 0.55, strcat(['$\bar{N}$=',num2str(nFits(i)*2, '%.0f')]));
%             t = text(9, 0.45, strcat(['$\bar{N}$=',num2str(nFits(i)*2, '%.0f')])); % only if both scales are combined
%             t = text(9, 0.45, strcat(['\textsf{$\bar{N}$=',num2str(nFits(i)*2, '%.0f}')])); 
%             t = text(9, 0.45, strcat(['$\bar{N}\mathsf{=',num2str(nFits(i)*2, '%.0f}$')])); 
            t = text(9, 0.45, strcat(['$\bar{N}\mathsf{=',num2str(nFits(i)*2, '%.0f}$')])); 
            set(t,'interpreter','latex','FontSize',11);
        end
    end
    
    sprintf('Done')
    if ~isempty(saveTag)
        saveas(gcf, sprintf('%scompareOrientations_scale%d_%s.png', savePath, scale, saveTag));
    end
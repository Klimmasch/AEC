%% This script compares multiple simulations and investiges their binocularity values
function eLifeBinocularities(saveName, bins)
    filePath = '/home/aecgroup/aecdata/Results/'; %'/SAB2018/' or 'eLifePaper'
    savePath = '/home/aecgroup/aecdata/Results/eLifePaper/plots/Binocularity/'
%     method = 'energy'
    method = 'responses'
    
    thresh = 0.2;
    useFreq = 1
%     filePath = '/home/mrklim/fiasAECGroup/aecdata/Results/'; %'/SAB2018/' or 'eLifePaper'
%     savePath = '/home/mrklim/fiasAECGroup/aecdata/Results/eLifePaper/plots/'

    if length(bins) == 1
        nbins = bins;
    else
        nbins = length(bins)-1;
    end

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
        },{ % horizontal only (33 px)
            {'eLifePaper/explFilterSizes/19-05-23_500000iter_1_fsize6std_filtBoth_34_prob_1_seed1'},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_2_fsize6std_filtBoth_34_prob_1_seed2'},
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_3_fsize6std_filtBoth_34_prob_1_seed3'},
            {'eLifePaper/explFilterSizes/18-10-01_500000iter_4_fsize6std_filtBoth_34_prob_1_seed4'},
            {'eLifePaper/explFilterSizes/19-05-24_500000iter_5_fsize6std_filtBoth_34_prob_1_seed5'},
        },{ % orthogonal left eye
            {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_46_prob_1_seed1'},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_2_fsize6std_filtBoth_46_prob_1_seed2'},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_3_fsize6std_filtBoth_46_prob_1_seed3'},
            {'eLifePaper/explFilterSizes/18-10-22_500000iter_4_fsize6std_filtBoth_46_prob_1_seed4'},
            {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_46_prob_1_seed5'},
        },{% monocular blur --> left eye only or both?
            {'eLifePaper/explFilterSizes/19-06-04_500000iter_1_fsize6std_filtBoth_45_prob_1_seed1'},
            {'eLifePaper/explFilterSizes/18-10-18_500000iter_2_fsize6std_filtBoth_45_prob_1_seed2'},
            {'eLifePaper/explFilterSizes/18-10-19_500000iter_3_fsize6std_filtBoth_45_prob_1_seed3'},
            {'eLifePaper/explFilterSizes/18-10-18_500000iter_4_fsize6std_filtBoth_45_prob_1_seed4'},
            {'eLifePaper/explFilterSizes/19-06-05_500000iter_5_fsize6std_filtBoth_45_prob_1_seed5'},
        },{ % strabismic case (5 degree)
%             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_5_seed2'},
%             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_5_seed3'},
%             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_5_seed4'},
%             % 10 degree strabism angle (fixed)
%             {'eLifePaper/strabism/19-05-29_500000iter_1_fixAllAt6m_filtB_29_strabAngle_10_seed1'},
%             {'eLifePaper/strabism/19-02-18_500000iter_2_fixAllAt6m_filtB_29_strabAngle_10_seed2'},
%             {'eLifePaper/strabism/19-02-18_500000iter_3_fixAllAt6m_filtB_29_strabAngle_10_seed3'},
%             {'eLifePaper/strabism/19-02-18_500000iter_4_fixAllAt6m_filtB_29_strabAngle_10_seed4'},
%             {'eLifePaper/strabism/19-05-29_500000iter_5_fixAllAt6m_filtB_29_strabAngle_10_seed5'},
            % 10 deg strabism (eye free to move)
            {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2'},
            {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3'},
            {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4'},
            {'eLifePaper/strabism/19-06-26_500000iter_1_learingActivefiltB_29_strabAngle_10_seed1'},
            {'eLifePaper/strabism/19-06-26_500000iter_5_learingActivefiltB_29_strabAngle_10_seed5'},

%             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2'},
%             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3'},
%             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4'},
        }
    };
    labelz = {'normal', 'vertical', 'horizontal', 'orthogonal', 'monocular', 'strabismic'};
%     labelz = {'Normal', 'Vertical', 'Horizontal', 'Orthogonal', 'Monocular', 'Strabismic'};



%     models = {
%         {   % normal case
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2'},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3'},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4'},
%         },{
%             % objSize = 2 m
%             {'eLifePaper/spreadBinocularity/18-11-05_500000iter_2_test_objSize=2m'},
%             {'eLifePaper/spreadBinocularity/18-11-06_500000iter_3_test_objSize=2m_s3'},
%             {'eLifePaper/spreadBinocularity/18-11-06_500000iter_4_test_objSize=2m_s4'},
%
%         },{ % objSize = 1.5 m
%             {'eLifePaper/spreadBinocularity/18-11-07_500000iter_2_objSize_1.5_patchSize_8_seed2'},
%             {'eLifePaper/spreadBinocularity/18-11-07_500000iter_3_objSize_1.5_patchSize_8_seed3'}, % not completed
%             {'eLifePaper/spreadBinocularity/18-11-07_500000iter_4_objSize_1.5_patchSize_8_seed4'}, % not completed
%         },{ % objSize = 1 m
%             {'eLifePaper/spreadBinocularity/18-11-08_500000iter_2_objSize_1.0_patchSize_8_seed2'},
%             {'eLifePaper/spreadBinocularity/18-11-08_500000iter_3_objSize_1.0_patchSize_8_seed3'},
% %             {'eLifePaper/spreadBinocularity/18-11-06_500000iter_2_test_objSize=1m'},   % please replace with another seed
%         }
%     }
%     labelz = {'objSize = 3m', 'objSize = 2m', 'objSize = 1.5m', 'objSize = 1m'};


%     models = {
%             {   % s = 0
%                 {'SAB2018/laplacianPolicy/18-07-26_500000iter_2_noLearning_initMethod_4_lapSigma_0.00_seed2'},
%                 {'SAB2018/laplacianPolicy/18-07-23_500000iter_3_noLearning_initMethod_4_lapSigma_0.00_seed3'},
%                 {'SAB2018/laplacianPolicy/18-07-21_500000iter_4_noLearning_initMethod_4_lapSigma_0.00_seed4'},
%             },{%
% %                 {'SAB2018/laplacianPolicy/18-07-26_500000iter_2_noLearning_initMethod_4_lapSigma_0.22_seed2'},
% %                 {'SAB2018/laplacianPolicy/18-07-23_500000iter_3_noLearning_initMethod_4_lapSigma_0.22_seed3'},
% %                 {'SAB2018/laplacianPolicy/18-07-21_500000iter_4_noLearning_initMethod_4_lapSigma_0.22_seed4'},
% %             },{
% %                 {'SAB2018/laplacianPolicy/18-07-26_500000iter_2_noLearning_initMethod_4_lapSigma_0.50_seed2'},
% %                 {'SAB2018/laplacianPolicy/18-07-23_500000iter_3_noLearning_initMethod_4_lapSigma_0.50_seed3'},
% %                 {'SAB2018/laplacianPolicy/18-07-21_500000iter_4_noLearning_initMethod_4_lapSigma_0.50_seed4'},
% %             },{
%                 {'SAB2018/laplacianPolicy/18-07-26_500000iter_2_noLearning_initMethod_4_lapSigma_1.00_seed2'},
%                 {'SAB2018/laplacianPolicy/18-07-23_500000iter_3_noLearning_initMethod_4_lapSigma_1.00_seed3'},
%                 {'SAB2018/laplacianPolicy/18-07-21_500000iter_4_noLearning_initMethod_4_lapSigma_1.00_seed4'},
%             },{
%                 {'SAB2018/laplacianPolicy/18-07-26_500000iter_2_noLearning_initMethod_4_lapSigma_2.00_seed2'},
%                 {'SAB2018/laplacianPolicy/18-07-23_500000iter_3_noLearning_initMethod_4_lapSigma_2.00_seed3'},
%                 {'SAB2018/laplacianPolicy/18-07-21_500000iter_4_noLearning_initMethod_4_lapSigma_2.00_seed4'},
%             },{
%                 {'SAB2018/laplacianPolicy/18-08-11_500000iter_2_noLearning_initMethod_4_lapSigma_5.00_seed2'},
%                 {'SAB2018/laplacianPolicy/18-08-13_500000iter_3_noLearning_initMethod_4_lapSigma_5.00_seed3'},
%                 {'SAB2018/laplacianPolicy/18-08-11_500000iter_4_noLearning_initMethod_4_lapSigma_5.00_seed4'},
%             },{
%                 {'SAB2018/laplacianPolicy/18-08-11_500000iter_2_noLearning_initMethod_4_lapSigma_10.00_seed2'},
%                 {'SAB2018/laplacianPolicy/18-08-13_500000iter_3_noLearning_initMethod_4_lapSigma_10.00_seed3'},
%                 {'SAB2018/laplacianPolicy/18-08-11_500000iter_4_noLearning_initMethod_4_lapSigma_10.00_seed4'},
% %             },{
% %                 {'SAB2018/laplacianPolicy/18-09-26_500000iter_2_noLearning_initMethod_4_lapSigma_15.00_seed2'},
% %                 {'SAB2018/laplacianPolicy/18-09-27_500000iter_3_noLearning_initMethod_4_lapSigma_15.00_seed3'},
% %                 {'SAB2018/laplacianPolicy/18-09-27_500000iter_4_noLearning_initMethod_4_lapSigma_15.00_seed4'},
%             },{
%                 {'SAB2018/laplacianPolicy/18-09-26_500000iter_2_noLearning_initMethod_4_lapSigma_20.00_seed2'},
%                 {'SAB2018/laplacianPolicy/18-09-27_500000iter_3_noLearning_initMethod_4_lapSigma_20.00_seed3'},
%                 {'SAB2018/laplacianPolicy/18-09-27_500000iter_4_noLearning_initMethod_4_lapSigma_20.00_seed4'},
%             }
%     };
%     % labelz = {'s = 0', 's = 0.2', 's = 0.5', 's = 1', 's = 2', 's = 5', 's = 10', 's = 15', 's = 20'};
%     labelz = {'s = 0', 's = 1', 's = 2', 's = 5', 's = 10', 's = 20'};


%     models = {
%         {   % normal case
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2'},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3'},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4'},
%         },{
%             % normal slant rotation distribution mean 0, std=10
%             {'eLifePaper/spreadBinocularity/18-11-17_500000iter_2_testSimV1_rotatePlane=-20,20deg_s2'},
%             {'eLifePaper/spreadBinocularity/18-11-19_500000iter_3_testSimV1_rotatePlane=-20,20deg_s3'},
%             {'eLifePaper/spreadBinocularity/18-11-17_500000iter_4_testSimV1_rotatePlane=-20,20deg_s4'},
%
% %         },{ % objSize = 1.5 m
% %             {'eLifePaper/spreadBinocularity/18-11-07_500000iter_2_objSize_1.5_patchSize_8_seed2'},
% %             {'eLifePaper/spreadBinocularity/18-11-07_500000iter_3_objSize_1.5_patchSize_8_seed3'}, % not completed
% %             {'eLifePaper/spreadBinocularity/18-11-07_500000iter_4_objSize_1.5_patchSize_8_seed4'}, % not completed
% %         },{ % objSize = 1 m
% %             {'eLifePaper/spreadBinocularity/18-11-08_500000iter_2_objSize_1.0_patchSize_8_seed2'},
% %             {'eLifePaper/spreadBinocularity/18-11-08_500000iter_3_objSize_1.0_patchSize_8_seed3'},
% % %             {'eLifePaper/spreadBinocularity/18-11-06_500000iter_2_test_objSize=1m'},   % please replace with another seed
%         },{   % increased slant rotation
%             {'eLifePaper/spreadBinocularity/18-11-30_500000iter_2_testSim2_rotatePlane_unif-80,80_s2'},
%             {'eLifePaper/spreadBinocularity/18-11-30_500000iter_3_testSimV2_rotatePlane=-80,80deg_s3'},
%             {'eLifePaper/spreadBinocularity/18-11-30_500000iter_4_testSimV2_rotatePlane=-80,80deg_s4'},
%         },{   % increased slant rotation + tilt
%             {'eLifePaper/spreadBinocularity/18-12-03_500000iter_2_yaw+tiltmaxYaw_80_maxTilt_45_seed2'},
%             {'eLifePaper/spreadBinocularity/18-12-03_500000iter_3_yaw+tiltmaxYaw_80_maxTilt_45_seed3'},
%             {'eLifePaper/spreadBinocularity/18-12-03_500000iter_4_yaw+tiltmaxYaw_80_maxTilt_45_seed4'},
%         }
%     };
%     labelz = {'no rotation', 'rotation: m=0, s=10', 'rotation: univ. [-80, 80]', 'rotation: univ. [-80, 80] +\ntilt: univ. [-45, 45]'};

%     models = {
%         {   % normal case
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2'},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3'},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4'},
%         },{
%             % reduced fixation time: 5 iterations
%             {'eLifePaper/spreadBinocularity/18-11-19_1000000iter_2_objS_3_interval_5_seed2'},
%             {'eLifePaper/spreadBinocularity/18-11-19_1000000iter_3_objS_3_interval_5_seed3'},
%             {'eLifePaper/spreadBinocularity/18-11-19_1000000iter_4_objS_3_interval_5_seed4'},
%
%         },{ % reduced fication time: 2 iterations
%             {'eLifePaper/spreadBinocularity/18-11-21_2500000iter_2_objS_3_interval_2_seed2'},
%             {'eLifePaper/spreadBinocularity/18-11-21_2500000iter_3_objS_3_interval_2_seed3'},
%             {'eLifePaper/spreadBinocularity/18-11-21_2500000iter_4_objS_3_interval_2_seed4'},
%         }
%     };
%     labelz = {'iter=10', 'iter=5', 'iter=2'};

%     models = {
%         {   % normal case
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2'},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3'},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4'},
%         },{ % 2 deg strab
%             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_2_seed2'},
%             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_2_seed3'},
%             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_2_seed4'},
%         },{ % 5 deg strab
%             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_5_seed2'},
%             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_5_seed3'},
%             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_5_seed4'},
%         },{ % 10 deg
%             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2'},
%             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3'},
%             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4'},
%         }
%     };
%     labelz = {'normal', 'strab. (2 deg)', 'strab. (5 deg)', 'strab. (10 deg)'};

%     models = {
%         {   % normal case
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2'},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3'},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4'},
%         },{ % 10 deg (adaptable vergence angle)
%             {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2'},
%             {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3'},
%             {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4'},
%         },{ % 10 deg (fixed vergence angle)
%             {'eLifePaper/strabism/18-12-17_500000iter_2_fixVergAngleAt6mfiltB_29_strabAngle_10_seed2'},
%             {'eLifePaper/strabism/18-12-17_500000iter_3_fixVergAngleAt6mfiltB_29_strabAngle_10_seed3'},
%             {'eLifePaper/strabism/18-12-17_500000iter_4_fixVergAngleAt6mfiltB_29_strabAngle_10_seed4'},
%         }
%     };
%     labelz = {'normal', 'strab. (10 deg)', 'fixed strab. (10 deg)'};

%     models = {
%         {   % strabismic case, fixed angle 10deg
%             {'eLifePaper/strabism/18-12-17_500000iter_2_fixVergAngleAt6mfiltB_29_strabAngle_10_seed2'},
%             {'eLifePaper/strabism/18-12-17_500000iter_3_fixVergAngleAt6mfiltB_29_strabAngle_10_seed3'},
%             {'eLifePaper/strabism/18-12-17_500000iter_4_fixVergAngleAt6mfiltB_29_strabAngle_10_seed4'},
%         },{ % fixed strabismus, only 2 BFs used for encoding
%             {'eLifePaper/strabism/19-01-07_500000iter_2_fixVergAngleAt6m_2BasisUsedfiltB_29_strabAngle_10_seed2'},
%             {'eLifePaper/strabism/19-01-07_500000iter_3_fixVergAngleAt6m_2BasisUsedfiltB_29_strabAngle_10_seed3'},
%             {'eLifePaper/strabism/19-01-07_500000iter_4_fixVergAngleAt6m_2BasisUsedfiltB_29_strabAngle_10_seed4'},
%         },{ % fixed strabismus, random BFs initialization
%             {'eLifePaper/strabism/19-01-07_500000iter_2_fixVergAngleAt6m_2BasisUsed_filtB_29_strabAngle_10_seed2'},
%             {'eLifePaper/strabism/19-01-07_500000iter_3_fixVergAngleAt6m_2BasisUsed_filtB_29_strabAngle_10_seed3'},
%             {'eLifePaper/strabism/19-01-07_500000iter_4_fixVergAngleAt6m_2BasisUsed_filtB_29_strabAngle_10_seed4'},
%         },{ % fixed strabismus increased strab angle
%             {'eLifePaper/strabism/19-01-13_500000iter_2_fixAt1m_objRange05-2filtB_29_strabAngle_15_seed2'},
%             {'eLifePaper/strabism/19-01-13_500000iter_3_fixAt1m_objRange05-2filtB_29_strabAngle_15_seed3'},
%             {'eLifePaper/strabism/19-01-13_500000iter_4_fixAt1m_objRange05-2filtB_29_strabAngle_15_seed4'},
%         }
%     };
%     labelz = {'strab. (10 deg)', 'strab. (10 deg), 2 BFs used', 'fixed strab. (10 deg), random BFs init', 'strab. (15 deg)'};

%     models = {
%                {   % normal case
%             {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1'},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2'},
%             {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3'},
%             {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4'},
%             {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5'},
%         },{ % no RL
%             {'eLifePaper/vergenceInfluence/19-05-15_500000iter_1_noLearning_initMethod_2_lapSigma_0_seed1'},
%             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_2_noLearninginitMethod_2_lapSigma_0_seed2'},
%             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_3_noLearninginitMethod_2_lapSigma_0_seed3'},
%             {'eLifePaper/vergenceInfluence/19-04-27_500000iter_4_noLearninginitMethod_2_lapSigma_0_seed4'},
%             {'eLifePaper/vergenceInfluence/19-05-15_500000iter_5_noLearning_initMethod_2_lapSigma_0_seed5'},
%         },{ % no RL & no Disparities
%             {'eLifePaper/vergenceInfluence/19-05-16_500000iter_1_noLearning_initMethod_4_lapSigma_0_seed1'},
%             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_2_noLearning_initMethod_4_lapSigma_0_seed2'},
%             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_3_noLearning_initMethod_4_lapSigma_0_seed3'},
%             {'eLifePaper/vergenceInfluence/19-04-29_500000iter_4_noLearning_initMethod_4_lapSigma_0_seed4'},
%             {'eLifePaper/vergenceInfluence/19-05-16_500000iter_5_noLearning_initMethod_4_lapSigma_0_seed5'},
%         }
%     }
%     labelz = {'normal', 'no RL', 'no RL & no disp.'};

    colors = [[77,175,74]; [55,126,184]; [152,78,163]; [255,127,0]; [228,26,28]; [166,86,40]; [118,42,131]]./256;
    
    %% first plot: compare all coarse scales
%     h = figure('position', [300 400 1400 400]);
    h = figure('position', [300 400 1400 800]);
%     set(gcf, 'defaulttextinterpreter','latex')
    
    means = zeros(2, length(labelz), nbins);
    stds = zeros(2, length(labelz), nbins);
    nUsed = zeros(2, length(labelz));
    
    for scale = 1 : 2

%         binocularityVals = {};
        
        for m = 1:length(models)
            nModels = length(models{m});
            if length(bins) == 1
                binocs = zeros(nModels, bins);
            else
                binocs = zeros(nModels, length(bins));
            end
            belowThresh = zeros(nModels, 1);
%             binVals = [];
            
            for i = 1 : nModels
                model = load(strcat(filePath, models{m}{i}{1}, '/model.mat'));
                model = model.model;

                % comparing the energy in left and right eye
                if strcmp(method, 'energy')
                    leftBasis = model.scModel{scale}.basis(1:end/2, :);
                    rightBasis = model.scModel{scale}.basis(end/2+1:end, :);

        %             binocularity = (sqrt(sum(leftBasis.^2)) - sqrt(sum(rightBasis.^2))) ./ (sqrt(sum(leftBasis.^2)) + sqrt(sum(rightBasis.^2)));
        %             binocularity = -binocularity;
                    binocularity = (sqrt(sum(rightBasis.^2)) - sqrt(sum(leftBasis.^2))) ./ (sqrt(sum(leftBasis.^2)) + sqrt(sum(rightBasis.^2)));
                    binocularity = squeeze(binocularity(1, :));

                % comparing responses to dominant stimuli
                elseif strcmp(method, 'responses')
                    Bases = model.scModel{scale}.basis;
                    if ~useFreq
                        paraL = load(strcat(model.savePath, '/modelAt500000/scModel', num2str(scale), "eyes", num2str(2), '.mat')); % 6 obsolete
                        paraR = load(strcat(model.savePath, '/modelAt500000/scModel', num2str(scale), "eyes", num2str(3), '.mat')); % 5 obsolete
                    else
                        paraL = load(strcat(model.savePath, '/modelAt500000/scModel', num2str(scale), "eyes", num2str(2), '_freq.mat'));
                        paraR = load(strcat(model.savePath, '/modelAt500000/scModel', num2str(scale), "eyes", num2str(3), '_freq.mat'));
                    end
            %         remove unnecessary entries
                    paraL.Fitted_Gabor_Filter_Parameters(:,6) = [];
                    paraR.Fitted_Gabor_Filter_Parameters(:,6) = [];

                    binocularity = [];
                    for b = 1:size(Bases, 2)
                        if ~useFreq
                            stimL = reshape(gabor(paraL.Fitted_Gabor_Filter_Parameters(b, :), 4), [1, 64]);
                            stimR = reshape(gabor(paraR.Fitted_Gabor_Filter_Parameters(b, :), 4), [1, 64]);
                        else
                            stimL = reshape(gaborFreq(paraL.Fitted_Gabor_Filter_Parameters(b, :), 4), [1, 64]);
                            stimR = reshape(gaborFreq(paraR.Fitted_Gabor_Filter_Parameters(b, :), 4), [1, 64]);
                        end
                        errL = paraL.Error(b);
                        errR = paraR.Error(b);
                        
                        if errL > thresh || errR > thresh
%                             sprintf('Base number %d is above threashold', b)
                            continue
                        end

                        basisL = Bases(1:end/2, b);
                        basisR = Bases(end/2+1:end, b);
                        
                        % normalize our Gabors
%                         stimL = bsxfun(@rdivide, stimL, sqrt(sum(stimL .^ 2)));
%                         stimR = bsxfun(@rdivide, stimR, sqrt(sum(stimR .^ 2)));

                        maxRespL = abs(stimL * basisL);
                        maxRespR = abs(stimR * basisR);

                        if maxRespL > maxRespR
            %                 respL = sum(sum(basisL .* stimL));
            %                 respR = sum(sum(basisR .* stimL));
                            respL = abs(stimL * basisL);
                            respR = abs(stimL * basisR);
%                             binocularity(end+1) = (respR - respL)/(respR + respL);
                        else
                            respL = abs(stimR * basisL);
                            respR = abs(stimR * basisR);
%                             binocularity(end+1) = (respR - respL)/(respR + respL);
                        end
                        binocularity(end+1) = (respR - respL)/(respR + respL);
%                         binocularity(end+1) = (respL - respR)/(respR + respL); % reverse sign
                    end
                end

%                 bins = [-1, -0.85, -0.5, -0.15, 0.15, 0.5, 0.85, 1];
%                 bins = [-1 : 2/7 : 1]; % try equidistant bins
%                 bins = [-1 : 2/21 : 1];
%                   binVals = [binVals; binocularity'];
                  [N, ~] = histcounts(binocularity, bins);
                  binocs(i, :) = binocs(i, :) + ((N./sum(N))*100);  % measure in percentage
                  belowThresh(i) = sum(N);
            end

            meanBinocs = mean(binocs, 1);
            stdBinocs = std(binocs, 1);
            
            means(scale, m, :) = meanBinocs;
            stds(scale, m, :) = stdBinocs;
            nUsed(scale, m) = mean(belowThresh);
            
%             binocularityVals{m} = binVals;
            
            % plot all in separate plots
            if scale == 1
                subplot(2, length(models), m);
                title(sprintf(labelz{m}));
            else
                subplot(2, length(models), length(models)+m);
            end
            hold on;
            bar(meanBinocs, 0.8);
            errorbar(meanBinocs, stdBinocs, '.black', 'CapSize', 4);
            grid minor;
            xlabel('Binocularity');
            if m == 1
                ylabel('Number of Bases [%]');
            end
            
            if scale == 1
                ylim([0 100]);
            elseif scale == 2
                ylim([0 70]);
            end
            xlim([0.5, bins(end)+0.5]);
            xticks([1, ceil(nbins/2), nbins]);
            xticklabels({'-1', '0', '1'});
%             xlim([0.5, 21.5]);
            set(gca, 'TickLabelInterpreter', 'latex')
        end
    end
    if ~isempty(saveName)
        saveas(h, sprintf('%ssingleBinocs_%s.png', savePath, saveName));
    %     saveas(h, sprintf('%scompareMoreBinocs_%s.png', savePath, saveName));
    end
    
    nCols = 2;
    if nCols == 1
        figure('position', [300 400 1400 300]); % single column
    elseif nCols == 2
        figure('position', [300 350 1000 625]); % 2 columns
    end
    
    for m = 1 : length(models)
        subplot(nCols, ceil(length(models)/nCols), m);
        
        me = reshape(means(:, m, :), [2, nbins])';
        st = reshape(stds(:, m, :), [2, nbins])';
        % bar(me, 'grouped')
        
%         colrs = {[ 0, 204, 0]./256, [ 255, 153, 0]./256}; %green and yellow
        
        colrs = {[27,158,119]./256, [ 217,95,2]./256}; %green and orange, color-blind safe
%         colrs = {[141,160,203]./356, [252,141,98]./356}; % purple and orange
        
        h = bar(me, 'hist', 'BarWidth', 0.6);
        
%         h(1).FaceColor = colrs{1};%'green'; 
%         h(2).FaceColor = colrs{2};%'yellow'; 
        
        % consistend color labeling version
        h(1).FaceColor = colors(m,:);
%         h(1).EdgeColor = colors(m,:);
        h(2).FaceColor = colors(m,:);
        h(2).FaceAlpha = 0.5;
%         h(2).FaceColor = [1,1,1];
%         h(2).EdgeColor = colors(m,:);
        
        hold on;
        title(labelz{m});
        if m < 4
            ylim([-0.4, 90]);
        else
            ylim([-0.2, 90]);
        end
        
        yticks([0, 30, 60, 90]);
        
        xticks([1,round(nbins/2), nbins]);
        xticklabels({-1,0,1});
        xlim([0.5, 7.5]);
        if m == 1 || m == 4
%             ylabel('Number of Bases [%]');
            ylabel('number of bases [%]');
        end
        if m == 5
%             xlabel('Binocularity Index');
            xlabel('binocularity index');
        end
%         if nCols == 1
%             grid on;
%         elseif nCols == 2
%             grid on;
%         end
        
        if m == 1%length(models)
            lg = legend('coarse', 'fine', 'AutoUpdate', 'off');
            legend('Location', 'northwest');
%             legend('Location', 'best');
%             legend('Location', 'bestoutside');
            legend('boxoff');
            lg.FontSize = 8;
            if nCols == 1
                lg.Position = [0.16    0.85    0.0    0.05];
            elseif nCols == 2
%                 lg.Position = [0.16    0.85    0.0    0.05];
            end
        end
        hold on;
        
        ngroups = nbins;
        nbars = 2;
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, me(:, i), st(:, i), 'k', 'linestyle', 'none', 'CapSize', 3);
%             errorbar(x, me(:, i), st(:, i), 'linestyle', 'none', 'CapSize', 3, 'Color', colors(m,:));
        end
        
        
%         if m < 4
%             t = text(4.3, 85, strcat(['$\bar{N}_c$=',num2str(nUsed(1, m), '%.0f')]));
%             set(t,'interpreter','Latex','FontSize',12);
%             t = text(4.3, 80, strcat(['$\bar{N}_f$=',num2str(nUsed(2, m), '%.0f')]));
%             set(t,'interpreter','Latex','FontSize',12);
%             t = text(4.6, 87.3, strcat(['$\bar{N}_c$=',num2str(nUsed(1, m), '%.0f')]));
%             set(t,'interpreter','Latex','FontSize',12);
%             t = text(4.6, 82.3, strcat(['$\bar{N}_f$=',num2str(nUsed(2, m), '%.0f')]));
%             set(t,'interpreter','Latex','FontSize',12);
        if nCols == 1
            t = text(4.6, 85.5, strcat(['$\bar{N}_c\mathsf{=',num2str(nUsed(1, m), '%.0f}$')]));
            set(t,'interpreter','Latex','FontSize',11);
            t = text(4.6, 80.4, strcat(['$\bar{N}_f\mathsf{=',num2str(nUsed(2, m), '%.0f}$')]));
            set(t,'interpreter','Latex','FontSize',11);
        elseif nCols == 2
            t = text(5, 81.5, strcat(['$\bar{N}_c\mathsf{=',num2str(nUsed(1, m), '%.0f}$')]));
            set(t,'interpreter','Latex','FontSize',11);
            t = text(5, 74.5, strcat(['$\bar{N}_f\mathsf{=',num2str(nUsed(2, m), '%.0f}$')]));
            set(t,'interpreter','Latex','FontSize',11);
        end
%         else
%             t = text(2, 80, strcat(['$\bar{N}_c$=',num2str(nUsed(1, m), '%.0f')]));
%             set(t,'interpreter','Latex','FontSize',12);
%             t = text(2, 75, strcat(['$\bar{N}_f$=',num2str(nUsed(2, m), '%.0f')]));
%             set(t,'interpreter','Latex','FontSize',12);
%         end
        
    end
    if ~isempty(saveName)
        saveas(gcf, sprintf('%scombinedBinocs_%s.png', savePath, saveName));
    %     saveas(h, sprintf('%scompareMoreBinocs_%s.png', savePath, saveName));
    end
end

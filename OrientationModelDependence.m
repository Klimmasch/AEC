% Script for analysing orientations from different models. 

path = '/home/aecgroup/aecdata/Results/SAB2018/';

threshold = 0.2;
orientation = 1; % first bin == vertical, bins =  0:15:165

% usage: {{model1s1, model1s2, ...},{model2s1, model2s2,..}, ...}
models = {
    %%%% Plot: 0deg-MD-wDisp-wLearning %%%
    
%     {
%     {'/normalInput/18-03-14_500000iter_2_normal2', 2, 2},
%     {'/normalInput/18-03-14_500000iter_3_normal3', 2, 2},
%     {'/normalInput/18-03-15_500000iter_4_normal4', 2, 2},
%     {'/normalInput/18-03-10_500000iter_5_normal5', 2, 2},
%     {'/normalInput/18-03-10_500000iter_6_normal6', 2, 2}
%     },
%    
%     {
%     {'/monocularDeprivation_diffProbs_local_s/18-03-14_500000iter_2_monDep_filtR3_prob01_s2', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-12_500000iter_3_monDep_filtR3_prob01', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-14_500000iter_4_monDep_filtR3_prob01_s4', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_5_monDep_filtR3_prob01_s5', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_6_monDep_filtR3_prob01_s6', 2, 2}
%     },
%     
%     {
%     {'/monocularDeprivation_diffProbs_local_s/18-03-10_500000iter_2_monDep_filtR3_prob025', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-15_500000iter_3_monDep_filtR3_prob025_s3', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-14_500000iter_4_monDep_filtR3_prob025_s4', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_5_monDep_filtR3_prob025_s5', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-17_500000iter_6_monDep_filtR3_prob025_s6', 2, 2}
%     },
%     
%     {
%     {'/monocularDeprivation_diffProbs_local_s/18-03-14_500000iter_2_monDep_filtR3_prob05_s2', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-12_500000iter_3_monDep_filtR3_prob05', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-15_500000iter_4_monDep_filtR3_prob05_s4', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_5_monDep_filtR3_prob05_s5', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-17_500000iter_6_monDep_filtR3_prob05_s6', 2, 2}
%     },
%     
%     {
%     {'/monocularDeprivation_diffProbs_local_s/18-03-10_500000iter_2_monDep_filtR3_prob075', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_3_monDep_filtR3_prob075_s3', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-15_500000iter_4_monDep_filtR3_prob075_s4', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-17_500000iter_5_monDep_filtR3_prob075_s5', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-18_500000iter_6_monDep_filtR3_prob075_s6', 2, 2}
%     },
%     
%     {
%     {'/monocularDeprivation_diffProbs_local_s/18-03-11_500000iter_2_monDep_filtR3_prob09', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_3_monDep_filtR3_prob09_s3', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-15_500000iter_4_monDep_filtR3_prob09_s4', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-17_500000iter_5_monDep_filtR3_prob09_s5', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-18_500000iter_6_monDep_filtR3_prob09_s6', 2, 2}
%     },
%     
%     {
%     {'/monocularDeprivation_diffProbs_local_s/18-03-15_500000iter_2_monDep_filtR3_prob1_s2', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-13_500000iter_3_monDep_filtR3_prob1', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_4_monDep_filtR3_prob1_s4', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-18_500000iter_5_monDep_filtR3_prob1_s5', 2, 2},
%     {'/monocularDeprivation_diffProbs_local_s/18-03-18_500000iter_6_monDep_filtR3_prob1_s6', 2, 2}
%     }
    
    %%% Plot: 0deg-MD-wDisp-w/oLearning %%%
    
%     {
%     {'/compareOrientedInputw-oVergence/18-05-25_500000iter_16_ActorLR=0_s16', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-05-29_500000iter_17_ActorLR=0_s17', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-05-30_500000iter_18_ActorLR=0_s18', 2, 2}
%     },
%     
%     {
%     {'/compareOrientedInputw-oVergence/18-05-28_500000iter_2_noVerg_monDep_filtR3_prob01_s2', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-01_500000iter_3_noVerg_monDep_filtR3_prob01_s3', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-05_500000iter_4_noVerg_monDep_filtR3_prob01_s4', 2, 2}
%     },
%     
%     {
%     {'/compareOrientedInputw-oVergence/18-05-28_500000iter_2_noVerg_monDep_filtR3_prob025_s2', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-02_500000iter_3_noVerg_monDep_filtR3_prob025_s3', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-06_500000iter_4_noVerg_monDep_filtR3_prob025_s4', 2, 2}
%     },
%    
%     {
%     {'/compareOrientedInputw-oVergence/18-05-27_500000iter_2_noVerg_monDep_filtR3_prob05_s2', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-05-31_500000iter_3_noVerg_monDep_filtR3_prob05_s3', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-04_500000iter_4_noVerg_monDep_filtR3_prob05_s4', 2, 2}
%     },
%    
%     {
%     {'/compareOrientedInputw-oVergence/18-05-29_500000iter_2_noVerg_monDep_filtR3_prob075_s2', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-02_500000iter_3_noVerg_monDep_filtR3_prob075_s3', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-06_500000iter_4_noVerg_monDep_filtR3_prob075_s4', 2, 2}
%     },
%     
%     {
%     {'/compareOrientedInputw-oVergence/18-05-28_500000iter_2_noVerg_monDep_filtR3_prob1_s2', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-01_500000iter_3_noVerg_monDep_filtR3_prob1_s3', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-05_500000iter_4_noVerg_monDep_filtR3_prob1_s4', 2, 2}
%     }
    
    %%% Plot: 0deg-MD-w/oDisp-w/oLearning %%%
    
%     {
%     {'/compareOrientedInputw-oVergence/18-05-25_500000iter_16_ActorLR=0_OD=3m_0disp_s16', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-05-29_500000iter_17_ActorLR=0_OD=3m_0disp_s17', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-05-30_500000iter_18_ActorLR=0_OD=3m_0disp_s18', 2, 2}
%     },
%     
%     {
%     {'/compareOrientedInputw-oVergence/18-05-26_500000iter_2_0disp_monDep_filtR3_prob01_s2', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-05-30_500000iter_3_0disp_monDep_filtR3_prob01_s3', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-03_500000iter_4_0disp_monDep_filtR3_prob01_s4', 2, 2}
%     },
%     
%     {
%     {'/compareOrientedInputw-oVergence/18-05-26_500000iter_2_0disp_monDep_filtR3_prob025_s2', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-05-30_500000iter_3_0disp_monDep_filtR3_prob025_s3', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-04_500000iter_4_0disp_monDep_filtR3_prob025_s4', 2, 2}
%     },
%    
%     {
%     {'/compareOrientedInputw-oVergence/18-05-25_500000iter_2_0disp_monDep_filtR3_prob05_s2', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-05-29_500000iter_3_0disp_monDep_filtR3_prob05_s3', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-02_500000iter_4_0disp_monDep_filtR3_prob05_s4', 2, 2}
%     },
%    
%     {
%     {'/compareOrientedInputw-oVergence/18-05-27_500000iter_2_0disp_monDep_filtR3_prob075_s2', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-05-31_500000iter_3_0disp_monDep_filtR3_prob075_s3', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-04_500000iter_4_0disp_monDep_filtR3_prob075_s4', 2, 2}
%     },
%     
%     {
%     {'/compareOrientedInputw-oVergence/18-05-26_500000iter_2_0disp_monDep_filtR3_prob1_s2', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-05-30_500000iter_3_0disp_monDep_filtR3_prob1_s3', 2, 2},
%     {'/compareOrientedInputw-oVergence/18-06-03_500000iter_4_0disp_monDep_filtR3_prob1_s4', 2, 2}
%     }
    
%     %%% Plot: normal condition comparison of wDisp wLearning, wDisp w/oLearning
%     %%% and w/oDisp w/oLearning
%     
%     {
%     {'/normalInput/18-03-16_500000iter_2_normal_gf33-01-01_2', 2, 1},
%     {'/normalInput/18-03-16_500000iter_3_normal_gf33-01-01_3', 2, 1},
%     {'/normalInput/18-03-17_500000iter_4_normal_gf33-01-01_4', 2, 1},
%     {'/normalInput/18-03-18_500000iter_5_normal_gf33-01-01_5', 2, 1},
%     {'/normalInput/18-03-18_500000iter_6_normal_gf33-01-01_6', 2, 1}
%     },
%     
%     {
%     {'/compareOrientedInputw-oVergence/18-05-25_500000iter_16_ActorLR=0_s16', 2, 1},
%     {'/compareOrientedInputw-oVergence/18-05-29_500000iter_17_ActorLR=0_s17', 2, 1},
%     {'/compareOrientedInputw-oVergence/18-05-30_500000iter_18_ActorLR=0_s18', 2, 1}
%     },
%     
%     {
%     {'/compareOrientedInputw-oVergence/18-05-25_500000iter_16_ActorLR=0_OD=3m_0disp_s16', 2, 1},
%     {'/compareOrientedInputw-oVergence/18-05-29_500000iter_17_ActorLR=0_OD=3m_0disp_s17', 2, 1},
%     {'/compareOrientedInputw-oVergence/18-05-30_500000iter_18_ActorLR=0_OD=3m_0disp_s18', 2, 1}
%     },
    

    %%% Plot: horizontal condition comparison of wDisp wlearning, wDisp w/oLearning
    %%% and w/oDisp w/oLearning
    
%     {
%     {'/compareOrientedInputw-oVergence/18-06-02_500000iter_16_horOnly_+vergence_s16', 2, 1},
%     {'/compareOrientedInputw-oVergence/18-06-03_500000iter_17_horOnly_+vergence_s17', 2, 1},
%     {'/compareOrientedInputw-oVergence/18-06-03_500000iter_18_horOnly_+vergence_s18', 2, 1}
%     },
    
%     {
%     {'/compareOrientedInputw-oVergence/18-06-01_500000iter_16_horOnly_ActorLR=0_s16', 2, 1},
%     {'/compareOrientedInputw-oVergence/18-06-02_500000iter_17_horOnly_ActorLR=0_s17', 2, 1},
%     {'/compareOrientedInputw-oVergence/18-06-03_500000iter_18_horOnly_ActorLR=0_s18', 2, 1}
%     },
    
%     {
%     {'/compareOrientedInputw-oVergence/18-06-01_500000iter_16_horOnly_ActorLR=0_OD=3m_0disp_s16', 2, 1},
%     {'/compareOrientedInputw-oVergence/18-06-02_500000iter_17_horOnly_ActorLR=0_OD=3m_0disp_s17', 2, 1},
%     {'/compareOrientedInputw-oVergence/18-06-03_500000iter_18_horOnly_ActorLR=0_OD=3m_0disp_s18', 2, 1}
%     },
%     
    
    }

bins = 1:size(models);
Ns = zeros(size(models, 1), 1); % mean relative nums
Ns_err = zeros(size(models, 1), 1); % std relative nums
bf_nums = zeros(size(models, 1), 1); % num of fits below threshold
ttest_sample1 = models{size(models, 1)};
ttest_sample2 = models{size(models, 1)};

for i = 1:size(models)
    tmp_N = 1:size(models{i}); % relative nums for each model
    for j = 1:length(tmp_N)
        para = load(strcat(path, models{i}{j}{1}, '/scModel', num2str(models{i}{j}{2}), "eyes", num2str(models{i}{j}{3}), '.mat'));
        Resnorm_Set = para.Error;
        Parameter_Set = para.Fitted_Gabor_Filter_Parameters;
        idx = find(Resnorm_Set < threshold);
        
        if (j == 2) % s4 num of used fits
            bf_nums(i) = length(idx);
        end

        bins_h = -7.5:15:172.5;
       
        N = histcounts(mod(Parameter_Set(idx,2)*180/pi+7.5, 180)-7.5,bins_h);
        tmp_N(j) = N(orientation)/sum(N);
    end
    
    if (i == 1)
        ttest_sample1 = tmp_N;
    end
    
    if (i == size(models, 1))
        ttest_sample2 = tmp_N;
    end
    
    Ns(i) = mean(tmp_N);
    Ns_err(i) = std(tmp_N);
end

hold on;
bar(bins, Ns*100, 1);
errorbar(bins, Ns*100, Ns_err*100, '.');
line([1,3],[43,43], "Color", "black")
line([1,1],[42,44], "Color", "black")
line([3,3],[42,44], "Color", "black")
text(3.95, 44, "*", "FontSize", 12)
hold off;

grid on
%xlabel('Normal condition')
ylabel('Percentage of Bases [%]')
xticks([1, 2, 3])%, 4, 5, 6, 7])
xticklabels({"D+L+", "D+L-", "D-L-"})%, 50, 75, 90, 100})
xlim([0.5, 3.5])
ylim([0 50])

title({"Vertical Orienation - Horizontal Condition","w/ Disparities w/ Learning","w/ Disparities w/o Learning","w/o Disparities w/o Learning"}, "FontSize", 10)

% t = text(1, 5, strcat("N_{s4}=", num2str(bf_nums(1))), 'FontSize', 12,'fontWeight','bold');
%     set(t,'Rotation',90);

[h,p] = ttest(ttest_sample1, ttest_sample2);

text(3, 45, strcat("p=", num2str(round(p, 3))),'FontSize', 10,'fontWeight','bold')

set(gca,'FontSize',15,'fontWeight','bold') %,'FontName','Courier')
%set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold') %,'FontName','Courier')

saveas(gcf, '0deg-HC-wDwL-wDwoL-woDwoL.png','png')

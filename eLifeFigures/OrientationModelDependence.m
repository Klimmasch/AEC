% Script for analysing orientations from different models.
% @param    subPlt can either be 'A' or 'B'
function OrientationModelDependence(saveTag)

    % filePath = '/home/aecgroup/aecdata/Results/SAB2018/';
    % filePath = '/home/mrklim/fiasAECGroup/aecdata/Results/SAB2018/';
    filePath = '/home/aecgroup/aecdata/Results/eLifePaper/';
    savePath = '/home/aecgroup/aecdata/Results/eLifePaper/plots/VergInfl/'

    threshold = 0.2;
    orientation = 1; % first bin == vertical, bins =  0:15:165

    for subPlt = 1:2
        if subPlt == 1
            scale = 2;  % choose between 1 and 2, exept for Laplacian Case where there is only one scale
            eyes = 2   % choose between 1 (both), 2 (left), 3 (right)
            % usage: {{model1s1, model1s2, ...},{model2s1, model2s2,..}, ...}
            models = {



                 %% Influence of vergence learning
                {   % normal case
                    {'explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1', scale, eyes},
                    {'explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2', scale, eyes},
                    {'explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3', scale, eyes},
                    {'explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4', scale, eyes},
                    {'explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5', scale, eyes},
        %             {'explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2', 2, eyes},
        %             {'explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3', 2, eyes},
        %             {'explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4', 2, eyes},
                },{ % no RL
                    {'vergenceInfluence/19-05-15_500000iter_1_noLearning_initMethod_2_lapSigma_0_seed1', scale, eyes},
                    {'vergenceInfluence/19-04-27_500000iter_2_noLearninginitMethod_2_lapSigma_0_seed2', scale, eyes},
                    {'vergenceInfluence/19-04-27_500000iter_3_noLearninginitMethod_2_lapSigma_0_seed3', scale, eyes},
                    {'vergenceInfluence/19-04-27_500000iter_4_noLearninginitMethod_2_lapSigma_0_seed4', scale, eyes},
                    {'vergenceInfluence/19-05-15_500000iter_5_noLearning_initMethod_2_lapSigma_0_seed5', scale, eyes},
        %             {'vergenceInfluence/19-04-27_500000iter_2_noLearninginitMethod_2_lapSigma_0_seed2', 2, eyes},
        %             {'vergenceInfluence/19-04-27_500000iter_3_noLearninginitMethod_2_lapSigma_0_seed3', 2, eyes},
        %             {'vergenceInfluence/19-04-27_500000iter_4_noLearninginitMethod_2_lapSigma_0_seed4', 2, eyes},
                },{ % no RL & no Disparities
                    {'vergenceInfluence/19-05-16_500000iter_1_noLearning_initMethod_4_lapSigma_0_seed1', scale, eyes},
                    {'vergenceInfluence/19-04-29_500000iter_2_noLearning_initMethod_4_lapSigma_0_seed2', scale, eyes},
                    {'vergenceInfluence/19-04-29_500000iter_3_noLearning_initMethod_4_lapSigma_0_seed3', scale, eyes},
                    {'vergenceInfluence/19-04-29_500000iter_4_noLearning_initMethod_4_lapSigma_0_seed4', scale, eyes},
                    {'vergenceInfluence/19-05-16_500000iter_5_noLearning_initMethod_4_lapSigma_0_seed5', scale, eyes},
        %             {'vergenceInfluence/19-04-29_500000iter_2_noLearning_initMethod_4_lapSigma_0_seed2', 2, eyes},
        %             {'vergenceInfluence/19-04-29_500000iter_3_noLearning_initMethod_4_lapSigma_0_seed3', 2, eyes},
        %             {'vergenceInfluence/19-04-29_500000iter_4_noLearning_initMethod_4_lapSigma_0_seed4', 2, eyes},
                }

                };

    %         names = {'normal case', 'no RL', {'no RL','no disp'}};
    %         names = {'normal case', 'no RL', 'no RL\nno disp'};
    %         names = {'normal case', 'no RL', '\begin{tabular}{c}no RL\\no disp\end{tabular}'};
%             names = {'normal', 'random disp.', 'zero disp.'};
            names = {'normal', 'random disp.', 'zero disp.'};
            names = cellfun(@(x) strrep(x,' ','\newline'), names,'UniformOutput',false);

            % nc baseline = 33.4

        elseif subPlt == 2 % Laplacian plot
            scale = 1;  % choose between 1 and 2, exept for Laplacian Case where there is only one scale
            eyes = 2;   % choose between 1 (both), 2 (left), 3 (right)
            % usage: {{model1s1, model1s2, ...},{model2s1, model2s2,..}, ...}
            models = {
                %%% Plot: Laplacian Disparity distribution, varying spread
            %     {
            %     {'/laplacianPolicy/18-07-26_500000iter_2_noLearning_initMethod_4_lapSigma_0.00_seed2', scale, eyes},
            %     {'/laplacianPolicy/18-07-23_500000iter_3_noLearning_initMethod_4_lapSigma_0.00_seed3', scale, eyes},
            %     {'/laplacianPolicy/18-07-21_500000iter_4_noLearning_initMethod_4_lapSigma_0.00_seed4', scale, eyes},
            %     },
            %     {
            %     {'/laplacianPolicy/18-07-26_500000iter_2_noLearning_initMethod_4_lapSigma_0.22_seed2', scale, eyes},
            %     {'/laplacianPolicy/18-07-23_500000iter_3_noLearning_initMethod_4_lapSigma_0.22_seed3', scale, eyes},
            %     {'/laplacianPolicy/18-07-21_500000iter_4_noLearning_initMethod_4_lapSigma_0.22_seed4', scale, eyes},
            %     },
            %     {
            %     {'/laplacianPolicy/18-07-26_500000iter_2_noLearning_initMethod_4_lapSigma_0.50_seed2', scale, eyes},
            %     {'/laplacianPolicy/18-07-23_500000iter_3_noLearning_initMethod_4_lapSigma_0.50_seed3', scale, eyes},
            %     {'/laplacianPolicy/18-07-21_500000iter_4_noLearning_initMethod_4_lapSigma_0.50_seed4', scale, eyes},
            %     },
            %     {
            %     {'/laplacianPolicy/18-07-26_500000iter_2_noLearning_initMethod_4_lapSigma_1.00_seed2', scale, eyes},
            %     {'/laplacianPolicy/18-07-23_500000iter_3_noLearning_initMethod_4_lapSigma_1.00_seed3', scale, eyes},
            %     {'/laplacianPolicy/18-07-21_500000iter_4_noLearning_initMethod_4_lapSigma_1.00_seed4', scale, eyes},
            %     },
            %     {
            %     {'/laplacianPolicy/18-07-26_500000iter_2_noLearning_initMethod_4_lapSigma_2.00_seed2', scale, eyes},
            %     {'/laplacianPolicy/18-07-23_500000iter_3_noLearning_initMethod_4_lapSigma_2.00_seed3', scale, eyes},
            %     {'/laplacianPolicy/18-07-21_500000iter_4_noLearning_initMethod_4_lapSigma_2.00_seed4', scale, eyes},
            %     }

            %% single (fine) scale model trained on Laplacian Policy
                    {   % elife
                        {'laplacianPolicy/19-06-06_500000iter_1_noLearn_fineScOnly_initMethod_4_lapSigma_0.00_seed1', 1, eyes},
                        {'laplacianPolicy/18-11-21_500000iter_2_noLearn_FineScOnly_initMethod_4_lapSigma_0.00_seed2', 1, eyes},
                        {'laplacianPolicy/18-11-22_500000iter_3_noLearn_FineScOnly_initMethod_4_lapSigma_0.00_seed3', 1, eyes},
                        {'laplacianPolicy/18-11-26_500000iter_4_noLearn_FineScOnly_initMethod_4_lapSigma_0.00_seed4', 1, eyes},
                        {'laplacianPolicy/19-06-06_500000iter_5_noLearn_fineScOnly_initMethod_4_lapSigma_0.00_seed5', 1, eyes},
                    },{%
                        {'laplacianPolicy/19-06-06_500000iter_1_noLearn_fineScOnly_initMethod_4_lapSigma_0.22_seed1', 1, eyes},
                        {'laplacianPolicy/18-11-21_500000iter_2_noLearn_FineScOnly_initMethod_4_lapSigma_0.22_seed2', 1, eyes},
                        {'laplacianPolicy/18-11-22_500000iter_3_noLearn_FineScOnly_initMethod_4_lapSigma_0.22_seed3', 1, eyes},
                        {'laplacianPolicy/18-11-26_500000iter_4_noLearn_FineScOnly_initMethod_4_lapSigma_0.22_seed4', 1, eyes},
                        {'laplacianPolicy/19-06-06_500000iter_5_noLearn_fineScOnly_initMethod_4_lapSigma_0.22_seed5', 1, eyes},
                    },{
                        {'laplacianPolicy/19-06-10_500000iter_1_noLearn_fineScOnly_initMethod_4_lapSigma_1.00_seed1', 1, eyes},
                        {'laplacianPolicy/18-11-21_500000iter_2_noLearn_FineScOnly_initMethod_4_lapSigma_1.00_seed2', 1, eyes},
                        {'laplacianPolicy/18-11-22_500000iter_3_noLearn_FineScOnly_initMethod_4_lapSigma_1.00_seed3', 1, eyes},
                        {'laplacianPolicy/18-11-26_500000iter_4_noLearn_FineScOnly_initMethod_4_lapSigma_1.00_seed4', 1, eyes},
                        {'laplacianPolicy/19-06-07_500000iter_5_noLearn_fineScOnly_initMethod_4_lapSigma_1.00_seed5', 1, eyes},
                    },{
                        {'laplacianPolicy/19-06-17_500000iter_1_noLearn_fineScOnly_initMethod_4_lapSigma_2.00_seed1', 1, eyes},
                        {'laplacianPolicy/18-11-21_500000iter_2_noLearn_FineScOnly_initMethod_4_lapSigma_2.00_seed2', 1, eyes},
                        {'laplacianPolicy/18-11-22_500000iter_3_noLearn_FineScOnly_initMethod_4_lapSigma_2.00_seed3', 1, eyes},
                        {'laplacianPolicy/18-11-26_500000iter_4_noLearn_FineScOnly_initMethod_4_lapSigma_2.00_seed4', 1, eyes},
                        {'laplacianPolicy/19-06-17_500000iter_5_noLearn_fineScOnly_initMethod_4_lapSigma_2.00_seed5', 1, eyes},
                    },{
                        {'laplacianPolicy/19-06-18_500000iter_1_noLearn_fineScOnly_initMethod_4_lapSigma_10.00_seed1', 1, eyes},
                        {'laplacianPolicy/18-11-21_500000iter_2_noLearn_FineScOnly_initMethod_4_lapSigma_10.00_seed2', 1, eyes},
                        {'laplacianPolicy/18-11-22_500000iter_3_noLearn_FineScOnly_initMethod_4_lapSigma_10.00_seed3', 1, eyes},
                        {'laplacianPolicy/18-11-26_500000iter_4_noLearn_FineScOnly_initMethod_4_lapSigma_10.00_seed4', 1, eyes},
                        {'laplacianPolicy/19-06-18_500000iter_5_noLearn_fineScOnly_initMethod_4_lapSigma_10.00_seed5', 1, eyes},
                    },{
                        {'laplacianPolicy/19-06-21_500000iter_1_noLearn_fineScOnly_initMethod_4_lapSigma_20.00_seed1', 1, eyes},
                        {'laplacianPolicy/18-11-21_500000iter_2_noLearn_FineScOnly_initMethod_4_lapSigma_20.00_seed2', 1, eyes},
                        {'laplacianPolicy/18-11-22_500000iter_3_noLearn_FineScOnly_initMethod_4_lapSigma_20.00_seed3', 1, eyes},
                        {'laplacianPolicy/18-11-26_500000iter_4_noLearn_FineScOnly_initMethod_4_lapSigma_20.00_seed4', 1, eyes},
                        {'laplacianPolicy/19-06-19_500000iter_5_noLearn_fineScOnly_initMethod_4_lapSigma_20.00_seed5', 1, eyes},
                    }
            };

            names = {'0', '0.2', '1', '2', '10', '20'};
        end

        useFreq = 1
        bins = 1:length(models);
        Ns = zeros(length(models), 1); % mean relative nums
        Ns_err = zeros(length(models), 1); % std relative nums
        bf_nums = zeros(length(models), 1); % num of fits below threshold
    %     ttest_sample1 = models{length(models)};
    %     ttest_sample2 = models{length(models)};
        ttest_samples = {};

        for i = 1:length(models)
            tmp_N = 1:length(models{i}); % relative nums for each model
            for j = 1:length(tmp_N)
                if ~useFreq
                    para = load(strcat(filePath, models{i}{j}{1}, '/modelAt500000/scModel', num2str(models{i}{j}{2}), "eyes", num2str(models{i}{j}{3}), '.mat'));
                else
                    para = load(strcat(filePath, models{i}{j}{1}, '/modelAt500000/scModel', num2str(models{i}{j}{2}), "eyes", num2str(models{i}{j}{3}), '_freq.mat'));
                end
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

    %         if (i == 1)
    %             ttest_sample1 = tmp_N;
    %         end
    %         if (i == size(models, 1))
    %             ttest_sample2 = tmp_N;
    %         end
            ttest_samples{end+1} = tmp_N;

            Ns(i) = mean(tmp_N);
            Ns_err(i) = std(tmp_N);
        end

        if subPlt == 1
    %         figure('Position', [675 545 277 349]); % for different vergenc regimes
%             figure('Position', [660 569 368 323]);
            figure('Position', [660 569 1000 360]);
            subplot(1,30,1:9);
        elseif subPlt == 2
    %         figure('Position', [680 556 560 420]); % for different laplace distributions
    %         figure('Position', [1009 543 551 349]);
%             figure('Position', [1008 569 730 323]);
            subplot(1,30,12:30);
        end
    %     set(gcf, 'defaulttextinterpreter','latex')
        hold on;
        bar(bins, Ns*100, 0.7);
        errorbar(bins, Ns*100, Ns_err*100, 'k', 'Capsize', 8, 'linestyle', 'none');
        % line([1,3],[43,43], "Color", "black")
        % line([1,1],[42,44], "Color", "black")
        % line([3,3],[42,44], "Color", "black")
        % text(3.95, 44, "*", "FontSize", 12)
        hold off;

    %     grid on
        %xlabel('Normal condition')
        if subPlt == 1
%             ylabel('vertically tuned bases [%]')
            ylabel('vertically tuned RFs [%]')
        end

        if subPlt == 2
    %         xlabel('Std. of Laplacian Disparity Distribution');
            xlabel('std. of Laplacian disparity distribution');
        end
    %     xlabel('Simulation');
        % xticks([1, 2, 3, 4, 5]);%, 6, 7])
        xticks(1:length(names))
        xticklabels(names);
%         xtickangle(45)
        % workaround for linebreak in xticklabels
    %     x = '\begin{tabular}{c} line1 \\ line2 \end{tabular}';
    %     set(gca,'xtick',1:3,'XTickLabel', x, 'TickLabelInterpreter','latex');
    %     set(gca, 'TickLabelInterpreter', 'latex')
        % xticklabels({"D+L+", "D+L-", "D-L-"})%, 50, 75, 90, 100})
        % xlim([0.5, 3.5])
        xlim([0.5, length(xticks) + 0.5])
        % ylim([0 50])
        ylim([0 40])
        yticks([0,10,20,30,40]);

        % title({"Vertical Orienation - Horizontal Condition","w/ Disparities w/ Learning","w/ Disparities w/o Learning","w/o Disparities w/o Learning"}, "FontSize", 10)
        % title({"Vertical Orienation - Horizontal Condition","w/ Disparities w/ Learning","w/ Disparities w/o Learning","w/o Disparities w/o Learning"}, "FontSize", 10)
        % title({"Number of vertically tuned neurons","for different Laplacian Distributions"}, "FontSize", 10)

        % t = text(1, 5, strcat("N_{s4}=", num2str(bf_nums(1))), 'FontSize', 12,'fontWeight','bold');
        %     set(t,'Rotation',90);

        if subPlt == 1
            [h1,p1] = ttest(ttest_samples{1}, ttest_samples{2});
            [h2,p2] = ttest(ttest_samples{2}, ttest_samples{3});

            sprintf('first comparision: %d', h1)
            sprintf('second comparision: %d', h2)

            lineHeight1 = 37;
            lineHeight2 = 34;
            smallLineHeight = 0.7;

            line([1,2],[lineHeight1, lineHeight1], 'Color', 'k');
            line([1,1],[lineHeight1+smallLineHeight, lineHeight1-smallLineHeight], 'Color', 'k');
            line([2,2],[lineHeight1+smallLineHeight, lineHeight1-smallLineHeight], 'Color', 'k');

        %     text(3, 45, strcat("p=", num2str(round(p, 3))),'FontSize', 10,'fontWeight','bold')
    %         text(1.15, lineHeight1+1, strcat("p=", num2str(round(p1, 3))));
            text(1.45, lineHeight1+1, strcat("*"), 'fontsize', 16);
        %     set(gca,'FontSize',15,'fontWeight','bold') %,'FontName','Courier')
            %set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold') %,'FontName','Courier')

            line([2,3],[lineHeight2, lineHeight2], 'Color', 'k');
            line([2,2],[lineHeight2+smallLineHeight, lineHeight2-smallLineHeight], 'Color', 'k');
            line([3,3],[lineHeight2+smallLineHeight, lineHeight2-smallLineHeight], 'Color', 'k');
        %     text(3, 45, strcat("p=", num2str(round(p, 3))),'FontSize', 10,'fontWeight','bold')
            text(2.45, lineHeight2+1, strcat("*"), 'fontsize', 16);
    %         text(2.15, lineHeight2+1, strcat("p=", num2str(round(p2, 3))));
        elseif subPlt == 2
            % todo: add horizontal line for NC for the laplacian plot
            NCmean = 33.4;
            line([0,7],[NCmean, NCmean], 'LineStyle', '--', 'Color', 'k');
    %         line([0,7],[NCmean, NCmean], 'LineStyle', '--');
        end
    end
    if ~isempty(saveTag)
        saveas(gcf, sprintf('%sorientationDependence_scale%d_%s.png', savePath, scale, saveTag));
    end

    %% old models
    %%%% Plot: 0deg-MD-wDisp-wLearning %%%
            % {
            % {'/normalInput/18-03-14_500000iter_2_normal2', 2, 2},
            % {'/normalInput/18-03-14_500000iter_3_normal3', 2, 2},
            % {'/normalInput/18-03-15_500000iter_4_normal4', 2, 2},
            % {'/normalInput/18-03-10_500000iter_5_normal5', 2, 2},
            % {'/normalInput/18-03-10_500000iter_6_normal6', 2, 2}
            % },
            % {
            % {'/monocularDeprivation_diffProbs_local_s/18-03-14_500000iter_2_monDep_filtR3_prob01_s2', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-12_500000iter_3_monDep_filtR3_prob01', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-14_500000iter_4_monDep_filtR3_prob01_s4', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_5_monDep_filtR3_prob01_s5', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_6_monDep_filtR3_prob01_s6', 2, 2}
            % },
            % {
            % {'/monocularDeprivation_diffProbs_local_s/18-03-10_500000iter_2_monDep_filtR3_prob025', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-15_500000iter_3_monDep_filtR3_prob025_s3', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-14_500000iter_4_monDep_filtR3_prob025_s4', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_5_monDep_filtR3_prob025_s5', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-17_500000iter_6_monDep_filtR3_prob025_s6', 2, 2}
            % },
            % {
            % {'/monocularDeprivation_diffProbs_local_s/18-03-14_500000iter_2_monDep_filtR3_prob05_s2', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-12_500000iter_3_monDep_filtR3_prob05', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-15_500000iter_4_monDep_filtR3_prob05_s4', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_5_monDep_filtR3_prob05_s5', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-17_500000iter_6_monDep_filtR3_prob05_s6', 2, 2}
            % },
            % {
            % {'/monocularDeprivation_diffProbs_local_s/18-03-10_500000iter_2_monDep_filtR3_prob075', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_3_monDep_filtR3_prob075_s3', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-15_500000iter_4_monDep_filtR3_prob075_s4', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-17_500000iter_5_monDep_filtR3_prob075_s5', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-18_500000iter_6_monDep_filtR3_prob075_s6', 2, 2}
            % },
            % {
            % {'/monocularDeprivation_diffProbs_local_s/18-03-11_500000iter_2_monDep_filtR3_prob09', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_3_monDep_filtR3_prob09_s3', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-15_500000iter_4_monDep_filtR3_prob09_s4', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-17_500000iter_5_monDep_filtR3_prob09_s5', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-18_500000iter_6_monDep_filtR3_prob09_s6', 2, 2}
            % },
            % {
            % {'/monocularDeprivation_diffProbs_local_s/18-03-15_500000iter_2_monDep_filtR3_prob1_s2', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-13_500000iter_3_monDep_filtR3_prob1', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-16_500000iter_4_monDep_filtR3_prob1_s4', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-18_500000iter_5_monDep_filtR3_prob1_s5', 2, 2},
            % {'/monocularDeprivation_diffProbs_local_s/18-03-18_500000iter_6_monDep_filtR3_prob1_s6', 2, 2}
            % }

            %%% Plot: 0deg-MD-wDisp-w/oLearning %%%
            % {
            % {'/compareOrientedInputw-oVergence/18-05-25_500000iter_16_ActorLR=0_s16', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-05-29_500000iter_17_ActorLR=0_s17', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-05-30_500000iter_18_ActorLR=0_s18', 2, 2}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-05-28_500000iter_2_noVerg_monDep_filtR3_prob01_s2', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-01_500000iter_3_noVerg_monDep_filtR3_prob01_s3', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-05_500000iter_4_noVerg_monDep_filtR3_prob01_s4', 2, 2}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-05-28_500000iter_2_noVerg_monDep_filtR3_prob025_s2', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-02_500000iter_3_noVerg_monDep_filtR3_prob025_s3', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-06_500000iter_4_noVerg_monDep_filtR3_prob025_s4', 2, 2}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-05-27_500000iter_2_noVerg_monDep_filtR3_prob05_s2', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-05-31_500000iter_3_noVerg_monDep_filtR3_prob05_s3', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-04_500000iter_4_noVerg_monDep_filtR3_prob05_s4', 2, 2}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-05-29_500000iter_2_noVerg_monDep_filtR3_prob075_s2', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-02_500000iter_3_noVerg_monDep_filtR3_prob075_s3', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-06_500000iter_4_noVerg_monDep_filtR3_prob075_s4', 2, 2}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-05-28_500000iter_2_noVerg_monDep_filtR3_prob1_s2', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-01_500000iter_3_noVerg_monDep_filtR3_prob1_s3', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-05_500000iter_4_noVerg_monDep_filtR3_prob1_s4', 2, 2}
            % }

            %%% Plot: 0deg-MD-w/oDisp-w/oLearning %%%
            % {
            % {'/compareOrientedInputw-oVergence/18-05-25_500000iter_16_ActorLR=0_OD=3m_0disp_s16', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-05-29_500000iter_17_ActorLR=0_OD=3m_0disp_s17', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-05-30_500000iter_18_ActorLR=0_OD=3m_0disp_s18', 2, 2}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-05-26_500000iter_2_0disp_monDep_filtR3_prob01_s2', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-05-30_500000iter_3_0disp_monDep_filtR3_prob01_s3', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-03_500000iter_4_0disp_monDep_filtR3_prob01_s4', 2, 2}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-05-26_500000iter_2_0disp_monDep_filtR3_prob025_s2', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-05-30_500000iter_3_0disp_monDep_filtR3_prob025_s3', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-04_500000iter_4_0disp_monDep_filtR3_prob025_s4', 2, 2}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-05-25_500000iter_2_0disp_monDep_filtR3_prob05_s2', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-05-29_500000iter_3_0disp_monDep_filtR3_prob05_s3', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-02_500000iter_4_0disp_monDep_filtR3_prob05_s4', 2, 2}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-05-27_500000iter_2_0disp_monDep_filtR3_prob075_s2', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-05-31_500000iter_3_0disp_monDep_filtR3_prob075_s3', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-04_500000iter_4_0disp_monDep_filtR3_prob075_s4', 2, 2}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-05-26_500000iter_2_0disp_monDep_filtR3_prob1_s2', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-05-30_500000iter_3_0disp_monDep_filtR3_prob1_s3', 2, 2},
            % {'/compareOrientedInputw-oVergence/18-06-03_500000iter_4_0disp_monDep_filtR3_prob1_s4', 2, 2}
            % }

            %%% Plot: normal condition comparison of wDisp wLearning, wDisp w/oLearning
            %%% and w/oDisp w/oLearning
        %     {
        %     {'/normalInput/18-03-16_500000iter_2_normal_gf33-01-01_2', 2, 1},
        %     {'/normalInput/18-03-16_500000iter_3_normal_gf33-01-01_3', 2, 1},
        %     {'/normalInput/18-03-17_500000iter_4_normal_gf33-01-01_4', 2, 1},
        %     {'/normalInput/18-03-18_500000iter_5_normal_gf33-01-01_5', 2, 1},
        %     {'/normalInput/18-03-18_500000iter_6_normal_gf33-01-01_6', 2, 1}
        %     },
        %     {
        %     {'/compareOrientedInputw-oVergence/18-05-25_500000iter_16_ActorLR=0_s16', 2, 1},
        %     {'/compareOrientedInputw-oVergence/18-05-29_500000iter_17_ActorLR=0_s17', 2, 1},
        %     {'/compareOrientedInputw-oVergence/18-05-30_500000iter_18_ActorLR=0_s18', 2, 1}
        %     },
        %     {
        %     {'/compareOrientedInputw-oVergence/18-05-25_500000iter_16_ActorLR=0_OD=3m_0disp_s16', 2, 1},
        %     {'/compareOrientedInputw-oVergence/18-05-29_500000iter_17_ActorLR=0_OD=3m_0disp_s17', 2, 1},
        %     {'/compareOrientedInputw-oVergence/18-05-30_500000iter_18_ActorLR=0_OD=3m_0disp_s18', 2, 1}
        %     },

            %%% Plot: horizontal condition comparison of wDisp wlearning, wDisp w/oLearning
            %%% and w/oDisp w/oLearning
            % {
            % {'/compareOrientedInputw-oVergence/18-06-02_500000iter_16_horOnly_+vergence_s16', 2, 1},
            % {'/compareOrientedInputw-oVergence/18-06-03_500000iter_17_horOnly_+vergence_s17', 2, 1},
            % {'/compareOrientedInputw-oVergence/18-06-03_500000iter_18_horOnly_+vergence_s18', 2, 1}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-06-01_500000iter_16_horOnly_ActorLR=0_s16', 2, 1},
            % {'/compareOrientedInputw-oVergence/18-06-02_500000iter_17_horOnly_ActorLR=0_s17', 2, 1},
            % {'/compareOrientedInputw-oVergence/18-06-03_500000iter_18_horOnly_ActorLR=0_s18', 2, 1}
            % },
            % {
            % {'/compareOrientedInputw-oVergence/18-06-01_500000iter_16_horOnly_ActorLR=0_OD=3m_0disp_s16', 2, 1},
            % {'/compareOrientedInputw-oVergence/18-06-02_500000iter_17_horOnly_ActorLR=0_OD=3m_0disp_s17', 2, 1},
            % {'/compareOrientedInputw-oVergence/18-06-03_500000iter_18_horOnly_ActorLR=0_OD=3m_0disp_s18', 2, 1}
            % },

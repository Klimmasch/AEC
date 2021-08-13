%% This script gathers the data from all accoding simulations and generates
%% an overview of the performance in the testing procedure

%% depicted are results from normal, vertical, horizontal, orthogonal,
%% strabismic and monocular rearing

%% TODO: change ylim([0.01, 2.5]) for testing performance!
function eLifeTestHealthyModels(saveName)
    filePath = '/home/aecgroup/aecdata/Results/'; %'/SAB2018/' or 'eLifePaper'
    savePath = '/home/aecgroup/aecdata/Results/eLifePaper/plots/Aniseikonia/'


    
    %% Modes for aniseikonia (increasing image size)
    models = {
        { % normal case
            {'eLifePaper/aniseikonia/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1/modelAt500000'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2/modelAt500000'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3/modelAt500000'};
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4/modelAt500000'};
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5/modelAt500000'};
%         },{ % ani 0.03
%             {'eLifePaper/aniseikonia/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1/modelAt500000_ani0.03'};
        },{ % ani 0.05
            {'eLifePaper/aniseikonia/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1/modelAt500000_ani0.05'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2/modelAt500000_ani0.05'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3/modelAt500000_ani0.05'};
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4/modelAt500000_ani0.05'};
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5/modelAt500000_ani0.05'};
        },{ % ani 0.1
            {'eLifePaper/aniseikonia/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1/modelAt500000_ani0.10'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2/modelAt500000_ani0.10'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3/modelAt500000_ani0.10'};
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4/modelAt500000_ani0.10'};
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5/modelAt500000_ani0.10'};
        },{ % ani 0.15
            {'eLifePaper/aniseikonia/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1/modelAt500000_ani0.15'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2/modelAt500000_ani0.15'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3/modelAt500000_ani0.15'};
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4/modelAt500000_ani0.15'};
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5/modelAt500000_ani0.15'};
        },{ % ani 0.20
            {'eLifePaper/aniseikonia/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1/modelAt500000_ani0.20'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2/modelAt500000_ani0.20'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3/modelAt500000_ani0.20'};
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4/modelAt500000_ani0.20'};
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5/modelAt500000_ani0.20'};
        },{ % ani 0.25
            {'eLifePaper/aniseikonia/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1/modelAt500000_ani0.25'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_2_fsize6std_filtBoth_29_prob_1_seed2/modelAt500000_ani0.25'};
            {'eLifePaper/explFilterSizes/18-09-28_500000iter_3_fsize6std_filtBoth_29_prob_1_seed3/modelAt500000_ani0.25'};
            {'eLifePaper/explFilterSizes/18-10-02_500000iter_4_fsize6std_filtBoth_29_prob_1_seed4/modelAt500000_ani0.25'};
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_5_fsize6std_filtBoth_29_prob_1_seed5/modelAt500000_ani0.25'};
%         },{ % ani 0.5
%             {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+0.50]_scLR_[0.2,+0.2]_seed1'};
%             {'eLifePaper/aniseikonia/20-11-10_500000iter_4_aniseikonia+_aniLR_[0-+0.50]_scLR_[0.2,+0.2]_seed4'};
%         },{ % ani 1
%             {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+1.00]_scLR_[0.2,+0.2]_seed1'},
%             {'eLifePaper/aniseikonia/'}
        }}
    names = {'normal', 'ani 0.05', 'ani 0.1', 'ani 0.15', 'ani 0.20', 'ani 0.25'}%, 'ani 0.5', 'ani 1'};
%     colors = [[35,132,67]; [120,198,121]; [194,230,153]; [255,255,204]]./256; % single hue
    colors = [[35,139,69]; [116,196,118]; [186,228,179]; [237,248,233]]./256; % single hue, more saturated?
    
%     colors = [[0,104,55]; [49,163,84]; [120,198,121]; [173,221,142]; [217,240,163]; [255,255,204]]./256; % multi hue
%     colors = [[0,109,44]; [49,163,84]; [116,196,118]; [161,217,155]; [199,233,192]; [237,248,233]; [242,252,240]]./256; % single hue 6
    colors = [[0,90,50]; [35,139,69]; [65,171,93]; [116,196,118]; [161,217,155]; [199,233,192]; [237,248,233]]./256; % single hue 7
    
%     colors = [[0,109,44]; [49,163,84]; [116,196,118]; [161,217,155]; [199,233,192]; [237,248,233]]./256; % single hue green this one


    
    trainTime = 500000;
    takeDataAtIteration = 10
%     rawDat = zeros(length(models), trainTime);
%     means = zeros(length(models), trainTime);
%     stds = zeros(length(models), trainTime);
    testRes = zeros(1960, length(models{1}));
    movingMedians = zeros(length(models{1}), trainTime/10);
%     movingIQRs = zeros(length(models{1}), trainTime/10);
    interval = 10;
    ind = interval : interval : trainTime;
    
    trainVals = {};
    testVals = {};
    means = zeros(length(models), length(models{1}));
    
    for i = 1:length(models)
        trainRes = zeros(length(models{i}), trainTime/interval);
        for m = 1:length(models{i})
            
            model = load([filePath, models{i}{m}{1}, '/model.mat'], 'model');
            model = model.model;
            trainRes(m, :) = abs(model.vergerr_hist(ind)');
%             if m == 1   % only extract testing data from one experiment
%                 try
                    testRes(:, m) = model.testResult3(:,takeDataAtIteration); 
%                 end
%             end
            means(i, m) = mean(abs(model.testResult3(:,takeDataAtIteration)));
            
%             for j = windowSize+1 : model.trainTime/model.interval
% %                 movingMedians(m,j) = median(abs(model.vergerr_hist(j - windowSize : j)));
%                 movingMedians(m, j) = median(abs(trainRes(m, j - windowSize : j)));
% %                 movingIQRs(m, j) = iqr(abs(trainRes(m, j - windowSize : j)));
%             end
        end
        
        trainVals{i} = trainRes;
        testVals{i} = testRes;
    end
    
    %% Test performance plot
    figure('Position', [680 561 500 400]);
    hold on
%     boxHandl = boxplot(abs(testRes));
%     boxHandl = boxplot(abs(testRes), 'Colors', colors);
%     set(boxHandl,{'linew'},{1});
%     ax = gca;
%     ax.XTickLabels = names;
%     xlabel('Input Regime');
    
    % setting the colors
%     a = get(get(ax,'children'),'children'); % creates an array of line elements
%     t = get(a, 'tag'); % creates an array of boxplot elements, in here you can find the boxes
%     for i = 1:length(models)
%         set(a(2*length(models)+i), 'Color', colors{i});
%     end
    
    
    x = [0,0.05,0.1,0.15,0.2,0.25];
    
%     % for single model
%     y = mean(abs(testRes), 1);
%     sd=std(abs(testRes), 0, 1);
    % for multiple models
    y = mean(means, 2);
    sd = std(means, 0, 2);
    
    % transform into arcsec
%     scaling_factor = 1/(3600/60); % arcmin, replace 60 by 1 for arcsec
%     scaling_factor = (28*y(1))/atand(1/257.34); % corrected arcsec
    scaling_factor = 1;
    y = y ./ scaling_factor;
    sd = sd ./ scaling_factor;
    
    yyaxis left
    for i = 1:length(x)
        e = errorbar(x(i),y(i),sd(i),-sd(i),'o');
        e.Color = colors(i,:);
        e.MarkerSize = 8;
        e.CapSize = 10;
        e.LineWidth = 2;
    end
    
    xlim([-0.01,0.26])
    ylim([0, 0.8])
    
%     f = fit(x',y','poly2','Lower',[0,0,0])
    f = fit(x',y,'poly2')
    ax = plot(f, 'b');
    ax.LineWidth = 2;
%     ax.set('Color',colors(1,:));

    ylabel(sprintf('stereo acuity [deg]'))
%     ylabel(sprintf('stereo acuity [arcmin]'))
%     ylabel(sprintf('corrected stereo acuity [arcsec]'))
    xlabel('degree of aniseikonia [%]');
    xticklabels({0,5,10,15,20,25})
%     ylabel('')
%     set(gca, 'TickLabelInterpreter', 'latex')
    ax = gca;
    set(ax, 'YColor', 'k')

    legend off
%     grid on;
%     box off;
    % add 1-px patch as background
    hp1 = patch([-0.01, 8, 8, -0.01], [atand(1/257.34) / scaling_factor, atand(1/257.34) / scaling_factor, 0, 0], ...
                [150 / 255, 150 / 255, 150 / 255], 'LineStyle', 'none', 'LineWidth', 2, 'FaceAlpha', 0.2);

    yyaxis right
    scaling_factor = 1/(3600); % transform deg to arcsec
    max_res = atand(1/257.34)/ scaling_factor; % width of 1 px in arcsec
    
    scaling_factor = scaling_factor / (28 / max_res); % scale our max res to the max res of humans
%     scaling_factor = 1/(3600/60); % arcmin
    y = y ./ scaling_factor;
    sd = sd ./ scaling_factor;
    
    for i = 1:length(x)
        e = errorbar(x(i),y(i),sd(i),-sd(i),'o');
        e.Color = colors(i,:);
        e.MarkerSize = 8;
        e.CapSize = 10;
        e.LineWidth = 2;
    end
    
    ylim([0, 0.8] / scaling_factor)
%     ylabel("stereo acuity [arcsec]")
    ylabel(sprintf("corrected stereo\nacuity [arcsec]"))
    
    ax = gca;
    set(ax, 'YColor', 'k')
    
%     hp1.Parent = ax;
%     uistack(hp1,'bottom');
%     saveas(gca, sprintf('%sHKTestPerformance_%s.png', savePath, saveName));
    
    % remove outliers
%     outl = findobj(boxHandl, 'tag', 'Outliers');
%     set(outl, 'Visible', 'off');
%     ylim([-0.01, 2.5]);
%     ylim([0, 3]);
    ax = gca();
    ax.FontSize=17;
    
    saveas(gca, sprintf('%sHKTestPerf_poly2fit_%s.png', savePath, saveName));


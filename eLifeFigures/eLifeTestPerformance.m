%% This script gathers the data from all accoding simulations and generates
%% an overview of the performance in the testing procedure

%% depicted are results from normal, vertical, horizontal, orthogonal,
%% strabismic and monocular rearing

%% TODO: change ylim([0.01, 2.5]) for testing performance!
%% ATTENTION: Models with corrected strabismic angle need special treatment,
%% see below.
function eLifeTestPerformance(saveName)
    filePath = '/home/aecgroup/aecdata/Results/'; %'/SAB2018/' or 'eLifePaper'
    savePath = '/home/aecgroup/aecdata/Results/eLifePaper/plots/'

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
            {'eLifePaper/strabism/18-10-24_500000iter_2_filtB_29_strabAngle_10_seed2'},
            {'eLifePaper/strabism/18-10-24_500000iter_3_filtB_29_strabAngle_10_seed3'},
            {'eLifePaper/strabism/18-10-23_500000iter_4_filtB_29_strabAngle_10_seed4'},
            {'eLifePaper/strabism/19-06-26_500000iter_1_learingActivefiltB_29_strabAngle_10_seed1'},
            {'eLifePaper/strabism/19-06-26_500000iter_5_learingActivefiltB_29_strabAngle_10_seed5'},
%             {'eLifePaper/strabism/'},
%             % 10 degree strabism angle + fix obj.Dist and vergAngle
% %             {'eLifePaper/strabism/19-02-18_500000iter_2_fixAllAt6m_filtB_29_strabAngle_10_seed2'},
% %             {'eLifePaper/strabism/19-02-18_500000iter_3_fixAllAt6m_filtB_29_strabAngle_10_seed3'},
% %             {'eLifePaper/strabism/19-02-18_500000iter_4_fixAllAt6m_filtB_29_strabAngle_10_seed4'},
        }
    };
    names = {'normal', 'vertical', 'horizontal', 'orthogonal', 'monocular', 'strabismic'};
%     names = {'Normal', 'Vertical', 'Horizontal', 'Orthogonal', 'Monocular', 'Strabismic'};

    %% Modes for decreasing aniseikonia
    models = {
        { % normal case
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1'},
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1'}
        },{ % ani 0.1
            {'eLifePaper/aniseikonia/20-09-11_500000iter_1_aniseikonia_aniLR_[0-0.10]_scLR_[0.2,0.2]_seed1'},
            {'eLifePaper/aniseikonia/20-09-11_500000iter_1_aniseikonia_aniLR_[0-0.10]_scLR_[0.2,0.2]_seed1'}
        },{ % ani 0.25
            {'eLifePaper/aniseikonia/20-09-11_500000iter_1_aniseikonia_aniLR_[0-0.25]_scLR_[0.2,0.2]_seed1'},
            {'eLifePaper/aniseikonia/20-09-11_500000iter_1_aniseikonia_aniLR_[0-0.25]_scLR_[0.2,0.2]_seed1'}
        },{ % ani 0.5
            {'eLifePaper/aniseikonia/20-09-11_500000iter_1_aniseikonia_aniLR_[0-0.50]_scLR_[0.2,0.2]_seed1'},
            {'eLifePaper/aniseikonia/20-09-11_500000iter_1_aniseikonia_aniLR_[0-0.50]_scLR_[0.2,0.2]_seed1'}
        },{ % ani 0.75
            {'eLifePaper/aniseikonia/20-09-15_500000iter_1_aniseikonia_aniLR_[0-0.75]_scLR_[0.2,0.2]_seed1'},
            {'eLifePaper/aniseikonia/20-09-15_500000iter_1_aniseikonia_aniLR_[0-0.75]_scLR_[0.2,0.2]_seed1'}
%         },{ % ani 1
%             {'eLifePaper/aniseikonia/'}
        },{ % ani 2
            {'eLifePaper/aniseikonia/20-09-15_500000iter_1_aniseikonia_aniLR_[0-2.00]_scLR_[0.2,0.2]_seed1'},
            {'eLifePaper/aniseikonia/20-09-15_500000iter_1_aniseikonia_aniLR_[0-2.00]_scLR_[0.2,0.2]_seed1'}
        }}
    names = {'normal', 'ani 0.1', 'ani 0.25', 'ani 0.5', 'ani 0.75', 'ani 2'};

    %% Modes for aniseikonia (increasing image size)
    models = {
        { % normal case
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1'},
            {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1'}
        },{ % ani 0.1
            {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+0.10]_scLR_[0.2,+0.2]_seed1'},
            {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+0.10]_scLR_[0.2,+0.2]_seed1'}
        },{ % ani 0.15
            {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+0.15]_scLR_[0.2,+0.2]_seed1'},
            {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+0.15]_scLR_[0.2,+0.2]_seed1'}
        },{ % ani 0.25
            {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+0.25]_scLR_[0.2,+0.2]_seed1'},
            {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+0.25]_scLR_[0.2,+0.2]_seed1'}
        },{ % ani 0.5
            {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+0.50]_scLR_[0.2,+0.2]_seed1'},
            {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+0.50]_scLR_[0.2,+0.2]_seed1'}
        },{ % ani 1
            {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+1.00]_scLR_[0.2,+0.2]_seed1'},
            {'eLifePaper/aniseikonia/20-10-16_500000iter_1_aniseikonia+_aniLR_[0-+1.00]_scLR_[0.2,+0.2]_seed1'}
        }}
    names = {'normal', 'ani 0.1', 'ani 0.15', 'ani 0.25', 'ani 0.5', 'ani 1'};
%     %% Modes for recovery of strabismic rearing
%     models = {
%         { % normal case
%             {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1'},
%             {'eLifePaper/explFilterSizes/19-05-18_500000iter_1_fsize6std_filtBoth_29_prob_1_seed1'}
%         },{ % trained with strabismic BFs, recovery with Sc learning rate 0.2
%             {'eLifePaper/inducedStrabism/20-08-31_500000iter_2_prelearnedBF-3deg-fixedStrab_test'},
%             {'eLifePaper/inducedStrabism/20-08-31_500000iter_2_prelearnedBF-3deg-fixedStrab_test'}
%         },{ % trained with strabismic BFs, recovery with Sc learning rate 0.02
%             {'eLifePaper/inducedStrabism/20-09-05_500000iter_2_prelearnedBF-3deg-fixedStrab_sc-eta_002-002'},
%             {'eLifePaper/inducedStrabism/20-09-05_500000iter_2_prelearnedBF-3deg-fixedStrab_sc-eta_002-002'}
%         },{ % trained with strabismic BFs, recovery with Sc learning rate 0.002
%             {'eLifePaper/inducedStrabism/20-09-16_500000iter_2_prelearnedBF-3deg-fixedStrab_sc-eta_0002-0002'},
%             {'eLifePaper/inducedStrabism/20-09-16_500000iter_2_prelearnedBF-3deg-fixedStrab_sc-eta_0002-0002'}
%         },{ % trained with strabismic BFs, recovery with Sc learning rate 0.0002
%             {'eLifePaper/inducedStrabism/20-09-18_500000iter_2_prelearnedBF-3deg-fixedStrab_sc-eta_00002-00002'},
%             {'eLifePaper/inducedStrabism/20-09-18_500000iter_2_prelearnedBF-3deg-fixedStrab_sc-eta_00002-00002'}
%         }}
%     names = {'normal', 'SC LR 0.2', 'SC LR 0.02', 'SC LR 0.002', 'SC LR 0.0002'};

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
%     names = {'normal', 'strab. (2 deg)', 'strab. (5 deg)', 'strab. (10 deg)'};



    trainTime = 500000;
    takeDataAtIteration = 10
%     rawDat = zeros(length(models), trainTime);
%     means = zeros(length(models), trainTime);
%     stds = zeros(length(models), trainTime);
    testRes = zeros(1960, length(models));
    movingMedians = zeros(length(models{1}), trainTime/10);
%     movingIQRs = zeros(length(models{1}), trainTime/10);
    interval = 10;
    ind = interval : interval : trainTime;
%     windowSize = 1000
    windowSize = 2000
%     colors = {[27,120,55], [90,174,97], [166,219,160], [231,212,232], [194,165,207], [153,112,171], [118,42,131]}; % shades of green into margenta
%     colors = {[77,175,74], [55,126,184], [152,78,163], [255,127,0], [228,26,28], [27,120,55]}; % green, blue, margenta, orange, red, dark green
%     colors = {[1,102,94], [90,180,172], [216,179,101], [140,81,10]};        % turquise, light turquise, light brown, brown
%     colors = {[69,117,180], [145,191,219], [252,141,89], [215,48,39]};    % dark blue, light blue, orange, red
%     colors = {[27,120,55], [127,191,123], [175,141,195], [118,42,131]};   % dark green, light green, light purple, dark purple !
%     colors = {[77,175,74], [55,126,184], [152,78,163], [228,26,28]};      % green, blue, purple, red
    colors = {[77,175,74], [55,126,184], [152,78,163], [255,127,0], [228,26,28], [27,120,55], [118,42,131]}; % green, blue, margenta, orange, red, dark green    plotHandles = {};
    % now for all 6 sets of data
    colors = {[77,175,74]./256, [55,126,184]./256, [152,78,163]./256, [255,127,0]./256, [228,26,28]./256, [166,86,40]./256, [118,42,131]./256};
    colors = [[77,175,74]./256; [55,126,184]./256; [152,78,163]./256; [255,127,0]./256; [228,26,28]./256; [166,86,40]./256; [118,42,131]./256];
    colors = [[77,175,74]; [55,126,184]; [152,78,163]; [255,127,0]; [228,26,28]; [166,86,40]; [118,42,131]]./256;

    trainVals = {};
    testVals = {};

    figure;
    hold on;
    for i = 1:length(models)
        trainRes = zeros(length(models{i}), trainTime/interval);
        for m = 1:length(models{i})
%             if ~(i==6) %% needs to be filled in for corrected strabismic models!
                model = load([filePath, models{i}{m}{1}, '/modelAt500000/model.mat'], 'model');
%             else % load model with corrected strabismic angle
%                 model = load([filePath, models{i}{m}{1}, '/modelAt500000_strabAngle0/model.mat'], 'model');
%             end
            model = model.model;
            trainRes(m, :) = abs(model.vergerr_hist(ind)');
            if m == 1   % only extract testing data from one experiment
%                 try
                    testRes(:, i) = model.testResult3(:,takeDataAtIteration);
%                 end
            end

            for j = windowSize+1 : model.trainTime/model.interval
%                 movingMedians(m,j) = median(abs(model.vergerr_hist(j - windowSize : j)));
                movingMedians(m, j) = median(abs(trainRes(m, j - windowSize : j)));
%                 movingIQRs(m, j) = iqr(abs(trainRes(m, j - windowSize : j)));
            end
        end

        trainVals{i} = trainRes;
        testVals{i} = testRes;

        %% plotting
        xDat = [(windowSize + 1) * 10 : 10 : length(movingMedians)*10]';

        means = mean(movingMedians, 1)';
        stds = std(movingMedians, 1)';

%         [hl, hp] = boundedline(xDat, movingMedians(i, windowSize+1:end), movingIQRs(i, windowSize+1:end), 'alpha');

        %% for more than two models
%         [hl, hp] = boundedline(xDat, means(windowSize+1:end), stds(windowSize+1:end), 'alpha');
%         hl.Color = colors(i, :);%./255;
%         hl.LineWidth = 1;
%         hp.FaceColor = hl.Color;
%         plotHandles{i} = hl;

        %% simpler version for less models
        hl = line(xDat,means(windowSize+1:end), 'Color', colors(i,:));
        plotHandles{i} = hl;

%         % alternative to circumvent segmentation fault in upper function
%         lo = means(windowSize+1:end) - stds(windowSize+1:end);
%         hi = means(windowSize+1:end) + stds(windowSize+1:end);
%         hp = patch([xDat; xDat(end:-1:1); xDat(1)], [lo; hi(end:-1:1); lo(1)], colors(i,:));
%         hold on;
%
% %         hl = plot(xDat, means(windowSize+1:end), 'Color', colors(i,:));
%         hl = line(xDat,means(windowSize+1:end), 'Color', colors(i,:));
%
%         set(hp, 'edgecolor', 'none', 'facealpha', 0.3);
% %         set(hl, 'color', 'r', 'marker', 'x');
%
%         plotHandles{i} = hl;

    end


    %% Train performance plot

%     grid on;
%     xlabel('Time');
    xlabel('training iteration');
    ylabel(sprintf('abs. vergence error [deg]'));

    %TODO: needs adjustment sometimes
%     legend({'1', names{1}, '2', names{2}, '3', names{3}, '4', names{4}, '5', names{5}, '6', names{6}}, 'Location', 'east');
    legend({'1', names{1}, '2', names{2}, '3', names{3}, '4', names{4}}, 'Location', 'east');
    [~,hObj] = legend([plotHandles{:}]);
    hL=findobj(hObj, 'type', 'line');
    set(hL, 'linewidth', 3);
    ylim([0.2, 1.3])


    ylim([0.2, 2.2])

% %     title('Only Vertical Input');
    if ~isempty(saveName)
        saveas(gca, sprintf('%sTrainPerf_abs_medians_%s.png', savePath, saveName));
    end

    %% Test performance plot
    figure('Position', [680 561 672 413]);
%     boxHandl = boxplot(abs(testRes));
    boxHandl = boxplot(abs(testRes), 'Colors', colors);
    set(boxHandl,{'linew'},{1});
    ax = gca;
    ax.XTickLabels = names;
%     xlabel('Input Regime');

    % setting the colors
%     a = get(get(ax,'children'),'children'); % creates an array of line elements
%     t = get(a, 'tag'); % creates an array of boxplot elements, in here you can find the boxes
%     for i = 1:length(models)
%         set(a(2*length(models)+i), 'Color', colors{i});
%     end


    ylabel(sprintf('abs. vergence error [deg]\nat testing'))
%     ylabel('')
%     set(gca, 'TickLabelInterpreter', 'latex')


%     grid on;
    box off;
    % add 1-px patch as background
    hp1 = patch([0, 8, 8, 0], [0.22, 0.22, 0, 0], ...
                [150 / 255, 150 / 255, 150 / 255], 'LineStyle', 'none', 'LineWidth', 2, 'FaceAlpha', 0.2);

    hp1.Parent = ax;
    uistack(hp1,'bottom');
    saveas(gca, sprintf('%sTestPerformanceAt%d_%s.png', savePath, takeDataAtIteration, saveName));

    % remove outliers
    outl = findobj(boxHandl, 'tag', 'Outliers');
    set(outl, 'Visible', 'off');
    ylim([-0.01, 2.5]);
    ylim([0, 3]);
    saveas(gca, sprintf('%sTestPerformanceAt%d_%s_noOutl.png', savePath, takeDataAtIteration, saveName));

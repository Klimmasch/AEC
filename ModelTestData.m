classdef ModelTestData < handle
    properties
        %model data
        trainTime;      %number of all training iterations
        interval;       %fixation interval
        recerr_hist;    %history of rec error
        disp_hist;      %history of disparity
        vergerr_hist;   %history of vergence error
        verge_actual;   %actual vergence angle
        verge_desired;  %desired vergence angle (output of RL)
        Z;              %object depth
        fixZ;           %depth of fixation
        % td_hist;        %history of td error
        % feature_hist;   %history of feature vector
        relCmd_hist;    %relative changes in muscle commands
        cmd_hist;       %history of absulute muscle commands
        metCost_hist;   %history of metabolic costs
    end

    methods
        function obj = ModelTestData(nIterations, interval)
            obj.trainTime = nIterations;
            obj.interval = interval;
            obj.recerr_hist = zeros(nIterations, 2);    %coarse and fine scale error
            obj.disp_hist = zeros(nIterations, 1);      %saved disparity values
            obj.vergerr_hist = zeros(nIterations, 1);   %saved vergence values
            obj.verge_actual = zeros(nIterations, 1);   %saved vergence values
            obj.verge_desired = zeros(nIterations, 1);  %saved vergence values
            obj.Z = zeros(nIterations, 1);              %saved object depth
            obj.fixZ = zeros(nIterations, 1);           %saved fixation depth

            % obj.td_hist = zeros(nIterations, 1);        %history of td error
            % obj.feature_hist = zeros(nIterations, 5);
            obj.relCmd_hist = zeros(nIterations, 2);
            obj.cmd_hist = zeros(nIterations, 2);
            obj.metCost_hist = zeros(nIterations, 1);
        end

        %% Plotting and saving graphs
        function testPlotSave(this, savePath)
            windowSize = 125;
            if (this.trainTime < windowSize * this.interval)
                windowSize = round(this.trainTime / this.interval / 5);
            end
            % only take the last value before the image/texture is changed
            ind = this.interval : this.interval : this.trainTime;

            %% Simple Moving Average Vergence Error
            vergerr = filter(ones(1, windowSize) / windowSize, 1, abs(this.vergerr_hist(ind)));

            figure;
            hold on;
            grid on;
            % Raw vergance error values
            plot(this.interval : this.interval : size(this.vergerr_hist), abs(this.vergerr_hist(ind)), ...
                 'color', [1, 0.549, 0], 'LineWidth', 0.8);

            % Simple Moving Average
            plot((windowSize + 1) * this.interval : this.interval : size(this.vergerr_hist), vergerr(windowSize + 1 : end), ...
                 'color', 'b', 'LineWidth', 1.3);

            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Vergence Error [deg]', 'FontSize', 12);
            title('Moving Average of the Vergence Error');
            legend('verg err', 'SMA(verg err)');
            plotpath = sprintf('%s/mvngAvgVergErr', savePath);
            saveas(gcf, plotpath, 'png');

            %% Root Mean Squared Error
            vergerr = this.vergerr_hist(ind);
            rmse = zeros(length(1 : windowSize : length(vergerr) - mod(length(vergerr), 2)), 1); %cut if odd length
            k = 1 : windowSize : length(vergerr) - mod(length(vergerr), 2);
            for i = 1 : length(rmse)
                rmse(i) = mean(vergerr(k(i) : k(i) + windowSize - 1) .^ 2);
            end
            rmse = sqrt(rmse);

            figure;
            hold on;
            grid on;
            plot(windowSize * this.interval : windowSize * this.interval : (length(vergerr) - mod(length(vergerr), 2)) * this.interval, ...
                 rmse, 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (windowSize=%d)', windowSize * this.interval), 'FontSize', 12);
            ylabel('RMSE Vergence Error [deg]', 'FontSize', 12);
            title('RMSE of the Vergence Error');
            plotpath = sprintf('%s/rmseVergErr', savePath);
            saveas(gcf, plotpath, 'png');

            %% Reconstruction Error
            recerr_L = filter(ones(1, windowSize) / windowSize, 1, this.recerr_hist(ind, 1));
            recerr_S = filter(ones(1, windowSize) / windowSize, 1, this.recerr_hist(ind, 2));

            figure;
            hold on;
            grid on;
            plot((windowSize + 1) * this.interval : this.interval : size(this.recerr_hist), recerr_L(windowSize + 1 : end), 'r', 'LineWidth', 1.3);
            plot((windowSize + 1) * this.interval : this.interval : size(this.recerr_hist), recerr_S(windowSize + 1 : end), 'b', 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Reconstruction Error [AU]', 'FontSize', 12);
            legend('Coarse', 'Fine');
            plotpath = sprintf('%s/recErr', savePath);
            saveas(gcf, plotpath, 'png');

            %% Vergence angle
            figure;
            hold on;
            grid on;
            plot(this.verge_desired, 'color', [0, 0.7255, 0.1765], 'LineWidth', 1.3);
            plot(this.verge_actual, 'b', 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Angle [deg]', 'FontSize', 12);
            legend('desired', 'actual');
            title('Vergence');
            plotpath = sprintf('%s/vergenceAngle', savePath);
            saveas(gcf, plotpath, 'png');

            %% Muscel graphs
            figure;
            hold on;
            grid on;
            subplot(3, 1, 1);
            plot(this.cmd_hist(:, 2), 'color', [rand, rand, rand], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            title('Absolute Muscel Commands');

            subplot(3, 1, 2);
            plot(this.relCmd_hist, 'color', [rand, rand, rand], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            title('Relative Muscel Commands');

            subplot(3, 1, 3);
            plot(this.metCost_hist, 'color', [rand, rand, rand], 'LineWidth', 1.3);
            xlabel(sprintf('Iteration # (interval=%d)', this.interval), 'FontSize', 12);
            ylabel('Value', 'FontSize', 12);
            title('Metabolic Costs');

            plotpath = sprintf('%s/muscelGraphs', savePath);
            saveas(gcf, plotpath, 'png');
        end
    end
end

% This script takes some files and generated the plots and makes some testing on them
function completePlotting(reinitRenderer)
    files = {...%'/home/klimmasch/projects/results/model_12-Apr-2016_17:41:14_100000_nonhomeo_1_Bestdiscrete_highRes_cluster/model.mat', ...
%         '/home/klimmasch/projects/results/model_20-May-2016_13:10:14_500000_nonhomeo_1_2m_newImplem_highResSmSc_noMF/model.mat', ...
%         '/home/klimmasch/projects/results/model_21-May-2016_19:37:02_500000_nonhomeo_1_2m_newImplem_newSmallScale_1pcMF/model.mat', ...
%         '/home/klimmasch/projects/results/model_21-May-2016_19:37:45_500000_nonhomeo_1_2m_newImplem_newSmallScale_5pcMF/model.mat', ...
%         '/home/klimmasch/projects/results/model_21-May-2016_19:38:32_500000_nonhomeo_1_2m_newImplem_newSmallScale_10pcMF/model.mat',...
%         '/home/klimmasch/projects/results/model_20-May-2016_13:04:15_500000_nonhomeo_1_2m_newImplem_highRes-halfSize-SmallScale/model.mat',...
%         '/home/klimmasch/projects/results/model_21-May-2016_19:39:29_500000_nonhomeo_1_2m_newImplem_newSmallScale_50pcMF/model.mat'};
        '/home/klimmasch/projects/results/model_23-May-2016_21:11:52_500000_nonhomeo_1_2m_newImplem_oldScalesHalfStride_0pcMF/model.mat', ...
        '/home/klimmasch/projects/results/model_23-May-2016_11:21:31_500000_nonhomeo_1_2m_newImplem_oldScalesHalfStride_1pcMF/model.mat', ...
        '/home/klimmasch/projects/results/model_23-May-2016_22:03:15_500000_nonhomeo_1_2m_newImplem_oldScalesHalfStride_5pcMF/model.mat', ...
        '/home/klimmasch/projects/results/model_23-May-2016_11:20:48_500000_nonhomeo_1_2m_newImplem_oldScaleHalfStride_10pcMF/model.mat',...
        '/home/klimmasch/projects/results/model_23-May-2016_11:25:30_500000_nonhomeo_1_2m_newImplem_oldScalesHalfStride_30pcMF/model.mat',...
        '/home/klimmasch/projects/results/model_23-May-2016_11:19:32_500000_nonhomeo_1_2m_newImplem_oldScales_50pcMF/model.mat'};
    %     '/home/lelais/Documents/MATLAB/results/model_20-Apr-2016_20:55:52_500000_nonhomeo_1_CACLAVarLu_stdR_var_10-5_lMFB_0.3_LR_1/model.mat', ...
    %     '/home/lelais/Documents/MATLAB/results/model_20-Apr-2016_21:00:24_500000_nonhomeo_1_CACLAVarLu_stdR_var_10-4_lMFB_0.3_LR_1/model.mat', ...
    %     '/home/lelais/Documents/MATLAB/results/model_19-Apr-2016_18:06:46_500000_nonhomeo_1_CACLAVarLu_stdR_var_10-5_lMFB_0.3_LR_1_winit_1_0.01_2/model.mat', ...
    %     '/home/lelais/Documents/MATLAB/results/model_19-Apr-2016_18:09:26_500000_nonhomeo_1_CACLAVarBp_stdR_var_10-5_lMFB_0.3_LR_1_winit_1_0.01_2/model.mat', ...
    %     '/home/lelais/Documents/MATLAB/results/model_20-Apr-2016_19:19:19_500000_nonhomeo_1_CACLAVarLu_stdR_var_10-5_lMFB_0.3_LR_1_winit_1_0.01_2_posActorInit/model.mat'};
        %,...
        %'/model.mat', ...
    %     };

    simulator = OpenEyeSim('create');
    if reinitRenderer
        simulator.reinitRenderer();
    else
        simulator.initRenderer();
    end
    
    stimulusRange = 33;
    vergRange = 7;
    objDists = [0.5:0.5:2];
    objRange = length(objDists);
    zeroVergInd = ceil(vergRange/2);

    meanOffset = zeros(objRange, length(files), stimulusRange);
    for i = 1:length(files)
        model = load(files{i});
        model = model.model;
%         model.allPlotSave;
%         testModel(model, 23, [0.5, 1, 1.5, 2], [-2 : 0.2 : 2], [20, 20], 1, 0, 1, 1);
        testModel2(model, 50, 0, 1, simulator, 0); % use more stimuli than in the textures file, just in case

        sprintf('######plotting and testing completed in %s #######', files{i});
%         close all;
%         nSamples = size(model.testResult3, 1);
%         if nSamples == 0
%             simulator = OpenEyeSim('create');
%             simulator.initRenderer();
%             simulator.reinitRenderer();
%             testModel2(model, 50, 0, 1, simulator, 0);
%         end
        nSamples = size(model.testResult3, 1)
        if nSamples ~= 924
            testModel2(model, 50, 0, 1, simulator, 0);
        end
        % extract all values that are generated with zero vergErr init
        tmpResult = []; % these are the values of the 33 input stimuli that started with 0 vergErr at the end of 10 iterations
        for j = 1:nSamples


            if mod(j, vergRange) == zeroVergInd % if we started at zero vergence Error
                tmpResult = [tmpResult, model.testResult3(j,end)]; %append the last value
            end

            if length(tmpResult) == stimulusRange
                objInd = ceil((j * objRange) / nSamples);
                meanOffset(objInd, i, :) = tmpResult;
                tmpResult = [];
            end
        end
    end
    %% todo: remove outliers from the plot and scale vergence error axis to the same range
    for k = 1:objRange
        figure; hold on; grid on; grid minor;
        [~, tmpSize1, tmpSize2] = size(meanOffset(k, :, :));
        boxplot(reshape(meanOffset(k, :, :), [tmpSize1, tmpSize2])', 'labels', [0 1 5 10 30 50]);
        title(sprintf('Offset from zero vergence error at %.1f m', objDists(k)));
        xlabel('lambda_{MF} in %');
        ylabel('Vergence Error in deg');
        saveas(gcf, sprintf('../results/biasesAt%.2f_OldScale.png', objDists(k)), 'png');

    end
end
% This script takes some files and generated the plots and makes some testing on them

files = {...%'/home/klimmasch/projects/results/model_12-Apr-2016_17:41:14_100000_nonhomeo_1_Bestdiscrete_highRes_cluster/model.mat', ...
    '/home/klimmasch/projects/results/model_19-Apr-2016_18:20:30_500000_nonhomeo_1_CACLAvarBp_vergRew_winit00017-001-004_var-5_alphas1-1/model.mat', ...
    '/home/klimmasch/projects/results/model_19-Apr-2016_18:20:37_500000_nonhomeo_1_CACLAvarLu_vergRew_winit00017-001-004_var-5_alphas1-1/model.mat', ...
    '/home/klimmasch/projects/results/model_19-Apr-2016_18:20:08_500000_nonhomeo_1_CACLAvarAl_vergRew_winit00017-001-004_var-5_alphas1-1/model.mat'};
    %,...
    %'/model.mat', ...
    %'/model.mat'};

for i = 1:length(files);
    model = load(files{i});
    model = model.model;
    model.allPlotSave;
    % testModel(model, randomizationSeed, objRange, vergRange, repeat, randStimuli, randObjRange, plotIt, saveTestResults)
    testModel(model, 23, [0.5, 1, 1.5, 2], [-2 : 0.2 : 2], [25, 25], 0, 0, 1, 1);
    
    sprintf('######plotting and testing completed in %s #######', files{i});
    close all;
end
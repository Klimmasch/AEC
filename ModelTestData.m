classdef ModelTestData < handle
    properties
        %model data
        trainTime;      %number of all training iterations
        recerr_hist;    %history of rec error
        disp_hist;      %history of disparity
        vergerr_hist;   %history of vergence error
        verge_actual;   %actual vergence angle
        verge_desired;  %desired vergence angle (output of RL)
        Z;              %object depth
        fixZ;           %depth of fixation
%         g_hist;         %history of nat gradient change
%         td_hist;        %history of td error
%         feature_hist;   %history of feature vector
        relCmd_hist;    %relative changes in muscle commands
        cmd_hist;       %history of absulute muscle commands
        metCost_hist;   %history of metabolic costs
    end

    methods
        function obj = ModelTestData(nIterations)
            obj.trainTime = nIterations;
            obj.recerr_hist = zeros(nIterations, 2);    %coarse and fine scale error
            obj.disp_hist = zeros(nIterations, 1);      %saved disparity values
            obj.vergerr_hist = zeros(nIterations, 1);   %saved vergence values
            obj.verge_actual = zeros(nIterations, 1);   %saved vergence values
            obj.verge_desired = zeros(nIterations, 1);  %saved vergence values
            obj.Z = zeros(nIterations, 1);              %saved object depth
            obj.fixZ = zeros(nIterations, 1);           %saved fixation depth

%             obj.g_hist = zeros(nIterations, 1);         %history of nat gradient change
%             obj.td_hist = zeros(nIterations, 1);        %history of td error
%             obj.feature_hist = zeros(nIterations, 5);
            obj.relCmd_hist = zeros(nIterations, 2);
            obj.cmd_hist = zeros(nIterations, 2);
            obj.metCost_hist = zeros(nIterations,1);
        end
    end
end

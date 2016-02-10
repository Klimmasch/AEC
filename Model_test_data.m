classdef Model_test_data < handle
    properties        
        %model data
        recerr_hist;    %history of rec error
        disp_hist;      %history of disparity
        vergerr_hist;   %history of vergence error
        verge_actual;   %actual vergence angle
        verge_desired;  %desired vergence angle (output of RL)
        Z;              %object depth
        fixZ;           %depth of fixation 
        g_hist;         %history of nat gradient change
        td_hist;        %history of td error
        feature_hist;   %history of feature vector
    end
    
    methods
        function obj = Model_test_data()
                                   
            obj.recerr_hist = [];      %coarse and fine scale error
            obj.disp_hist = [];        %saved disparity values
            obj.vergerr_hist = [];     %saved vergence values
            obj.verge_actual = [];     %saved vergence values
            obj.verge_desired = [];    %saved vergence values
            obj.Z =[];                 %saved object depth
            obj.fixZ = [];             %saved fixation depth
            
            obj.g_hist = [];         %history of nat gradient change
            obj.td_hist = [];        %history of td error
            obj.feature_hist = [];
            
        end
    end
end
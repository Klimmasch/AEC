function model = config(learned_file,texture_file,train_time,sparse_coding)

Basis_S = [];
Basis_L = [];
weights = [];
weights_hist = cell(2,1);

%load parameters of existing model
configpath = sprintf('config/%s',learned_file);
if(~isempty(learned_file))
    load(configpath);
    Basis_S = model.scmodel_Small.Basis;
    Basis_L = model.scmodel_Large.Basis;
    
    weights = model.rlmodel.Weights; %actor-critic NN weights
    weights{3} = model.rlmodel.J;   
    weights{4} = model.rlmodel.g;   %nat gradient
        
    %weight history saved
    if ~isempty(model.rlmodel.Weights_hist)
        weights_hist = model.rlmodel.Weights_hist;
    else
        weights_hist = model.rlmodel.Weights;
    end            
end

loadBasis = 0;
loadweights = 0;
%setup parameters for sparse coding - FINE (small) SCALE 
Basis_num_used = 10;%number of basis used to encode in sparse mode
Basis_size = 128;%size of each (binocular) base vector (200)
Basis_num = 288;%total basis number (128)
eta = 0.5;%learning rate (0.01) 
Temperature = 0.01;%temperature in softmax
Dsratio = 2;        %downsampling ratio (target resolution = 8x8)
PARAMSC_S = {Basis_num_used,Basis_size,Basis_num,eta,Temperature,Dsratio,Basis_S,loadBasis};

%setup parameters for sparse coding - COARSE (large) SCALE 
Basis_num_used = 10;%number of basis used to encode in sparse mode
Basis_size = 128;%size of each (binocular) base vector (200)
Basis_num = 288;%total basis number (128)
eta = 0.5;%learning rate (0.01) 
Temperature = 0.01;%temperature in softmax
Dsratio = 8;        %downsampling ratio (target resolution = 8x8)
PARAMSC_L = {Basis_num_used,Basis_size,Basis_num,eta,Temperature,Dsratio,Basis_L,loadBasis};

PARAMSC = {PARAMSC_L,PARAMSC_S};

%setup parameters for reinforcement learning
Action = [-8 -4 -2 -1 -0.5 0 0.5 1 2 4 8];%vergence angles

alpha_v = 0.4;%learning rate to update the value function (0.05)
alpha_n = 0.4;%learning rate of natural policy gradient (0.05)
alpha_p = 0.4;%learning rate to update the policy function (1)
xi = 0.3;       %discount factor
gamma = 0.01;   %learning rate to update cumulative value; (0.01)
Temperature = 1;%temperature in softmax function in policy network
S0 = PARAMSC_L{3} + PARAMSC_S{3};%number of neurons in the input layer (Small + Large scale)
weight_range = [0.4,0.05];%maximum initial weight
lambda = 0.01;   %reguralization factor (0.01)
PARAMRL = {Action,alpha_v,alpha_n,alpha_p,xi,gamma,Temperature,lambda,S0,weight_range,loadweights,weights,weights_hist};

Interval = 10;%period to change a new environment for the eye (10)

PARAMModel = {learned_file, texture_file, train_time, sparse_coding, Interval};

PARAM = {PARAMModel,PARAMSC,PARAMRL};
model = Model(PARAM);   %create instance of object

function model = config(learnedFile, textureFile, trainTime, sparseCodingType)

basisSmall = [];
basisLarge = [];
weights = [];
weightsHist = cell(2, 1);

%load parameters of existing model
configpath = sprintf('config/%s', learnedFile);
if (~isempty(learnedFile))
    load(configpath);
    basisSmall = model.scmodel_Small.Basis;
    basisLarge = model.scmodel_Large.Basis;

    weights = model.rlmodel.Weights;    %actor-critic NN weights
    weights{3} = model.rlmodel.J;
    weights{4} = model.rlmodel.g;       %nat gradient

    %weight history saved
    if (~isempty(model.rlmodel.Weights_hist))
        weightsHist = model.rlmodel.Weights_hist;
    else
        weightsHist = model.rlmodel.Weights;
    end
end
loadBasis = uint8(~isempty(learnedFile));
loadweights = uint8(~isempty(learnedFile));

%setup parameters for sparse coding - FINE (small) SCALE
Basis_num_used = 10;    %number of basis used to encode in sparse mode
Basis_size = 128;       %size of each (binocular) base vector (200)
Basis_num = 288;        %total basis number (128)
eta = 0.2;              %learning rate | origin 0.01 | Lukas 0.1 | Alex P 0.5 | Chong 0.2
Temperature = 0.01;     %temperature in softmax | origin 0.01
Dsratio = 2;            %downsampling ratio (target resolution = 8x8)
PARAMSC_S = {Basis_num_used, Basis_size, Basis_num, eta, Temperature, Dsratio, basisSmall, loadBasis};

%setup parameters for sparse coding - COARSE (large) SCALE
Basis_num_used = 10;    %number of basis used to encode in sparse mode
Basis_size = 128;       %size of each (binocular) base vector (200)
Basis_num = 288;        %total basis number (128)
eta = 0.2;              %learning rate | origin 0.01 | Lukas 0.1 | Alex P 0.5
Temperature = 0.01;     %temperature in softmax | origin 0.01
Dsratio = 8;            %downsampling ratio (target resolution = 8x8)
PARAMSC_L = {Basis_num_used, Basis_size, Basis_num, eta, Temperature, Dsratio, basisLarge, loadBasis};

PARAMSC = {PARAMSC_L, PARAMSC_S};

%setup parameters for reinforcement learning
Action = [-8 -4 -2 -1 -0.5 0 0.5 1 2 4 8]; %vergence angles (discrete policy)

alpha_v = 0.9;                      %learning rate to update the value function | origin 0.05 | Chong 1 | Lukas 0.9 | Alex P 0.4
alpha_n = 0.025;                    %learning rate of natural policy gradient | origin 0.05 | Chong 0.025 | Lukas 0.1 | Alex P 0.4
alpha_p = 0.002;                    %learning rate to update the policy function | origin 1 | Chong 0.002 | Lukas 0.01 | Alex P 0.4
xi = 0.3;                           %discount factor | origin 0.3 | Alex P 0.3
gamma = 0.01;                       %learning rate to update cumulative value | origin 1
Temperature = 1e-5;                 %temperature in softmax function in policy network | origin 1
                                    %if policy is continuous, this value
                                    %serves as variance for the actor
% S0 = PARAMSC_L{3} + PARAMSC_S{3} + 2;%number of neurons in the input layer (Small + Large scale + Muscle activities)
S0 = PARAMSC_L{3} + PARAMSC_S{3};   %number of neurons in the input layer (Small + Large scale)
weight_range = [0.05, 0.1];         %maximum initial weight | origin [0.4, 0.05]
lambda = 0.01;                      %reguralization factor | origin 0.01
continuous = uint8(1);              %indicates if the policy is discrete or continuous
PARAMRL = {Action, alpha_v, alpha_n, alpha_p, xi, gamma, Temperature, lambda, S0, weight_range, loadweights, weights, weightsHist, continuous};

interval = 10;                      %period to change a new environment for the eye | origin 10
lambdaMuscleFB = 0.0001;            %factor of musscle activity feedback to RL feature vector

% Reward function parameters, i.e. their proportions to the reward function
% R elem [-2, 0]
lambdaMet = 1;                      %metabolic costs factor
lambdaRec = 23.8;                   %reconstruction error factor
lambdaV = 0.1;                      %value networks input->output weights factor | L1 norm 0.1 | L2 norm 0.04
lambdaP1 = 0.1;                     %policy networks input->hidden weights factor | L1 norm 0.1 | L2 norm 3.2
lambdaP2 = 1.2;                     %policy networks hidden->output weights factor | L1 norm 1.2 | L2 norm 430.8
PARAMModel = {learnedFile, textureFile, trainTime, sparseCodingType, interval, lambdaMuscleFB, lambdaMet, lambdaRec, lambdaV, lambdaP1, lambdaP2};

PARAM = {PARAMModel, PARAMSC, PARAMRL};
model = Model(PARAM);

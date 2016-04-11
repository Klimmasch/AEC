function model = config(learnedFile, textureFile, trainTime, sparseCodingType)

basisSmall = [];
basisLarge = [];
weights = [];
weightsHist = cell(2, 1);

% DEPRECATED
% TODO: update model reload
%%% Loading existing model
configpath = sprintf('config/%s', learnedFile);
if (~isempty(learnedFile))
    sprintf('Model reloading function is DEPRECATED and therefore currently not fully supported!');

    load(configpath);
    basisSmall = model.scmodel_Small.Basis;
    basisLarge = model.scmodel_Large.Basis;

    % weights = model.rlmodel.Weights;    %actor-critic NN weights
    % weights{3} = model.rlmodel.J;
    % weights{4} = model.rlmodel.g;       %nat gradient

    % % weight history saved
    % if (~isempty(model.rlmodel.Weights_hist))
    %     weightsHist = model.rlmodel.Weights_hist;
    % else
    %     weightsHist = model.rlmodel.Weights;
    % end
end
loadBasis = uint8(~isempty(learnedFile));
loadweights = uint8(~isempty(learnedFile));

%%% Sparce Coding parameters
% FINE (small) SCALE
Basis_num_used = 10;    %number of basis used to encode in sparse mode
Basis_size = 128;       %size of each (binocular) base vector (200)
Basis_num_fine = 288;   %total basis number (128)
eta = 0.2;              %learning rate | origin 0.01 | Lukas 0.1 | Alex P 0.5 | Chong 0.2
Temperature = 0.01;     %temperature in softmax | origin 0.01
Dsratio = 2;            %downsampling ratio (target resolution = 8x8)
PARAMSC_S = {Basis_num_used, Basis_size, Basis_num_fine, eta, Temperature, Dsratio, basisSmall, loadBasis};

% COARSE (large) SCALE
Basis_num_used = 10;    %number of basis used to encode in sparse mode
Basis_size = 128;       %size of each (binocular) base vector (200)
Basis_num_coarse = 288; %total basis number (128)
eta = 0.2;              %learning rate | origin 0.01 | Lukas 0.1 | Alex P 0.5
Temperature = 0.01;     %temperature in softmax | origin 0.01
Dsratio = 8;            %downsampling ratio (target resolution = 8x8)
PARAMSC_L = {Basis_num_used, Basis_size, Basis_num_coarse, eta, Temperature, Dsratio, basisLarge, loadBasis};

PARAMSC = {PARAMSC_L, PARAMSC_S};

%%% Reinforcement Learning parameters
% Critic and Actor implementation rlFlavour = [Critic, Actor]
% 0 = CNGAC             Critic Continuous Natural-Gradient Actor-Critc [Chong's implementation]
% 1 = CRG               Critic Continuous Regular Gradient
% 2 = CACLA             Critic Continuous Actor Critic Learning Automaton
%
% 0 = CNGAC             Actor Continuous Natural-Gradient Actor-Critc [Chong's implementation]
% 1 = CRG               Actor Continuous Regular Gradient [linear network topology]
% 2 = CACLALin          Actor Continuous Actor Critic Learning Automaton [linear network topology]
% 3 = CACLAVarLin       Actor Continuous Actor Critic Learning Automaton with (delta std) * update [linear network topology]
% 4 = CACLA             Actor Continuous Actor Critic Learning Automaton
% 5 = CACLAVar          Actor Continuous Actor Critic Learning Automaton with (delta std) * update
% 6 = CNGACFI           Actor Continuous Natural-Gradient Actor-Critc with Fisher Information matrix TODO: unsupported yet
rlFlavour = [uint8(2), uint8(5)];

continuous = uint8(1);                          %indicates if the policy is discrete or continuous
% Action = [-8, -4, -2, -1, -0.5, 0, ...        %vergence angles (discrete policy)
%           0.5, 1, 2, 4, 8];
Action = [-8, -4, -2, -1, -0.5, -0.2, -0.1, ... %vergence angles (discrete policy) enable half pixel resolution
          0, 0.1, 0.2, 0.5, 1, 2, 4, 8];
alpha_v = 0.9;                                  %learning rate to update the value function | origin 0.05 | Chong 1 | Lukas 0.9 | Alex P 0.4
alpha_n = 0.025;                                %learning rate of natural policy gradient | origin 0.05 | Chong 0.025 | Lukas 0.1 | Alex P 0.4
alpha_p = 0.1;                                  %learning rate to update the policy function | origin 1 | Chong 0.002 | Lukas 0.01 | Alex P 0.4 | linear 0.002
xi = 0.3;                                       %discount factor | origin 0.3 | Alex P 0.3
gamma = 0.3;                                    %learning rate to update cumulative value | origin 1
varianceRange = [1e-4, 1e-5];                   %variance of action output [start, end]
varDec = -(log(2) * trainTime) / log(varianceRange(2) / varianceRange(1)); %action variance decay factor
if continuous
    inputDim = PARAMSC_L{3} + PARAMSC_S{3} + 1; %number of neurons in the input layer (Small + Large scale + Muscle activities)
else
    inputDim = PARAMSC_L{3} + PARAMSC_S{3};     % only small + large scale basis function inputs in discrete case
    varianceRange = 1;
end
hiddenDim = 50;                                     %number of neurons in the hidden layer
weight_range = [1 / inputDim, ...                   %maximum initial weight [critic_ji, actor_ji, actor_kj]
                50 / (inputDim * hiddenDim), ...    %origin [0.05, 0.4, 0.4] | Lukas [0.1, 0.05, 0.05]
                1 / hiddenDim];                     %linear [1, 1, 1]
lambda = 0.01;                                      %reguralization factor | origin 0.01
deltaVar = 1;                                       %TD error variance tracking/approximating (CACLAVar)
eta = 0.001;                                        %TD error variance variance scaling factor (CACLAVar)
fiScale = 1e-5;                                     %scaling factor of Fisher Information matrix (CNGFI)
PARAMRL = {Action, alpha_v, alpha_n, alpha_p, xi, gamma, varianceRange, lambda, inputDim, weight_range, ...
           loadweights, weights, weightsHist, continuous, deltaVar, eta, fiScale, rlFlavour, varDec, hiddenDim};

%%% Model parameters
% Camera parameters
% offset = 0;               %vertical offset between left and right (0 in the Simulator!!!)
focalLength = 257.34;       %focal length [px]
baseline = 0.056;           %interocular distance

% Object distance to eyes [m]
objDistMin = 0.5;
objDistMax = 2;

muscleInitMin = 0.00807;    %minimal initial muscle innervation
muscleInitMax = 0.07186;    %maximal --"--

interval = 10;              %period for changing the stimulus for the eye | origin 10
lambdaMuscleFB = 1.0722;    %factor of muscle activity feedback to RL feature vector
                            %Proportion MF/feature:
                            % 0.5% = 0.0179 | 1% = 0.0357 | 5% = 0.1787 | 10% = 0.3574
                            % 20% = 0.7148 | 30% = 1.0722 | 40% = 1.4296 | 50% = 1.7871 | 100% = 3.5741

% Reward function parameters, i.e. their proportions to the reward function
% R elem [-2, 0]
lambdaRec = 4.929;          %reconstruction error factor | 4.929
lambdaMet = 0;              %metabolic costs factor | 0.204
lambdaV = 7.0282e-04;       %value networks input->output weights factor | L1 norm 7.0282e-04
lambdaP1 = 0.019;           %policy networks input->hidden weights factor | L1 norm 0.019
lambdaP2 = 0.309;           %policy networks hidden->output weights factor | L1 norm 0.309
PARAMModel = {learnedFile, textureFile, trainTime, sparseCodingType, focalLength, baseline, ...
              objDistMin, objDistMax, muscleInitMin, muscleInitMax, interval, ...
              lambdaMuscleFB, lambdaMet, lambdaRec, lambdaV, lambdaP1, lambdaP2};

PARAM = {PARAMModel, PARAMSC, PARAMRL};
model = Model(PARAM);

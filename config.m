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
% 0 = Chong's implementation
% 1 = CRG       Critic Continuous Regular Gradient
% 2 = CACLA     Critic Continuous Actor Critic Learning Automaton
%
% 0 = Chong's implementation
% 1 = CRG       Actor Continuous Regular Gradient
% 2 = CACLA     Actor Continuous Actor Critic Learning Automaton
% 3 = CACLAVar  Actor Continuous Actor Critic Learning Automaton with (delta std) * update
% 4 = CNGFI     Actor Continuous Natural Gradient with Fisher Information matrix TODO: unsupported yet
rlFlavour = [2, 3];

Action = [-8 -4 -2 -1 -0.5 0 0.5 1 2 4 8]; %vergence angles (discrete policy)
alpha_v = 0.9;                          %learning rate to update the value function | origin 0.05 | Chong 1 | Lukas 0.9 | Alex P 0.4
alpha_n = 0.025;                        %learning rate of natural policy gradient | origin 0.05 | Chong 0.025 | Lukas 0.1 | Alex P 0.4
alpha_p = 0.002;                        %learning rate to update the policy function | origin 1 | Chong 0.002 | Lukas 0.01 | Alex P 0.4
xi = 0.3;                               %discount factor | origin 0.3 | Alex P 0.3
gamma = 0.3;                            %learning rate to update cumulative value | origin 1
variance = 1e-5;                        %variance of action output
                                        %if policy is continuous, this value
                                        %serves as variance for the actor
S0 = PARAMSC_L{3} + PARAMSC_S{3} + 1;   %number of neurons in the input layer (Small + Large scale + Muscle activities)
weight_range = [0.0017, 0.0017, 0.05];  %maximum initial weight [critic_ji, actor_ji, actor_kj] | origin [0.05, 0.4, 0.4] | Lukas [0.1, 0.05, 0.05]
lambda = 0.01;                          %reguralization factor | origin 0.01
continuous = uint8(1);                  %indicates if the policy is discrete or continuous
deltaVar = 1;                           %delta variance tracking (CACLAVar)
eta = 0.001;                            %delta vairance scaling factor (CACLAVar)
fiScale = 1e-5;                         %scaling factor of Fisher Information matrix (CNGFI)
PARAMRL = {Action, alpha_v, alpha_n, alpha_p, xi, gamma, variance, lambda, S0, weight_range, ...
           loadweights, weights, weightsHist, continuous, deltaVar, eta, fiScale, rlFlavour};

%%% Model parameters
% Camera parameters
% offset = 0;           %vertical offset between left and right (0 in the Simulator!!!)
focalLength = 257.34;   %focal length [px]
baseline = 0.056;       %interocular distance

% Object distance to eyes [m]
objDistMin = 0.5;
objDistMax = 2;

muscleInitMin = 0.00807;    %minimal initial muscle innervation
muscleInitMax = 0.07186;    %maximal --"--

interval = 10;                          %period to change a new environment for the eye | origin 10
lambdaMuscleFB = 1.0722;                %factor of muscle activity feedback to RL feature vector
                                        %Proportion MF/feature | 0.5% = 0.0179 | 1% = 0.0357 | 5% = 0.1787 | 10% = 0.3574
                                        % 30% = 1.0722 | 50% = 1.7871 | 100% = 3.5741

% Reward function parameters, i.e. their proportions to the reward function
% R elem [-2, 0]
lambdaRec = 4.929;                      %reconstruction error factor | 4.929
lambdaMet = 0;                          %metabolic costs factor | 0.204
lambdaV = 7.0282e-04;                   %value networks input->output weights factor | L1 norm 7.0282e-04
lambdaP1 = 0.019;                       %policy networks input->hidden weights factor | L1 norm 0.019
lambdaP2 = 0.309;                       %policy networks hidden->output weights factor | L1 norm 0.309
PARAMModel = {learnedFile, textureFile, trainTime, sparseCodingType, focalLength, baseline, ...
              objDistMin, objDistMax, muscleInitMin, muscleInitMax, interval, ...
              lambdaMuscleFB, lambdaMet, lambdaRec, lambdaV, lambdaP1, lambdaP2};

PARAM = {PARAMModel, PARAMSC, PARAMRL};
model = Model(PARAM);

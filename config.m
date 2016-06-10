function model = config(textureFile, trainTime, sparseCodingType)

weights = [];
weightsHist = cell(2, 1);
loadBasis = uint8(0);
loadweights = uint8(0);

%%% Sparce Coding parameters
% Scales := [coarse, less_coarse, ..., fine], i.e. [peripheral vision, ..., central vision]
nBasis = [288, 288];            % total number of basis | origin [288, 288]
nBasisUsed = [10, 10];          % number of basis used to encode in sparse mode | origin [10, 10]
basisSize = [128, 128];         % size of each (binocular) base vector: patchSize * patchSize * 2 (left + right eye) | origin [128, 128]
eta = [0.2, 0.2];               % learning rate [origin 0.01 | Lukas 0.1 | Alex P 0.5, origin 0.01 | Lukas 0.1 | Alex P 0.5 | Chong 0.2]
temperature = [0.01, 0.01];     % temperature in softmax | origin 0.01

PARAMSC = {nBasis, nBasisUsed, basisSize, eta, temperature};

%%% Reinforcement Learning parameters
% Critic and Actor implementation rlFlavour = [Critic, Actor]
% 0 = CACLA             Critic Continuous Actor Critic Learning Automaton
% 1 = CRG               Critic Continuous Regular Gradient
%
% 0 = CACLAVarLu        Actor Continuous Actor Critic Learning Automaton with (delta std) * update [Lukas's interpretation of CACLA appoach]
% 1 = CACLAVarAl        Actor Continuous Actor Critic Learning Automaton with (delta std) * update [Alex L's interpretation of CACLA appoach]
% 2 = CACLAVarBp        Actor Continuous Actor Critic Learning Automaton with (delta std) * update [CACLA appoach with std. Backpropagation]
% 3 = CACLAVar          Actor Continuous Actor Critic Learning Automaton with (delta std) * update
% 4 = CACLAVar2         Actor Continuous Actor Critic Learning Automaton with (delta std) * update [non-linear output layer]
% # = CNGACFI           Actor Continuous Natural-Gradient Actor-Critc with Fisher Information matrix TODO: unsupported yet
rlFlavour = [uint8(0), uint8(0)];

continuous = uint8(1);                               % indicates if the policy is discrete (= 0) or continuous (= 1)
% actionSpace = [-8, -4, -2, -1, -0.5, 0, ...        % vergence angles (discrete policy)
%                 0.5, 1, 2, 4, 8];
actionSpace = [-8, -4, -2, -1, -0.5, -0.2, -0.1, ... % vergence angles (discrete policy) enable half pixel resolution
                0, 0.1, 0.2, 0.5, 1, 2, 4, 8];
alpha_v = 1;                                         % learning rate to update the value function | origin 0.05 | Chong 1 | Lukas 0.9 | Alex P 0.4
alpha_n = 0.025;                                     % learning rate of natural policy gradient | origin 0.05 | Chong 0.025 | Lukas 0.1 | Alex P 0.4
alpha_p = 1;                                         % learning rate to update the policy function | origin 1 | Chong 0.002 | Lukas 0.01 | Alex P 0.4 | linear 0.002
xi = 0.3;                                            % discount factor | origin 0.3 | Alex P 0.3
gamma = 0.3;                                         % learning rate to update cumulative value | origin 1

varianceRange = [1e-5, 1e-5];                        % variance of action output, i.e. variance of Gaussian policy [training_start, training_end]
                                                     % corresponds to softMax temperature in discrete RL models
if (size(varianceRange, 2) == 1 || varianceRange(1) <= varianceRange(2))
    varDec = 0; % no variance decay
else
    varDec = -(log(2) * trainTime) / log(varianceRange(2) / varianceRange(1)); % action variance decay factor
end

outputDim = 2;                                      % number of neurons in the output layer and amount of eye muscles
if (continuous == 1)
    inputDim = sum(PARAMSC{1}) + outputDim;         % number of neurons in the input layer (Small + Large scale + Muscle activities)
else
    inputDim = sum(PARAMSC{1});                     % only small + large scale basis function inputs in discrete case
    varianceRange = 1;
end
hiddenDim = 50;                                     % number of neurons in the hidden layer
dimensions = [inputDim, hiddenDim, outputDim];

weight_range = [1 / inputDim, ...                   % maximum initial weight [critic_ji, actor_ji, actor_kj]
                0.01, ...                           % origin [0.05, 0.4, 0.4] | Lukas [0.1, 0.05, 0.05] | AL 100 / (inputDim * hiddenDim)
                2 / hiddenDim * outputDim];         % linear [1/inputDim, 1/inputDim, -]
lambda = 0.01;                                      % reguralization factor | origin 0.01
deltaVar = 1;                                       % TD error variance tracking/approximating (CACLAVar)
eta = 0.001;                                        % TD error variance scaling factor (CACLAVar)
fiScale = 1e-5;                                     % scaling factor of Fisher Information matrix (CNGFI)

PARAMRL = {actionSpace, alpha_v, alpha_n, alpha_p, xi, gamma, varianceRange, lambda, dimensions, weight_range, ...
           loadweights, weights, weightsHist, continuous, deltaVar, eta, fiScale, rlFlavour, varDec};

%%% Model parameters
% Image processing constants
patchSize = 8;                                  % patch size [pixel] of one basis functon, i.e. "receptive field" size | origin 8
% [peripheral vision, ..., central vision]
dsRatio = [8, 1];                               % downsampling ratio, i.e. how many pixels in original image
                                                % correspond to how many pixels in downsampled image | origin [8, 2]
pxFieldOfViewOrig = [128, 40];                  % fields of view in original image [pixel] | origin [128, 80]
pxFieldOfView = pxFieldOfViewOrig ./ dsRatio;   % fields of view in downsampled image [pixel] (previously called fovea)
                                                % pxFieldOfView = FieldOfView in original image [pixel] / dsRatio
stride = [patchSize / 2, patchSize / 2];        % image patch strides | origin [1, patchSize / 2]

% Sanity check of parameter values
if (~all(diff(pxFieldOfViewOrig) < 0))
    sprintf('pxFieldOfViewOrig must contain decreasing values,\ndue to convention it must hold [peripheral vision, ..., central vision].')
    return;
% elseif (mod(pxFieldOfView(2 : end), dsRatio(1 : end - 1) ./ dsRatio(2 : end)))
%     sprintf('pxFieldOfView(scale + 1) / (dsRatio(scale) / dsRatio(scale + 1)) must be an integer')
%     return;
% elseif (mod(pxFieldOfView(1 : end - 1) - pxFieldOfView(2 : end) ./ dsRatio(1 : end - 1) ./ dsRatio(2 : end), 2))
%     sprintf('pxFieldOfView(scale) - (pxFieldOfView(scale + 1) / (dsRatio(scale) / dsRatio(scale + 1))) must be an even integer')
%     return;
end

% Camera parameters
% offset = 0;               % vertical offset between left and right (0 in the iCub Simulator!)
focalLength = 257.34;       % focal length [px]
baseline = 0.056;           % interocular distance

% Object distance to eyes [m]
objDistMin = 0.5;
objDistMax = 2;

% Fixation distance [m]
% used for eye fixation initialization
fixDistMin = 0.3379;
fixDistMax = 3.2219;

% muscle initialization: correspond now to the minimum and maximum distance
% the eyes should be looking at
muscleInitMin = 0.00807;       % minimal initial muscle innervation orig: 0.00807 corr. to vergAngleMin | 0 corr. to 1 deg
muscleInitMax = 0.07186;       % maximal --"--, orig: 0.07186 corr. to vergAngleMax | 0.1 corrs. to 12.7 deg

interval = 10;              % period for changing the stimulus for the eyes | origin 10
lambdaMuscleFB = 0.0357;    % factor of muscle activity feedback to RL feature vector
                            % Proportion MF/feature:
                            % 0.5% = 0.0179 | 1% = 0.0357 | 5% = 0.1787 | 10% = 0.3574
                            % 20% = 0.7148 | 30% = 1.0722 | 40% = 1.4296 | 50% = 1.7871 | 100% = 3.5741

% Reward function parameters, i.e. their proportions to the reward function
% R elem [-2, 0]
lambdaRec = 4.929;          % reconstruction error factor | privious 77.12% = 4.929 | 100% = 6.391
lambdaMet = 0.012;          % metabolic costs factor | privious 12.75% =  0.204 | 10% = 0.161 | 5% = 0.081 | 1% = 0.016 | 0.75% = 0.012 | 0.5% = 0.008
lambdaV = 7.0282e-04;       % value networks input->output weights factor | L1 norm 7.0282e-04
lambdaP1 = 0.019;           % policy networks input->hidden weights factor | L1 norm 0.019
lambdaP2 = 0.309;           % policy networks hidden->output weights factor | L1 norm 0.309
PARAMModel = {textureFile, trainTime, sparseCodingType, focalLength, baseline, ...
              objDistMin, objDistMax, muscleInitMin, muscleInitMax, interval, ...
              lambdaMuscleFB, lambdaMet, lambdaRec, lambdaV, lambdaP1, lambdaP2, ...
              patchSize, pxFieldOfView, dsRatio, stride, fixDistMin, fixDistMax};

PARAM = {PARAMModel, PARAMSC, PARAMRL};
model = Model(PARAM);

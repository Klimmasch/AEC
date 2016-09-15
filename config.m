function model = config(textureFile, trainTime, sparseCodingType)

%%% Model parameters
% Image processing constants
patchSize = 8;                                  % patch size [pixel] of one basis functon, i.e. "receptive field" size | origin 8

% [peripheral vision, ..., central vision]
dsRatio = [4, 1];                               % downsampling ratio, i.e. how many pixels in original image
                                                % correspond to how many pixels in downsampled image | origin [8, 2]
pxFieldOfViewOrig = [128, 40];                  % fields of view in original image [pixel] | origin [128, 80]
pxFieldOfView = pxFieldOfViewOrig ./ dsRatio;   % fields of view in downsampled image [pixel] (previously called fovea)
                                                % pxFieldOfView = FieldOfView in original image [pixel] / dsRatio
stride = [patchSize / 2, patchSize / 2];        % image patch strides | origin [1, patchSize / 2]
cutout = uint8(0);                              % Manages whether cutout procedure is applied [1] or not [0]
overlap = [0];                                  % Overlap between the different layers measured in units of FINE scale - works only in conjunction with cutout

% Sanity check of parameter values
if (~all(diff(pxFieldOfViewOrig) < 0))
    sprintf('pxFieldOfViewOrig must contain decreasing values,\ndue to convention it must hold [peripheral vision, ..., central vision].')
    return;
% TODO: needs to be inserted and checked when model.preprocessImageCutout() is integrated
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
baseline = 0.056;           % interocular distance [m]

% Object distance to eyes [m]
objDistMin = 0.5; % origin 0.5
objDistMax = 6;   % origin 2

% Fixation distance [m]
% used for eye fixation initialization
fixDistMin = 0.3379;
fixDistMax = 6; %3.2219 for objDistMax = 2m

% Muscle initialization [%]: correspond now to the minimum and maximum distance
% the eyes should be looking at. [lateral rectus, medial rectus]
% some correspondances (distance: [lateral, medial] activation):
% 0.5m : [0, 0.0726], 1.5m : [0, 0.0166], 1.5m-2deg : [0, 0.0441], 2m : [0, 0.0089], 3m : [0, 0.0011], 3.22m : [0, 0],
% 4m : [0.0027, 0], 6m : [0.0064, 0], 10m : [0.0093, 0], Inf : [0.0136, 0]
muscleInitMin = [0, 0];             % minimal initial muscle innervation orig: 0.00807 corr. to vergAngleMin | 0 corr. to 1 deg
muscleInitMax = [0.0064, 0.0166];   % maximal --"--, orig: 0.07186 corr. to vergAngleMax | 0.1 corrs. to 12.7 deg

% period for changing the stimulus for the eyes | origin 10
interval = 10;

%%% Factor of Muscle activity feedback to RL feature vector
% here the % is in respect to the feature signal,
% i.e. 100% feedback = mean(sum(total muscle force)) ~= mean(sum(feature_vetor(1 : PARAMSC{1})))
% 0.1% = 0.036 | 0.5% = 0.0179 | 1% = 0.0357 | 2.5% = 0.0893 | 5% = 0.1787 | 7.5% = 0.2678
% 10% = 0.3574 | 20% = 0.7148 | 30% = 1.0722 | 40% = 1.4296 | 50% = 1.7871 | 100% = 3.5741
%
%%%
% here the % corresponds to the mean of all single components of the feature vector
% mean(mean(model.feature_hist(:, 1 : 576), 2)) = 0.0059
% mean(sum(model.feature_hist(:, 1 : 576), 2)) = 3.3721
% mean(model.cmd_hist) = [0.0229, 0.0465], 0.0465 / 0.0059 = 1 / lambdaMuscleFB
%
% 1% = 0.001268817 | 5% = 0.0063
% 10% = 0.0127 | 20% = 0.0254 | 30% = 0.0381 | 40% = 0.0508 | 50% = 0.0634
% 60% = 0.0761 | 70% = 0.0888 | 80% = 0.1015 | 90% = 0.1142 | 100% = 0.1269
% 150% = 0.1903 | 200% = 0.2538 | 250% = 0.3172 | 300% = 0.3806
%
%%%
lambdaMuscleFB = 0.1269;

%%% Reward function parameters, i.e. their "proportions" to the whole reward signal
%%% Reconstruction error factor
% the % is in respect to the reconstruction error, i.e. 100% X means signal X is as strong as
% 6.391 * mean reconstruction error on average, whereby 6.391 * mean reconstruction error ~= 1
% pure 15.647% = 1 | privious 77.12% = 4.929 | 100% = 6.391 | set to 1 for simplicity
%%
lambdaRec = 1;

% due to recError reduction at dsRatio = [8, 2], lambdaRec needs to be scaled accordingly
% this is a temporary solution, you should recalculate everything if you introduce 3 scales!
if (dsRatio(end) > 1)
    lambdaRec = lambdaRec * 1.7706;
end

%%% Metabolic costs factor
% per_cent = mean(model.metCost_hist) * lambdaMet * lambdaRec / mean(recError) * lambdaRec
% privious 12.75% =  0.204 | 10% = 0.161 | 5% = 0.081 | 1% = 0.016 | 0.8% = 0.0128 | 0.75% = 0.012 | 0.6% = 0.0096 | 0.5% = 0.008
%%% NEW
% sum(sum(model.recerr_hist))/size(model.recerr_hist, 1) = 0.2997 @ dsRatio = [8, 1]
% sum(sum(model.recerr_hist))/size(model.recerr_hist, 1) = 0.1693 @ dsRatio = [8, 2]
% mean(model.metCost_hist) = 1.0380
%
% 0.5% = 0.0014435 | 0.25% = 0.00072175 | 0.1% = 0.0002887
% 1% = 0.0029 | 2% = 0.0058 | 3% = 0.0087 | 4% = 0.0114 | 5% = 0.0144
% 6% = 0.0173 | 7% = 0.0203 | 8% = 0.0232 | 9% = 0.0261 | 10% = 0.0289
% 15% = 0.0435 | 25% = 0.0722 | 30% = 0.0866 | 50% = 0.1443 | 100% = 0.2887
lambdaMet = 0;

% due to the dependancy of mean(model.metCost_hist) * lambdaMet * lambdaRec / mean(recError) * lambdaRec = x%
% lambdaMet needs to be scaled accordingly
lambdaMet = lambdaMet * lambdaRec;

PARAMModel = {textureFile, trainTime, sparseCodingType, focalLength, baseline, ...
              objDistMin, objDistMax, muscleInitMin, muscleInitMax, interval, ...
              lambdaMuscleFB, lambdaRec, lambdaMet, patchSize, pxFieldOfView, ...
              dsRatio, stride, fixDistMin, fixDistMax, overlap, cutout};

%%% Sparce Coding parameters
% Scales := [coarse, less_coarse, ..., fine], i.e. [peripheral vision, ..., central vision]
nBasis = [400, 400];                                        % total number of basis | origin [288, 288]
nBasisUsed = [10, 10];                                      % number of basis used to encode in sparse mode | origin [10, 10]
basisSize = [(patchSize ^ 2) * 2, (patchSize ^ 2) * 2];     % size of each (binocular) base vector: patchSize * patchSize * 2 (left + right eye) | origin [128, 128]
eta = [0.2, 0.2];                                           % learning rate [origin 0.01 | Lukas 0.1 | Alex P 0.5, origin 0.01 | Lukas 0.1 | Alex P 0.5 | Chong 0.2]
temperature = [0.01, 0.01];                                 % temperature in softmax | origin 0.01

% consistency check
if ((length(pxFieldOfViewOrig) ~= length(dsRatio)) ...
 || (length(stride) ~= length(dsRatio)) ...
 || (length(nBasis) ~= length(dsRatio)) ...
 || (length(nBasisUsed) ~= length(dsRatio)) ...
 || (length(basisSize) ~= length(dsRatio)) ...
 || (length(eta) ~= length(dsRatio)) ...
 || (length(temperature) ~= length(dsRatio)))
    sprintf('Error: For usage of #%d scales all respective SC and model parameters needs to have length %d.', ...
            length(dsRatio), length(dsRatio))
    return;
end

if (length(dsRatio) > 1 && (length(overlap) ~= length(dsRatio) - 1))
    sprintf('Error: For usage of #%d scales overlap needs to have length %d.', length(dsRatio) - 1)
    return;
end

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
alpha_p = 0.5;                                       % learning rate to update the policy function | origin 1 | Chong 0.002 | Lukas 0.01 | Alex P 0.4 | linear 0.002
xi = 0.3;                                            % discount factor | origin 0.3 | Alex P 0.3
gamma = 0.3;                                         % learning rate to update cumulative value | origin 1
regularizer = 1 - 1e-3;                              % actor weight regularization via factorial downscaling

varianceRange = [1e-5, 1e-5];                        % variance of action output, i.e. variance of Gaussian policy [training_start, training_end]
                                                     % corresponds to softMax temperature in discrete RL models
if (length(varianceRange) == 1 || varianceRange(1) == varianceRange(2))
    varDec = 0; % no variance decay
elseif (varianceRange(1) < varianceRange(2))
    sprintf('Error: It must hold varianceRange(1) >= varianceRange(2)')
    return;
else
    varDec = -(log(2) * trainTime) / log(varianceRange(2) / varianceRange(1)); % action variance decay factor
end

outputDim = 2;                                      % number of neurons in the output layer and amount of eye muscles
if (continuous == 1)
    inputDim = sum(PARAMSC{1}) + outputDim;         % number of neurons in the input layer (Small + Large scale + Muscle activities)
else
    inputDim = sum(PARAMSC{1});                     % only small + large scale basis function inputs in discrete models
    varianceRange = 1;
    outputDim = 1;                                  % only one delta angle output in discrete models
end
hiddenDim = 50;                                     % number of neurons in the hidden layer
dimensions = [inputDim, hiddenDim, outputDim];

weight_range = [1 / inputDim, ...                   % maximum initial weight [critic_ji, actor_ji, actor_kj]
                1 / (inputDim * hiddenDim), ...     % origin [0.05, 0.4, 0.4] | Lukas [0.1, 0.05, 0.05] | AL 100 / (inputDim * hiddenDim)
                2 / (hiddenDim * outputDim)];       % linear [1/inputDim, 1/inputDim, -]
lambda = 0.01;                                      % reguralization factor | origin 0.01
deltaVar = 1;                                       % TD error variance tracking/approximating (CACLAVar)
eta = 0.001;                                        % TD error variance scaling factor (CACLAVar)
fiScale = 1e-5;                                     % scaling factor of Fisher Information matrix (CNGFI)

PARAMRL = {actionSpace, alpha_v, alpha_n, alpha_p, xi, gamma, varianceRange, lambda, dimensions, weight_range, ...
           continuous, deltaVar, eta, fiScale, rlFlavour, varDec, regularizer};

PARAM = {PARAMModel, PARAMSC, PARAMRL};
model = Model(PARAM);

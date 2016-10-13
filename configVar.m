% Model object configuration and generation
% with support of variable parameter vector
function model = configVar(myvarargin)

% --------------------
% Experiment paramters
% --------------------

% stimulus file name
[found, textureFile, myvarargin] = parseparam(myvarargin, 'textureFile');
if (~found)
    textureFile = {'Textures_mcgillManMade40.mat', 'Textures_mcgillManMade100.mat'};
end

% training duration
[found, trainTime, myvarargin] = parseparam(myvarargin, 'trainTime');
if (~found)
    trainTime = 1000;
end

if (~isscalar(trainTime) || trainTime < 1)
    error('trainTime must be scalar > 0');
end

% points in time of intermediate test procedure during training
[found, testAt, myvarargin] = parseparam(myvarargin, 'testAt');
if (~found)
    testAt = [1000000 : 1000000 : trainTime];
end

% sparse coding type
[found, sparseCodingType, myvarargin] = parseparam(myvarargin, 'sparseCodingType');
if (~found)
    sparseCodingType = uint8(0);
end

if (~isscalar(sparseCodingType) || sparseCodingType < 0 || sparseCodingType > 1)
    error('sparseCodingType must be scalar in {0, 1}');
end

% ----------------
% Model parameters
% ----------------

% period for changing the stimulus for the eyes | origin 10
[found, interval, myvarargin] = parseparam(myvarargin, 'interval');
if (~found)
    interval = 10;
end

if (~isscalar(interval) || interval < 1)
    error('interval must be scalar >= 1');
end

%%% Image processing constants

% patch size [pixel] of one basis functon, i.e. "receptive field" size | origin 8
[found, patchSize, myvarargin] = parseparam(myvarargin, 'patchSize');
if (~found)
    patchSize = 8;
end

if (~isscalar(patchSize) || patchSize < 4)
    error('patchSize must be scalar >= 4');
end

% [peripheral, intermediate ..., central vision]
% downsampling ratio, i.e. how many pixels in original image
% correspond to how many pixels in downsampled image | origin [8, 2]
[found, dsRatio, myvarargin] = parseparam(myvarargin, 'dsRatio');
if (~found)
    dsRatio = [4, 1];
end

% fields of view in original image [pixel] | origin [128, 80]
[found, pxFieldOfViewOrig, myvarargin] = parseparam(myvarargin, 'pxFieldOfViewOrig');
if (~found)
    pxFieldOfViewOrig = [128, 40];
end

% fields of view in downsampled image [pixel] (previously called fovea)
% pxFieldOfView = FieldOfView in original image [pixel] / dsRatio
[found, pxFieldOfView, myvarargin] = parseparam(myvarargin, 'pxFieldOfView');
if (~found)
    pxFieldOfView = pxFieldOfViewOrig ./ dsRatio;
end

if (~all(diff(pxFieldOfViewOrig) < 0))
    error('pxFieldOfViewOrig must contain decreasing values,\ndue to convention it must hold [peripheral, intermediate ..., central vision].');
% TODO: needs to be inserted and checked when model.preprocessImageCutout() is integrated
% elseif (mod(pxFieldOfView(2 : end), dsRatio(1 : end - 1) ./ dsRatio(2 : end)))
%     sprintf('pxFieldOfView(scale + 1) / (dsRatio(scale) / dsRatio(scale + 1)) must be an integer')
%     return;
% elseif (mod(pxFieldOfView(1 : end - 1) - pxFieldOfView(2 : end) ./ dsRatio(1 : end - 1) ./ dsRatio(2 : end), 2))
%     sprintf('pxFieldOfView(scale) - (pxFieldOfView(scale + 1) / (dsRatio(scale) / dsRatio(scale + 1))) must be an even integer')
%     return;
end

% image patch strides | origin [1, patchSize / 2]
[found, stride, myvarargin] = parseparam(myvarargin, 'stride');
if (~found)
    stride = [patchSize / 2, patchSize / 2];
end

% flag indicates whether cutout procedure is applied [1] or not [0]
[found, cutout, myvarargin] = parseparam(myvarargin, 'cutout');
if (~found)
    cutout = uint8(0);
end

if (~isscalar(cutout) || cutout < 0 ||  cutout > 1)
    error('cutout must be scalar in {0, 1}');
end

% overlap between the different layers measured in units of FINE scale,
% works only in conjunction with cutout
[found, overlap, myvarargin] = parseparam(myvarargin, 'overlap');
if (~found)
    overlap = 0;
end

if (length(dsRatio) > 1 && (length(overlap) ~= length(dsRatio) - 1))
    error('For usage of #%d scales overlap needs to have length %d.', length(dsRatio) - 1);
end

%%% Camera parameters

% vertical offset between left and right (0 in the iCub Simulator!)
% [found, offset, myvarargin] = parseparam(myvarargin, 'offset');
% if (~found)
%     offset = 0;
% end

% focal length [px]
[found, focalLength, myvarargin] = parseparam(myvarargin, 'focalLength');
if (~found)
    focalLength = 257.34;
end

% interocular distance [m]
[found, baseline, myvarargin] = parseparam(myvarargin, 'baseline');
if (~found)
    baseline = 0.056;
end

% Object distance to eyes [m]
[found, objDistMin, myvarargin] = parseparam(myvarargin, 'objDistMin');
if (~found)
    objDistMin = 0.5; % origin 0.5
end

[found, objDistMax, myvarargin] = parseparam(myvarargin, 'objDistMax');
if (~found)
    objDistMax = 6;   % origin 2
end

% Fixation distance [m]
% used for eye fixation initialization
[found, fixDistMin, myvarargin] = parseparam(myvarargin, 'fixDistMin');
if (~found)
    fixDistMin = 0.3379;
end

[found, fixDistMax, myvarargin] = parseparam(myvarargin, 'fixDistMax');
if (~found)
    fixDistMax = 6; % 3.2219 for objDistMax = 2m
end

% Muscle initialization [%]: correspond now to the minimum and maximum distance
% the eyes should be looking at. [lateral rectus, medial rectus]
% some correspondances (distance: [lateral, medial] activation):
% 0.5m : [0, 0.0726], 1.5m : [0, 0.0166], 1.5m-2deg : [0, 0.0441], 2m : [0, 0.0089], 3m : [0, 0.0011], 3.22m : [0, 0],
% 4m : [0.0027, 0], 6m : [0.0064, 0], 10m : [0.0093, 0], Inf : [0.0136, 0]

% minimal initial muscle innervation
% orig: 0.00807 corr. to vergAngleMin | 0 corr. to 1 deg
[found, muscleInitMin, myvarargin] = parseparam(myvarargin, 'muscleInitMin');
if (~found)
    muscleInitMin = [0, 0];
end

% maximal initial muscle innervation
% orig: 0.07186 corr. to vergAngleMax | 0.1 corrs. to 12.7 deg
[found, muscleInitMax, myvarargin] = parseparam(myvarargin, 'muscleInitMax');
if (~found)
    muscleInitMax = [0.0064, 0.0166];
end

%%% Factor of Muscle activity feedback in RL feature vector
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
[found, lambdaMuscleFB, myvarargin] = parseparam(myvarargin, 'lambdaMuscleFB');
if (~found)
    lambdaMuscleFB = 0.1269;
end

%%% Reward function parameters, i.e. their "proportions" to the whole reward signal

%%% Reconstruction error factor
% the % is in respect to the reconstruction error, i.e. 100% X means signal X is as strong as
% 6.391 * mean reconstruction error on average, whereby 6.391 * mean reconstruction error ~= 1
% pure 15.647% = 1 | privious 77.12% = 4.929 | 100% = 6.391 | set to 1 for simplicity
%
[found, lambdaRec, myvarargin] = parseparam(myvarargin, 'lambdaRec');
if (~found)
    lambdaRec = 1;
end

% due to recError reduction at dsRatio = [8, 2], lambdaRec needs to be scaled accordingly
% this is a temporary solution, you should recalculate everything if you introduce 3 scales!
if (dsRatio(end) > 1)
    lambdaRec = lambdaRec * 1.7706;
end

%%% Metabolic costs factor (range)
% per_cent = mean(model.metCost_hist) * model.lambdaMet / model.lambdaRec
%
% mean(model.metCost_hist) = 1.0380
% meanR=mean(sum(model.recerr_hist,2)) = 0.1038
% mean(model.metCost_hist) = 0.5727
% 0.5% = 0.0014435 | 0.25% = 0.00072175 | 0.1% = 0.0002887
% 1% = 0.0029 | 2% = 0.0058 | 3% = 0.0087 | 4% = 0.0114 | 5% = 0.0144
% 6% = 0.0173 | 7% = 0.0203 | 8% = 0.0232 | 9% = 0.0261 | 10% = 0.0289 | 12.5% = 0.0360875
% 15% = 0.0435 | 17.5% = 0.0505225 | 20% = 0.05774 | 25% = 0.0722 | 30% = 0.0866 | 50% = 0.1443 | 100% = 0.2887
%%%%%
%
% 0% = 0, 5.2632% = 0.0085, 10.5263% = 0.0171, 15.7895% = 0.0256, 21.0526% = 0.0341, 26.3158% = 0.0427,
% 31.5789% = 0.0512, 36.8421% = 0.0598, 42.1053% = 0.0683, 47.3684% = 0.0768, 52.6316% = 0.0854,
% 57.8947% = 0.0939, 63.1579% = 0.1024, 68.4211% = 0.1110, 73.6842% = 0.1195, 78.9474% = 0.1281,
% 84.2105% = 0.1366, 89.4737% = 0.1451, 94.7368% = 0.1537, 100.0000% = 0.1622
%
[found, metCostRange, myvarargin] = parseparam(myvarargin, 'metCostRange');
if (~found)
    metCostRange = [0, 0];
end

% due to the dependancy of mean(model.metCost_hist) * metCostRange * lambdaRec / mean(recError) * lambdaRec = x%
% metCostRange needs to be scaled accordingly
metCostRange = metCostRange .* lambdaRec; % #hack

if (length(metCostRange) == 1 || metCostRange(1) == metCostRange(2))
    metCostDec = 0; % no decay
elseif (metCostRange(1) < metCostRange(2))
    error('It must hold metCostRange(1) >= metCostRange(2)');
else
    % metCostDec = -(log(2) * trainTime) / log(metCostRange(2) / metCostRange(1)); % metCost decay factor
    metCostDec = metCostRange(1) - metCostRange(2);
end

PARAMModel = {textureFile, trainTime, testAt, sparseCodingType, focalLength, baseline, ...
              objDistMin, objDistMax, muscleInitMin, muscleInitMax, interval, ...
              lambdaMuscleFB, lambdaRec, metCostRange, patchSize, pxFieldOfView, ...
              dsRatio, stride, fixDistMin, fixDistMax, overlap, cutout, metCostDec};

% ------------------------
% Sparce Coding parameters
% ------------------------

% Scales := [coarse, less_coarse, ..., fine], i.e. [peripheral vision, ..., central vision]

% total number of basis | origin [288, 288]
[found, nBasis, myvarargin] = parseparam(myvarargin, 'nBasis');
if (~found)
    nBasis = [400, 400];
end

% number of basis used to encode in sparse mode | origin [10, 10]
[found, nBasisUsed, myvarargin] = parseparam(myvarargin, 'nBasisUsed');
if (~found)
    nBasisUsed = [10, 10];
end

% size of each (binocular) base vector: patchSize * patchSize * 2 (left + right eye) | origin [128, 128]
[found, basisSize, myvarargin] = parseparam(myvarargin, 'basisSize');
if (~found)
    basisSize = [(patchSize ^ 2) * 2, (patchSize ^ 2) * 2];
end

% SC learning rate(s)
% origin 0.01 | Lukas 0.1 | Alex P 0.5, origin 0.01 | Lukas 0.1 | Alex P 0.5 | Chong 0.2
[found, sc_eta, myvarargin] = parseparam(myvarargin, 'sc_eta');
if (~found)
    sc_eta = [0.2, 0.2];
end

% temperature in softmax | origin 0.01
[found, temperature, myvarargin] = parseparam(myvarargin, 'temperature');
if (~found)
    temperature = [0.01, 0.01];
end

if ((length(pxFieldOfViewOrig) ~= length(dsRatio)) ...
 || (length(stride) ~= length(dsRatio)) ...
 || (length(nBasis) ~= length(dsRatio)) ...
 || (length(nBasisUsed) ~= length(dsRatio)) ...
 || (length(basisSize) ~= length(dsRatio)) ...
 || (length(sc_eta) ~= length(dsRatio)) ...
 || (length(temperature) ~= length(dsRatio)))
    error('For usage of #%d scales all respective SC and model parameters needs to have length %d.', ...
          length(dsRatio), length(dsRatio));
end

PARAMSC = {nBasis, nBasisUsed, basisSize, sc_eta, temperature};

% ---------------------------------
% Reinforcement Learning parameters
% ---------------------------------

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
[found, rlFlavour, myvarargin] = parseparam(myvarargin, 'rlFlavour');
if (~found)
    rlFlavour = [uint8(0), uint8(0)];
end

if (length(rlFlavour) ~= 2)
    error('It must hold length(rlFlavour) = 2, as rlFlavour = [Critic, Actor]');
end

if (rlFlavour(1) < 0 || rlFlavour(1) > 1)
    error('rlFlavour(1) must be scalar in {0, 1}');
end

if (rlFlavour(2) < 0 || rlFlavour(2) > 4)
    error('rlFlavour(2) must be scalar in {0, 1, ..., 4}');
end

% indicates if the policy is discrete (= 0) or continuous (= 1)
[found, continuous, myvarargin] = parseparam(myvarargin, 'continuous');
if (~found)
    continuous = uint8(1);
end

% vergence angles (discrete policy) enable half pixel resolution
[found, actionSpace, myvarargin] = parseparam(myvarargin, 'actionSpace');
if (~found)
    actionSpace = [-8, -4, -2, -1, -0.5, -0.2, -0.1, ...
                    0, 0.1, 0.2, 0.5, 1, 2, 4, 8];
end

% Critic learning rate (value function)
% origin 0.05 | Chong 1 | Lukas 0.9 | Alex P 0.4
[found, alpha_v, myvarargin] = parseparam(myvarargin, 'alpha_v');
if (~found)
    alpha_v = 0.75;
end

% CRG Critic discount factor
% origin 0.3 | Alex P 0.3
[found, xi, myvarargin] = parseparam(myvarargin, 'xi');
if (~found)
    xi = 0.3;
end

% CACLA Critic discount factor
% origin 1
% TODO: fuse both discounting factors
[found, gamma, myvarargin] = parseparam(myvarargin, 'gamma');
if (~found)
    gamma = 0.3;
end

% Actor learning rate of natural policy gradient
% origin 0.05 | Chong 0.025 | Lukas 0.1 | Alex P 0.4
[found, alpha_n, myvarargin] = parseparam(myvarargin, 'alpha_n');
if (~found)
    alpha_n = 0.025;
end

% Actor learning rate of Gaussean policy
% origin 1 | Chong 0.002 | Lukas 0.01 | Alex P 0.4 | linear 0.002
[found, alpha_p, myvarargin] = parseparam(myvarargin, 'alpha_p');
if (~found)
    alpha_p = 0.5;
end

% Actor weight regularization via factorial downscaling
[found, regularizer, myvarargin] = parseparam(myvarargin, 'regularizer');
if (~found)
    regularizer = 1 - 1e-3;
end

% variance of action output, i.e. variance of Gaussian policy [training_start, training_end]
% corresponds to softMax temperature in discrete RL models
[found, varianceRange, myvarargin] = parseparam(myvarargin, 'varianceRange');
if (~found)
    varianceRange = [1e-5, 1e-5];
end

if (length(varianceRange) == 1 || varianceRange(1) == varianceRange(2))
    % no variance decay
    varDec = 0;
elseif (varianceRange(1) < varianceRange(2))
    sprintf('Error: It must hold varianceRange(1) >= varianceRange(2)')
    return;
else
    % varDec = -(log(2) * trainTime) / log(varianceRange(2) / varianceRange(1)); % action variance decay factor
    varDec = varianceRange(1) - varianceRange(2);
end

% Actor's number of neurons in the output layer and amount of eye muscles
[found, outputDim, myvarargin] = parseparam(myvarargin, 'outputDim');
if (~found)
    outputDim = 2;
end

% Critic's and Actor's number of neurons in the input layer (Small + Large scale + Muscle activities)
[found, inputDim, myvarargin] = parseparam(myvarargin, 'inputDim');
if (~found)
    if (continuous == 1)
        inputDim = sum(PARAMSC{1}) + outputDim; % number of neurons in the input layer (Small + Large scale + Muscle activities)
    else
        inputDim = sum(PARAMSC{1});             % only small + large scale basis function inputs in discrete models
        varianceRange = 1;
        outputDim = 1;                          % only one delta angle output in discrete models
    end
end

% Actor's number of neurons in the hidden layer
[found, hiddenDim, myvarargin] = parseparam(myvarargin, 'hiddenDim');
if (~found)
    hiddenDim = 50;
end

% all network dimensions at once
[found, dimensions, myvarargin] = parseparam(myvarargin, 'dimensions');
if (~found)
    dimensions = [inputDim, hiddenDim, outputDim];
end

% initial network weights
[found, weight_range, myvarargin] = parseparam(myvarargin, 'weight_range');
if (~found)
    weight_range = [1 / inputDim, ...                   % maximum initial weight [critic_ji, actor_ji, actor_kj]
                    1 / (inputDim * hiddenDim), ...     % origin [0.05, 0.4, 0.4] | Lukas [0.1, 0.05, 0.05] | AL 100 / (inputDim * hiddenDim)
                    2 / (hiddenDim * outputDim)];       % linear [1/inputDim, 1/inputDim, -]
end

% Actor's reguralization factor
% origin 0.01
% TODO: DEPRICATED because it became obsolet, cleanup needed
[found, lambda, myvarargin] = parseparam(myvarargin, 'lambda');
if (~found)
    lambda = 0.01;
end

% TD error variance tracking/approximating (CACLAVar)
[found, deltaVar, myvarargin] = parseparam(myvarargin, 'deltaVar');
if (~found)
    deltaVar = 1;
end

% TD error variance scaling factor (CACLAVar)
[found, rl_eta, myvarargin] = parseparam(myvarargin, 'rl_eta');
if (~found)
    rl_eta = 0.001;
end

% scaling factor of Fisher Information matrix (CNGFI)
[found, fiScale, myvarargin] = parseparam(myvarargin, 'fiScale');
if (~found)
    fiScale = 1e-5;
end

PARAMRL = {actionSpace, alpha_v, alpha_n, alpha_p, xi, gamma, varianceRange, lambda, dimensions, weight_range, ...
           continuous, deltaVar, rl_eta, fiScale, rlFlavour, varDec, regularizer};

PARAM = {PARAMModel, PARAMSC, PARAMRL};
model = Model(PARAM);

% -------------------------
% Parameter vector handling
% -------------------------

% Parses parameter vector and prunes it on successful find
% @param paramVector:    parameter vector
% @param param:          searched paramter name
%
% return found:         [0, 1] success flag
% return val:           value of parameter
% return paramVector:   updated parameter vector
function [found, val, paramVector] = parseparam(paramVector, param)

% search param name
isvar = cellfun(@(x) ischar(x) && strcmpi(x, param), paramVector);

if (sum(isvar) > 1)
    error('Parameters can only be passed once');
end

if (any(isvar))
    found = true;
    idx = find(isvar);
    val = paramVector{idx + 1};
    paramVector([idx, idx + 1]) = []; % prune param vector
else
    found = false;
    val = [];
end
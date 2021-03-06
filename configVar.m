% Model object configuration and generation with support of variable
% parameter vector
%
% @param varParamArray: variable parameter cell array
%                       {'paramName1', value1, 'paramName2', value2, ...}
% return model:         handle (pointer) to generated model instance
function model = configVar(varParamArray)

% ---------------------
% Experiment parameters
% ---------------------

% save parameter vector for documentation (model.varParamArray)
inputParams = varParamArray;

% stimulus file name
[found, textureFile, varParamArray] = parseparam(varParamArray, 'textureFile');
if (~found)
    textureFile = {'Textures_mcgillManMade40.mat', 'Textures_mcgillManMade100.mat'};
end

% training duration
[found, trainTime, varParamArray] = parseparam(varParamArray, 'trainTime');
if (~found)
    trainTime = 500000;
    % trainTime = 1000000; % for training with metabolic costs
end

if (~isscalar(trainTime) || trainTime < 1)
    error('trainTime must be scalar > 0');
end

% points in time of intermediate test procedure during training
[found, testAt, varParamArray] = parseparam(varParamArray, 'testAt');
if (~found)
    testAt = [250000 : 250000 : trainTime];
    % testAt = [500000 : 500000 : trainTime];
    % testAt = [100000 : 100000 : trainTime]; % for training without metCosts
end

% sparse coding type
[found, sparseCodingType, varParamArray] = parseparam(varParamArray, 'sparseCodingType');
if (~found)
    sparseCodingType = uint8(0);
end

if (~isscalar(sparseCodingType) || sparseCodingType < 0 || sparseCodingType > 1)
    error('sparseCodingType must be scalar in {0, 1}');
end

% ----------------
% Model parameters
% ----------------

% period for changing the stimulus for the eyes at training | origin 10
[found, interval, varParamArray] = parseparam(varParamArray, 'interval');
if (~found)
    interval = 10;
end

if (~isscalar(interval) || interval < 1)
    error('interval must be scalar >= 1');
end

% period for changing the stimulus for the eyes at testing | origin 10
[found, testInterval, varParamArray] = parseparam(varParamArray, 'testInterval');
if (~found)
    testInterval = interval * 2;
end

if (~isscalar(testInterval) || testInterval < 2)
    error('testInterval must be scalar >= 2');
end

%%% Image processing variables

% patch size [pixel] of one basis functon, i.e. "receptive field" size | origin 8
[found, patchSize, varParamArray] = parseparam(varParamArray, 'patchSize');
if (~found)
    patchSize = 8;
end

if (~isscalar(patchSize) || patchSize < 4)
    error('patchSize must be scalar >= 4');
end

% [peripheral, intermediate ..., central vision]
% downsampling ratio, i.e. how many pixels in original image
% correspond to one pixels in downsampled image | origin [8, 2]
[found, dsRatio, varParamArray] = parseparam(varParamArray, 'dsRatio');
if (~found)
    dsRatio = [4, 1];
end

% fields of view in original image [pixel] | origin [128, 80]
[found, pxFieldOfViewOrig, varParamArray] = parseparam(varParamArray, 'pxFieldOfViewOrig');
if (~found)
    pxFieldOfViewOrig = [128, 40];
end

% fields of view in downsampled image [pixel] (previously called fovea)
% pxFieldOfView = FieldOfView in original image [pixel] / dsRatio
[found, pxFieldOfView, varParamArray] = parseparam(varParamArray, 'pxFieldOfView');
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
[found, stride, varParamArray] = parseparam(varParamArray, 'stride');
if (~found)
    stride = [patchSize / 2, patchSize / 2];
end

% flag indicates whether cutout procedure is applied [1] or not [0]
[found, cutout, varParamArray] = parseparam(varParamArray, 'cutout');
if (~found)
    cutout = uint8(0);
end

if (~isscalar(cutout) || cutout < 0 ||  cutout > 1)
    error('cutout must be scalar in {0, 1}');
end

% overlap between the different layers measured in units of FINE scale,
% works only in conjunction with cutout
[found, overlap, varParamArray] = parseparam(varParamArray, 'overlap');
if (~found)
    overlap = 0;
end

if (length(dsRatio) > 1 && (length(overlap) ~= length(dsRatio) - 1))
    error('For usage of #%d scales overlap needs to have length %d.', length(dsRatio) - 1);
end

[found, whitening, varParamArray] = parseparam(varParamArray, 'whitening');
if (~found)
    whitening = 0;
end

%%% Camera parameters

% vertical offset between left and right (0 in the iCub Simulator!)
% [found, offset, varParamArray] = parseparam(varParamArray, 'offset');
% if (~found)
%     offset = 0;
% end

% focal length [px]
[found, focalLength, varParamArray] = parseparam(varParamArray, 'focalLength');
if (~found)
    focalLength = 257.34;
end

% interocular distance [m]
[found, baseline, varParamArray] = parseparam(varParamArray, 'baseline');
if (~found)
    baseline = 0.056;
end

% Object distance to eyes [m]
[found, objDistMin, varParamArray] = parseparam(varParamArray, 'objDistMin');
if (~found)
    objDistMin = 0.5; % origin 0.5
end

[found, objDistMax, varParamArray] = parseparam(varParamArray, 'objDistMax');
if (~found)
    objDistMax = 6;   % origin 2
end

% Fixation distance [m]
% used for eye fixation initialization
[found, fixDistMin, varParamArray] = parseparam(varParamArray, 'fixDistMin');
if (~found)
    fixDistMin = 0.3379;
end

[found, fixDistMax, varParamArray] = parseparam(varParamArray, 'fixDistMax');
if (~found)
    fixDistMax = 10; % 3.2219 for objDistMax = 2m
end

% object plane size where the stimuli are projected on
[found, objSize, varParamArray] = parseparam(varParamArray, 'objSize');
if (~found)
    objSize = 3; % 3 m fills coarse and fine scale at 6 m
end

% maximal yaw angle of object plane
[found, maxYaw, varParamArray] = parseparam(varParamArray, 'maxYaw');
if (~found)
    maxYaw = 0;
end

% maximal tilt angle of object plane
[found, maxTilt, varParamArray] = parseparam(varParamArray, 'maxTilt');
if (~found)
    maxTilt = 0;
end

% maximal roll angle of object plane
[found, maxRoll, varParamArray] = parseparam(varParamArray, 'maxRoll');
if (~found)
    maxRoll = 0;
end

% Muscle initialization [%]: correspond now to the minimum and maximum distance
% the eyes should be looking at. [lateral rectus, medial rectus]
% some correspondances (distance: [lateral, medial] activation):
% 0.5m : [0, 0.0726], 1.5m : [0, 0.0166], 1.5m-2deg : [0, 0.0441], 2m : [0, 0.0089], 3m : [0, 0.0011], 3.22m : [0, 0],
% 4m : [0.0027, 0], 6m : [0.0064, 0], 10m : [0.0093, 0], Inf : [0.0136, 0]

% minimal initial muscle innervation
% orig: 0.00807 corr. to vergAngleMin | 0 corr. to 1 deg
[found, muscleInitMin, varParamArray] = parseparam(varParamArray, 'muscleInitMin');
if (~found)
    muscleInitMin = [0, 0];
end

% maximal initial muscle innervation
% orig: 0.07186 corr. to vergAngleMax | 0.1 corrs. to 12.7 deg
[found, muscleInitMax, varParamArray] = parseparam(varParamArray, 'muscleInitMax');
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
% standard before normalization: 0.1269
[found, lambdaMuscleFB, varParamArray] = parseparam(varParamArray, 'lambdaMuscleFB');
if (~found)
    lambdaMuscleFB = 1;
end

%%% Reward function parameters, i.e. their "proportions" to the whole reward signal

%%% Reconstruction error factor
% the % is in respect to the reconstruction error, i.e. 100% X means signal X is as strong as
% 6.391 * mean reconstruction error on average, whereby 6.391 * mean reconstruction error ~= 1
% pure 15.647% = 1 | privious 77.12% = 4.929 | 100% = 6.391 | set to 1 for simplicity
%
[found, lambdaRec, varParamArray] = parseparam(varParamArray, 'lambdaRec');
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
[found, metCostRange, varParamArray] = parseparam(varParamArray, 'metCostRange');
if (~found)
    metCostRange = [0, 0];
end

% due to the dependancy of mean(model.metCost_hist) * metCostRange * lambdaRec / mean(recError) * lambdaRec = x%
% metCostRange needs to be scaled accordingly
metCostRange = metCostRange .* lambdaRec; % #hack

% Metabolic costs decay factor
if (length(metCostRange) == 1 || metCostRange(1) == metCostRange(2))
    metCostDec = 0; % no decay
elseif (metCostRange(1) < metCostRange(2))
    error('It must hold metCostRange(1) >= metCostRange(2)');
else
    % metCostDec = -(log(2) * trainTime) / log(metCostRange(2) / metCostRange(1));  % exponential decay factor
    metCostDec = metCostRange(1) - metCostRange(2);                                 % linear decay factor
end

%%% muscle initialization / reset method
% 0 ('simple')      : use random muscle commands out of [0, 0.1] for lateral and [0,
%                     0.2] for medial rectus
% 1 ('advanced')    : first, uniformly draw a object distance, then set one
%                     muscle actity to 0 and calculate the other muscle activity
% 2 ('moreAdvanced'): first, uniformly draw a object distance, then uniformly draw
%                     a point from all possible mucle commands that fixate this distance
% 3 ('perturbed')   : take the last command and add a random vector with radius
%                     uniformly drawn from [0, 0.02] (max 2 % muscle activity for both muscles)
% 4 ('laplacian')   : draw fixation distances that produce laplacian distributed disparities
%%%
[found, initMethod, varParamArray] = parseparam(varParamArray, 'initMethod');
if (~found)
    initMethod = 2;
end

% Keep(0) or normalize(1) feature fector by z-transform
[found, normFeatVect, varParamArray] = parseparam(varParamArray, 'normFeatVect');
if (~found)
    normFeatVect = 1;
end

if ((normFeatVect ~= 0) && (normFeatVect ~= 1))
    error('normFeatVect must be a scalar in {0, 1}');
end

% Desired standard deviation of each z-transformed variable
[found, desiredStdZT, varParamArray] = parseparam(varParamArray, 'desiredStdZT');
if (~found)
    desiredStdZT = 0.02;
end

if (desiredStdZT <= 0)
    error('desiredStdZT must be a scalar in ]0, inf[');
end

[found, lapSig, varParamArray] = parseparam(varParamArray, 'lapSig');
if (~found)
    lapSig = 1;     % if initMethod == 4, this produces disparities mostly between -2 and 2 deg
    % lapSig = 0.5;   % mostly between [-1, 1]
end

[found, strabAngle, varParamArray] = parseparam(varParamArray, 'strabAngle');
if (~found)
    strabAngle = 0;
end

%%% filters for left and right images to simulate altered rearing conditions
%   and the probabilities at which they should occur.
%   Leave empty if no filters should be applied.
[found, filterLeft, varParamArray] = parseparam(varParamArray, 'filterLeft');
if (~found)
    filterLeft = [];
    filterRight = []; % hotfix for certain parameters
elseif filterLeft == 1
    % filterLeft = orientedGaussian(9,9,0.1); % leaves only vertical edges
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(9,9,0.1); % leaves only vertical edges
elseif filterLeft == 2
    % filterLeft = orientedGaussian(9,0.1,9); % leaves only horizontal edges
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(9,0.1,9);
elseif filterLeft == 3
    %filterLeft = orientedGaussian(240,240,240); % blurrs the whole image (monocular deprivation)
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(240,240,240);
elseif filterLeft == 4
    %filterLeft = orientedGaussian(17,17,0.1);
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(17,17,0.1);
elseif filterLeft == 5
    %filterLeft = orientedGaussian(17,0.1,17);
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(17,0.1,17);
elseif filterLeft == 6
    % filterLeft = orientedGaussian(33,33,0.1);
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,33,0.1); % only verticals
elseif filterLeft == 7
    % filterLeft = orientedGaussian(33,0.1,33);
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,0.1,33); % only horizontals
elseif filterLeft == 8
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,0.1,0.1);
elseif filterLeft == 9
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,0.2,0.2);
elseif filterLeft == 10
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,0.5,0.5);
elseif filterLeft == 11
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,1,1);
elseif filterLeft == 12
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,2,2);
elseif filterLeft == 13
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,5,5);
elseif filterLeft == 14
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(100,1,1);
elseif filterLeft == 15
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(100,5,5);
elseif filterLeft == 16
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(100,10,10);
elseif filterLeft == 17
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(100,50,50);
elseif filterLeft == 18
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(100,100,100);
elseif filterLeft == 19
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(240,0.1,33); % only horizontals
elseif filterLeft == 20
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(240,0.1,100); % only horizontals
elseif filterLeft == 21
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(240,0.1,240); % only horizontals
elseif filterLeft == 22
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(240,0.1,1000); % only horizontals
elseif filterLeft == 23
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,0.1,6); % only horizontals
elseif filterLeft == 24
      filterLeft = {};
      [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,0.1,10); % only horizontals
elseif filterLeft == 25
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(33,0.1,1); % only horizontals
%%%%%%%%%%%%%%%% Journal Paper Params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filterLeft == 26 % vertical
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(240,240,0.1); % only verticals
    filterRight = filterLeft;
elseif filterLeft == 27 % horizontal
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(240,0.1,240); % only horizontals
    filterRight = filterLeft;
elseif filterLeft == 28 % orthogonal
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(240,240,0.1); % only verticals
    filterRight = {};
    [filterRight{1}, filterRight{2}] = orientedGaussianVectors(240,0.1,240); % only horizontals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filterLeft == 29 % normal case
    % add a little blurr in both eyes to make up for the non-zero distortion
    % in the other cases
    sdd = 0.1;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(10,0.1,0.1);
    filterRight = filterLeft;
%%%%%%%%%%%%%%%%%%%%%%%%% vertical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filterLeft == 30 % vertical 33
    sdd = 33;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,sdd,0.1);
    filterRight = filterLeft;
elseif filterLeft == 31 %vertical 100
    sdd = 100;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,sdd,0.1);
    filterRight = filterLeft;
elseif filterLeft == 32 % vertical 240
    sdd = 240;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,sdd,0.1);
    filterRight = filterLeft;
elseif filterLeft == 33 %vertical 500
    sdd = 500;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,sdd,0.1);
    filterRight = filterLeft;
%%%%%%%%%%%%%%%%%%%%%%%% horizontal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filterLeft == 34 % horizontal 33
    sdd = 33;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,0.1,sdd);
    filterRight = filterLeft;
elseif filterLeft == 35 % horizontal 100
    sdd = 100;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,0.1,sdd);
    filterRight = filterLeft;
elseif filterLeft == 36 % horizontal 240
    sdd = 240;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,0.1,sdd);
    filterRight = filterLeft;
elseif filterLeft == 37 % horizontal 500
    sdd = 500;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,0.1,sdd);
    filterRight = filterLeft;
%%%%%%%%%%%%%%%%%%% finer grained filter size %%%%%%%%%%%%%%%%%%%%%%
elseif filterLeft == 38 % vertical 9
    sdd = 9;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,sdd,0.1);
    filterRight = filterLeft;
elseif filterLeft == 39 % vertical 17
    sdd = 17;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,sdd,0.1);
    filterRight = filterLeft;
elseif filterLeft == 40 % horizontal 9
    sdd = 9;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,0.1,sdd);
    filterRight = filterLeft;
elseif filterLeft == 41 % horizontal 17
    sdd = 17;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,0.1,sdd);
    filterRight = filterLeft;
%%%%%%%%%%%%%%%%%%% orthogonal filters %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filterLeft == 42
    sdd = 33;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,0.1,sdd);
    filterRight = {};
    [filterRight{1}, filterRight{2}] = orientedGaussianVectors(6*sdd,sdd,0.1);
elseif filterLeft == 43
    sdd = 17;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,0.1,sdd);
    filterRight = {};
    [filterRight{1}, filterRight{2}] = orientedGaussianVectors(6*sdd,sdd,0.1);
elseif filterLeft == 44
    sdd = 100;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,0.1,sdd);
    filterRight = {};
    [filterRight{1}, filterRight{2}] = orientedGaussianVectors(6*sdd,sdd,0.1);
%%%%%%%%%%%%%%%%%% monocular filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filterLeft == 45
    sdd = 33;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,0.1,0.1);
    filterRight = {};
    [filterRight{1}, filterRight{2}] = orientedGaussianVectors(6*sdd,sdd,sdd);
%%%%%%%%%%%%%%%%%%% orthogonal filters %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filterLeft == 46
    sdd = 33;
    filterLeft = {};
    [filterLeft{1}, filterLeft{2}] = orientedGaussianVectors(6*sdd,sdd,0.1);
    filterRight = {};
    [filterRight{1}, filterRight{2}] = orientedGaussianVectors(6*sdd,0.1,sdd);
end

[found, filterLeftProb, varParamArray] = parseparam(varParamArray, 'filterLeftProb');
if (~found)
    filterLeftProb = 1;
end

% just a workaround to test multiple different same filters for left and
% right eyes
% filterRight = filterLeft;
% filterRightProb = filterLeftProb;

% [found, filterRight, varParamArray] = parseparam(varParamArray, 'filterRight');
% if (~found)
%     filterRight = filterLeft;
% elseif filterRight == 1
%     % filterRight = orientedGaussian(240,240,240);
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(240,33,33);
%     % [filterRight{1}, filterRight{2}] = orientedGaussianVectors(240,1000,1000);
% elseif filterRight == 2
%     % filterRight = orientedGaussian(240,240,240);
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(240,240,240);
%     % [filterRight{1}, filterRight{2}] = orientedGaussianVectors(240,1000,1000);
% elseif filterRight == 3
%     % filterRight = orientedGaussian(240,240,240);
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(240,1000,1000);
%     % [filterRight{1}, filterRight{2}] = orientedGaussianVectors(240,1000,1000);
% end

[found, filterRightProb, varParamArray] = parseparam(varParamArray, 'filterRightProb');
if (~found)
    filterRightProb = filterLeftProb;
end

% [found, filterRight, varParamArray] = parseparam(varParamArray, 'filterRight');
% if (~found)
%     filterRight = [];
% end
% if filterRight == 1
%     % filterRight = orientedGaussian(9,9,0.1);
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(9,9,0.1);
% elseif filterRight == 2
%     % filterRight = orientedGaussian(9,0.1,9);
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(9,0.1,9);
% elseif filterRight == 3
%     % filterRight = orientedGaussian(240,240,240);
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(240,240,240);
% elseif filterRight == 4
%     % filterRight = orientedGaussian(17,17,0.1);
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(17,17,0.1);
% elseif filterRight == 5
%     % filterRight = orientedGaussian(17,0.1,17);
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(17,0.1,17);
% elseif filterRight == 6
%     % filterRight = orientedGaussian(33,33,0.1);
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(33,33,0.1);
% elseif filterRight == 7
%     % filterRight = orientedGaussian(33,0.1,33);
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(33,0.1,33);
% elseif filterRight == 8
%     filterRight = {};
%     [filterRight{1}, filterRight{2}] = orientedGaussianVectors(33,0.1,0.1);
% end
%
% [found, filterRightProb, varParamArray] = parseparam(varParamArray, 'filterRightProb');
% if (~found)
%     filterRightProb = 1;
% end

PARAMModel = {textureFile, trainTime, testAt, sparseCodingType, focalLength, baseline, ...
              objDistMin, objDistMax, muscleInitMin, muscleInitMax, interval, ...
              lambdaMuscleFB, lambdaRec, metCostRange, patchSize, pxFieldOfView, ...
              dsRatio, stride, fixDistMin, fixDistMax, overlap, cutout, metCostDec, ...
              initMethod, inputParams, normFeatVect, desiredStdZT, testInterval, ...
              filterLeft, filterLeftProb, filterRight, filterRightProb, ...
              whitening, lapSig, strabAngle, objSize, maxYaw, maxTilt, maxRoll};

% ------------------------
% Sparse Coding parameters
% ------------------------

% Scales := [coarse, less_coarse, ..., fine], i.e. [peripheral vision, ..., central vision]

% total number of basis | origin [288, 288]
[found, nBasis, varParamArray] = parseparam(varParamArray, 'nBasis');
if (~found)
    nBasis = [400, 400];
end

% number of basis used to encode in sparse mode | origin [10, 10]
[found, nBasisUsed, varParamArray] = parseparam(varParamArray, 'nBasisUsed');
if (~found)
    nBasisUsed = [10, 10];
end

% size of each (binocular) base vector: patchSize * patchSize * 2 (left + right eye) | origin [128, 128]
[found, basisSize, varParamArray] = parseparam(varParamArray, 'basisSize');
if (~found)
    basisSize = [(patchSize ^ 2) * 2, (patchSize ^ 2) * 2];
end

% SC learning rate(s)
% origin 0.01 | Lukas 0.1 | Alex P 0.5, origin 0.01 | Lukas 0.1 | Alex P 0.5 | Chong 0.2
[found, sc_eta, varParamArray] = parseparam(varParamArray, 'sc_eta');
if (~found)
    sc_eta = [0.2, 0.2];
end

% temperature in softmax | origin 0.01
[found, temperature, varParamArray] = parseparam(varParamArray, 'temperature');
if (~found)
    temperature = [0.01, 0.01];
end

% initialization of basis functions
[found, BFinit, varParamArray] = parseparam(varParamArray, 'BFinit');
if (~found)
    BFinit = [2, 2];     % initialization with non-aligned Gabor wavelets
    % BFinit = [1, 1];   % initialization with random white noise
end

% fitting of basis functions
[found, fitFreq, varParamArray] = parseparam(varParamArray, 'BFfitFreq');
if (~found)
    BFfitFreq = 1;     % fitting of basis functions uses frequency instead of wavelength
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

PARAMSC = {nBasis, nBasisUsed, basisSize, sc_eta, temperature, BFinit, [BFfitFreq, BFfitFreq]};

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
[found, rlFlavour, varParamArray] = parseparam(varParamArray, 'rlFlavour');
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
[found, continuous, varParamArray] = parseparam(varParamArray, 'continuous');
if (~found)
    continuous = uint8(1);
end

% vergence angles (discrete policy) enable half pixel resolution
[found, actionSpace, varParamArray] = parseparam(varParamArray, 'actionSpace');
if (~found)
    actionSpace = [-8, -4, -2, -1, -0.5, -0.2, -0.1, ...
                    0, 0.1, 0.2, 0.5, 1, 2, 4, 8];
end

% Critic learning rate (value function)
% origin 0.05 | Chong 1 | Lukas 0.9 | Alex P 0.4
% [found, alpha_v, varParamArray] = parseparam(varParamArray, 'alpha_v');
% if (~found)
%     alpha_v = 0.75;
% end

% Critic learning rate range (value function)
% origin 0.05 | Chong 1 | Lukas 0.9 | Alex P 0.4
[found, criticLRRange, varParamArray] = parseparam(varParamArray, 'criticLRRange');
if (~found)
    criticLRRange = [0.75];
end

% Critic learning rate decay factor
if (length(criticLRRange) == 1 || criticLRRange(1) == criticLRRange(2))
    critLRDec = 0; % no variance decay
elseif (criticLRRange(1) < criticLRRange(2))
    error('It must hold criticLRRange(1) >= criticLRRange(2)');
else
    % critLRDec = -(log(2) * trainTime) / log(criticLRRange(2) / criticLRRange(1)); % exponential decay factor
    critLRDec = criticLRRange(1) - criticLRRange(2);                                % linear decay factor
end

% CRG Critic discount factor
% origin 0.3 | Alex P 0.3
[found, xi, varParamArray] = parseparam(varParamArray, 'xi');
if (~found)
    xi = 0.3;
end

% CACLA Critic discount factor
% origin 1
% TODO: fuse both discounting factors
[found, gamma, varParamArray] = parseparam(varParamArray, 'gamma');
if (~found)
    gamma = 0.3;
end

% Actor learning rate of natural policy gradient
% origin 0.05 | Chong 0.025 | Lukas 0.1 | Alex P 0.4
[found, alpha_n, varParamArray] = parseparam(varParamArray, 'alpha_n');
if (~found)
    alpha_n = 0.025;
end

% % Actor learning rate of Gaussean policy
% % origin 1 | Chong 0.002 | Lukas 0.01 | Alex P 0.4 | linear 0.002
% [found, alpha_p, varParamArray] = parseparam(varParamArray, 'alpha_p');
% if (~found)
%     alpha_p = 0.5;
% end

% Actor learning rate of Gaussean policy
% origin 1 | Chong 0.002 | Lukas 0.01 | Alex P 0.4 | linear 0.002
[found, actorLRRange, varParamArray] = parseparam(varParamArray, 'actorLRRange');
if (~found)
    actorLRRange = [0.5, 0];
end

% Actor learning rate decay factor
if (length(actorLRRange) == 1 || actorLRRange(1) == actorLRRange(2))
    actLRDec = 0; % no variance decay
elseif (actorLRRange(1) < actorLRRange(2))
    error('It must hold actorLRRange(1) >= actorLRRange(2)');
else
    % actLRDec = -(log(2) * trainTime) / log(actorLRRange(2) / actorLRRange(1)); % exponential decay factor
    actLRDec = actorLRRange(1) - actorLRRange(2);                                % linear decay factor
end

% Actor weight regularization via factorial downscaling
% how it works: actor.wp_ji = (1 - (actor.actorRegularizer * actor.learnRate)) * actor.wp_ji;
% previously 1e-4 | 1e-3 / actorLRRange(1); ensures a regularization factor of 1-1e-3 at the beginning of the simulation.
[found, actorRegularizer, varParamArray] = parseparam(varParamArray, 'actorRegularizer');
if (~found)
    actorRegularizer = 1e-5;
end

% critic weight regularization via factorial downscaling
% works as in the actor
[found, criticRegularizer, varParamArray] = parseparam(varParamArray, 'criticRegularizer');
if (~found)
    criticRegularizer = 0;
end

% Variance of action output, i.e. variance of Gaussian policy [training_start, training_end]
% corresponds to softMax temperature in discrete RL models
[found, varianceRange, varParamArray] = parseparam(varParamArray, 'varianceRange');
if (~found)
    varianceRange = [1e-5, 1e-5];
end

% Action variance decay factor
if (length(varianceRange) == 1 || varianceRange(1) == varianceRange(2))
    % no variance decay
    varDec = 0;
elseif (varianceRange(1) < varianceRange(2))
    error('It must hold varianceRange(1) >= varianceRange(2)');
else
    % varDec = -(log(2) * trainTime) / log(varianceRange(2) / varianceRange(1)); % exponential decay factor
    varDec = varianceRange(1) - varianceRange(2);                                % linear decay factor
end

% Actor`s number of neurons in the output layer and amount of eye muscles
[found, outputDim, varParamArray] = parseparam(varParamArray, 'outputDim');
if (~found)
    outputDim = 2;
end

% Constant bias in input vector/feature vector (0 := off, >0 := take that value as bias)
[found, bias, varParamArray] = parseparam(varParamArray, 'bias');
if (~found)
    bias = desiredStdZT; % 0 before, 0.02 in curr. approach
end

% Critic`s and Actor`s number of neurons in the input layer (Small + Large scale + Muscle activities)
% note: if you specify a bias and define the input dimension, remember to
% add the bias to that!
[found, inputDim, varParamArray] = parseparam(varParamArray, 'inputDim');
if (~found)
    if (continuous == 1)
        inputDim = sum(PARAMSC{1}) + outputDim; % number of neurons in the input layer (Small + Large scale + Muscle activities) 802
    else
        inputDim = sum(PARAMSC{1});             % only small + large scale basis function inputs in discrete models
        varianceRange = 1;
        outputDim = 1;                          % only one delta angle output in discrete models
    end
    if (bias > 0)
        inputDim = inputDim + 1;
    end
end

% Actor`s number of neurons in the hidden layer
[found, hiddenDim, varParamArray] = parseparam(varParamArray, 'hiddenDim');
if (~found)
    hiddenDim = 50;
end

% all network dimensions at once
[found, dimensions, varParamArray] = parseparam(varParamArray, 'dimensions');
if (~found)
    dimensions = [inputDim, hiddenDim, outputDim];
end

% initial network weights
[found, weight_range, varParamArray] = parseparam(varParamArray, 'weight_range');
if (~found)
    weight_range = [1 / inputDim, ...                   % maximum initial weight [critic_ji, actor_ji, actor_kj]
                    1 / (inputDim * hiddenDim), ...     % origin [0.05, 0.4, 0.4] | Lukas [0.1, 0.05, 0.05] | AL 100 / (inputDim * hiddenDim)
                    1 / (hiddenDim * outputDim)];       % linear [1/inputDim, 1/inputDim, -]
end

% Actor`s reguralization factor
% origin 0.01
% TODO: DEPRICATED because it became obsolet, cleanup needed
[found, lambda, varParamArray] = parseparam(varParamArray, 'lambda');
if (~found)
    lambda = 0.01;
end

% TD error variance tracking/approximating (CACLAVar)
[found, deltaVar, varParamArray] = parseparam(varParamArray, 'deltaVar');
if (~found)
    deltaVar = 1;
end

% TD error variance scaling factor (CACLAVar)
[found, rl_eta, varParamArray] = parseparam(varParamArray, 'rl_eta');
if (~found)
    rl_eta = 0.001;
end

% scaling factor of Fisher Information matrix (CNGFI)
[found, fiScale, varParamArray] = parseparam(varParamArray, 'fiScale');
if (~found)
    fiScale = 1e-5;
end

PARAMRL = {actionSpace, criticLRRange, alpha_n, actorLRRange, xi, gamma, varianceRange, lambda, dimensions, weight_range, ...
           continuous, deltaVar, rl_eta, fiScale, rlFlavour, varDec, actorRegularizer, critLRDec, actLRDec, bias, criticRegularizer};

if (~isempty(varParamArray))
    error('varParamArray contains unrecognized elements, p.e. %s', varParamArray{1});
end

PARAM = {PARAMModel, PARAMSC, PARAMRL};
model = Model(PARAM);

% -------------------------
% Parameter vector handling
% -------------------------

% Parses parameter vector and prunes it on successful find
% @param paramVector:   parameter vector
% @param param:         searched paramter name
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

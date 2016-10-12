if isempty(gcp('nocreate'))
    mPool = parpool(); % creates a parallel pool with a default number of workers that chould be specified by the cluster this script in running on
end

% upperBounds = [1e-4, 1e-5];
% lowerBounds = [5e-6, 1e-6];
% params = combvec(upperBounds, lowerBounds);
% nParams = length(params)

varNames = ['varDec', 'alpha_p'];
var1 = [[1e-4, 5e-6]; [1e-4, 1e-6]; [1e-5, 5e-6]; [1e-5, 1e-6]];	% variance decay boundaries
var2 = [0.5; 0.3];													% for example actor learnig rate 

nParams = length(var1) * length(var2); 


%% small test
% results = zeros(2, nParams);
% 
% parfor (ind = 1 : nParams)
% 	results(:, ind) = params(:, ind)
% end
% 
% results

%% parameter section
nIters = 2000000;   
rSeed = 1;                      % default
% textureFiles = {'Textures_mcgillManMade40.mat', 'Textures_mcgillManMade100.mat'};
% textureFiles = {'mcGillTest2.mat', 'mcGillTest1.mat'}; % test files containing less images

% simulator = prepareSimulator(textureFiles); % at first, we try to use a shared simulator for all threads
% additional params, that should be set in the config
paramIdentifiers = [];
params = [];

%% main loop, TODO: test if renderers are not interfering (anaglyphs ...)
parfor ind = 1 : nParams
    OES2Muscles(nIters, rSeed, paramIdentifiers, params(ind), sprintf('varDec%g--%g', params(ind)));
end
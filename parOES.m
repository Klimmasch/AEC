if isempty(gcp('nocreate'))
    mPool = parpool(); % creates a parallel pool with a default number of workers that chould be specified by the cluster this script in running on
end

% upperBounds = [1e-4, 1e-5];
% lowerBounds = [5e-6, 1e-6];
% params = combvec(upperBounds, lowerBounds);
% nParams = length(params)

% varNames = {'varDec', 'alpha_p'};
% var1 = {[1e-4, 5e-6]; [1e-4, 1e-6]};	% variance decay boundaries
% var2 = {0.5; 0.3};													% for example actor learnig rate 

varNames = {'interval', 'alpha_v'};
var1 = {10, 20, 50, 100};
var2 = {0.75};

nParams = length(var1) * length(var2); 
PARAMS = cell(nParams, 2);				%column vector, each row contains a set of params for one run

dummy = 1;
for i = 1 : length(var1)
	for j = 1 : length(var2)
		PARAMS(dummy, :) = {cell2mat(var1(i)), cell2mat(var2(j))};
		dummy = dummy + 1;
	end
end

%% small test
% results = zeros(2, nParams);
% 
% parfor (ind = 1 : nParams)
% 	results(:, ind) = params(:, ind)
% end
% 
% results

%% parameter section
nIters = 500000; % 2000000;   
rSeed = 1;                      % default
% textureFiles = {'Textures_mcgillManMade40.mat', 'Textures_mcgillManMade100.mat'};
% textureFiles = {'mcGillTest2.mat', 'mcGillTest1.mat'}; % test files containing less images

% simulator = prepareSimulator(textureFiles); % at first, we try to use a shared simulator for all threads
% additional params, that should be set in the config
% paramIdentifiers = ['alpha_v'];
% params = [0.75];

%% main loop, TODO: test if renderers are not interfering (anaglyphs ...)
parfor ind = 1 : nParams
    OES2Muscles(nIters, rSeed, PARAMS(ind), sprintf('testCluster_interval%d', var1{ind})); % sprintf('varDec%g--%g', PARAMS{ind})
end
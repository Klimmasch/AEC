if isempty(gcp('nocreate'))
    mPool = parpool(); % creates a parallel pool with a default number of workers that chould be specified by the cluster this script in running on
end

upperBounds = [1e-4, 1e-5];
lowerBounds = [5e-6, 1e-6];

params = combvec(upperBounds, lowerBounds);
nParams = length(params)

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

simulator = prepareSimulator(); % at first, we try to use a shared simulator for all threads

%% main loop
parfor ind = 1 : nParams
    
    OES2Muscles(nIters, rSeed, simulator, 'paramIdentifier', params(ind), sprintf('varDec%g--%g', params(ind));
end



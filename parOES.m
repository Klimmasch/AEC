%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                                                                                 %%%%%%%%%%%
%%%%%%%%%%%                 parallel Matlab simulation with Open Eye Simulator              %%%%%%%%%%%
%%%%%%%%%%%                         Author: Lukas Klimmasch                                 %%%%%%%%%%%
%%%%%%%%%%%                                                                                 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Usage:  This scprit creates a pool of parallel workers and distributes          %%%%%%%%%%%
%%%%%%%%%%%         a number of jobs with different parameter combinations provided         %%%%%%%%%%%
%%%%%%%%%%%         in the first part of the file.                                          %%%%%%%%%%%
%%%%%%%%%%%         It is possible to start a series of experiments with all                %%%%%%%%%%%
%%%%%%%%%%%         combinations of two variables that can take different ranges.           %%%%%%%%%%%
%%%%%%%%%%%         To faciliate readability of resulting file names, one can               %%%%%%%%%%%
%%%%%%%%%%%         specify descriptors of the parameters that will be used to              %%%%%%%%%%%
%%%%%%%%%%%         compose an according file name.                                         %%%%%%%%%%%
%%%%%%%%%%%         A general parameter section provides the opportunity to choose          %%%%%%%%%%%
%%%%%%%%%%%         the number of iterations and the seed for each simulation run.          %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

function parOES(nWorkers)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters that should be explored %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% varNames = {'interval', 'alpha_v'};
% var1 = {10, 20, 50, 100};
% var2 = {0.75};

% varNames = {'regularizer', 'actorLRRange'};
% var1 = {1e-2, 1e-3, 1e-4};
% var2 = {[0.5], [0.5, 0], [1, 0]};

% varDescr = {'regul', 'actorLR'};
% var1Descr = {'1e-2', '5e-2', '1e-3', '5e-3', '1e-4'};
% var2Descr = {'0.5', '0.5To0', '1To0'};

% varNames = {'gamma', 'interval'};
% var1 = {0.1, 0.3, 0.9};
% var2 = {10, 50, 100};

varNames = {'criticLRRange', 'actorLRRange'};
var1 = {[1, 1], [1, 0], [0.75, 0.75], [0.75, 0], [0.5, 0.5], [0.5, 0], [0.25, 0.25], [0.25, 0]};
var2 = {[1, 1], [1, 0], [0.75, 0.75], [0.75, 0], [0.5, 0.5], [0.5, 0], [0.25, 0.25], [0.25, 0]};

%% descriptive parameter names used in folder name
% varDescr = {'cDiscout', 'interval'};
% var1Descr = {'01', '03', '09'};
% var2Descr = {'10', '50', '100'};

varDescr = {'CriticLR', 'ActorLR'};
var1Descr = {'[1,1]', '[1,0]', '[0.75,0.75]', '[0.75,0]', '[0.5,0.5]', '[0.5,0]', '[0.25,0.25]', '[0.25,0]'};
var2Descr = {'[1,1]', '[1,0]', '[0.75,0.75]', '[0.75,0]', '[0.5,0.5]', '[0.5,0]', '[0.25,0.25]', '[0.25,0]'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% general parameter section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nIters = 2000000;               % number of iterations
rSeed = 1;                      % random seed
folderName = 'CriticLR vs ActorLR';
% folderName = 'Regularizer vs Actor Learning Rate' % no ';' intended.
% folderName = 'Discount Factor vs Interval' % no ';' intended.
%TODO enable definition of other parameters that are not changed.
% standardParams = {'textureFile', textureFiles, 'trainTime', trainTime, 'testAt', testAt, 'sparseCodingType', sparseCodingType};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% staring here, the rest is done automatically and should - in gerneral - not be altered  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myCluster = parcluster('local');
myCluster.NumWorkers = nWorkers;
% saveProfile(myCluster);

% create a parallel worker pool if there is none
if (isempty(gcp('nocreate')))
    mPool = parpool(nWorkers);
end

nParams = length(var1) * length(var2);
paramValues = cell(nParams, 2); % column vector, each row contains a set of params for one run

iter = 1;
for i = 1 : length(var1)
    for j = 1 : length(var2)
        paramValues(iter, :) = {cell2mat(var1(i)), cell2mat(var2(j))};
        paramStrings(iter, :) = {cell2mat(var1Descr(i)), cell2mat(var2Descr(j))};
        iter = iter + 1;
    end
end

%% small test
% results = zeros(2, nParams);
%
% parfor (ind = 1 : nParams)
%   results(:, ind) = params(:, ind)
% end
%
% results

%% if no other descriptors are given, take original values
if isempty(varDescr)
    varDescr = varNames;
end

if isempty(var1Descr)
    var1Descr = var1;
end

if isempty(var2Descr)
    var2Descr = var2;
end

% textureFiles = {'Textures_mcgillManMade40.mat', 'Textures_mcgillManMade100.mat'};
% textureFiles = {'mcGillTest2.mat', 'mcGillTest1.mat'}; % test files containing less images
% simulator = prepareSimulator(textureFiles); % at first, we try to use a shared simulator for all threads

%% main loop
parfor ind = 1 : nParams
    sprintf('%s=[%4.2f,%4.2f], %s=[%4.2f,%4.2f]', varNames{1}, paramValues{ind, 1}(1), paramValues{ind, 1}(2), varNames{2}, paramValues{ind, 2}(1), paramValues{ind, 2}(2))

    % OES2Muscles(nIters, rSeed, 1, paramValues(ind), sprintf('%s%s_%s%s', varDescr{1}, paramStrings{ind, 1}, varDescr{2}, paramStrings{ind, 2})); % sprintf('varDec%g--%g', paramValues{ind})

    % OES2Muscles(nIters, rSeed, 1, ...
    %             {varNames{1}, paramValues{ind, 1}, varNames{2}, paramValues{ind, 2}}, ...
    %             sprintf('cluster_%s_%4.2f_%s_%.0f', varDescr{1}, paramValues{ind, 1}, varDescr{2}, paramValues{ind, 2}));

    OES2Muscles(nIters, rSeed, 1, ...
                {varNames{1}, paramValues{ind, 1}, varNames{2}, paramValues{ind, 2}}, ...
                folderName, sprintf('cluster_%s_%s_%s_%s', varDescr{1}, paramStrings{ind, 1}, varDescr{2}, paramStrings{ind, 2}));
end

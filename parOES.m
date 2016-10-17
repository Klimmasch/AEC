%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%																					%%%%%%%%%%%
%%%%%%%%%%%					parallel Matlab simulation with Open Eye Simulator 				%%%%%%%%%%%
%%%%%%%%%%%							Author: Lukas Klimmasch 								%%%%%%%%%%%
%%%%%%%%%%%																					%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Usage:	This scprit creates a pool of parallel workers and distributes			%%%%%%%%%%%
%%%%%%%%%%% 		a number of jobs with different parameter combinations provided			%%%%%%%%%%%
%%%%%%%%%%%			in the first part of the file.											%%%%%%%%%%%
%%%%%%%%%%%			It is possible to start a series of experiments with all 				%%%%%%%%%%%
%%%%%%%%%%%			combinations of two variables that can take different ranges.			%%%%%%%%%%%
%%%%%%%%%%%			To faciliate readability of resulting file names, one can 				%%%%%%%%%%%
%%%%%%%%%%%			specify descriptors of the parameters that will be used to 				%%%%%%%%%%%
%%%%%%%%%%%			compose an according file name.											%%%%%%%%%%%
%%%%%%%%%%%			A general parameter section provides the opportunity to choose			%%%%%%%%%%%
%%%%%%%%%%%			the number of iterations and the seed for each simulation run.			%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters that should be explored %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% varNames = {'interval', 'alpha_v'};
% var1 = {10, 20, 50, 100};
% var2 = {0.75};

% varNames = {'regularizer', 'actorLRRange'};
% var1 = {1e-2, 1e-3, 1e-4};
% var2 = {[0.5], [0.5, 0], [1, 0]};

% varDescr = {'regul', 'actorLR'};
% var1Descr = {'1e-2', '1e-3', '1e-4'};
% var2Descr = {'0.5', '0.5To0', '1To0'};

varNames = {'gamma', 'interval'};
var1 = {0.1, 0.3, 0.9};
var2 = {10, 50, 100};

varDescr = {'cDiscout', 'interval'};
var1Descr = {'01', '03', '09'};
var2Descr = {'10', '50', '100'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% general parameter section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nIters = 2000000;   			% number of iterations
rSeed = 1;                      % random seed
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% staring here, the rest is done automatically and should - in gerneral - not be altered 	%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(gcp('nocreate'))
    mPool = parpool(); % creates a parallel pool with a default number of workers that chould be specified by the cluster this script in running on
end

nParams = length(var1) * length(var2); 
PARAMS = cell(nParams, 2);				%column vector, each row contains a set of params for one run

dummy = 1;
for i = 1 : length(var1)
	for j = 1 : length(var2)
		PARAMS(dummy, :) = {cell2mat(var1(i)), cell2mat(var2(j))};
		dummy = dummy + 1;
	end
end

dummy = 1;
for i = 1 : length(var1)
	for j = 1 : length(var2)
		PARAMSDescr(dummy, :) = {cell2mat(var1Descr(i)), cell2mat(var2Descr(j))};
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
% additional params, that should be set in the config
% paramIdentifiers = ['alpha_v'];
% params = [0.75];

%% main loop, TODO: test if renderers are not interfering (anaglyphs ...)
parfor ind = 1 : nParams
    OES2Muscles(nIters, rSeed, PARAMS(ind), sprintf('%s%s_%s%s', varDescr{1}, PARAMSDescr{ind, 1}, varDescr{2}, PARAMSDescr{ind, 2})); % sprintf('varDec%g--%g', PARAMS{ind})
end
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

% varNames = {'gamma', 'interval'};
% var1 = {0.1, 0.3, 0.9};
% var2 = {10, 50, 100};
% varDescr = {'cDiscout', 'interval'};
% numberFormatVar1 = '%1.1f';
% numberFormatVar2 = '%d';
% experimentDirName = 'Discount Factor vs Interval' % no ';' intended.

% varNames = {'gamma', 'metCostRange'};
% var1 = {0.1, 0.3, 0.6, 0.9};
% var2 = {[0], [0.01], [0.025], [0.05]};
% varDescr = {'gamma', 'metCost'};
% numberFormatVar1 = '%1.1f';
% numberFormatVar2 = '[%1.3f]';
% experimentDirName = 'GammaVsMetCostsBias0.02'

% varNames = {'gamma', 'metCostRange'};
% var1 = {0.01, 0.07, 0.1, 0.2};
% var2 = {[0], [0.0017], [0.0033], [0.005]};
% varDescr = {'gamma', 'metCost'};
% numberFormatVar1 = '%1.2f';
% numberFormatVar2 = '[%1.4f]';
% experimentDirName = 'GammaVsMetCostsBiasFineGrain'

% varNames = {'actorLRRange', 'varianceRange'};
% var1 = {[1,1], [0.5,0.5], [1, 0], [0.5, 0]};
% var2 = {[5e-5,5e-5], [5e-5,0], [1e-5,1e-5], [1e-5,0]};
% varDescr = {'actorLRRange', 'varianceRange'};
% numberFormatVar1 = '[%1.1f-%1.1f]';
% numberFormatVar2 = '[%1.0e-%1.0e]';
% experimentDirName = 'actorLR_vs_varRange_norm_feat_mio' % no ';' intended.

% varNames = {'regularizer', 'actorLRRange'};
% var1 = {1e-4, 1e-5, 1e-6};    % {1e-2, 1e-3, 1e-4};
% var2 = {[1,1], [0.5,0.5], [1, 0], [0.5, 0]};
% varDescr = {'regul', 'actorLR'};
% numberFormatVar1 = '%1.0e';
% numberFormatVar2 = '[%1.2f-%1.2f]';

% experimentDirName = 'actorLR_vs_regul_norm_feat' % no ';' intended.

% varNames = {'criticLRRange', 'actorLRRange'};
% var1 = {[1, 1], [1, 0], [0.75, 0.75], [0.75, 0], [0.5, 0.5], [0.5, 0], [0.25, 0.25], [0.25, 0]};
% var2 = {[1, 1], [1, 0], [0.75, 0.75], [0.75, 0], [0.5, 0.5], [0.5, 0], [0.25, 0.25], [0.25, 0]};

% var1Descr = {'[1,1]', '[1,0]', '[0.75,0.75]', '[0.75,0]', '[0.5,0.5]', '[0.5,0]', '[0.25,0.25]', '[0.25,0]'};
% var2Descr = {'[1,1]', '[1,0]', '[0.75,0.75]', '[0.75,0]', '[0.5,0.5]', '[0.5,0]', '[0.25,0.25]', '[0.25,0]'};

% varDescr = {'CriticLR', 'ActorLR'};
% numberFormatVar1 = '[%1.2f-%1.2f]';
% numberFormatVar2 = '[%1.2f-%1.2f]';

% experimentDirName = 'CriticLR vs ActorLR';

% varNames = {'criticLRRange', 'varianceRange'};
% var1 = {[1, 0], [0.75, 0.75]};
% var2 = {[1e-4, 0], [5e-5, 0], [1e-5, 0], [5e-6, 0], [1e-4, 5e-6], [5e-5, 5e-6], [1e-5, 5e-6], [5e-6, 5e-6], [5e-5, 5e-5], [1e-5, 1e-5], [5e-6, 5e-6]};

% varDescr = {'critic', 'varRange'};
% numberFormatVar1 = '[%1.2f-%1.2f]';
% numberFormatVar2 = '[%1.0e-%1.0e]';

% experimentDirName = 'varDec_new'

% varNames = {'actorLRRange', 'varianceRange'};
% var1 = {[1, 1], [0.5, 0.5], [1, 0], [0.5, 0]};
% var2 = {[1e-4, 1e-4], [1e-5, 1e-5], [1e-4, 1e-6], [1e-5, 1e-6], [1e-4, 0], [1e-5, 0]};

% varDescr = {'actor', 'varRange'};
% numberFormatVar1 = '[%1.1f-%1.1f]';
% numberFormatVar2 = '[%1.0e-%1.0e]';

% experimentDirName = 'steplength_actorVsVariance_reg1e-5'

% varNames = {'actorLRRange', 'varianceRange'};
% var1 = {[1, 1], [0.5, 0.5], [1, 0], [0.5, 0]};
% % var2 = {[5e-5, 5e-5], [5e-5, 0],[1e-5, 1e-5], [1e-5, 0], [4e-5, 4e-5], [4e-5, 0],[3e-5, 3e-5], [3e-5, 0], [2e-5, 2e-5], [2e-5, 0]};
% var2 = {[4e-5, 4e-5], [4e-5, 0],[3e-5, 3e-5], [3e-5, 0], [2e-5, 2e-5], [2e-5, 0], [0.5e-5, 0.5e-5], [0.5e-5, 0]};
% varDescr = {'actor', 'varRange'};
% numberFormatVar1 = '[%1.1f-%1.1f]';
% numberFormatVar2 = '[%1.0e-%1.0e]';
%
% experimentDirName = 'steplength_actorVsVariance_reg1e-5_1mio'

% varNames = {'criticLRRange', 'gamma'};
% var1 = {[1, 1], [0.5, 0.5], [1, 0], [0.5, 0]};
% var2 = {0.1, 0.3, 0.6, 0.9};
% varDescr = {'actor', 'varRange'};
% numberFormatVar1 = '[%1.1f-%1.1f]';
% numberFormatVar2 = '[%1.1f]';

% experimentDirName = 'wobbling_gammaVsCriticLR_1mio'

% varNames = {'gamma', 'metCostRange'};
% var1 = {0.1, 0.3, 0.6, 0.9};
% % var2 = {[0], [0.01], [0.025], [0.05]};
% var2 = {[0], [0.01], [0.05], [0.1]};
% varDescr = {'gamma', 'metCost'};
% numberFormatVar1 = '%1.1f';
% numberFormatVar2 = '[%1.2f]';
% experimentDirName = 'GammaVsMetCosts_0,5mio'
% experimentDirName = 'KillTheBugs'


% varNames = {'criticLRRange', 'metCostRange'};
% var1 = {[1, 1], [1, 0], [0.5, 0.5], [0.5, 0]};
% % var2 = {[0], [0.01], [0.025], [0.05]};
% var2 = {[0.025, 0.025], [0.025, 0], [0.05, 0.05], [0.05, 0], [0, 0]};
% varDescr = {'critLR', 'metCost'};
% numberFormatVar1 = '[%1.1f-%1.1f]';
% numberFormatVar2 = '[%1.3f-%1.3f]';
% % 12 runs: 1992 mb per cpu on autumnchat
% experimentDirName = 'CritLRVsMetCosts_1mio'

% var1Descr = {'[1,1]', '[1,0]', '[0.75,0.75]', '[0.75,0]', '[0.5,0.5]', '[0.5,0]', '[0.25,0.25]', '[0.25,0]'};
% var2Descr = {'[1,1]', '[1,0]', '[0.75,0.75]', '[0.75,0]', '[0.5,0.5]', '[0.5,0]', '[0.25,0.25]', '[0.25,0]'};
% varDescr = {'CriticLR', 'ActorLR'};
% numberFormatVar1 = '[%1.2f-%1.2f]';
% numberFormatVar2 = '[%1.2f-%1.2f]';
% experimentDirName = 'CriticLR vs ActorLR';

% varNames = {'regularizer', 'desiredStdZT'};
% var1 = {5e-5, 1e-5, 5e-6, 1e-6};
% var2 = {0.005, 0.0075, 0.01, 0.015, 0.02, 0.025};
% var1Descr = {'5e-5', '1e-5', '5e-6', '1e-6'};
% var2Descr = {'0.005', '0.0075', '0.01', '0.015', '0.02', '0.025'};
% varDescr = {'regularizer', 'desiredStdZT'};
% experimentDirName = 'Regularizer_vs_desiredStdZT_short';

% varNames = {'lambdaMuscleFB', 'desiredStdZT'};
% var1 = {0, 0.5, 1, 2};
% var2 = {0.005, 0.0075, 0.01, 0.015, 0.02};
% var1Descr = {'0', '0.5', '1', '2'};
% var2Descr = {'0.005', '0.0075', '0.01', '0.015', '0.02'};
% varDescr = {'lambdaMuscleFB', 'desiredStdZT'};
% experimentDirName = 'lambdaMuscleFB_vs_desiredStdZT_seed2';

% varNames = {'gamma', 'metCostRange'};
% var1 = {0.01, 0.05, 0.1};
% % var2 = {[0], [0.01], [0.03], [0.05]};
% var2 = {[0], [0.05], [0.075], [0.1]};
% varDescr = {'gamma', 'metCost'};
% numberFormatVar1 = '%1.2f';
% numberFormatVar2 = '[%1.3f]';
% experimentDirName = 'GammaVsMetCosts_FineGrainNoBias'

% varNames = {'filterLeft', 'filterLeftProb'};
% % var1 = {8, 9, 10, 11, 12, 13};
% var1 = {13, 12, 11, 10, 9, 8};
% % var2 = {0.1, 0.25, 0.5, 0.75, 0.9, 1};
% var2 = {1};
% varDescr = {'filtL', 'filtLProb'};
% numberFormatVar1 = '%d';
% numberFormatVar2 = '%d';
% experimentDirName = 'uniformBlurrInput'

varNames = {'filterLeft', 'filterLeftProb'};
% var1 = {26, 27, 28};
var1 = {34};% {45, 46};%, 30, 31, 32, 33, 34, 35, 36, 37};
% var1 = fliplr(var1);
% var2 = {0.1, 0.25, 0.5, 0.75, 0.9, 1};
var2 = {1};
varDescr = {'filtBoth', 'prob'};
numberFormatVar1 = '%d';
numberFormatVar2 = '%d';
% experimentDirName = 'edgeDeprivNew'
experimentDirName = 'explFilterSizes'
globalParams = {};
globalName = 'fsize6std_longer_'

% varNames = {'filterRight', 'filterRightProb'};
% % var1 = {8, 9, 10, 11, 12, 13};
% var1 = {2};
% var2 = {0, 1};
% varDescr = {'MDfilterR', 'filtProb'};
% numberFormatVar1 = '%d';
% numberFormatVar2 = '%d';
% experimentDirName = 'disp+verg_influence'
%
% globalParams = {'objDistMin', 3, 'objDistMax', 3, 'fixDistMin', 3, 'fixDistMax', 3, ...
%  'initMethod', 1, 'actorLRRange', [0], 'varianceRange', [0], 'weight_range', [0,0,0]};
% globalName = 'OD3_initMet1_ALR0_Var0_weights0_'

% varNames = {'initMethod', 'lapSig'};
% % var1 = {8, 9, 10, 11, 12, 13};
% var1 = {4};
% % var2 = {0, 0.22, 0.5, 1, 2};
% var2 = {200};%, 200};
% varDescr = {'initMethod', 'lapSigma'};
% numberFormatVar1 = '%d';
% numberFormatVar2 = '%0.2f';
% experimentDirName = 'laplacianPolicy'
%
% globalParams = {'actorLRRange', [0], 'varianceRange', [0], 'weight_range', [0,0,0], ...
%                 'nBasis', [400], 'nBasisUsed', [10], 'basisSize', [128], ...
%                 'sc_eta', [0.2], 'temperature', [0.01], 'dsRatio', [1], ...
%                 'pxFieldOfViewOrig', [40], 'pxFieldOfView', [40], 'stride', [4]};
% globalName = 'noLearn_fineScOnly_'

% varNames = {'initMethod', 'lapSig'};
% var1 = {4};%, 4};
% var2 = {0};
% varDescr = {'initMethod', 'lapSigma'};
% numberFormatVar1 = '%d';
% numberFormatVar2 = '%d';
% experimentDirName = 'vergenceInfluence'
% globalParams = {'actorLRRange', [0], 'varianceRange', [0], 'weight_range', [0,0,0], 'filterLeft', 29};
% globalName = 'noLearning_'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% general parameter section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nIters = 2000000;               % number of iterations
rSeeds = [1, 2, 3]     % random seeds
% rSeeds = [2];
% rSeeds = [1, 5]       % random seeds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% staring here, the rest is done automatically and should - in gerneral - not be altered  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nParams = length(var1) * length(var2);
paramValues = cell(nParams, 2); % column vector, each row contains a set of params for one run

iter = 1;
for i = 1 : length(var1)
    for j = 1 : length(var2)
        paramValues(iter, :) = {cell2mat(var1(i)), cell2mat(var2(j))};
        % paramStrings(iter, :) = {cell2mat(var1Descr(i)), cell2mat(var2Descr(j))}; % temporary solution for not renaming all folders ...
        paramStrings(iter, :) = {num2str(var1{i}, numberFormatVar1), num2str(var2{j}, numberFormatVar2)};
        iter = iter + 1;
    end
end

%% if no other descriptors are given, take original values
if isempty(varDescr)
    varDescr = varNames;
end

myCluster = parcluster('local');
myCluster.NumWorkers = nWorkers;
% saveProfile(myCluster);

% create a parallel worker pool if there is none
if (isempty(gcp('nocreate')))
    mPool = parpool(nWorkers);
end

% textureFiles = {'Textures_mcgillManMade40.mat', 'Textures_mcgillManMade100.mat'};
% textureFiles = {'mcGillTest2.mat', 'mcGillTest1.mat'}; % test files containing less images
% simulator = prepareSimulator(textureFiles); % at first, we try to use a shared simulator for all threads

%% main loop
if length(rSeeds) == 1
    parfor ind = 1 : nParams
    % for ind = 1 : nParams
        % sprintf('%s=[%4.2f,%4.2f], %s=[%4.2f,%4.2f]', varNames{1}, paramValues{ind, 1}(1), paramValues{ind, 1}(2), varNames{2}, paramValues{ind, 2}(1), paramValues{ind, 2}(2))

        sprintf('%s_%s_%s_%s', varDescr{1}, paramStrings{ind, 1}, varDescr{2}, paramStrings{ind, 2})
        OES2Muscles(nIters, rSeeds(1), 1, ...
                    {varNames{1}, paramValues{ind, 1}, varNames{2}, paramValues{ind, 2}, globalParams{:}}, ...
                    experimentDirName, sprintf('%s%s_%s_%s_%s_seed%d', globalName, varDescr{1}, paramStrings{ind, 1}, varDescr{2}, paramStrings{ind, 2}, rSeeds(1)));
    end
else
    parfor ind1 = 1 : length(rSeeds)
        for ind2 = 1 : nParams
        % for ind2 = 1 : length(rSeeds)
            sprintf('%s_%s_%s_%s', varDescr{1}, paramStrings{ind2, 1}, varDescr{2}, paramStrings{ind2, 2})
            OES2Muscles(nIters, rSeeds(ind1), 1, ...
                        {varNames{1}, paramValues{ind2, 1}, varNames{2}, paramValues{ind2, 2}, globalParams{:}}, ...
                        experimentDirName, sprintf('%s%s_%s_%s_%s_seed%d', globalName, varDescr{1}, paramStrings{ind2, 1}, varDescr{2}, paramStrings{ind2, 2}, rSeeds(ind1)));
        end
    end
end

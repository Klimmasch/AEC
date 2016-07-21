% This script prepares a training and a testing set from all images present in the passed directory.
% The generated train and test matricies are saved into the config folder.
%
% @param stimFolder         image folder
% @param nTestFiles         either #testStimuli or proportion of testStimuli,
%                           i.e. nTestFiles elem. [0, 2, #allStimuli - 1] and [0, 1[
% @param stimulusSetName    name of stimulus set
function generateTrainTestSets(stimFolder, nTestFiles, stimulusSetName)

files = dir(sprintf('%s/*.bmp', stimFolder));
nFiles = length(files);

% savety check
if ((nTestFiles == 1) || (nTestFiles >= nFiles))
    sprintf('Error: invalid value for nTestFiles!\nnTestFiles elem. [0, 2, #allStimuli - 1] and [0, 1[')
    return;
end

% special case: use all stimuli in a single set
if (nTestFiles <= 0)
    stimulusSet = cell(nFiles, 1);
    for i = 1 : nFiles
        stimulusSet{i} = sprintf('../aec/Textures/%s/%s', stimulusSetName, files(i).name);
    end
    save(sprintf('../config/Textures_%sAll.mat', stimulusSetName), 'stimulusSet');
    sprintf('Added %d files to the stimulus set', nFiles)
    return;
end

if (nTestFiles < 1)
    nTestFiles = floor(nFiles * nTestFiles);
end

shuffledFiles = files(randperm(nFiles));

% add nTestFiles images to the testing set
testSet = cell(nTestFiles, 1);
for i = 1 : nTestFiles
    testSet{i} = sprintf('../aec/Textures/%s/%s', stimulusSetName, shuffledFiles(i).name);
end
save(sprintf('../config/Textures_%sTest.mat', stimulusSetName), 'testSet');

% add remaining stimuli to training set
trainSet = cell(nFiles - nTestFiles, 1);
j = 1;
for i = nTestFiles + 1 : nFiles
    trainSet{j} = sprintf('../aec/Textures/%s/%s', stimulusSetName, shuffledFiles(i).name);
    j = j + 1;
end
save(sprintf('../config/Textures_%sTrain.mat', stimulusSetName), 'trainSet');

sprintf('Added %d files to the testing set and %d to the trainig set', nTestFiles, nFiles - nTestFiles)

%for savety:
% for i = 1:length(testSet)
%     for j = 1:length(texture)
%         if strcmp(testSet{i}, texture{j})
%             sprintf('double occurence: %s', texture(j))
%         end
%     end
% end
% display('done')

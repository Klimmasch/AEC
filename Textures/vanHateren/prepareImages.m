% This script prepares a training and a testing set from all images that
% are saves in the van Hateren folder and saves it in the accoding files in
% the config folder

files = dir('*.bmp');
nFiles = length(files);
nTestFiles = 50;

shuffledFiles = files(randperm(nFiles));
texture = {};

% add nTestFiles images to the testing set
for i = 1:nTestFiles
    texture{i} = ['../aec/Textures/vanHateren/' shuffledFiles(i).name];
end
texture = texture';
save('../../config/Textures_vanHaterenTest.mat', 'texture')

texture = {};
% add all others to training set
for i = 1:(nFiles - nTestFiles)
    texture{i} = ['../aec/Textures/vanHateren/' shuffledFiles(nTestFiles + i).name];
end
texture = texture';
save('../../config/Textures_vanHaterenTrain.mat', 'texture')

sprintf('Added %d files to the testing Set and %d to the trainig set', nTestFiles, nFiles - nTestFiles) 
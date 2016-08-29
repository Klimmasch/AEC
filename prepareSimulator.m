%%% Script for instantiation of a simulator that can be used for testing etc.
%%% Load will be the same images that are used for training and testing. 
function simulator = prepareSimulator()

    textureFiles = {'Textures_mcgillManMade40.mat', 'Textures_mcgillManMade100.mat'};
%     textureFiles = {'mcGillTest2.mat', 'mcGillTest1.mat'}; % test files containing less images
    
    % simulator = OpenEyeSim('create');
    simulator = OpenEyeSimV5('create');

    simulator.initRenderer();
%     simulator.reinitRenderer();
    
    % load Textures
    nTextures = 0;
    tmpTexInd = 1;
    for i = 1:length(textureFiles)
        texture = load(sprintf('config/%s', textureFiles{i}));
        texture = texture.texture;
        nTextures = nTextures + length(texture);
        if i == 1
            nTestTextures = nTextures;      % save number of test textures for proper indexing later
        end
        textureInd = 1;
        for j = tmpTexInd : nTextures
            simulator.add_texture(j, texture{textureInd});
            textureInd = textureInd + 1;
        end
        tmpTexInd = tmpTexInd + nTextures;
    end
    sprintf('%d textures were added to the buffer from the training and testing files', nTextures)
    sprintf('The first %d textures are the testing set.', nTestTextures)
    
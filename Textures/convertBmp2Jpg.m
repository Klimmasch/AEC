
folder = '../config/'
matFiles = dir(strcat(folder, '*mcgill*(jpg).mat')) % I made safety copies of all files and appended '(jpg)' 
nMats = length(matFiles);

for mFile = 1:nMats
    fName = strcat(folder, matFiles(mFile).name)
    file = load(fName)
    nImages = length(file.texture)
    
    for texture = 1:nImages
        file.texture{texture} = regexprep(file.texture{texture}, '.bmp', '.jpg')
    end
    texture = file.texture;
    save(fName, 'texture')
end
% mex command to build the C files

% dynamic linking
% mexFile=' OpenEyeSim_.cc'; % stable renderer
mexFile=' OpenEyeSimV5_.cc'; % experimental version
if mislocked(mexFile)
    munlock(mexFile)
end

% get paths
% ATTENTION: this enviroment variable needs to be defined by
% export OPENSIM_INSTALL_DIR='/directory/to/OpenSimInstall'
oSimDir = getenv('OPENSIM_INSTALL_DIR');
oSimLibDir = sprintf(' \''%s/lib/lib', oSimDir);

oSimIncludeDir1 = [oSimDir '/sdk/include'];
oSimIncludeDir2 = [oSimDir '/sdk/include/SimTK/include'];
mexIncludeDir1 = ['include/mexplus'];
mexIncludeDir2 = ['include'];
include = sprintf(' -I\''%s\'' -I\''%s\'' -I\''%s\'' -I\''%s\''', oSimIncludeDir1, oSimIncludeDir2, mexIncludeDir1, mexIncludeDir2);

%set shared libs that are used
libs=[oSimLibDir 'osimSimulation.so'''...
    oSimLibDir 'osimCommon.so'''...
    oSimLibDir 'osimActuators.so'''...
    oSimLibDir 'osimTools.so'''...
    oSimLibDir 'osimAnalyses.so'''...
    oSimLibDir 'SimTKcommon.so'''...
    oSimLibDir 'SimTKsimbody.so'''...
    oSimLibDir 'SOIL.so'''...
    oSimLibDir 'SimTKmath.so'''];
    %oSimLibDir 'SimTKlapack.so'''];
disp(oSimLibDir)
disp(libs)
eval(['mex -v -lGL -lglut -lGLU -DGL_GLEXT_PROTOTYPES -L/usr/lib/ -llapack -lblas' mexFile libs include]);
sprintf('Compilation done!')

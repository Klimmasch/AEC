Active Efficient Coding - Setup instructions
============================================

Before you are able to run any new experiments with the AEC environment,
you need to compile the OpenEyeSimRendererV3, which is used to generate
visual scenes in a virtual environment. For that you execute:

$ setupOESRendererV3.sh

If you get a success message, everything is set up properly.

In case the OpenEyeSim_.cc is modified, the OESRendererV3 needs to be
recompiled again. For that you execute:

$ recompileOESRendererV3.sh

In case you want to use the depricated checkEnvironment renderer, you
need to copy all library objects back to the repositories root directory:

$ cp OpenSimInstall/lib/*.so* .


#!/bin/bash
#
# Script for setting up and compiling OpenEyeSim Renderer V3
#
echo ""
echo "OpenEyeSim Renderer V3 Setup"
echo "============================"

# Export neccessary environment variables
echo -n "Exporting neccessary environment variables..."
TMP="$(echo $(pwd)/OpenSimInstall)"
export OPENSIM_INSTALL_DIR=$TMP
export LD_LIBRARY_PATH=$TMP/lib
echo "export LD_LIBRARY_PATH=$TMP/lib" >> ~/.bashrc
echo "done."

# Compile MEX object
echo -n "Compiling MEX object..."
LOG="./err.log"
matlab -nodisplay -r "makeOpenEyeSim; quit" | tee -a $LOG
if [ $? -eq 0 ]
then
    exit 0
else
    echo "Error occured! See $LOG for details."
    exit 1
fi

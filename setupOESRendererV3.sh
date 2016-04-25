#!/bin/bash
#
# Script for setting up and compiling OpenEyeSim Renderer V3
#
echo ""
echo "OpenEyeSim Renderer V3 Setup"
echo "============================"
echo -n "Exporting neccessary environment variables..."

# Export neccessary environment variables
TMP="$(echo $(pwd)/OpenSimInstall)"
echo "export OPENSIM_INSTALL_DIR='$TMP'" >> ~/.bashrc
echo "export LD_LIBRARY_PATH='$TMP/lib'" >> ~/.bashrc
echo "done."

# Compile MEX object
LOG="./err.log"
echo -n "Compiling MEX object..."
matlab -nodisplay -r "makeOpenEyeSim; quit" >> $LOG 2>&1
if [ $? -eq 0 ]
then
    echo "done."
else
    echo "Error occured! See $LOG for details."
    exit 1
fi

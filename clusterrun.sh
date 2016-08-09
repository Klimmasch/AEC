#!/bin/bash

# login
srun --pty -p sleuths --exclusive -w marple bash

# display generation
X :0 2>/dev/null &
export DISPLAY=:0
unset SESSION_MANAGER

# neccessary  path variable
REPO_DIR=$(echo $(pwd))
export LD_LIBRARY_PATH="$REPO_DIR/OpenSimInstall/lib"

# start experiment
matlab -nodisplay -r "$1"

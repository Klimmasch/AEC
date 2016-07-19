Active Efficient Coding - Setup instructions
============================================

Before you are able to run any new experiments with the AEC environment, you
need to compile the OpenEyeSimRenderer, which is used to generate visual
scenes in a virtual environment. For that you execute:

$ setupOESRenderer.sh

If you get a success message, everything is set up properly. If you get error
messages, p.e. if makeOpenEyeSim couldn't be located, try executing
makeOpenEyeSim within the Matlab environment (within the GUI). The script was
tested within GNOME, but within other desktop environments the afore mentioned
errors can occur.

In case the OpenEyeSim_.cc is modified, the OESRenderer needs to be recompiled
again. For that you execute once more:

$ setupOESRenderer.sh

In case you want to use the depricated checkEnvironment renderer, you need to
copy all library objects back to the repositories root directory:

$ cp OpenSimInstall/lib/*.so* .


Experiment execution instructions
=================================

There are currently three experimental situations supported:
1) Discrete action space reinforcement learner models

2) Continuous action space reinforcement learner models with a medial rectus eye
   muscle, further single/1-muscle models

3) Continuous action space reinforcement learner models with medial rectus and
   lateral rectus eye muscles, further 2-muscles models

Depending on the situtation, you need to execute a different experimental/main
script, which you can of course use as an inspiration to implement your own
experimental/main script.

Situation) Main script
----------------------
1) OESDiscrete.m
2) OES1Muscle.m
3) OES2Muscles.m

You either execute one of those directly within Matlab or in case you want to do
it within the shell:

$ matlab -nodisplay -r "OES2Muscles(with, respective, parameter, passing)"


The algorithm works roughly in the following way:

a) OES2Muscles calls config.m
   In config.m you need to define all model parameters beforehand

b) config.m calls Model.m to create a model object instance

c) Model.m calls SparseCoding*.m to create the necessary sparse coder object
   instances

d) Model.m calls ReinforcementLearning*.m to create respective RL object
   instances

   d.2) in case of *-muscles models, ReinforcementLearningCont.m calls
        *Actor*.m and *Critic*.m to create respective actor and critic instances

e) Model.m initiates/preallocates some necessary constants and variables and
   finalizes the model object creation

f) OES*.m script continues by backing up all used scritps, classes and config
   files into ../results/timestamp_user_defined_name_of_experiment/

g) OES*.m script conducts the experiment and saves all resulting figures of the
   training and of course the model file into the directory mentioned in f)

h) OES*.m script executes testModelContinuous.m to conduct the testing procedure
   whereby all figures generated are also saved into the folder mentioned in f)

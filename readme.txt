This simulation was used in the following article: Neymotin SA,
  McDougal RA, Sherif MA, Fall CP, Hines ML, Lytton WW.  Neuronal
  calcium wave propagation varies with changes in endoplasmic
  reticulum parameters: a computer model.  Neural Computation 2015 (in
  press).

The code in this folder generates Fig. 2 (basic calcium wave) and
Fig. 11 (varying AMPA synapse stimulation parameters and view effects
on the calcium wave).

The simulations were tested/developed on LINUX systems, but may run on
Microsoft Windows or Mac OS.

To run the demo, you will need the NEURON simulator (available at
http://www.neuron.yale.edu) compiled with python enabled. To draw the
output you will need to have Matplotlib installed (
http://matplotlib.org/ ).

Instructions to run the model 
- unzip the file
- cd ca1dDemo
- nrnivmodl (compiles the NMODL files)

The nrnivmodl command will produce an architecture-dependent folder
with a script called special.  On 64 bit systems the folder is x86_64.

Note that these simulations will run with NEURON's variable time step
activated, in order to reduce the time it takes to run the simulation.
However, the simulation can still take a long time to run, depending
on your hardware setup. Therefore, the code is setup to save
simulation data to a folder called /data within the simulation
directory.

----------------------------------------------------------------------

Running and plotting the baseline calcium wave figure (Fig. 2):

# run the following code in a terminal from within the directory
  containing the model files:
python
from batch import *
baseRun() # run and save a file of numpy arrays
baseDraw()  # should show a plot of fig. 2 in the article

The figure shows Ca2+ wave propagation with baseline
parameters. Elevated IP3 stimulus placed at mid-dendrite (500 um on
y-axis) after 2 s past start of simulation. The plot on the left
depicts cytosolic [Ca2+] showing a wave of increased
concentration. The plot on the right of ER [Ca2+ ] shows a mirror
image wave of decreased concentration as Ca2+ is released to cytosol.

----------------------------------------------------------------------

Running and plotting Fig. 11 (electrochemical model which shows effect
of AMPA receptor stimulation on release of calcium from ER): # you
will be running 2 simulations with 2 different configuration files
(one for no AMPA receptor stimulation and one for stimulation of AMPA
receptor with 150 inputs). Both configuration files insert a set of
ion channels in the dendritic section.

# this first simulation takes less time to run (~75 seconds on a Xeon
  E5/Core i7 Integrated Memory Controller processor):
python -i cawave.py AMPA0.cfg

# the following simulation took ~54 minutes when run using the same processor:
python -i cawave.py AMPA150.cfg

To plot the output from these simulations, run the following (after
you exit from the previous simulations):
python
from batch import *
execfile('plot_fig11.py')

The figure shows electrical stimulation with increased number of AMPA
activations enhancing Ca2+ waves induced by IP3 (2.5 mM at 7 s). The
left column shows control simulation: Ca2+ wave with no AMPA inputs
prior to the IP3 stimulus. Middle column shows Ca2+ wave with train of
150 AMPA inputs (onset: 3 s; interspike interval: 25 ms) prior to the
IP3 stimulus. The column on the right is a comparison of voltage
(top), ER Ca2+ (middle), and cytosolic Ca2+ (bottom) in control
(black) and simulation with 150 AMPA inputs (red).


NOTE: If you are interested in running the simulation for a shorter
period of time, you can modify the configuration files AMPA0.cfg and
AMPA150.cfg. Change the value of tstop under [run] to the time you're
interested in (in milliseconds). However, the resulting figure will
not be identical to Fig. 11 in the publication.
----------------------------------------------------------------------

For questions/comments email: 
 mohamed dot sherif dot md at gmail dot com
 or
 samn at neurosim dot downstate dot edu
 or
 robert dot mcdougal at yale dot edu

20160915 This updated version from the Lytton lab allows their models
which contain misc.mod and misc.h to compile on the mac.

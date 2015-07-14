#Takes cluster sizes from a trajectory and saves them as a .txt file
import markov as mk

t0=0 #Initial time of interest in ps
tf=500000 #Final time of interest in ps
dt = 5 #Minimum time spacing between .xtc windows in ps (nstxtcout*dt)
xtc = 'md_whole_restart.xtc' #.xtc file name
tpr = 'md_dummy.tpr' #.tpr file name
cutoff=.5 #minimum distance between any two atoms for system to be considered a cluster
ats=29 #atoms per molecule (for parsing the .gro file into separate molecules; all atoms in the gro file need to be part of molecules of the same size for cluster analysis. Currently working on fixing this.)
outputFileName = 'mdrun2_sizes.txt'

#---------------------------------------------------------#

tlist = range(t0,tf+1, dt)

mk.sizerun(tlist, xtc, tpr, cutoff, ats, fnm=outputFileName)

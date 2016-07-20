#Takes cluster sizes from a trajectory and saves them as a .txt file
import imp
mk = imp.load_source('mk','/home/mansbac2/coarsegraining/code/markov/markov.py')

t0=0 #Initial time of interest in ps
tf=TT #Final time of interest in ps
dt = 50 #Minimum time spacing between .xtc windows in ps (nstxtcout*dt)
xtc = 'md_noW.xtc' #.xtc file name
tpr = 'md_dummy.tpr' #.tpr file name
cutoff=.5 #minimum distance between any two atoms for system to be considered a cluster
ats=29 #atoms per molecule (for parsing the .gro file into separate molecules; all atoms in the gro file need to be part of molecules of the same size for cluster analysis. Currently working on fixing this.)
outputFileName = 'mdrun_sizes.dat'

#---------------------------------------------------------#

tlist = range(t0,tf+1, dt)

mk.sizerun(tlist, xtc, tpr, cutoff, ats, fnm=outputFileName)

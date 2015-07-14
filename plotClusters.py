#Plots cluster sizes stored in a .txt file output by savesizes.py.
import markov as mk, numpy as np

t0=0 #Initial time of interest in ps
tf=218 #Final time of interest in ps
trange=2500000 #Time window to use for blocking in ps
dt=10 #Time between Markov steps
dtmin=5 #Minimum time between windows (nstxtcout*dt in the .mdp file)
factor = 1 #Scaling of time
molecules = 343 #Number of molecuels in the simulation

f1='mdrun1_sizes.txt' #Name of output file of savesizes.py
#f2='longrun2_sizes.txt'
#f3='longrun3_sizes.txt'
flist=[f1]
plotTimeScale = 1000 #Scales units of time to ns as is.

#-------------------------------------------------------------#

t0 = t0*factor
tf = tf*factor
trange=trange*factor
dt = dt*factor
dtmin = dtmin*factor

tlist = range(t0,tf+1, dt)

tlist=np.array(tlist)

mk.plotdata(flist, tlist/plotTimeScale, molecules, t0=t0, dt=dt, dtmin=dtmin, tf=tf, merlist=[1,2,3,4], exe=0) #Plots cluster sizes of 1, 2, 3, and 4. Will not crash if no 4-mers were produced during the simulation, but will crash if fewer or more cluster sizes are listed here.
mk.plotdata(flist, tlist/plotTimeScale, molecules, t0=t0, dt=dt, dtmin=dtmin, tf=tf, merlist=[5,6,7,8], exe=0) #Plots cluster sizes of 5, 6, 7, and 8

mk.pqintfull(flist, dt, 4*dt, t0=t0, tf=tf, dtmin=dtmin)

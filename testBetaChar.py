import numpy as np,imp,time#,markov as mk
mk = imp.load_source('mk','/home/rachael/Analysis_and_run_code/coarsegraining/markov/markov.py')
betafname = "betaCounts_short_1.dat"
t0 = 400000
tf = 400000
dt = 500
xtc = '../sizedata/md_noW.xtc'
tpr = '../sizedata/md_dummy.tpr'
cutoff = 0.5
ats = 29
mols = 14
start = time.clock()
tlist = range(t0,tf+1,dt)
tlist = [tf]
betaf = open(betafname,'w')
#write out indices to one file
#write out to second file for each time frame a list of eigenvalues + Rg + asphericity for all clusters (so these lines will change lengths as Nclusts changes)
#write out to third file for each time frame a list of principal eigenvectors for all clusters (lines will also change lengths)
for t in tlist:
	betacount = mk.betaClust(t,xtc,tpr,'tempc.gro',cutoff,ats,True)
	#print inds
	#print eigvals
	#print eigvecs
	#print Rhs
	lineb = ''
	print("len(betacount) = ",len(betacount))
	for i in range(len(betacount)):
		lineb+=str(t)+' '
		for j in range(ats):
			#linds[j+i*ats] = inds[i]
			lineb += str(betacount[i][j])+' '
		lineb+='\n'
	betaf.write(lineb)	
betaf.close()

#eigS =mk.clustMorph(tf,xtc,tpr,'temp.gro',cutoff,ats,boxL)
end = time.clock()
print "time elapsed: ",end-start
#print len(inds)
#print inds

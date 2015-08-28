import numpy as np,imp,time#,markov as mk
mk = imp.load_source('mk','/home/rachael/Analysis_and_run_code/coarsegraining/markov/markov.py')
indsfname = "tempInds_1.dat"
scalarfname = "clustScalars_short_2.dat"
vecfname = "tempVecs_1.dat"
hydrfname = "tempRh_1.dat"
mRfname = "mass_v_Rg_short_2.dat"
t0 = 0
tf = 500000
dt = 500
xtc = 'md_noW_restart.xtc'
tpr = 'md_dummy.tpr'
cutoff = 0.5
ats = 29
mols = 343
start = time.clock()
tlist = range(t0,tf+1,dt)
#tlist = [t0]
indsf = open(indsfname,'w')
scalarf = open(scalarfname,'w')
vecf = open(vecfname,'w')
hydrf = open(hydrfname,'w')
mRf = open(mRfname,'w')
boxL = [33.90642,33.90642,33.90642]
#linds = np.zeros(ats*mols)
#write out indices to one file
#write out to second file for each time frame a list of eigenvalues + Rg + asphericity for all clusters (so these lines will change lengths as Nclusts changes)
#write out to third file for each time frame a list of principal eigenvectors for all clusters (lines will also change lengths)
for t in tlist:
	(inds,eigvals,eigvecs,Rhs,ms,Rgs) = mk.clustMorph(t,xtc,tpr,'tempc.gro',cutoff,ats,True)
	#print inds
	#print eigvals
	#print eigvecs
	#print Rhs
	linei = ''
	lines = ''
	linev = ''
	liner = ''
	linem = ''
	for i in range(len(inds)):
		linei+=str(t)+' '
		for j in range(ats):
			#linds[j+i*ats] = inds[i]
			linei += str(inds[i])+' '
		linei+='\n'
	for i in range(0,len(eigvals),5):
		#for j in range(5):
		lines += str(t)+' '+str(eigvals[i])+' '+ str(eigvals[i+1])+' ' + str(eigvals[i+2])+' '+str(eigvals[i+3])+' '+str(eigvals[i+4])+'\n'
	for i in range(len(ms)):
	
		#linem += str(ms[i])
		#linem += ' '
		#linem += str(eigvals[4*i+3])
		#linem += '\n'
		linem += str(ms[i])+' '+str(Rhs[i])+' '+str(Rgs[i])+'\n'
		#linem += '\n'		
	for i in range(len(Rhs)):
		liner += str(t)+' '+str(Rhs[i])+'\n'
	for i in range(len(eigvecs)):
		#for j in range(3):
		#for k in range(3):
		linev += str(t)+' '+str(eigvecs[i])+'\n'
	indsf.write(linei)
	scalarf.write(lines)
	vecf.write(linev)
	hydrf.write(liner)
	mRf.write(linem)	
hydrf.close()
indsf.close()
scalarf.close()
vecf.close()

#eigS =mk.clustMorph(tf,xtc,tpr,'temp.gro',cutoff,ats,boxL)
end = time.clock()
print "time elapsed: ",end-start
#print len(inds)
#print inds

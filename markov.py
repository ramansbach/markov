import numpy as np, time, os, numpy.linalg as npla, matplotlib.pyplot as plt, scipy.linalg as spla, time, scipy.stats as spst, datatoolbox as dt, string, matplotlib.ticker as mtick
from scipy import weave
from operator import sub, div
params = {'axes.titlesize':70,
		'axes.labelsize':60,
		'text.fontsize':60,
		'font.size':50,
		'lines.markersize':6,
		'lines.linewidth':4,
		'text.usetex':True,
		'xtick.major.pad':7,
		'ytick.major.pad':18}
plt.rcParams.update(params)

#Takes a .gro file and returns a 1d list of positions. If the gro file format changes, this will break
def readGro(fName): 
	with open(fName, 'r') as myF:
		myLns = myF.read().splitlines()

	return np.array([[float(myLns[i][20:].split()[0]), float(myLns[i][20:].split()[1]), float(myLns[i][20:].split()[2])] for i in range(2, len(myLns)-1)]).flatten()

#Gets the positions of all atoms in trajectory trj, run from tpr file tpr, at time t, written to output gro file outGro.
def getPos(t, trj, tpr, outGro):
	os.system('echo 0 | gmx trjconv -f ' + trj + ' -o ' + outGro + ' -b ' + str(t) + ' -e ' + str(t) + ' -s ' + tpr)
	return readGro(outGro)

#Gets the positions of all atoms from a given run from the separated .gro files produced by trjconv.
def getpos2(ind, grobase):
	return readGro(grobase + str(ind) + '.gro')

#Gets the neighbors of atom ind from a list of potential indices potentialInds (each of which are indices of the list of all peptides listed in peplist). Being neighbors is defined as two peptides having any two atoms separated by less than cutoff. ats is the number of atoms per peptide in peplist. Summary: ind: index to check, cutoff: distance that defines neighbors (as separation of any two atoms), peplist: list of all atoms, potentialInds: potential neighbors, ats: atoms per peptide.
def getNeigh(ind, cutoff, peplist, potentialInds, ats):
	ret = []

	cutsq = cutoff**2
	support = '#include <math.h>'

	code = """

	int i, j;
	return_val = 0;
	for(i=0; i<Npep1[0]/3; i++){
		for(j=0; j<Npep2[0]/3; j++){
			if ((pow(pep1[3*i]-pep2[3*j],2) + pow(pep1[3*i+1]-pep2[3*j+1],2) + pow(pep1[3*i+2]-pep2[3*j+2],2)) < cutsq){
				return_val = 1;
				break;
			}
		}
		if(return_val == 1)
			break;
	}
			"""
	pep1 = peplist[ind*3*ats:(ind+1)*3*ats] #Assumes all peptides have ats atoms in them. The 3 is for 3 coords in 3 dimensions.
	for i in range(len(potentialInds)):
		pep2 = peplist[potentialInds[i]*3*ats:(potentialInds[i]+1)*3*ats]
		test = weave.inline(code,['pep1', 'pep2', 'cutsq'], support_code = support, libraries = ['m'])
		if test == 1:
			ret.append(potentialInds[i])
	
	return ret

#Returns an array of arrays. Each inner array is a list of peptide indices that are in the same cluster. A cluster is defined as the largest list of molecules for which each molecule in the list is neighbors either directly or indirectly (neighbors of neighbors of neighbors etc...) neighbors with each other molecule. ind is the atom to check the cluster of, cutoff is the minimum distance that defines neighbors, peplist is a list of all atoms in the simulation, potentialInds is the list of indices that could possibly be neighbors with ind, ats is the number of atoms per peptide (this is important for dividing up pepList and must be constant for all peptides in pepList), and printMe is a boolean that will cause the immediate neighbors of ind to be printed if it is true (more for debuging and checking things).
def getClust(ind, cutoff, pepList, potentialInds, ats, printMe):
	neighInds = getNeigh(ind, cutoff, pepList, potentialInds, ats)
	if printMe:
		print("Neighbors of " + str(ind) + " are found to have indices of: ")
		print(neighInds)

	for neighInd in neighInds:
		potentialInds.remove(neighInd)
	for neighInd in neighInds:
		vals = getClust(neighInd, cutoff, pepList, potentialInds, ats, printMe)
		if len(vals) > 0:
			neighInds += vals

	return neighInds

#Returns a numpy array with the cluster size that each peptide is part of (eg. [2, 2, 3, 3, 1, 3, ...] means peptides 0 and 1 are part of dimers, peptides 2, 3, and 5 are part of trimers, peptide 4 is a monomer, etc...). t is the time frame of trajectory xtc to look at from run tpr, outgro is the temporary gro file to write to when getting atom positions, cutoff is the minimum distance between two atoms in a peptide to define the peptides as being part of the same cluster, and ats is the number of atoms in a peptide.
def sizeprof(t, xtc, tpr, outgro, cutoff, ats):
	peps = getPos(t, xtc, tpr, outgro)
	os.system('rm ' + outgro)
	pots = range(len(peps)/3/ats)
	sizes = np.zeros(len(peps)/3/ats)
	while len(pots) > 0:
		init = pots[0]
		pots.remove(init)
		clust = getClust(init, cutoff, peps, pots, ats, False) + [init] #getClust changes pots during execution.
		for at in clust:
			sizes[at] = len(clust)
			
	return sizes
			
#Should be the same as sizeprof except it's used in the framework of using trjconv to separate all frames into .gro files an then perform analysis on all of them. When all frames are wanted, this is faster than the other way.
def sizeprof2(ind, grobase, cutoff, ats):
	peps = getpos2(ind, grobase)
	pots = range(len(peps)/3/ats)
	sizes = np.zeros(len(peps)/3/ats)
	while len(pots) > 0:
		init = pots[0]
		pots.remove(init)
		clust = getClust(init, cutoff, peps, pots, ats, False) + [init] #getClust changes pots during execution.
		for at in clust:
			sizes[at] = len(clust)
	return sizes

#Given a list of cluster sizes (e.g. [1, 1, 5, 2, 2, 1, ...] corresponds to a monomer, a monomer, a pentamer, a dimer, a dimer, a monomer, etc...) returns a histogram of those cluster sizes (e.g. [50, 20, 10, ...] corresponds to 50 monomers, 20 dimers, 10 trimers, etc...)
def getHist(myList):
	return(np.array([[i, myList.count(i)*i] for i in range(1, 64)]))

#Calculates the un-normalized transition matrix for transitions from one cluster size to another at each step. sizeray is a list of arrays with each array containing the cluster size that monomer i (numbered from zero) was in at that given timestep (so each array occurs at a different time in the simualtion) where i is the index of a given element in the array (eg. [[1,1,1], [2,2,1]] would mean all three peptides were monomers initially but at step 2 peptides 0 and 1 dimerized). inc is at what increments these transitions should be kept (so if inc=5 then tranistions from sizeray[0] to sizeray[5], sizeray[5] to sizeray[10], etc... will be the transitions observed). el0 is the first element to look at and elf is the last element to look at. The return is a matrix whose rows correspond to sizes transitioned from and columns correspond to sizes transitioned to, and whose elements correspond to the number of times such a transition occured over the course of the simulation. 
def sizechanges(sizeray, inc, el0=0, elf=-1):
	if elf==-1:
		elf = len(sizeray)
	ceil = max([max(sizeray[i]) for i in range(len(sizeray))])
	dclusts = []
	for i in range(el0, elf-inc, inc):
		T=np.zeros((ceil,ceil))
		for j in range(len(sizeray[i])):
			T[sizeray[i][j]-1][sizeray[i+inc][j]-1] += 1
		dclusts.append(T)

	return dclusts #dclusts is a list of cluster changes between each timestep with the columns of each element corresponding to the state transitioned to and the rows corresponding to the states transitioned from (e.g. [[100, 20, 1, 0, ...], [0, 2, 2, 0,...], ...] would be one elemnet of dclusts and corresponds to 100 monomers remaining monomers, 20 monomers becoming dimers, 1 monomer and dimer forming a trimer, and 1 dimer staying a dimer).

#Same as sizechanges, except all steps are checked at the given increment rather than only increment steps (so given inc=5, transitions between elements 0 and 5, 1 and 6, 2 and 7, etc... are checked instead of just transitions between elements 0 and 5, 5 and 10, 10 and 15, etc...).
def sizechanges2(sizeray, inc, el0=0, elf=-1):
	if elf==-1:
		elf=len(inc)
	ceil = max([max(sizeray[i]) for i in range(len(sizeray))])
	dclusts = []
	for i in range(el0, elf):
		if i+inc < len(sizeray):
			T=np.zeros((ceil,ceil))
			for j in range(len(sizeray[i])):
				T[sizeray[i][j]-1][sizeray[i+inc][j]-1] += 1
			dclusts.append(T)

	return dclusts

#Expands dclusts (see sizechanges) so that all matrices are the same size (in the return of sizechanges each dclust will have a size that is no larger than the largest cluster formed up until that point; the return of this will make all matrices be as large as the largest cluster formed over the entire simulation).
def expdclusts(dclusts):
	maxl = max([len(dclusts[i]) for i in range(len(dclusts))])
	ret = []
	for i in range(len(dclusts)):
		ret.append(np.zeros((maxl,maxl)))
	for i in range(len(dclusts)):
		for j in range(len(dclusts[i])):
			for k in range(len(dclusts[i])):
				ret[i][j][k] = dclusts[i][j][k]

	return ret

#Normalizes a transition cluster by row (so rows of the input should correspond to state transitioned from and columns of the input should correspond to state transitioned to).
def normclust(clust):
	return(np.array([clust[i]/(float(sum(clust[i])) if sum(clust[i]) > 0 else 1) for i in range(len(clust))]))

#Saves dclusts, which is a series of differences in clusters, each of which is an array whose elements correspond to the number of times in a given step a cluster of size [row-number] transitioned to a cluster of size [column-number]. It is saved to a text file such that columns are delimited by spaces, rows by \n, and matrices by an additional \n. This is not as efficient as it could be but is more readable.
def savedclusts(dclusts, fnm = 'out.txt'):
	outstr = ''
	for i in range(len(dclusts)):
		for j in range(len(dclusts[i])):
			for k in range(len(dclusts[i][j])):
				outstr += '%s ' % int(dclusts[i][j][k])
			outstr +='\n'
		outstr += '\n'
	
	with open(fnm, 'w') as myf:
		myf.write(outstr)

#Returns a list of clusters sizes at each timestep in tlist. Each clustersize is a list of sizes 
def getsizes(tlist, xtc, tpr, outgro, cutoff, ats):
	return [sizeprof(tlist[i], xtc, tpr, outgro, cutoff, ats) for i in range(len(tlist))]

#Same as getsizes, but opperates assuming .xtc file has already been divided into .gro files by trjconv.
def getsizes2(grobase, cutoff, ats):
	fs = os.listdir('.')
	try:
		fmax = max([float(fs[i][len(grobase):fs[i].index('.')]) for i in range(len(fs)) if ((fs[i].find('.gro') > 0) and (fs[i].find(grobase) == 0))])
		return([sizeprof2(i, grobase, cutoff, ats) for i in range(int(fmax))])
	except ValueError:
		print('Error in file names. File %s*.gro other than what was produced by trjconv.')
		fmax = int(raw_input('To continue, delete th(is/ese) offending file(s) and run again. Alternately, input largest file number output by trjconv: '))
		return([sizeprof2(i, grobase, cutoff, ats) for i in range(fmax)])

#Saves a list of cluster sizes for all time steps in tlist from run trajectory xtc.
def sizerun(tlist, xtc, tpr, cutoff, ats, fnm='out.txt', outgro = 'temp.gro'):
	savedclusts([getsizes(tlist, xtc, tpr, outgro, cutoff, ats)], fnm=fnm)

#Same as sizerun, but divides the .xtc file into separate .gro files with trjconv. This should be faster when all frames are needed, but does not have the option of using fewer than all.
def sizerun2(xtc, tpr, cutoff, ats, fnm='out.txt', outgro = 'temp.gro', ind=''):

	grobase = outgro[:outgro.index('.')]
	if ind == '':
		os.system('trjconv -s %s -f %s -o %s -sep' % (tpr, xtc, outgro))
	else:
		os.system('trjconv -s %s -f %s -o %s -n %s -sep' % (tpr, xtc, outgro, ind))

	savedclusts([getsizes2(grobase, cutoff, ats)], fnm=fnm)
	os.system('rm %s*.gro' % grobase)
	return

#transfix loops through the size lists in fnm and eliminates any size changes that exist for less than stepmin steps (and so are transient).
def transfix(fnm, stepmin=5, ofnm = 'out.txt'):
	dat = readdclusts(fnm=fnm)[0]
	for i in range(len(dat[0])):
		check = dat[0][i]
		newcheck = -1
		checkind = []
		checkind2 = []
		for j in range(1, len(dat)):
			if dat[j][i] != check: #Cluster size of a given peptide may have changed
				if newcheck == -1: #This is a new change we are investigating.
					newcheck = dat[j][i]
					checkind.append(j)
				else:
					if dat[j][i] == newcheck: #The cluster size has persisted
						checkind.append(j)
						if len(checkind) > stepmin: #It has persisted for long enough to count
							check = newcheck #It's the new norm.
							checkind = []
							newcheck = -1
							if len(checkind2) > 0: #any sizes in checkind2 were transient
								for ind in checkind2:
									dat[ind][i] = check
								checkind2 = []
					else: #While investigating one cluster size change, another has happened
						checkind2.append(j) #Look at this secondary change (this will count even if this is a fourth different cluster size (ie 1 1 2 3 4 will have two counted for checkind2). This should not happen frequently and when it does, the system should have equilibrated into the steady state cluster size after stepmin steps.
						if len(checkind2) > stepmin: #This secondary investigation proved to be correct (ie 1 1 2 3 2 3 3 3 3 3 3 3, so the 2's will be replaced with 3's)
							for ind in checkind: #checkind must have been a short time transition
								dat[ind][i] = dat[j][i]
							check = dat[j][i] #it's the new norm
							newcheck = -1
							checkind = []
							checkind2 = []
			else: #cluster size is the same as it was
				if len(checkind) > 0: #All other changes must have been transient. This won't quite be accurate in the rare case of 1 1 2 2 2 2 1 2 2 2 2 1 2 2 2 2, etc... Eventually it should settle one way or another.
					for ind in checkind:
						dat[ind][i] = check
					checkind = []
					newcheck = -1
				if len(checkind2) > 0:
					for ind in checkind2:
						dat[ind][i] = check
					checkind2 = []

	savedclusts([dat], fnm=ofnm)

	return

	
#The reverse operation of savedclusts. It reads a file written by savedclusts and returns information of the forms [np.array(mat0), np.array(mat1), np.array(mat2), ...]
def readdclusts(fnm = 'out.txt'):
	with open(fnm, 'r') as myf:
		mylns = myf.read().splitlines()

	ret = []
	tempray = []
	it=0
	for i in range(len(mylns)):
		if mylns[i] == '':
			ret.append(np.array(tempray))
			tempray=[]
		else:
			tempray.append([])
			for elm in mylns[i].split():
				tempray[-1].append(int(elm))

	return ret

#Converts a file of clusters to a file of sizes. Sizes are much faster to work with so this was done to convert them.
def clusts2sizes(iname, oname):
	fullclusts=readdclusts(iname)
	atmax = max([max(fullclusts[0][i]) for i in range(len(fullclusts[0]))])
	fullsizes = np.zeros((len(fullclusts), atmax+1))
	for i in range(len(fullclusts)):
		clusts = fullclusts[i]
		for j in range(len(clusts)):
			clust = clusts[j]
			for k in range(len(clust)):
				fullsizes[i][clust[k]] = len(clust)

	savedclusts([fullsizes], fnm=oname)

#Gets the infinitesimal generator matrix (or as close to it as possible given a minimum timestep dt) from a cluster change matrix (rows correspond to states transitioned out and and columns states transitioned to; see sizechanges). dt is the time difference between steps.
def getQ(dclust, dt=1.):
	pclust = normclust(dclust)
	Q = np.zeros((len(pclust),len(pclust)))

	for i in range(len(pclust)):
		for j in range(len(pclust)):
			Q[i][j] = pclust[i][j]/dt

	for i in range(len(Q)):
		Q[i][i] = -sum([Q[i][j] for j in range(len(Q)) if j != i])

	return Q

#Gets the infinitesimal generator matrix using Inamura 2.11 from a list of cluster change matrices dclusts (for each element: rows correspond to states transitioned out of and columns states transitioned to; see sizechanges). dt is the time difference between steps.
def getQ2(dclusts, dt=1.):
	htimes = np.array([sum([dclusts[i][k][j] for k in range(len(dclusts[0][0])) for i in range(len(dclusts))])*dt for j in range(len(dclusts[0]))]) #This should be a 1d array of holding times in each state. Should this be weight or number? Since I am counting N with number I will count rate with number too. As long as I am consistent it should be effectively the same (just a global rescaling different).
	for i in range(len(htimes)):
		if htimes[i] == 0:
			htimes[i] = 1
	dclust = sumclust(dclusts)
	#print("dclust:")
	#print(dclust)
	Qoff = [[dclust[i][j]/float(htimes[i]) if i!=j else 0 for j in range(len(dclust[i]))] for i in range(len(dclust))]
	#print("before:")
	#print(Qoff)
	for i in range(len(htimes)):
		if htimes[i] <= dt*(i+1):
			#print(i)
			#print(htimes)
			for j in range(len(Qoff)):
				Qoff[i][j] = 0.0
				Qoff[j][i] = 0.0
	#print("after:")
	#print(Qoff)

	#Qoff = np.array([[dclust[i][j]/htimes[i] for j in range(len(dclust[i]))] for i in range(len(dclust))]) #Not sure if it is faster to have the if or the extra calculations. Not a bottleneck xxx so don't worry too much.
	#return np.array([[Qoff[i][j] if i!=j else if i > 0 if i < len(Qoff)-1 sum(Qoff[i][:i] + Qoff[i][i+1:]) else sum(Qoff[:i]) else sum(Qoff[i:]) for i in range(len(dclusts))] for j in range(len(dclusts[0]))]) #Can't think of a good way to compress things like this and still be ok on the edges. It is a mess anyway.
	for i in range(len(Qoff)):
		if i > 0:
			if i < len(Qoff)-1:
				Qoff[i][i] = -sum(Qoff[i][:i] + Qoff[i][i+1:])
			else:
				Qoff[i][i] = -sum(Qoff[i][:-1])

		else:
			Qoff[i][i] = -sum(Qoff[i][1:])
	
	return np.array(Qoff)


#Mutiplies num copies of mat together
def matpow(mat, num):
	ret = mat
	for i in range(1, num):
		ret = np.dot(ret, mat)

	return ret
	
#Given a transition rate matrix Q and an initial concentration or number x, integrates forward in time by exponentiating Qt and multiplying row vector x0 from the right (so xi = x0 * e^(Q*ti)).
def qint(Q, x0, ts):
	ret = [x0]
	for t in ts:
		ret.append(np.dot(x0, spla.expm(Q*t)))

	return ret

#Integrates x0 forward in time using a discrete transition matrix formed from dclusts.
def pint(dclusts, x0, its):
	totdclust = sumclust(dclusts)
	P=normclust(totdclust)
	if len(x0) > len(P):
		x0=x0[:len(P)]
	hists=[x0]

	for i in range(its-1):
		hists.append(np.dot(hists[-1], P))

	return [hists, P]

#Given a list of cluster changes dclusts, returns the sum of all these changes (e.g. dclusts==[[[1,0],[1,1]], [[-1, 1],[1,1]]] would return [[0,1],[2,2]].
def sumclust(dclusts):
	return(np.array([[sum([dclusts[i][k][j] for i in range(len(dclusts))]) for j in range(len(dclusts[0]))] for k in range(len(dclusts[0]))]))

#Plots a series of histograms from different runs on top of each other for comparison. dclusts is a list containing a dclust (matrix of cluster changes containing the count of clustered that transitioned from a size given by the row to a size given by the column) for each run and tss is a list containing a time series for each run.
def comphists(dclusts, tss):
	simhists = [[np.array([sum(dclusts[k][i][j]) for j in range(len(dclusts[k][i]))]) for i in range(len(dclusts[k]))] for k in range(len(dclusts))]
	cols = ['b', 'g', 'r', 'c', 'm']
	for k in range(len(dclusts)):
		for j in range(min(5, len(simhists[k][0]))):
			plt.plot(tss[k], np.array([simhists[k][m][j] for m in range(len(simhists[k]))]), color = cols[j], label="%s-mer %s" % (j+1, k))

	plt.legend(loc=0)
	plt.show()

def comphists2(fnms, ts):
	fullevos = [getClustEvo(readdclusts(fnm=fnm)) for fnm in fnms]
	num = 125.
	ts=np.array(ts)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.spines['top'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

	ax1 = fig.add_subplot(221)
	ax2 = fig.add_subplot(222)
	ax3 = fig.add_subplot(223)
	ax4 = fig.add_subplot(224)

	marks = ['-', '-', '-', 'o', 'v']
	colors = ['indianred', 'goldenrod', 'lightseagreen', 'cornflowerblue']
	for j in range(3):
		ax1.plot(ts/1000., [fullevos[j][0][i][0]/num for i in range(len(fullevos[0][0]))], marks[j], color = colors[j], label = 'Simulation %s' % (j+1), markeredgecolor = colors[j], markersize=8)
		ax2.plot(ts/1000., [fullevos[j][0][i][1]/num for i in range(len(fullevos[0][0]))], marks[j], color = colors[j], label = 'Simulation %s' % (j+1), markeredgecolor = colors[j], markersize=8)
		ax3.plot(ts/1000., [fullevos[j][0][i][2]/num for i in range(len(fullevos[0][0]))], marks[j], color = colors[j], label = 'Simulation %s' % (j+1), markeredgecolor = colors[j], markersize=8)
		ax4.plot(ts/1000., [fullevos[j][0][i][3]/num for i in range(len(fullevos[0][0]))], marks[j], color = colors[j], label = 'Simulation %s' % (j+1), markeredgecolor = colors[j], markersize=8)

	plt.legend(loc=(-.6,1.5), fontsize=20)
	plt.text(-.27, .97, 'a)', transform=ax1.transAxes, fontsize=60)
	plt.text(-.30, .97, 'b)', transform=ax2.transAxes, fontsize=60)
	plt.text(-.27, .97, 'c)', transform=ax3.transAxes, fontsize=60)
	plt.text(-.30, .97, 'd)', transform=ax4.transAxes, fontsize=60)
	ax.text(.20, -.12, 'Scaled simulation time in ns', transform=ax.transAxes, fontsize=60)
	ax.set_ylabel('Mass fraction')
	ax.yaxis.labelpad = 70
	fig.subplots_adjust(wspace=.5)

	plt.show()

#Calculates the msd between a simulated (directly from simulation) histogram and an integrated (from integration forward in time from initial conditions using transition probabilities) histogram. Only uses times in intts, and computes it for all cluster sizes in mers (ie [1,2,3,4,6] would compute things based on monomers, dimers, trimers, tetramers, and hexamers).
def dev(simts, simhist, intts, inthist, mers):
	ret = []

	for mer in mers:
		ret.append(np.mean([(simhist[simts.index(intts[i])][mer-1] - inthist[i][mer-1])**2 for i in range(len(intts))]))

	return ret

#Turns a list of q matricies into a matrix of changes over time.
def qseries(qlist, dclist):
	qout = []
	for i in range(len(qlist[0])):
		row = []
		for j in range(len(qlist[0])):
			row.append([qlist[k][i][j] if (sum(dclist[k][i]) > 0)  else 'na' for k in range(len(qlist))])

		qout.append(row)

	return qout

#Returns a series of Q matrices each obtained by blocking the data in size file fnm into time blocks of size tblock, with t ranging from t0 to tf, dtmin being the minimum stepsize that can be taken in the file, dt the step size taken to obtain the data (needs to be enough larger than dt so that we are above the Markove mixing time), and sc is the value by which the Q elements are scaled.
def getqser(fnm, t0=0, tf=350000, tblock=70000, dt=700, dtmin=70, sc=1000.):
	qs = []
	dclist = []
	dat=readdclusts(fnm=fnm)
	for t in range(t0, tf-tblock+1, tblock):
		dclusts = sizechanges(dat[0], dt/dtmin, el0=t/dtmin, elf=(t+tblock)/dtmin)
		qs.append(getQ2(dclusts, dt=dt)*sc)
		dclist.append(sumclust(dclusts))

	return(qseries(qs, dclist))

#Plots a list of Q matricies by element given in qsers over times tsrange. 
def plotqser(qsers, tsrange):
	ax1=plt.subplot2grid((4,4),(0,0), colspan=4, rowspan=4)
	colors=['blue', 'green', 'red']
	tshift = (tsrange[1]-tsrange[0])/2.
	for i in range(4):
		for j in range(4):
			ax=plt.subplot2grid((4,4),(i,j))
			for n in range(len(qsers)):
				if ((i < len(qsers[n])) & (j < len(qsers[n]))):
					ts=[tsrange[k]+tshift for k in range(len(tsrange)) if qsers[n][i][j][k] !='na']
					ys=[qsers[n][i][j][m] for m in range(len(tsrange)) if qsers[n][i][j][m] != 'na']
					ax.plot(ts, ys, label='Run %s' % (n+1), color=colors[n])
					ax.set_xlim([tsrange[0],tsrange[-1]+2*tshift])
				if ((i==0) & (j==0)):
					plt.legend(loc=(4.62,-1.7), fontsize = 25)

			ax.tick_params(axis='both', labelsize=20, pad=1)

	ax1.set_xlabel('Time in ns')
	ax1.set_ylabel('Transition Rate over next 10 ns in 1/ns')
	plt.text(.40, -.12, 'Time in ns', transform=ax1.transAxes, fontsize=50)
	plt.text(-.15, 1.02, 'Transition Rate over next 10 ns in 1/ns', transform=ax1.transAxes, fontsize=50, rotation='vertical')
	plt.text(-.09, .18, 'From 4-mer', transform=ax1.transAxes, fontsize=33, rotation='vertical')
	plt.text(-.09, .46, 'From 3-mer', transform=ax1.transAxes, fontsize=33, rotation='vertical')
	plt.text(-.09, .72, 'From 2-mer', transform=ax1.transAxes, fontsize=33, rotation='vertical')
	plt.text(-.09, .98, 'From 1-mer', transform=ax1.transAxes, fontsize=33, rotation='vertical')
	plt.text(.06, -.06, 'To 1-mer', transform=ax1.transAxes, fontsize=33)
	plt.text(.32, -.06, 'To 2-mer', transform=ax1.transAxes, fontsize=33)
	plt.text(.58, -.06, 'To 3-mer', transform=ax1.transAxes, fontsize=33)
	plt.text(.85, -.06, 'To 4-mer', transform=ax1.transAxes, fontsize=33)
	
	plt.show()

#Gets an average rate matrix by blocking and averaging over blocks and files for all files listed in flist.
def qoutshort(flist, dt=700, dtmin=70, tblock=17500, t0=0, tf=350000, sc=1000):
	qs = []
	qsize=[]
	boolray=[]
	for f in flist:
		dat = readdclusts(fnm=f)[0]
		for t in np.arange(t0, tf-tblock+1, tblock):
			dclusts = sizechanges(dat, int(dt/dtmin), el0=int(t/dtmin), elf=int((t+tblock)/dtmin))
			tot = sumclust(dclusts)
			boolray.append([((sum(tot[i])>0) | (sum(tot[:,i])>0)) for i in range(len(tot))])
			qs.append(getQ2(dclusts, dt=dt)*sc)

		qsize.append(len(qs[-1]))

	clustmax = max(qsize)
	means =[]
	errs = []
	dat = []
	for i in range(clustmax):
		dattemp = []
		for j in range(clustmax):
			stat=[]
			for k in range(len(qs)):
				if ((i < len(qs[k])) & (j < len(qs[k]))):
					if boolray[k][i]: #Makes sure we aren't artificially lowering transition rates that haven't had a chance to transition. This could be wrong though: transitions into a higher size from a lower size should be zero when they are zero
						stat.append(qs[k][i][j])

			if len(stat) > 0:
				dattemp.append([np.mean(stat), np.std(stat, ddof=1)])
			else:
				dattemp.append([0,0])

		dat.append(dattemp)

	for i in range(len(dat)):
		for j in range(len(dat[i])):
			if i == j:
				num = -sum([dat[i][m][0] for m in range(len(dat[i])) if m != i])
				stdev = (sum([dat[i][m][1]**2 for m in range(len(dat[m])) if m != i]) + np.std([dat[i][m][0] for m in range(len(dat[i])) if m!= i], ddof=1)**2/(len(dat)-1))**.5

				dat[i][j] = [num, stdev]

	return dat

#Same as qoutshort except gets the standard deviation rather than the standard error.
def qoutshort_std(flist, dt=700, dtmin=70, tblock=17500, t0=0, tf=350000, sc=1000):
	qs = []
	qsize=[]
	boolray=[]
	for f in flist:
		dat = readdclusts(fnm=f)[0]
		for t in np.arange(t0, tf-tblock+1, tblock):
			dclusts = sizechanges(dat, int(dt/dtmin), el0=int(t/dtmin), elf=int((t+tblock)/dtmin))
			tot = sumclust(dclusts)
			boolray.append([((sum(tot[i])>0) | (sum(tot[:,i])>0)) for i in range(len(tot))])
			qs.append(getQ2(dclusts, dt=dt)*sc)

		qsize.append(len(qs[-1]))

	clustmax = max(qsize)
	means =[]
	errs = []
	dat = []
	for i in range(clustmax):
		dattemp = []
		for j in range(clustmax):
			stat=[]
			for k in range(len(qs)):
				if ((i < len(qs[k])) & (j < len(qs[k]))):
					if boolray[k][i]: #Makes sure we aren't artificially lowering transition rates that haven't had a chance to transition. This could be wrong though: transitions into a higher size from a lower size should be zero when they are zero
						stat.append(qs[k][i][j])

			if len(stat) > 0:
				dattemp.append([np.mean(stat), np.std(stat, ddof=1)])
			else:
				dattemp.append([0,0])

		dat.append(dattemp)

	for i in range(len(dat)):
		for j in range(len(dat[i])):
			if i == j:
				num = -sum([dat[i][m][0] for m in range(len(dat[i])) if m != i])
				stdev = (sum([dat[i][m][1]**2 for m in range(len(dat[m])) if m != i]) + np.std([dat[i][m][0] for m in range(len(dat[i])) if m!= i], ddof=1)**2)**.5

				dat[i][j] = [num, stdev]

	return dat

#Given the output of quoutshort, saves said output to a .txt file in a nice .tex format.
def saveqout(dat, fnm='q.txt', decs = 1):
	out = '\\begin{table}[!ht]\n\t\\caption{}\n\t\\begingroup\n\t\\fontsize{8}{12}\n\t\\newcolumntype{d}[1]{D{.}{.}{#1}}\n\t\\noindent\\textsc{\n\t\t\\begin{tabular}{|l||%s}\n\t\t\t\\hhline{|-%s|}\n\t\t\t & %s \\\\\\hhline{|=::%s|}\n' % ('d{1}rr|'*len(dat), '---'*len(dat), '& '.join([('\\multicolumn{3}{c|}{\\textbf{To %s-mer}} ' % (i+1)) for i in range(len(dat))]), '==='*len(dat))
	for i in range(len(dat)):
		out+=('\t\t\t\\textbf{From %s-mer} & ' % (i+1)) + '& '.join([('%.1f&$\\pm$&%.1f ' % (dat[i][j][0], dat[i][j][1])) for j in range(len(dat[i]))]) + ' \\\\\n'

	out+='\t\t\t\\hhline{|-%s|}\n\t\t\\end{tabular}\n\t}\n\t\\endgroup\n\t\\label{}\n\\end{table}' % ('---'*len(dat))
	with open(fnm, 'w') as myf:
		myf.write(out)

	return
	
#Given dat, the output of qoutshort_std, thist method calculates the slowest relaxation rate of the system by forming a matrix by drawing off diagonal elements from a gaussian distribution with the same meand and stdev as given by dat. Negative values are rejected for off diagonal elements so the distribution is technically not normal. its is the tau used in calculating the matrix. its is the number of times this process is repeated. The rate constant and its error given by these randomly drawn matrices is returned.
def eigerr(dat, dt=700, its=100):
	rlist = []
	for t in range(its):
		mat = np.zeros((len(dat), len(dat)))
		for i in range(len(dat)):
			for j in range(len(dat)):
				if i != j:
					if dat[i][j][1] != 0:
						mat[i][j] = srand(dat[i][j][0], dat[i][j][1])
					else:
						mat[i][j] = dat[i][j][0]

		for i in range(len(mat)):
			mat[i][i] = -sum([mat[i][j] for j in range(len(mat[i])) if i != j])

		rlist.append(gettime(spla.expm(dt*mat),dt))

	med = rlist[len(rlist)/2]
	slist = sorted(rlist)
	fixedlist2=np.array(slist[len(slist)/40:-len(slist)/40])
	fixedlist = np.array([rlist[i] for i in range(len(rlist)) if rlist[i] < 10*med])
	print("95 percent confidence interval: (%s, %s)" % (fixedlist2[0], fixedlist2[-1]))
	print("Mean of everything: %s" % np.mean(rlist))
	return(np.mean(fixedlist), np.std(fixedlist, ddof=1))

#Draws a shifted random number from the normal distribution with mean mu and stdev sig. It is shifted so that numbers below zero and above 2*mu are rejected.
def srand(mu, sig):
	num = np.random.normal(mu, sig)
	while ((num < 0) | (num > 2*mu)):
		num = np.random.normal(mu, sig)

	return num

def plotqt(flist, dt=700, dtmin=70, tblock=1750, t0=0, tf=35000, sc=1000, ticks=12, tscale=1000.):
	qs = []
	qsize=[]
	boolray=[]
	ts = np.arange(t0, tf-tblock+1, tblock)/tscale
	for f in flist:
		dat = readdclusts(fnm=f)[0]
		qtemp = []
		booltemp=[]
		for t in np.arange(t0, tf-tblock+1, tblock):
			dclusts = sizechanges(dat, int(dt/dtmin), el0=int(t/dtmin), elf=int((t+tblock)/dtmin))
			tot = sumclust(dclusts)
			booltemp.append([((sum(tot[i])>0) | (sum(tot[:,i])>0)) for i in range(len(tot))]) #If clusters of that size existed at this time then it is true.
			qtemp.append(getQ2(dclusts, dt=dt)*sc)

		qsize.append(len(qtemp[-1]))
		qs.append(qtemp)
		boolray.append(booltemp)

	clustmax = max(qsize)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.axis('off')
	axs = [fig.add_subplot(clustmax, clustmax, i+1) for i in range(clustmax**2)]
	boolsize = []
	for m in range(len(qs)):
		for k in range(len(ts)):
			while len(boolray[m][k]) < clustmax:
				boolray[m][k].append(False)

			qs[m][k] = np.concatenate((qs[m][k], np.zeros((len(qs[m][k]),clustmax-len(qs[m][k])), dtype=qs[m][k].dtype)), axis=1)
	
	for i in range(clustmax):
		for j in range(clustmax):
			dat = [np.mean([qs[m][k][i][j] for m in range(len(qs)) if (boolray[m][k][i])]) for k in range(len(ts))] #Check that we don't go larger than qs[simnumber] ever formed. boolray[k][m][i] checks the i component.
			try:
				errs = [np.std([qs[m][k][i][j] for m in range(len(qs)) if boolray[m][k][i]], ddof=1) for k in range(len(ts))]
			except Exception:
				errs=[0]*len(dat)
			axs[i*clustmax+j].errorbar(ts, dat, yerr=errs, color='black')

	for myax in axs:
		#myax.tick_params(axis='both', labelsize=10, pad=1)
		myax.tick_params(axis='both', labelsize=17, pad=1)
		myax.locator_params(axis='both', nbins=ticks)
		myax.set_xlim(ts[0], ts[-1])
	
	for i in range(clustmax):
		#axs[i*clustmax].text(-.45-(i%2)*.4, 1.2, 'From %s-mer' % (i+1), fontsize=20, rotation='vertical', transform=axs[i*clustmax].transAxes)
		#axs[clustmax*(clustmax-1) + i].text(.1, -.5, 'To %s-mer' % (i+1), fontsize=20, transform=axs[clustmax*(clustmax-1)+i].transAxes) #Technically the transform is the only one that matters.
		axs[i*clustmax].text(-.35-(i%2)*.2, .9, 'From %s-mer' % (i+1), fontsize=30, rotation='vertical', transform=axs[i*clustmax].transAxes)
		axs[clustmax*(clustmax-1) + i].text(.24, -.32, 'To %s-mer' % (i+1), fontsize=30, transform=axs[clustmax*(clustmax-1)+i].transAxes) #Technically the transform is the only one that matters.



	ax.text(.4, -.12, 'Time in ns', transform=ax.transAxes, fontsize=60)
	ax.text(-.15, .7, 'Weight Fraction', rotation='vertical', transform=ax.transAxes, fontsize=60)
	plt.show()

#Rather than plotting all Qs like plotqser, plots the average and standard deviation of the list. Also performs a two-tailed p-test, testing the null hypothesis that the slope of the data is zero. The looping and data access are terribly done.
def plotqserav(qsers, tsrange):

	#ax1=plt.subplot2grid((8,8),(0,0), colspan=8, rowspan=8)
	tshift = (tsrange[1]-tsrange[0])/2.
	xray = [[.7, .7, .7, .07], [.7, .7, .7, .7], [.7, .7, .7, .7], [.07, .7, .07, .07]]
	yray = [[.1, .8, .8, .8], [.8, .1, .8, .8], [.8, .8, .1, .8], [.8, .8, .8, .8]]
	qav = []
	for i in range(8):
		qadd=[]
		for j in range(8):
			#ax=plt.subplot2grid((8,8),(i,j))
			runys=[]
			runerrs=[]
			statlist=[]
			for t in range(len(tsrange)):
				stat=[]
				for n in range(len(qsers)):
					if ((i < len(qsers[n])) & (j < len(qsers[n]))):
						if qsers[n][i][j][t] != 'na':
							stat.append(qsers[n][i][j][t])
							statlist.append(qsers[n][i][j][t])

				if len(stat) > 0:
					runys.append(np.mean(stat))
					runerrs.append(np.std(stat, ddof=1))
				else:
					runys.append('na')
					runerrs.append('na')

			ts=[tsrange[k]+tshift for k in range(len(tsrange)) if runys[k] !='na']
			ys=[runys[k] for k in range(len(runys)) if runys[k] != 'na']
			errs=[runerrs[k] for k in range(len(runys)) if runys[k] != 'na']
			if len(ts) > 2:
				m, b, r, p, err = spst.linregress(ts, ys)
				#plt.text(xray[i][j], yray[i][j], 'p=%.3f' % p, transform=ax.transAxes, fontsize=20)

			#ax.errorbar(ts, ys, yerr=errs, label='Run %s' % (n+1), color='black')
			#ax.set_xlim([tsrange[0],tsrange[-1]+2*tshift])
			#ax.tick_params(axis='both', labelsize=20, pad=1)
			qaverr = 2*np.std(statlist, ddof=1)/np.sqrt(len(statlist))
			if abs(qaverr) > 0:
				msd = int(np.floor(np.log10(qaverr)))
				print("From %s to %s" % (i, j))
				print(statlist)
			else:
				msd=1
				print("From %s to %s" % (i, j))
				print(statlist)
			qadd.append('%s$\\pm$%s' % (round(np.mean(statlist), -msd), round(qaverr, -msd)))
		qav.append(qadd)

	#ax1.set_xlabel('Time in ns')
	#ax1.set_ylabel('Transition Rate over next 10 ns in 1/ns')
	#plt.text(.40, -.12, 'Time in ns', transform=ax1.transAxes, fontsize=50)
	#plt.text(-.15, 1.02, 'Transition Rate over next 10 ns in 1/ns', transform=ax1.transAxes, fontsize=50, rotation='vertical')
	#plt.text(-.09, .18, 'From 4-mer', transform=ax1.transAxes, fontsize=33, rotation='vertical')
	#plt.text(-.09, .46, 'From 3-mer', transform=ax1.transAxes, fontsize=33, rotation='vertical')
	#plt.text(-.09, .72, 'From 2-mer', transform=ax1.transAxes, fontsize=33, rotation='vertical')
	#plt.text(-.09, .98, 'From 1-mer', transform=ax1.transAxes, fontsize=33, rotation='vertical')
	#plt.text(.06, -.06, 'To 1-mer', transform=ax1.transAxes, fontsize=33)
	#plt.text(.32, -.06, 'To 2-mer', transform=ax1.transAxes, fontsize=33)
	#plt.text(.58, -.06, 'To 3-mer', transform=ax1.transAxes, fontsize=33)
	#plt.text(.85, -.06, 'To 4-mer', transform=ax1.transAxes, fontsize=33)
	
	#plt.show()
	return(qav)

def tautest(dtlist, tlist, fnm, plot=False):
	tf = tlist[-1]
	dtmin = tlist[1]-tlist[0]
	sizes = readdclusts(fnm=fnm)[0]
	print(sizechanges(sizes, 1))
	print(1/0)
	fullevo = getClustEvo(fullclusts) #returns a list of histograms at each time step and a list of cluster size differences between steps (dclusts).
	x0 = fullevo[0][0] #The [0][0] is getting the first histogram ([0] gives all histograms lists, [0] give the first histogram list in the structure [[1, 125], [2, 0], ...].
	for i in range(1, len(fullevo[0])):
		if fullevo[0][i][0] != x0[0]:
			t0 = tlist[i-1] #Sets t0 to be the time immediately before the first cluster is formed.
			break
	t0=0
	msds=[]
	Ps = []

	for dt in dtlist:
		evo = getClustEvo([fullclusts[i] for i in range(t0/dtmin, tf/dtmin, dt/dtmin)])
		partts = [tlist[i] for i in range(t0/dtmin, tf/dtmin, dt/dtmin)]
		inthist, intP =pint(evo[1], x0, len(partts)) #The [1] in evo is for getting the dclusts rather than the histograms at each step. 

		if plot:
			plt.plot(partts, [inthist[i][0] for i in range(len(inthist))], label='t=%s' % dt)
		
		msds.append(dev(tlist, fullevo[0], partts, inthist, [1,2,3]))
		Ps.append(intP)

	if plot:
		plt.legend(loc=0, fontsize=20)
		plt.plot(tlist, [fullevo[0][i][0] for i in range(len(fullevo[0]))], label='sim')
		plt.show()
	return(msds, Ps)

def tautest(flist, dtlist, tlist, plot=False):

	for fnm in flist:
		print('notfinished')

#Plots 4 subplots, each corresponding to a given agggregate size, at different values of dt given by dtlist. tf is the final time to check
def tauplot(dtlist, tf, tlist, fnm):
	dtmin = tlist[1]-tlist[0]
	tlist = np.array(tlist)
	fullclusts = sizechanges(readdclusts(fnm=fnm)[0], (tlist[1]-tlist[0]), el0=tlist[0]/(tlist[1]-tlist[0]), elf=tlist[-1]/(tlist[1]-tlist[0]))
	num = float(sum([len(fullclusts[0][i]) for i in range(len(fullclusts[0]))]))
	fullevo = getClustEvo(fullclusts)
	x0 = fullevo[0][0] #The [0][0] is getting the first histogram ([0] gives all histograms lists, [0] give the first histogram list in the structure [[1, 125], [2, 0], ...].
	for i in range(1, len(fullevo[0])):
		if fullevo[0][i][0] != x0[0]:
			t0 = tlist[i-1]
			break
	t0=0
	msds=[]
	Ps = []
	ihists = []

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.spines['top'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

	ax1 = fig.add_subplot(221)
	ax2 = fig.add_subplot(222)
	ax3 = fig.add_subplot(223)
	ax4 = fig.add_subplot(224)

	marks = ['-', '-', '-', '^', '<', '>']
	colors = ['indianred', 'goldenrod', 'lightseagreen']
	for j in range(len(dtlist)):
		dt = dtlist[j]
		evo = getClustEvo([fullclusts[i] for i in range(t0/dtmin, tf/dtmin, dt/dtmin)])
		partts = np.array([tlist[i] for i in range(t0/dtmin, tf/dtmin, dt/dtmin)])
		inthist, intP =pint(evo[1], x0, len(partts)) #The [1] in evo is for getting the dclusts rather than the histograms at each step. 
		ihists.append(inthist)
		Ps.append(intP)
		ax1.plot(partts/1000., [inthist[i][0]/num for i in range(len(inthist))], marks[j%len(marks)], color = colors[j%len(colors)], label = r'$\tau$ = %s ps' % dt)
		ax2.plot(partts/1000., [inthist[i][1]/num for i in range(len(inthist))], marks[j%len(marks)], color = colors[j%len(colors)], label = r'$\tau$ = %s ps' % dt)
		ax3.plot(partts/1000., [inthist[i][2]/num for i in range(len(inthist))], marks[j%len(marks)], color = colors[j%len(colors)], label = r'$\tau$ = %s ps' % dt)
		ax4.plot(partts/1000., [inthist[i][3]/num for i in range(len(inthist))], marks[j%len(marks)], color = colors[j%len(colors)], label = r'$\tau$ = %s ps' % dt)

	
	ax1.plot(tlist/1000., [fullevo[0][i][0]/num for i in range(len(fullevo[0]))], label='Simulation', color='black')
	ax2.plot(tlist/1000., [fullevo[0][i][1]/num for i in range(len(fullevo[0]))], label='Simulation', color='black')
	ax3.plot(tlist/1000., [fullevo[0][i][2]/num for i in range(len(fullevo[0]))], label='Simulation', color='black')
	ax4.plot(tlist/1000., [fullevo[0][i][3]/num for i in range(len(fullevo[0]))], label='Simulation', color='black')
	plt.legend(loc=(-.6,1.5), fontsize=20)
	
	plt.text(-.27, .97, 'a)', transform=ax1.transAxes, fontsize=60)
	plt.text(-.30, .97, 'b)', transform=ax2.transAxes, fontsize=60)
	plt.text(-.27, .97, 'c)', transform=ax3.transAxes, fontsize=60)
	plt.text(-.30, .97, 'd)', transform=ax4.transAxes, fontsize=60)
	ax.text(.20, -.12, 'Scaled simulation time in ns', transform=ax.transAxes, fontsize=60)
	ax.set_ylabel('Mass fraction')
	ax.yaxis.labelpad = 70
	fig.subplots_adjust(wspace=.5)

	plt.show()

	return(msds, Ps)

#Plots the numbers of clusters of size mer vs time ts over the couse of a simulation. dclusts is the usual differences in cluster sizes over time, ax is the axis on which to plot, and label is the label to give this plot.
def plotsim(dclusts, ts, ax, mer, label="", line="-", color="blue", markersize=15):
	if mer <= len(dclusts[0]):
		if label == "":
			label = str(mer) + '-mer sim'
		hists = [np.array([sum(dclusts[i][j]) for j in range(len(dclusts[i]))]) for i in range(len(dclusts))]
		hists.append(np.array([sum(dclusts[-1][:,j]) for j in range(len(dclusts[-1][:,j]))]))
		ax.plot(ts, [hists[i][mer-1] for i in range(len(hists))], line, label=label, color=color, markeredgecolor=color, markersize=markersize)
	else:
		ax.plot(ts, np.zeros(len(ts)), line, label=label, color=color, markeredgecolor=color, markersize=markersize)

def plotint(dclusts, ts, ax, mer, label="", lines=['-','-'], colors=['black', 'black'], markersizes=[15,15]):
	if mer <= len(dclusts[0]):
		if label == "":
			label = str(mer) + '-mer int'
		x0 = np.array([sum(dclusts[0][j]) for j in range(len(dclusts[0]))])
		Q = getQ2(dclusts, dt=ts[1]-ts[0])

		hists = [np.array([sum(dclusts[i][j]) for j in range(len(dclusts[i]))]) for i in range(len(dclusts))]
		hists.append(np.array([sum(dclusts[-1][:,j]) for j in range(len(dclusts[-1][:,j]))]))
		ax.plot(ts, [hists[i][mer-1] for i in range(len(hists))], lines[0], label=labels[0], color=colors[0], markeredgecolor=colors[0], markersize=markersizes[0])

		ax.plot(ts, [hists[i][mer-1] for i in range(len(hists))], label=labels[0])
		ax.plot(ts, spla.expm(xxx))

	else:
		ax.plot(ts, np.zeros(len(ts)), lines[0], label=labels[0], color=colors[0], markeredgecolor=colors[0], markersize=markersizes[0])

def plotdata(flist, tlist, norm, exe=0, merlist=[1,2,3,4], dt=700, dtmin=70, t0=0, tf=350000, xlims= []):
	
	dclustslist = [sizechanges(readdclusts(fnm=f)[0], int(dt/dtmin), el0=int(t0/dtmin), elf=int(tf/dtmin)) for f in flist]
	nclustslist = [[dclust/norm for dclust in dclusts] for dclusts in dclustslist]
	
	fig=plt.figure()
	ax1=plt.subplot2grid((2,4),(0,0), colspan=2, rowspan=2)
	ax2=plt.subplot2grid((2,4),(0,2))
	ax3=plt.subplot2grid((2,4),(0,3))
	ax4=plt.subplot2grid((2,4),(1,2))
	ax5=plt.subplot2grid((2,4),(1,3))
	
	colors=['blue', 'green', 'red', 'purple', 'orange']
	lines=['-', '--', '|', 'x']
	msizes=[15, 15, 15, 8]

	plotsim(nclustslist[exe], tlist, ax1, merlist[0], label="%s-mer" % merlist[0], color='black', line='-')
	plotsim(nclustslist[exe], tlist, ax1, merlist[1], label="%s-mer" % merlist[1], color='black', line=':')
	plotsim(nclustslist[exe], tlist, ax1, merlist[2], label="%s-mer" % merlist[2], color='black', line='-.')
	plotsim(nclustslist[exe], tlist, ax1, merlist[3], label="%s-mer" % merlist[3], color='black', line='--')
	ax1.legend(loc=0, fontsize=35)
	
	for i in range(len(flist)):
		plotsim(nclustslist[i], tlist, ax2, merlist[0], label="Run "+str(i+1), color=colors[i%5])
		plotsim(nclustslist[i], tlist, ax3, merlist[1], label="Run "+str(i+1), color=colors[i%5])
		plotsim(nclustslist[i], tlist, ax4, merlist[2], label="Run "+str(i+1), color=colors[i%5])
		plotsim(nclustslist[i], tlist, ax5, merlist[3], label="Run "+str(i+1), color=colors[i%5])

	plt.legend(loc=0, fontsize=25)
	ax1.set_ylabel('Mass Fraction')
	plt.text(.8, -.12, 'Time in ns', transform=ax1.transAxes, fontsize=60)

	plt.text(0, 1.02, 'a)', transform=ax1.transAxes, fontsize=60)
	plt.text(0, 1.05, 'b)', transform=ax2.transAxes, fontsize=35)
	plt.text(0, 1.05, 'c)', transform=ax3.transAxes, fontsize=35)
	plt.text(0, 1.05, 'd)', transform=ax4.transAxes, fontsize=35)
	plt.text(0, 1.05, 'e)', transform=ax5.transAxes, fontsize=35)
	
	ax2.tick_params(axis='both', labelsize=25, pad=1)
	ax3.tick_params(axis='both', labelsize=25, pad=1)
	ax4.tick_params(axis='both', labelsize=25, pad=1)
	ax5.tick_params(axis='both', labelsize=25, pad=1)

	if xlims != []:
		ax1.set_xlim(xlims)
		ax2.set_xlim(xlims)
		ax3.set_xlim(xlims)
		ax4.set_xlim(xlims)
		ax5.set_xlim(xlims)

	ax2.set_ylim([0,1])
	ax3.set_ylim([0,1])
	ax4.set_ylim([0,1])
	ax5.set_ylim([0,1])
	
	plt.show()

#Deletes all rows and columns of a matrix for which the row is entirely zero (so in [[1,0,1,0],[0,0,0,0],[2,0,2,0],[0,0,0,0]] the matrix [[1,1],[2,2]] will be returned). 
def delzero(mat):
	bads = []
	dim = len(mat)
	for i in range(dim):
		bad = True
		for j in range(dim):
			if mat[i][j] != 0:
				bad=False
				break
		if bad:
			bads.append(i)

	if bads == []:
		return(mat)

	ret = np.zeros((dim-len(bads), dim-len(bads)))
	iret=0
	for i in range(len(mat)):
		if i in bads:
			continue
		jret=0
		for j in range(len(mat)):
			if j in bads:
				continue
			ret[iret][jret] = mat[i][j]
			jret+=1
		iret+=1

	return ret


#Prints eigenvalues and eigenvectors of a matrix in a readable format.
def printeig(P):
	eigsyst = npla.eig(P)
	for i in range(len(eigsyst[0])):
		print("val: %s" % eigsyst[0][i])
		print(eigsyst[1][i])

#Gets the largest eigenvalue that is less than one of a matrix.
def geteig(P, tol=10**-10):
	eigs = npla.eigvals(np.transpose(P))
	eigs.sort()
	m=1

	while True:
		if m > len(eigs):
			print("No non-one eigenvalue found")
			raise Exception("Matrix may be very nearly the identity matrix.")
		if isinstance(eigs[-m], float):
			if eigs[-m] + tol < 1:
				#if m > 2:
					#print("Repeat 1 eigenvalue found")
					#printeig(P)
				return(eigs[-m])

		elif eigs[-m].imag < tol:
			if eigs[-m].real + tol < 1:
				#if m > 2:
					#print("Repeat 1 or complex eigenvalue found")
					#printeig(P)
					#print(P)
				return(eigs[-m].real)
		elif abs(eigs[-m]) + tol < 1:
			#if m > 2:
				#print("Repeat 1 or complex eigenvalue found")
				#printeig(P)
				#print(P)
			return(abs(eigs[-m]))

		m+=1

def getrate(P, tau):
	return(-np.log(geteig(P))/tau)

def gettime(P, tau):
	return(-tau/np.log(geteig(P)))

#Given a list of size files flist and a list of times ts (as well as perhaps a dt at which to calculate the integration, a dtmin for the minimum dt in the file, tf for the final time, and t0 fo the initial time), plots the mass fraction of the cluster sizes as a function of time as well as the predicted cluster sizes based on the Q matrix derived from the whole simulation.
def qintplt(flist, ts, dt=100, dtmin=10, tf=50000, t0=0, xlims = [], ptstride=1, mersper=4):
	dclustslist = [sizechanges(readdclusts(fnm=f)[0], int(dt/dtmin), el0=int(t0/dtmin), elf=int(tf/dtmin)) for f in flist]
	simlist = []
	maxlen = max([len(dclustslist[i][0]) for i in range(len(dclustslist))])
	x0 = [sum(dclustslist[0][0][i]) for i in range(len(dclustslist[0][0]))]
	x0 = np.array(x0 + [0]*(maxlen-len(x0)))
	norm = float(sum(x0))
	for dclusts in dclustslist:
		totdclust = sumclust(dclusts)
	
		simhist = [np.array([sum(dclusts[i][j]) for j in range(len(dclusts[i]))]) for i in range(len(dclusts))]
		simhist.append(np.array([sum(dclusts[-1][j,:]) for j in range(len(dclusts[-1]))]))
		simlist.append(simhist)

	Q=getQ2([[[np.sum([dclustslist[k][m][i][j] for k in range(len(dclustslist)) if ((i < len(dclustslist[k][m])) & (j < len(dclustslist[k][m])))]) for j in range(maxlen)] for i in range(maxlen)] for m in range(len(dclustslist[0]))], dt= ts[1]-ts[0])

	inthist = [np.dot(x0, spla.expm(Q*ts[i])) for i in range(len(ts))]

	cols = ['b', 'g', 'r', 'c', 'k', 'm', 'LightCoral', 'Goldenrod', 'LawnGreen', 'SteelBlue']
	lines = ['-', '--', '-.', ':']
	fig = plt.figure()
	ax = fig.add_subplot(111)

	for i in range(maxlen/mersper+ (1 if maxlen%mersper > 0 else 0)):
		for j in range(min(mersper, maxlen-mersper*i)):
			dat = [[simlist[k][m][mersper*i+j]/norm for k in range(len(simlist)) if (mersper*i+j) < len(simlist[k][m])] for m in range(len(simlist[0]))]
			tempy=[np.mean(dat[m]) for m in range(len(simlist[0]))]
			temperr=[np.std(dat[m], ddof=1) for m in range(len(simlist[0]))]
			plt.errorbar([ts[z] for z in range(0, len(ts), ptstride)], [tempy[z] for z in range(0, len(ts), ptstride)], yerr=[temperr[z] for z in range(0, len(ts), ptstride)], color=cols[j], label="%s-mer" % (mersper*i+j+1))
			#plt.plot(ts/1000., np.array([inthist[k][mersper*i+j]/norm for k in range(len(inthist))]), '--', color = cols[j], label="%s-mer Markov dt=%s" % (mersper*i+j+1, ts[1]-ts[0]))
			plt.plot(ts, np.array([inthist[k][mersper*i+j]/norm for k in range(len(inthist))]), '--', color = cols[j])


		plt.text(.4, -.12, 'Time in ns', transform=ax.transAxes, fontsize=60)
		plt.ylabel('Mass Fraction')
		plt.legend(loc=0, fontsize=25)
		if xlims != []:
			plt.xlim(xlims[0], xlims[1])
		plt.ylim(ymin=0)
		plt.show()

#Takes a list of dclusts and sums it so that all transitions are added up (ie dclusts1 + dclusts2 + dclusts3 etc... with standard matrix addition).
def sumdclustlist(dclustslist):
	maxlen = max([len(dclustslist[i][0]) for i in range(len(dclustslist))])
	return [[[sum([dclustslist[k][m][i][j] for k in range(len(dclustslist)) if ((i < len(dclustslist[k][0])) & (j < len(dclustslist[k][0])))]) for j in range(maxlen)] for i in range(maxlen)] for m in range(len(dclustslist[0]))]

#Given a list of files containing sizes, plots mass fraction of cluster sizes as a function of time along with the initial conditions integrated forward using CTMC with dt0 as tau and DTMC with dt0 and dt1 (to show Chap-Kol. eqn is satisfied).
def pqintfull(flist, dt0, dt1, t0=0, tf=350000, dtmin=70, tscale=1000., xlims = [], ptstride=1, dim = [2,2], mers=[], yticks=12, digs=2, legsize=40):
	dclustslist=[sizechanges(readdclusts(fnm=flist[i])[0], int(dt0/dtmin), el0=int(t0/dtmin), elf=int(tf/dtmin)) for i in range(len(flist))]
	ts = np.array(np.arange(t0, tf, dt0))
	dclustslist2=[sizechanges(readdclusts(fnm=flist[i])[0], int(dt1/dtmin), el0=int(t0/dtmin), elf=int(tf/dtmin)) for i in range(len(flist))]
	ts2= np.array(np.arange(t0, tf, dt1))
	simlist = []
	maxlen = max([len(dclustslist[i][0]) for i in range(len(dclustslist))])
	x0 = [sum(dclustslist[0][0][i]) for i in range(len(dclustslist[0][0]))]
	x0 = np.array(x0 + [0]*(maxlen-len(x0)))
	norm = sum(x0)
	if mers == []:
		mers = range(1, maxlen+1)
	for dclusts in dclustslist:
		totdclust = sumclust(dclusts)
	
		simhist = [np.array([sum(dclusts[i][j]) for j in range(len(dclusts[i]))]) for i in range(len(dclusts))]
		simhist.append(np.array([sum(dclusts[-1][j,:]) for j in range(len(dclusts[-1]))]))
		simlist.append(simhist)

	Q=getQ2([[[np.sum([dclustslist[k][m][i][j] for k in range(len(dclustslist)) if ((i < len(dclustslist[k][m])) & (j < len(dclustslist[k][m])))]) for j in range(maxlen)] for i in range(maxlen)] for m in range(len(dclustslist[0]))], dt= ts[1]-ts[0])
	tot = sumdclustlist(dclustslist)
	P1 = normclust(sumclust(tot))
	tot = sumdclustlist(dclustslist2)
	P2 = normclust(sumclust(tot))
	inthist1 = [np.dot(x0, spla.expm(Q*ts[i])) for i in range(len(ts))]
	inthist2 = [x0]
	inthist3 = [x0]
	for i in range(1, len(ts)):
		inthist2.append(np.dot(inthist2[-1], P1))
		if i % 4 == 0:
			inthist3.append(np.dot(inthist3[-1], P2))

	cols = ['b', 'g', 'r', 'c', 'm']
	lines = ['-', '--', '-.', ':']

	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.axis('off')
	axs = [fig.add_subplot(dim[0], dim[1], i+1) for i in range(dim[0]*dim[1])]
	for i in range(dim[0]):
		for j in range(dim[1]):
			num = dim[1]*i + j
			if num < len(mers):
				dat = [[simlist[k][m][mers[num]-1]/norm for k in range(len(simlist)) if (mers[num]-1) < len(simlist[k][m])] for m in range(len(simlist[0]))]
				tempts = ts/tscale
				tempys = [np.mean(dat[m]) for m in range(len(simlist[0]))]
				temperrs = [np.std(dat[m], ddof=1) for m in range(len(simlist[0]))]
				axs[num].errorbar([(ts/tscale)[z] for z in range(0, len(ts), ptstride)], [tempys[z] for z in range(0, len(ts), ptstride)], yerr=[temperrs[z] for z in range(0, len(ts), ptstride)], color='k', label='Simulation')
				axs[num].plot(ts/tscale, np.array([inthist1[k][mers[num]-1]/norm for k in range(len(inthist1))]), '-.', color = 'b', label=r'CTMC ($\tau$=%s ps)' % int(dt0))
				axs[num].plot(ts/tscale, np.array([inthist2[k][mers[num]-1]/norm for k in range(len(inthist2))]), '--', color = 'r', label=r'DTMC ($\tau$=%s ps)' % int(dt0))
				axs[num].plot([ts[m]/tscale for m in range(len(ts)) if m%4==0], np.array([inthist3[k][mers[num]-1]/norm for k in range(len(inthist3))]), '--', color = 'g', label=r'DTMC ($\tau$=%s ps)' % int(dt1)) #What is the m%4 about?
				axs[num].text(-.18+dim[1]*.02, .97, string.ascii_lowercase[num:num+1] + ')', transform=axs[num].transAxes, fontsize=40-5*dim[1])


	ax.text(.4, -.12, 'Time in ns', transform=ax.transAxes, fontsize=60)
	ax.text(-.15, .7, 'Weight Fraction', rotation='vertical', transform=ax.transAxes, fontsize=60)
	
	for myax in axs:
		myax.tick_params(axis='both', labelsize=40-6*dim[1], pad=1)
		myax.locator_params(axis='y', nbins=yticks)

	if xlims != []:
		for myax in axs:
			myax.set_xlim(xlims)
	
	"""
	for myax in axs:
		ys = myax.get_yticks()
		myax.set_yticks(np.linspace(ys[0], ys[-1], yticks))
		myax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
	"""
	
	axs[0].set_ylim([0,1])

	plt.legend(loc=0, fontsize=legsize)
	plt.show()

#pqintfull before I started tinkering with mers and plot dims
def pqintfull_old(flist, dt0, dt1, t0=0, tf=350000, dtmin=70, tscale=1000., xlims = [], ptstride=1):
	dclustslist=[sizechanges(readdclusts(fnm=flist[i])[0], int(dt0/dtmin), el0=int(t0/dtmin), elf=int(tf/dtmin)) for i in range(len(flist))]
	ts = np.array(np.arange(t0, tf, dt0))
	dclustslist2=[sizechanges(readdclusts(fnm=flist[i])[0], int(dt1/dtmin), el0=int(t0/dtmin), elf=int(tf/dtmin)) for i in range(len(flist))]
	ts2= np.array(np.arange(t0, tf, dt1))
	simlist = []
	maxlen = max([len(dclustslist[i][0]) for i in range(len(dclustslist))])
	x0 = [sum(dclustslist[0][0][i]) for i in range(len(dclustslist[0][0]))]
	x0 = np.array(x0 + [0]*(maxlen-len(x0)))
	norm = sum(x0)
	for dclusts in dclustslist:
		totdclust = sumclust(dclusts)
	
		simhist = [np.array([sum(dclusts[i][j]) for j in range(len(dclusts[i]))]) for i in range(len(dclusts))]
		simhist.append(np.array([sum(dclusts[-1][j,:]) for j in range(len(dclusts[-1]))]))
		simlist.append(simhist)

	Q=getQ2([[[np.sum([dclustslist[k][m][i][j] for k in range(len(dclustslist)) if ((i < len(dclustslist[k][m])) & (j < len(dclustslist[k][m])))]) for j in range(maxlen)] for i in range(maxlen)] for m in range(len(dclustslist[0]))], dt= ts[1]-ts[0])
	tot = sumdclustlist(dclustslist)
	P1 = normclust(sumclust(tot))
	tot = sumdclustlist(dclustslist2)
	P2 = normclust(sumclust(tot))
	inthist1 = [np.dot(x0, spla.expm(Q*ts[i])) for i in range(len(ts))]
	inthist2 = [x0]
	inthist3 = [x0]
	for i in range(1, len(ts)):
		inthist2.append(np.dot(inthist2[-1], P1))
		if i % 4 == 0:
			inthist3.append(np.dot(inthist3[-1], P2))

	cols = ['b', 'g', 'r', 'c', 'm']
	lines = ['-', '--', '-.', ':']

	for i in range(maxlen/4+ (1 if maxlen%4 > 0 else 0)):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.axis('off')
		ax1 = fig.add_subplot(221)
		ax2 = fig.add_subplot(222)
		ax3 = fig.add_subplot(223)
		ax4 = fig.add_subplot(224)
		axs = [ax1, ax2, ax3, ax4]
		for j in range(min(4, maxlen-4*i)):
			dat = [[simlist[k][m][4*i+j]/norm for k in range(len(simlist)) if (4*i+j) < len(simlist[k][m])] for m in range(len(simlist[0]))]
			tempts = ts/tscale
			tempys = [np.mean(dat[m]) for m in range(len(simlist[0]))]
			temperrs = [np.std(dat[m], ddof=1) for m in range(len(simlist[0]))]
			axs[j].errorbar([(ts/tscale)[z] for z in range(0, len(ts), ptstride)], [tempys[z] for z in range(0, len(ts), ptstride)], yerr=[temperrs[z] for z in range(0, len(ts), ptstride)], color='k', label='Simulation')
			#plt.plot(ts/tscale, np.array([inthist[k][4*i+j]/norm for k in range(len(inthist))]), '--', color = cols[j], label="%s-mer Markov dt=%s" % (4*i+j+1, ts[1]-ts[0]))
			axs[j].plot(ts/tscale, np.array([inthist1[k][4*i+j]/norm for k in range(len(inthist1))]), '-.', color = 'b', label=r'CTMC ($\tau$=%s ps)' % int(dt0))
			axs[j].plot(ts/tscale, np.array([inthist2[k][4*i+j]/norm for k in range(len(inthist2))]), '--', color = 'r', label=r'DTMC ($\tau$=%s ps)' % int(dt0))
			axs[j].plot([ts[m]/tscale for m in range(len(ts)) if m%4==0], np.array([inthist3[k][4*i+j]/norm for k in range(len(inthist3))]), '--', color = 'g', label=r'DTMC ($\tau$=%s ps)' % int(dt1))


		plt.text(.4, -.12, 'Time in ns', transform=ax.transAxes, fontsize=60)
		plt.text(-.15, .7, 'Weight Fraction', rotation='vertical', transform=ax.transAxes, fontsize=60)
	
		plt.text(-.18, .97, 'a)', transform=ax1.transAxes, fontsize=40)
		plt.text(-.18, .97, 'b)', transform=ax2.transAxes, fontsize=40)
		plt.text(-.18, .97, 'c)', transform=ax3.transAxes, fontsize=40)
		plt.text(-.18, .97, 'd)', transform=ax4.transAxes, fontsize=40)

		ax1.tick_params(axis='both', labelsize=40, pad=1)
		ax2.tick_params(axis='both', labelsize=40, pad=1)
		ax3.tick_params(axis='both', labelsize=40, pad=1)
		ax4.tick_params(axis='both', labelsize=40, pad=1)

		if xlims != []:
			ax1.set_xlim(xlims)
			ax2.set_xlim(xlims)
			ax3.set_xlim(xlims)
			ax4.set_xlim(xlims)
		ax1.set_ylim([0,1])

		plt.legend(loc=0, fontsize=20)
		plt.show()

#Test me. 
def qintfull3(dclustslist, ts):
	for dclusts in dclustslist:
		totdclust = sumclust(dclusts)
		Q=getQ(totdclust, dt=(ts[1]-ts[0]))
		print(Q)
	
		x0 = np.array([sum(dclusts[0][i]) for i in range(len(dclusts[0]))])
		simhist = [np.array([sum(dclusts[i][j]) for j in range(len(dclusts[i]))]) for i in range(len(dclusts))]
		simhist.append(np.array([sum(dclusts[-1][j,:]) for j in range(len(dclusts[-1]))]))
		inthist = [np.dot(x0, spla.expm(Q*ts[i])) for i in range(len(ts))]

		cols = ['b', 'g', 'r', 'c', 'm']
		lines = ['-', '--', '-.', ':']

		for i in range(len(simhist[0])/4+ (1 if len(simhist[0])%4 > 0 else 0)):
			for j in range(min(4, len(simhist[0])-4*i)):
				plt.plot(ts, np.array([simhist[k][4*i+j]/norm for k in range(len(simhist))]), '--', color=cols[j], label="%s-mer simulation" % (4*i+j+1))
				plt.plot(ts, np.array([inthist[k][4*i+j]/norm for k in range(len(inthist))]), color = cols[j], label="%s-mer Markov dt=10" % (4*i+j+1))


	plt.show()

if __name__ == '__main__':
	ats = 146
	cutoff = .5
	myxtc = 'sd40ns.xtc'
	mytpr = 'sd100.tpr'
	tmax = 40000
	dt=5000
	ts = range(0,tmax+1,dt)
	


	dclust_2 = cls.readclust(fnm = 'sd_2_40ns_clust10.txt')
	dclust_2_100 = cls.readclust(fnm = 'sd_2_40ns_clust100.txt')
	dclust_2_500 = cls.readclust(fnm = 'sd_2_40ns_clust500.txt')
	dclust_3 = cls.readclust(fnm = 'sd_3_40ns_clust10.txt')
	dclust_3_100 = cls.readclust(fnm = 'sd_3_40ns_clust100.txt')
	dclust_3_500 = cls.readclust(fnm = 'sd_3_40ns_clust500.txt')
	ts_2 = range(10,tmax+1,10)
	ts_2_100 = range(10,tmax+1,100)
	ts_2_500 = range(10,tmax+1,500)

	print("run 1, dt=100")
	print(npla.eigvals(cls.normclust(cls.sumclust(dclust100))))
	print(npla.eigvals(spla.expm(cls.getQ(cls.sumclust(dclust100), dt=100)*100)))
	print("run 1, dt=500")
	print(npla.eigvals(cls.normclust(cls.sumclust(dclust500))))
	print(npla.eigvals(spla.expm(cls.getQ(cls.sumclust(dclust500), dt=500)*500)))

	print("run 2, dt=100")
	print(npla.eigvals(cls.normclust(cls.sumclust(dclust_2_100))))
	print(npla.eigvals(spla.expm(cls.getQ(cls.sumclust(dclust_2_100), dt=100)*100)))
	print("run 2, dt=500")
	print(npla.eigvals(cls.normclust(cls.sumclust(dclust_2_500))))
	print(npla.eigvals(spla.expm(cls.getQ(cls.sumclust(dclust_2_500), dt=500)*500)))

	#cls.comphist([dclust100, dclust_2_100, dclust_3_100], [ts100, ts_2_100, ts_2_100])

	cls.qintfull2(dclust100, ts100, dclust500, ts500)

	print(1/0)

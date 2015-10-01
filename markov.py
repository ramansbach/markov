import numpy as np, time, random, os, subprocess,numpy.linalg as npla, matplotlib.pyplot as plt, scipy.linalg as spla, time, scipy.stats as spst
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
	boxL1 = float(myLns[len(myLns)-1].split()[0])
	boxL2 = float(myLns[len(myLns)-1].split()[1])
	boxL3 = float(myLns[len(myLns)-1].split()[2])
	return (np.array([[float(myLns[i][20:].split()[0]), float(myLns[i][20:].split()[1]), float(myLns[i][20:].split()[2])] for i in range(2, len(myLns)-1)]).flatten(),np.array([boxL1,boxL2,boxL3]))

#Takes a 1d list of positions, a filename, and box vectors and writes out a gro file
#molStructure is a list defining the name,type, and number of each atom (assume all molecules are the same)
def writeGro(fName,poslist,boxV,molStructure):
	f = open(fName,'w')
	atomno = len(molStructure)/2
	molno = len(poslist)/(3*atomno)
	f.write("This is a k-mer library file\n")
	f.write(str(len(poslist)/3)+"\n")
	#print atomno
	#print molno
	#print np.size(poslist)
	atomind = 1
	for m in range(molno):
		for a in range(atomno):
			
			pos = poslist[3*m*atomno+3*a:3*m*atomno+3*a+3]
			aname = molStructure[2*(a % atomno)]
			atype = molStructure[2*(a % atomno)+1]
			#ano = molStructure[2*(a % atomno)+2]
			#print pos			
			line = "%8s%7s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" % (aname,atype,atomind,pos[0],pos[1],pos[2],0.0,0.0,0.0)
			
			f.write(line+"\n")
			atomind+=1
	boxline = "%f %f %f \n" % (boxV[0],boxV[1],boxV[2])
	f.write(boxline)
	f.close()
			


#Gets the positions of all atoms in trajectory trj, run from tpr file tpr, at time t, written to output gro file outGro.
def getPos(t, trj, tpr, outGro):
	os.system('echo 0 | trjconv -f ' + trj + ' -o ' + outGro + ' -b ' + str(t) + ' -e ' + str(t) + ' -s ' + tpr)
	return readGro(outGro)[0]

#same as before except also returns box-length variable
def getPosB(t,trj,tpr,outGro):
	os.system('echo 0 | trjconv -f ' + trj + ' -o ' + outGro + ' -b ' + str(t) + ' -e ' + str(t) + ' -s ' + tpr)
	return readGro(outGro)

#same as before except also returns box-length variable and makes whole pbcs
def getPosBWhole(t,trj,tpr,outGro):
	os.system('echo 0 | trjconv -f ' + trj + ' -o ' + outGro + ' -b ' + str(t) + ' -e ' + str(t) + ' -s ' + tpr + ' -pbc whole')
	return readGro(outGro)

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

#Gets the neighbors of atom ind from a list of potential indices potentialInds (each of which are indices of the list of all peptides listed in peplist). Being neighbors is defined as two peptides having any two atoms separated by less than cutoff. ats is the number of atoms per peptide in peplist. Summary: ind: index to check, cutoff: distance that defines neighbors (as separation of any two atoms), peplist: list of all atoms, potentialInds: potential neighbors, ats: atoms per peptide. This version assumes periodic boundary conditions
def getNeighPBC(ind, cutoff, peplist, potentialInds, ats,boxlx,boxly,boxlz):
	ret = []

	cutsq = cutoff**2
	support = '#include <math.h>'

	code = """

	int i, j;
	return_val = 0;
	double x,y,z;
	for(i=0; i<Npep1[0]/3; i++){
		for(j=0; j<Npep2[0]/3; j++){
			x = pep1[3*i]-pep2[3*j];
			y = pep1[3*i+1]-pep2[3*j+1];
			z = pep1[3*i+2]-pep2[3*j+2];
			x = x - ((double) boxlx)*round(x/((double) boxlx));
			y = y - ((double) boxly)*round(y/((double) boxly));
			z = z - ((double) boxlz)*round(z/((double) boxlz));
			if (x*x + y*y + z*z < cutsq){
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
		test = weave.inline(code,['pep1', 'pep2', 'cutsq','boxlx','boxly','boxlz'], support_code = support, libraries = ['m'])
		if test == 1:
			ret.append(potentialInds[i])
	
	return ret

#Returns an array of arrays. Each inner array is a list of peptide indices that are in the same cluster. A cluster is defined as the largest list of molecules for which each molecule in the list is neighbors either directly or indirectly (neighbors of neighbors of neighbors etc...) neighbors with each other molecule. ind is the atom to check the cluster of, cutoff is the minimum distance that defines neighbors, peplist is a list of all atoms in the simulation, potentialInds is the list of indices that could possibly be neighbors with ind, ats is the number of atoms per peptide (this is important for dividing up pepList and must be constant for all peptides in pepList), and printMe is a boolean that will cause the immediate neighbors of ind to be printed if it is true (more for debuging and checking things).
def getClust(ind, cutoff, pepList, potentialInds, ats, printMe):
	#start = time.clock()
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
	#end = time.clock();
	#print end - start;
	return neighInds

#Returns an array of arrays. Each inner array is a list of peptide indices that are in the same cluster. A cluster is defined as the largest list of molecules for which each molecule in the list is neighbors either directly or indirectly (neighbors of neighbors of neighbors etc...) neighbors with each other molecule. ind is the atom to check the cluster of, cutoff is the minimum distance that defines neighbors, peplist is a list of all atoms in the simulation, potentialInds is the list of indices that could possibly be neighbors with ind, ats is the number of atoms per peptide (this is important for dividing up pepList and must be constant for all peptides in pepList), and printMe is a boolean that will cause the immediate neighbors of ind to be printed if it is true (more for debuging and checking things). Assumes PBC
def getClustPBC(ind, cutoff, pepList, potentialInds, ats, boxlx,boxly,boxlz,printMe):
	#start = time.clock()
	neighInds = getNeighPBC(ind, cutoff, pepList, potentialInds, ats,boxlx,boxly,boxlz)
	if printMe:
		print("Neighbors of " + str(ind) + " are found to have indices of: ")
		print(neighInds)

	for neighInd in neighInds:
		potentialInds.remove(neighInd)
	for neighInd in neighInds:
		vals = getClust(neighInd, cutoff, pepList, potentialInds, ats, printMe)
		if len(vals) > 0:
			neighInds += vals
	#end = time.clock();
	#print end - start;
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
		clust = getClust(init, cutoff, peps, pots, ats, False) + [init]
		for at in clust:
			sizes[at] = len(clust)
			
	return sizes

def clustInds(t,xtc,tpr,outgro,cutoff,ats,rm=True):
	#return the index of which clusters each peptide is in
	#start = time.clock();
	peps = getPos(t,xtc,tpr,outgro)
	if rm:
		os.system('rm '+outgro)
	pots = range(len(peps)/3/ats)
	inds = np.zeros(len(peps)/3/ats)
	ind = 1
	while len(pots) > 0:
		init = pots[0]
		pots.remove(init)
		clusts = getClustPBC(init,cutoff,peps,pots,ats,False) + [init]
		for clust in clusts:
			inds[clust] = ind
		ind+=1
	#end = time.clock()
	#t = end - start
	#print "time: ", t
	return inds

def hydroRadClust(t,xtc,tpr,outgro,cutoff,ats,rm=True):
	#find and return the hydrodynamic radii of all clusters at some timestep t
	Rhs = np.array([])
	(peps,box_length) = getPosB(t,xtc,tpr,outgro)
	if rm:
		os.system('rm '+outgro)
	pots = range(len(peps)/3/ats)
	inds = np.zeros(len(peps)/3/ats)
	ind = 1
	eigvals = np.array([])
	eigvecs = np.array([])
	while len(pots) > 0:
		init = pots[0]
		pots.remove(init)
		clusts = getClust(init,cutoff,peps,pots,ats,False) + [init]
		#clusts is a list of peptides that are found in the cluster
		#each index in clusts corresponds to the indices index*ats*3:(index+1)*ats*3 in peps
		pepList = np.zeros(len(clusts)*ats*3)
		curr = 0
		for clust in clusts:
			inds[clust] = ind
			
			pepList[curr*ats*3:(curr+1)*ats*3]=peps[clust*ats*3:(clust+1)*ats*3]
			curr+=1
		Rh = hydroRad(pepList)
		Rhs = np.append(Rhs,Rh)
	return Rhs
		
def writeClustLibrary(t,xtc,tpr,outgro,cutoff,ats,molStruct,rm=True):
	#write out a library of .gro files of clusters
	(peps,box_length) = getPosB(t,xtc,tpr,outgro)
	if rm:
		os.system('rm '+outgro)
	pots = range(len(peps)/3/ats)
	inds = np.zeros(len(peps)/3/ats)
	ind = 1
	clustinds = {}
	while len(pots) > 0:
		init = pots[0]
		pots.remove(init)
		clusts = getClustPBC(init,cutoff,peps,pots,ats,box_length[0],box_length[1],box_length[2],False) + [init]
		#clusts is a list of peptides that are found in the cluster
		#each index in clusts corresponds to the indices index*ats*3:(index+1)*ats*3 in peps
		pepList = np.zeros(len(clusts)*ats*3)
		curr = 0
		for clust in clusts:
			pepList[curr*ats*3:(curr+1)*ats*3]=peps[clust*ats*3:(clust+1)*ats*3]
			curr+=1
		mer = len(clusts)
		if mer in clustinds:
			clustinds[mer]+=1
		else:
			clustinds[mer] = 1
		fName = str(mer)+"mer_"+str(clustinds[mer])+".gro"
		
		writeGro(fName,pepList,box_length,molStruct)
		
def clustMakeup(distrib,nmols):
	#given a distribution of cluster sizes, output a list consisting of numbers of k-mers of up to nmols that fits the distribution reasonably well but also the correct number of molecules
	#let distrib be a numpy array
	nclusts = 0.0
	for i in range(len(distrib)):
		nclusts += distrib[i]*i
	nclusts = nmols/nclusts
	ndistrib = np.zeros(np.size(distrib))
	for i in range(len(ndistrib)):
		ndistrib[i] = round(distrib[i]*nclusts)
	nmols = 0
	nclusts = 0
	for i in range(len(ndistrib)):
		nmols += i*ndistrib[i]
		nclusts += ndistrib[i]
	return (ndistrib,nmols,nclusts)

def prepInitConds(distrib,nmols,libraryns,boxsize):
	#given a desired distribution of k-mers and a desired number of molecules, prepares a (non-solvated) .gro file with k-mers using library files
	#libraryns is a list of how many configurations are available for k-mers as library files
	(ndistrib,nmols,nclusts) = clustMakeup(distrib,nmols)
	
	#os.system('echo 0 | trjconv -f ' + trj + ' -o ' + outGro + ' -b ' + str(t) + ' -e ' + str(t) + ' -s ' + tpr)
	#cmd = 'export GMX_MAXBACKUP=-1'
	#e1 = subprocess.Popen(cmd,stdout=subprocess.PIPE)
	#e1.wait()
	ind = 0
	for k in range(len(ndistrib)):
		n = ndistrib[k]
		for i in range(int(n)):
			#randomly choose library file
			#call gromacs to insert it into box
			nfiles = libraryns[k]
			filei = random.randint(1,nfiles)
			fnamei = str(k+1)+"mer_"+str(filei)+".gro"
			if ( (k==0) and (i==0)):
				#call editconf
				command = ['editconf','-f',fnamei,'-o','init_'+str(ind)+'.gro','-bt','triclinic','-box',str(boxsize[0]),str(boxsize[1]),str(boxsize[2])]
				ind+=1
			else:
				#call genconf
				command = ['genbox','-cp','init_'+str(ind-1)+'.gro','-o','init_'+str(ind)+'.gro','-ci',fnamei,'-nmol','1','-try','50']	
				ind+=1
			p = subprocess.Popen(command,stdout=subprocess.PIPE)
			p.wait()
	for j in range(ind-1):
		os.system('rm init_'+str(j)+'.gro')
	#cmd = 'export GMX_MAXBACKUP=99'
	#e2 = subprocess.Popen(cmd,stdout=subprocess.PIPE)
	#e2.wait()
	#print "under construction"

def betaCharacter(posList,ats,cutoff,bbs,bblist):
	#characterize the "beta-sheetness" of a given cluster based on how many consecutive beads have bonds
	N = int(len(posList)/3)
	cutoff2 = cutoff*cutoff
	#aMat = np.zeros([N,N]) #adjacency matrix for cluster
	betaBonds = np.zeros(bbs)
	print "len(betabonds) = ",len(betaBonds)
	support = '#include <math.h>'
	code = """
	 double d;
	int ** aMat;
	int beta;
	int pbeta;
    int mj,mc;
         int i,j;
	aMat = (int**) malloc((N/ats)*bbs*sizeof(int*));
    for (int k = 0; k < (N/ats)*bbs; k++){
        aMat[k] = (int *) malloc((N/ats)*bbs*sizeof(int));
    }

	for (int m = 0; m < (N/ats)*bbs; m++){
		mc = m/((int) bbs);
		for (int n = 0; n < m; n++){
			mj = n/((int) bbs);
			if (mj!=mc){
				i = int(bblist[m % bbs])+int(ats)*(int(m)/int(bbs));
				j = int(bblist[n % bbs])+int(ats)*(int(n)/int(bbs));
				d = (posList[3*i]-posList[3*j])*(posList[3*i]-posList[3*j])+(posList[3*i+1]-posList[3*j+1])*(posList[3*i+1]-posList[3*j+1])+(posList[3*i+2]-posList[3*j+2])*(posList[3*i+2]-posList[3*j+2]);
				if (d < ((double) cutoff2)){
					aMat[m][n] = 1;
				}
				else{
					aMat[m][n] = 0;
				}
			}
		}

	}

	int totbonds = 0;
	for (int i = 1; i < (N/ats)*bbs; i++){
		beta = 0;
		for (int j = i; j < (N/ats)*bbs; j++){
			if (aMat[j][j-i] == 1){
				beta++;
			}
			else{
			    if (beta > 0){
				betaBonds[beta-1] = betaBonds[beta-1]+1.0;
				totbonds++;
			    }
				beta = 0;
				//pbeta = 0;
			}
		}
        if (beta > 0){
            betaBonds[beta-1] = betaBonds[beta-1]+1.0;
	    totbonds++;
        }

	}
	
	free(aMat);
	"""
	weave.inline(code,['posList','N','ats','cutoff2','betaBonds','bbs','bblist'],support_code = support,libraries=['m'])
	print "len(betabonds) is ",len(betaBonds)
    	return betaBonds

def betaClust(t,xtc,tpr,outgro,cutoff,ats,bbs,bblist,rm=True):
	#return a set of counts (normalized by the total number) of "beta-bonds" or appearing contiguous segments of bonds between atoms
	#like a beta-ladder kind of thing for each cluster
	(peps,box_length) = getPosB(t,xtc,tpr,outgro)
	#print box_length
	if rm:
		os.system('rm '+outgro)
	pots = range(len(peps)/3/ats)
	inds = np.zeros(len(peps)/3/ats)
	ind = 1
	betachars = np.empty((0,bbs),float)
	
	while len(pots) > 0:
		
		init = pots[0]
		pots.remove(init)
		clusts = getClust(init,cutoff,peps,pots,ats,False) + [init]
		#clusts is a list of peptides that are found in the cluster
		#each index in clusts corresponds to the indices index*ats*3:(index+1)*ats*3 in peps
		pepList = np.zeros(len(clusts)*ats*3)
		curr = 0
		#mass = len(clusts);
		for clust in clusts:
			pepList[curr*ats*3:(curr+1)*ats*3]=peps[clust*ats*3:(clust+1)*ats*3]
			curr+=1
		betabonds = betaCharacter(pepList,ats,cutoff,bbs,bblist)
		#print("found betabonds: ",betabonds)
		betachars = np.append(betachars,np.array([betabonds]),axis=0)
		#print(betachars)
	#print(len(betachars))
	return betachars
		

def clustMorph(t,xtc,tpr,outgro,cutoff,ats,rm=True):
	#return the index of which clusters each peptide is in
	#also the eigenvectors and eigenvalues of the gyration tensor of each cluster
	#also the hydrodynamic radius
	#start = time.clock();
	(peps,box_length) = getPosB(t,xtc,tpr,outgro)
	#print box_length
	if rm:
		os.system('rm '+outgro)
	pots = range(len(peps)/3/ats)
	inds = np.zeros(len(peps)/3/ats)
	ind = 1
	eigvals = np.array([])
	eigvecs = np.array([])
	ms = np.array([])
	Rhs = np.array([])
	Rgs = np.array([])
	while len(pots) > 0:
		init = pots[0]
		pots.remove(init)
		clusts = getClust(init,cutoff,peps,pots,ats,False) + [init]
		#clusts is a list of peptides that are found in the cluster
		#each index in clusts corresponds to the indices index*ats*3:(index+1)*ats*3 in peps
		pepList = np.zeros(len(clusts)*ats*3)
		curr = 0
		#mass = len(clusts);
		for clust in clusts:
			inds[clust] = ind
			
			pepList[curr*ats*3:(curr+1)*ats*3]=peps[clust*ats*3:(clust+1)*ats*3]
			curr+=1
		gyrationTensor = gyrTens(pepList,box_length)
		eigstuff = np.linalg.eig(gyrationTensor)
		ind+=1
		eigMorph = np.zeros(5)
		eigval = np.sort(eigstuff[0])
		eigMorph[0] = eigval[0]
		eigMorph[1] = eigval[1]
		eigMorph[2] = eigval[2]
		eigMorph[3] = eigMorph[0]+eigMorph[1]+eigMorph[2] #Rg
		eigMorph[4] = 1.5*eigMorph[2]-0.5*eigMorph[3]
		
		eigvals = np.append(eigvals,eigMorph)
		#eigvecs = np.append(eigvecs,eigstuff[1])
		Rh = hydroRad(pepList)
		Rhs = np.append(Rhs,Rh)
		mstuff = np.zeros(3)
		mass = float(len(clusts))
		#print Rh
		#print mass
		Rgs = np.append(Rgs,eigMorph[3])
		ms = np.append(ms,mass)
		#ms = np.append(ms,np.array([mass,eigMorph[3],Rh]))
	#end = time.clock()
	#t = end - start
	#print "time: ", t
	return (inds,eigvals,eigvecs,Rhs,ms,Rgs)


def gyrTens(posList,box_length):
	#compute the full gyration tensor and return it
	gT = np.zeros([3,3])
	for i in range(3):
		for j in range(3):
			gT[i][j] = gyrTensxyC(posList,i,j,box_length[i],box_length[j])
		
	return gT


	

def gyrTensxy(posList,x,y,boxlx,boxly):
    #given a list of atom positions for a cluster, find the gyration tensor entry x,y
    gxy = 0
    N = len(posList)/3
    for R in range(N):
        for S in range(N):
            V = posList[3*R+x]-posList[3*S+x]
            V = V - boxlx*round(V/boxlx)
            U = posList[3*R+y]-posList[3*S+y]
            U = U - boxly*round(U/boxly)
            #print U,V
            gxy = gxy + V*U
    gxy = gxy/(2*N**2)
    return gxy

def hydroRad(posList):
	#given a list of atom positions for a cluster, find the average theoretical hydrodynamic radius
	support = '#include <math.h>'
	code = """
	return_val = 0;
	int i,j;
	double Rh = 0;
	int npep = NposList[0]/3;
	for (i = 0; i < npep;i++){
		for (j=0; j < npep;j++){
			if (i!=j){
				Rh += 1.0/sqrt((posList[3*i]-posList[3*j])*(posList[3*i]-posList[3*j])+(posList[3*i+1]-posList[3*j+1])*(posList[3*i+1]-posList[3*j+1])+(posList[3*i+2]-posList[3*j+2])*(posList[3*i+2]-posList[3*j+2]));
			}
		}	
	}
	Rh = Rh/(npep*npep);
	Rh = 1.0/Rh;
	return_val = Rh;
	"""
	Rh = weave.inline(code,['posList'],support_code = support,libraries=['m'])
	return Rh

def gyrTensxyC(posList,x,y,boxlx,boxly):
    #same as gyrTensxy except using inline C code to test for efficiency
    support = '#include <math.h>'
    code = """
    return_val = 0;
    int R,S;
    double V,U;
    double gxy = 0;
    for (R = 0; R<NposList[0]/3;R++){
        for (S = 0; S<NposList[0]/3;S++){
            V = posList[3*R+x] - posList[3*S+x];
            V = V - boxlx*round(V/boxlx);
            U = posList[3*R+y] - posList[3*S+y];
            U = U - boxly*round(U/boxly);
            gxy = gxy+U*V;
        }
    }
    gxy = gxy/(2*(NposList[0]/3)*(NposList[0]/3));
    //return_val = Py::new_reference_to(Py::Double(gxy));
    return_val = gxy;
    """
    gxy = weave.inline(code,['posList','x','y','boxlx','boxly'],support_code = support,libraries=['m'])
    return gxy
		

#Given a list of cluster sizes (e.g. [1, 1, 5, 2, 2, 1, ...] corresponds to a monomer, a monomer, a pentamer, a dimer, a dimer, a monomer, etc...) returns a histogram of those cluster sizes (e.g. [50, 20, 10, ...] corresponds to 50 monomers, 20 dimers, 10 trimers, etc...)
def getHist(myList):
	return(np.array([[i, myList.count(i)*i] for i in range(1, 64)]))

#Calculates the un-normalized transition matrix for transitions from one cluster size to another at each step. sizeray is a list of arrays with each array containing the cluster size that monomer i (numbered from zero) was in at that given timestep (so each array occurs at a different time in the simualtion) where i is the index of a given element in the array (eg. [[1,1,1], [2,2,1]] would mean all three peptides were monomers initially but at step 2 peptides 0 and 1 dimerized). inc is at what increments these transitions should be kept (so if inc=5 then tranistions from sizeray[0] to sizeray[5], sizeray[5] to sizeray[10], etc... will be the transitions observed). el0 is the first element to look at and elf is the last element to look at. The return is a matrix whose rows correspond to sizes transitioned from and columns correspond to sizes transitioned to, and whose elements correspond to the number of times such a transition occured over the course of the simulation. 
def sizechanges(sizeray, inc, el0=0, elf=-1):
	if elf==-1:
		elf = len(sizeray)
	ceil = max([max(sizeray[i]) for i in range(len(sizeray))])
	dclusts = []
	for i in range(el0, elf-inc+1, inc):
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

#Saves a list of cluster sizes for all time steps in tlist from run trajectory xtc.
def sizerun(tlist, xtc, tpr, cutoff, ats, fnm='out.txt', outgro = 'temp.gro'):
	savedclusts([getsizes(tlist, xtc, tpr, outgro, cutoff, ats)], fnm=fnm)
	
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
def qout(flist, dt=700, dtmin=70, tblock=17500, t0=0, tf=350000, sc=1000):
	dclist = []
	qs = []
	qsize=[]
	boolray=[]
	for f in flist:
		dat = readdclusts(fnm=f)[0]
		for t in range(t0, tf-tblock+1, tblock):
			dclusts = sizechanges(dat, dt/dtmin, el0=t/dtmin, elf=(t+tblock)/dtmin)
			tot = sumclust(dclusts)
			boolray.append([((sum(tot[i])>0) | (sum(tot[:,i])>0)) for i in range(len(tot))])
			qs.append(getQ2(dclusts, dt=dt)*sc)
			#print(qs[-1])

		qsize.append(len(qs[-1]))

	clustmax = max(qsize)
	means =[]
	errs = []
	ret = []
	#tempret=[]
	dat = []
	for i in range(clustmax):
		temp = []
		#temp2 = []
		dattemp = []
		for j in range(clustmax):
			stat=[]
			for k in range(len(qs)):
				if ((i < len(qs[k])) & (j < len(qs[k]))):
					if boolray[k][i]:
						stat.append(qs[k][i][j])

			if len(stat) > 0:
				#print("i: %s j: %s " % (i, j))
				#print(stat)
				if np.std(stat) > 0:
					msd = int(np.floor(np.log10(np.std(stat)*2/np.sqrt(len(stat)))))-2
				else:
					msd=1
				temp.append(r"%s$\pm$%s" % (round(np.mean(stat), -msd), round(np.std(stat)*2/np.sqrt(len(stat)), -msd)))
				#temp2.append(r"%s$\pm$%s" % (np.mean(stat), np.std(stat)*2/np.sqrt(len(stat))))
				dattemp.append([np.mean(stat), np.std(stat)/np.sqrt(len(stat))])
			else:
				temp.append("NA")
				dattemp.append("NA")
		ret.append(temp)
		#tempret.append(temp2)
		dat.append(dattemp)

	while ret[-1][-1] == "NA": #This fixes the case when a cluster was temporarily formed during the run but did not last for long enough to be significant.
		ret = [[ret[i][j] for j in range(len(ret[i])-1)] for i in range(len(ret)-1)]
		dat = [[dat[i][j] for j in range(len(dat[i])-1)] for i in range(len(dat)-1)]
		#tempret = [[tempret[i][j] for j in range(len(tempret[i])-1)] for i in range(len(tempret)-1)]
				
	ret2=[]
	Q=[]
	for i in range(len(dat)):
		rettemp=[]
		qtemp=[]
		for j in range(len(dat[i])):
			if i != j:
				if dat[i][j][1] > 0:
					msd = int(np.floor(np.log10(dat[i][j][1])))-2
				else:
					msd=1
				#rettemp.append(r"%s$\pm$%s" % (round(dat[i][j][0], -msd), round(dat[i][j][1], -msd)))
				rettemp.append(r"%.1f$\pm$%.1f" % (dat[i][j][0], dat[i][j][1]))
				qtemp.append(dat[i][j][0])
			else:
				num = -sum([dat[i][m][0] for m in range(len(dat[i])) if m != i])
				stdev = (sum([dat[i][m][1]**2 for m in range(len(dat[m])) if m != i]) + np.std([dat[i][m][0] for m in range(len(dat[i])) if m!= i])**2/(len(dat)-1))**.5
				if stdev > 0:
					msd = int(np.floor(np.log10(stdev)))-1
				else:
					msd=1
				#rettemp.append(r"%s$\pm$%s" % (round(num, -msd), round(stdev, -msd)))
				rettemp.append(r"%.1f$\pm$%.1f" % (num, stdev))
				qtemp.append(num)

		ret2.append(rettemp)
		Q.append(qtemp)

	#temp.append(r"%s$\pm$%s" % (round(np.mean(stat), -msd), round(np.std(stat)*2/np.sqrt(len(stat)), -msd)))
	#print(tempret)
	print(geteig(spla.expm(700*[[Q[i][j] for j in range(len(Q[i]))] for i in range(len(Q))])))

	return ret2

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
					runerrs.append(np.std(stat))
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
			qaverr = 2*np.std(statlist)/np.sqrt(len(statlist))
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

#Plots a series of clusters given by a list of .txt files listed in flist against times tlist. norm is the number of molecules, exe is the example plot to overlay on the left side while all other plots will be showed together divided up by cluster sizes on the right size. merlist is a list of agregate sizes to plot.
def plotdata(flist, tlist, norm, exe=0, merlist=[1,2,3,4], dt=100, dtmin=10, t0=0, tf=40000):
	
	dclustslist = [sizechanges(readdclusts(fnm=f)[0], dt/dtmin, el0=t0/dtmin, elf=tf/dtmin) for f in flist]
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
def geteig(P):
	eigs = npla.eigvals(np.transpose(P))
	eigs.sort()
	m=1

	while True:
		if m > len(eigs):
			print("No non-one eigenvalue found")
			raise Exception("Matrix may be very nearly the identity matrix.")
		if isinstance(eigs[-m], float):
			if eigs[-m] + .0000001 < 1:
				if m > 2:
					print("Repeat 1 eigenvalue found")
					printeig(P)
				return(eigs[-m])

		elif eigs[-m].imag < .0000001:
			if eigs[-m].real + .0000001 < 1:
				if m > 2:
					print("Repeat 1 or complex eigenvalue found")
					printeig(P)
					print(P)
				return(eigs[-m].real)

		m+=1

def getrate(P, tau):
	return(-np.log(geteig(P)/tau))

#Given dclustslist (a list of serieses of matrices counting the number of times a rownumber-mer tranistioned to a columnnumber-mer) plots the histograms of cluster sizes as well as those predicted by integrating forward in time using continuous and discrete Markov processes as models.
def qintfull(dclustslist, ts):
	simlist = []
	maxlen = max([len(dclustslist[i][0]) for i in range(len(dclustslist))])
	x0 = [sum(dclustslist[0][0][i]) for i in range(len(dclustslist[0][0]))]
	x0 = np.array(x0 + [0]*(maxlen-len(x0)))
	for dclusts in dclustslist:
		totdclust = sumclust(dclusts)
	
		simhist = [np.array([sum(dclusts[i][j]) for j in range(len(dclusts[i]))]) for i in range(len(dclusts))]
		simhist.append(np.array([sum(dclusts[-1][j,:]) for j in range(len(dclusts[-1]))]))
		simlist.append(simhist)

	Q=getQ2([[[np.sum([dclustslist[k][m][i][j] for k in range(len(dclustslist)) if ((i < len(dclustslist[k][m])) & (j < len(dclustslist[k][m])))]) for j in range(maxlen)] for i in range(maxlen)] for m in range(len(dclustslist[0]))], dt= ts[1]-ts[0])
	print(Q)
	print(npla.eigvals(spla.expm(Q*700)))
	inthist = [np.dot(x0, spla.expm(Q*ts[i])) for i in range(len(ts))]

	cols = ['b', 'g', 'r', 'c', 'm', 'k', 'LightCoral', 'Goldenrod', 'LawnGreen', 'SteelBlue']
	lines = ['-', '--', '-.', ':']
	fig = plt.figure()
	ax = fig.add_subplot(111)

	for i in range(maxlen/4+ (1 if maxlen%4 > 0 else 0)):
		for j in range(min(4, maxlen-4*i)):
			dat = [[simlist[k][m][4*i+j]/125. for k in range(len(simlist)) if (4*i+j) < len(simlist[k][m])] for m in range(len(simlist[0]))]
			plt.errorbar(ts/1000., [np.mean(dat[m]) for m in range(len(simlist[0]))], yerr=[np.std(dat[m])/np.sqrt(len(dat[m])) for m in range(len(simlist[0]))], color=cols[j], label="%s-mer" % (4*i+j+1))
			#plt.plot(ts/1000., np.array([inthist[k][4*i+j]/125. for k in range(len(inthist))]), '--', color = cols[j], label="%s-mer Markov dt=%s" % (4*i+j+1, ts[1]-ts[0]))
			plt.plot(ts/1000., np.array([inthist[k][4*i+j]/125. for k in range(len(inthist))]), '--', color = cols[j])


		plt.text(.4, -.12, 'Time in ns', transform=ax.transAxes, fontsize=60)
		plt.ylabel('Mass Fraction')
		plt.legend(loc=0, fontsize=25)
		plt.show()

def sumdclustlist(dclustslist):
	maxlen = max([len(dclustslist[i][0]) for i in range(len(dclustslist))])
	return [[[sum([dclustslist[k][m][i][j] for k in range(len(dclustslist)) if ((i < len(dclustslist[k][0])) & (j < len(dclustslist[k][0])))]) for j in range(maxlen)] for i in range(maxlen)] for m in range(len(dclustslist[0]))]

def pqintfull(flist, dt0, dt1, t0=0, tf=350000, dtmin=70):
	dclustslist=[sizechanges(readdclusts(fnm=flist[i])[0], dt0/dtmin, el0=t0/dtmin, elf=tf/dtmin) for i in range(len(flist))]
	ts = np.array(range(t0, tf+1, dt0))
	dclustslist2=[sizechanges(readdclusts(fnm=flist[i])[0], dt1/dtmin, el0=t0/dtmin, elf=tf/dtmin) for i in range(len(flist))]
	ts2= np.array(range(t0, tf+1, dt1))
	simlist = []
	maxlen = max([len(dclustslist[i][0]) for i in range(len(dclustslist))])
	x0 = [sum(dclustslist[0][0][i]) for i in range(len(dclustslist[0][0]))]
	x0 = np.array(x0 + [0]*(maxlen-len(x0)))
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
			dat = [[simlist[k][m][4*i+j]/64. for k in range(len(simlist)) if (4*i+j) < len(simlist[k][m])] for m in range(len(simlist[0]))]
			axs[j].errorbar(ts/1000., [np.mean(dat[m]) for m in range(len(simlist[0]))], yerr=[np.std(dat[m])/np.sqrt(len(dat[m])) for m in range(len(simlist[0]))], color='k', label='Simulation')
			#plt.plot(ts/1000., np.array([inthist[k][4*i+j]/64. for k in range(len(inthist))]), '--', color = cols[j], label="%s-mer Markov dt=%s" % (4*i+j+1, ts[1]-ts[0]))
			axs[j].plot(ts/1000., np.array([inthist1[k][4*i+j]/64. for k in range(len(inthist1))]), '-.', color = 'b', label='Q')
			axs[j].plot(ts/1000., np.array([inthist2[k][4*i+j]/64. for k in range(len(inthist2))]), '--', color = 'r', label=r'P($\tau$)')
			axs[j].plot([ts[m]/1000. for m in range(len(ts)) if m%4==0], np.array([inthist3[k][4*i+j]/64. for k in range(len(inthist3))]), '--', color = 'g', label=r'P($4\tau$)')


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

		plt.legend(loc=0, fontsize=30)
		plt.show()

def qintfull4(flist, t0=0, tf=350000, dt=700, dtmin=70): #Should probably add tblock at some point and get Q via blocking.
	tlist = np.array(range(t0, tf+1, dt))
	dclustslist = []
	qlist = []
	for fnm in flist:
		dclustslist.append(sizechanges(readdclusts(fnm=fnm)[0], dt/dtmin, el0=t0/dtmin, elf=tf/dtmin))
		qlist.append(getQ2(dclustslist[-1], dt=dt))

	print(getQ2(sumdclustlist(dclustslist), dt=dt))
	print(qlist)
	

	print(1/0)
	


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
				plt.plot(ts, np.array([simhist[k][4*i+j]/64. for k in range(len(simhist))]), '--', color=cols[j], label="%s-mer simulation" % (4*i+j+1))
				plt.plot(ts, np.array([inthist[k][4*i+j]/64. for k in range(len(inthist))]), color = cols[j], label="%s-mer Markov dt=10" % (4*i+j+1))


	plt.show()
#Same as quintfull but with the addition that the system is also integrated forward in time using data from a second set of dclusts at a different discrete time. If a Markov process is an accurate model of the system then the continuous from dclusts, the discrete from dclusts, the discrete from dclusts2, and the actual data from dclusts should all approximately match.
def qintfull2(dclusts, ts, dclusts2, ts2):
	totdclust = sumclust(dclusts)
	totdclust2 = sumclust(dclusts2)
	P=normclust(totdclust)
	P2=normclust(totdclust2)
	Q=getQ(totdclust, dt=(ts[1]-ts[0]))

	x0 = np.array([sum(dclusts[0][i]) for i in range(len(dclusts[0]))])
	simhist = [np.array([sum(dclusts[i][j]) for j in range(len(dclusts[i]))]) for i in range(len(dclusts))]
	inthist = [np.dot(x0, spla.expm(Q*ts[i])) for i in range(len(ts))]

	simhist2=[x0]
	for i in range(len(ts)-1):
		simhist2.append(np.dot(simhist2[-1], P))

	simhist3=[x0]
	for i in range(len(ts2)-1):
		simhist3.append(np.dot(simhist3[-1], P2))

	cols = ['b', 'g', 'r', 'c', 'm']

	for i in range(len(simhist[0])/5+ (1 if len(simhist[0])%5 > 0 else 0)):
		for j in range(min(5, len(simhist[0])-5*i-1)):
			plt.plot(ts, np.array([simhist[k][5*i+j] for k in range(len(simhist))]), color = cols[j], label="%s-mer simulation" % (5*i+j+1))
			plt.plot(ts, np.array([inthist[k][5*i+j] for k in range(len(inthist))]), '--', color = cols[j], label="%s-mer continuous Markov" % (5*i+j+1))
			plt.plot(ts, np.array([simhist2[k][5*i+j] for k in range(len(simhist2))]), 'x', color = cols[j], label="%s-mer Markov dt=%s" % (5*i+j+1, ts[1]-ts[0]))
			plt.plot(ts2, np.array([simhist3[k][5*i+j] for k in range(len(simhist3))]), 'o', color = cols[j], label="%s-mer Markov dt=%s" % (5*i+j+1, ts2[1]-ts2[0]))

		plt.legend(loc=1)
		plt.show()

#No different the 2 at the moment.
def qintfull3(dclusts, ts, dclusts2, ts2):
	totdclust = sumclust(dclusts)
	totdclust2 = sumclust(dclusts2)
	P=normclust(totdclust)
	P2=normclust(totdclust2)
	Q=getQ(totdclust, dt=(ts[1]-ts[0]))

	x0 = np.array([sum(dclusts[0][i]) for i in range(len(dclusts[0]))])
	simhist = [np.array([sum(dclusts[i][j]) for j in range(len(dclusts[i]))]) for i in range(len(dclusts))]
	inthist = [np.dot(x0, spla.expm(Q*ts[i])) for i in range(len(ts))]

	simhist2=[x0]
	for i in range(len(ts)-1):
		simhist2.append(np.dot(simhist2[-1], P))

	simhist3=[x0]
	for i in range(len(ts2)-1):
		simhist3.append(np.dot(simhist3[-1], P2))

	cols = ['b', 'g', 'r', 'c', 'm']

	for i in range(len(simhist[0])/5+ (1 if len(simhist[0])%5 > 0 else 0)):
		for j in range(min(5, len(simhist[0])-5*i-1)):
			plt.plot(ts, np.array([simhist[k][5*i+j] for k in range(len(simhist))]), color = cols[j], label="%s-mer simulation" % (5*i+j+1))
			plt.plot(ts, np.array([inthist[k][5*i+j] for k in range(len(inthist))]), '--', color = cols[j], label="%s-mer continuous Markov" % (5*i+j+1))
			plt.plot(ts, np.array([simhist2[k][5*i+j] for k in range(len(simhist2))]), 'x', color = cols[j], label="%s-mer Markov dt=%s" % (5*i+j+1, ts[1]-ts[0]))
			plt.plot(ts2, np.array([simhist3[k][5*i+j] for k in range(len(simhist3))]), 'o', color = cols[j], label="%s-mer Markov dt=%s" % (5*i+j+1, ts2[1]-ts2[0]))

		plt.legend(loc=1)
		plt.show()

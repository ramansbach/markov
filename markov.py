import numpy as np, time, random, os, subprocess,numpy.linalg as npla, matplotlib.pyplot as plt, scipy.linalg as spla, time, scipy.stats as spst, scipy.spatial.distance as spsd

from scipy import weave
from operator import sub, div
#params = {'axes.titlesize':70,
#		'axes.labelsize':60,
#		'text.fontsize':60,
#		'font.size':50,
#		'lines.markersize':6,
#		'lines.linewidth':4,
#		'text.usetex':True,
#		'xtick.major.pad':7,
#		'ytick.major.pad':18}
#plt.rcParams.update(params)
#Takes a .gro file and returns a 1d list of positions, box lengths, and a 1d list of whether the molecule does
 #or does not contain a particular residue (used to tell if it's charged or uncharged)
#resInfo is a tuple ('resname',pos), where position is the index within the molecule of the residue

def recluster1(ind,coms,cutoff,box):
	
	currind = ind
	#mols = range(np.shape(coms,1)))
	inds = np.zeros(np.shape(mols))
	masses = np.zeros(np.shape(mols))
	while len(mols) > 0:
		currmol = mols[0]
		mols = mols.remove(currmol)
		clust = getClust1(currmol,mols,coms,cutoff,box)
		for c in clust:
			mols = mols.remove(c)
			inds[c] = currind
			masses[c] = len(clust)
		currind += 1
	return (inds,masses)
		

def readGroQ(fName,resInfo,ats):
        with open(fName, 'r') as myF:
                myLns = myF.read().splitlines()
        boxL1 = float(myLns[len(myLns)-1].split()[0])
        boxL2 = float(myLns[len(myLns)-1].split()[1])
        boxL3 = float(myLns[len(myLns)-1].split()[2])
        resBool = np.zeros([(len(myLns)-2)/ats,1])

        j=0
    
        for i in range(resInfo[1]+2,len(myLns)-1,ats):
            #print(myLns[i][5:8])
            if myLns[i][5:8] == resInfo[0]:
                #print(myLns[i][5:8])
                #print i
                resBool[j]=1.0
            j+=1
                #print(i/ats)
        return (np.array([[float(myLns[i][20:].split()[0]), float(myLns[i][20:].split()[1]), float(myLns[i][20:].split()[2])] for i in range(2, len(myLns)-1)]).flatten(),np.array([boxL1,boxL2,boxL3]),resBool)

def getPosQ(t,trj,tpr,outGro,resInfo,ats):
        os.system('echo 0 | trjconv -f ' + trj + ' -o ' + outGro + ' -b ' + str(t) + ' -e ' + str(t) + ' -s ' + tpr)
        #subprocess.check_output([])        
        return readGroQ(outGro,resInfo,ats)

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
def writeGro(snapno,poslist,boxV,qbools,molStructureQ=["1PAQ","BB","1PAQ","SC1","2PHE","BB","2PHE","SC1","2PHE","SC2","2PHE","SC3","3ALA","BB","4GLY","BB","5COA","BB","6COP","BB","6COP","SC1","6COP","BB","7VIN","BB","8BEN","BB","8BEN","SC1","8BEN","BB","9VIN","BB","10COA","BB","11COP","BB","11COP","SC1","11COP","BB","12GLY","BB","13ALA","BB","14PHE","BB","14PHE","SC1","14PHE","SC2","14PHE","SC3","15PAQ","BB","15PAQ","SC1"],molStructure=["1PAS","BB","1PAS","SC1","2PHE","BB","2PHE","SC1","2PHE","SC2","2PHE","SC3","3ALA","BB","4GLY","BB","5COA","BB","6COP","BB","6COP","SC1","6COP","BB","7VIN","BB","8BEN","BB","8BEN","SC1","8BEN","BB","9VIN","BB","10COA","BB","11COP","BB","11COP","SC1","11COP","BB","12GLY","BB","13ALA","BB","14PHE","BB","14PHE","SC1","14PHE","SC2","14PHE","SC3","15PAS","BB","15PAS","SC1"]
):
	atomno = len(molStructure)/2
	molno = len(poslist)/(3*atomno) 
	qno = int(np.sum(qbools)) 
	nqno = int(molno-qno)
	fName = 'Q'+str(qno)+'NQ'+str(nqno)+'mer_'+str(snapno)+'.gro' 
	tName = 'Q'+str(qno)+'NQ'+str(nqno)+'mer_'+str(snapno)+'.top'
	t = open(tName,'w')
 	t.write("#include \"martini.itp\"\n")
 	t.write("#include \"Protein.itp\"\n")
 	t.write("#include \"Protein_NP.itp\"\n")  
 	t.write("[ system ]\n\n")
 	t.write(";name\n")
 	t.write("Martini system for k-mer library\n\n")
   	t.write("[ molecules ]\n\n")
     	t.write(";name\tnumber\n")
     	t.write("Protein\t{}\n".format(nqno))
 	t.write("ProteinQ\t{}\n".format(qno))      
	t.close()
	f = open(fName,'w')
	f.write("This is a k-mer library file\n")
	f.write(str(len(poslist)/3)+"\n")
	#print atomno
	#print molno
	#print np.size(poslist)
	atomind = 1
     
	for m in range(molno):
		for a in range(atomno):
			
			pos = poslist[3*m*atomno+3*a:3*m*atomno+3*a+3]
			if qbools[m]:
				aname = molStructureQ[2*(a % atomno)]
			else:
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
	#os.system('echo 0 | trjconv -f ' + trj + ' -o ' + outGro + ' -b ' + str(t) + ' -e ' + str(t) + ' -s ' + tpr)
	
     p1 = subprocess.Popen(["trjconv","-f",trj,"-o",outGro,"-b",str(t),"-e",str(t),"-s",tpr],stdin = subprocess.PIPE)
     p1.communicate("0")
     return readGro(outGro)

#same as before except assumes it's given a grofile
def getPosGro(groname):
	return readGro(groname)

#same as before except assumes we have the .gro files already
def getPosB2(ind, grobase):
	return readGro(grobase + str(ind) + '.gro')

#same as before except also returns box-length variable and makes whole pbcs
def getPosBWhole(t,trj,tpr,outGro):
	os.system('echo 0 | trjconv -f ' + trj + ' -o ' + outGro + ' -b ' + str(t) + ' -e ' + str(t) + ' -s ' + tpr + ' -pbc whole')
	return readGro(outGro)

#Gets the neighbors of atom ind from a list of potential indices potentialInds (each of which are indices of the list of all peptides listed in peplist). Being neighbors is defined as two peptides having any two atoms separated by less than cutoff. ats is the number of atoms per peptide in peplist. Summary: ind: index to check, cutoff: distance that defines neighbors (as separation of any two atoms), peplist: list of all atoms, potentialInds: potential neighbors, ats: atoms per peptide.
def getNeighContact(ind, cutoff, peplist, potentialInds, ats):
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
 
 #python version for debugging Gets the neighbors of atom ind from a list of potential indices potentialInds (each of which are indices of the list of all peptides listed in peplist). Being neighbors is defined as two peptides having any two atoms separated by less than cutoff. ats is the number of atoms per peptide in peplist. Summary: ind: index to check, cutoff: distance that defines neighbors (as separation of any two atoms), peplist: list of all atoms, potentialInds: potential neighbors, ats: atoms per peptide.
def getNeighContactPy(ind, cutoff, peplist, potentialInds, ats):
	ret = []

	cutsq = cutoff**2

	pep1 = peplist[ind*3*ats:(ind+1)*3*ats] #Assumes all peptides have ats atoms in them. The 3 is for 3 coords in 3 dimensions.
	for i in range(len(potentialInds)):
		pep2 = peplist[potentialInds[i]*3*ats:(potentialInds[i]+1)*3*ats]
  		test = 0
  		for k in range(len(pep1)/3):
			for l in range(len(pep2)/3):
				if (pep1[3*k]-pep2[3*l])* (pep1[3*k]-pep2[3*l]) + (pep1[3*k+1]-pep2[3*l+1])* (pep1[3*k+1]-pep2[3*l+1]) + (pep1[3*k+2]-pep2[3*l+2])* (pep1[3*k+2]-pep2[3*l+2])< cutsq:
					test = 1
					break
			if test == 1:
				break
		if test == 1:
			ret.append(potentialInds[i])
	
	return ret
 

#Gets the neighbors of atom ind from a list of potential indices potentialInds (each of which are indices of the list of all peptides listed in peplist). Being neighbors is defined as two peptides having any two atoms separated by less than cutoff. ats is the number of atoms per peptide in peplist. Summary: ind: index to check, cutoff: distance that defines neighbors (as separation of any two atoms), peplist: list of all atoms, potentialInds: potential neighbors, ats: atoms per peptide.
def getNeigh(ind, cutoff, peplist, potentialInds, ats, metric = 'contact'):
	ret = []
	if metric == 'contact':     
		ret = getNeighContact(ind,cutoff,peplist,potentialInds,ats)
	elif metric == 'optical':     
		ret = getNeighContact(ind,cutoff,peplist,potentialInds,ats) 
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
def getClust(ind, cutoff, pepList, potentialInds, ats, printMe, metric = 'contact'):
	#start = time.clock()
	neighInds = getNeigh(ind, cutoff, pepList, potentialInds, ats, metric)
	if printMe:
		print("Neighbors of " + str(ind) + " are found to have indices of: ")
		print(neighInds)

	for neighInd in neighInds:
		potentialInds.remove(neighInd)
	for neighInd in neighInds:
		vals = getClust(neighInd, cutoff, pepList, potentialInds, ats, printMe, metric)
		if len(vals) > 0:
			neighInds += vals
	#end = time.clock();
	#print end - start;
	return neighInds
 
 #Returns an array of arrays. Each inner array is a list of peptide indices that are in the same cluster. A cluster is defined as the largest list of molecules for which each molecule in the list is neighbors either directly or indirectly (neighbors of neighbors of neighbors etc...) neighbors with each other molecule. ind is the atom to check the cluster of, cutoff is the minimum distance that defines neighbors, peplist is a list of all atoms in the simulation, potentialInds is the list of indices that could possibly be neighbors with ind, ats is the number of atoms per peptide (this is important for dividing up pepList and must be constant for all peptides in pepList), and printMe is a boolean that will cause the immediate neighbors of ind to be printed if it is true (more for debuging and checking things).
def getClustTest(ind, cutoff, pepList, potentialInds, ats, printMe):
	#start = time.clock()
	neighInds = getNeighTest(ind, cutoff, pepList, potentialInds, ats)
	if printMe:
		print("Neighbors of " + str(ind) + " are found to have indices of: ")
		print(neighInds)

	for neighInd in neighInds:
		potentialInds.remove(neighInd)
	if ind == 3:
		2
	for neighInd in neighInds:
		vals = getClustTest(neighInd, cutoff, pepList, potentialInds, ats, printMe)
		if len(vals) > 0:
			neighInds += vals
	#end = time.clock();
	#print end - start;
	return neighInds
 
def getNeighTest(ind, cutoff, peplist, potentialInds, ats):
	ret = []

	cutsq = cutoff
#	support = '#include <math.h>'
#
#	code = """
#
#	int i, j;
#	return_val = 0;
#	for(i=0; i<Npep1[0]/3; i++){
#		for(j=0; j<Npep2[0]/3; j++){
#			if ((pow(pep1[3*i]-pep2[3*j],2) + pow(pep1[3*i+1]-pep2[3*j+1],2) + pow(pep1[3*i+2]-pep2[3*j+2],2)) < cutsq){
#				return_val = 1;
#				break;
#			}
#		}
#		if(return_val == 1)
#			break;
#	}
#			"""
#	pep1 = peplist[ind*3*ats:(ind+1)*3*ats] #Assumes all peptides have ats atoms in them. The 3 is for 3 coords in 3 dimensions.
#	for i in range(len(potentialInds)):
#		pep2 = peplist[potentialInds[i]*3*ats:(potentialInds[i]+1)*3*ats]
#		test = weave.inline(code,['pep1', 'pep2', 'cutsq'], support_code = support, libraries = ['m'])
#		if test == 1:
#			ret.append(potentialInds[i])
	for p in potentialInds:
		d = abs(peplist[ind]-peplist[p])
		if d < cutsq:
			ret.append(p)
 
	return ret

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
		clusts = getClust(init,cutoff,peps,pots,ats,False) + [init]

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
		
def writeClustLibrary(t,xtc,tpr,outgro,cutoff,ats,molStruct,resInfo,rm=True):
	#write out a library of .gro files of clusters
	(peps,box_length,qbool) = getPosQ(t,xtc,tpr,outgro,resInfo,ats)
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
		qboolc = np.zeros(len(clusts))
		curr = 0
		for clust in clusts:
			pepList[curr*ats*3:(curr+1)*ats*3]=peps[clust*ats*3:(clust+1)*ats*3]

   			qboolc[curr] = qbool[clust]
			curr+=1
		pepList = fixPBC(pepList,box_length,ats,cutoff)
		mer = len(clusts)
		if mer in clustinds:
			clustinds[mer]+=1
		else:
			clustinds[mer] = 1
		
		
		writeGro(clustinds[mer],pepList,box_length,qboolc)
  


#use the aromarecluster file which has group index and peptide index in it, to put together a list of clusters as defined in that file and then write those clusters out separately like writeClustLibrary
def writeSubClusterLibrary(aromaFile,t,xtc,tpr,outgro,resInfo,ats,cutoff):
	a = open(aromaFile,'r')
	lines = a.readline()
	a.close()
	clusters = dict()
	for line in lines:
		spline = line.split()
		ind = int(spline[1]) 
		pepind = int(spline[30])
		if clusters.has_key(ind):
			j[ind].append(pepind)
		else:
			clusters[ind] = [pepind]

	(peps,box_length,qbool) = getPosQ(t,xtc,tpr,outgro,resInfo,ats)
	for key in clusters.keys():
		clustinds = clusters[key]
		pepList = np.zeros(3*len(clustinds))
		for c in range(len(clustinds)):
			pepList[3*c:(3*c+3)] = peps[3*clustinds[c]:(3*clustinds[c]+3)]
		pepList = fixPBC(pepList,box_length,ats,cutoff) 
		
		
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
	print ndistrib
	for i in range(len(ndistrib)):
		nmols += (i+1)*ndistrib[i]
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
	print "nmols is ", nmols
	for k in range(len(ndistrib)-1,-1,-1):
		n = ndistrib[k]
		for i in range(int(n)):
			#randomly choose library file
			#call gromacs to insert it into box
			nfiles = libraryns[k]
			filei = random.randint(1,nfiles)
			fnamei = str(k+1)+"mer_"+str(filei)+".gro"
			f = open("confs.gro","w")
			if ( (k==4) and (i==0)):
				#call editconf
				command = ['editconf','-f',fnamei,'-o','init_'+str(ind)+'.gro','-bt','triclinic','-box',str(boxsize[0]),str(boxsize[1]),str(boxsize[2])]
				p = subprocess.call(command)
				#(outdata,errdata) = p.communicate()
				ind+=1
				f.write("configuration "+fnamei)
			else:
				tries = 200
				trying=-1
				while trying==-1:
					#call genconf
					f.write("configuration "+fnamei)
					#print "trying ",tries
					#f = open("temp.txt","w")
					command = ['genbox','-cp','init_'+str(ind-1)+'.gro','-o','init_'+str(ind)+'.gro','-ci',fnamei,'-nmol','1','-try',str(tries)]	
					
					p = subprocess.call(command)
					#(outdata,errdata) = p.communicate()
					#p.wait()
					#f.close()
					#f2 = open("temp.txt","r")
					#outdata = f2.read()
					#trying = outdata.find('Added 1 molecules')
					#f2.close()
					#tries+=20
					ind+=1
					trying = 1
			
			#print "START: ",outdata
			#print "END: ",errdata
	f.close()
	for j in range(ind-1):
		os.system('rm init_'+str(j)+'.gro')
	#cmd = 'export GMX_MAXBACKUP=99'
	#e2 = subprocess.Popen(cmd,stdout=subprocess.PIPE)
	#e2.wait()
	#print "under construction"
def aromaCross(arom,molInds,posList,boxL):
	
	b11 = molInds[arom[0]*3:arom[0]*3+3]
	b12 = molInds[arom[1]*3:arom[1]*3+3]
	b13 = molInds[arom[2]*3:arom[2]*3+3]
	b11p = posList[b11]
	b12p = posList[b12]
	b13p = posList[b13]
	db1 = b13p - b12p
	db2 = b12p - b11p
	for i in range(3):
		db1[i] = db1[i] - boxL[i]*round(db1[i]/boxL[i])
		db2[i] = db2[i] - boxL[i]*round(db2[i]/boxL[i])
	aromaPlane = np.cross(db1,db2)
	return aromaPlane

def aromaticData(posList,ats,beadList,beadMasses,aromBeads,aromMass,boxLs):
	#characterize the COM location and orientation of the aromatic rings in the backbone (beads 10, 11, 12, 13, 14, 15, 16, 17, 21, 20, 19 for DFAG)
	#right now assumes molecules have been made whole and that there are three aromatic rings in the backbone that we are being given
	#orient is the two vectors pointing from the first to the second aromatic ring and from the second to the third
	#also return the three vectors of the planes of the aromatic rings as aromaPlanes
	N = int(len(posList)/3) #number of atoms
	nmols = N/ats
	comsPBC = np.zeros([nmols,3])
	orients = np.zeros([nmols,6])
	aromaPlanes = np.zeros([nmols,9])
	#loop over each group of aromatic rings and calculate COM and orient
	#molinds = indices of each molecule: 0:3*ats-1, 3*ats:2(3*ats)-1,...,(N-1)(3*ats):N(3*ats)-1
	#indices of each aromatic group: molinds[3*beadList:3*beadList+2]
	for i in range(nmols):
		molInds = range(i*(3*ats),(i+1)*(3*ats))
		#aromInds = molinds[3*beadList:3*beadList+2]
		#endInds = molinds[3*endList:3*endList+2]
		comPBC = np.array([0.0,0.0,0.0])
		M = 0.0
		j = 0
		beadInd0 = molInds[beadList[0]*3:beadList[0]*3+3]
		b0 = posList[beadInd0]
		#print "b0 ",b0
		m0 = beadMasses[0]
		for bead in beadList:
			beadInds = molInds[bead*3:bead*3+3]
			b = posList[beadInds]
			m = beadMasses[j]
			j+=1
			db = b - b0
			db = db - boxLs * np.round(db/boxLs)
			comPBC += m*db
			M+=m
		
		comPBC = comPBC/M + b0
		arom1 = aromBeads[0:3]
		arom2 = aromBeads[3:6]
		arom3 = aromBeads[6:9]
		coma1 = np.array([0.0,0.0,0.0])
		coma2 = np.array([0.0,0.0,0.0])
		coma3 = np.array([0.0,0.0,0.0])
		coma1 = aromaCOM(arom1,posList,aromMass[0:3],molInds,boxLs)
		coma2 = aromaCOM(arom2,posList,aromMass[3:6],molInds,boxLs)
		coma3 = aromaCOM(arom3,posList,aromMass[6:9],molInds,boxLs)
		aromaPlane1 = aromaCross(arom1,molInds,posList,boxLs)
		aromaPlane2 = aromaCross(arom2,molInds,posList,boxLs)
		aromaPlane3 = aromaCross(arom3,molInds,posList,boxLs)
		#normalize vectors
		aromaPlane1 = aromaPlane1/np.sqrt(sum(aromaPlane1**2))
		aromaPlane2 = aromaPlane2/np.sqrt(sum(aromaPlane2**2))
		aromaPlane3 = aromaPlane3/np.sqrt(sum(aromaPlane3**2))
		comaa = coma3 - coma2
		comaa = comaa - boxLs*np.round(comaa/boxLs)
		comab = coma2 - coma1
		comab = comab - boxLs*np.round(comab/boxLs)
		orient = np.concatenate((comaa,comab),axis=0)
		#put comPBC inside box if it isn't
		#print "comPBC is ",comPBC
		comPBC = inBox(comPBC,boxLs)
		#print "comPBC is now ",comPBC
		comsPBC[i,:] = comPBC
		orients[i,:] = orient
		aromaPlanes[i,:] = np.concatenate((aromaPlane1,aromaPlane2,aromaPlane3),axis=0)
	return (comsPBC,orients,aromaPlanes)

def inBox(r,boxLs):
	for i in range(3):
		r[i] = r[i] - boxLs[i]*((int(r[i])/int(boxLs[i])))
	return r

def aromaCOM(arom,posList,aromMass,molInds,boxL):
	b0 = molInds[arom[0]*3:arom[0]*3+3]
	b0p = posList[b0]
	b1 = molInds[arom[1]*3:arom[1]*3+3]
	b1p = posList[b1]
	b2 = molInds[arom[2]*3:arom[2]*3+3]
	b2p = posList[b2]
	db1 = b1p - b0p
	db2 = b2p - b0p 
	#print db1
	#print db2
	#print boxL
	db1 = db1 - boxL*np.round(db1/boxL)
	db2 = db2 - boxL*np.round(db2/boxL)
	coma = (aromMass[1]*db1 + aromMass[2]*db2)/(aromMass[0]+aromMass[1]+aromMass[2]) + b0p
	return coma

def aromaticSnap(infname,cutoff,ats,beadList,beadMasses,aromBeads,aromMass):
	(peps,box_length) = getPosGro(infname)
	pots = range(len(peps)/3/ats)
	inds = np.zeros(len(peps)/3/ats)
	ind = 1
	achars = np.empty((0,30),float)
	cNum = 0
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
		(coms,orients,aromaPlanes,aromaCOMs) = aromaticDataC(pepList,ats,beadList,beadMasses,aromBeads,aromMass,box_length)
	
		cNums = cNum*np.ones([len(coms),1])
		cSizes = len(clusts)*np.ones([len(coms),1])
		clustList = np.array(clusts)
		clustList = np.reshape(clustList,[len(coms),1])
		try:
			adata = np.concatenate((cNums,cSizes,coms,orients,aromaPlanes,aromaCOMs,clustList),axis=1)
			achars = np.concatenate((achars,adata),axis=0)
		except:
			np.shape(clustList)      
		cNum+=1
	
	return achars
 
def getCOM(poslist,masslist):
    #return com coordinates of a single molecule
    com = np.zeros(3)
    N = len(poslist)
    support = '#include <math.h>'
    code = """
    double M = 0.0;
    double rx = 0.0;
    double ry = 0.0;
    double rz = 0.0;
    for (int i = 0; i < N/3; i++){
            rx += double(masslist[i])*double(poslist[3*i]);
            ry += double(masslist[i])*double(poslist[3*i+1]);
            rz += double(masslist[i])*double(poslist[3*i+2]);
            M += double(masslist[i]);
    }
    rx /= M;
    ry /= M;
    rz /= M;
    com[0] = rx;
    com[1] = ry;
    com[2] = rz;
    """
    weave.inline(code,['poslist','N','com','masslist'],support_code = support,libraries=['m'])
    return com    
    
    
def getCOMpy(poslist,masslist):
    #return com coordinates of a single molecule, written in python
    N = len(poslist)
    com = np.zeros(3)
    rx = 0.
    ry = 0.
    rz = 0.
    for i in range(N/3):
        rx += masslist[i]*poslist[3*i]
        ry += masslist[i]*poslist[3*i+1]
        rz += masslist[i]*poslist[3*i+2]
    M = np.sum(masslist)
    rx /= M
    ry /= M
    rz /= M
    com[0] = rx
    com[1] = ry
    com[2] = rz
    return com
    
def getCOMnumpy(poslist,masslist):
    #return com coordinates of a single molecule, written in python using numpy
    poslist3 = poslist.reshape([len(poslist)/3,3])
    com = np.dot(masslist,poslist3)/masslist.sum()
    return com
    
def getCOMs(poslist,masslist,ats):
    #return the center of mass of each molecule in a position list of peptides
    N = len(poslist)/3/ats #total number of molecules
    comlocs = np.zeros(3*N)
    for i in range(N):
        Rcom = getCOM(poslist[3*ats*i:3*ats*i+3*ats],masslist)
        comlocs[3*i:3*i+3] = Rcom
    return comlocs
'''    
def getCOMsWRONG(poslist,masslist,ats):
    #return the center of mass of each molecule in a position list of peptides
    N = len(poslist)/3/ats #total number of molecules
    comlocs = np.zeros(3*N)
    for i in range(N):
        Rcom = getCOMWRONG(poslist[3*ats*i:3*ats*i+3*ats],masslist)
        comlocs[3*i:3*i+3] = Rcom
    return comlocs
    
def getCOMspy(poslist,masslist,ats):
    #return the center of mass of each molecule in a position list of peptides
    N = len(poslist)/3/ats #total number of molecules
    comlocs = np.zeros(3*N)
    for i in range(N):
        Rcom = getCOMpy(poslist[3*ats*i:3*ats*i+3*ats],masslist)
        comlocs[3*i:3*i+3] = Rcom
    return comlocs
    
def getCOMsnumpy(poslist,masslist,ats):
    #return the center of mass of each molecule in a position list of peptides
    N = len(poslist)/3/ats #total number of molecules
    comlocs = np.zeros(3*N)
    for i in range(N):
        Rcom = getCOMnumpy(poslist[3*ats*i:3*ats*i+3*ats],masslist)
        comlocs[3*i:3*i+3] = Rcom
    return comlocs
'''
    
def corrDim(comsnap,emax,estep):
    #calculate C(eps) from 0 to emax with steps of estep, where C(eps) is 
    #the correlation sum on the set of center-of-mass distances given as comsnap
    
    distsq = spsd.pdist(comsnap,'sqeuclidean')

    epss = np.arange(estep,emax+estep,estep) 
    ce = np.zeros(len(epss))

    support = '#include <math.h>'
    code = """
    for (int i = 0; i < Ndistsq[0]; i++){
        for (int k = 0; k < Nepss[0]; k++){
            if (distsq[i] <= epss[k]*epss[k]){
                ce[k]++;
                
            }
        }    
    }
    """
    weave.inline(code,['emax','estep','distsq','ce','epss'],support_code=support,libraries=['m'])
    N = len(comsnap)
    ce = ce/(N*(N-1))
#    for i in range(Ndsq):
#        
#        for k in range(len(epss)):
#            if distsq[i] < epss[k]*epss[k]:
#                ce[k]+=1

    return (epss,ce)
    
def corrDimP(comsnap,emax,estep):
    #calculate C(eps) from 0 to emax with steps of estep, where C(eps) is 
    #the correlation sum on the set of center-of-mass distances given as comsnap
    
    distsq = spsd.pdist(comsnap,'sqeuclidean')

    epss = np.arange(estep,emax+estep,estep) 
    ce = np.zeros(len(epss))

    #support = '#include <math.h>'
    #code = 
    """
    for (int i = 0; i < Ndistsq[0]; i++){
        for (int k = 0; k < Nepss[0]; k++){
            if (distsq[i] <= epss[k]*epss[k]){
                ce[k]++;
                
            }
        }    
    }
    """
   # weave.inline(code,['emax','estep','distsq','ce','epss'],support_code=support,libraries=['m'])
    N = len(comsnap)
    ce = ce/(N*(N-1))
    for i in range(len(distsq)):
        
        for k in range(len(epss)):
            if distsq[i] <= epss[k]*epss[k]:
                ce[k]+=1

    return (epss,ce)
 
def corrcalc(t,xtc,tpr,outgro,ats,emax,estep,masslist,fnbase,rm=True):
    #first calculate and write out the C(e) function for a full snapshot
    (peps,box_length) = getPosB(t,xtc,tpr,outgro)

    if rm:
        os.system('rm '+outgro)
    comsnap = getCOMs(peps,masslist,ats) #get the center of mass locations of the molecules in the box
    comsnapr = comsnap.reshape([len(comsnap)/3,3])
    (epsnap,cdsnap) = corrDim(comsnapr,emax,estep) #get a matrix of e in row 1 and C(e) in row 2
    f = open(fnbase+'_'+str(t/1000)+'_snap.dat','w')
    for i in range(len(cdsnap)):
        f.write('{0}\t{1}\n'.format(epsnap[i],cdsnap[i]))
    f.close()

def corrcalcClust(t,xtc,tpr,outgro,ats,emax,estep,masslist,fnbase,cutoff,rm=True):
    #then calculate and write out the C(e) function for each cluster
    (peps,box_length) = getPosB(t,xtc,tpr,outgro)
	#print box_length
    if rm:
        os.system('rm '+outgro)
    pots = range(len(peps)/3/ats)
    ind = 0
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
        pepList = fixPBC(pepList,box_length,ats,cutoff)
        #if len(pepList)/3/ats == 1:
         #   print "help"
        pepcoms = getCOMs(pepList,masslist,ats)
        pepcomsr = pepcoms.reshape([len(pepcoms)/3,3])
        (epc,cdc) = corrDim(pepcomsr,emax,estep)
        f = open(fnbase+'_'+str(t/1000)+'_'+str(len(pepList)/3/ats)+'_c_'+str(ind)+'.dat','w')
        for i in range(len(cdc)):
            f.write('{0}\t{1}\n'.format(epc[i],cdc[i]))
        f.close()
        ind+=1
    
def testcorrcalc(infile,outfile,emax,estep):
    #make sure corrcalc does the same thing that Jiang's test code does
    f = open(infile)
    lines = f.readlines()
    f.close()
    comsnap = np.zeros([len(lines),3])
    l = 0
    for line in lines:
        spline = line.split()
        comsnap[l,0] = float(spline[0])
        comsnap[l,1] = float(spline[1])
        comsnap[l,2] = float(spline[2])
        l+=1
    (epstest,cdtest) = corrDim(comsnap,emax,estep)
    o = open(outfile,'w')
    for i in range(len(cdtest)):
        o.write('{0}\t{1}\n'.format(epstest[i],cdtest[i]))
    o.close()

def aromaticClust(t,xtc,tpr,outgro,cutoff,ats,beadList,beadMasses,aromBeads,aromMass,rm=True):
	(peps,box_length) = getPosB(t,xtc,tpr,outgro)
	#print box_length
	if rm:
		os.system('rm '+outgro)
	pots = range(len(peps)/3/ats)
	inds = np.zeros(len(peps)/3/ats)
	ind = 1
	achars = np.empty((0,30),float)
	cNum = 0
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
		(coms,orients,aromaPlanes,aromaCOMs) = aromaticDataC(pepList,ats,beadList,beadMasses,aromBeads,aromMass,box_length)
	
		cNums = cNum*np.ones([len(coms),1])
		cSizes = len(clusts)*np.ones([len(coms),1])
		clustList = np.array(clusts)
		clustList = np.reshape(clustList,[len(coms),1])
		try:
			adata = np.concatenate((cNums,cSizes,coms,orients,aromaPlanes,aromaCOMs,clustList),axis=1)
			achars = np.concatenate((achars,adata),axis=0)
		except:
			np.shape(clustList)      
		cNum+=1
	
	return achars

def aromaticDataC(posList,ats,beadList,beadMasses,aromBeads,aromMass,boxLs):
	coms = aromaticCOMC(posList,ats,beadList,beadMasses,boxLs)
	(orients,aromaPlanes,aromaCOMs) = aromaticPlaneC(posList,ats,aromBeads,aromMass,boxLs)
	return (coms,orients,aromaPlanes,aromaCOMs)

def aromaticCOMC(posList,ats,beadList,beadMasses,boxLs):
	#characterize the COM location and orientation of the aromatic rings in the backbone (beads 10, 11, 12, 13, 14, 15, 16, 17, 21, 20, 19 for DFAG)
	#right now assumes molecules have been made whole
	N = int(len(posList)/3) #number of atoms
	nmols = N/ats
	coms = np.zeros(3*nmols)
	#loop over each group of aromatic rings and calculate COM and orient
	#molinds = indices of each molecule: 0:3*ats-1, 3*ats:2(3*ats)-1,...,(N-1)(3*ats):N(3*ats)-1
	#indices of each aromatic group: molinds[3*beadList:3*beadList+2]
	support = '#include <math.h>'
	code = """
	int * molInds;
	molInds = (int *) malloc(sizeof(int)*3*ats);
	double com [3];
	com[0] = 0.0;
	com[1] = 0.0;
	com[2] = 0.0;


	double M,bx,by,bz,m,m0,b0x,b0y,b0z,dbx,dby,dbz;
	int bead,beadInd0x,beadInd0y,beadInd0z,beadIndx,beadIndy,beadIndz,end1x,end1y,end1z,end2x,end2y,end2z;
	for (int i = 0; i < nmols; i++){
		for (int k = 0; k < 3*ats; k++){
			molInds[k] = i*3*ats+k;
		}
		for (int k = 0; k < 3; k++){
			com[k] = 0.0;
		}
		M = 0.0;
		beadInd0x = molInds[beadList[0]*3];
		beadInd0y = molInds[beadList[0]*3+1];
		beadInd0z = molInds[beadList[0]*3+2];
		b0x = posList[beadInd0x];
		b0y = posList[beadInd0y];
		b0z = posList[beadInd0z];
		for (int j = 0; j < NbeadList[0]; j++){
			bead = beadList[j];
			beadIndx = molInds[bead*3];
			beadIndy = molInds[bead*3+1];
			beadIndz = molInds[bead*3+2];
			bx = posList[beadIndx];
			by = posList[beadIndy];
			bz = posList[beadIndz];
			dbx = bx - b0x;
			dby = by - b0y;
			dbz = bz - b0z;
			dbx = dbx - double(boxLs[0])*round(dbx/double(boxLs[0]));
			dby = dby - double(boxLs[1])*round(dby/double(boxLs[1]));
			dbz = dbz - double(boxLs[2])*round(dbz/double(boxLs[2]));
			m = beadMasses[j];
			com[0] += m*dbx;
			com[1] += m*dby;
			com[2] += m*dbz;
			M+=m;
		}
		com[0] = com[0]/M+b0x;
		com[1] = com[1]/M+b0y;
		com[2] = com[2]/M+b0z;
		coms[3*i] = com[0];
		coms[3*i+1] = com[1];
		coms[3*i+2] = com[2];
	}
	free(molInds);
	"""
	weave.inline(code,['posList','ats','nmols','beadList','beadMasses','coms','boxLs'],support_code = support,libraries=['m'])
	coms = np.reshape(coms,[nmols,3])
	return coms

def aromaticPlaneC(posList,ats,aromBeads,aromMass,boxLs):
    N = int(len(posList)/3)
    nmols = N/ats
    orients = np.zeros(6*nmols)
    aromaPlanes = np.zeros(9*nmols)
    aromaCOMs = np.zeros(9*nmols)
    support = '#include <math.h>'
    code = """
    int * molInds;
    molInds = (int *) malloc(sizeof(int)*3*ats);
    double orient [6];
    for (int i = 0; i < 6; i++){
        orient[i] = 0.0;
    }
    
    double comas [9];
    int b0 [9];
    int b1 [9];
    int b2 [9];
    double b0p [9];
    double b1p [9];
    double b2p [9];
    double db1 [9];
    double db2 [9];
    double dc1 [9];
    double dc2 [9];
    double norm;
    for (int i = 0; i < 9; i++){
        b0[i] = 0;
        b1[i] = 0;
        b2[i] = 0;
        b0p[i] = 0.0;
        b1p[i] = 0.0;
        b2p[i] = 0.0;
        db1[i] = 0.0;
        db2[i] = 0.0;
        dc1[i] = 0.0;
        dc2[i] = 0.0;
        comas[i] = 0.0;
    }
    
    
    for (int i = 0; i < nmols; i++){
        for (int k = 0; k < 3*ats; k++){
            molInds[k] = i*3*int(ats)+k;
        }
        for (int c = 0; c < 3; c++){
        for (int j = 0; j < 3; j++){
            b0[j+3*c] = molInds[int(aromBeads[0+3*c])*3+j];
            b1[j+3*c] = molInds[int(aromBeads[1+3*c])*3+j];
            b2[j+3*c] = molInds[int(aromBeads[2+3*c])*3+j];
            b0p[j+3*c] = posList[b0[j+3*c]];
            b1p[j+3*c] = posList[b1[j+3*c]];
            b2p[j+3*c] = posList[b2[j+3*c]];
            db1[j+3*c] = b1p[j+3*c]-b0p[j+3*c];
            db2[j+3*c] = b2p[j+3*c]-b0p[j+3*c];
            db1[j+3*c] = db1[j+3*c] - double(boxLs[j])*round(db1[j+3*c]/double(boxLs[j]));
            db2[j+3*c] = db2[j+3*c] - double(boxLs[j])*round(db2[j+3*c]/double(boxLs[j]));
            dc1[j+3*c] = b2p[j+3*c] - b1p[j+3*c];
            dc2[j+3*c] = b1p[j+3*c] - b0p[j+3*c];
            dc1[j+3*c] = dc1[j+3*c] - double(boxLs[j])*round(dc1[j+3*c]/double(boxLs[j]));
            dc2[j+3*c] = dc2[j+3*c] - double(boxLs[j])*round(dc2[j+3*c]/double(boxLs[j]));
            comas[j+3*c] = (aromMass[1+3*c]*db1[j+3*c]+aromMass[2+3*c]*db2[j+3*c])/(aromMass[0+3*c]+aromMass[1+3*c]+aromMass[2+3*c])+b0p[j+3*c];
	    aromaCOMs[9*i+j+3*c] = comas[j+3*c];
            
        }
            aromaPlanes[9*i+3*c] = dc1[3*c+1]*dc2[3*c+2] - dc1[3*c+2]*dc2[3*c+1];
            aromaPlanes[9*i+3*c+1] = -dc1[3*c]*dc2[3*c+2] + dc1[3*c+2]*dc2[3*c];
            aromaPlanes[9*i+3*c+2] = dc1[3*c]*dc2[3*c+1] - dc1[3*c+1]*dc2[3*c];
            norm = aromaPlanes[9*i+3*c]*aromaPlanes[9*i+3*c]+aromaPlanes[9*i+3*c+1]*aromaPlanes[9*i+3*c+1]+aromaPlanes[9*i+3*c+2]*aromaPlanes[9*i+3*c+2];
            norm = sqrt(norm);
            aromaPlanes[9*i+3*c] = aromaPlanes[9*i+3*c]/norm;
            aromaPlanes[9*i+3*c+1] = aromaPlanes[9*i+3*c+1]/norm;
            aromaPlanes[9*i+3*c+2] = aromaPlanes[9*i+3*c+2]/norm;
        }
        orient[3] = comas[6] - comas[3] - double(boxLs[0])*round((comas[6]-comas[3])/double(boxLs[0]));
        orient[4] = comas[7] - comas[4] - double(boxLs[1])*round((comas[7]-comas[4])/double(boxLs[1]));
        orient[5] = comas[8] - comas[5] - double(boxLs[2])*round((comas[8]-comas[5])/double(boxLs[2]));
        orient[0] = comas[3] - comas[0] - double(boxLs[0])*round((comas[3]-comas[0])/double(boxLs[0]));
        orient[1] = comas[4] - comas[1] - double(boxLs[1])*round((comas[4]-comas[1])/double(boxLs[1]));
        orient[2] = comas[5] - comas[2] - double(boxLs[2])*round((comas[5]-comas[2])/double(boxLs[2]));
        for (int o = 0; o < 6; o++){
            orients[6*i+o] = orient[o];
        }

    }

    """
    weave.inline(code,['posList','ats','nmols','aromBeads','aromMass','orients','boxLs','aromaPlanes','aromaCOMs'],support_code = support,libraries=['m'])
    orients = np.reshape(orients,[nmols,6])
    aromaPlanes = np.reshape(aromaPlanes,[nmols,9])
    aromaCOMs = np.reshape(aromaCOMs,[nmols,9])
    return (orients,aromaPlanes,aromaCOMs)

def betaCharacter(posList,ats,cutoff,bbs,bblist):
	#characterize the "beta-sheetness" of a given cluster based on how many consecutive beads have bonds
	N = int(len(posList)/3)
	cutoff2 = cutoff*cutoff
	#aMat = np.zeros([N,N]) #adjacency matrix for cluster
	betaBonds = np.zeros(bbs)
	#print "len(betabonds) = ",len(betaBonds)
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
	for (int k = 0; k < (N/ats)*bbs; k++){
        free(aMat[k]);
    	}
	free(aMat);
	"""
	weave.inline(code,['posList','N','ats','cutoff2','betaBonds','bbs','bblist'],support_code = support,libraries=['m'])
	#print "len(betabonds) is ",len(betaBonds)
    	return betaBonds

def betaClust(t,grobase,cutoff,ats,bbs,bblist,ind,rm=True):
	#return a set of counts (normalized by the total number) of "beta-bonds" or appearing contiguous segments of bonds between atoms
	#like a beta-ladder kind of thing for each cluster
	(peps,box_length) = getPosB2(ind,grobase)
	#print box_length
	if rm:
		os.system('rm '+outgro)
	pots = range(len(peps)/3/ats)
	inds = np.zeros(len(peps)/3/ats)
	ind = 1
	betachars = np.empty((0,bbs),float)
        ms = np.array([])	
	while len(pots) > 0:
		
		init = pots[0]
		pots.remove(init)
		clusts = getClust(init,cutoff,peps,pots,ats,False) + [init]
		#clusts is a list of peptides that are found in the cluster
		#each index in clusts corresponds to the indices index*ats*3:(index+1)*ats*3 in peps
		pepList = np.zeros(len(clusts)*ats*3)
		curr = 0
		mass = float(len(clusts))
		for clust in clusts:
			pepList[curr*ats*3:(curr+1)*ats*3]=peps[clust*ats*3:(clust+1)*ats*3]
			curr+=1
		betabonds = betaCharacter(pepList,ats,cutoff,bbs,bblist)
		#print("found betabonds: ",betabonds)
		betachars = np.append(betachars,np.array([betabonds]),axis=0)
		#print(betachars)
		ms = np.append(ms,mass)
	#print(len(betachars))
	return (ms,betachars)
		

def fixPBC(peps,box,ats,cutoff):
	#return positions fixed across PBCs for calculation of structural metrics like Rh and Rg
	#create the list of fixed positions
	fixedXYZ = peps.copy()
	potInds = range(1,len(peps)/(ats*3))
	#the first ats*3 coordinates are the coordinates of the first atom
	fixedXYZ[0:3*ats] = fixCoords(fixedXYZ[0:3*ats].copy(),fixedXYZ[0:3].copy(),box)
	correctInds = [0]
	while len(correctInds) > 0:
		atom = correctInds.pop()
		neighs = getNeigh(atom,cutoff,peps,potInds,ats)
		for n in neighs:

			potInds.remove(n)
			correctInds.append(n)
			fixedXYZ[3*ats*n:3*ats*(n+1)] = fixCoords(fixedXYZ[3*ats*n:3*ats*(n+1)].copy(),fixedXYZ[3*atom*ats:3*atom*ats+3].copy(),box)
	return fixedXYZ

def fixCoords(pos,posinit,box):
	#fix all coords based on the initial coordinate and the periodic boundary conditions
	for i in range(len(pos)/3):
		dr = pos[3*i:3*i+3] - posinit
		dr = dr - box*np.round(dr/box)
		pos[3*i:3*i+3] = dr + posinit
	return pos
'''
def e2edist(pepList):
    #find the largest distance between any two beads in the cluster
    

def end2end(t,xtc,tpr,outgro,cutoff,ats,rm=True):
    #run through all clusters and return cluster size and largest distance between
    #any two beads in the cluster
    (peps,box_length) = getPosB(t,xtc,ptr,outgro)
    if rm:
        os.system('rm '+outgro)
    pots = range(len(peps)/3/ats)
    inds = np.zeros(len(peps)/3/ats)
    ind = 1
    e2es = np.array([])
    ms = np.array([])
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
        pepList = fixPBC(pepList,box_length,ats,cutoff)
        ee = e2edist(pepList)
        mass = float(len(clust))
        ms = np.append(ms,mass)
        e2es = np.append(e2es,ee)
    return (ms,e2es)
'''

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
		pepList = fixPBC(pepList,box_length,ats,cutoff)
		#writeGro('mol_'+str(ind)+'.gro',pepList,box_length)
		gyrationTensor = gyrTens(pepList,box_length)
		eigstuff = np.linalg.eig(gyrationTensor)
		ind+=1
		eigMorph = np.zeros(5)
		eigval = np.sort(eigstuff[0])
		eigMorph[0] = eigval[0]
		eigMorph[1] = eigval[1]
		eigMorph[2] = eigval[2]
		eigMorph[3] = eigMorph[0]+eigMorph[1]+eigMorph[2] #Rg^2
		eigMorph[4] = 1.5*eigMorph[2]-0.5*eigMorph[3]
		
		eigvals = np.append(eigvals,eigMorph)
		eigOrder = np.argsort(eigstuff[0])
		eigvec = eigstuff[1]
		eigvect = np.append(eigvec[:,eigOrder[2]],eigvec[:,eigOrder[1]])
		eigvect = np.append(eigvect,eigvec[:,eigOrder[0]])
		eigvecs = np.append(eigvecs,eigvect)
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
            V = V - double(boxlx)*round(double(V)/double(boxlx));
            U = posList[3*R+y] - posList[3*S+y];
            U = U - double(boxly)*round(double(U)/double(boxly));
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
  
if __name__ == '__main__':
#    comsnap = np.array([1.,0.,0.,2.,0.,0.,3.,0.,0.,4.,0.,0.])
#    emax = 3.5
#    estep = 0.5
#    ce = corrDim(comsnap,emax,estep)
#    print ce
    emax = 42
    estep = 0.1
    folder = '/home/rachael/cluster/data/mansbac2/coarsegraining/MD/DFAG/378mols_run4/4_md/'
    ti = 400000
    tf = 400000
    dt = 400000
    xtc = folder+'sizedata/md_noW.xtc'
    tpr = folder+'sizedata/md_dummy.tpr'
    ats = 29
    fnbase = folder+'corrdimdata/mcdim'
    masslist = np.ones(ats)*45.
    m72 = [0,1,2,6,7,21,22,23,27,28]
    masslist[m72] = 72.
    outgro = 'temp.gro'
    cutoff=0.5
    for t in range(ti,tf+dt,dt):
        print t
        corrcalc(t,xtc,tpr,outgro,ats,emax,estep,masslist,fnbase,rm=True)
        #corrcalcClust(t,xtc,tpr,outgro,ats,emax,estep,masslist,fnbase,cutoff,rm=True)

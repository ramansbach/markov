import numpy as np, imp,os
mk = imp.load_source('mk','/home/rachael/Analysis_and_run_code/coarsegraining/markov/markov.py')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('text',usetex=False)
import subprocess,timeit
#use visualization to check
#data: cluster size + charged/uncharged molecules + COM + distance of ends from COM (distribution of cluster end distance from COM?)
#first get 
#markov fixPBC will make whole
#need a new readGro function that returns a list of 1s and 0s for charged and uncharged molecules based on the residue name

#Takes a .gro file and returns a 1d list of positions, box lengths, and a 1d list of whether the molecule does
 #or does not contain a particular residue (used to tell if it's charged or uncharged)
#resInfo is a tuple ('resname',pos), where position is the index within the molecule of the residue
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


#returns the following data
#data: cluster size + charged/uncharged molecules + COM + distance of ends from COM (distribution of cluster end distance from COM?)
#resInfo is ('PAQ',0) for DFAG,masses is a list of the beads in .gro file order mass in amu (either 72 or 46)
#endInds are the indices in the .gro file containing the 4 charged beads (in DFAG, that is [0,1,27,28]
def clustStruct(t,xtc,tpr,outgro,cutoff,resInfo,ats,masses,endInds,bins,elongRat=8,elongMols=400,rm=True):
    (peps,box_length,qbool) = getPosQ(t,xtc,tpr,outgro,resInfo,ats)
    if rm:
            os.system('rm '+outgro)
    pots = range(len(peps)/3/ats)
    inds = np.zeros(len(peps)/3/ats)
    ind = 1
    hQs = np.zeros([len(pots),bins]) #each row is a histogram which is the count of charged molecules found within the spatial region from the center--may be the sum of slices for elongated molecules
    #each row signifies a different sized cluster
    hnQs = np.zeros([len(pots),bins]) #same as previous, except for noncharged
    rQs = np.zeros([len(pots),1])
    rQsnorm = np.zeros([len(pots),1])
    rnQs = np.zeros([len(pots),1])
    rnQsnorm = np.zeros([len(pots),1])
    while len(pots) > 0:
        init = pots[0]
        pots.remove(init)
        clust = mk.getClust(init,cutoff,peps,pots,ats,False)+[init]
        pepList = np.zeros(len(clust)*ats*3)
        curr = 0
        #clust is a list of monomers that are found in the cluster
        #each index in clusts corresponds to the indices index*ats*3:(index+1)*ats*3 in peps
        for c in clust:
            inds[c] = ind
            pepList[curr*ats*3:(curr+1)*ats*3]=peps[c*ats*3:(c+1)*ats*3] #1d list of positions in cluster
            curr+=1
        pepList = mk.fixPBC(pepList,box_length,ats,cutoff)
        com = clustCOM(pepList,masses)
        nQs = np.sum(qbool[clust]) #number of charged molecules in cluster
        nMs = len(clust) #number of molecules in cluster
        dEs = endDist(com,pepList,ats,endInds)
        #pep3 = np.reshape(pepList,[len(pepList)/3,3])
        #fig=plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(pep3[:,0],pep3[:,1],pep3[:,2],c='b',marker='.')
        #ax.scatter(com[0],com[1],com[2],c='r',marker='^')
        
#        fig2 = plt.figure()
 #       ax2 = fig2.add_subplot(111, projection='3d')
        qinds = getQinds(qbool[clust],len(endInds)*3)
        dEq = dEs[qinds==1]
        dEnq= dEs[qinds==0]
        dEqr = dEq.reshape(len(dEq)/3,3)
        dEnqr = dEnq.reshape(len(dEnq)/3,3)
        #dEq = dEs[qinds==1,:]
        #dEnq = dEs[qinds==0,:]
  #      ax2.scatter(dEqr[:,0],dEqr[:,1],dEqr[:,2],c='r',marker='o')
   #     ax2.scatter(dEnqr[:,0],dEnqr[:,1],dEnqr[:,2],c='b',marker='o')
    #    ax2.scatter(0,0,0,c='g',marker='^')
        gTens = mk.gyrTens(pepList,box_length)
        
        gStuff = np.linalg.eig(gTens)
        gVals = gStuff[0]
        gMoms = gStuff[1]
        rat = max(gVals)/min(gVals)
    
        if rat > elongRat and nMs > elongMols:
            #assume that we have a sinuous kind of molecule which requires special treatment
            h = elongHist(gVals,gMoms,dEs,bins)
        else:
            (rq,hq) = radialHist(dEq,bins)
            (rnq,hnq) = radialHist(dEnq,bins)
        hQs[len(clust),:] += hq
        rQs[len(clust)] += rq #track mean distance from COM
        rQsnorm[len(clust)]+=1 #normalization for mean (let's not do standard deviation yet)
        rnQs[len(clust)]+= rnq
        rnQsnorm[len(clust)]+=1
        hnQs[len(clust),:] += hnq
        #fig3 = plt.figure()
        #ax3 = fig3.add_subplot(111)
        #ax3.scatter(dEs[qinds==0,0],dEs[qinds==0,1],c='b',marker='o')
        #ax3.scatter(dEs[qinds==1,0],dEs[qinds==1,1],c='r',marker='o')
        #ax3.scatter(0,0,0,c='g',marker='^')
        #ax.set_xlabel('X Label')
        #ax.set_ylabel('Y Label')
        #ax.set_zlabel('Z Label')
        plt.show()
        #if ind > 4:
        #    pots = []
        
        ind+=1
        '''
        for i in range(len(rQs)):
        rQs[i] /= rQsnorm[i]
        rnQs[i] /= rnQsnorm[i]
        '''
    return (rQs,rQsnorm,rnQs,rnQsnorm,hQs,hnQs)

def elongHist(gVals,gMoms,pepList,bins):
    print("under construction")

def radialHist(beadDists,bins):
    #construct a histogram of the radial distances from the coms of the beads
    #compute radial distances from (x,y,z) points in relation to COM and track the max distance

    h = np.zeros([1,bins])
    rmax = 0
    rs = np.zeros([len(beadDists)/3])
    for i in range(len(beadDists)/3):
        r = np.sqrt(beadDists[3*i]*beadDists[3*i]+beadDists[3*i+1]*beadDists[3*i+1]+beadDists[3*i+2]*beadDists[3*i+2])
        rs[i] = r
        if r > rmax:
            rmax = r
    (h,bE) = np.histogram(rs,bins=bins,range=(0.0,rmax))
    
    #return histogram
    return (np.mean(rs),h)
        
def getQinds(qbool,l):
    #All this does is multiply qbool so that it has an entry for each end bead rather than just for the molecules
    qinds = np.zeros(l*len(qbool))
    for i in range(len(qbool)):
        if qbool[i]==1:
            qinds[l*i:l*i+l] = 1.0
    return qinds    

            
def clustCOM(pepList,masslist):
	'''
	[x1,y1,z1,x2,y2,z2,...,xn,yn,zn] pepList > len should be multiple of ats
	for x1,y1,z1, corresponding mass is masses[(ind % (3*ats))/3] <== integer division
	'''
	comx = 0
	comy = 0
	comz = 0
	M = np.sum(masslist)*len(pepList)/(3*len(masslist))
	for i in range(len(pepList)/3):
		xi = 3*i
		yi = 3*i+1
		zi = 3*i+2
		mi = masslist[i%len(masslist)]
		comx += mi*pepList[xi]
		comy += mi*pepList[yi]
		comz += mi*pepList[zi]
	comx /= M
	comy /= M
	comz /= M
	return (comx,comy,comz)
 

def endDist(com,pepList,ats,endInds):
    '''
    goes through the list, checking the (x-com,y-com,z-com) locations for the 4 charged beads in each peptide and returning in a 3d list
    need separate lists for charged beads and non-charged beads    
    '''
    molno = len(pepList)/(3*ats)
    dEs = np.zeros([molno*len(endInds)*3,1])
    for i in range(molno):
        k=0
        for j in endInds:
            ix = 3*ats*i+3*j
            iy = ix+1
            iz = iy+1
            x = pepList[ix] - com[0]
            y = pepList[iy] - com[1]
            z = pepList[iz] - com[2]
            dEs[3*len(endInds)*i+k] = x
            dEs[3*len(endInds)*i+k+1] = y
            dEs[3*len(endInds)*i+k+2] = z
            k+=3
    return dEs

if __name__ == "__main__":
    xtc = 'md_noW.xtc'
    tpr = 'md_dummy.tpr'
    t = 200000
    outgro = 'temp.gro'
    masses = 46*np.ones([29,1])
    linds = [0,1,2,6,7,21,22,27,28]
    masses[linds] = 72
    cutoff = 0.5
    resInfo = ('PAQ',0)
    ats=29
    endInds = [0,27]
    bins=10
    (rqs,rqsnorm,rnqs,rnqsnorm,hqs,hnqs) = clustStruct(t,xtc,tpr,outgro,cutoff,resInfo,ats,masses,endInds,bins)
    
    '''
    ===============
    =Visualization=
    ===============
    for i in range(50):
        if np.sum(hqs[i,:])!=0:
            fig = plt.figure()
            plt.bar(range(10),hqs[i,:]/np.sum(hqs[i,:]),color='r')
            plt.bar(range(10),hnqs[i,:]/np.sum(hnqs[i,:]),color='b')
            plt.title("size {}".format(i))
    plt.figure()
    plt.scatter(range(len(rqs)),rqs,color='r')
    plt.scatter(range(len(rnqs)),rnqs,color='b')
            
    plt.show()
    '''
        

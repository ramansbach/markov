import numpy as np, imp,os
try:
	mk = imp.load_source('mk','/home/mansbac2/coarsegraining/code/markov/markov.py')
except IOError:
	mk = imp.load_source('mk','/home/rachael/cluster/home/mansbac2/coarsegraining/code/markov/markov.py')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('text',usetex=False)
import subprocess,timeit
#check the number of lines in a file
def checklines(fname):
	f = open(fname,"r")
	nl = 0
	while True:
		line = f.readline()
		if line=="":
			break
		nl+=1
	return nl

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

def getPosQ(t,trj,tpr,outGro,resInfo,ats,writeGro=True):
        if writeGro:
            os.system('echo 0 | trjconv -f ' + trj + ' -o ' + outGro + ' -b ' + str(t) + ' -e ' + str(t) + ' -s ' + tpr)
        #subprocess.check_output([])        
        return readGroQ(outGro,resInfo,ats)

#function that takes a 1d list of peptide locations, the number of atoms in the
#peptide, and list of 1s and 0s where 1 is charged and 0 is uncharged, corresponding
#to molecules that are charged and uncharged, and endInds gives the locations
#of the termini in the molecule
def endToEnd(pepList,ats,endInds,qbool):
    molno = len(pepList)/(3*ats)
    qnum = np.sum(qbool)
    nqnum = len(qbool)-qnum
    deeQ = np.zeros(qnum)
    deeNQ = np.zeros(nqnum)
    indq = 0
    indnq = 0
    for i in range(molno):
    
        ix1 = 3*ats*i+3*endInds[0]
        iy1 = ix1+1
        iz1 = iy1+1
        x1 = pepList[ix1]
        y1 = pepList[iy1]
        z1 = pepList[iz1]
        
        ix2 = 3*ats*i+3*endInds[1]
        iy2 = ix2+1
        iz2 = iy2+1
        x2 = pepList[ix2]
        y2 = pepList[iy2]
        z2 = pepList[iz2]
        
        dee = np.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
        if dee > 6:
            print dee
        if qbool[i]:
            deeQ[indq] = dee
            indq+=1
        else:
            deeNQ[indnq] = dee
            indnq+=1
            
            
    return (deeQ,deeNQ)

#gets all the end to end distances of the monomers in a cluster (simply the distance between termini beads)
#write out t clustSize end to end dist for each monomer, then each monomer is tagged with cluster size
#and time and also write out different file for Q and NQ molecules
def endToEndMono(t,xtc,tpr,outgro,cutoff,resInfo,ats,masses,endInds,outfq,outfnq):
    (peps,box_length,qbool) = getPosQ(t,xtc,tpr,outgro,resInfo,ats,writeGro=True)
    os.system('rm '+outgro)
    pots = range(len(peps)/3/ats)
    inds = np.zeros(len(peps)/3/ats)
    ind = 1
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
        if ind == 8:
            print curr
            
        pepList = mk.fixPBC(pepList,box_length,ats,cutoff)
        (deeQ,deeNQ) = endToEnd(pepList,ats,endInds,qbool[clust])
        for i in range(len(deeQ)):
            if deeQ[i] > 6:
                print deeQ
            
            outfq.write('{0} {1} {2}\n'.format(t,len(clust),deeQ[i]))
        for i in range(len(deeNQ)):
            if deeNQ[i] > 6:
                print deeNQ
            outfnq.write('{0} {1} {2}\n'.format(t,len(clust),deeNQ[i]))
        ind+=1
    return
        
#returns the following data
#data: cluster size + charged/uncharged molecules + COM + distance of ends from COM (distribution of cluster end distance from COM?)
#resInfo is ('PAQ',0) for DFAG,masses is a list of the beads in .gro file order mass in amu (either 72 or 46)
#endInds are the indices in the .gro file containing the 4 charged beads (in DFAG, that is [0,1,27,28]
#outf is a file open for writing passed into the function
def clustStruct(t,xtc,tpr,outgro,cutoff,resInfo,ats,masses,endInds,bins,outf,outfrs,raboutz=False,writeGro=True,elongRat=8,elongMols=400,rm=True):
    (peps,box_length,qbool) = getPosQ(t,xtc,tpr,outgro,resInfo,ats,writeGro=writeGro)

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
        indmax = np.argmax(gStuff[0])
        gVals = gStuff[0]
        gVals = np.sort(gVals)
        gMoms = gStuff[1]
        gPrinc = gMoms[:,indmax]
        g2 = gMoms[:,np.argsort(gVals)[1]] #second principal eigenvector
        rat = gVals[2]/gVals[0]
    
        if raboutz:
            #find the average perpendicular distance from the z' axis
            #where the z' axis is the principal moment of the gyration tensor
                        
            (rq,hq,rsq1) = radialHistZ(dEq,bins,gPrinc)
            (rnq,hnq,rsnq1) = radialHistZ(dEnq,bins,gPrinc)
            (rq2,hq2,rsq2) = radialHistZ(dEq,bins,g2)
            (rnq2,hnq2,rsnq2) = radialHistZ(dEnq,bins,g2)

        else:
            (rq,hq,rsq) = radialHist(dEq,bins)
            (rnq,hnq,rsnq) = radialHist(dEnq,bins)
	nC = len(clust)
	nQ = np.sum(qbool[clust])
	outf.write("\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(t,nC,nQ,float(nQ)/(float(nC)),gVals[0],gVals[1],gVals[2],rat,rq,rnq))
	for w in range(bins):
		outf.write("\t{}".format(hq[w]))
	for w in range(bins):
		outf.write("\t{}".format(hnq[w]))
        if raboutz:
                for w in range(len(rsq1)):
                        outfrs.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(t,nC,rsq1[w],rsq2[w],1))
                for w in range(len(rsnq1)):
                        outfrs.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(t,nC,rsnq1[w],rsnq2[w],0))
        else:
                for w in range(len(rsq)):
		        outfrs.write("{0}\t{1}\t{2}\t{3}\n".format(t,nC,rsq[w],1))#1 for charged
 	        for w in range(len(rsnq)):
		        outfrs.write("{0}\t{1}\t{2}\t{3}\n".format(t,nC,rsnq[w],0))#0 for charged        
	outf.write("\n")
     
        hQs[len(clust),:] += hq
        rQs[len(clust)] += rq #track mean distance from COM
        rQsnorm[len(clust)]+=1 #normalization for mean (let's not do standard deviation yet)
        rnQs[len(clust)]+= rnq
        rnQsnorm[len(clust)]+=1
        hnQs[len(clust),:] += hnq
        '''
        Error Checking Visualization
        
        fig3 = plt.figure()
        dEqr = np.reshape(dEq,[len(dEq)/3,3])
        dEnqr = np.reshape(dEnq,[len(dEnq)/3,3])
        ax3 = fig3.add_subplot(111,projection='3d')
        ax3.scatter(dEnqr[:,0],dEnqr[:,1],dEnqr[:,2],c='b',marker='o')
        ax3.scatter(dEqr[:,0],dEqr[:,1],dEqr[:,2],c='r',marker='o')
        ax3.scatter(0,0,0,c='g',marker='^')
        ax3.quiver(0,0,0,gPrinc[0],gPrinc[1],gPrinc[2])
        #print "rq is {}".format(rq)
        #print "rnq is {}".format(rnq)
        fig = plt.figure()
        plt.subplot(211)
        plt.bar(range(bins),hq/float(np.sum(hq)),color='r')
        plt.subplot(212)
        plt.bar(range(bins),hnq/float(np.sum(hnq)),color='b')
        #ax.set_xlabel('X Label')
        #ax.set_ylabel('Y Label')
        #ax.set_zlabel('Z Label')
        plt.show()
        if ind > 1:
            pots = []
        '''
        ind+=1
        '''
        for i in range(len(rQs)):
        rQs[i] /= rQsnorm[i]
        rnQs[i] /= rnQsnorm[i]
        '''
    if rm:
        os.system('rm '+outgro)
    return (rQs,rQsnorm,rnQs,rnQsnorm,hQs,hnQs)

def elongHist(gVals,gMoms,pepList,bins):
    print("under construction")
    
def radialHistZ(beadDists,bins,gPrinc):
    #the same as radialHist, except it doesn't compute the distance from the com
    #instead, it computes the perpendicular distance from the principal
    #moment of the gyration tensor
    #com is at (0,0,0) by construction
    h = np.zeros([1,bins])
    rzs = np.zeros([len(beadDists)/3])
    zx = gPrinc[0]
    zy = gPrinc[1]
    zz = gPrinc[2]
    rzmax = 0
    for i in range(len(beadDists)/3):
        vx = beadDists[3*i]
        vy = beadDists[3*i+1]
        vz = beadDists[3*i+2]
        dot = vx*zx+vy*zy+vz*zz
        r = np.sqrt((vx-dot*zx)*(vx-dot*zx)+(vy-dot*zy)*(vy-dot*zy)+(vz-dot*zz)*(vz-dot*zz))
        rzs[i] = r
        if r > rzmax:
            rzmax = r
    (h,bE) = np.histogram(rzs,bins=bins,range=(0.0,rzmax))
    
    return (np.mean(rzs),h,rzs)

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
    return (np.mean(rs),h,rs)
        
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
    '''
    ===============
    ===DEBUGGING===
    ===============
    '''
    testarray = np.array([0.0,0.0,0.0,1.0,0.0,0.0,2.0,0.0,0.0,3.0,0.0,0.0,4.0,0.0,0.0,1.0,-2.0,0.0,1.5,-1.5,0.0,1.75,-1.4,0.0,1.9,-1.5,0.0,2.0,-2.0,0.0,1.0,-4.0,0.0,1.5,-4.5,0.0,2.0,-5.0,0.0,2.5,-4.5,0.0,3.0,-4.0,0.0])
    testqbool = np.array([1,1,1])
    endInds = [0,4]
    ats=5
    (tQ,tNQ) = endToEnd(testarray,ats,endInds,testqbool)
    print tQ
    print tNQ
    '''
    ==================
    ===Sanity Check===
    ==================
    testDists = np.array([0.0,7.0/(4*np.sqrt(2)),9.0/(4*np.sqrt(2)),0.0,9.0/(4*np.sqrt(2)),7.0/(4*np.sqrt(2))])
    testV = np.array([0.0,1.0/np.sqrt(2),1.0/np.sqrt(2)])
    print radialHistZ(testDists,10,testV)    
  
    xtc = 'md_noW.xtc'
    tpr = 'md_dummy.tpr'
    t = 200000
    outgro = 'temp.gro'
    outf = open('test.txt','w')
    masses = 46*np.ones([29,1])
    linds = [0,1,2,6,7,21,22,27,28]
    masses[linds] = 72
    cutoff = 0.5
    resInfo = ('PAQ',0)
    ats=29
    endInds = [0,27]
    bins=10
    (rqs,rqsnorm,rnqs,rnqsnorm,hqs,hnqs) = clustStruct(t,xtc,tpr,outgro,cutoff,resInfo,ats,masses,endInds,bins,outf,raboutz=True)
    
    
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
        

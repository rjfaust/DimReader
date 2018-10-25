import numpy as np
import scipy
import pickle
import math
import sys
import json
import copy
import blockPowerIter
import multiprocessing

def load_data(filename):
    f = open(filename, "rb")
    m = pickle.load(f)
    d = pickle.load(f)
    return (m, d)

def calcSimMat(pPoints,sigma):
    n = len(pPoints)
    m = len(pPoints[0])

    simMat = [0] * n

    for i in range(n):
        simMat[i] = [0] * n

    for i in range(n):
        j = i
        while j < n:
            mag = np.inner(np.subtract(pPoints[i], pPoints[j]), np.subtract(pPoints[i], pPoints[j]))
            sij = pow(math.e, -mag / pow(sigma, 2))
            simMat[i][j] = sij
            simMat[j][i] = sij

            j = j + 1
    # import pylab
    # pylab.matshow(simMat,cmap="Reds")
    # pylab.colorbar()
    # pylab.show()
    return simMat


def calcExpandedLs(simMat,d):
    n=len(simMat)

    ls = np.zeros((n * d, n * d))

    for i in range(n):
        for j in range(n):
            if i == j:
                val = sum(simMat[i]) - simMat[i][j]
            else:
                val = -simMat[i][j]
            m = np.eye(d) * val
            ls[i * d:(i + 1) * d, j * d:(j + 1) * d] = m
    return ls


def calcLs(simMat,d):
    ls = []
    # ls = []
    n=len(simMat)

    for i in range(n):
        ls.append([])
        for j in range(n):
                if i==j:
                    ls[i].append(sum(simMat[i])-simMat[i][j])
                else:
                    ls[i].append(-simMat[i][j])
    return ls

def findPert(mat,ls,penalty):
    n = len(mat[0])
    m1=np.dot(np.transpose(mat),mat)
    m2=np.multiply(penalty,ls)
    result = np.dot(np.transpose(mat),mat)-np.multiply(penalty,ls)
    vals,vecs = scipy.linalg.eigh(m1)
    # print("BTB: ",vals)
    vals,vecs = scipy.linalg.eigh(ls)
    # print("LS: ",vals)
    posCorrection = max(vals)*penalty*np.eye(len(ls))


    vals, vecs = scipy.linalg.eigh(result)
    # print(vals)

    result = np.add(posCorrection,result)
    # vals, vecs= scipy.linalg.eigh(np.dot(np.transpose(blocks),blocks))
    # print(vals)
    # vals,vecs = scipy.linalg.eigh(-np.multiply(penalty,ls))
    # print(vals)
    vals, vecs = scipy.linalg.eigh(result)
    # print(vals)
    # vals, vecs = scipy.linalg.eig(result)
    # v1=np.dot(m1,vecs[:,-1])
    # v2=np.dot(m2,vecs[:,-1])

    return vals,vecs

def findPertFromBlocks(blocks,ls,penalty,which,d):
    simPen = np.multiply(penalty,ls)
    bT =np.array(blocks).transpose((0,2,1))
    bTb = []
    n = len(blocks)
    for i in range(n):
        bTb.append(np.dot(bT[i],blocks[i]))

    def blockDot(v):
        resultV = np.zeros(n*d)
        for i in range(n):
            vSect = v[i * d:(i + 1) * d]
            sect = np.dot(bTb[i], vSect)
            resultV[i * d:(i + 1) * d] = sect
        return resultV

    v = np.random.randint(0,100,n*d)

    def lsDot(v):
        print("Mult ls")
        #
        # resultV = np.zeros(n*d)
        # identRow = np.eye(d)
        # for i in range(n-1):
        #     identRow = np.hstack((identRow,np.eye(d)))
        # for i in range(n):
        #     row = np.repeat(ls[i],d)
        #     v2 = np.multiply(row,v)
        #     resultV[i * d:(i + 1) * d] = np.dot(identRow,v2)
        #
        # return resultV
        procs = []
        cpus = multiprocessing.cpu_count()
        resultV = multiprocessing.Array('d', range(n * d))
        # for i in range(n * d):
        #     resultV[i] = 0
        chunksize = int(np.floor(n / cpus))
        for i in range(cpus):
            minI = chunksize * i
            if i < cpus - 1:
                maxI = chunksize * (i + 1)
            else:
                maxI = n

            procs.append(multiprocessing.Process(target=lsDotFunc,
                                                 args=(n, d, ls,v, minI, maxI, resultV)))

        for proc in procs:
            proc.start()

        for proc in procs:
            proc.join()

        result = []
        for i in range(n * d):
            result.append(resultV[i])
        return result

    print("init")
    powerIter = blockPowerIter.BlockPowerIteration(n * d, lsDot,1)
    print("finding largest eig of ls")
    eigVects, eigVals = powerIter.calculateEig()
    correction = penalty*eigVals[0]
    # correction = 0

    def matDot(v):
        print("mult mat")
        # resultV = np.zeros(n*d)
        # matDotFunc(n,d,bTb, simPen,correction, v, 0, n,resultV)
        procs = []
        cpus = multiprocessing.cpu_count()
        resultV = multiprocessing.Array('d',range(n*d))
        for i in range(n*d):
            resultV[i]=0
        chunksize = int(np.floor(n/cpus))
        for i in range(cpus):
            minI = chunksize *i
            if i < cpus -1:
                maxI = chunksize*(i+1)
            else:
                maxI = n

            procs.append(multiprocessing.Process(target=matDotFunc, args=(n,d,bTb,simPen,correction,v,minI,maxI,resultV)))

        for proc in procs:
            proc.start()

        for proc in procs:
            proc.join()


        result = []
        for i in range(n*d):
            result.append(resultV[i])
        return result

    print("finding eig")
    powerIter = blockPowerIter.BlockPowerIteration(n*d,matDot,which+1)
    eigVects,eigVals = powerIter.calculateEig()
    print(eigVals)

    return eigVals,eigVects

def lsDotFunc(n,d,ls,v,minI, maxI,resultV):
    ident = np.eye(d)
    identRow = np.tile(ident,(1, n))
    # identRow = np.zeros((d,n*d))
    # for i in range(n):
    #     identRow[:,i*d:(i+1)*d] = np.eye(d)
    for i in range(minI, maxI):
        # print("MIN: ", minI, ", ", i , " Max:", maxI)
        row = np.repeat(ls[i],d)
        v2 = np.multiply(row,v)
        resultV[i * d:(i + 1) * d] = np.dot(identRow,v2)

def matDotFunc(n,d,bTb, simPen,correction, v, minI, maxI, resultV):
    ident = np.eye(d)
    identRow= np.tile(ident,(1,n))
    # identRow = np.zeros((d,n*d))
    # for i in range(n):
    #     identRow[:,i*d:(i+1)*d] = np.eye(d)
    for i in range(minI, maxI):
        ident= np.eye(d)
        identRow[:,i * d:(i + 1) * d] = np.add(bTb[i], np.add(-simPen[i][i]*ident,correction*ident))
        simRow = -simPen[i]
        simRow[i]=1
        v2 = np.multiply(np.repeat(simRow,d),v)
        temp = np.dot(identRow,v2)
        # temp[:,i * d:(i + 1) * d] = np.add(temp[:,i * d:(i + 1) * d],correction*np.eye(d))
        resultV[i * d:(i + 1) * d] =temp
        identRow[:,i * d:(i + 1) * d] = ident

        # for j in range(n):
        #     if i == j:
        #         chunk = np.subtract(bTb[i], simPen[i][j] * np.eye(d))
        #         chunk = np.add(chunk,correction*np.eye(d))
        #         vSect = v[j * d:(j + 1) * d]
        #         sect = np.dot(chunk, vSect)
        #         resultV[i * d:(i + 1) * d] = np.add(resultV[i * d:(i + 1) * d], sect)
        #     else:
        #         chunk = (-1 * simPen[i][j]) * np.eye(d)
        #         vSect = v[j * d:(j + 1) * d]
        #         sect = np.dot(chunk, vSect)
        #         resultV[i * d:(i + 1) * d] = np.add(resultV[i * d:(i + 1) * d], sect)




def extract_block_submatrix(m, n, d, k):
    return scipy.array(list(m[k*x:k*x+k, d*x:d*x+d] for x in range(n)))


def recover_pert(mat,data,sigma,penalty,which):
    n = len(data.points)
    k = len(data.points[0])
    d = len(data.origPoints[0])
    pPoints = data.points
    # print("calc simMat")
    simMat = calcSimMat(pPoints,sigma)
    # print("calc ls")
    ls = calcExpandedLs(simMat,d)
    # blocks = extract_block_submatrix(mat,n,d,k)
    # print("finding pert")
    vals,vecs=findPert(mat,ls,penalty)
    frontInd = 0
    backInd = len(vals)-1
    biggerInd =0
    count =-1
    # while count < which:
    #     m = max(abs(vals[frontInd]),abs(vals[backInd]))
    #     if m==abs(vals[frontInd]):
    #         biggerInd =frontInd
    #         frontInd += 1
    #     else:
    #         biggerInd  == backInd
    #         backInd-=1
    #     count +=1
    # print(biggerInd)
    return vecs[:,(len(vals)-1)-which]





def recoverPertFromBlocks(blocks,data,sigma,penalty,which):
    n = len(data.points)
    k = len(data.points[0])
    d = len(data.origPoints[0])
    pPoints = data.points
    print("calc simMat")
    simMat = calcSimMat(pPoints, sigma)
    print("calc ls")
    ls = calcLs(simMat, d)
    print("finding pert")
    vals,vecs = findPertFromBlocks(blocks,ls,penalty,which,d)
    return vecs[which]

def createPlot(data,mat,k,pert,scale=1):
    import pylab

    points = scipy.array(data.points)
    colors = [[0.0, 0.0, 0.0] for i in range(len(points))]
    if(data.classes is not None):
        uniqClass = set(data.classes)
        colors = np.array(colors)
        for i, name in enumerate(uniqClass):
            pts = [i for i, x in enumerate(data.classes) if x == name]
            colors[pts]=[float(data.colors[i][0])/255.0,float(data.colors[i][1])/255.0,float(data.colors[i][2])/255.0]
    elif(data.colors is not None):
        colors = np.divide(data.colors,255)
    else:
    #    data.colors = data.calcColors()
        if data.colors is not None:
            colors = np.divide(data.colors,255)

    # scale = 1
    for i in range(150):
        proj =  np.dot(mat[k*i:k*i+k],pert)
        pylab.plot([points[i, 0], points[i, 0] + proj[0] * scale],
                   [points[i, 1], points[i, 1] + proj[1] * scale], 'k-')

    pylab.scatter(points[:, 0], points[:, 1], c=colors)
    pylab.show()




def createPlotFromBlock(data,blocks,k,pert,scale=1):
    import pylab

    points = scipy.array(data.points)
    colors = [[0.0, 0.0, 0.0] for i in range(len(points))]
    if (data.classes is not None):
        uniqClass = set(data.classes)
        colors = np.array(colors)
        for i, name in enumerate(uniqClass):
            pts = [i for i, x in enumerate(data.classes) if x == name]
            colors[pts] = [float(data.colors[i][0]) / 255.0, float(data.colors[i][1]) / 255.0,
                           float(data.colors[i][2]) / 255.0]
    elif (data.colors is not None):
        colors = data.colors
    # scale = 500

    n = len(data.points)

    pylab.scatter(points[:, 0], points[:, 1], c=colors)
    pylab.show()


def createMatrixFromBlocks(blocks,n,d):
    mat =[]
    for i in range(n):
        xrow = blocks[i][0]
        yrow = blocks[i][0]
        if (i > 0):
            zeros = np.zeros(i * d)
            xrow = np.hstack((zeros, xrow))
            yrow = np.hstack((zeros, yrow))
        zeros = np.zeros((n - (i + 1)) * d)
        xrow = np.hstack((xrow, zeros))
        yrow = np.hstack((yrow, zeros))
        mat.append(xrow)
        mat.append(yrow)
    return mat

def writePert(pert,sigma,penalty,inFile,which):
    inFile =inFile.replace(".pkl",'')
    filename = inFile + "_pert_s"+str(sigma) + "_l" + str(penalty)  + "_"+str(which)+".csv"
    f = open(filename, "w")
    for row in pert:
        line = ""
        for i in range(len(row) - 1):
            line = line + str(row[i]) + ","
        line = line + str(row[len(row) - 1]) + "\n"
        f.write(line)
    f.close()





if __name__=="__main__":
    # sys.argv=[None, "Matrices/Iris/IrisTsneMat.pkl","10","50","0"]
    # sys.argv = [None,"IrisTsneBlocks.pkl",10,30,0]
    if(len(sys.argv)<=3):
        print("Invalid input")
        print("Expected Input: filename sigma (similarity parameter) lambda (similarity penalty) which vector (opt)")
        exit(0)
    if(len(sys.argv)==5):
        which= int(sys.argv[4])
    else:
        which=0
    mat,data = load_data(sys.argv[1])
    n = len(data.points)
    k = len(data.points[0])
    d = len(data.origPoints[0])

    s = np.shape(mat[0][0])

    if s and d== len(mat[0][0]):
        pert = recoverPertFromBlocks(mat,data,float(sys.argv[2]),float(sys.argv[3]),which)
        print("creatingMat")
        # mat = createMatrixFromBlocks(mat,n,d)
    elif len(mat[0])!=n*d:
        blocks = np.reshape(mat,(n,k,d))
        pert = recoverPertFromBlocks(blocks,data,float(sys.argv[2]),float(sys.argv[3]),which)
    else:
    # pylab.subplot(2, 1, 1)
        pert = recover_pert(mat,data,float(sys.argv[2]),float(sys.argv[3]),which)
        print("Perturbation:", np.reshape(pert, (n, d)))
        print(list(pert).index(max(pert)))
        print(sum(np.power(pert,2)))

        createPlot(data,mat,k,pert)
    writePert(np.reshape(pert,(n,d)),float(sys.argv[2]),float(sys.argv[3]),sys.argv[1],which)



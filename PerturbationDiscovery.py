import numpy as np
import math
import sys
import json
import blockPowerIteration
import datetime
import time
import multiprocessing
from DimReader import runTangentMapProjection

def load_data(filename):
    f = open(filename, "r")
    tMap = json.loads(f.read())
    return  tMap

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

    return simMat

def calcLs(simMat):
    ls = []
    n=len(simMat)

    for i in range(n):
        ls.append([])
        for j in range(n):
                if i==j:
                    ls[i].append(sum(simMat[i])-simMat[i][j])
                else:
                    ls[i].append(-simMat[i][j])
    return ls

def findPert(tangents,ls,penalty,which,d):
    simPen = np.multiply(penalty,ls)
    bT =np.array(tangents).transpose((0,2,1))
    bTb = []
    n = len(tangents)
    for i in range(n):
        bTb.append(np.dot(bT[i],tangents[i]))


    v = np.random.randint(0,100,n*d)

    def lsDot(v):
        procs = []
        cpus = multiprocessing.cpu_count()
        resultV = multiprocessing.Array('d', range(n * d))

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

    powerIter = blockPowerIteration.BlockPowerIteration(n * d, lsDot,1)
    print("finding largest eig of ls")
    eigVects, eigVals = powerIter.calculateEig()
    correction = penalty*eigVals[0]


    def matDot(v):
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
    powerIter = blockPowerIteration.BlockPowerIteration(n*d,matDot,which+1)
    eigVects,eigVals = powerIter.calculateEig()
    print(eigVals)

    return eigVals,eigVects

def lsDotFunc(n,d,ls,v,minI, maxI,resultV):
    ident = np.eye(d)
    identRow = np.tile(ident,(1, n))

    for i in range(minI, maxI):

        row = np.repeat(ls[i],d)
        v2 = np.multiply(row,v)
        resultV[i * d:(i + 1) * d] = np.dot(identRow,v2)

def matDotFunc(n,d,bTb, simPen,correction, v, minI, maxI, resultV):
    ident = np.eye(d)
    identRow= np.tile(ident,(1,n))

    for i in range(minI, maxI):
        ident= np.eye(d)
        identRow[:,i * d:(i + 1) * d] = np.add(bTb[i], np.add(-simPen[i][i]*ident,correction*ident))
        simRow = -simPen[i]
        simRow[i]=1
        v2 = np.multiply(np.repeat(simRow,d),v)
        temp = np.dot(identRow,v2)

        resultV[i * d:(i + 1) * d] =temp
        identRow[:,i * d:(i + 1) * d] = ident


def recover_pert(tMap,sigma,penalty, which):
    projPoints = []
    tangents = []

    for pt in tMap:
        projPoints.append(pt["range"])
        tangents.append(pt["tangent"])
    d = len(tMap[0]["domain"])
    simMat = calcSimMat(projPoints,sigma)
    ls = calcLs(simMat)
    vals,vecs = findPert(tangents,ls,penalty,which,d)
    return vecs[which]


def writePert(pert,sigma,penalty,inFile,which):
    inFile =inFile[:inFile.rfind(".")]
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
    if(len(sys.argv)<=3):
        print("Invalid input")
        print("Expected Input: [tangent map file] [sigma (similarity parameter)] [lambda (similarity penalty)] [which vector (optional)]")
        exit(0)
    if(len(sys.argv)==5):
        which= int(sys.argv[4])
    else:
        which=0

    tMap = load_data(sys.argv[1])
    pert = recover_pert(tMap,float(sys.argv[2]), float(sys.argv[3]), which)
    n = len(tMap)
    d = len(tMap[0]["domain"])
    writePert(np.reshape(pert,(n,d)),float(sys.argv[2]),float(sys.argv[3]),sys.argv[1],which)
    date = str(datetime.datetime.fromtimestamp(time.time())).replace(" ", "_")
    date = date.split(".")[0]

    inputFile = sys.argv[1]

    prefix = inputFile[:inputFile.rfind(".")] + "_Tangent-Map_" + date

    runTangentMapProjection(tMap,[np.reshape(pert,(n,d)).tolist()],prefix)



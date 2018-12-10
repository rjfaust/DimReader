import sys
import DualNum
import csv
import numpy as np
import tSNE
import multiprocessing
import datetime
import time
import Grid
import json

class ProjectionRunner:
    def __init__(self,projection,params=None):
        self.params = params
        self.projection = projection
        self.firstRun = False

    def calculateValues(self, points, perturbations=None):


        self.points = points
        self.origPoints = points
        self.resultVect = [0] * (len(self.points) * 2)

        n = len(points)
        self.dualNumPts = DualNum.DualNum(np.array(points), np.zeros(np.shape(points)))

        self.perturb = perturbations
        if not self.firstRun:
            #initial run of tsne to get seed parameters
            p = self.projection(points, self.params)
            p.run()
            self.firstRun = True
            self.runParams= p.getExecParams()


        if (n != 0):
            procs = []
            cpus = multiprocessing.cpu_count()
            xOutArray = multiprocessing.Array('d', range(n))
            yOutArray = multiprocessing.Array('d', range(n))
            outDotArray = multiprocessing.Array('d', range(2 * n))

            if (cpus > n):
                cpus = 1
            chunksize = int(np.floor(float(n) / cpus))
            for i in range(cpus):
                minI = chunksize * i

                if (i < cpus - 1):
                    maxI = chunksize * (i + 1)
                else:
                    maxI = n

                procs.append(multiprocessing.Process(target=self.loopFunc,
                                                     args=(self.dualNumPts, minI, maxI, outDotArray,
                                                            self.runParams, xOutArray, yOutArray)))

            for proc in procs:
                proc.start()

            for proc in procs:
                proc.join()


            points = []
            self.resultVect = [0] * (2 * n)
            for i in range(n):
                self.resultVect[i] = outDotArray[i]
                self.resultVect[i + n] = outDotArray[i + n]
                points.append([xOutArray[i], yOutArray[i]])

        self.points = points

    def loopFunc(self, pts, minI, maxI, dotArr, runParams, xPts, yPts):
        for i in range(minI, maxI):

            if self.perturb is None:
                pts.dot[i][self.axis] = 1
            else:
                pts.dot[i] = self.perturb[i]
            m = len(self.origPoints[0])
            p = self.projection(self.dualNumPts,self.params)
            p.setExecParams(runParams)


            results = p.run()

            # print("Min: ", minI, "I: ", i, "max: ", maxI)

            dotArr[2 * i] = results[i][0].dot
            dotArr[2 * i + 1] = results[i][1].dot
            n = len(pts.val)

            if i == 0:
                for j in range(len(self.points)):
                    xPts[j] = results[j][0].val
                    yPts[j] = results[j][1].val
            pts.dot[i] = np.zeros(m)



projections = ["tsne", "Tangent-Map"]
projectionClasses=[tSNE.tSNE,None.TangentMapProjection]
projectionParamOpts = [tSNE.paramOpts,[]]


def readFile(filename):
    read = csv.reader(open(filename, 'rt'))

    points = []
    firstLine = next(read)
    headers = []
    rowDat = []
    head = False
    for i in range(0, len(firstLine)):
        try:
            rowDat.append(float(firstLine[i]))
        except:
            head = True
            break
    if head:
        headers = firstLine
    else:
        points.append(rowDat)

    for row in read:
        rowDat = []
        for i in range(0, len(row)):
            try:
                rowDat.append(float(row[i]))
            except:
                print("invalid data type - must be numeric")
                exit(0)
        points.append(rowDat)
    return points

def readTangentMap(filename):
    f = open(filename,"r")
    tMap = json.loads(f.read())
    f.close()

    return tMap

def generateGrid(points):
    gridSize = 10
    points = np.array(points)
    xmax = max(points[:, 0])
    xmin = min(points[:, 0])
    ymax = max(points[:, 1])
    ymin = min(points[:, 1])

    gridCoord = []
    if (xmin == xmax):
        xmax += 1
        xmin -= 1
    if (ymin == ymax):
        ymin -= 1
        ymax += 1
    if ymax > xmax:
        xmax = ymax
    else:
        ymax = xmax
    if ymin < xmin:
        xmin = ymin
    else:
        ymin = xmin
    yrange = ymax - ymin
    xrange = xmax - xmin

    xstep = float(xrange) / (gridSize - 1)
    ystep = float(yrange) / (gridSize - 1)

    for i in range(gridSize + 1):
        gridCoord.append([])
        for j in range(gridSize + 1):
            gridCoord[i].append([(xmin - xstep / 2.0) + xstep * j, (ymax + ystep / 2.0) - ystep * i])
    return gridCoord

def calcGrid(points,dVects):#date, grid, gridCoord, ind):
    gridCoord = generateGrid(points)
    g = Grid.Grid(points, dVects, gridCoord)

    grid = g.calcGridPoints()

    n = len(gridCoord)
    m = len(gridCoord[0])


    gridRows = []
    for row in range(n - 1):
        for col in range(m - 1, ):
            NW = grid[m * row + col]
            NE = grid[m * row + (col + 1)]
            SW = grid[m * (row + 1) + col]
            SE = grid[m * (row + 1) + (col + 1)]
            Row = (n - 2) - row
            Col = col
            d = {"NW":str(NW),
                 "NE":str(NE),
                 "SW":str(SW) ,
                 "SE":str(SE) ,
                 "Row":str(Row) ,
                 "Col":str(Col)}
            gridRows.append(d)
    return gridRows



def runProjection(projection, points, perturbations, parameters,filePrefix):



    pr = ProjectionRunner(projection,parameters)

    n = len(points)


    pts = []
    print("Running DimReader")
    for i in range(len(perturbations)):
        pert = perturbations[i]
        pr.calculateValues(np.array(points), np.array(pert))

        derivVect = pr.resultVect
        projPts = pr.points
        data = []
        for j in range(n):
            data.append({"domain":points[j],
                         "range": projPts[j],
                         "inputPert": pert[j],
                         "outputPert":derivVect[j]
                        })
        if len(perturbations)>1:
            fileName = filePrefix +"_"+str(i)+".dimreader"
        else:
            fileName = filePrefix + ".dimreader"

        output  ={"points":data}
        grid = calcGrid(points,derivVect)
        output.update({"grid":grid})
        f = open(fileName, "w")
        f.write(json.dumps(output))
        f.close()

def tmProject(tMap,pert):
    n = len(pert)
    result = []
    for i in range(n):
        proj = np.dot(tMap[i]["tangent"],pert[i])
        result.append(proj[0])
        result.append(proj[1])

    return result

def runTangentMapProjection(tMap,perts,prefix):

    print("Running DimReader")
    n= len(tMap)
    for i in range(len(perts)):
        pert = perts[i]

        derivVect = tmProject(tMap,pert)

        data = []
        points = []
        for j in range(n):
            data.append({"domain":tMap[j]["domain"],
                         "range": tMap[j]["range"],
                         "inputPert": pert[j],
                         "outputPert":derivVect[j]
                        })
            points.append(tMap[j]["domain"])
        if len(perts)>1:
            fileName = prefix +"_"+str(i)+".dimreader"
        else:
            fileName = prefix + ".dimreader"

        output  ={"points":data}


        grid = calcGrid(points,derivVect)
        output.update({"grid":grid})
        f = open(fileName, "w")
        f.write(json.dumps(output))
        f.close()

if __name__ == "__main__":
    if (len(sys.argv) >= 4):
        inputFile = sys.argv[1]
        perturbFile = sys.argv[2]
        projection = sys.argv[3]

        if str.lower(projection) not in map(str.lower, projections):
            print("Invalid Projection")
            print("Projection Options:")
            for i,opt in enumerate(projections):
                print("\t" + opt)
                print("\t\tOptional Parameters: "+str(projectionParamOpts[i]))
            exit(0)

        projInd = list(map(str.lower, projections)).index(str.lower(projection))


        date = str(datetime.datetime.fromtimestamp(time.time())).replace(" ", "_")
        date = date.split(".")[0]

        prefix = inputFile[:inputFile.rfind(".")] + "_" + projections[projInd] + "_" + date

        if str.lower(projection) != "tangent-map":
            inputPts = readFile(inputFile)

            if str.lower(perturbFile) == "all":
                perturbVects = []
                n, m = np.shape(inputPts)
                for i in range(m):
                    currPert = np.zeros((n, m))
                    for j in range(n):
                        currPert[j][i] = 1
                    perturbVects.append(currPert.tolist())
            else:
                perturbVects = [readFile(perturbFile)]


            if (len(sys.argv)>4):
                params = []
                for i in range(4, len(sys.argv)):
                    params = [sys.argv[i]]
            else:
                params = []


            runProjection(projectionClasses[projInd], inputPts, perturbVects,params,prefix)
            print("Output File Prefix: ",prefix)

        else:
            tMap  = readTangentMap(inputFile)
            if str.lower(perturbFile) == "all":
                perturbVects = []
                n = len(tMap)
                m  = len(tMap[0]["domain"])
                for i in range(m):
                    currPert = np.zeros((n, m))
                    for j in range(n):
                        currPert[j][i] = 1
                    perturbVects.append(currPert.tolist())
            else:
                perturbVects = [readFile(perturbFile)]

            runTangentMapProjection(tMap, perturbVects, prefix)
            print("Output File Prefix: ", prefix)


    else:
        print("DimReader.py [input file] [perturbation file] [Projection] [optional parameters (in order)]")
        print("-To perturb all dimensions, perturbation file = all")
        print("-When using the Tangent Map as the projection, the input file = tangent map file.")
        print("Projection Options:")
        for i,opt in enumerate(projections):
            print("\t" + opt)
            print("\t\tOptional Parameters: " + str(projectionParamOpts[i]))

        exit(0)



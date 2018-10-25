import sys
import DualNum
import csv
import numpy as np
import tSNE
import multiprocessing
import datetime
import time
import Grid

class Projection:
    def __init__(self,projection,params=None):
        self.params = params
        self.projection = projection
        self.firstRun = False
        self.resultVect = [0] * (len(self.points) * 2)

    def calculateValues(self, points=None, perturbations=None):

        if points is None:
            points = self.origPoints
        else:
            self.points = points
            self.origPoints = points
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
            p = projection(self.dualNumPts,self.params)
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
projectionClasses=[tSNE.tSNE]


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

    #    gridPoints = []
    #    for i in range(n):
    #        gridPoints.append([])
    #        for j in range(m):
    #            gridPoints[i].append(vVector[m * i + j])

    # if ind is not None:
    #     gridF = open(date + "_grid" + str(ind) + ".csv", "w")
    # else:
    #     gridF = open(date + "_grid.csv", "w")
    # gridF.write("NW,NE,SW,SE,Row,Col\n")
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


    # gridF.close()

def runProjection(projection, points, perturbations):
    # derivVects = []
    # for pert in perturbations:
    #     projection.calculateValues(np.array(points), np.array(pert))
    #     derivVects.append(projection.resultVect)
    #     projPts = projection.points

    date = str(datetime.datetime.fromtimestamp(time.time())).replace(" ", "_")
    date = date.split(".")[0]
    fileName = date + "_output.csv"

    f = open(fileName, "w")

    n = len(points)

    pts = []
    for i in range(perturbations):
        pert = perturbations[i]
        projection.calculateValues(np.array(points), np.array(pert))

        # derivVects.append(projection.resultVect)
        derivVect = projection.resultVect
        projPts = projection.points
        data = []
        for j in range(n):
            data.append({"Domain":points[j].tolist(),
                         "Range": projPts[i].tolist(),
                         "OutputPert":derivVect[i].tolistt()
                        }
                        )

    headers = "ProjectedX,ProjectedY"
    if (len(perturbations) > 1):
        for i in range(len(perturbations)):
            headers += ",dx" + str(i) + ",dy" + str(i)
        headers += "\n"
    else:
        headers += ",dx,dy\n"
    f.write(headers)
    for i in range(n):
        row = str(projPts[i][0]) + "," + str(projPts[i][1])
        for j in range(len(perturbations)):
            row += "," + str(derivVects[j][2 * i]) + "," + str(derivVects[j][2 * i + 1])
        row += "\n"
        f.write(row)
    f.close()



if __name__ == "__main__":
    if (len(sys.argv) >= 4):
        inputFile = sys.argv[1]
        perturbFile = sys.argv[2]
        projection = sys.argv[3]

        if str.lower(projection) not in map(str.lower, projections):
            print("Invalid Projection")
            print("Projection Options:")
            for opt in projections:
                print("\t" + opt)
            exit(0)

        projInd = list(map(str.lower, projections)).index(str.lower(projection))
        if (projInd < 5):
            inputPts = readFile(inputFile)
        else:
            projection = projectionClasses[projInd]()
            projection.loadMat(inputFile)
            inputPts = projection.origPoints

        if str.lower(perturbFile) == "all":
            perturbVects = []
            n, m = np.shape(inputPts)
            for i in range(m):
                currPert = np.zeros((n, m))
                for j in range(n):
                    currPert[j][i] = 1
                perturbVects.append(currPert)
        else:
            perturbVects = [readFile(perturbFile)]

        points = DualNum.DualNum(inputPts, perturbVects)

        if (projInd < 5):
            if (len(sys.argv) == 5):
                projection = projectionClasses[projInd](float(sys.argv[4]))
            else:
                projection = projectionClasses[projInd]()

        runProjection(projection, inputPts, perturbVects)

    else:
        print("DimReaderScript [input file] [perturbation file] [Projection] [optional parameter]")
        print("For all dimension perturbations, perturbation file = all")
        print("Projection Options:")
        for opt in projections:
            print("\t" + opt)

        exit(0)



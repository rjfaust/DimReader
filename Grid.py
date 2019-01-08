from numpy import linalg
import numpy as np

class Grid():
    def __init__(self,points,resultVect,gridCoords):
        self.points = points
        #currently assumes the points given in a 2D array with the [n][0] being the smallest x and y and
        # [n][m] the largest (x values increase from left to right, y values from bottom to top)
        self.gridCoords = gridCoords
        self.resultVect= resultVect

    def addAvgNeighbors(self,c):
        print("avg neighb init")
        for i in range(self.nrow* self.ncol):
            c.append([])
            self.resultVect.append(0)
            for j in range(self.nrow * self.ncol):
                c[2*len(self.points)+i].append(0)

        # calculate the weights for the average of the neighbors for each vi
        for i in range(len(c[0])):
            cInd = 2 * len(self.points) + i
            c[cInd][i] = 1
            if i % self.ncol == 0:
                # it is on the left edge (no neighbors to the left)
                if i < self.ncol:
                    # it is the top left corner
                    c[cInd][i + 1] = -.5
                    c[cInd][i +self.ncol] = -.5

                elif i >= self.ncol * (self.nrow - 1):
                    # it is the bottom left corner
                    c[cInd][i + 1] = -.5
                    c[cInd][i - self.ncol] = -.5

                else:
                    # interior point on the left edge
                    c[cInd][i - self.ncol] = -1.0 / 3.0
                    c[cInd][i + 1] = -1.0 / 3.0
                    c[cInd][i + self.ncol] = -1.0 / 3.0

            elif (i + 1) % self.ncol == 0:
                # point is on the right edge (no neighbors to the right)
                if i < self.ncol:
                    # it is the top right corner
                    c[cInd][i - 1] = -.5
                    c[cInd][i + self.ncol] = -.5
                elif i >= self.ncol * (self.nrow - 1):
                    # it is the bottom right corner
                    c[cInd][i - 1] = -.5
                    c[cInd][i - self.ncol] = -.5

                else:
                    # interior point on the right edge
                    c[cInd][i - self.ncol] = -1.0 / 3.0
                    c[cInd][i - 1] = -1.0 / 3.0
                    c[cInd][i + self.ncol] = -1.0 / 3.0

            else:
                if i < self.ncol:
                    # point is an interior point on the top row
                    c[cInd][i + self.ncol] = -1.0 / 3.0
                    c[cInd][i - 1] = -1.0 / 3.0
                    c[cInd][i + 1] = -1.0 / 3.0

                elif i >= self.ncol * (self.nrow - 1):
                    c[cInd][i - self.ncol] = -1.0 / 3.0
                    c[cInd][i - 1] = -1.0 / 3.0
                    c[cInd][i + 1] = -1.0 / 3.0
                else:
                    c[cInd][i + self.ncol] = -.25
                    c[cInd][i - self.ncol] = -.25
                    c[cInd][i - 1] = -.25
                    c[cInd][i + 1] = -.25

    def addAvgNeighborsInner(self,c):
        if(self.nrow>2 and self.ncol>2):

            # calculate the weights for the average of the neighbors for each vi
            for i in range(((self.nrow) * (self.ncol))):

                if(i % self.ncol != 0 and (i + 1) % self.ncol != 0 and i<self.ncol * (self.nrow - 1) and i>self.ncol):
                    c.append([])
                    self.resultVect.append(0)
                    cInd = len(c)-1
                    for j in range(self.nrow * self.ncol):
                        c[cInd].append(0)

                    c[cInd][i] = 1
                    #print(str(int(i / self.nrow)) + ", " + str(i % self.nrow))
                    c[cInd][i + self.ncol] = -.25
                    c[cInd][i - self.ncol] = -.25
                    c[cInd][i - 1] = -.25
                    c[cInd][i + 1] = -.25


    def calcGridPoints(self):

        #assumes grid is rectangular (equal number of squares in each row)
        self.nrow = len(self.gridCoords)
        self.ncol = len(self.gridCoords[0])
        c = []
        #initialize c (first n rows are for x derivatives, the second n are for y derivative, the rest
        #are for the average neighbors of each corner)
        print("inits")
        for i in range(2*len(self.points)):
            c.append([])
            for j in range(self.nrow*self.ncol):
                c[i].append(0)


        corners = []

        #iterate over each square
        for i in range(self.nrow-1):
            corners.append([])
            for j in range(self.ncol-1):
                corners[i].append([self.gridCoords[i][j+1],self.gridCoords[i][j],self.gridCoords[i+1][j],self.gridCoords[i+1][j+1]])

        print("corners")
        for k in range(len(self.points)):
            # print(k)
            p = self.points[k]
            i = 0
            j =0
            found = False
            #max iterations is max(num rows, num cols)
            while(not found):
                #find the maxX,minX, miny and maxY of the corners
                maxX = corners[i][j][0][0]
                minX = corners[i][j][0][0]
                maxY= corners[i][j][0][1]
                minY= corners[i][j][0][1]
                for l in range(1,3):
                    if(corners[i][j][l][0]>maxX):
                        maxX= corners[i][j][l][0]
                    elif (corners[i][j][l][0] < minX):
                        minX = corners[i][j][l][0]
                    if (corners[i][j][l][1] > maxY):
                        maxY = corners[i][j][l][1]
                    if (corners[i][j][l][1] < minY):
                        minY = corners[i][j][l][1]

                #check if it is in the current square
                if (p[0]<=maxX or abs(p[0]-maxX)<=pow(10,-5) )and (p[0]>=minX or abs(p[0]-minX)<=pow(10,-5)) and (p[1]>=minY or abs(p[1]-minY)<=pow(10,-5))and (p[1]<=maxY or abs(p[1]-maxY)<=pow(10,-5)):

                    xscaled = (p[0]-corners[i][j][2][0])/(corners[i][j][3][0]-corners[i][j][2][0])
                    yscaled = (p[1]-corners[i][j][2][1])/(corners[i][j][1][1]-corners[i][j][2][1])


                    if (xscaled >= yscaled):
                        # point is in lower triangle
                        weightYv3=0
                        weightYv2 = -1.0
                        weightYv1 = 1.0

                        weightXv3 = -1.0
                        weightXv2 = 1.0
                        weightXv1 = 0

                        c[2*k][(i+1)*self.ncol+j+1] = weightXv2
                        c[2*k+1][(i+1)*self.ncol+j+1] = weightYv2

                    else:
                        weightYv1 = 0
                        weightYv3 = -1.0
                        weightYv2 = 1.0

                        weightXv3 = 0
                        weightXv2 = -1.0
                        weightXv1 = 1.0

                        c[2*k][i * self.ncol + j] = weightXv2
                        c[2*k+1][i*self.ncol+j] = weightYv2
                    c[2*k][i * self.ncol + (j + 1)] = weightXv1
                    c[2*k][(i + 1) * self.ncol + j] = weightXv3
                    c[2*k+1][i * self.ncol + (j + 1)] = weightYv1
                    c[2*k+1][(i + 1) * self.ncol + j] = weightYv3
                    found = True
                    #if it is of smaller y value, add to row counter
                elif p[1]<minY:# or p[1]>maxY:
                    # if it is of greater x value, add to col counter
                    if(p[0]>maxX):# or p[0]<minX):
                        j+=1
                    i += 1
                # if it is of greater x value, add to col counter
                elif (p[0]>maxX):# or p[0]<minX):
                    j += 1

        self.addAvgNeighbors(c)

        #Perform matrix calculations to get V0,...,Vn
        c = np.array(c)
        resultVect = np.array(self.resultVect)

        self.vVect = linalg.lstsq(c,resultVect,rcond=0.000001)

        return self.vVect[0]





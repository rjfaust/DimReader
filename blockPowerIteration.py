import numpy as np
import math


class BlockPowerIteration():
    def __init__(self,n,multFunc,numVects,initB=None):
        if(initB==None):

            self.initB = [1.0 / math.sqrt(n)] * n

        self.numVects = numVects
        self.multFunc = multFunc
        self.n = n

    def calculateEig(self):
        eigVects = [0] * self.numVects
        eigVals = [0] * self.numVects


        for i in range(self.numVects):
            b = self.initB

            oldB = [0] * self.n


            while sum(abs(np.power(b, 2) - np.power(oldB, 2))) > pow(10, -5):

                oldB = b

                Ab = self.multFunc(b)
                norm = sum(np.power(Ab, 2))

                norm = np.power(norm, .5)
                b = Ab / norm

                newV = b
                for v in range(0, i):
                    newV = newV - eigVects[v] * sum(b * eigVects[v])
                b = newV


            eigVects[i] = b
            ab = self.multFunc(b)

            val = 0

            for j in range(len(b)):
                if (abs(b[j]) > pow(10, -5)):
                    val = ab[j] / b[j]

            if (val == 0):

                eigVects[i] = np.array([1 / math.sqrt(self.n)] * self.n)
                self.initB = [0] * len(self.A)
                for j in range(self.n):
                    n = sum(range(1, self.n + 1))
                    if (j % 2 == 0):
                        self.initB[j] = math.sqrt((j + 1) / n)
                    else:
                        self.initB[self.n - j - 1] = math.sqrt((j + 1) / n)

            eigVals[i] = val

        return eigVects, eigVals


class BlockInversePowerIteration():
    def __init__(self, n, multFunc, numVects, initB=None, singular=False):
        if (initB == None):
            self.initB = [1.0 / math.sqrt(n)] * n
            self.initB = [self.initB] * numVects
            self.dNum = False

            self.initGiven = False
        else:
            self.initB = initB
            self.initGiven = True


        self.n = n
        self.multFunc = multFunc
        self.numVects = numVects
        self.singular = singular
        self.oldParam = [0] * self.numVects

    def calculateEig(self):
        eigVects = [0] * self.numVects
        eigVals = [0] * self.numVects

        for i in range(self.numVects):

            b = self.initB[i]

            oldB = [0] * self.n

            oldB = np.array(oldB)


            while sum(abs(np.power(b, 2) - np.power(oldB, 2))) > pow(10, -9):

                oldB = b

                self.oldParam[i] = oldB



                c = BlockConjugateGradient(self.n,self.multFunc, b)

                b = c.solve(np.array([0] * self.n))

                newV = b
                for v in range(0, i):
                    newV = newV - np.dot(eigVects[v],np.dot(b, eigVects[v]))
                    # newV= newV-eigVects[v] * sum(b * eigVects[v])

                b = newV
                norm = pow(sum(pow(b, 2)), .5)

                if (abs(norm) > 0):
                    b = b / norm
                else:
                    b = b



            eigVects[i] = b
            ab = self.multFunc(b)
            val = 0

            for p in range(self.n):
                if (abs(b[p]) > pow(10, -5)):
                    val = ab[p] / b[p]
            if (val == 0):
                eigVects[i] = np.array([1 / math.sqrt(n)] * n)
                if (not self.initGiven):
                    for l in range(self.numVects):
                        for j in range(self.n):
                            k = sum(range(1, self.n + 1))
                            if (j % 2 == 0):
                                self.initB[l][j] = math.sqrt((j + 1) / k)
                            else:
                                self.initB[l][self.n - j - 1] = math.sqrt((j + 1) / k)

            eigVals[i] = val

        return eigVects, eigVals


class BlockConjugateGradient():
    def __init__(self, n, multFunc, b, singular=False):
        self.n =n
        self.multFunc = multFunc
        self.b = b
        self.singular = singular


    def solve(self, initX):
        X = initX


        r = np.subtract(self.b, self.multFunc(X))
        rT = np.transpose(r)

        P = r
        for i in range(self.n):

            PT = P.T
            rsold = rT.dot(r)
            AP = self.multFunc(P)
            PTAP = PT.dot(AP)

            # Allowed? What does it mean for PTAP to be zero?
            if (round(PTAP, 6) == 0):
                break
            a = rsold / PTAP



            X = X + a * (P)
            r = r - a * AP
            rT = r.T
            rsnew = rT.dot(r)




            val = math.sqrt(rsnew)

            if (val < pow(10, -10)):
                break

            B = rsnew / rsold


            P = r + B * (P)

        return X



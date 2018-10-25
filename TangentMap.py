import numpy as np
import sys
import DimReader

def recoverMatrix2(matToRec, n, m):
    origMat = []
    for i in range(m * n):
        origMat.append([])
    for i in range(m):
        print("RUN: ", i)
        vect = []
        for j in range(n):
            p = np.zeros(m)
            p[i] = 1
            vect.append(p)
        v = matToRec.calculateAllValues(perturbations=vect)
        print("\n\n\n", np.shape(v), "\n\n\n")
        for j in range(n):
            origMat[j * m + i] = v[j]
    # origMat.append(v)

    return np.transpose(origMat)


def recoverMatrix3(matToRec, n, m):
    origMat = []
    for i in range(m * n):
        origMat.append([])
    for i in range(m):
        print("RUN: ", i)
        vect = []
        p = np.zeros(m)
        p[i] = 1
        vect = np.tile(p, n)
        # for j in range(n):
        #     vect.append(p)
        v = matToRec.dot(vect)
        print("\n\n\n", np.shape(v), "\n\n\n")
        for j in range(n):
            p = np.zeros(2 * n)
            p[2 * j] = v[2 * j]
            p[2 * j + 1] = v[2 * j + 1]
            origMat[j * m + i] = p
    # origMat.append(v)

    return np.transpose(origMat)


def recoverMatrix(matToRec, n, m):
    origMat = []
    for i in range(m):
        print("RUN: ", i)
        vect = np.zeros(m)
        vect[i] = 1
        v = matToRec.dot(vect)
        origMat.append(v)

    return np.transpose(origMat)


def recoverBlocksOnly(matToRec, n, m):
    origMat = []
    for i in range(m):
        print("RUN: ", i)
        vect = []
        p = np.zeros(m)
        p[i] = 1
        vect = np.tile(p, n)
        # for j in range(n):
        #     p = np.zeros(m)
        #     p[i] = 1
        #     vect.append(p)
        v = matToRec.dot(vect)
        print("\n\n\n", np.shape(v), "\n\n\n")
        # for j in range(n):
        #     p = np.zeros(2 * n)
        #     p[2 * j] = v[2 * j]
        #     p[2 * j + 1] = v[2 * j + 1]
        origMat.append(v)

    return np.transpose(origMat)



if __name__ == "__main__":

    if (len(sys.argv) >= 4):
        inputFile = sys.argv[1]
        projection = sys.argv[3]


        if str.lower(projection) not in map(str.lower, DimReader.projections):
            print("Invalid Projection")
            print("Projection Options:")
            for opt in DimReader.projections:
                print("\t" + opt)
            exit(0)

        projInd = list(map(str.lower, DimReader.projections)).index(str.lower(projection))
        if (projInd < 5):
            inputPts = DimReader.readFile(inputFile)
        else:

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







    else:
        print("DimReaderScript [input file] [perturbation file] [Projection] [optional parameter]")
        print("For all dimension perturbations, perturbation file = all")
        print("Projection Options:")
        for opt in projections:
            print("\t" + opt)

        exit(0)
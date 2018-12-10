import numpy as np
import sys
import DimReader
import json

projections = DimReader.projections
projectionParamOpts = DimReader.projectionParamOpts
projectionClasses = DimReader.projectionClasses

def calcTangentMap(points, projection,params):
    n,m = np.shape(points)
    tMap = []
    for i in range(n):
        if type(points[i]) ==type(np.array([])):
            p = points[i].tolist()
        else:
            p = points[i]

        pt = {"domain":p,
              "range":[0,0],
              "tangent":np.zeros((2,m)).tolist()
              }
        tMap.append(pt)
    print("Progress: 0%")
    pj = DimReader.ProjectionRunner(projection, params)
    for i in range(m):
        currPert = np.zeros((n, m))
        for j in range(n):
            currPert[j][i] = 1


        pj.calculateValues(np.array(points),currPert)

        outPerts = pj.resultVect
        projPts = pj.points
        if i==0:
            for j in range(n):
                tMap[j]["range"]= projPts[j]

        for j in range(n):
            tMap[j]["tangent"][0][i] = outPerts[2*j]
            tMap[j]["tangent"][1][i] = outPerts[2*j+1]

        print("Progress: ", int(((i+1)/m)*10000)/100.0 ,"%")
    return tMap




if __name__ == "__main__":

    if (len(sys.argv) >= 3):
        inputFile = sys.argv[1]
        projection = sys.argv[2]


        if str.lower(projection) not in map(str.lower, DimReader.projections) and str.lower(projection) != "tangent-map":
            print("Invalid Projection")
            print("Projection Options:")
            for opt in DimReader.projections:
                if opt != "Tangent-Map":
                    print("\t" + opt)
            exit(0)

        projInd = list(map(str.lower, DimReader.projections)).index(str.lower(projection))

        inputPts = DimReader.readFile(inputFile)

        if (len(sys.argv)>3):
            params = []
            for i in range(3, len(sys.argv)):
                params = [sys.argv[i]]
        else:
            params = []

        tMap = calcTangentMap(inputPts,projectionClasses[projInd],params)

        fName = inputFile[:inputFile.rfind(".")]+"_TangentMap_"+projections[projInd] +".tmap"

        f = open(fName, "w")
        f.write(json.dumps(tMap))
        f.close()






    else:
        print("DimReaderScript [input file] [Projection] [optional parameters]")
        print("For all dimension perturbations, perturbation file = all")
        print("Projection Options:")
        for opt in projections:
            if opt!="Tangent-Map":
                print("\t" + opt)

        exit(0)
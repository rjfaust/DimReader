# DimReader

## Prerequisites
### Python 3
DimReader was created for Python3.  It may work on Python2 but this has not been tested. 

## Overview

This repository provides the source code for DimReader.  The full paper can be found at https://rjfaust.github.io/DimReader.pdf

DimReader offers 3 methods:
    
    1. DimReader: The user provides a perturbation for each point and it calculates the effect of perturbingi
    each point in the specified way.

    2. Tangent Map Recovery: The calculation of the entire tangent map for a given dataset and projection.

    3. Perturbation Discovery: The calculation of the perturbation that changes the projection the most 
    (requires the tangent map).

## Running DimReader
DimReader can be run with the command

```
python DimReader.py dataFile.csv perturbationFile.csv projection [optional parameters]
```

### Data File
 Should be in CSV format and all data should be numeric.

### Perturbation File
 A csv where each line is the perturbation for the corresponding point.  In a row (perturbation), there should be a value for each dimension of the original data.  For example, if the data point is [1,2,3,4] the corresponding row of the perturbation file might look like "1,0,0,0".

If you would like run DimReader once for each dimension (perturbing the data points in that dimension), enter "all" instead of a perturbation file.

### Projections


Currently, only two types of projections are supported:

    1. tSNE 
    Paremters: Perplexity (default = 20), 
        Number of iterations (default 1000), 
        Number of Dimensions to project to (default 2)

    2. Projection from the Tangent Map 


### Optional Parameters
Some projections have user tunable paramters. If that is the case, those parameters must be entered in the order that they are shown in the description (it will read them in in that order).  For example, for tSNE, if we want to specify the number of iterations, we first have to enter the perplexity. 




## Running Tangent Map Recovery
The command to recover the tangent map is

```
python TangentMap.py dataFile projection [optional parameters]
```

All of these command parameters are the same as those used in the DimReader command (above). 


## Running Perturbation Discovery
Perturbation Discovery finds the perturbation that changes the projection the most under the constraints:

    1. The perturbation for each point should be the one that maximizes the change in the projected point

    2. Points that are projected to similar places should be perturbed in similar ways



The command to run perturbation discovery is 

```
python PertrubationDiscovery.py tangentMapFile.tmap sigma lambda [optional which_vector]
```

### Tangent Map File
The tangent map file should be a file created by runing the Tangent Map Recovery (discussed above)

### Sigma
Sigma is a similarity parameter that is used to decide how similar two projected points are.  The similarity equation is

![Simliarity Equation](/images/similarity.png  "Similarity Equation")

Sigma, in the denominator of the exponent, determines how close projected points have to be to be considered similar. Further details are given in the paper.

### Lambda
Lambda is a parameter to determine how much smoothing we want.  Details of how this parameter is used is given in the paper. 

### Which Vector
The perturbation discovery boils down to an eigenvalue problem where the first eiegenvector is the perturbation that changes the projection the most. However, other eigenvectors may contain interesting information so the user can specify which perturbation vector they want to see (from 0 to n-1). By default, this is 0. 





 

# README

## Description
These are the mono-solver and poly-solver using quadrature based moment method (QMOM) for the droplets of different sizes moving with the same or different velocities. These are developed based on OpenFOAM v2106 (ESI-OpenCFD, 2021. ESI OpenCFD Release OpenFOAM® v2106. https://www.openfoam.com/news/main-news/openfoam-v2106) and OpenQBMM v7 (Passalacqua, A., Heylmun, J., Icardi, M., Madadi, E., Bachant, P., Hu, X., Weaver, J., 2021. OpenQBMM/OpenQBMM: OpenQBMM 7.0.0 for OpenFOAM v2106. Zenodo. doi:10.5281/ZENODO.5039574). 

**monodisperseDropletNOIFoam**: is for the droplets moving with the same velocity

**polydisperseDropletNOIFoam**: is for the droplets moving with the different velocities

## Installation

+ Install OpenFOAM v2106 and OpenQBMM v7 first
+ Copy the folder `myMultiphase` to the folder `OpenQBMM_7.0.0/applications/solvers`
+ Copy the folders `myPDFTransportModels` and `myPopulationBalanceModels` to the folder `OpenQBMM_7.0.0/src/quadratureMethods`
+ Use `wmake` to compile first in folder `myPDFTransportModels`
+ Use `wmake` to compile then in folder `myPopulationBalanceModels`
+ Use `Allwmake` to compile in folder `myMultiphase`
+ Modify the following file `myMultiphase\twoPhaseSystem\phaseModels\phaseModel\phaseModels.C`: delete the commenting symbols `/*` at line 96 and `*/` at line 172
+ Use `Allwmake` to compile again in folder `myMultiphase`

## Usage

### monodisperseDropletNOIFoam

#### Freely falling droplets

Please find the setups in the testing case of freely falling droplets in folder `testing_cases/freelyFallingDroplets_QMOM`, and run the simulation as follows:

+ blockMesh
+ setFields
+ monodisperseDropletNOIFoam

### polydisperseDropletNOIFoam

Please find the setups in the testing case of freely falling droplets in folder `testing_cases/aerosolsTransmissionChen_QMOM`, and run the simulation as follows:

+ blockMesh
+ polydisperseDropletNOIFoam

**Note**: To use the RNG $k-\epsilon$ model, please add the following sentences in the file `multiphaseCompressibleTurbulenceModels.C` located at `OpenFOAM-v2106/src/phaseSystemModels/reactingEuler/multiphaseSystem/turbulence` and compile it

```
#include "RNGkEpsilon.H"
makeRASModel(RNGkEpsilon);
```


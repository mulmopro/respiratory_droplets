# README

## Description
These are the solvers using the pseudo-single phase model (PSPM) for the small-size droplets or aerosols that are passively transported by the airflow. These are developed based on OpenFOAM v2106 (ESI-OpenCFD, 2021. ESI OpenCFD Release OpenFOAM® v2106. https://www.openfoam.com/news/main-news/openfoam-v2106). 

**buoyantSimpleCFoam**: is a steady solver which solves an additional transport equation for a scalar droplet concentration.

**multiSctBuoyantSimpleFoam**: is a steady solver which solves the species transport equations. The species mass fraction will affect the thermal physical properties of the gas.

## Installation

+ Install OpenFOAM v2106 first
+ Copy the folder `buoyantSimpleCFoam` and `multiSctBuoyantSimpleFoam` to the folder `OpenFOAM-v2106/applications/solvers/heatTransfer`
+ enter the folder `buoyantSimpleCFoam`, and compile with the command `wmake`
+ enter the folder `multiSctBuoyantSimpleFoam`, and compile with the command `wmake`

## Usage

### buoyantSimpleCFoam

Please find the setups in the testing case of transmission of tracer particles in an office under mixing ventilation in folder `testing_cases/officeMVLiu2022_PSPM`, and run the simulation as follows:

+ create mesh using `snappyHexMesh`
```
blockMesh
surfaceFeatureExtract
snnappyHexMesh -overwrite
```

+ run simulation parallelly

```
decomposePar
mpirun -np 8 buoyantSimpleCFoam -parallel > log.buoyantSimpleCFoam
```

+ post-process the results to obtain velocity (U), temperature (T), and CO2 concentration (CO2) along the sampling lines

```
reconstructPar
buoyantSimpleCFoam -postProcess -func "sample"
```

### multiSctBuoyantSimpleFoam

Please find the setups in the testing case of transmission of tracer gas in an office under dispalcement ventilation in folder `testing_cases/officeDVTian2019_PSPM`, and run the simulation as follows:

+ create mesh using `snappyHexMesh`
```
blockMesh
surfaceFeatureExtract
snnappyHexMesh -overwrite
```

+ run simulation parallelly

```
decomposePar
mpirun -np 8 multiSctBuoyantSimpleFoam -parallel > log.multiSctBuoyantSimpleFoam
```

+ post-process the results to obtain velocity (U), temperature (T), and CO2 concentration (CO2) along the sampling lines

```
reconstructPar
multiSctBuoyantSimpleFoam -postProcess -func "sample"
```






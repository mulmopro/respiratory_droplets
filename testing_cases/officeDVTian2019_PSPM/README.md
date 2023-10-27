# Transmission of tracer gas in an office under displacement ventilation

## Description

In the experiment of Tian et al. [^1], velocity, temperature, and CO2 concentration were measured in an office with displacement ventilation (DV), mixing ventilation (MV), and stratum ventilation (SV).

This testing case is to reproduce the case of DV with source located in armpit (z=1.0 m). The solver `multiSctBuoyantSimpleFoam` was adopted to solve the species transport equations for water vapor and CO2.

## Run

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

[^1]: Tian, X., Li, B., Ma, Y., Liu, D., Li, Y., & Cheng, Y. (2019). Experimental study of local thermal comfort and ventilation performance for mixing, displacement and stratum ventilation in an office. Sustainable Cities and Society, 50, 101630.
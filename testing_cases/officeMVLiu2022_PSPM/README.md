# Transmission of tracer particles in an office under mixing ventilation

## Description

In the experiment of Liu et al. [^1], velocity, temperature, and particle concentration were measured in an office with displacement ventilation (DV) and mixing ventilation (MV).

This testing case is to reproduce the case of MV. The solver `buoyantSimpleCFoam` was adopted to solve the concentration transport equations for the tracer particles.

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
mpirun -np 8 buoyantSimpleCFoam -parallel > log.buoyantSimpleCFoam
```

+ post-process the results to obtain velocity (U), temperature (T), and particle concentration (C) along the sampling lines

```
reconstructPar
buoyantSimpleCFoam -postProcess -func "sample"
```

[^1]: Liu, S., Koupriyanov, M., Paskaruk, D., Fediuk, G., & Chen, Q. (2022). Investigation of airborne particle exposure in an office with mixing and displacement ventilation. Sustainable Cities and Society, 79, 103718.

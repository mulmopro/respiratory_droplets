# README

## Description
This is the SM-solver using sectional method for the droplets of different sizes moving with the same or different velocities. This is developed based on the OpenFOAM v8 (OpenFOAM-Foundation, 2020. OpenFOAM 8 Released. https://openfoam.org/release/8/). 

A new `phaseSystem` - `interfaceCompositionPhaseChangePopulationBalanceMultiphaseSystem` is added to the standard solver `multiphaseEulerFoam`, which couples the phase system for phase change with that for population balance equation.

## Installation

+ Install OpenFOAM v8 first
+ Copy the folder `PBEInterfaceCompositionPhaseChangePhaseSystem` to the folder `OpenFOAM-8/applications/solvers/multiphase/multiphaseEulerFoam/phaseSystems/PhaseSystems`
+ Copy the file `multiphaseSystems.C` to the folder `OpenFOAM-8/applications/solvers/multiphase/multiphaseEulerFoam/multiphaseEulerFoam/multiphaseSystems` (please make backup first)
+ Enter the folder `OpenFOAM-8/applications/solvers/multiphase/multiphaseEulerFoam/multiphaseEulerFoam/multiphaseSystems` and compile it using `wmake`

## Usage

Please find the setups in the testing case of freely falling droplets in folder `testing_cases/freelyFallingDroplets_SM`, and run the simulation as follows:

+ blockMesh
+ setFields
+ multiphaseEulerFoam

## Postprocess

The minimum and maximum value of droplets diameter (`d.liquid`) will be extracted during the simulation and stored in folder `freelyFallingDroplets_SM/postProcessing/fieldMinMax`.


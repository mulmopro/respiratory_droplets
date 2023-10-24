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

### multiSctBuoyantSimpleFoam




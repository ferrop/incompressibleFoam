# incompressibleFoam

incompressibleFoam is a single phase solver (incompressibleFoam) with 2 forms of Rhie & Chow interpolation and 2 pressure formulations

## Rhie & Chow interpolation
### Consistent momentum interpolation (CMI)
The CMI allows to have time-step and relaxation factors independant results. This approach is very robust and avoid pressure-velocity decoupling. 
However the method is dissipative. 
In the CMI we use the flux at face rather than an interpolation of the velocity from cell centre to face. 
### Non-consistent momentum interpolation (NCMI)

## Pressure formulation
### Standard form
The solved pressure is the physical pressure (divided by the density)
### Corrected form
The solved pressure is the corrected pressure $p_c$ such as :
$
p = p^0 + p_c
$

## Installation

# How to compile ?
The code compiles for OpenFOAM ESI 2412+ (wwww.openfoam.com)
After loading OpenFOAM environnement 
```bash
wmake all

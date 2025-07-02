# incompressibleFoam

incompressibleFoam is a single phase solver (incompressibleFoam) with two forms of the Rhie & Chow interpolation and two pressure formulations.
The solver features are :
- Seven time schemes, from steady-state to BDF and (E)SDIRK up to third order, have been coded.
- Two forms of Rhie & Chow interpolation (CMI and NCMI)
- Two forms of pressure (standard and corrected)
- 2nd order field extrapolation
- Correct handling of pressure relaxation

Limitations :
- Moving grids are currently not supported
- The proposed time schemes are only available for ddt(U)
- If a turbulence model is used, the standard OpenFOAM schemes have to be used for the turbulence variables (ddt(k), ddt(omega) etc..). The implementation supports an (E)SDIRK scheme for ddt(U) and an OpenFOAM's standard scheme for turbulence. 
  However, this is rather unorthodox and only steadyState, localEuler, Euler, backward or CrankNicolson should be used with turbulence modeling.

## Contact
Paulin FERRO : paulin.ferro@g-met.fr

# References
P. Ferro, P. Landel, C. Landrodie, and M. Pescheux. incompressiblefoam: a new time consistent framework with bdf
and dirk integration schemes. 2024. URL https://arxiv.org/abs/2411.08688.

## Available $\frac{d\textbf{u}}{dt}$ schemes
### Transient schemes
 | Name | Description |
 |----------------|-------------|
 | _Euler_ | 1st order BDF1  |
 | _backward_ | 2nd order BDF2  |
 | _BDF3_ | 3d order BDF3 |

 | Name | Description |
 |----------------|-------------|
 | _CrankNicolson_ | 1 stage, 2nd order (E)SDIRK |
 | _ESDIRK23_1_ | 2 stages, 3d order (E)SDIRK |
 | _ESDIRK23_2_ | 2 stages, 3d order (E)SDIRK |
 | _DIRK22_ | 2 stages, 2nd order SDIRK |
 | _DIRK33_ | 3 stages, 3d order SDIRK |
 | _DIRK43_ | 4 stages, 3d order SDIRK |

With the _CrankNicolson_ scheme, the off-centering coefficient is controlled as follows :
```cpp
    ddt(U)          CrankNicolson;
    ocCoeff         1.0;
````
If ocCoeff = 0, _CrankNicolson_ scheme is equivalent to a pure 1st order implicit Euler scheme. The full 2nd order scheme
is obtained with ocCoeff = 1.

### Steady-state schemes
 | Name | Description |
 |----------------|-------------|
 | _steadyState_ | steady scheme for SIMPLE mode calculations |
 | _localEuler_ | LTS schemes for pseudo-transient approaches |

For _steadyState_, the correct settings are :
- in the fvSolutions file: set nOuterCorrectors and nCorrectors to 1.
- for the residuals definition, the fields have to be completed with "Final" notation.
- for the linear solvers use also the notation with "Final".

Note that convergence based on _residualControl_ is currently not supported for steadyState scheme.
```cpp
solvers
{
    ...
    pFinal
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-08;
        relTol          0.01;
    }
    ...
}
PIMPLE
{
   ...
   nOuterCorrectors          1;
   nCorrectors               1;
   ...
}
relaxationFactors
{
    fields
    {
        pFinal     0.5;
    }
    equations
    {
        UFinal     0.5;
        kFinal     0.3;
        omegaFinal 0.3;
    }
}
````
## Rhie & Chow interpolation
### Consistent momentum interpolation (CMI)
The CMI allows to avoid the dependency of time-step and relaxation factors. This approach is very robust and avoids pressure-velocity decoupling. 
However, it is dissipative. In the CMI, the old-time flux is used at face rather than an interpolation of the old-time velocity from cell center to face. The CMI is activated using the _consistentRhieChow_ entry (default value true).
```cpp
PIMPLE
{
   ...
    consistentRhieChow       true;
   ...
}
````
### Non-consistent momentum interpolation (NCMI)
The NCMI is less dissipative than the CMI. However, pressure-velocity decoupling (checkerboarding) might occurs at small time steps and the solution is time step and relaxation factors sensitive.
In the NCMI, old-time flux are calculated by linear interpolation of old-time velocities from cell center to face. 

## Pressure formulation
### Standard form
The solved pressure is the physical pressure $p$ (divided by the density). The standard pressure form is controlled by _pressureCorrectionForm_ entry (default value false).
```cpp
PIMPLE
{
   ...
    pressureCorrectionForm    false;
   ...
}
````
### Corrected form
The solved pressure is the corrected pressure $p_c$ such as :
$p = p^0 + p_c$ with $p^0$ the old-time or last iteration pressure (if steady scheme). With NCMI the corrected form will likely leads to strong checkerboard effects.

```cpp
PIMPLE
{
   ...
    pressureCorrectionForm    true;
   ...
}
````
Note that only fixedValue and zeroGradient boundary conditions are supported for $p$ if the pressure corrected form is chosen.
## Field extrapolation
This option allows to start the first temporal iteration with fields (flux, velocity and pressure) calculated using a 2nd order extrapolation. 
This option is particulary interesting for PISO transient simulations. The activation of extrapolation will likely constrain the maximum Courant number.

$\boldsymbol{u}^{n+1}_P = \frac{d t}{d t_0}(\boldsymbol{u}^{n}_P-\boldsymbol{u}^{n-1}_P)+\boldsymbol{u}^{n}_P$

This option is not relevant for _steadyState_ scheme.

The OpenFOAM behavior is a first order extrapolation.
# Installation

## How to compile ?
The code compiles for OpenFOAM ESI 2412+ (wwww.openfoam.com)
After loading OpenFOAM environnement : 
```bash
cd incompressibleFoam
wmake all
```
# Licence
incompressibleFoam is published under the GNU GPL Version 3 license, which can be found in the Licence file.

# Disclaimer
The use of this code is at your own risk. The authors, contributors, and copyright holders of this project are not liable for any direct or indirect damages resulting from the use or inability to use this code.
This code is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and non-infringement. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in an action of contract, tort, or otherwise, arising from, out of, or in connection with the code or the use or other dealings in the code.
Users are encouraged to review the code carefully and ensure it is suitable for their needs before using it. Users must also ensure they comply with all applicable laws and regulations in their jurisdiction when using this code.
By using this code, you acknowledge that you have read this disclaimer, understood its contents, and agree to its terms.

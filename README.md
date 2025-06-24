# incompressibleFoam

incompressibleFoam is a single phase solver (incompressibleFoam) with two forms of the Rhie & Chow interpolation and two pressure formulations.
The solver features are :
- Seven time schemes, from steady-state to BDF and (E)SDIRK up to third order, have been coded.
- Two forms of Rhie & Chow interpolation (CMI and NCMI)
- Two forms of pressure (standard and corrected)
- 2nd order field extrapolation
- Correct handling of pressure relaxation

# References
P. Ferro, P. Landel, C. Landrodie, and M. Pescheux. incompressiblefoam: a new time consistent framework with bdf
and dirk integration schemes. 2024. URL https://arxiv.org/abs/2411.08688.

## Available time schemes
### Transient schemes
 | Name | Description | Name | Description |
 |----------------|-------------|-------------|-------------|
 | Euler | 1st order BDF1  | DIRK22 | 2 stages 2nd order SDIRK |
 | backward | 2nd order BDF2  | DIRK33 | 3 stages 3d order SDIRK |
 | BDF3 | 3d order BDF3  | DIRK43 | 4 stages 3d order SDIRK |
 | CrankNicolson | 2nd order (E)SDIRK  | ESDIRK23_1 | 2 stages 3d order (E)SDIRK |
 | ESDIRK23_2 | 2 stages 3d order (E)SDIRK |  | |

With the Crank-Nicolson scheme, the off-centering coefficient is controlled as follows :
```cpp
    ddt(U)          CrankNicolson;
    ocCoeff         1.0;
````
If ocCoeff = 0, CrankNicolson is equivalent to a pure 1st order implicit Euler scheme. The full 2nd order scheme
is obtained with ocCoeff = 1.

### Steady-state schemes
 | Name | Description |
 |----------------|-------------|
 | steadyState | steady scheme for SIMPLE mode calculations |
 | localEuler | LTS schemes for pseudo-transient approaches |

For steadyState, the correct settings in the fvSolutions file, is to set nOuterCorrectors and nCorrectors to 1 and to control relaxation with the final fields.
fieldsExtrapolation has to be set to false. Note that convergence based on residualControl is currently not available with the steadyState scheme.
```cpp
PIMPLE
{
   ...
   fieldsExtrapolation      false;
   nOuterCorrectors          1;
   nCorrectors               1;
   ...
}
relaxationFactors
{
    fields
    {
        pFinal    0.5; //<--relaxation factor entry for the steadyState scheme
    }
    equations
    {
        UFinal    0.5; //<--relaxation factor entry for the steadyState scheme
    }
}
````
## Rhie & Chow interpolation
### Consistent momentum interpolation (CMI)
The CMI allows to avoid the dependency of time-step and relaxation factors. This approach is very robust and avoids pressure-velocity decoupling. 
However, it is dissipative. 
In the CMI, the old-time flux is used at face rather than an interpolation of the old-time velocity from cell center to face. The CMI is activated using the following consistentRhieChow entry (default value true).
```cpp
PIMPLE
{
   ...
    consistentRhieChow       true;
   ...
}
````
### Non-consistent momentum interpolation (NCMI)
The NCMI is less dissipative than the CMI. However, pressure-velocity decoupling (checkerboarding) might occurs at small time steps.
In the NCMI, old-time flux are calculated by linear interpolation of old-time velocities from cell center to face. 

## Pressure formulation
### Standard form
The solved pressure is the physical pressure $p$ (divided by the density). The standard pressure form is controlled by pressureCorrectionForm entry (default value false).
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
## Field extrapolations
This option allows to start the first temporal iteration with fields (flux, velocity and pressure) calculated using a 2nd order extrapolation. 
This option is particulary interesting for PISO transient simulations. The activation of extrapolation will likely constrains the maximal Courant number.

$\boldsymbol{u}^{n+1}_P = \frac{d t}{d t_0}(\boldsymbol{u}^{n}_P-\boldsymbol{u}^{n-1}_P)+\boldsymbol{u}^{n}_P$

The OpenFOAM behavior is a first order extrapolation.
# Installation

## How to compile ?
The code compiles for OpenFOAM ESI 2412+ (wwww.openfoam.com)
After loading OpenFOAM environnement : 
```bash
wmake all
```
# Licence
incompressibleFoam is published under the GNU GPL Version 3 license, which can be found in the LICENSE file.

# Disclaimer
The use of this code is at your own risk. The authors, contributors, and copyright holders of this project are not liable for any direct or indirect damages resulting from the use or inability to use this code.
This code is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and non-infringement. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in an action of contract, tort, or otherwise, arising from, out of, or in connection with the code or the use or other dealings in the code.
Users are encouraged to review the code carefully and ensure it is suitable for their needs before using it. Users must also ensure they comply with all applicable laws and regulations in their jurisdiction when using this code.
By using this code, you acknowledge that you have read this disclaimer, understood its contents, and agree to its terms.

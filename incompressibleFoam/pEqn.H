// We build all necessary fields before Rhie & Chow interpolation
volScalarField as( convDiffUEqn.A() );

if (ddtUEqn.hasDiag()) //Necessary check for steadyState scheme
{
    at = ddtUEqn.A();
}

volScalarField atild( at + as*a[i][i] );
surfaceScalarField    atf(linearInterpolate(at));
surfaceScalarField atildf(linearInterpolate(atild));
surfaceScalarField rAUfeq("rAUfeq", ci*alphaU / atildf);

volVectorField convDiffH( a[i][i]*convDiffUEqn.H() );

//Rhie & Chow interpolation 1/2 part:
//Convection/diffusion flux and other sources
surfaceScalarField phiRC
(
     "phiRC",
      alphaU * fvc::flux(convDiffH) / atildf
    + alphaU * sumAijRUfj / atildf
);

if(pressureCorrectionForm)
{
    if (consistentRhieChow)
    {
        phiRC -= fvc::snGrad(p.oldTime(),"fvc::snGrad(p)")*mesh.magSf()*rAUfeq;
    }
    else
    {
        phiRC -= fvc::flux(fvc::grad(p.oldTime(),"fvc::grad(p)"))*rAUfeq;
    }

}
//Rhie & Chow interpolation 2/2 part:
//Old time and previous iteration contribution
surfaceScalarField phiOldAndRelax
(
     timeMethod->phiOldAndRelax(U,phi,atf,atildf,alphaU,consistentRhieChow)
);

phiRC += phiOldAndRelax;

//Here we constrain the flux of Rhie & Chow interpolation
//to be consistent with velocity boundary conditions
//In standard OF solvers this is done by constrainHbyA function.
timeMethod->constrainFlux(phiRC,U,pressureCorrectionForm ? pc : p, MRF);
MRF.makeRelative(phiRC);

if (pressureCorrectionForm ? pc.needReference() : p.needReference())
{
    adjustPhi(phiRC, U, pressureCorrectionForm ? pc : p);
}

// Update the pressure BCs to ensure flux consistency
// This function is needed only if a fixedFluxPressure BC is used
constrainPressure(pressureCorrectionForm ? pc : p, U, phiRC, rAUfeq, MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUfeq, pressureCorrectionForm ? pc : p,"laplacian(rAUf,p)") == fvc::div(phiRC)
    );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiRC - pEqn.flux();
    }
}
if (pimple.finalIter())
{
    #include "continuityErrs.H"
}

//We update the cell centre velocity field with all contributions
U = timeMethod->updateU(i,U,sumAijRUj,convDiffH,ddtUEqn,atild,p,alphaU);
U.correctBoundaryConditions();
fvOptions.correct(U);

// Explicitly relax pressure
if(pressureCorrectionForm)
{
    p = p.oldTime() + alphaP*pc;
}
else
{
    p.relax(alphaP);
}

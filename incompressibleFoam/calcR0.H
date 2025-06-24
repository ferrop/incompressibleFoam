Info<< "Explicit step " << endl;

MRF.correctBoundaryVelocity(U);
MRF.makeRelative(phi);

fvVectorMatrix convDiffUEqn
(
     fvm::div(phi, U)
   + turbulence->divDevReff(U)
   + MRF.DDt(U)
);

MRF.makeAbsolute(phi);

RU[0]()  = - convDiffUEqn.A()*U + convDiffUEqn.H() ;

if (consistentRhieChow)
{
    RUf[0]() = - linearInterpolate(convDiffUEqn.A())*phi + fvc::flux(convDiffUEqn.H()) ;
}
else
{
    RUf[0]() = fvc::flux(RU[0]());
}

MRF.correctBoundaryVelocity(U);

// Solve the Momentum equation
fvVectorMatrix ddtUEqn
(
    timeMethod->ddt(U)
);

fvVectorMatrix convDiffUEqn
(
    fvm::div(phi, U)
  + turbulence->divDevReff(U)
  + MRF.DDt(U)
     ==
    fvOptions(U)
);

fvOptions.constrain(convDiffUEqn);

fvVectorMatrix UEqn
(
    ddtUEqn
  + a[i][i] * convDiffUEqn
);

if (pimple.momentumPredictor())
{
    UEqn.relax(alphaU);

    solve
    (
        UEqn == - ci*fvc::grad(p) + sumAijRUj
    );

    fvOptions.correct(U);
}

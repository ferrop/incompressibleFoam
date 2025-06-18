scalar alphaU(1);
scalar alphaP(1);

if (pimple.finalIter())
{
    alphaU = mesh.equationRelaxationFactor(U.name() + "Final");
    alphaP = mesh.fieldRelaxationFactor(p.name() + "Final");
}
else
{
    alphaU = mesh.equationRelaxationFactor(U.name());
    alphaP = mesh.fieldRelaxationFactor(p.name());
}

Info << "Relaxing UEqn with factor : " << alphaU << endl;
Info << "Relaxing p with factor : " << alphaP << endl;

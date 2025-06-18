List<scalarList> a = timeMethod->a();
int nRKloop = a.size();

List<autoPtr<volVectorField>>     RU (nRKloop);
List<autoPtr<surfaceScalarField>> RUf (nRKloop);

forAll (RU, iField)
{
    word nameRi ("RU" + name(iField));

    RU[iField].reset
    (
        new volVectorField
        (
            IOobject
            (
                nameRi,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("0",U.dimensions()/dimTime,Zero)/*,
            fvPatchFieldBase::zeroGradientType()*/
        )
    );
}

forAll (RUf, iField)
{
    word nameRi ("R" + name(iField));

    RUf[iField].reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                nameRi,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0",phi.dimensions()/dimTime,Zero)
        )
    );
    RUf[iField]().setOriented(true);
}

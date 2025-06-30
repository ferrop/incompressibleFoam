Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating pressure correction field\n" << endl;
volScalarField pc
(
    IOobject
    (
        "pc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0",p.dimensions(),0.0),
    p.boundaryField().types()
);

//The pressure correction boundary conditions is 0 when fixedValue is used for p
typename GeometricField<scalar, fvPatchField, volMesh>::Boundary pcBF = pc.boundaryField();
forAll(pc.boundaryField(), patchI)
{
    if( isType<fixedValueFvPatchScalarField>(pc.boundaryField()[patchI]))
    {
        pcBF[patchI] == scalarField(pcBF.size(),scalar(0.0));
    }
}

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
mesh.setFluxRequired(pc.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

autoPtr<timeSchemes> timeMethod
(
    Foam::timeSchemes::New(mesh)
);

IOdictionary fvSchemes
(
    IOobject
    (
        "fvSchemes",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
      )
);
word timeScheme = fvSchemes.subDict("ddtSchemes").get<word>("ddt(U)");

volScalarField at
(
    IOobject
    (
        "at",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0",dimless/dimTime,Zero)
);

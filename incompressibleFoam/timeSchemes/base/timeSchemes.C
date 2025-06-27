// ************************************************************************* //

#include "timeSchemes.H"
#include "fixedFluxExtrapolatedPressureFvPatchScalarField.H"
#include "EulerDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeSchemes, 0);
    defineRunTimeSelectionTable
    (
        timeSchemes,
        components
    );
}

Foam::timeSchemes::timeSchemes
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
           "fvSchemes",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeSchemes::~timeSchemes()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar Foam::timeSchemes::c
(
    int i
)
{
    scalar ci(0);

    for (int j = 0; j < (i+1) ; j++)
    {
        ci += this->a()[i][j];
    }
    return ci;
}

volVectorField Foam::timeSchemes::sumAijRj
(
    int i,
    const List<autoPtr<volVectorField>>& RiPtr
)
{
    volVectorField sumAijRj
    (
        IOobject
        (
            "sumAijRj",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0",RiPtr[0]().dimensions(),Zero)
    );

    for (int j = 0; j < i ; j++)
    {
        sumAijRj += this->a()[i][j]*RiPtr[j]();
    }

    return sumAijRj;
}

volScalarField Foam::timeSchemes::sumAijRj
(
    int i,
    const List<autoPtr<volScalarField>>& RiPtr
)
{
    volScalarField sumAijRj
    (
        IOobject
        (
            "sumAijRj",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0",RiPtr[0]().dimensions(),Zero)
    );

    for (int j = 0; j < i ; j++)
    {
        sumAijRj += this->a()[i][j]*RiPtr[j]();
    }

    return sumAijRj;
}
surfaceVectorField Foam::timeSchemes::sumAijRfj
(
    int i,
    const List<autoPtr<surfaceVectorField>>& RfiPtr,
    bool orient
)
{
    surfaceVectorField sumAijRfj
    (
        IOobject
        (
            "sumAijRfj",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0",RfiPtr[0]().dimensions(),Zero)
    );

    sumAijRfj.setOriented(orient);

    for (int j = 0; j < i ; j++)
    {
        sumAijRfj += this->a()[i][j]*RfiPtr[j]();
    }

    return sumAijRfj;
}
surfaceScalarField Foam::timeSchemes::sumAijRfj
(
    int i,
    const List<autoPtr<surfaceScalarField>>& RfiPtr,
    bool orient
)
{
    surfaceScalarField sumAijRfj
    (
        IOobject
        (
            "sumAijRfj",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0",RfiPtr[0]().dimensions(),Zero)
    );

    sumAijRfj.setOriented(orient);

    for (int j = 0; j < i ; j++)
    {
        sumAijRfj += this->a()[i][j]*RfiPtr[j]();
    }

    return sumAijRfj;
}

fvVectorMatrix Foam::timeSchemes::ddt
(
    const volVectorField& vf
)
{
    return fv::EulerDdtScheme<vector>(mesh_).fvmDdt(vf);
}

fvScalarMatrix Foam::timeSchemes::ddt
(
    const volScalarField& vf
)
{
    return fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(vf);
}

surfaceScalarField Foam::timeSchemes::phiOldAndRelax
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& atf,
    const surfaceScalarField& atildf,
    const scalar& alphaU,
    const bool& consistentRhieChow
)
{
    if (consistentRhieChow)
    {
        return
        (
                   alphaU *atf*phi.oldTime() / atildf
            + (1.0-alphaU)*phi.prevIter()
        );
    }
    else
    {
        return
        (
                   alphaU *atf*fvc::flux(U.oldTime()) / atildf
            + (1.0-alphaU)*phi.prevIter()
        );
    }
}

void Foam::timeSchemes::constrainFlux
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& p,
    const IOMRFZoneList& MRF
)
{
    // This is an adaptation of constrainHbyA fonction
    // We constrain the flux of HbyA instead of HbyA

    typename GeometricField<scalar, fvsPatchField, surfaceMesh>::Boundary& phib = phi.boundaryFieldRef();
    typename GeometricField<vector, fvsPatchField, surfaceMesh>::Boundary Sfb = mesh_.Sf().boundaryField();
    typename GeometricField<vector, fvPatchField, volMesh>::Boundary Ub = U.boundaryField();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
        !U.boundaryField()[patchi].assignable()
        && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
        (
            p.boundaryField()[patchi]
        )
        )
        {
            phib[patchi] = (Ub[patchi] & Sfb[patchi]);
        }
    }
}

void Foam::timeSchemes::extrapolateFields
(
    int i,
    surfaceScalarField& phi,
    volVectorField& U,
    volScalarField& p
)
{
    if (mesh_.time().timeIndex() > 2 ) // only extrapolate if data are available
    {
        scalar dt_dt0(this->c(i)*mesh_.time().deltaTValue()/mesh_.time().deltaT0Value());
        phi = dt_dt0*(phi.oldTime() - phi.oldTime().oldTime()) + phi.oldTime() ;
        U   = dt_dt0*(  U.oldTime() -   U.oldTime().oldTime()) +   U.oldTime() ;
        p   = dt_dt0*(  p.oldTime() -   p.oldTime().oldTime()) +   p.oldTime() ;
    }
    else
    {
        phi = phi.oldTime() ;
        U   = U.oldTime()   ;
        p   = p.oldTime()   ;
    }
}

volVectorField Foam::timeSchemes::updateU
(
    const int& i,
    const volVectorField& U,
    const volVectorField& sumAijRUj,
    const volVectorField& convDiffH,
    const fvVectorMatrix& ddtUEqn,
    const volScalarField& atild,
    const volScalarField& p,
    const scalar& alphaU
)
{
    return
    (
         alphaU
         *
         (
            ddtUEqn.H()
          - fvc::grad(p)*this->c(i)
          + sumAijRUj
          + convDiffH
         ) / atild
         + (1.0-alphaU)*U.prevIter()
    );
};

bool Foam::timeSchemes::read()
{
    if (regIOobject::read())
    {
        return true;
    }

    return false;
}

// ************************************************************************* //

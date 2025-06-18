
#include "backward.H"
#include "addToRunTimeSelectionTable.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(backward, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        backward,
        components
    );

}

Foam::backward::backward
(
    const fvMesh& mesh
)
:
    timeSchemes(mesh)
{
    a_.setSize(1);
    a_[0].setSize(1);

    a_[0][0] = 1;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

List<scalarList> Foam::backward::a
()
{
    return(a_);
}

volVectorField Foam::backward::sumAijRj
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
            RiPtr[0]().mesh().time().timeName(),
            RiPtr[0]().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        RiPtr[0]().mesh(),
        dimensionedVector("0",RiPtr[0]().dimensions(),Zero)
    );

    return sumAijRj;
}

volScalarField Foam::backward::sumAijRj
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
            RiPtr[0]().mesh().time().timeName(),
            RiPtr[0]().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        RiPtr[0]().mesh(),
        dimensionedScalar("0",RiPtr[0]().dimensions(),Zero)
    );

    return sumAijRj;
}

surfaceScalarField Foam::backward::sumAijRfj
(
    int i,
    const List<autoPtr<surfaceScalarField>>& RfiPtr
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

    sumAijRfj.setOriented(true);

    return sumAijRfj;
}

fvVectorMatrix Foam::backward::ddt
(
    const volVectorField& vf
)
{
    return fv::backwardDdtScheme<vector>(mesh_).fvmDdt(vf);
}

fvScalarMatrix Foam::backward::ddt
(
    const volScalarField& vf
)
{
    return fv::backwardDdtScheme<scalar>(mesh_).fvmDdt(vf);
}

surfaceScalarField Foam::backward::phiOldAndRelax
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& atf,
    const surfaceScalarField& atildf,
    const scalar& alphaU,
    const bool& consistentRhieChow
)
{
    scalar deltaT = backward::deltaT_();
    scalar deltaT0 = backward::deltaT0_(mesh_);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (consistentRhieChow)
    {
        return
        (
              alphaU *atf*phi.oldTime()/atildf*coefft0/coefft
            - alphaU *atf*phi.oldTime().oldTime()/atildf*coefft00/coefft
            + (1.0-alphaU)*phi.prevIter()
        );
    }
    else
    {
        return
        (
               alphaU *atf*fvc::flux(U.oldTime())/atildf*coefft0/coefft
             - alphaU *atf*fvc::flux(U.oldTime().oldTime())/atildf*coefft00/coefft
            + (1.0-alphaU)*phi.prevIter()
        );
    }
}

scalar backward::deltaT_() const
{
    return mesh_.time().deltaTValue();
}

scalar backward::deltaT0_() const
{
    return mesh_.time().deltaT0Value();
}

scalar backward::deltaT0_(const fvMesh& mesh) const
{
    if (mesh_.time().timeIndex() < 2)
    {
        return GREAT;
    }
    else
    {
        return deltaT0_();
    }
}

void Foam::backward::message(int i)
{
    Info<<"TIME SCHEME : backward"<<endl;
}

bool Foam::backward::read()
{
    if (read())
    {
        return true;
    }

    return false;
}

// ************************************************************************* //

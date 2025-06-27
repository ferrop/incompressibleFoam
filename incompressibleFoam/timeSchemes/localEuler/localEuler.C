// ************************************************************************* //

#include "localEuler.H"
#include "localEulerDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(localEuler, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        localEuler,
        components
    );

}

Foam::localEuler::localEuler
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

List<scalarList> Foam::localEuler::a
()
{
    return(a_);
}

volVectorField Foam::localEuler::sumAijRj
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

volScalarField Foam::localEuler::sumAijRj
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

surfaceScalarField Foam::localEuler::sumAijRfj
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

fvVectorMatrix Foam::localEuler::ddt
(
    const volVectorField& vf
)
{
    return fv::localEulerDdtScheme<vector>(mesh_).fvmDdt(vf);
}

fvScalarMatrix Foam::localEuler::ddt
(
    const volScalarField& vf
)
{
    return fv::localEulerDdtScheme<scalar>(mesh_).fvmDdt(vf);
}

surfaceScalarField Foam::localEuler::phiOldAndRelax
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& atf,
    const surfaceScalarField& atildf,
    const scalar& alphaU,
    const bool& consistentRhieChow
)
{
    volScalarField rDeltaT = fv::localEulerDdt::localRDeltaT(mesh_) ;
    if (consistentRhieChow)
    {
        return
        (
                   alphaU *fvc::interpolate(rDeltaT)*phi.oldTime() / atildf
            + (1.0-alphaU)*phi.prevIter()
        );
    }
    else
    {
        return
        (
                   alphaU *fvc::interpolate(rDeltaT)*fvc::flux(U.oldTime()) / atildf
            + (1.0-alphaU)*phi.prevIter()
        );
    }
}

void Foam::localEuler::message(int i)
{
    Info<<"TIME SCHEME : localEuler"<<endl;
}

bool Foam::localEuler::read()
{
    if (read())
    {
        return true;
    }

    return false;
}
// ************************************************************************* //

// ************************************************************************* //

#include "Euler.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Euler, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        Euler,
        components
    );

}

Foam::Euler::Euler
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

List<scalarList> Foam::Euler::a
()
{
    return(a_);
}

volVectorField Foam::Euler::sumAijRj
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

volScalarField Foam::Euler::sumAijRj
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

surfaceScalarField Foam::Euler::sumAijRfj
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

void Foam::Euler::message(int i)
{
    Info<<"TIME SCHEME : Euler"<<endl;
}

bool Foam::Euler::read()
{
    if (read())
    {
        return true;
    }

    return false;
}
// ************************************************************************* //

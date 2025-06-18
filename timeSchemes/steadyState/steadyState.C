// ************************************************************************* //

#include "steadyState.H"
#include "addToRunTimeSelectionTable.H"
#include "steadyStateDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(steadyState, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        steadyState,
        components
    );

}

Foam::steadyState::steadyState
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

List<scalarList> Foam::steadyState::a
()
{
    return(a_);
}

volVectorField Foam::steadyState::sumAijRj
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

volScalarField Foam::steadyState::sumAijRj
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

surfaceScalarField Foam::steadyState::sumAijRfj
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

fvVectorMatrix Foam::steadyState::ddt
(
    const volVectorField& vf
)
{
    return fv::steadyStateDdtScheme<vector>(mesh_).fvmDdt(vf);
}

fvScalarMatrix Foam::steadyState::ddt
(
    const volScalarField& vf
)
{
    return fv::steadyStateDdtScheme<scalar>(mesh_).fvmDdt(vf);
}

surfaceScalarField Foam::steadyState::phiOldAndRelax
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& atf,
    const surfaceScalarField& atildf,
    const scalar& alphaU,
    const word& fluxMethod
)
{
    return
    (
         (1.0-alphaU)*phi.prevIter()
    );
}

void Foam::steadyState::message(int i)
{
    const dictionary& fvSolution = mesh_.thisDb().lookupObject<IOdictionary>("fvSolution");
    dictionary pimpleDict(fvSolution.subDict("PIMPLE"));
    label nCorrPIMPLE(pimpleDict.lookupOrDefault<label>("nOuterCorrectors",1));
    label nCorrPISO(pimpleDict.lookupOrDefault<label>("nCorrectors",1));

    Info<<"TIME SCHEME : steadyState"<<endl;

    if (nCorrPIMPLE > 1)
    {
        Info<<"WARNING : with steadyState use nOuterCorrectors 1 ;"<<endl;
    }
    if (nCorrPISO > 1)
    {
        Info<<"WARNING : with steadyState use nCorrectors 1 ;"<<endl;
    }
}

bool Foam::steadyState::read()
{
    if (read())
    {
        return true;
    }

    return false;
}
// ************************************************************************* //

// ************************************************************************* //

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::timeSchemes> Foam::timeSchemes::New
(
     const fvMesh& mesh
)
{
    const IOdictionary dict
    (
        IOobject
        (
            "fvSchemes",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false // Do not register
          )
    );

    ITstream& is = dict.subDict("ddtSchemes").lookup("ddt(U)");
    word timeSchemes;
    is >> timeSchemes;

    Info<< "Time scheme : " << timeSchemes << endl;

    auto cstrIter = componentsConstructorTablePtr_->cfind(timeSchemes);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            timeSchemes,
            "timeSchemes",
            timeSchemes,
            *componentsConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    return
        autoPtr<Foam::timeSchemes>
        (
            cstrIter()(mesh)
        );
}

// ************************************************************************* //

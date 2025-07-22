// ************************************************************************* //

#include "CrankNicolson.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CrankNicolson, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        CrankNicolson,
        components
    );

}

Foam::CrankNicolson::CrankNicolson
(
    const fvMesh& mesh
)
:
    timeSchemes(mesh)
{
    ITstream& is = optionalSubDict("ddtSchemes").lookup("ddt(U)");
    word dummy;
    scalar ocCoeff;
    is >> dummy;
    if (!is.eof())
    {
        is >> ocCoeff;
    }
    if (ocCoeff < 0 || ocCoeff > 1)
    {
        FatalIOErrorInFunction(*this)
            << "Off-centreing coefficient = " << ocCoeff
            << " should be >= 0 and <= 1"
            << exit(FatalIOError);
    }

    a_.setSize(2);
    a_[0].setSize(1);
    a_[1].setSize(2);

    a_[0][0] = 0.0 ;
    a_[1][0] = ocCoeff/(ocCoeff + 1.0); a_[1][1] = 1.0/(ocCoeff + 1.0);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

List<scalarList> Foam::CrankNicolson::a
()
{
     return(a_);
}

void Foam::CrankNicolson::message(int i)
{
    Info<<"TIME SCHEME : CrankNicolson"<<" : ITERATION "<<i+1<<"/"<<a_.size()<<endl;
}

bool Foam::CrankNicolson::read()
{
    if (read())
    {
        return true;
    }

    return false;
}
// ************************************************************************* //

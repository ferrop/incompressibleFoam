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
    timeSchemes(mesh),
    ocCoeff_(optionalSubDict("ddtSchemes").get<scalar>("ocCoeff"))
{
    if (ocCoeff_ < 0 || ocCoeff_ > 1)
    {
        FatalIOErrorInFunction(*this)
            << "Off-centreing coefficient = " << ocCoeff_
            << " should be >= 0 and <= 1"
            << exit(FatalIOError);
    }

    a_.setSize(2);
    a_[0].setSize(1);
    a_[1].setSize(2);

    a_[0][0] = 0.0 ;
    a_[1][0] = ocCoeff_/(ocCoeff_ + 1.0); a_[1][1] = 1.0/(ocCoeff_ + 1.0);
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
        subDict("ddtSchemes").readEntry("ocCoeff", ocCoeff_);
        return true;
    }

    return false;
}
// ************************************************************************* //

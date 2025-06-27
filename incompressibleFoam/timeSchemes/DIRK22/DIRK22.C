// ************************************************************************* //

#include "DIRK22.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DIRK22, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        DIRK22,
        components
    );

}

Foam::DIRK22::DIRK22
(
    const fvMesh& mesh
)
:
    timeSchemes(mesh)
{
    a_.setSize(2);
    a_[0].setSize(1);
    a_[1].setSize(2);

    const scalar gamma(1.0-Foam::sqrt(2.0)/2.0);

    a_[0][0] = gamma;
    a_[1][0] = 1.0-gamma;     a_[1][1] = gamma;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

List<scalarList> Foam::DIRK22::a
()
{
    return(a_);
}

void Foam::DIRK22::message(int i)
{
    Info<<"TIME SCHEME : DIRK22"<<" : ITERATION "<<i+1<<"/"<<a_.size()<<endl;
}

bool Foam::DIRK22::read()
{
    if (read())
    {
        return true;
    }

    return false;
}

// ************************************************************************* //

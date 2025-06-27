// ************************************************************************* //

#include "EDIRK23_1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EDIRK23_1, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        EDIRK23_1,
        components
    );

}

Foam::EDIRK23_1::EDIRK23_1
(
    const fvMesh& mesh
)
:
    timeSchemes(mesh)
{
    a_.setSize(3);
    a_[0].setSize(1);
    a_[1].setSize(2);
    a_[2].setSize(3);

    a_[0][0] = 0.0;
    a_[1][0] = 7./12.;   a_[1][1] = 7./12.;
    a_[2][0] = 5./14.;   a_[2][1] = -6./7.;    a_[2][2] = 3./2.;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

List<scalarList> Foam::EDIRK23_1::a
()
{
    return(a_);
}

void Foam::EDIRK23_1::message(int i)
{
    Info<<"TIME SCHEME : EDIRK23_1"<<" : ITERATION "<<i+1<<"/"<<a_.size()<<endl;
}

bool Foam::EDIRK23_1::read()
{
    if (read())
    {
        return true;
    }

    return false;
}
// ************************************************************************* //

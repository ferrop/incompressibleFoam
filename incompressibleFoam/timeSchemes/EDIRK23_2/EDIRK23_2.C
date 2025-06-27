// ************************************************************************* //

#include "EDIRK23_2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EDIRK23_2, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        EDIRK23_2,
        components
    );

}

Foam::EDIRK23_2::EDIRK23_2
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
    a_[1][0] = 3./4.;    a_[1][1] = 3./4.;
    a_[2][0] = 7./18.;   a_[2][1] = -4./18.;    a_[2][2] = 15./18.;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

List<scalarList> Foam::EDIRK23_2::a
()
{
    return(a_);
}

void Foam::EDIRK23_2::message(int i)
{
    Info<<"TIME SCHEME : EDIRK23_2"<<" : ITERATION "<<i+1<<"/"<<a_.size()<<endl;
}

bool Foam::EDIRK23_2::read()
{
    if (read())
    {
        return true;
    }

    return false;
}
// ************************************************************************* //

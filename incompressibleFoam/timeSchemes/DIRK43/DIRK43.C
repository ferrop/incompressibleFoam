// ************************************************************************* //

#include "DIRK43.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DIRK43, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        DIRK43,
        components
    );
}

Foam::DIRK43::DIRK43
(
    const fvMesh& mesh
)
:
    timeSchemes(mesh)
{
    a_.setSize(4);
    a_[0].setSize(1);
    a_[1].setSize(2);
    a_[2].setSize(3);
    a_[3].setSize(4);

    a_[0][0] = 0.5;

    a_[1][0] = 1.0/6.0;
    a_[1][1] = 0.5;

    a_[2][0] = -0.5;
    a_[2][1] = 0.5;
    a_[2][2] = 0.5;

    a_[3][0] = 3.0/2.0;
    a_[3][1] = -3.0/2.0;
    a_[3][2] = 0.5;
    a_[3][3] = 0.5;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

List<scalarList> Foam::DIRK43::a
()
{
    return(a_);
}

void Foam::DIRK43::message(int i)
{
    Info<<"TIME SCHEME : DIRK43"<<" : ITERATION "<<i+1<<"/"<<a_.size()<<endl;
}

bool Foam::DIRK43::read()
{
    if (read())
    {
        return true;
    }

    return false;
}
// ************************************************************************* //

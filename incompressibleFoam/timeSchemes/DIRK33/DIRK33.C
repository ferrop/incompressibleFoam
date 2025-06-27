// ************************************************************************* //

#include "DIRK33.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DIRK33, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        DIRK33,
        components
    );

}

Foam::DIRK33::DIRK33
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

    scalar gamma(0.43586652150845899941601945);
    a_[0][0] = gamma;

    a_[1][0] = (1+gamma)/2. - gamma;
    a_[1][1] = gamma;

    a_[2][0] = -((6*gamma*gamma-16.*gamma+1.)/4.);
    a_[2][1] =  ((6*gamma*gamma-20.*gamma+5.)/4.);
    a_[2][2] = gamma;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

List<scalarList> Foam::DIRK33::a
()
{
    return(a_);
}

void Foam::DIRK33::message(int i)
{
    Info<<"TIME SCHEME : DIRK33"<<" : ITERATION "<<i+1<<"/"<<a_.size()<<endl;
}

bool Foam::DIRK33::read()
{
    if (read())
    {
        return true;
    }

    return false;
}
// ************************************************************************* //

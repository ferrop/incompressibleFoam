// ************************************************************************* //

#include "BDF3.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BDF3, 0);

    addToRunTimeSelectionTable
    (
        timeSchemes,
        BDF3,
        components
    );

}

Foam::BDF3::BDF3
(
    const fvMesh& mesh
)
:
    timeSchemes(mesh)
{
/*
    a_.setSize(1);
    a_[0].setSize(1);
    a_[0][0] = 1;
*/
/*
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
*/

    a_.setSize(2);
    a_[0].setSize(1);
    a_[1].setSize(2);

    const scalar gamma(1.0-Foam::sqrt(2.0)/2.0);
    a_[0][0] = gamma;
    a_[1][0] = 1.0-gamma;     a_[1][1] = gamma;


}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

List<scalarList> Foam::BDF3::a
()
{
    if (mesh_.time().timeIndex() < 3)
    {
        return(a_);
    }
    else
    {
        a_.setSize(1);
        a_[0].setSize(1);
        a_[0][0] = 1;

        return(a_);
    }
}

volVectorField Foam::BDF3::sumAijRj
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

    if (mesh_.time().timeIndex() < 3)
    {
        for (int j = 0; j < i ; j++)
        {
            sumAijRj += this->a()[i][j]*RiPtr[j]();
        }
    }

    return sumAijRj;
}

volScalarField Foam::BDF3::sumAijRj
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

    if (mesh_.time().timeIndex() < 3)
    {
        for (int j = 0; j < i ; j++)
        {
            sumAijRj += this->a()[i][j]*RiPtr[j]();
        }
    }
    return sumAijRj;
}

surfaceScalarField Foam::BDF3::sumAijRfj
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

    if (mesh_.time().timeIndex() < 3)
    {
        for (int j = 0; j < i ; j++)
        {
            sumAijRfj += this->a()[i][j]*RfiPtr[j]();
        }
    }
    return sumAijRfj;
}

fvVectorMatrix Foam::BDF3::ddt
(
    const volVectorField& vf
)
{
    scalar a1(double(11./6.));
    scalar a2(double(18./6.));
    scalar a3(double(9./6.));
    scalar a4(double(2./6.));

    if (mesh_.time().timeIndex() < 3)
    {
        a1 = 1.0;
        a2 = 1.0;
        a3 = 0.0;
        a4 = 0.0;
    }

    tmp<fvVectorMatrix> tfvm
    (
        new fvVectorMatrix
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvVectorMatrix& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/deltaT_();

    fvm.diag() = (a1*rDeltaT)*mesh_.V();

    if (mesh_.moving())
    {/*
        fvm.source() = rDeltaT*
        (
            a2*vf.oldTime().primitiveField()*this->V().oldTime()
          - a3*vf.oldTime().oldTime().primitiveField()*this->V().oldTime().oldTime()
          + a4*vf.oldTime().oldTime().oldTime().primitiveField()*this->V().oldTime().oldTime().oldTime()
        );*/
    }
    else
    {
        fvm.source() = rDeltaT*mesh_.V()*
        (
            a2*vf.oldTime().primitiveField()
          - a3*vf.oldTime().oldTime().primitiveField()
          + a4*vf.oldTime().oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}

fvScalarMatrix Foam::BDF3::ddt
(
    const volScalarField& vf
)
{
    scalar a1(double(11./6.));
    scalar a2(double(18./6.));
    scalar a3(double(9./6.));
    scalar a4(double(2./6.));

    if (mesh_.time().timeIndex() < 3)
    {
        a1 = 1.0;
        a2 = 1.0;
        a3 = 0.0;
        a4 = 0.0;
    }

    tmp<fvScalarMatrix> tfvm
    (
        new fvScalarMatrix
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvScalarMatrix& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/deltaT_();

    fvm.diag() = (a1*rDeltaT)*mesh_.V();
    fvm.source() = rDeltaT*mesh_.V()*
    (
         a2*vf.oldTime().primitiveField()
       - a3*vf.oldTime().oldTime().primitiveField()
       + a4*vf.oldTime().oldTime().oldTime().primitiveField()
    );

    return tfvm;
}

surfaceScalarField Foam::BDF3::phiOldAndRelax
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& atf,
    const surfaceScalarField& atildf,
    const scalar& alphaU,
    const bool& consistentRhieChow
)
{
    scalar a1(double(11./6.));
    scalar a2(double(18./6.));
    scalar a3(double(9./6.));
    scalar a4(double(2./6.));

    if (mesh_.time().timeIndex() < 3)
    {
        a1 = 1.0;
        a2 = 1.0;
        a3 = 0.0;
        a4 = 0.0;
    }

    if (consistentRhieChow)
    {
        return
        (
              alphaU *atf*phi.oldTime()/atildf*(a2/a1)
            - alphaU *atf*phi.oldTime().oldTime()/atildf*(a3/a1)
            + alphaU *atf*phi.oldTime().oldTime().oldTime()/atildf*(a4/a1)
            + (1.0-alphaU)*phi.prevIter()
        );
    }
    else
    {
        return
        (
              alphaU *atf*fvc::flux(U.oldTime())/atildf*(a2/a1)
            - alphaU *atf*fvc::flux(U.oldTime().oldTime())/atildf*(a3/a1)
            + alphaU *atf*fvc::flux(U.oldTime().oldTime().oldTime())/atildf*(a4/a1)
            + (1.0-alphaU)*phi.prevIter()
        );
    }
}

scalar BDF3::deltaT_() const
{
    return mesh_.time().deltaTValue();
}

scalar BDF3::deltaT0_() const
{
    return mesh_.time().deltaT0Value();
}

scalar BDF3::deltaT0_(const fvMesh& mesh) const
{
    if (mesh_.time().timeIndex() < 3)
    {
        return GREAT;
    }
    else
    {
        return deltaT0_();
    }
}

void Foam::BDF3::message(int i)
{
    if (mesh_.time().timeIndex() < 3)
    {
        Info<<"Starting with DIRK22"<<" : ITERATION "<<i+1<<"/"<<a_.size()<<endl;
    }
    else
    {
        Info<<"TIME SCHEME : BDF3"<<endl;
    }
}

bool Foam::BDF3::read()
{
    if (read())
    {
        return true;
    }

    return false;
}

// ************************************************************************* //

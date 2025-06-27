/* -------------------------------------------------------------------------- //

Class
    Foam::BDF3

Description
    Third order implicit BDF3 scheme

SourceFiles
    BDF3.C

\*---------------------------------------------------------------------------*/

#ifndef BDF3_H
#define BDF3_H

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class BDF3 Declaration
\*---------------------------------------------------------------------------*/

class BDF3
:
    public timeSchemes
{
    // Private Member Functions

        List<scalarList> a_;

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        scalar deltaT0_(const fvMesh& mesh) const;

public:

    // Constructors
    TypeName("BDF3");

        //- Construct from components
        BDF3
        (
            const fvMesh& mesh
        );

    // Member Functions
        static autoPtr<BDF3> New
        (
            const fvMesh& mesh
        );

        //- Give message for the consol
        virtual void message(int i) ;

        //- Return Butcher table
        virtual List<scalarList> a();

        //- Return true if the RK scheme has explicit step
        virtual bool explicitFirstStep()
        {
            return false;
        }
        //- Return the sum aijRj
        virtual volVectorField sumAijRj
        (
            int i,
            const List<autoPtr<volVectorField>>& RiPtr
        );

        virtual volScalarField sumAijRj
        (
            int i,
            const List<autoPtr<volScalarField>>& RiPtr
        );

        virtual surfaceScalarField sumAijRfj
        (
            int i,
            const List<autoPtr<surfaceScalarField>>& RfiPtr
        );

        //- Return the derivative of psi
        virtual fvVectorMatrix ddt
        (
            const volVectorField& vf
        );

        virtual fvScalarMatrix ddt
        (
            const volScalarField& vf
        );


        //- Return the old contribution of fluxes
        // for Rie Chow interpolation

        virtual surfaceScalarField phiOldAndRelax
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const surfaceScalarField& atf,
            const surfaceScalarField& atildf,
            const scalar& alphaU,
            const bool& consistentRhieChow
        );
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

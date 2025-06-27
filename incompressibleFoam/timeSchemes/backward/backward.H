/* -------------------------------------------------------------------------- //

Class
    Foam::backward

Description
    2d order implicit Euler scheme

SourceFiles
    backward.C

\*---------------------------------------------------------------------------*/

#ifndef backward_H
#define backward_H

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class backward Declaration
\*---------------------------------------------------------------------------*/

class backward
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
    TypeName("backward");

        //- Construct from components
        backward
        (
            const fvMesh& mesh
        );

    // Member Functions
        static autoPtr<backward> New
        (
            const fvMesh& mesh
        );

        //- Give message for the consol
        virtual void message(int i);

        //- Return Butcher table
        virtual List<scalarList> a();

        //- Return true if the RK scheme has explicit step
        virtual bool explicitFirstStep() { return false; }

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

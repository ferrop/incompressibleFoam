/* -------------------------------------------------------------------------- //

Class
    Foam::localEuler

Description
    First order implicit localEuler scheme

SourceFiles
    localEuler.C

\*---------------------------------------------------------------------------*/

#ifndef localEuler_H
#define localEuler_H

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class localEuler Declaration
\*---------------------------------------------------------------------------*/

class localEuler
:
    public timeSchemes
{
    // Private Member Functions

        List<scalarList> a_;

public:

    // Constructors
    TypeName("localEuler");

        //- Construct from components
        localEuler
        (
            const fvMesh& mesh
        );

    // Member Functions
        static autoPtr<localEuler> New
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

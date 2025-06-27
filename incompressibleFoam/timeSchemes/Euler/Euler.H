/* -------------------------------------------------------------------------- //

Class
    Foam::Euler

Description
    First order implicit Euler scheme

SourceFiles
    Euler.C

\*---------------------------------------------------------------------------*/

#ifndef Euler_H
#define Euler_H

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class Euler Declaration
\*---------------------------------------------------------------------------*/

class Euler
:
    public timeSchemes
{
    // Private Member Functions

        List<scalarList> a_;

public:

    // Constructors
    TypeName("Euler");

        //- Construct from components
        Euler
        (
            const fvMesh& mesh
        );

    // Member Functions
        static autoPtr<Euler> New
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
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/* -------------------------------------------------------------------------- //

Class
    Foam::DIRK33

Description
    Third order DIRK33 scheme

SourceFiles
    DIRK33.C

\*---------------------------------------------------------------------------*/

#ifndef DIRK33_H
#define DIRK33_H

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class DIRK33 Declaration
\*---------------------------------------------------------------------------*/

class DIRK33
:
    public timeSchemes
{
    // Private Member Functions

        List<scalarList> a_;

public:

    // Constructors
    TypeName("DIRK33");

        //- Construct from components
        DIRK33
        (
            const fvMesh& mesh
        );

    // Member Functions
        static autoPtr<DIRK33> New
        (
            const fvMesh& mesh
        );

        //- Give message for the consol
        virtual void message(int i);

        //- Return Butcher table
        virtual List<scalarList> a();

        //- Return true if the RK scheme has explicit step
        virtual bool explicitFirstStep() { return false; }

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

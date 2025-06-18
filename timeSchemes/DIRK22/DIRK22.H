/* -------------------------------------------------------------------------- //

Class
    Foam::DIRK22

Description
    2nd order DIRK22 scheme

SourceFiles
    DIRK22.C

\*---------------------------------------------------------------------------*/

#ifndef DIRK22_H
#define DIRK22_H

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class DIRK22 Declaration
\*---------------------------------------------------------------------------*/

class DIRK22
:
    public timeSchemes
{
    // Private Member Functions

        List<scalarList> a_;

public:

    // Constructors
    TypeName("DIRK22");

        //- Construct from components
        DIRK22
        (
            const fvMesh& mesh
        );

    // Member Functions
        static autoPtr<DIRK22> New
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

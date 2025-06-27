/* -------------------------------------------------------------------------- //

Class
    Foam::DIRK43

Description
    Third order DIRK43 scheme

SourceFiles
    DIRK43.C

\*---------------------------------------------------------------------------*/

#ifndef DIRK43_H
#define DIRK43_H

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class DIRK43 Declaration
\*---------------------------------------------------------------------------*/

class DIRK43
:
    public timeSchemes
{
    // Private Member Functions

        List<scalarList> a_;

public:

    // Constructors
    TypeName("DIRK43");

        //- Construct from components
        DIRK43
        (
            const fvMesh& mesh
        );

    // Member Functions
        static autoPtr<DIRK43> New
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

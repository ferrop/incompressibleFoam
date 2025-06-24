/* -------------------------------------------------------------------------- //

Class
    Foam::EDIRK23_2

Description
    third order two stages EDIRK23_2 scheme

SourceFiles
    EDIRK23_2.C

\*---------------------------------------------------------------------------*/

#ifndef EDIRK23_2_H
#define EDIRK23_2_H

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class EDIRK23_2 Declaration
\*---------------------------------------------------------------------------*/

class EDIRK23_2
:
    public timeSchemes
{
    // Private Member Functions

        List<scalarList> a_;

public:

    // Constructors
    TypeName("EDIRK23_2");

        //- Construct from components
        EDIRK23_2
        (
            const fvMesh& mesh
        );

    // Member Functions
        static autoPtr<EDIRK23_2> New
        (
            const fvMesh& mesh
        );

        //- Give message for the consol
        virtual void message(int i);

        //- Return Butcher table
        virtual List<scalarList> a();

        //- Return true if the RK scheme has explicit step
        virtual bool explicitFirstStep() { return true; }

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

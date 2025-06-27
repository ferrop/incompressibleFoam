/* -------------------------------------------------------------------------- //

Class
    Foam::EDIRK23_1

Description
    third order two stages EDIRK23_1 scheme

SourceFiles
    EDIRK23_1.C

\*---------------------------------------------------------------------------*/

#ifndef EDIRK23_1_H
#define EDIRK23_1_H

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class EDIRK23_1 Declaration
\*---------------------------------------------------------------------------*/

class EDIRK23_1
:
    public timeSchemes
{
    // Private Member Functions

        List<scalarList> a_;

public:

    // Constructors
    TypeName("EDIRK23_1");

        //- Construct from components
        EDIRK23_1
        (
            const fvMesh& mesh
        );

    // Member Functions
        static autoPtr<EDIRK23_1> New
        (
            const fvMesh& mesh
        );

        //- Give message for the consol
        virtual void message(int i);

        //- Return Butcher table
        virtual List<scalarList> a();

        //- Return true if the RK scheme has explicit step
        virtual bool explicitFirstStep()
        {
            if (mesh_.time().timeIndex() < 2)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

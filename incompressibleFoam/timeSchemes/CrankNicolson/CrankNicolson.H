/* -------------------------------------------------------------------------- //

Class
    Foam::CrankNicolson

Description
    2nd order impplicit CrankNicolson scheme

SourceFiles
    CrankNicolson.C

\*---------------------------------------------------------------------------*/

#ifndef CrankNicolson_H
#define CrankNicolson_H

#include "timeSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class CrankNicolson Declaration
\*---------------------------------------------------------------------------*/

class CrankNicolson
:
    public timeSchemes
{
    // Private Member Functions

        List<scalarList> a_;

public:

    // Constructors
    TypeName("CrankNicolson");

        //- Construct from components
        CrankNicolson
        (
            const fvMesh& mesh
        );

    // Member Functions
        static autoPtr<CrankNicolson> New
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
            return true;
        }
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\

Class
    Foam::timeSchemes

Description
    Class for time schemes

SourceFiles
    timeSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef timeSchemes_H
#define timeSchemes_H

#include "runTimeSelectionTables.H"
#include "fvCFD.H"
#include "IOMRFZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class timeSchemes Declaration
\*---------------------------------------------------------------------------*/

class timeSchemes
:
    public IOdictionary
{
    // Private data
protected:

        const fvMesh& mesh_;

    // Private Member Functions

        //- No copy construct
        timeSchemes
        (
            const timeSchemes&
        ) = delete;

        //- No copy assignment
        void operator=(const timeSchemes&) = delete;


public:

    // Constructors
    TypeName("timeSchemes");

        //- Construct from components
        timeSchemes
        (
            const fvMesh& mesh
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            timeSchemes,
            components,
            (const fvMesh& mesh),
            (mesh)
        );

    //- Destructor
    virtual ~timeSchemes();

    // Member Functions
        static autoPtr<timeSchemes> New
        (
            const fvMesh& mesh
        );

        //- Correction of flux for consistence with
        // U boundary conditions
        void constrainFlux
        (
            surfaceScalarField& phi,
            const volVectorField& U,
            const volScalarField& p,
            const IOMRFZoneList& MRF
        );

        //- 2nd order extrapolation
        void extrapolateFields
        (
            int i,
            surfaceScalarField& phi,
            volVectorField& U,
            volScalarField& p
        );

        //- Return Butcher table
        virtual List<scalarList> a() = 0;

        //- Return c(i) abscissa coefficient
        scalar c
        (
            int i
        );

        //- Give message for the consol
        virtual void message(int i) = 0;

        //- Return true if the RK scheme has explicit step
        virtual bool explicitFirstStep() = 0;

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

        //- Return the sum aijRfj
        virtual surfaceScalarField sumAijRfj
        (
            int i,
            const List<autoPtr<surfaceScalarField>>& RfiPtr,
            bool orient
        );

        //- Return the sum aijRfj
        virtual surfaceVectorField sumAijRfj
        (
            int i,
            const List<autoPtr<surfaceVectorField>>& RfiPtr,
            bool orient
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

        //- Velocity Update
        virtual volVectorField updateU
        (
            const int& i,
            const volVectorField& U,
            const volVectorField& sumAijRUj,
            const volVectorField& convDiffH,
            const fvVectorMatrix& ddtUEqn,
            const volScalarField& atild,
            const volScalarField& p,
            const scalar& alphaU
        );
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

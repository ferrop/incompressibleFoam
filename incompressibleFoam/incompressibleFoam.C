/* -------------------------------------------------------------------------- //

Application
    incompressibleFoam.C

Description
    Steady and Transient solver for incompressible,
    turbulent flow of Newtonian fluids
    with BDF and DIRK integration methods.

    Two momentum interpolations (Rhie & Chow interpolation) are available
    as well as two pressure formulations (standard and corrected)

Authors :
Paulin FERRO, Pierre-Etienne MEILLER, Paul LANDEL, Carla  LANDRODIE and Marc PESCHEUX
SARL G-MET Technologies, www.g-met.fr,
63 rue d'Hyères, Six-Fours-Les-Plages, FRANCE

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "fvcSmooth.H"
#include "timeSchemes.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"

    turbulence->validate();

    if (!LTS && (timeScheme != "steadyState"))
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else if (timeScheme != "steadyState")
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;
        scalar dT(runTime.deltaTValue());

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "createR.H"

        for (int i = 0; i < nRKloop ; i++)
        {
            timeMethod->message(i);

            #include "initializeVariable.H"

            //Let's start Pressure-Velocity coupling !
            if ( (i == 0 && timeMethod->explicitFirstStep()) )
            {
                //For ESDIRK (CrankNicolson, ESDIRK23...)
                //we need a empty iteration for calculating R and Rf
                #include "calcR0.H"
            }
            else
            {
                while (pimple.loop())
                {
                    p.storePrevIter();
                    phi.storePrevIter();
                    U.storePrevIter();
                    MRF.makeRelative(phi);

                    #include "relaxationCoefficient.H"
                    #include "UEqn.H"

                    while (pimple.correct())
                    {
                        #include "pEqn.H"
                    }

                    laminarTransport.correct();
                    {
                       // For (E)SDIRK, we need to modify
                       // the advection time to be consistent
                       // with the RK abscissae
                       runTime.setDeltaT(ci*dT);
                       turbulence->correct();
                       runTime.setDeltaT(dT);
                    }

                    MRF.makeAbsolute(phi);
                    // Calculate R and Rf after PU update
                    #include "updateR.H"
                }
            }
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

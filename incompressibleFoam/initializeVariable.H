//If PIMPLE is used (pimple.nCorrPIMPLE() > 1) :
//we decide to reset the variables at the begining of a Runge-Kutta loop.
//It is more consistent than using the available fields = previous RK loop
//especially in case of PISO mode.
//If PISO is used (pimple.nCorrPIMPLE() = 1) :
//2nd order flux, velocity and pression extrapolations is recommended
// fieldsExtrapolation = true
//This improves time accuracy

bool consistentRhieChow =
(
    pimple.dict().lookupOrDefault<bool>("consistentRhieChow", true)
);
bool fieldsExtrapolation =
(
    pimple.dict().lookupOrDefault<bool>("fieldsExtrapolation", false)
);
bool pressureCorrectionForm =
(
    pimple.dict().lookupOrDefault<bool>("pressureCorrectionForm", false)
);

scalar ci(timeMethod->c(i));

if (fieldsExtrapolation)
{
    timeMethod->extrapolateFields(i,phi,U,p);
}
else
{
    phi = phi.oldTime() ;
    U   = U.oldTime()   ;
    p   = p.oldTime()   ;
}
if(pressureCorrectionForm)
{
    pc == 0*pc;
}
// Centre et face source term calculation for DIRK
volVectorField      sumAijRUj(timeMethod->sumAijRj(i,RU));
surfaceScalarField sumAijRUfj(timeMethod->sumAijRfj(i,RUf,true));

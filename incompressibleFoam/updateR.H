if (pimple.finalIter())
{
    RU[i]()  = - convDiffUEqn.A()*U + convDiffUEqn.H() ;

    if (consistentRhieChow)
    {
        RUf[i]() = - linearInterpolate(convDiffUEqn.A())*phi + fvc::flux(convDiffUEqn.H()) ;
    }
    else
    {
        RUf[i]() = fvc::flux(RU[i]());
    }
}

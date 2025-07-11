if (!pimple.flow())
{
    if (pimple.models())
    {
        fvModels.correct();
    }

    if (pimple.thermophysics())
    {
        #include "YEqn.H"
        #include "EEqn.H"
    }
}
else
{
    if (!mesh.steady() && pimples.firstPimpleIter())
    {
        #include "rhoEqn.H"
    }

    if (pimple.models())
    {
        fvModels.correct();
    }

    #include "UEqn.H"

    if (pimple.thermophysics())
    {
        #include "YEqn.H"
        #include "EEqn.H"
    }

    // --- PISO loop
    while (pimple.correct())
    {
        #include "pEqn.H"
    }

    if (pimples.pimpleTurbCorr(i))
    {
        turbulence.correct();
        thermophysicalTransport.correct();
    }

    if (!mesh.steady() && pimples.finalPimpleIter())
    {
        rho = thermo.rho();
    }

    enthalpy = thermo.he();
}

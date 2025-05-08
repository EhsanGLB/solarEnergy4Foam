
if(typeFluid[i] == "forcedConvection")
{
    #include "forcedConvection/UEqnFC.H"

    for (int corr = 0; corr < nCorr; corr++)
    {
        #include "forcedConvection/pEqnFC.H"
    }
}


if(typeFluid[i] == "naturalConvection")
{
    #include "naturalConvection/UEqnNC.H"

    for (int corr = 0; corr < nCorr; corr++)
    {
        #include "naturalConvection/pEqnNC.H"
    }
}


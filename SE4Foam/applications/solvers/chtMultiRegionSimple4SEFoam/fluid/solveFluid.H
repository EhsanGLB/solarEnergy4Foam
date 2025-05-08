//  Pressure-velocity SIMPLE corrector

    p.storePrevIter();


if(typeFluid[i] == "forcedConvection")
{
    #include "forcedConvection/UEqnFC.H"
    #include "forcedConvection/pEqnFC.H"
    #include "forcedConvection/TEqnFC.H"
}


if(typeFluid[i] == "naturalConvection")
{
    #include "naturalConvection/UEqnNC.H"
    #include "naturalConvection/pEqnNC.H"
    #include "naturalConvection/TEqnNC.H"
}

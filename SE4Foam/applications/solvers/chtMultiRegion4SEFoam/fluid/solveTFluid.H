
if(typeFluid[i] == "forcedConvection")
{
    #include "forcedConvection/TEqnFC.H"
}


if(typeFluid[i] == "naturalConvection")
{
    #include "naturalConvection/TEqnNC.H"
}



#include <iostream>
#include "CoolPropLib.h"
int main(){


    double Tc = 300.15;
    double Pc = 1e5;
    double rhoc = PropsSI("D", "T", Tc, "P", Pc, "INCOMP::TVP1");
    double kappac = PropsSI("L", "T", Tc, "P", Pc, "INCOMP::TVP1");
    double muc = PropsSI("V", "T", Tc, "P", Pc, "INCOMP::TVP1");
    double Cpc = PropsSI("C", "T", Tc, "P", Pc, "INCOMP::TVP1");

    std::cout << "rhoc: " << rhoc << std::endl;
    std::cout << "kappac: " << kappac << std::endl;
    std::cout << "muc: " << muc << std::endl;
    std::cout << "Cpc: " << Cpc << std::endl;
}


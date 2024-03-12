
    rho = fluidPtr->rho();
    kappa = fluidPtr->kappa();
    Cp = fluidPtr->Cp();
    mu = fluidPtr->mu();

    turbulence->correct();

    kappaEff = porosity*((turbulence->nut()/Prt)*rho*Cp + kappa) + (1-porosity)*kappaPM;


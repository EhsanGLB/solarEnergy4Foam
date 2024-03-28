
    rho = fluidPtr->rho();
    kappa = fluidPtr->kappa();
    Cp = fluidPtr->Cp();
    mu = fluidPtr->mu();

    turbulence->correct();

    kappat = (turbulence->nut()/Prt)*rho*Cp;
    kappat.correctBoundaryConditions();
    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;


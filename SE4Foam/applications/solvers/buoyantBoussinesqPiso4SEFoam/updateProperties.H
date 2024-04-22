
    rho = fluidPtr->rho();
    kappa = fluidPtr->kappa();
    Cp = fluidPtr->Cp();
    mu = fluidPtr->mu();

    turbulence->correct();
    kappat = (turbulence->nut()/Prt)*rho*Cp;
    kappat.correctBoundaryConditions();

    rhoEff = porosity*rho + (1-porosity)*rhoPM;
    CpEff = (porosity*rho*Cp + (1-porosity)*rhoPM*CpPM)/rhoEff;
    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;


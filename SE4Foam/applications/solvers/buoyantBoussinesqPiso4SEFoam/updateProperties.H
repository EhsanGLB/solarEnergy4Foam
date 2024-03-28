
    rho = fluidPtr->rho();
    kappa = fluidPtr->kappa();
    Cp = fluidPtr->Cp();
    mu = fluidPtr->mu();

    turbulence->correct();

    rhoEff = porosity*rho + (1-porosity)*rhoPM;
    CpEff = (porosity*rho*Cp + (1-porosity)*rhoPM*CpPM)/rhoEff;
    kappat = (turbulence->nut()/Prt)*rho*Cp;
    kappat.correctBoundaryConditions();
    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;


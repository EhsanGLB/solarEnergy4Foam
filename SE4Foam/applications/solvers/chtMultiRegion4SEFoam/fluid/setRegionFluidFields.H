    const fvMesh& mesh = fluidRegions[i];

    volScalarField& T = TFluid[i];
    volScalarField& p = pFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];


    volScalarField& rho = rhoFluid[i];
    volScalarField& Cp = CpFluid[i];
    volScalarField& kappa = kappaFluid[i];
    volScalarField& mu = muFluid[i];
    volScalarField& beta = betaFluid[i];
    word& type = typeFluid[i];
    fluid& fluidProp = fluidPropertiesFluid[i];

    rho = fluidProp.rho();
    kappa = fluidProp.kappa();
    Cp = fluidProp.Cp();
    mu = fluidProp.mu();
    rho.correctBoundaryConditions();
    Cp.correctBoundaryConditions();
    kappa.correctBoundaryConditions();
    mu.correctBoundaryConditions();

    if(type == "naturalConvection")
    {
        dimensionedScalar& betaf = betafFluid[i];
        beta = betaf;
        beta.correctBoundaryConditions();
    }

    mesh.schemesDict().setFluxRequired(p.name());


    incompressible::turbulenceModel& turb = turbulence[i];
    dimensionedScalar& Prt = PrtFluid[i];
    turb.correct();
    volScalarField& kappat = kappatFluid[i];
    kappat = (turb.nut()/Prt)*rho*Cp;
    kappat.correctBoundaryConditions();


    porousZones& pZones = pZonesFluid[i];
    volScalarField& porosity = porosityFluid[i];
    volScalarField& rhoPM = rhoPMFluid[i];
    volScalarField& CpPM = CpPMFluid[i];
    volScalarField& kappaPM = kappaPMFluid[i];
    volScalarField& hsPM = hsPMFluid[i];


    volScalarField& rhoEff = rhoEffFluid[i];
    volScalarField& CpEff = CpEffFluid[i];
    volScalarField& kappaEff = kappaEffFluid[i];
    rhoEff = porosity*rho + (1-porosity)*rhoPM;
    CpEff = (porosity*rho*Cp + (1-porosity)*rhoPM*CpPM)/rhoEff;
    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;
    rhoEff.correctBoundaryConditions();
    CpEff.correctBoundaryConditions();
    kappaEff.correctBoundaryConditions();


    radiation::radiationModel& rad = radiation[i];

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
    dimensionedScalar& Prt = PrtFluid[i];
    incompressible::turbulenceModel& turb = turbulence[i];
    fluid& fluidProp = fluidPropertiesFluid[i];

    rho = fluidProp.rho();
    kappa = fluidProp.kappa();
    Cp = fluidProp.Cp();
    mu = fluidProp.mu();

    if(type == "naturalConvection")
    {
        dimensionedScalar& betaf = betafFluid[i];
        beta = betaf;
    }

    mesh.schemesDict().setFluxRequired(p.name());

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];
    turb.correct();

    radiation::radiationModel& rad = radiation[i];

    porousZones& pZones = pZonesFluid[i];
    volScalarField& porosity = porosityFluid[i];
    volScalarField& rhoPM = rhoPMFluid[i];
    volScalarField& CpPM = CpPMFluid[i];
    volScalarField& kappaPM = kappaPMFluid[i];

    volScalarField& kappaEff = kappaEffFluid[i];

    kappaEff = porosity*((turb.nut()/Prt) * rho * Cp + kappa) + (1-porosity)*kappaPM;

    // Initialise fluid field pointer lists
    PtrList<volScalarField> TFluid(fluidRegions.size());
    PtrList<volScalarField> pFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    

    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to TFluid\n" << endl;
        TFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to pFluid\n" << endl;
        pFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(UFluid[i]) & fluidRegions[i].Sf()
            )
        );
    }

//------------------- SE4Foam ------------------//
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volScalarField> CpFluid(fluidRegions.size());
    PtrList<volScalarField> kappaFluid(fluidRegions.size());
    PtrList<volScalarField> muFluid(fluidRegions.size());
    PtrList<volScalarField> betaFluid(fluidRegions.size());
    PtrList<singlePhaseTransportModel> laminarTransport(fluidRegions.size());
    PtrList<word> typeFluid(fluidRegions.size());
    PtrList<fluid> fluidPropertiesFluid(fluidRegions.size());


    forAll(fluidRegions, i)
    {
        Info<< "    Adding to rhoFluid\n" << endl;
        rhoFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("rho", dimDensity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to CpFluid\n" << endl;
        CpFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Cp",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("Cp", dimSpecificHeatCapacity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to kappaFluid\n" << endl;
        kappaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappa",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("kappa", dimThermalConductivity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to muFluid\n" << endl;
        muFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "mu",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("mu", dimMass/dimLength/dimTime, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to betaFluid\n" << endl;
        betaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "beta",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("beta", dimless/dimTemperature, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to laminarTransport\n" << endl;
        laminarTransport.set
        (
            i,
            new singlePhaseTransportModel
            (
                UFluid[i],
                phiFluid[i]
            )
        );

        typeFluid.set
        (
            i,
            new word (laminarTransport[i].lookup("type"))
        );

        Info<< "    Adding to fluidProperties\n" << endl;
        fluidPropertiesFluid.set
        (
            i,
            new fluid
            (
                laminarTransport[i],
                UFluid[i],
                pFluid[i],
                TFluid[i]
            )
        );
    }

//------------------- SE4Foam -> Base fluid ------------------//
    // forced convection
    PtrList<volScalarField> rAUFluid(fluidRegions.size());
    PtrList<pimpleControl> pimpleFluid(fluidRegions.size());


    // natural convection
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<dimensionedScalar> betafFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TRefFluid(fluidRegions.size());
    PtrList<volScalarField> rhokFluid(fluidRegions.size());


    forAll(fluidRegions, i)
    {
      // forced convection
      if(typeFluid[i] == "forcedConvection")
      {
        Info<< "Reading field rAU if present\n" << endl;
        rAUFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rAU",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                runTime.deltaT(),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        Info<< "    Adding to pimpleFluid\n" << endl;
        pimpleFluid.set
        (
            i,
            new pimpleControl
            (
                fluidRegions[i],
                "PIMPLE"
            )
        );
      }


      // natural convection
      if(typeFluid[i] == "naturalConvection")
      {
        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        if(laminarTransport[i].found("beta"))
        {
            betafFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("beta"))
            );
        }

        if(laminarTransport[i].found("TRef"))
        {
            TRefFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("TRef"))
            );
        }

        Info<< "    Adding to rhokFluid\n" << endl;
        rhokFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhok",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("rhok", dimless, SMALL),
                "zeroGradient"
            )
        );
      }
    }

//------------------- SE4Foam -> Turbulence ------------------//
    PtrList<volScalarField> kappatFluid(fluidRegions.size());
    PtrList<incompressible::turbulenceModel> turbulence(fluidRegions.size());
    PtrList<dimensionedScalar> PrtFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        Info<< "    Adding to kappat\n" << endl;
        kappatFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappat",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("kappat", dimThermalConductivity, SMALL)
            )
        );

        Info<< "    Adding to turbulence\n" << endl;
        turbulence.set
        (
            i,
            autoPtr<incompressible::turbulenceModel>
            (
                incompressible::turbulenceModel::New
                (
                    UFluid[i],
                    phiFluid[i],
                    laminarTransport[i]
                )
            ).ptr()
        );

        PrtFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("Prt"))
        );
    }


//------------------- SE4Foam -> Porous Media ------------------//
    PtrList<porousZones> pZonesFluid(fluidRegions.size());
    PtrList<volScalarField> porosityFluid(fluidRegions.size());
    PtrList<volScalarField> rhoPMFluid(fluidRegions.size());
    PtrList<volScalarField> CpPMFluid(fluidRegions.size());
    PtrList<volScalarField> kappaPMFluid(fluidRegions.size());
    PtrList<volScalarField> hsPMFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        Info<< "    Adding to pZones\n" << endl;
        pZonesFluid.set
        (
            i,
            new porousZones
            (
                fluidRegions[i]
            )
        );

        Info<< "    Adding to porosity\n" << endl;
        porosityFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "porosity",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("porosity", dimless, 1.0),
                "zeroGradient"
            )
        );

        Info<< "    Adding to rhoPM\n" << endl;
        rhoPMFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhoPM",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("rhoPM", dimDensity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to CpPM\n" << endl;
        CpPMFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "CpPM",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("CpPM", dimSpecificHeatCapacity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to kappaPM\n" << endl;
        kappaPMFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappaPM",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("kappaPM", dimThermalConductivity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to hsPM\n" << endl;
        hsPMFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "hsPM",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("hsPM", dimensionSet(1, -1, -3, 0, 0, 0, 0), SMALL),
                "zeroGradient"
            )
        );

        forAll( pZonesFluid[i], zoneI )
        {
            const label & zoneId( pZonesFluid[i][zoneI].zoneId() );
            const labelList & cells(fluidRegions[i].cellZones()[zoneId]);
            const scalar & zonePorosity( pZonesFluid[i][zoneI].porosity() );

            const dictionary& thermalPropertiesPM = pZonesFluid[i][zoneI].dict().subDict("thermalProperties");
            dimensionedScalar zoneRho(thermalPropertiesPM.lookup("rho"));
            dimensionedScalar zoneCp(thermalPropertiesPM.lookup("Cp"));
            dimensionedScalar zoneKappa(thermalPropertiesPM.lookup("kappa"));
            dimensionedScalar zoneHs(thermalPropertiesPM.lookup("hs"));

            forAll( cells, cellI )
            {
                porosityFluid[i][cells[cellI]] = zonePorosity;
                rhoPMFluid[i][cells[cellI]] = zoneRho.value();
                CpPMFluid[i][cells[cellI]] = zoneCp.value();
                kappaPMFluid[i][cells[cellI]] = zoneKappa.value();
                hsPMFluid[cells[cellI]] = zoneHs.value();
            }
        }

        rhoPMFluid[i].correctBoundaryConditions();
        CpPMFluid[i].correctBoundaryConditions();
        kappaPMFluid[i].correctBoundaryConditions();
        hsPMFluid[i].correctBoundaryConditions();
    }

//------------------- SE4Foam -> Effective thermal properties ------------------//
    PtrList<volScalarField> rhoEffFluid(fluidRegions.size());
    PtrList<volScalarField> CpEffFluid(fluidRegions.size());
    PtrList<volScalarField> kappaEffFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        Info<< "    Adding to rhoEff\n" << endl;
        rhoEffFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhoEff",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("rhoEff", dimDensity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to CpEff\n" << endl;
        CpEffFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "CpEff",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("CpEff", dimSpecificHeatCapacity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to kappaEff\n" << endl;
        kappaEffFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappaEff",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("kappaEff", dimThermalConductivity, SMALL),
                "zeroGradient"
            )
        );

    }

//------------------- SE4Foam -> Radiation ------------------//
    PtrList<radiation::radiationModel> radiation(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        Info<< "    Adding to radiation\n" << endl;
        radiation.set
        (
            i,
            radiation::radiationModel::New
            (
                TFluid[i]
            )
        );
    }

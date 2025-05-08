    Info<< "Reading thermophysical properties\n" << endl;

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//------------------- SE4Foam ------------------//
    Info<< "Constructing field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field Cp\n" << endl;
    volScalarField Cp
    (
        IOobject
        (
            "Cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Cp", dimSpecificHeatCapacity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field kappa\n" << endl;
    volScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field mu\n" << endl;
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mu", dimMass/dimLength/dimTime, SMALL),
        "zeroGradient"
    );

#   include "createPhi.H"

    Info<< "Reading transport Properties\n" << endl;
    singlePhaseTransportModel laminarTransport(U, phi);
    dimensionedScalar beta(laminarTransport.lookup("beta"));
    dimensionedScalar TRef(laminarTransport.lookup("TRef"));

    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        1.0 - beta*(T - TRef)
    );

    autoPtr<fluid> fluidPtr(fluid::New(laminarTransport, U, p, T));

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p.name());

//------------------- SE4Foam -> Turbulence ------------------//
    Info<< "Constructing field kappat\n" << endl;
    volScalarField kappat
    (
        IOobject
        (
            "kappat",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappat", dimThermalConductivity, SMALL)
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    dimensionedScalar Prt(laminarTransport.lookup("Prt"));

//------------------- SE4Foam -> Porous Media ------------------//
    porousZones pZones(mesh);

    IOdictionary porousZonesDict
    (
        IOobject
        (
            "porousZones",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    volScalarField porosity
    (
        IOobject
        (
            "porosity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("porosity", dimless, 1.0),
        "zeroGradient"
    );

    volScalarField rhoPM
    (
        IOobject
        (
            "rhoPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoPM", dimDensity, SMALL),
        "zeroGradient"
    );

    volScalarField CpPM
    (
        IOobject
        (
            "CpPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("CpPM", dimSpecificHeatCapacity, SMALL),
        "zeroGradient"
    );

    volScalarField kappaPM
    (
        IOobject
        (
            "kappaPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappaPM", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

    volScalarField hsPM
    (
        IOobject
        (
            "hsPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("hsPM", dimensionSet(1, -1, -3, 0, 0, 0, 0), SMALL),
        "zeroGradient"
    );

    forAll( pZones, zoneI )
    {
        const label & zoneId( pZones[zoneI].zoneId() );
        const labelList & cells(mesh.cellZones()[zoneId]);
        const scalar & zonePorosity( pZones[zoneI].porosity() );

        const dictionary& thermalPropertiesPM = pZones[zoneI].dict().subDict("thermalProperties");
        dimensionedScalar zoneRho(thermalPropertiesPM.lookup("rho"));
        dimensionedScalar zoneCp(thermalPropertiesPM.lookup("Cp"));
        dimensionedScalar zoneKappa(thermalPropertiesPM.lookup("kappa"));
        dimensionedScalar zoneHs(thermalPropertiesPM.lookup("hs"));

        forAll( cells, cellI )
        {
            porosity[cells[cellI]] = zonePorosity;
            rhoPM[cells[cellI]] = zoneRho.value();
            CpPM[cells[cellI]] = zoneCp.value();
            kappaPM[cells[cellI]] = zoneKappa.value();
            hsPM[cells[cellI]] = zoneHs.value();
        }
    }

    rhoPM.correctBoundaryConditions();
    CpPM.correctBoundaryConditions();
    kappaPM.correctBoundaryConditions();
    hsPM.correctBoundaryConditions();

//------------------- SE4Foam -> Effective thermal properties ------------------//
    Info<< "Constructing field kappaEff\n" << endl;
    volScalarField kappaEff
    (
        IOobject
        (
            "kappaEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappaEff", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

//------------------- SE4Foam -> Radiation ------------------//
    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(T)
    );

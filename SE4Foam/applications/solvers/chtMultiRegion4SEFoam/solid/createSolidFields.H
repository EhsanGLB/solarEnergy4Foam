    // Initialise solid field pointer lists
    PtrList<volScalarField> rhos(solidRegions.size());
    PtrList<volScalarField> cps(solidRegions.size());
    PtrList<volScalarField> rhosCps(solidRegions.size());
    PtrList<volScalarField> Kappas(solidRegions.size());
    PtrList<volScalarField> Ts(solidRegions.size());

    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding to rhos\n" << endl;
        rhos.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        Info<< "    Adding to cps\n" << endl;
        cps.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "cp",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        rhosCps.set
        (
            i,
            new volScalarField("rhosCps", rhos[i]*cps[i])
        );

        Info<< "    Adding to Kappas\n" << endl;
        Kappas.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Kappa",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        Info<< "    Adding to Ts\n" << endl;
        Ts.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );
    }

//------------------- SE4Foam -> Heat source ------------------//
    PtrList<IOdictionary> heatSource(solidRegions.size());
    PtrList<dimensionedScalar> hsSolid(solidRegions.size());

    forAll(solidRegions, i)
    {
        heatSource.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "heatSource",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        if(heatSource[i].found("hs"))
        {
            hsSolid.set
            (
                i,
                new dimensionedScalar (heatSource[i].lookup("hs"))
            );
        }
        else
        {
            hsSolid.set
            (
                i,
                new dimensionedScalar ("hs", dimThermalConductivity*dimTemperature/dimArea, SMALL)
            );
        }
    }

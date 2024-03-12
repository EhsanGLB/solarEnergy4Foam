
    forAll(T.internalField(), i)
    {
        scalar Ti = T.internalField()[i];
        rho.internalField()[i] = PropsSI("D", "T", Ti, "P", PRef, fluidName);
        Cp.internalField()[i] = PropsSI("C", "T", Ti, "P", PRef, fluidName);
        kappa.internalField()[i] = PropsSI("L", "T", Ti, "P", PRef, fluidName);
        mu.internalField()[i] = PropsSI("V", "T", Ti, "P", PRef, fluidName);
    }

    forAll(T.boundaryField(), i)
    {
	forAll(T.boundaryField()[i], j)
	{
            scalar Ti = T.boundaryField()[i][j];
            rho.boundaryField()[i][j] = PropsSI("D", "T", Ti, "P", PRef, fluidName);
            Cp.boundaryField()[i][j] = PropsSI("C", "T", Ti, "P", PRef, fluidName);
            kappa.boundaryField()[i][j] = PropsSI("L", "T", Ti, "P", PRef, fluidName);
            mu.boundaryField()[i][j] = PropsSI("V", "T", Ti, "P", PRef, fluidName);
	}
    }

    turbulence->correct();
    kappaEff = porosity*((turbulence->nut()/Prt)*rho*Cp + kappa) + (1-porosity)*kappaPM;

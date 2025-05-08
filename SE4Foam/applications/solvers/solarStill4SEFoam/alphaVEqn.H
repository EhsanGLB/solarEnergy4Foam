{
    fvScalarMatrix alphaVEqn
    (
        fvm::div(phi, alphaV)
      + fvm::SuSp(-fvc::div(phi), alphaV)
      - fvm::laplacian(DV, alphaV)
    );

    alphaVEqn.relax();
    alphaVEqn.solve();

    rhok = 1.0 - beta*(T - TRef) - betaV*(alphaV - alphaVRef);
}

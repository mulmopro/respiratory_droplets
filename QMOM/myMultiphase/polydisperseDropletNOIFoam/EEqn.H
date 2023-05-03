{
    fluid.correctEnergyTransport();

    autoPtr<twoPhaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    twoPhaseSystem::heatTransferTable&
        heatTransfer = heatTransferPtr();
/*
    Info<< "fluid.heatTransfer1()" << average(fluid.heatTransfer1()->D()) << endl;
    Info<< "fluid.heatTransfer2()" << average(fluid.heatTransfer2()->D()) << endl;

    Info<< "*heatTransfer[phase1.name()]" << average(heatTransfer[phase1.name()]->D()) << endl;
    Info<< "*heatTransfer[phase2.name()]" << average(heatTransfer[phase2.name()]->D()) << endl;

    Info<< "fluid.heatTransfer1()" << average(fluid.heatTransfer1()->source()) << endl;
    Info<< "*heatTransfer[phase1.name()]" << average(heatTransfer[phase1.name()]->source()) << endl;
*/

    
    volScalarField& he1 = thermo1.he();
    volScalarField& he2 = thermo2.he();

    tmp<volScalarField> tK1(phase1.K());
    const volScalarField& K1(tK1());

    tmp<volScalarField> tK2(phase2.K());
    const volScalarField& K2(tK2());
/*
    Info<< "K1: " << average(K1) << endl;
    Info<< "phaseK1: " << average(phaseK1) << endl;
    Info<< "K2: " << average(K2) << endl;
    Info<< "phaseK2: " << average(phaseK2) << endl;
*/
    //volScalarField Cpv1("Cpv1", thermo1.Cpv());
    //volScalarField Cpv2("Cpv2", thermo2.Cpv());
    
    //volScalarField Kh(fluid.Kh());

    fvScalarMatrix E1Eqn
    (
        fvm::ddt(alpha1, rho1, he1) + fvm::div(alphaRhoPhi1, he1)
      - fvm::Sp(contErr1, he1)

      + fvc::ddt(alpha1, rho1, K1) + fvc::div(alphaRhoPhi1, K1)
      - contErr1*K1
      + (
            he1.name() == thermo1.phasePropertyName("e")
          ? fvc::div(fvc::absolute(alphaPhi1, alpha1, U1), p)
          + p*fvc::ddt(alpha1)
          : -alpha1*dpdt
        )

      - fvm::laplacian
        (
            fvc::interpolate(alpha1)
           *fvc::interpolate(thermo1.alphaEff(phase1.turbulence().mut())),
            he1
        )
    );

    E1Eqn.relax();

    E1Eqn -=
    (
        //fluid.heatTransfer1()
        *heatTransfer[phase1.name()]
      + alpha1*rho1*(U1&g)
      + fvOptions(alpha1, rho1, he1)
    );

    fvScalarMatrix E2Eqn
    (
        fvm::ddt(alpha2, rho2, he2) + fvm::div(alphaRhoPhi2, he2)
      - fvm::Sp(contErr2, he2)

      + fvc::ddt(alpha2, rho2, K2) + fvc::div(alphaRhoPhi2, K2)
      - contErr2*K2
      + (
            he2.name() == thermo2.phasePropertyName("e")
          ? fvc::div(fvc::absolute(alphaPhi2, alpha2, U2), p)
          + p*fvc::ddt(alpha1)
          : -alpha2*dpdt
        )

      - fvm::laplacian
        (
            fvc::interpolate(alpha2)
           *fvc::interpolate(thermo2.alphaEff(phase2.turbulence().mut())),
            he2
        )
    );

    E2Eqn.relax();

    E2Eqn -=
    (
        //fluid.heatTransfer2()
        *heatTransfer[phase2.name()]
      + alpha2*rho2*(U2&g)
      + fvOptions(alpha2, rho2, he2)
    );

    fvOptions.constrain(E1Eqn);
    E1Eqn.solve();
    fvOptions.correct(he1);

    fvOptions.constrain(E2Eqn);
    E2Eqn.solve();
    fvOptions.correct(he2);

    phase1.correctThermo();
    thermo1.correct();

    Info<< "min " << thermo1.T().name()
        << " " << min(thermo1.T()).value() << endl;
    
    phase2.correctThermo();
    thermo2.correct();
    
    Info<< "min " << thermo2.T().name()
        << " " << min(thermo2.T()).value() << endl;
}

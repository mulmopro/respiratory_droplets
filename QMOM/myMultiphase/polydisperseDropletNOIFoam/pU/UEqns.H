Info<< "Constructing momentum equations" << endl;

MRF.correctBoundaryVelocity(U1);
MRF.correctBoundaryVelocity(U2);
MRF.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

volScalarField Kd(fluid.Kd());
{
    volScalarField Vm(fluid.Vm());

    const volScalarField dmdt(fluid.dmdt());
    const volScalarField dmdt12(posPart(dmdt));
    const volScalarField dmdt21(negPart(dmdt));

    //Info<< "Average dmdt: " << average(dmdt) << endl;
    //    << "Average dmdt21: " << average(dmdt21) << endl
    //    << "Average dmdt12: " << average(dmdt12) << endl;

    {
        U1Eqn =
        (
            fvm::ddt(alpha1, rho1, U1) + fvm::div(alphaRhoPhi1, U1)
          - fvm::Sp(contErr1, U1)
          + MRF.DDt(alpha1*rho1 + Vm, U1)
          + fluid.divDevRhoReff1()
         ==
          - Vm
           *(
                fvm::ddt(U1)
              + fvm::div(phi1, U1)
              - fvm::Sp(fvc::div(phi1), U1)
              - DDtU2
            )
         + fvOptions(alpha1, rho1, U1)
        );
        //U1Eqn += dmdt21*U2 - fvm::Sp(dmdt21, U1);
        //U2Eqn += dmdt12*U2 - fvm::Sp(dmdt12, U1);
        U1Eqn -= dmdt21*U2 + fvm::Sp(dmdt12, U1);
        U1Eqn.relax();
        U1Eqn += fvm::Sp(Kd, U1);
        fvOptions.constrain(U1Eqn);
        U1.correctBoundaryConditions();
        fvOptions.correct(U1);
    }

    {
        U2Eqn =
        (
            fvm::ddt(alpha2, rho2, U2) + fvm::div(alphaRhoPhi2, U2)
          - fvm::Sp(contErr2, U2)
          + MRF.DDt(alpha2*rho2 + Vm, U2)
          + fluid.divDevRhoReff2()
         ==
          - Vm
           *(
                fvm::ddt(U2)
              + fvm::div(phi2, U2)
              - fvm::Sp(fvc::div(phi2), U2)
              - DDtU1
            )
          + fvOptions(alpha2, rho2, U2)
        );
        //U2Eqn += dmdt21*U1 - fvm::Sp(dmdt21, U2);
        //U2Eqn += dmdt12*U1 - fvm::Sp(dmdt12, U2);
        U2Eqn += dmdt12*U1 + fvm::Sp(dmdt21, U2);
        U2Eqn.relax();
        U2Eqn += fvm::Sp(Kd, U2);
        fvOptions.constrain(U2Eqn);
        U2.correctBoundaryConditions();
        fvOptions.correct(U2);
    }
}

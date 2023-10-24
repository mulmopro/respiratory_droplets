/*
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,C)")
    )
);
*/
{
    fvScalarMatrix CEqn
    (
        //fvm::ddt(rho, C)
        fvm::div(phi, C)
        - fvm::laplacian(turbulence->muEff()/Sc, C)
        ==
        //reaction->R(Yi)
        fvOptions(rho, C)
    );

    CEqn.relax();

    fvOptions.constrain(CEqn);

    CEqn.solve(mesh.solver("C"));

    fvOptions.correct(C);
}

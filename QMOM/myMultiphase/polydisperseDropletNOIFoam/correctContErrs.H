contErr1 =
    fvc::ddt(alpha1, rho1) + fvc::div(alphaRhoPhi1);

contErr2 =
    fvc::ddt(alpha2, rho2) + fvc::div(alphaRhoPhi2);

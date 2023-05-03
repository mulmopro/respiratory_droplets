/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "myEvapULPNormDBCoBCPolydispersePhaseModel.H"
#include "fvMatrix.H"
#include "twoPhaseSystem.H"
#include "fixedValueFvPatchFields.H"
#include "cyclicFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "emptyFvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "momentFieldSets.H"
#include "vectorList.H"
#include "fixedFaceFvPatchScalarField.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::updateVelocity()
{
    // Correct mean velocity using the new velocity moments
    /*
    U_ =
        quadrature_.velocityMoments()[1]
       /Foam::max
        (
            quadrature_.moments()[1],
            residualAlpha_*rho()
        );
    U_.correctBoundaryConditions();

    phiPtr_() = fvc::flux(U_);
    alphaPhi_ = fvc::interpolate(*this)*phiPtr_();
    alphaRhoPhi_ = fvc::interpolate(rho())*alphaPhi_;
    */
   
    volVectorField& U = this->U();
    surfaceScalarField& phi = this->phi();
    surfaceScalarField& alphaPhi = this->alphaPhi();
    surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi();

    //***
    //Info<< "after phi" << endl;
    /*
    label facei = 0;
    Info<< "In velocity correction, before update" << endl;
    Info<< "phi.boundaryField()[0][" << facei << "] = " << phi.boundaryField()[0][facei] << endl;
    facei = 83;
    Info<< "phi.boundaryField()[2][" << facei << "] = " << phi.boundaryField()[2][facei] << endl;
    facei = 20077;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 20078;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19839;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    */
    //***
    /*
    U =
        quadrature_.velocityMoments()[1]
       /Foam::max
        (
            quadrature_.moments()[1],
            this->residualAlpha()*this->rho()
        );
    */
    U = dimensionedVector("zero", dimVelocity, Foam::vector(0.0, 0.0, 0.0));
    forAll(Us_, nodei)
    {
        U += alphas_[nodei]*Us_[nodei];
    }
    U /= Foam::max((*this), this->residualAlpha());
    U.correctBoundaryConditions();

    phi = fvc::flux(this->U());
    alphaPhi = fvc::interpolate(*this)*phi;
    alphaRhoPhi = fvc::interpolate(this->rho())*alphaPhi;
    //***
    //Info<< "after phi" << endl;
    /*
    facei = 0;
    Info<< "In velocity correction, after update" << endl;
    Info<< "phi.boundaryField()[0][" << facei << "] = " << phi.boundaryField()[0][facei] << endl;
    facei = 83;
    Info<< "phi.boundaryField()[2][" << facei << "] = " << phi.boundaryField()[2][facei] << endl;
    facei = 20077;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 20078;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19839;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    */
    //***
}


template<class BasePhaseModel>
Foam::scalar Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::coalescenceSource
(
    const label celli,
    const label momentOrder
)
{
    scalar cSource = 0.0;
    if (!coalescence_)
    {
        return cSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volScalarNode& node1 = nodes[nodei];
        scalar weight1 = node1.primaryWeight()[celli];
        scalar abscissa1 = Foam::max(node1.primaryAbscissae()[0][celli], SMALL);
        scalar n1 = node1.n(celli, weight1, abscissa1);
        scalar d1 = node1.d(celli, abscissa1);

        forAll(nodes, nodej)
        {
            const volScalarNode& node2 = nodes[nodej];
            scalar weight2 = node2.primaryWeight()[celli];
            scalar abscissa2 = Foam::max(node2.primaryAbscissae()[0][celli], SMALL);

            scalar n2 = node2.n(celli, weight2, abscissa2);
            scalar d2 = node2.d(celli, abscissa2);
            vector Ur = Us_[nodei][celli] - Us_[nodej][celli];

            //- Diameter is used to calculate the coalesence kernel in place
            //  of the abscissa, handled inside kernel
            cSource +=
                0.5*n1
               *(
                    n2
                   *(
                        pow // Birth
                        (
                            (abscissa1) + (abscissa2),
                            momentOrder
                        )
                      - pow(abscissa1, momentOrder)
                      - pow(abscissa2, momentOrder)
                    )
                )*coalescenceKernel_->Ka(d1, d2, Ur, celli);
        }
    }
    return cSource;
}


template<class BasePhaseModel>
Foam::vector Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::coalescenceSourceU
(
    const label celli,
    const label momentOrder
)
{
    vector cSource = Zero;
    if (!coalescence_ || momentOrder == 1)
    {
        return cSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volScalarNode& node1 = nodes[nodei];
        scalar weight1 = node1.primaryWeight()[celli];
        scalar abscissa1 = Foam::max(node1.primaryAbscissae()[0][celli], SMALL);
        scalar n1 = node1.n(celli, weight1, abscissa1);
        scalar d1 = node1.d(celli, abscissa1);

        forAll(nodes, nodej)
        {
            const volScalarNode& node2 = nodes[nodej];
            scalar weight2 = node2.primaryWeight()[celli];
            scalar abscissa2 = Foam::max(node2.primaryAbscissae()[0][celli], SMALL);

            scalar n2 = node2.n(celli, weight2, abscissa2);
            scalar d2 = node2.d(celli, abscissa2);
            vector Ur = Us_[nodei][celli] - Us_[nodej][celli];

            //- Diameter is used to calculate the coalesence kernel in place
            //  of the abscissa, handled inside kernel
            cSource +=
                0.5*n1*Us_[nodei][celli]
               *(
                    n2
                   *(
                        pow // Birth
                        (
                            (abscissa1) + (abscissa2),
                            momentOrder
                        )
                      - pow(abscissa1, momentOrder)
                      - pow(abscissa2, momentOrder)
                    )*coalescenceKernel_->Ka(d1, d2, Ur, celli)
                );
        }
    }

    return cmptMultiply(cSource, validDirections_);
}


template<class BasePhaseModel>
Foam::scalar Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::breakupSource
(
    const label celli,
    const label momentOrder
)
{
    scalar bSource = 0.0;
    if (!breakup_)
    {
        return bSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volScalarNode& node = nodes[nodei];
        scalar weight = node.primaryWeight()[celli];
        scalar abscissa = Foam::max(node.primaryAbscissae()[0][celli], SMALL);
        scalar d = node.d(celli, abscissa);
        scalar n = node.n(celli, weight, abscissa);

        //- Diameter is used to calculate the breakup kernel in place
        //  of the abscissa
        bSource +=
            n
           *breakupKernel_->Kb(d, celli)
           *breakupKernel_->massNodeSource  //Birth and death
            (
                abscissa,
                momentOrder
            );
    }
    return bSource;
}


template<class BasePhaseModel>
Foam::vector Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::breakupSourceU
(
    const label celli,
    const label momentOrder
)
{
    vector bSource = Zero;
    if (!breakup_  || momentOrder == 1)
    {
        return bSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    forAll(nodes, nodei)
    {
        const volScalarNode& node = nodes[nodei];
        scalar weight = node.primaryWeight()[celli];
        scalar abscissa = Foam::max(node.primaryAbscissae()[0][celli], SMALL);
        scalar d = node.d(celli, abscissa);
        scalar n = node.n(celli, weight, abscissa);

        //- Diameter is used to calculate the breakup kernel in place
        //  of the abscissa
        bSource +=
            n*Us_[nodei][celli]
           *breakupKernel_->Kb(d, celli)
           *breakupKernel_->massNodeSource  //Birth and death
            (
                abscissa,
                momentOrder
            );
    }

    return cmptMultiply(bSource, validDirections_);
}


template<class BasePhaseModel>
void Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::solveSourceOde()
{
    if (!coalescence_ && !breakup_)
    {
        return;
    }

    coalescenceKernel_->preUpdate();
    breakupKernel_->preUpdate();

    volScalarMomentFieldSet& moments = quadrature_.moments();
    label nMoments = quadrature_.nMoments();
    PtrList<volVectorField>& Ups = quadrature_.velocityMoments();
    label nVelocityMoments = Ups.size();

    scalar globalDt = moments[0].mesh().time().deltaT().value();

    if (!solveOde_)
    {
        forAll(moments[0], celli)
        {
            forAll(moments, mi)
            {
                moments[mi][celli] +=
                    globalDt
                   *(
                        coalescenceSource(celli, mi)
                        + breakupSource(celli, mi)
                    );
            }

            forAll(Ups, mi)
            {
                Ups[mi][celli] +=
                    globalDt
                   *(
                        coalescenceSourceU(celli, mi)
                      + breakupSourceU(celli, mi)
                    );
            }
        }
        return;
    }


    Info << "Solving source terms in realizable ODE solver." << endl;

    forAll(moments[0], celli)
    {
        if ((*this)[celli] < 0.01)
        {
            continue;
        }

        scalarField oldMoments(nMoments, Zero);
        vectorField oldUps(nVelocityMoments, Zero);
        forAll(oldMoments, mi)
        {
            oldMoments[mi] = moments[mi][celli];
        }
        forAll(oldUps, mi)
        {
            oldUps[mi] = Ups[mi][celli];
        }

        //- Local time
        scalar localT = 0.0;

        // Initialize the local step
        scalar localDt = localDt_[celli];

        // Initialize RK parameters
        scalarField k1(nMoments_, Zero);
        scalarField k2(nMoments_, Zero);
        scalarField k3(nMoments_, Zero);

        vectorField k1U(nVelocityMoments, Zero);
        vectorField k2U(nVelocityMoments, Zero);
        vectorField k3U(nVelocityMoments, Zero);

        // Flag to indicate if the time step is complete
        bool timeComplete = false;

        // Check realizability of intermediate moment sets
        bool realizableUpdate1 = true;
        bool realizableUpdate2 = true;
        bool realizableUpdate3 = true;

        scalarList momentsSecondStep(nMoments_, 0.0);
        scalar error = 0.0;

        while (!timeComplete)
        {
            bool nullSource = true;
            do
            {
                nullSource = true;

                // First intermediate update
                forAll(oldMoments, mi)
                {
                    k1[mi] =
                        localDt
                       *(
                            coalescenceSource(celli, mi)
                          + breakupSource(celli, mi)
                        );
                    moments[mi][celli] = oldMoments[mi] + k1[mi];
                    nullSource =
                    (
                        mag(k1[mi]) < pow(SMALL, mi)
                     && nullSource
                    );
                }

                if (nullSource)
                {
                    break;
                }

                forAll(oldUps, mi)
                {
                    k1U[mi] =
                        localDt
                       *(
                            coalescenceSourceU(celli, mi)
                          + breakupSourceU(celli, mi)
                        );
                    Ups[mi][celli] = oldUps[mi] + k1U[mi];
                }
                realizableUpdate1 =
                        quadrature_.updateAllLocalQuadrature(celli, false);

                // Second moment update
                forAll(oldMoments, mi)
                {
                    k2[mi] =
                        localDt
                       *(
                            coalescenceSource(celli, mi)
                          + breakupSource(celli, mi)
                        );
                    moments[mi][celli] = oldMoments[mi] + (k1[mi] + k2[mi])/4.0;
                    momentsSecondStep[mi] = moments[mi][celli];
                }
                forAll(oldUps, mi)
                {
                    k2U[mi] =
                        localDt
                       *(
                            coalescenceSourceU(celli, mi)
                          + breakupSourceU(celli, mi)
                        );
                    Ups[mi][celli] = oldUps[mi] + (k1U[mi] + k2U[mi])/4.0;
                }
                realizableUpdate2 =
                        quadrature_.updateAllLocalQuadrature(celli, false);

                // Third moment update
                forAll(oldMoments, mi)
                {
                    k3[mi] =
                        localDt
                       *(
                            coalescenceSource(celli, mi)
                          + breakupSource(celli, mi)
                        );
                    moments[mi][celli] =
                        oldMoments[mi] + (k1[mi] + k2[mi] + 4.0*k3[mi])/6.0;
                }
                forAll(oldUps, mi)
                {
                    k3U[mi] =
                        localDt
                       *(
                            coalescenceSourceU(celli, mi)
                          + breakupSourceU(celli, mi)
                        );
                    Ups[mi][celli] =
                        oldUps[mi] + (k1U[mi] + k2U[mi] + 4.0*k3U[mi])/6.0;
                }
                realizableUpdate3 =
                        quadrature_.updateAllLocalQuadrature(celli, false);

                if
                (
                    !realizableUpdate1
                 || !realizableUpdate2
                 || !realizableUpdate3
                )
                {
                    Info << "Not realizable" << endl;

                    forAll(oldMoments, mi)
                    {
                        moments[mi][celli] = oldMoments[mi];
                    }
                    forAll(oldUps, mi)
                    {
                        Ups[mi][celli] = oldUps[mi];
                    }

                    // Updating local quadrature with old moments
                    quadrature_.updateAllLocalQuadrature(celli);

                    localDt /= 2.0;

                    if (localDt < minLocalDt_)
                    {
                        FatalErrorInFunction
                            << "Reached minimum local step in realizable ODE"
                            << nl
                            << "    solver. Cannot ensure realizability." << nl
                            << abort(FatalError);
                    }
                }
            }
            while
            (
                !realizableUpdate1
             || !realizableUpdate2
             || !realizableUpdate3
            );

            if (nullSource)
            {
                timeComplete = true;
                localT = 0.0;
                break;
            }

            for (label mi = 0; mi < nMoments; mi++)
            {
                scalar scalei =
                   Foam::max
                    (
                        mag(momentsSecondStep[mi]), mag(oldMoments[mi])
                    )*RTol_;

                if (scalei > 0)
                {
                    error +=
                        sqr
                        (
                            (momentsSecondStep[mi] - moments[mi][celli])/scalei
                        );
                }
            }
            error = sqrt(error/nMoments);
            if (error < SMALL)
            {
                timeComplete = true;
                localT = 0.0;
                break;
            }
            else if (error < 1)
            {
                localT += localDt;
                localDt_[celli] = localDt;

                localDt *= Foam::min(facMax_, Foam::max(facMin_, fac_/pow(error, 1.0/3.0)));

                scalar maxLocalDt = Foam::max(globalDt - localT, 0.0);
                localDt = Foam::min(maxLocalDt, localDt);

                forAll(oldMoments, mi)
                {
                    oldMoments[mi] = moments[mi][celli];
                }
                forAll(oldUps, mi)
                {
                    oldUps[mi] = Ups[mi][celli];
                }

                if (localDt == 0.0)
                {
                    timeComplete = true;
                    localT = 0.0;
                    break;
                }
            }
            else
            {
                localDt *= Foam::min(1.0, Foam::max(facMin_, fac_/pow(error, 1.0/3.0)));

                forAll(oldMoments, mi)
                {
                    moments[mi][celli] = oldMoments[mi];
                }
                forAll(oldUps, mi)
                {
                    Ups[mi][celli] = oldUps[mi];
                }
                quadrature_.updateAllLocalQuadrature(celli);

                if (localDt < minLocalDt_)
                {
                    Info<<" cell "<<celli
                        <<" error: " <<error
                        <<" dt: "<<localDt<<endl;
                    FatalErrorInFunction
                        << "Reached minimum local step in realizable ODE"
                        << nl
                        << "    solver. Cannot ensure realizability." << nl
                        << abort(FatalError);
                }
            }
        }
    }

    quadrature_.updateAllMoments();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::myEvapULPNormDBCoBCPolydispersePhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseProperties, phaseName),
    pbeDict_
    (
        IOobject
        (
            "populationBalanceProperties",
            fluid.mesh().time().constant(),
            fluid.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    solveOde_(pbeDict_.lookupOrDefault("ode", false)),
    coalescence_(pbeDict_.lookup("coalescence")),
    breakup_(pbeDict_.lookup("breakup")),
    quadrature_(phaseName, fluid.mesh(), "RPlus"),
    nNodes_(quadrature_.nodes().size()),
    nMoments_(quadrature_.nMoments()),
    alphas_(nNodes_),
    dotMsLim_(nNodes_),
    Us_(quadrature_.velocities()),
    Vs_(nNodes_),
    ds_(nNodes_),
    maxD_("maxD", dimLength, this->phaseDict()),
    minD_("minD", dimLength, this->phaseDict()),
    W_("W", dimless, this->phaseDict()),
    Xi_("Xi", dimless, this->phaseDict()),
    minLocalDt_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("minLocalDt"))),
    localDt_(this->size(), fluid.mesh().time().deltaT().value()/10.0),
    ATol_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("ATol"))),
    RTol_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("RTol"))),
    facMax_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("facMax"))),
    facMin_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("facMin"))),
    fac_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("fac"))),
    validDirections_
    (
        (vector(this->fluid().mesh().solutionD()) + vector(1.0, 1.0, 1.0))/2.0
    )
{
    this->d_.writeOpt() = IOobject::AUTO_WRITE;

    const volVectorField& U(this->U());
    /*
    wordList phiTypes
    (
        U.boundaryField().size(),
        calculatedFvPatchScalarField::typeName
    );

    forAll(U.boundaryField(), i)
    {
        if
        (
            isA<fixedValueFvPatchVectorField>(U.boundaryField()[i])
         || isA<slipFvPatchVectorField>(U.boundaryField()[i])
         || isA<partialSlipFvPatchVectorField>(U.boundaryField()[i])
        )
        {
            phiTypes[i] = fixedValueFvsPatchScalarField::typeName;
        }
    }
    */

    forAll(alphas_, nodei)
    {
        alphas_.set
        (
            nodei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "alpha",
                        IOobject::groupName
                        (
                            this->name(),
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar("alpha", dimless, 0.0)
            )
        );

        dotMsLim_.set
        (
            nodei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "dotMsLim",
                        IOobject::groupName
                        (
                            this->name(),
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar("dotMsLim", dimMass/dimTime/dimVolume, 0.0)
            )
        );

        Vs_.set
        (
            nodei,
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "V",
                        IOobject::groupName
                        (
                            this->name(),
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid.mesh(),
                dimensionedVector("zeroV", dimVelocity, Zero)
            )
        );

        ds_.set
        (
            nodei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "d",
                        IOobject::groupName
                        (
                            this->name(),
                            Foam::name(nodei)
                        )
                    ),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluid.mesh(),
                minD_
            )
        );
    }

    // Set alpha value based on moments
    //volScalarField(*this) == quadrature_.moments()[1]/this->rho();
    //volScalarField(*this) == quadrature_.moments()[3]*Foam::constant::mathematical::pi/6.0;
    volScalarField(*this) == quadrature_.moments()[3]*(W_*pow(Xi_, 3.0))*Foam::constant::mathematical::pi/6.0;
    correct();

    const dictionary& pimpleDict =
        this->fluid().mesh().solutionDict().subDict("PIMPLE");
    label nCorrectors = pimpleDict.lookupOrDefault<label>("nFluxCorrectors", 0);
    if (nCorrectors > 0)
    {
        word patchName
        (
            pimpleDict.lookupOrDefault
            (
                "corrPatch",
                U.boundaryField()[0].patch().name()
            )
        );
        wordList boundaries(U.boundaryField().size(), "zeroGradient");
        forAll(boundaries, patchi)
        {
            if (U.boundaryField()[patchi].patch().name() == patchName)
            {
                //boundaries[patchi] = "fixedFace";
                boundaries[patchi] = "fixedValue";
            }
        }

        corr_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "corr",
                    this->fluid().mesh().time().timeName(),
                    this->fluid().mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->fluid().mesh(),
                dimensionedScalar("0", sqr(dimLength)/dimTime, 0.0),
                boundaries
            )
        );
    }
}


template<class BasePhaseModel>
void Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::setModels()
{
    coalescenceKernel_.set
    (
        new populationBalanceSubModels::aggregationKernels::coalescence
        (
            pbeDict_.subDict("coalescenceKernel"),
            this->fluid().mesh()
        )
    );
    breakupKernel_ =
    (
        Foam::populationBalanceSubModels::breakupKernel::New
        (
            pbeDict_.subDict("breakupKernel"),
            this->fluid().mesh()
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::~myEvapULPNormDBCoBCPolydispersePhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::correct()
{

    volScalarField& d = this->d();
    if (nNodes_ == 1)
    {
        alphas_[0] = *this;
        ds_[0] =
            Foam::min
            (
                Foam::max
                (
                    quadrature_.nodes()[0].primaryAbscissae()[0]*Xi_,
                    minD_
                ),
                maxD_
            );
        d = ds_[0];
        return;
    }
    else
    {
        d = dimensionedScalar("zero", dimLength, 0.0);
        volScalarField scale
        (
            (*this)
           /Foam::max
            (
                quadrature_.moments()[3],
                this->residualAlpha()*6.0/Foam::constant::mathematical::pi/(W_*pow(Xi_, 3.0))
            )
        );

        forAll(quadrature_.nodes(), nodei)
        {
            const volScalarNode& node = quadrature_.nodes()[nodei];

            // Set alpha values such that the moment.1 is equal to the bounded
            // alpha
            alphas_[nodei] =
                node.primaryWeight()*(
                    pow
                    (
                        node.primaryAbscissae()[0],
                        3
                    )
                )*scale;

            alphas_[nodei].max(0);
            alphas_[nodei].min(1);

            //  Calculate bubble diameter based on bubble mass (abscissa)
            ds_[nodei] =
                Foam::min
                (
                    Foam::max
                    (
                        node.primaryAbscissae()[0]*Xi_,
                        minD_
                    ),
                    maxD_
                );
            d += alphas_[nodei]*ds_[nodei];
        }

        d /= Foam::max((*this), this->residualAlpha());
        d.max(minD_);
    }
}


template<class BasePhaseModel>
void Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::relativeTransport()
{
    // Do nothing if only mean is used
    if (nNodes_ == 1)
    {
        return;
    }

    Info<< "Transporting moments based on relative flux" << endl;

    quadrature_.interpolateNodes();
    const PtrList<surfaceScalarNode>& nodesOwn = quadrature_.nodesOwn();
    const PtrList<surfaceScalarNode>& nodesNei = quadrature_.nodesNei();

    // ***
    /*
        label celli = 11760;
        Info<< "Before solving Mp, at celli: " << celli << endl;
        //Info<< "m0>1 and m1<0, at celli:" << celli << endl;
        Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
            << ", " << quadrature_.moments()[1][celli]
            << ", " << quadrature_.moments()[2][celli]
            << ", " << quadrature_.moments()[3][celli]
            << endl;
        */
    // ***

    // Transport moments with relative flux
    forAll(quadrature_.moments(), mEqni)
    {
        volScalarField& m = quadrature_.moments()[mEqni];
        dimensionedScalar zeroPhi("zero", this->phiPtr_().dimensions(), 0.0);

        // Create total flux field so that the individual fluxes can be
        // summed together
        volScalarField relativeDivVp
        (
            IOobject
            (
                "relativeDivVp",
                this->fluid().mesh().time().timeName(),
                this->fluid().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->fluid().mesh(),
            dimensionedScalar("zero", m.dimensions()/dimTime, 0.0)
        );

        for (label nodei = 0; nodei < nNodes_; nodei++)
        {
            surfaceScalarField phiv("phiv", fvc::flux(Vs_[nodei]));

            // Calculate size moment flux
            surfaceScalarField rFluxVp
            (
                nodesNei[nodei].primaryWeight()
               *(
                    pow
                    (
                        nodesNei[nodei].primaryAbscissae()[0],
                        mEqni
                    )
                )*Foam::min(phiv, zeroPhi)
              + nodesOwn[nodei].primaryWeight()
               *pow
                (
                    nodesOwn[nodei].primaryAbscissae()[0],
                    mEqni
                )*Foam::max(phiv, zeroPhi)
            );

            relativeDivVp += fvc::surfaceIntegrate(rFluxVp);

            // ***
            /*
            Info<< "Relative flux for nodei: " << nodei << endl 
                << "facei = 23184: " << phiv[23184] << endl
                << "facei = 23422: " << phiv[23422] << endl
                << "facei = 23423: " << phiv[23423] << endl
                << endl;
            Info<< "relativeDivVp for nodei: " << nodei << endl 
                << "celli = 11760: " << relativeDivVp[11760] << endl;

            Info<< "nodesnei for nodei: " << nodei << endl
                << "facei = 23184: w = " << nodesNei[nodei].primaryWeight()[23184] 
                << "xi = " << nodesNei[nodei].primaryAbscissae()[0][23184] << endl
                << "facei = 23422: w = " << nodesNei[nodei].primaryWeight()[23422] 
                << "xi = " << nodesNei[nodei].primaryAbscissae()[0][23422] << endl
                << "facei = 23423: w = " << nodesNei[nodei].primaryWeight()[23423] 
                << "xi = " << nodesNei[nodei].primaryAbscissae()[0][23423] << endl
                << endl;
            Info<< "nodesown for nodei: " << nodei << endl
                << "facei = 23184: w = " << nodesOwn[nodei].primaryWeight()[23184] 
                << "xi = " << nodesOwn[nodei].primaryAbscissae()[0][23184] << endl
                << "facei = 23422: w = " << nodesOwn[nodei].primaryWeight()[23422] 
                << "xi = " << nodesOwn[nodei].primaryAbscissae()[0][23422] << endl
                << "facei = 23423: w = " << nodesOwn[nodei].primaryWeight()[23423] 
                << "xi = " << nodesOwn[nodei].primaryAbscissae()[0][23423] << endl
                << endl;
            */
            // ***
        }

        // Solve relative size moment transport equation
        fvScalarMatrix mEqn
        (
            fvm::ddt(m)
          + relativeDivVp
        );

        // ***
        //label celli = 11760;
        /*
        Info<< "Before solving Mp" << endl;
        //Info<< "m0>1 and m1<0, at celli:" << celli << endl;
        Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
            << ", " << quadrature_.moments()[1][celli]
            << ", " << quadrature_.moments()[2][celli]
            << ", " << quadrature_.moments()[3][celli]
            << endl;
        */
        // ***

        mEqn.relax();
        mEqn.solve();

        // ****
        /*
        celli = 11760;
        Info<< "mi: " << mEqni << ", relativeDivVp: " << relativeDivVp[celli] << endl;
        */
        /*
        celli = 11760;
        Info<< "After solving Mp" << endl;
        Info<< "mi: " << mEqni << ", relativeDivVp: " << relativeDivVp[celli] << endl;
        Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
            << ", " << quadrature_.moments()[1][celli]
            << ", " << quadrature_.moments()[2][celli]
            << ", " << quadrature_.moments()[3][celli]
            << endl;
        */
        forAll (quadrature_.moments()[0], celli)
        {
            if (quadrature_.moments()[0][celli] > 1.0)
            {
                if (quadrature_.moments()[1][celli] < 0.0)
                {
                    Info<< "After solving Mp but before Up" << endl;
                    Info<< "m0>1 and m1<0, at celli:" << celli << endl;
                    Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
                        << ", " << quadrature_.moments()[1][celli]
                        << ", " << quadrature_.moments()[2][celli]
                        << ", " << quadrature_.moments()[3][celli]
                        << endl;
                    Info<< "mi: " << mEqni << ", relativeDivVp: " << relativeDivVp[celli] << endl;
                }
            }
        }
        // ****
    }

    // ****
    /*
        celli = 11760;
        Info<< "After solving Mp, at celli: " << celli << endl;
        //Info<< "mi: " << mEqni << ", relativeDivVp: " << relativeDivVp[celli] << endl;
        Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
            << ", " << quadrature_.moments()[1][celli]
            << ", " << quadrature_.moments()[2][celli]
            << ", " << quadrature_.moments()[3][celli]
            << endl; 
    */   
    // ***

    forAll(quadrature_.velocityMoments(), mEqni)
    {
        volVectorField& Up = quadrature_.velocityMoments()[mEqni];
        dimensionedScalar zeroPhi("zero", this->phiPtr_().dimensions(), 0.0);

        // Create total flux field so that the individual fluxes can be
        // summed together
        volVectorField relativeDivPp
        (
            IOobject
            (
                "relativeDivPp",
                this->fluid().mesh().time().timeName(),
                this->fluid().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->fluid().mesh(),
            dimensionedVector("zero", Up.dimensions()/dimTime, Zero)
        );

        for (label nodei = 0; nodei < nNodes_; nodei++)
        {
            surfaceScalarField phiv("phiv", fvc::flux(Vs_[nodei]));

            // Calculate velocity moment flux
            surfaceVectorField rFluxPp
            (
                "rFluxPp",
                quadrature_.velocitiesNei()[nodei]
               *nodesNei[nodei].primaryWeight()
               *(
                    pow
                    (
                        nodesNei[nodei].primaryAbscissae()[0],
                        mEqni
                    )
                )*Foam::min(phiv, zeroPhi)
              + quadrature_.velocitiesOwn()[nodei]
               *nodesOwn[nodei].primaryWeight()
               *pow
                (
                    nodesOwn[nodei].primaryAbscissae()[0],
                    mEqni
                )*Foam::max(phiv, zeroPhi)
            );

            relativeDivPp += fvc::surfaceIntegrate(rFluxPp);
        }

        // Solve relative velocity moment transport equation
        fvVectorMatrix UpEqn
        (
            fvm::ddt(Up)
          + relativeDivPp
        );

        UpEqn.relax();
        UpEqn.solve();
    }

    // ****
    forAll (quadrature_.moments()[0], celli)
    {
        if (quadrature_.moments()[0][celli] > 1.0)
        {
            if (quadrature_.moments()[1][celli] < 0.0)
            {
                Info<< "After solving Up" << endl;
                Info<< "m0>1 and m1<0, at celli:" << celli << endl;
                Info<< "M0 and M1 and M2: " << quadrature_.moments()[0][celli] 
                    << ", " << quadrature_.moments()[1][celli]
                    << ", " << quadrature_.moments()[2][celli]
                    << endl;
            }
        }
    }
    // ****

    quadrature_.updateAllQuadrature();
    correct();

    // calc dotMsLim
    PtrList<volScalarNode>& nodes = quadrature_.nodes();
    label sizeIndex = nodes[0].sizeIndex();
    //bool lengthBased = nodes[0].lengthBased();
    //- GlobalDt
    scalar globalDt = quadrature_.moments()[0].mesh().time().deltaT().value();
    //- Minimum abscissa
    scalar minAbscissa = minD_.value() / Xi_.value();

    //scalar minM3 = this->residualAlpha().value() * 6.0 / Foam::constant::mathematical::pi / (W_.value()*pow(Xi_.value(), 3));

    forAll(nodes, pNodeI)
    {
        forAll(quadrature_.moments()[1], celli)
        {
            //if (quadrature_.moments()[3][celli] >= minM3)
            if (nodes[pNodeI].primaryAbscissae()[sizeIndex][celli] >= minAbscissa)
            {
                volScalarNode& node = nodes[pNodeI];
                scalar& abscissa = node.primaryAbscissae()[sizeIndex][celli];
                scalar& weight = node.primaryWeight()[celli];

                dotMsLim_[pNodeI][celli] = W_.value()*pow(Xi_.value(), 3) * this->rho()[celli] * Foam::constant::mathematical::pi * weight * pow(abscissa, 2) * (abscissa - minAbscissa) / (2.0 * globalDt);
            }
        }
    }
}


template<class BasePhaseModel>
void Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::averageTransport
(
    const PtrList<fvVectorMatrix>& AEqns
)
{
    // Evaporation
    //Info<< "Evaporation in averageTransport" << endl;

    PtrList<volScalarNode>& nodes = quadrature_.nodes();
    label sizeIndex = nodes[0].sizeIndex();
    //bool lengthBased = nodes[0].lengthBased();
    //- GlobalDt
    scalar globalDt = quadrature_.moments()[0].mesh().time().deltaT().value();
    //- Minimum abscissa
    scalar minAbscissa = minD_.value() / Xi_.value();

    scalar minM3 = this->residualAlpha().value() * 6.0 / Foam::constant::mathematical::pi / (W_.value()*pow(Xi_.value(), 3));

    //const volScalarField& dvdtf = quadrature_.moments()[0].mesh().lookupObject<volScalarField>("dvdtf");

    //const volScalarField& alpha_liquid = quadrature_.moments()[0].mesh().lookupObject<volScalarField>("alpha.liquid");

    //dimensionedScalar MSmall("mSmall", inv(dimArea), 1e-8);
    //dimensionedScalar ASmall("ASmall", dimless, 1e-8);

    forAll(nodes, pNodeI)
    {
        const volScalarField& dvdtf = quadrature_.moments()[0].mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                "dvdtf",
                Foam::name(pNodeI)
            )
        );
        volScalarField drift = 2 * dvdtf / (Foam::constant::mathematical::pi) / (W_*pow(Xi_, 3));

        forAll(quadrature_.moments()[1], celli)
        {
            if (quadrature_.moments()[3][celli] >= minM3)
            {
                volScalarNode& node = nodes[pNodeI];
                scalar& abscissa = node.primaryAbscissae()[sizeIndex][celli];
                scalar& weight = node.primaryWeight()[celli];

                if (abscissa > minAbscissa)
                {
                    abscissa = abscissa - drift[celli] / (weight * abscissa * abscissa) * globalDt;
                }

                if (abscissa <= minAbscissa)
                {
                    abscissa = minAbscissa;
                }
            }
        }
    }
    quadrature_.updateAllMoments();
    quadrature_.updateAllQuadrature();


    // Correct mean flux
    const PtrList<surfaceScalarNode>& nodesOwn = quadrature_.nodesOwn();
    const PtrList<surfaceScalarNode>& nodesNei = quadrature_.nodesNei();
    //***
    //Info<< "before zerophi" << endl;
    //***
    dimensionedScalar zeroPhi("zero", this->phiPtr_().dimensions(), 0.0);
    //***
    //Info<< "before phi" << endl;
    //***
    surfaceScalarField phi(this->phiPtr_());

    //***
    //Info<< "after phi" << endl;
    /*
    label facei = 0;
    Info<< "Before correction" << endl;
    Info<< "phi.boundaryField()[0][" << facei << "] = " << phi.boundaryField()[0][facei] << endl;
    facei = 83;
    Info<< "phi.boundaryField()[2][" << facei << "] = " << phi.boundaryField()[2][facei] << endl;
    facei = 20076;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 20077;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19839;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19601;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19837;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19839;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19840;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    */
    //***

    //Info<< "Before corr" <<  endl;

    const dictionary& pimpleDict =
        this->fluid().mesh().solutionDict().subDict("PIMPLE");
    label nCorrectors = pimpleDict.lookupOrDefault<label>("nFluxCorrectors", 0);
    scalar phicLim = pimpleDict.lookupOrDefault<scalar>("phicLim", 1.0);
    if (corr_.valid())
    {
        volScalarField& corr = corr_.ref();

        for (label i = 0; i < nCorrectors; i++)
        {
            //- Update interpolated nodes since upwinding direction
            //  may have changed
            quadrature_.interpolateNodes();
            volScalarField meanM1Flux
            (
                IOobject
                (
                    "meanM1Flux",
                    this->fluid().mesh().time().timeName(),
                    this->fluid().mesh()
                ),
                this->fluid().mesh(),
                dimensionedScalar("zero", inv(dimTime), Zero)
            );

            for (label nodei = 0; nodei < nNodes_; nodei++)
            {
                meanM1Flux +=
                    fvc::surfaceIntegrate
                    (
                        nodesNei[nodei].primaryWeight()
                       *(
                            pow
                            (
                                nodesNei[nodei].primaryAbscissae()[0],
                                3.0
                            )
                       )
                       *Foam::min(phi, zeroPhi)
                      + nodesOwn[nodei].primaryWeight()
                       *(
                            pow
                            (
                                nodesOwn[nodei].primaryAbscissae()[0],
                                3.0
                            )
                       )
                       *Foam::max(phi, zeroPhi)
                    );
            }

            dimensionedScalar minCoeff
            (
                dimensionedScalar::lookupOrDefault
                (
                    "minCoeff",
                    pimpleDict,
                    dimensionSet(0,0,0,0,0,0,0),
                    this->residualAlpha().value()
                )
            );

            fvScalarMatrix corrEqn
            (
                ((*this)*6.0/Foam::constant::mathematical::pi/(W_*pow(Xi_, 3)) - quadrature_.moments()[3])
               /this->fluid().mesh().time().deltaT()
              + meanM1Flux
              + fvm::laplacian
                (
                    Foam::max(quadrature_.moments()[3], minCoeff),
                    corr,
                    "laplacian(" + quadrature_.moments()[3].name() + ",corr)"
                )
            );
            corrEqn.relax();
            corrEqn.solve();
            surfaceScalarField phic(fvc::snGrad(corr)*this->fluid().mesh().magSf());
            forAll(phi, facei)
            {
                if (phic[facei] >= 0.0)
                {
                    phic[facei] = min(phic[facei], phicLim*mag(phi[facei]));
                }
                else
                {
                    phic[facei] = max(phic[facei], -phicLim*mag(phi[facei]));
                }
            }
            phi += phic;

            //***
            /*
            facei = 34535;
            Info<< "snGrad(corr)[" << facei << "] = " << fvc::snGrad(corr)[facei] << endl;
            facei = 34500;
            Info<< "snGrad(corr)[" << facei << "] = " << fvc::snGrad(corr)[facei] << endl;
            */
            //***
        }
    }
    quadrature_.interpolateNodes();

    //Info<< "After corr" <<  endl;

    // ****
    /*
        label celli = 11760;
        Info<< "In average Transport, before solving Mp, at celli: " << celli << endl;
        //Info<< "mi: " << mEqni << ", relativeDivVp: " << relativeDivVp[celli] << endl;
        Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
            << ", " << quadrature_.moments()[1][celli]
            << ", " << quadrature_.moments()[2][celli]
            << ", " << quadrature_.moments()[3][celli]
            << endl;    
    */
    //****
    /*
    label celli = 9961;
    Info<< "In average Transport, before solving Mp, at celli: " << celli << endl;
    Info<< "alpha[" << celli << "] = " << (*this)[celli] << endl;
    Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
        << ", " << quadrature_.moments()[1][celli]
        << ", " << quadrature_.moments()[2][celli]
        << ", " << quadrature_.moments()[3][celli]
        << endl;
    //Info<< "corr[" << celli << "] = " << corr_[celli] << endl;
    facei = 19601;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl
        << "nodesNei[0].primaryWeight()[" << facei << "] = " << nodesNei[0].primaryWeight()[facei] << endl
        << "nodesNei[0].primaryAbscissae()[0][" << facei << "] = " << nodesNei[0].primaryAbscissae()[0][facei] << endl
        << "nodesNei[1].primaryWeight()[" << facei << "] = " << nodesNei[1].primaryWeight()[facei] << endl
        << "nodesNei[1].primaryAbscissae()[0][" << facei << "] = " << nodesNei[1].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[0].primaryWeight()[" << facei << "] = " << nodesOwn[0].primaryWeight()[facei] << endl
        << "nodesOwn[0].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[0].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[1].primaryWeight()[" << facei << "] = " << nodesOwn[1].primaryWeight()[facei] << endl
        << "nodesOwn[1].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[1].primaryAbscissae()[0][facei]<< endl;

    facei = 19837;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl
        << "nodesNei[0].primaryWeight()[" << facei << "] = " << nodesNei[0].primaryWeight()[facei] << endl
        << "nodesNei[0].primaryAbscissae()[0][" << facei << "] = " << nodesNei[0].primaryAbscissae()[0][facei] << endl
        << "nodesNei[1].primaryWeight()[" << facei << "] = " << nodesNei[1].primaryWeight()[facei] << endl
        << "nodesNei[1].primaryAbscissae()[0][" << facei << "] = " << nodesNei[1].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[0].primaryWeight()[" << facei << "] = " << nodesOwn[0].primaryWeight()[facei] << endl
        << "nodesOwn[0].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[0].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[1].primaryWeight()[" << facei << "] = " << nodesOwn[1].primaryWeight()[facei] << endl
        << "nodesOwn[1].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[1].primaryAbscissae()[0][facei]<< endl;

    facei = 19839;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl
        << "nodesNei[0].primaryWeight()[" << facei << "] = " << nodesNei[0].primaryWeight()[facei] << endl
        << "nodesNei[0].primaryAbscissae()[0][" << facei << "] = " << nodesNei[0].primaryAbscissae()[0][facei] << endl
        << "nodesNei[1].primaryWeight()[" << facei << "] = " << nodesNei[1].primaryWeight()[facei] << endl
        << "nodesNei[1].primaryAbscissae()[0][" << facei << "] = " << nodesNei[1].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[0].primaryWeight()[" << facei << "] = " << nodesOwn[0].primaryWeight()[facei] << endl
        << "nodesOwn[0].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[0].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[1].primaryWeight()[" << facei << "] = " << nodesOwn[1].primaryWeight()[facei] << endl
        << "nodesOwn[1].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[1].primaryAbscissae()[0][facei]<< endl;

    facei = 19840;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl
        << "nodesNei[0].primaryWeight()[" << facei << "] = " << nodesNei[0].primaryWeight()[facei] << endl
        << "nodesNei[0].primaryAbscissae()[0][" << facei << "] = " << nodesNei[0].primaryAbscissae()[0][facei] << endl
        << "nodesNei[1].primaryWeight()[" << facei << "] = " << nodesNei[1].primaryWeight()[facei] << endl
        << "nodesNei[1].primaryAbscissae()[0][" << facei << "] = " << nodesNei[1].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[0].primaryWeight()[" << facei << "] = " << nodesOwn[0].primaryWeight()[facei] << endl
        << "nodesOwn[0].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[0].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[1].primaryWeight()[" << facei << "] = " << nodesOwn[1].primaryWeight()[facei] << endl
        << "nodesOwn[1].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[1].primaryAbscissae()[0][facei]<< endl;


    celli = 10080;
    Info<< "In average Transport, before solving Mp, at celli: " << celli << endl;
    Info<< "alpha[" << celli << "] = " << (*this)[celli] << endl;
    Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
        << ", " << quadrature_.moments()[1][celli]
        << ", " << quadrature_.moments()[2][celli]
        << ", " << quadrature_.moments()[3][celli]
        << endl;
    //Info<< "corr[" << celli << "] = " << corr_[celli] << endl;
    facei = 19839;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl
        << "nodesNei[0].primaryWeight()[" << facei << "] = " << nodesNei[0].primaryWeight()[facei] << endl
        << "nodesNei[0].primaryAbscissae()[0][" << facei << "] = " << nodesNei[0].primaryAbscissae()[0][facei] << endl
        << "nodesNei[1].primaryWeight()[" << facei << "] = " << nodesNei[1].primaryWeight()[facei] << endl
        << "nodesNei[1].primaryAbscissae()[0][" << facei << "] = " << nodesNei[1].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[0].primaryWeight()[" << facei << "] = " << nodesOwn[0].primaryWeight()[facei] << endl
        << "nodesOwn[0].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[0].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[1].primaryWeight()[" << facei << "] = " << nodesOwn[1].primaryWeight()[facei] << endl
        << "nodesOwn[1].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[1].primaryAbscissae()[0][facei]<< endl;

    facei = 20076;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl
        << "nodesNei[0].primaryWeight()[" << facei << "] = " << nodesNei[0].primaryWeight()[facei] << endl
        << "nodesNei[0].primaryAbscissae()[0][" << facei << "] = " << nodesNei[0].primaryAbscissae()[0][facei] << endl
        << "nodesNei[1].primaryWeight()[" << facei << "] = " << nodesNei[1].primaryWeight()[facei] << endl
        << "nodesNei[1].primaryAbscissae()[0][" << facei << "] = " << nodesNei[1].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[0].primaryWeight()[" << facei << "] = " << nodesOwn[0].primaryWeight()[facei] << endl
        << "nodesOwn[0].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[0].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[1].primaryWeight()[" << facei << "] = " << nodesOwn[1].primaryWeight()[facei] << endl
        << "nodesOwn[1].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[1].primaryAbscissae()[0][facei]<< endl;

    facei = 20077;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl
        << "nodesNei[0].primaryWeight()[" << facei << "] = " << nodesNei[0].primaryWeight()[facei] << endl
        << "nodesNei[0].primaryAbscissae()[0][" << facei << "] = " << nodesNei[0].primaryAbscissae()[0][facei] << endl
        << "nodesNei[1].primaryWeight()[" << facei << "] = " << nodesNei[1].primaryWeight()[facei] << endl
        << "nodesNei[1].primaryAbscissae()[0][" << facei << "] = " << nodesNei[1].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[0].primaryWeight()[" << facei << "] = " << nodesOwn[0].primaryWeight()[facei] << endl
        << "nodesOwn[0].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[0].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[1].primaryWeight()[" << facei << "] = " << nodesOwn[1].primaryWeight()[facei] << endl
        << "nodesOwn[1].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[1].primaryAbscissae()[0][facei]<< endl;

    facei = 0;
    Info<< "phi[" << facei << "] = " << phi.boundaryField()[0][facei] << endl
        << "nodesNei[0].primaryWeight()[" << facei << "] = " << nodesNei[0].primaryWeight().boundaryField()[0][facei] << endl
        << "nodesNei[0].primaryAbscissae()[0][" << facei << "] = " << nodesNei[0].primaryAbscissae()[0].boundaryField()[0][facei] << endl
        << "nodesNei[1].primaryWeight()[" << facei << "] = " << nodesNei[1].primaryWeight().boundaryField()[0][facei] << endl
        << "nodesNei[1].primaryAbscissae()[0][" << facei << "] = " << nodesNei[1].primaryAbscissae()[0].boundaryField()[0][facei] << endl
        << "nodesOwn[0].primaryWeight()[" << facei << "] = " << nodesOwn[0].primaryWeight().boundaryField()[0][facei] << endl
        << "nodesOwn[0].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[0].primaryAbscissae()[0].boundaryField()[0][facei] << endl
        << "nodesOwn[1].primaryWeight()[" << facei << "] = " << nodesOwn[1].primaryWeight().boundaryField()[0][facei] << endl
        << "nodesOwn[1].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[1].primaryAbscissae()[0].boundaryField()[0][facei]<< endl;

    celli = 9960;
    Info<< "In average Transport, before solving Mp, at celli: " << celli << endl;
    Info<< "alpha[" << celli << "] = " << (*this)[celli] << endl;
    Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
        << ", " << quadrature_.moments()[1][celli]
        << ", " << quadrature_.moments()[2][celli]
        << ", " << quadrature_.moments()[3][celli]
        << endl;
    //Info<< "corr[" << celli << "] = " << corr_[celli] << endl;
    facei = 19599;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl
        << "nodesNei[0].primaryWeight()[" << facei << "] = " << nodesNei[0].primaryWeight()[facei] << endl
        << "nodesNei[0].primaryAbscissae()[0][" << facei << "] = " << nodesNei[0].primaryAbscissae()[0][facei] << endl
        << "nodesNei[1].primaryWeight()[" << facei << "] = " << nodesNei[1].primaryWeight()[facei] << endl
        << "nodesNei[1].primaryAbscissae()[0][" << facei << "] = " << nodesNei[1].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[0].primaryWeight()[" << facei << "] = " << nodesOwn[0].primaryWeight()[facei] << endl
        << "nodesOwn[0].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[0].primaryAbscissae()[0][facei] << endl
        << "nodesOwn[1].primaryWeight()[" << facei << "] = " << nodesOwn[1].primaryWeight()[facei] << endl
        << "nodesOwn[1].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[1].primaryAbscissae()[0][facei]<< endl;

    facei = 83;
    Info<< "phi[" << facei << "] = " << phi.boundaryField()[2][facei] << endl
        << "nodesNei[0].primaryWeight()[" << facei << "] = " << nodesNei[0].primaryWeight().boundaryField()[2][facei] << endl
        << "nodesNei[0].primaryAbscissae()[0][" << facei << "] = " << nodesNei[0].primaryAbscissae()[0].boundaryField()[2][facei] << endl
        << "nodesNei[1].primaryWeight()[" << facei << "] = " << nodesNei[1].primaryWeight().boundaryField()[2][facei] << endl
        << "nodesNei[1].primaryAbscissae()[0][" << facei << "] = " << nodesNei[1].primaryAbscissae()[0].boundaryField()[2][facei] << endl
        << "nodesOwn[0].primaryWeight()[" << facei << "] = " << nodesOwn[0].primaryWeight().boundaryField()[2][facei] << endl
        << "nodesOwn[0].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[0].primaryAbscissae()[0].boundaryField()[2][facei] << endl
        << "nodesOwn[1].primaryWeight()[" << facei << "] = " << nodesOwn[1].primaryWeight().boundaryField()[2][facei] << endl
        << "nodesOwn[1].primaryAbscissae()[0][" << facei << "] = " << nodesOwn[1].primaryAbscissae()[0].boundaryField()[2][facei]<< endl;


    facei = 0;
    Info<< "phi.boundaryField()[0][" << facei << "] = " << phi.boundaryField()[0][facei] << endl;
    facei = 83;
    Info<< "phi.boundaryField()[2][" << facei << "] = " << phi.boundaryField()[2][facei] << endl;
    facei = 20076;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 20077;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19838;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19601;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19837;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19839;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    facei = 19840;
    Info<< "phi[" << facei << "] = " << phi[facei] << endl;
    */
    //****
    // ***

    // Mean moment advection
    Info<< "Transporting moments with average velocity" << endl;
    forAll(quadrature_.moments(), mEqni)
    {
        volScalarField& m = quadrature_.moments()[mEqni];

        volScalarField meanDivUbMp
        (
            IOobject
            (
                "meanDivUbMp",
                this->fluid().mesh().time().timeName(),
                this->fluid().mesh()
            ),
            this->fluid().mesh(),
            dimensionedScalar("zero", m.dimensions()/dimTime, Zero)
        );

        for (label nodei = 0; nodei < nNodes_; nodei++)
        {
            // Update average size moment flux
            surfaceScalarField aFluxMp
            (
                "aFluxMp",
                nodesNei[nodei].primaryWeight()
               *(
                    pow
                    (
                        nodesNei[nodei].primaryAbscissae()[0],
                        mEqni
                    )
                )*Foam::min(phi, zeroPhi)
              + nodesOwn[nodei].primaryWeight()
               *pow
                (
                    nodesOwn[nodei].primaryAbscissae()[0],
                    mEqni
                )*Foam::max(phi, zeroPhi)
            );

            meanDivUbMp += fvc::surfaceIntegrate(aFluxMp);
        }
        // ****
        /*
        celli = 9961;
        Info<< "meanDivUbMp[" << celli << "] = " << meanDivUbMp[celli] << endl;
        celli = 10080;
        Info<< "meanDivUbMp[" << celli << "] = " << meanDivUbMp[celli] << endl;
        celli = 10081;
        Info<< "meanDivUbMp[" << celli << "] = " << meanDivUbMp[celli] << endl;
        */
        // ****
        // Solve average size moment transport
        fvScalarMatrix mEqn
        (
            fvm::ddt(m)
          - fvc::ddt(m)
          + meanDivUbMp
        );
        mEqn.relax();
        mEqn.solve();
    }

    // ****
    /*
        celli = 11760;
        Info<< "In average Transport, after solving Mp, at celli: " << celli << endl;
        //Info<< "mi: " << mEqni << ", relativeDivVp: " << relativeDivVp[celli] << endl;
        Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
            << ", " << quadrature_.moments()[1][celli]
            << ", " << quadrature_.moments()[2][celli]
            << ", " << quadrature_.moments()[3][celli]
            << endl; 
    */   
    //****
    /*
    celli = 9961;
    Info<< "In average Transport, after solving Mp, at celli: " << celli << endl;
    Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
        << ", " << quadrature_.moments()[1][celli]
        << ", " << quadrature_.moments()[2][celli]
        << ", " << quadrature_.moments()[3][celli]
        << endl;

    celli = 10080;
    Info<< "In average Transport, after solving Mp, at celli: " << celli << endl;
    Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
        << ", " << quadrature_.moments()[1][celli]
        << ", " << quadrature_.moments()[2][celli]
        << ", " << quadrature_.moments()[3][celli]
        << endl;

    celli = 9960;
    Info<< "In average Transport, after solving Mp, at celli: " << celli << endl;
    Info<< "M0 and M1 and M2 and M3: " << quadrature_.moments()[0][celli] 
        << ", " << quadrature_.moments()[1][celli]
        << ", " << quadrature_.moments()[2][celli]
        << ", " << quadrature_.moments()[3][celli]
        << endl;
        */
    //****
    // ***
    
    if (nNodes_ == 1)
    {
        //const volScalarField& U = this->U();
        forAll(quadrature_.velocityMoments(), mi)
        {
            quadrature_.velocityMoments()[mi] = this->U()*quadrature_.moments()[mi];
            quadrature_.velocityMoments()[mi].correctBoundaryConditions();
        }

        quadrature_.updateQuadrature();
        Us_[0] = this->U();
        return;
    }

    forAll(quadrature_.velocityMoments(), mEqni)
    {
        volVectorField& Up = quadrature_.velocityMoments()[mEqni];

        volVectorField meanDivUbUp
        (
            IOobject
            (
                "meanDivUbUp",
                this->fluid().mesh().time().timeName(),
                this->fluid().mesh()
            ),
            this->fluid().mesh(),
            dimensionedVector("zero", Up.dimensions()/dimTime, Zero)
        );

        for (label nodei = 0; nodei < nNodes_; nodei++)
        {
            // Update average velocity moment flux
            surfaceVectorField aFluxUp
            (
                "aFluxUp",
                quadrature_.velocitiesNei()[nodei]
               *nodesNei[nodei].primaryWeight()
               *(
                    pow
                    (
                        nodesNei[nodei].primaryAbscissae()[0],
                        mEqni
                    )
                )*Foam::min(phi, zeroPhi)
              + quadrature_.velocitiesOwn()[nodei]
               *nodesOwn[nodei].primaryWeight()
               *pow
                (
                    nodesOwn[nodei].primaryAbscissae()[0],
                    mEqni
                )*Foam::max(phi, zeroPhi)
            );

            meanDivUbUp += fvc::surfaceIntegrate(aFluxUp);
        }

        // Solve average velocity moment transport Equation
        fvVectorMatrix UpEqn
        (
            fvm::ddt(Up)
          - fvc::ddt(Up)
          + meanDivUbUp
        );

        UpEqn.relax();
        UpEqn.solve();
    }

    //****
    forAll (quadrature_.moments()[0], celli)
    {
        if (quadrature_.moments()[0][celli] > 1e-8)
        {
            if (quadrature_.moments()[1][celli] < 0.0)
            {
                Info<< "m0>1 and m1<0, at celli:" << celli << endl;
                Info<< "M0 and M1 and M2: " << quadrature_.moments()[0][celli] 
                    << ", " << quadrature_.moments()[1][celli]
                    << ", " << quadrature_.moments()[2][celli]
                    << endl;
            }
        }
    }
    //****


    quadrature_.updateAllQuadrature();

    // Solve for velocity abscissa directly since the momentum exchange
    // terms do not change the mass
    Info << "Solving for velocity abscissae" << endl;

    forAll(Us_, nodei)
    {
        //  Colisional time, forces velocities towards mean in the case of
        //  high volume fractions
        //  Could be replaced by radial distribution function
        /*
        volScalarField tauC
        (
            "tauC",
            Foam::max
            (
                (0.5 + 0.5*tanh(((*this) - 0.63)/0.01))*GREAT,
                this->residualAlpha()
            )
        );
        tauC.dimensions().reset(dimDensity/dimTime);
        */

        // *** mass transfer
        const volScalarField& dvdtf = quadrature_.moments()[0].mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                "dvdtf",
                Foam::name(nodei)
            )
        );
        volScalarField dmdt = dvdtf * this->rho();
        volScalarField dmdt12(posPart(dmdt));
        volScalarField dmdt21(negPart(dmdt));


        volScalarField alphaRhoi(alphas_[nodei]*this->rho());
        //volScalarField& U = this->U();
        // Solve for velocities using acceleration terms
        fvVectorMatrix UsEqn
        (
            alphaRhoi*fvm::ddt(Us_[nodei])
          - alphaRhoi*fvc::ddt(Us_[nodei])
          //+ fvm::Sp(tauC, Us_[nodei])
         ==
            AEqns[nodei]
          //+ tauC*this->U()
        );

        UsEqn -= dmdt21*this->otherPhase().U() + fvm::Sp(dmdt12, Us_[nodei]);

        UsEqn.relax();
        UsEqn.solve();
        //***
        Us_[nodei].correctBoundaryConditions();
        //***
    }
    quadrature_.updateAllMoments();

    // Update moments with breakup and coalescence sources
    solveSourceOde();

    //- Update mean velocity
    this->updateVelocity();

    // Update deviation velocity
    forAll(Vs_, nodei)
    {
        Vs_[nodei] = Us_[nodei] - this->U();
    }
}


template<class BasePhaseModel>
bool Foam::myEvapULPNormDBCoBCPolydispersePhaseModel<BasePhaseModel>::read(const bool readOK)
{
    bool read = false;
    if (readOK)
    {
        maxD_.readIfPresent(this->phaseDict());
        minD_.readIfPresent(this->phaseDict());
        read = true;
    }

    if (pbeDict_.modified())
    {
        const dictionary& odeDict(pbeDict_.subDict("odeCoeffs"));
        pbeDict_.lookup("coalescence") >> coalescence_;
        pbeDict_.lookup("breakup") >> breakup_;
        odeDict.lookup("minLocalDt") >> minLocalDt_;
        odeDict.lookup("ATol") >> ATol_;
        odeDict.lookup("RTol") >> RTol_;
        odeDict.lookup("facMax") >> facMax_;
        odeDict.lookup("facMin") >> facMin_;
        odeDict.lookup("fac") >> fac_;
        read = true;
    }

    return read || readOK;
}


// ************************************************************************* //

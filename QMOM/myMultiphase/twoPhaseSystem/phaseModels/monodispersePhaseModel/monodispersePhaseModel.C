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

#include "monodispersePhaseModel.H"
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::monodispersePhaseModel<BasePhaseModel>::monodispersePhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseProperties, phaseName)
    /*
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
    */
    //solveOde_(pbeDict_.lookupOrDefault("ode", false)),
    //coalescence_(pbeDict_.lookup("coalescence")),
    //breakup_(pbeDict_.lookup("breakup")),
    //quadrature_(phaseName, fluid.mesh(), "RPlus"),
    //nNodes_(quadrature_.nodes().size()),
    //nMoments_(quadrature_.nMoments()),
    //alphas_(nNodes_),
    //Us_(quadrature_.velocities()),
    //Vs_(nNodes_),
    //ds_(nNodes_),
    //maxD_("maxD", dimLength, this->phaseDict()),
    //minD_("minD", dimLength, this->phaseDict()),
    //minLocalDt_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("minLocalDt"))),
    //localDt_(this->size(), fluid.mesh().time().deltaT().value()/10.0),
    //ATol_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("ATol"))),
    //RTol_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("RTol"))),
    //facMax_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("facMax"))),
    //facMin_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("facMin"))),
    //fac_(readScalar(pbeDict_.subDict("odeCoeffs").lookup("fac"))),
    /*
    validDirections_
    (
        (vector(this->fluid().mesh().solutionD()) + vector(1.0, 1.0, 1.0))/2.0
    )
    */
{
/*
    populationBalance_ = populationBalanceModel::New
    (
        "populationBalance", pbeDict_, this->phi()
    );
*/
    //this->d_.writeOpt() = IOobject::AUTO_WRITE;
/*
    const volVectorField& U(this->U());
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
    //correct();
*/
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::monodispersePhaseModel<BasePhaseModel>::~monodispersePhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
template<class BasePhaseModel>
void Foam::monodispersePhaseModel<BasePhaseModel>::correct()
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
                    Foam::pow
                    (
                        quadrature_.nodes()[0].primaryAbscissae()[0]*6.0
                       /(this->rho()*Foam::constant::mathematical::pi)
                      + dimensionedScalar("smallVolume", dimVolume, SMALL),
                        1.0/3.0
                    ),
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
                quadrature_.moments()[1]/this->rho(),
                this->residualAlpha()
            )
        );

        forAll(quadrature_.nodes(), nodei)
        {
            const volScalarNode& node = quadrature_.nodes()[nodei];

            // Set alpha values such that the moment.1 is equal to the bounded
            // alpha
            alphas_[nodei] =
                node.primaryWeight()*node.primaryAbscissae()[0]/this->rho()*scale;
            alphas_[nodei].max(0);
            alphas_[nodei].min(1);

            //  Calculate bubble diameter based on bubble mass (abscissa)
            ds_[nodei] =
                Foam::min
                (
                    Foam::max
                    (
                        Foam::pow
                        (
                            node.primaryAbscissae()[0]*6.0
                           /(this->rho()*Foam::constant::mathematical::pi)
                          + dimensionedScalar("smallVolume", dimVolume, SMALL),
                            1.0/3.0
                        ),
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
*/
/*
template<class BasePhaseModel>
bool Foam::monodispersePhaseModel<BasePhaseModel>::read(const bool readOK)
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
*/


// ************************************************************************* //

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

#include "pbeSMDLimitedDBasedM3PhaseModel.H"
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
Foam::pbeSMDLimitedDBasedM3PhaseModel<BasePhaseModel>::pbeSMDLimitedDBasedM3PhaseModel
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
    dsm_
    (
        IOobject
        (
            "dsm",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar(dimLength, Zero)
    ),
    maxD_("maxD", dimLength, this->phaseDict()),
    minD_("minD", dimLength, this->phaseDict())
    //maxD_("maxD", dimless, this->phaseDict()),
    //minD_("minD", dimless, this->phaseDict())
{

    populationBalance_ = populationBalanceModel::New
    (
        "populationBalance", pbeDict_, this->phi()
    );
    /*
    Info<< "dsm: " << average(dsm_) << endl;
    dsm_ = this->d_;
    Info<< "dsm: " << average(dsm_) << endl;
    */
    volScalarField& alpha = *this;
    alpha = volumeFraction();
    correctDsm();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::pbeSMDLimitedDBasedM3PhaseModel<BasePhaseModel>::~pbeSMDLimitedDBasedM3PhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BasePhaseModel>
void Foam::pbeSMDLimitedDBasedM3PhaseModel<BasePhaseModel>::correctDsm()
{
    const scalarQuadratureApproximation& quadrature(populationBalance_->quadrature());
    //const PtrList<volScalarNode>& nodes = quadrature.nodes();
    //label sizeIndex = nodes[0].sizeIndex();

    //scalar d = 0.0;

    dsm_ = this->d_;

    //const volScalarField& alpha = *this;
    dimensionedScalar Zero("Zero", dimless, 0);

    forAll(quadrature.moments()[1], celli)
    {
        //if ((*this)[celli]> this->residualAlpha())
        if (pos0((*this)[celli] - this->residualAlpha()) > Zero)
        {
            dsm_[celli] = quadrature.moments()[3][celli] / max(quadrature.moments()[2][celli], VSMALL); 
        }
    }

    // ****
    //Info<< "dsm_: " << dsm_ << endl;
    // ****

    dsm_.max(minD_);
    dsm_.min(maxD_);

    Info<< "max(dsm): " << max(dsm_) << endl
        << "min(dsm): " << min(dsm_)
        << endl;

    return;
}


template<class BasePhaseModel>
void Foam::pbeSMDLimitedDBasedM3PhaseModel<BasePhaseModel>::solve()
{
    populationBalance_->solve();

    volScalarField& alpha = *this;
    alpha = volumeFraction();

    correctDsm();
}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField> Foam::pbeSMDLimitedDBasedM3PhaseModel<BasePhaseModel>::volumeFraction() const
{
    const scalarQuadratureApproximation& quadrature(populationBalance_->quadrature());

    tmp<volScalarField> alpha
    (
        quadrature.moments()[3] * Foam::constant::mathematical::pi /6.0
    );

    return alpha;
}


// ************************************************************************* //

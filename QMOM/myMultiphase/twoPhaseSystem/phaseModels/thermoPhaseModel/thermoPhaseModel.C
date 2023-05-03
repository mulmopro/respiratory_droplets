/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "thermoPhaseModel.H"

#include "twoPhaseSystem.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class ThermoType>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::thermoPhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseProperties, phaseName),
    thermo_(ThermoType::New(fluid.mesh(), this->name()))
{
    thermo_->validate
    (
        IOobject::groupName(phaseModel::typeName, this->name()),
        "h",
        "e"
    );

    turbulence_ =
    PhaseCompressibleTurbulenceModel<BasePhaseModel>::New
    (
        *this,
        thermo_->rho(),
        this->U(),
        this->alphaRhoPhi(),
        this->phi(),
        *this
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel, class ThermoType>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::~thermoPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
template<class BasePhaseModel, class ThermoType>
Foam::PhaseCompressibleTurbulenceModel<Foam::thermoPhaseModel<BasePhaseModel, ThermoType>>&
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::turbulence()
{
    return turbulence_();
}


template<class BasePhaseModel, class ThermoType>
const Foam::PhaseCompressibleTurbulenceModel<Foam::thermoPhaseModel<BasePhaseModel, ThermoType>>&
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::turbulence() const
{
    return turbulence_();
}
*/
/*
template<class BasePhaseModel, class ThermoType>
Foam::PhaseCompressibleTurbulenceModel<Foam::phaseModel>&
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::turbulence()
{
    return turbulence_();
}


template<class BasePhaseModel, class ThermoType>
const Foam::PhaseCompressibleTurbulenceModel<Foam::phaseModel>&
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::turbulence() const
{
    return turbulence_();
}
*/

template<class BasePhaseModel, class ThermoType>
Foam::PhaseCompressibleTurbulenceModel<BasePhaseModel>&
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::turbulence()
{
    return turbulence_();
}


template<class BasePhaseModel, class ThermoType>
const Foam::PhaseCompressibleTurbulenceModel<BasePhaseModel>&
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::turbulence() const
{
    return turbulence_();
}

/*
template<class BasePhaseModel, class ThermoType>
bool Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::compressible() const
{
    return !thermo_().incompressible();
}
*/

template<class BasePhaseModel, class ThermoType>
void Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::correctEnergyTransport()
{
    BasePhaseModel::correctEnergyTransport();

    turbulence_->correctEnergyTransport();
}


template<class BasePhaseModel, class ThermoType>
const Foam::rhoThermo&
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::thermo() const
{
    return thermo_();
}


template<class BasePhaseModel, class ThermoType>
Foam::rhoThermo&
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::thermoRef()
{
    return thermo_();
}

/*
template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::rho() const
{
    return thermo_->rho();
}
*/

template<class BasePhaseModel, class ThermoType>
const Foam::volScalarField&
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::rho() const
{
    return thermo_->rho();
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::Cp() const
{
    return thermo_->Cp();
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::mu() const
{
    return thermo_->mu();
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::mu
(
    const label patchi
) const
{
    return thermo_->mu(patchi);
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::nu() const
{
    return thermo_->nu();
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::nu
(
    const label patchi
) const
{
    return thermo_->nu(patchi);
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::kappa() const
{
    return thermo_->kappa();
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::kappa
(
    const label patchi
) const
{
    return thermo_->kappa(patchi);
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::alphahe() const
{
    return thermo_->alphahe();
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::alphahe
(
    const label patchi
) const
{
    return thermo_->alphahe(patchi);
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::muEff() const
{
    return turbulence_->muEff();
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::kappaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->kappaEff(alphat);
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->kappaEff(alphat, patchi);
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::alpha() const
{
    return thermo_->alpha();
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::alpha
(
    const label patchi
) const
{
    return thermo_->alpha(patchi);
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::alphaEff
(
    const volScalarField& alphat
) const
{
    return thermo_->alphaEff(alphat);
}


template<class BasePhaseModel, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::thermoPhaseModel<BasePhaseModel, ThermoType>::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return thermo_->alphaEff(alphat, patchi);
}


// ************************************************************************* //

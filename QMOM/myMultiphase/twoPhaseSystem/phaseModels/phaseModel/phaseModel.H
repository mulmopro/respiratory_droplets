/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added additional return functions so that class can
                            be extended to polydisperse
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

Class
    Foam::phaseModel

Description
    Base class for a moving phase model. Functions are made so that the class
    can be abstracted to a polydisperse phase.

SourceFiles
    phaseModel.C
    newPhaseModel.C
    phaseModels.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "transportModel.H"
#include "rhoThermo.H"
#include "runTimeSelectionTables.H"
//#include "myMonoKineticQuadratureApproximation.H"

//#include "phaseCompressibleTurbulenceModelFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseSystem;

template<class Phase>
class PhaseCompressibleTurbulenceModel;


class phaseModel
:
    public volScalarField,
    public transportModel
{
protected:
    // Protected data

        //- Reference to mesh
        const twoPhaseSystem& fluid_;

        //- Name of phase
        word name_;

        //- Phase dictionary
        dictionary phaseDict_;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha_;

        //- Optional maximum phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        //- Thermophysical properties
        //autoPtr<rhoThermo> thermo_;

        //- Mean velocity
        volVectorField U_;

        //- Volumetric flux of the phase
        surfaceScalarField alphaPhi_;

        //- Mass flux of the phase
        surfaceScalarField alphaRhoPhi_;

        //- Volumetric flux of the phase
        autoPtr<surfaceScalarField> phiPtr_;

        //- Diameter (alpha weighted diameter for polydisperse)
        volScalarField d_;

        //- Turbulence model
        //autoPtr<PhaseCompressibleTurbulenceModel<phaseModel>> turbulence_;

        //- Bool to use B&G bubble pressure based viscosity
        bool BGviscosity_;

        //- Kinetic Energy
        mutable tmp<volScalarField> K_;


public:

    //- Runtime type information
    ClassName("phaseModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseModel,
            dictionary,
            (
                const twoPhaseSystem& fluid,
                const dictionary& dict,
                const word& phaseName
            ),
            (fluid, dict, phaseName)
        );

    // Constructors
        phaseModel
        (
            const twoPhaseSystem& fluid,
            const dictionary& dict,
            const word& phaseName
        );

         //- Return clone
        autoPtr<phaseModel> clone() const;


    // Selectors

        static autoPtr<phaseModel> New
        (
            const twoPhaseSystem& fluid,
            const dictionary& dict,
            const word& phaseName
        );


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        //- Set population balance models
        virtual void setModels()
        {}

        //- Return the name of this phase
        const word& name() const
        {
            return name_;
        }

        //- Return the dict
        const dictionary& phaseDict() const
        {
            return phaseDict_;
        }

        //- Return if B&G viscosity model is used
        bool BGviscosity() const
        {
            return BGviscosity_;
        }

        //- Return the twoPhaseSystem
        const twoPhaseSystem& fluid() const
        {
            return fluid_;
        }

        //- Return the other phase in this two-phase system
        const phaseModel& otherPhase() const;


        //- Return the number of nodes
        virtual label nNodes() const
        {
            return 1;
        }

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        const dimensionedScalar& residualAlpha() const
        {
            return residualAlpha_;
        }

        //- Constant access to alpha field for nodei
        virtual const volScalarField& alphas(const label nodei) const
        {
            return *this;
        }

        //- Non-constant access to alpha field for nodei
        virtual volScalarField& alphas(const label nodei)
        {
            return *this;
        }

        //- Return the mean diameter
        const volScalarField& d() const
        {
            return d_;
        }

        //- Return the non-const mean diameter
        volScalarField& d()
        {
            return d_;
        }

        //- Return diameter of nodei
        virtual const volScalarField& ds(const label nodei) const
        {
            return d_;
        }

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

/*
        //- Return the turbulence model
        const PhaseCompressibleTurbulenceModel<phaseModel>& turbulence() const;

        //- Return non-const access to the turbulence model
        //  for correction
        PhaseCompressibleTurbulenceModel<phaseModel>& turbulence();

        //- Return the thermophysical model
        const rhoThermo& thermo() const
        {
            return thermo_();
        }

        //- Return non-const access to the thermophysical model
        //  for correction
        rhoThermo& thermo()
        {
            return thermo_();
        }

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return thermo_->nu();
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return thermo_->nu(patchi);
        }

        //- Return the laminar dynamic viscosity
        tmp<volScalarField> mu() const
        {
            return thermo_->mu();
        }

        //- Return the laminar dynamic viscosity for patch
        tmp<scalarField> mu(const label patchi) const
        {
            return thermo_->mu(patchi);
        }

        //- Return the thermal conductivity on a patch
        tmp<scalarField> kappa(const label patchi) const
        {
            return thermo_->kappa(patchi);
        }

        //- Return the thermal conductivity
        tmp<volScalarField> kappa() const
        {
            return thermo_->kappa();
        }

        //- Return the laminar thermal conductivity
        tmp<volScalarField> kappaEff
        (
            const volScalarField& alphat
        ) const
        {
            return thermo_->kappaEff(alphat);
        }

        //- Return the laminar thermal conductivity on a patch
        tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const
        {
            return thermo_->kappaEff(alphat, patchi);
        }

        //- Return the laminar thermal diffusivity for enthalpy
        tmp<volScalarField> alpha() const
        {
            return thermo_->alpha();
        }

        //- Return the laminar thermal diffusivity for enthalpy on a patch
        tmp<scalarField> alpha(const label patchi) const
        {
            return thermo_->alpha(patchi);
        }

        //- Return the effective thermal diffusivity for enthalpy
        tmp<volScalarField> alphaEff
        (
            const volScalarField& alphat
        ) const
        {
            return thermo_->alphaEff(alphat);
        }

        //- Return the effective thermal diffusivity for enthalpy on a patch
        tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const
        {
            return thermo_->alphaEff(alphat, patchi);
        }

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        tmp<scalarField> alphahe(const label patchi) const
        {
            return thermo_->alphahe(patchi);
        }

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        tmp<volScalarField> alphahe() const
        {
            return thermo_->alphahe();
        }

        //- Return the specific heat capacity
        tmp<volScalarField> Cp() const
        {
            return thermo_->Cp();
        }

        //- Return the density
        const volScalarField& rho() const
        {
            return thermo_->rho();
        }
*/

    //- Turbulence

        //- Return the turbulence model
        virtual const PhaseCompressibleTurbulenceModel<phaseModel>& turbulence() const = 0;

        //- Return non-const access to the turbulence model
        //  for correction
        virtual PhaseCompressibleTurbulenceModel<phaseModel>& turbulence() = 0;

    //- Thermo

        //- Return the thermophysical model
        virtual const rhoThermo& thermo() const = 0;

        //- Access the thermophysical model
        virtual rhoThermo& thermoRef() = 0;

        //- Return the density field
        //virtual tmp<volScalarField> rho() const = 0;

        //- Return the density
        virtual const volScalarField& rho() const = 0;

    // Species

        //- Return whether the phase is pure (i.e., not multi-component)
        virtual bool pure() const = 0;

        //- Return the species fraction equation
        virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi) = 0;

        //- Return the species mass fractions
        virtual const PtrList<volScalarField>& Y() const = 0;

        //- Return a species mass fraction by name
        virtual const volScalarField& Y(const word& name) const = 0;

        //- Access the species mass fractions
        virtual PtrList<volScalarField>& YRef() = 0;

        //- Return the active species mass fractions
        virtual const UPtrList<volScalarField>& YActive() const = 0;

        //- Access the active species mass fractions
        virtual UPtrList<volScalarField>& YActiveRef() = 0;

    // Transport

        //- Return the laminar dynamic viscosity
        virtual tmp<volScalarField> mu() const = 0;

        //- Return the laminar dynamic viscosity on a patch
        virtual tmp<scalarField> mu(const label patchi) const = 0;

        //- Return the laminar kinematic viscosity
        virtual tmp<volScalarField> nu() const = 0;

        //- Return the laminar kinematic viscosity on a patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;

        //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
        virtual tmp<volScalarField> alpha() const = 0;

        //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alpha(const label patchi) const = 0;

        //- Thermal diffusivity for temperature of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappa() const = 0;

        //- Thermal diffusivity for temperature of mixture
        //  for patch [J/m/s/K]
        virtual tmp<scalarField> kappa(const label patchi) const = 0;

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        virtual tmp<volScalarField> alphahe() const = 0;

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alphahe(const label patchi) const = 0;

        //- Return the effective dynamic viscosity
        virtual tmp<volScalarField> muEff() const = 0;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappaEff
        (
            const volScalarField& alphat
        ) const = 0;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [J/m/s/K]
        virtual tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const = 0;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff
        (
            const volScalarField& alphat
        ) const = 0;

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const = 0;

        //- Specific heat capacity
        virtual tmp<volScalarField> Cp() const = 0;

        //- Constant access to the mean velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Non-const access to the mean velocity
        volVectorField& U()
        {
            return U_;
        }

        //- Return the velocity
        /*
        tmp<volVectorField> U() const
        {
            return U_;
        }
        */

        //- Constant access to the velocity for nodei
        virtual const volVectorField& Us(const label nodei) const
        {
            return U_;
        }

        //- Non-const access to the velocity for nodei
        virtual volVectorField& Us(const label nodei)
        {
            return U_;
        }

        //- Deviation of velocity nodei from the mean
        virtual tmp<volVectorField> Vs(const label nodei) const;

        //- Constant access to the volumetric flux
        const surfaceScalarField& phi() const
        {
            return phiPtr_();
        }

        //- Non-const access to the volumetric flux
        surfaceScalarField& phi()
        {
            return phiPtr_();
        }

        //- Return the volumetric flux of the phase
        const surfaceScalarField& alphaPhi() const
        {
            return alphaPhi_;
        }

        //- Return non-const access to the volumetric flux of the phase
        surfaceScalarField& alphaPhi()
        {
            return alphaPhi_;
        }

        //- Return the mass flux of the phase
        const surfaceScalarField& alphaRhoPhi() const
        {
            return alphaRhoPhi_;
        }

        //- Return non-const access to the mass flux of the phase
        surfaceScalarField& alphaRhoPhi()
        {
            return alphaRhoPhi_;
        }

        //- Ensure that the flux at inflow BCs is preserved
        void correctInflowOutflow(surfaceScalarField& alphaPhi) const;

        //- Correct the phase properties
        virtual void correct();

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Correct the energy transport
        virtual void correctEnergyTransport();

        //- Relative transport of moments, not used in standard phase model
        virtual void relativeTransport()
        {
            return;
        }

        //- Average transport of first size moment
        virtual void transportAlpha()
        {
            return;
        }

        //- Average transport of moments, not used in standard phase model
        virtual void averageTransport(const PtrList<fvVectorMatrix>& AEqns)
        {
            return;
        }

        //- Solve change in size moments based on breakup and coalescence
        virtual void solveBreakupCoalescence()
        {
            return;
        }

        //- Read phase properties dictionary
        virtual bool read()
        {
            return true;
        }

        //- Read phase properties dictionary
        virtual bool read(const bool readOK)
        {
            return false;
        }

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Return the phase kinetic energy
        virtual tmp<volScalarField> K() const;

        //- solve PBE
        virtual void solve() 
        {
            return;
        }

        //virtual const myMonoKineticQuadratureApproximation& quadrature() const;

        virtual const volScalarField& dotMsLim(const label nodei) const
        {
            return d_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

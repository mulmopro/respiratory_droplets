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

#include "OneResistanceHeatTransferPhaseSystem.H"
#include "fvmSup.H"
#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
OneResistanceHeatTransferPhaseSystem
(
    const fvMesh& mesh, 
    const dimensionedVector& g
)
:
    BasePhaseSystem(mesh, g)
{
    heatTransfer_.set
    (
        new BlendedInterfacialModel<heatTransferModel>
        (
            this->lookup("heatTransfer"),
            (
                this->blendingMethods().found("heatTransfer")
              ? this->blendingMethods()["heatTransfer"]
              : this->blendingMethods()["default"]
            ),
            this->pair(),
            this->pair1In2(),
            this->pair2In1()
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
~OneResistanceHeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
/*
template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr
    (
        new phaseSystem::heatTransferTable()
    );

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        eqns.set
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    // Heat transfer across the interface
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const volScalarField K(heatTransferModelIter()->K());

        const phasePair& pair(this->phasePairs_[heatTransferModelIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const volScalarField& he(phase.thermo().he());
            volScalarField Cpv(phase.thermo().Cpv());

            *eqns[phase.name()] +=
                K*(otherPhase.thermo().T() - phase.thermo().T() + he/Cpv)
              - fvm::Sp(K/Cpv, he);
        }
    }

    // Source term due to mass transfer
    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        // Note that the phase heEqn contains a continuity error term, which
        // implicitly adds a mass transfer term of fvm::Sp(dmdt, he). These
        // additions do not include this term.

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        *eqns[phase1.name()] +=
            dmdt21*he2 - fvm::Sp(dmdt21, he1) + dmdt21*(K2 - K1);

        *eqns[phase2.name()] -=
            dmdt12*he1 - fvm::Sp(dmdt12, he2) + dmdt12*(K1 - K2);
    }

    return eqnsPtr;
}
*/

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
Kh(const label nodei) const
{
    return heatTransfer_->K(nodei, 0);
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::Kh() const
{
    tmp<volScalarField> tKh
    (
        new volScalarField
        (
            IOobject
            (
                "Kh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar
            (
                "Kh",
                dimensionSet(1, -1, -3, -1, 0),
                0.0
            )
        )
    );
    for (label nodei = 0; nodei < this->phase1().nNodes(); nodei++)
    {
        tKh.ref() += Kh(nodei);
    }
    return tKh;
}



template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer1()
{
    rhoThermo& thermo1 = this->phase1().thermoRef();
    rhoThermo& thermo2 = this->phase2().thermoRef();

    volScalarField& he1 = thermo1.he();
    //volScalarField& he2 = thermo2.he();

    volScalarField Cpv1("Cpv1", thermo1.Cpv());
   // volScalarField Cpv2("Cpv2", thermo2.Cpv());

    volScalarField K = Kh();
    //*****
    //Info<< "After obtaining K1 and K2: " << endl;
    //Info<< "K average: " << average(K) << endl;
    //Info<< "heat transfer model: " << heatTransfer_ << endl; 
    //*****
    return 
    (
        K * (thermo2.T() - thermo1.T()) 
      + K * he1 / Cpv1
      - fvm::Sp(K/Cpv1, he1)
    );
}


template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer2()
{
    rhoThermo& thermo1 = this->phase1().thermoRef();
    rhoThermo& thermo2 = this->phase2().thermoRef();

    //volScalarField& he1 = thermo1.he();
    volScalarField& he2 = thermo2.he();

    //volScalarField Cpv1("Cpv1", thermo1.Cpv());
    volScalarField Cpv2("Cpv2", thermo2.Cpv());

    volScalarField K = Kh();

    return 
    (
        K * (thermo1.T() - thermo2.T()) 
      + K * he2 / Cpv2
      - fvm::Sp(K/Cpv2, he2)
    );
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::twoPhaseSystem::heatTransferTable>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    autoPtr<twoPhaseSystem::heatTransferTable> eqnsPtr
    (
        new twoPhaseSystem::heatTransferTable()
    );

    twoPhaseSystem::heatTransferTable& eqns = eqnsPtr();

    eqns.set
    (
        this->phase1().name(),
        new fvScalarMatrix(this->phase1().thermo().he(), dimEnergy/dimTime)
    );

    eqns.set
    (
        this->phase2().name(),
        new fvScalarMatrix(this->phase2().thermo().he(), dimEnergy/dimTime)
    );

    const rhoThermo& thermo1 = this->phase1().thermo();
    const rhoThermo& thermo2 = this->phase2().thermo();

    const volScalarField& he1 = thermo1.he();
    const volScalarField& he2 = thermo2.he();

    const volScalarField Cpv1("Cpv1", thermo1.Cpv());
    const volScalarField Cpv2("Cpv2", thermo2.Cpv());

    volScalarField K = Kh();

    *eqns[this->phase1().name()] += K * (thermo2.T() - thermo1.T()) 
                                + K * he1 / Cpv1
                                - fvm::Sp(K/Cpv1, he1);

    *eqns[this->phase2().name()] += K * (thermo1.T() - thermo2.T()) 
                                + K * he2 / Cpv2
                                - fvm::Sp(K/Cpv2, he2);
    
    return eqnsPtr;
}


template<class BasePhaseSystem>
bool Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "TwoResistanceHeatTransferPhaseSystem.H"

//#include "BlendedInterfacialModel.H"
#include "heatTransferModel.H"

//#include "HashPtrTable.H"

#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvMatrix.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
TwoResistanceHeatTransferPhaseSystem
(
    const fvMesh& mesh, 
    const dimensionedVector& g
)
:
    BasePhaseSystem(mesh, g),
    Tf_
    (
        IOobject
        (
            "Tf",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("0", dimTemperature, 0)
    )
{
    heatTransfer1_.set
    (
        new BlendedInterfacialModel<heatTransferModel>
        (
            this->lookup(IOobject::groupName("heatTransfer", this->phase1().name())),
            //this->lookup("heatTransfer.liquid"),
            //this->lookup("heatTransfer"),
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

    heatTransfer2_.set
    (
        new BlendedInterfacialModel<heatTransferModel>
        (
            this->lookup(IOobject::groupName("heatTransfer", this->phase2().name())),
            //this->lookup("heatTransfer.gas"),
            //this->lookup("heatTransfer"),
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

    const volScalarField& T1(this->phase1().thermo().T());
    const volScalarField& T2(this->phase2().thermo().T());

    volScalarField K1 = Kh1();
    volScalarField K2 = Kh2();
    dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);

    Tf_ = (K1*T1 + K2*T2)/max(K1 + K2, HSmall);
    /*
    this->generatePairsAndSubModels
    (
        "heatTransfer",
        heatTransferModels_,
        false
    );

    // Check that models have been specified on both sides of the interfaces
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
    
        const phasePair& pair = this->phasePairs_[heatTransferModelIter.key()];

        if (!heatTransferModels_[pair].first())
        {
            FatalErrorInFunction
                << "A heat transfer model for the " << pair.phase1().name()
                << " side of the " << pair << " pair is not specified"
                << exit(FatalError);
        }
        if (!heatTransferModels_[pair].second())
        {
            FatalErrorInFunction
                << "A heat transfer model for the " << pair.phase2().name()
                << " side of the " << pair << " pair is not specified"
                << exit(FatalError);
        }
    }
    
    // Calculate initial Tf-s as if there is no mass transfer
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[heatTransferModelIter.key()];

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& T1(phase1.thermo().T());
        const volScalarField& T2(phase2.thermo().T());

        volScalarField H1(heatTransferModels_[pair].first()->K());
        volScalarField H2(heatTransferModels_[pair].second()->K());
        dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);

        Tf_.set
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Tf", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                (H1*T1 + H2*T2)/max(H1 + H2, HSmall)
            )
        );
    }
    */
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
~TwoResistanceHeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
/*
template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
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

    // Heat transfer with the interface
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair
        (
            this->phasePairs_[heatTransferModelIter.key()]
        );

        const volScalarField& Tf(*Tf_[pair]);

        const Pair<tmp<volScalarField>> Ks
        (
            heatTransferModelIter().first()->K(),
            heatTransferModelIter().second()->K()
        );

        forAllConstIter(phasePair, pair, iter)
        {
            const phaseModel& phase = iter();

            const volScalarField& he(phase.thermo().he());
            const volScalarField Cpv(phase.thermo().Cpv());
            const volScalarField& K(Ks[iter.index()]);

            *eqns[phase.name()] +=
                K*(Tf - phase.thermo().T())
              + K/Cpv*he - fvm::Sp(K/Cpv, he);
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

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        *eqns[phase1.name()] += - fvm::Sp(dmdt21, he1) + dmdt21*(K2 - K1);

        *eqns[phase2.name()] -= - fvm::Sp(dmdt12, he2) + dmdt12*(K1 - K2);

        if (this->heatTransferModels_.found(phasePairIter.key()))
        {
            const volScalarField& Tf(*Tf_[pair]);

            *eqns[phase1.name()] +=
                dmdt21*phase1.thermo().he(phase1.thermo().p(), Tf);

            *eqns[phase2.name()] -=
                dmdt12*phase2.thermo().he(phase2.thermo().p(), Tf);
        }
        else
        {
            *eqns[phase1.name()] += dmdt21*he2;

            *eqns[phase2.name()] -= dmdt12*he1;
        }
    }

    return eqnsPtr;
}
*/

template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
correctEnergyTransport()
{
    BasePhaseSystem::correctEnergyTransport();

    correctInterfaceThermo();
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
correctInterfaceThermo()
{
    const volScalarField& T1(this->phase1().thermo().T());
    const volScalarField& T2(this->phase2().thermo().T());

    volScalarField K1 = Kh1();
    volScalarField K2 = Kh2();
    dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);

    //const volScalarField& p(this->phase1.thermo().p());

    //const volScalarField L
    //(
    //    this->phase1.thermo().he(p, Tf) - this->phase2.thermo().he(p, Tf)
    //);

    //const volScalarField dmdt(this->dmdt());

    //Tf = (H1*T1 + H2*T2 + dmdt*L)/max(H1 + H2, HSmall);
    Tf_ = (K1*T1 + K2*T2)/max(K1 + K2, HSmall);

    Info<< "Tf"
        << ": min = " << min(Tf_.primitiveField())
        << ", mean = " << average(Tf_.primitiveField())
        << ", max = " << max(Tf_.primitiveField())
        << endl;
}

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
Kh1(const label nodei) const
{
    return heatTransfer1_->K(nodei, 0);
}

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
Kh2(const label nodei) const
{
    return heatTransfer2_->K(nodei, 0);
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::Kh1() const
{
    tmp<volScalarField> tKh
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("Kh", this->phase1().name()),
                //"Kh1",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar
            (
                IOobject::groupName("Kh", this->phase1().name()),
                dimensionSet(1, -1, -3, -1, 0),
                0.0
            )
        )
    );
    for (label nodei = 0; nodei < this->phase1().nNodes(); nodei++)
    {
        tKh.ref() += Kh1(nodei);
    }
    return tKh;
}

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::Kh2() const
{
    tmp<volScalarField> tKh
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("Kh", this->phase2().name()),
                //"Kh2",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar
            (
                IOobject::groupName("Kh", this->phase2().name()),
                dimensionSet(1, -1, -3, -1, 0),
                0.0
            )
        )
    );
    for (label nodei = 0; nodei < this->phase1().nNodes(); nodei++)
    {
        tKh.ref() += Kh2(nodei);
    }
    return tKh;
}



template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer1()
{
    rhoThermo& thermo1 = this->phase1().thermoRef();
    rhoThermo& thermo2 = this->phase2().thermoRef();

    volScalarField& he1 = thermo1.he();
    //volScalarField& he2 = thermo2.he();

    volScalarField Cpv1("Cpv1", thermo1.Cpv());
   // volScalarField Cpv2("Cpv2", thermo2.Cpv());

    //*****
    //Info<< "Before obtaining K1 and K2: " << endl;
    //*****

    volScalarField K1 = Kh1();
    volScalarField K2 = Kh2();
    dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);

    //*****
    //Info<< "After obtaining K1 and K2: " << endl;
    //Info<< "K1 average: " << average(K1) << endl
    //    << "K2 average: " << average(K2) << endl; 
    //Info<< "heat transfer model 1: " << heatTransfer1_ << endl 
    //    << "heat transfer model 2: " << heatTransfer2_ << endl; 
    //Info<< IOobject::groupName("heatTransfer", this->phase1().name()) << endl;
    //Info<< IOobject::groupName("heatTransfer", this->phase2().name()) << endl;
    //*****

    return 
    (
        (K1 * K2)/max(K1 + K2, HSmall) * (thermo2.T() - thermo1.T()) 
      + K1 * he1 / Cpv1
      - fvm::Sp(K1/Cpv1, he1)
    );
}


template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer2()
{
    rhoThermo& thermo1 = this->phase1().thermoRef();
    rhoThermo& thermo2 = this->phase2().thermoRef();

    //volScalarField& he1 = thermo1.he();
    volScalarField& he2 = thermo2.he();

    //volScalarField Cpv1("Cpv1", thermo1.Cpv());
    volScalarField Cpv2("Cpv2", thermo2.Cpv());

    volScalarField K1 = Kh1();
    volScalarField K2 = Kh2();
    dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);

    return 
    (
        (K1 * K2)/max(K1 + K2, HSmall) * (thermo1.T() - thermo2.T()) 
      + K2 * he2 / Cpv2
      - fvm::Sp(K2/Cpv2, he2)
    );
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::twoPhaseSystem::heatTransferTable>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
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

    volScalarField H1 = Kh1();
    volScalarField H2 = Kh2();
    //dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);
    //volScalarField HEff = (H1 * H2)/max(H1 + H2, HSmall);

    /*
    *eqns[this->phase1().name()] += HEff * (thermo2.T() - thermo1.T()) 
                                + H1 * he1 / Cpv1
                                - fvm::Sp(H1/Cpv1, he1);

    *eqns[this->phase2().name()] += HEff * (thermo1.T() - thermo2.T()) 
                                + H2 * he2 / Cpv2
                                - fvm::Sp(H2/Cpv2, he2);
    */
    *eqns[this->phase1().name()] += H1 * (Tf_ - thermo1.T()) 
                                + H1 * he1 / Cpv1
                                - fvm::Sp(H1/Cpv1, he1);

    *eqns[this->phase2().name()] += H2 * (Tf_ - thermo2.T()) 
                                + H2 * he2 / Cpv2
                                - fvm::Sp(H2/Cpv2, he2);
    
    return eqnsPtr;
}


template<class BasePhaseSystem>
bool Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::read()
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

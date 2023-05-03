/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "QMOMInterfaceCompositionPhaseChangePhaseSystem.H"
//#include "interfaceCompositionModel.H"
//#include "massTransferModel.H"
#include "../interfacialCompositionModels/interfaceCompositionModels/interfaceCompositionModel_1/interfaceCompositionModel.H"
#include "../interfacialCompositionModels/massTransferModels/massTransferModel/massTransferModel.H"
#include "fvmSup.H"
#include "fvc.H"
#include "fvm.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::iDmdt
() const
{
    tmp<volScalarField> tIDmdt = twoPhaseSystem::dmdt();

    //const phasePair& pair = this->pair();

    for
    (
        const word& member
        : interfaceCompositionModel2_->species()
    )
    {
        tIDmdt.ref() +=
            (
                //*(*iDmdtSu_[pair])[member]
                //+ *(*iDmdtSp_[pair])[member]*this->phase2().Y(member)

                *iDmdtSu_[member]
                + *iDmdtSp_[member]*this->phase2().Y(member)
            );
    }

    return tIDmdt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
QMOMInterfaceCompositionPhaseChangePhaseSystem
(
    const fvMesh& mesh, 
    const dimensionedVector& g
)
:
    BasePhaseSystem(mesh, g),
    nInterfaceCorrectors_
    (
        this->template getOrDefault<label>("nInterfaceCorrectors", 1)
    ),
    dvdtf_(this->phase1().nNodes())
{
    massTransferModel2_.set
    (
        new BlendedInterfacialModel<massTransferModel>
        (
            this->lookup(IOobject::groupName("massTransfer", this->phase2().name())),
            //this->lookup("heatTransfer.liquid"),
            //this->lookup("heatTransfer"),
            (
                this->blendingMethods().found("massTransfer")
              ? this->blendingMethods()["massTransfer"]
              : this->blendingMethods()["default"]
            ),
            this->pair(),
            this->pair1In2(),
            this->pair2In1()
        )
    );
    /*
    Info<< *this << endl;
    Info<< this->lookup(IOobject::groupName("interfaceComposition", this->phase2().name())) << endl;
    //Info<< this->lookup(IOobject::groupName("interfaceComposition", this->phase2().name()))[this->pair()] << endl;
    Info<< this->lookup(IOobject::groupName("massTransfer", this->phase2().name())) << endl;
    */
    //Info<< "before defining iCDictTable as dictTable" << endl;
    const phasePair::dictTable& iCDictTable(this->lookup(IOobject::groupName("interfaceComposition", this->phase2().name())));
    //Info<< "after defining iCDictTable" << endl;
    //Info<< "Dictionary: " << iCDictTable << endl;
    /*
    if (iCDictTable.found(this->pair2In1()))
    {
        Info<< "Dictionary for 2 in 1: " << iCDictTable[this->pair2In1()] << endl;
    }
    if (iCDictTable.found(this->pair1In2()))
    {
        Info<< "Dictionary for 1 in 2: " << iCDictTable[this->pair1In2()] << endl;
    }
    */
    Info<< this->pair() << endl;
    Info<< this->pair().phase1().name() << endl
        << this->pair().phase1().name() << endl;
/*
    Info<< "before defining iCDict as dictionary" << endl;
    const dictionary& iCDict(this->lookup(IOobject::groupName("interfaceComposition", this->phase2().name())));
    Info<< "after defining iCDict" << endl;
    Info<< "Dictionary: " << iCDict << endl;
*/
    /*
    if (iCDict.found(this->pair2In1()))
    {
        Info<< "Dictionary for 2in 1: " << iCDict[this->pair2In1()] << endl;
    }
    */

    interfaceCompositionModel2_.set
    (
        interfaceCompositionModel::New
        (
            //this->lookup(IOobject::groupName("interfaceComposition", this->phase2().name())),
            //this->subDict(IOobject::groupName("interfaceComposition", this->phase2().name())),
            //this->found(IOobject::groupName("interfaceComposition", this->phase2().name())),
            //this->pair()
            iCDictTable[this->pair1In2()],
            this->pair1In2()
        ).ptr()
    );

    // Generate mass transfer fields, initially assumed to be zero
    for (const word& member : interfaceCompositionModel2_->species())
    {
        //iDmdtSu_[pair]->set
        iDmdtSu_.set
        (
            member,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("iDmdtSu", this->pair().name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime)
            )
        );

        //iDmdtSp_[pair]->set
        iDmdtSp_.set
        (
            member,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("iDmdtSp", this->pair().name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime)
            )
        );
    }

    forAll(dvdtf_, nodei)
    {
        dvdtf_.set
        (
            nodei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "dvdtf",
                        Foam::name(nodei)
                        /*
                        IOobject::groupName
                        (
                            this->pair().name(),
                            Foam::name(nodei)
                        )
                        */
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimless/dimTime)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
~QMOMInterfaceCompositionPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdt
() const
{
    //Info<< "interfaceComposition: dmdt!" << endl;
    return BasePhaseSystem::dmdt() + this->iDmdt();
    //this->iDmdt();
    //return BasePhaseSystem::dmdt();
}

/*
template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const interfaceCompositionModel& compositionModel =
            interfaceCompositionModelIter();

        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        for (const word& member : compositionModel.species())
        {
            const volScalarField iDmdt
            (
                *(*iDmdtSu_[pair])[member]
              + *(*iDmdtSp_[pair])[member]*phase.Y(member)
            );

            this->addField(phase, "dmdt", iDmdt, dmdts);
            this->addField(otherPhase, "dmdt", - iDmdt, dmdts);
        }
    }

    return dmdts;
}
*/

template<class BasePhaseSystem>
Foam::autoPtr<Foam::twoPhaseSystem::massTransferTable>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
massTransfer() const
{
    autoPtr<twoPhaseSystem::massTransferTable> eqnsPtr =
        BasePhaseSystem::massTransfer();

    twoPhaseSystem::massTransferTable& eqns = eqnsPtr();

    for
    (
        const word& member
        : interfaceCompositionModel2_->species()
    )
    {
        const word name(IOobject::groupName(member, this->phase2().name()));

        *eqns[name] = *iDmdtSu_[member]
                     + fvm::Sp(*iDmdtSp_[member], this->phase2().Y(member));
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correct()
{
    BasePhaseSystem::correct();

    const volScalarField& Tf = this->Tf_;

    //const phasePair& pair = this->pair();

    for
    (
        const word& member
        : interfaceCompositionModel2_->species()
    )
    {
        const volScalarField KD(Km2()*interfaceCompositionModel2_->D(member));

        const volScalarField Yf(interfaceCompositionModel2_->Yf(member, Tf));

        //*(*iDmdtSu_[pair])[member] = this->phase1().rho()*KD*Yf;
        //*(*iDmdtSp_[pair])[member] = - this->phase1().rho()*KD;
        *iDmdtSu_[member] = this->phase2().rho()*KD*Yf;
        *iDmdtSp_[member] = - this->phase2().rho()*KD;
    }
}

template<class BasePhaseSystem>
void Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correctInterfaceThermo()
{
    // This loop solves for the interface temperatures, Tf, and updates the
    // interface composition models.
    //
    // The rate of heat transfer to the interface must equal the latent heat
    // consumed at the interface, i.e.:
    //
    // H1*(T1 - Tf) + H2*(T2 - Tf) == mDotL
    //                             == K*rho*(Yfi - Yi)*Li
    //
    // Yfi is likely to be a strong non-linear (typically exponential) function
    // of Tf, so the solution for the temperature is newton-accelerated

    const volScalarField H1(this->Kh1());
    const volScalarField H2(this->Kh2());
    const dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);

    volScalarField& Tf = this->Tf_;

    for (label i = 0; i < nInterfaceCorrectors_; ++ i)
    {
        volScalarField mDotL
        (
            IOobject
            (
                "mDotL",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimVolume/dimTime)
        );
        volScalarField mDotLPrime
        (
            IOobject
            (
                "mDotLPrime",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar(mDotL.dimensions()/dimTemperature)
        );

        // Add latent heats from forward and backward models

        interfaceCompositionModel2_->addMDotL
        (
            Km2(),
            Tf,
            mDotL,
            mDotLPrime
        );

        // Update the interface temperature by applying one step of newton's
        // method to the interface relation
        Tf -=
            (
                H1*(Tf - this->pair().phase1().thermo().T())
                + H2*(Tf - this->pair().phase2().thermo().T())
                + mDotL
            )
            /(
                max(H1 + H2 + mDotLPrime, HSmall)
            );

        Tf.correctBoundaryConditions();

        // Update the interface compositions
        interfaceCompositionModel2_->update(Tf);
    }

    forAll(dvdtf_, nodei)
    {
        //volScalarField tdmdt = twoPhaseSystem::dmdt();
        tmp<volScalarField>tdmdt = twoPhaseSystem::dmdt();
        for
        (
            const word& member
            : interfaceCompositionModel2_->species()
        )
        {
            const volScalarField KD(Km2(nodei)*interfaceCompositionModel2_->D(member));

            const volScalarField Yf(interfaceCompositionModel2_->Yf(member, Tf));

            //*(*iDmdtSu_[pair])[member] = this->phase1().rho()*KD*Yf;
            //*(*iDmdtSp_[pair])[member] = - this->phase1().rho()*KD;
            //*iDmdtSu_[member] = this->phase1().rho()*KD*Yf;
            //*iDmdtSp_[member] = - this->phase1().rho()*KD;
            //tdmdt += this->phase2().rho()*KD*(Yf - this->phase2().Y(member));
            tdmdt.ref() += this->phase2().rho()*KD*(Yf - this->phase2().Y(member));
        }

        dvdtf_[nodei] = tdmdt.ref() / this->phase1().rho();
        Info<< "Nodei: " << nodei << endl;
    }

/*
    Info<< "dvdtf" 
        << ": min = " << min(dvdtf_.primitiveField())
        << ", mean = " << average(dvdtf_.primitiveField())
        << ", max = " << max(dvdtf_.primitiveField())
        << endl;
*/
}


template<class BasePhaseSystem>
bool Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::read()
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


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
Km2(const label nodei) const
{
    return massTransferModel2_->K(nodei, 0);
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::Km2() const
{
    tmp<volScalarField> tKm
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("Km", this->phase2().name()),
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
                IOobject::groupName("Km", this->phase2().name()),
                dimensionSet(0, -2, 0, 0, 0),
                0.0
            )
        )
    );
    for (label nodei = 0; nodei < this->phase1().nNodes(); nodei++)
    {
        tKm.ref() += Km2(nodei);
    }
    return tKm;
}

/*
template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
heatTransfer1()
{
    //tmp<fvScalarMatrix> hT = BasePhaseSystem::heatTransfer1();
    
    rhoThermo& thermo1 = this->phase1().thermoRef();
    rhoThermo& thermo2 = this->phase2().thermoRef();

    volScalarField& he1 = thermo1.he();
    //volScalarField& he2 = thermo2.he();

    volScalarField Cpv1("Cpv1", thermo1.Cpv());
   // volScalarField Cpv2("Cpv2", thermo2.Cpv());

    // *****
    //Info<< "Before obtaining K1 and K2: " << endl;
    // *****

    volScalarField H1 = this->Kh1();
    volScalarField H2 = this->Kh2();
    volScalarField H1Fac(H1/(H1 + H2));

    dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);

    const volScalarField& Tf(this->Tf_);
    const volScalarField hef1(thermo1.he(thermo1.p(), Tf));
    const volScalarField hef2(thermo2.he(thermo2.p(), Tf));
    const volScalarField hc1(thermo1.hc());
    const volScalarField hc2(thermo2.hc());

    // Kinetic energy
    const volScalarField K1(this->phase1().K());
    const volScalarField K2(this->phase2().K());

    for
    (
        const word& member
        : interfaceCompositionModel2_->species()
    )
    {
        
    }

    return 
    (
        (H1 * H2)/max(H1 + H2, HSmall) * (thermo2.T() - thermo1.T()) 
      + H1 * he1 / Cpv1
      - fvm::Sp(H1/Cpv1, he1)
    );
}


template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
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
*/

template<class BasePhaseSystem>
Foam::autoPtr<Foam::twoPhaseSystem::heatTransferTable>
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    //Info<< "heatTransfer()" << endl;

    const phaseModel& phase1 = this->phase1();
    const phaseModel& phase2 = this->phase2();

    autoPtr<twoPhaseSystem::heatTransferTable> eqnsPtr
    (
        new twoPhaseSystem::heatTransferTable()
    );

    twoPhaseSystem::heatTransferTable& eqns = eqnsPtr();

    //Info<< "heatTransfer(): set eqns" << endl;

    eqns.set
    (
        phase1.name(),
        new fvScalarMatrix(phase1.thermo().he(), dimEnergy/dimTime)
    );

    eqns.set
    (
        phase2.name(),
        new fvScalarMatrix(phase2.thermo().he(), dimEnergy/dimTime)
    );

    //Info<< "heatTransfer(): calculating" << endl;

    const rhoThermo& thermo1 = phase1.thermo();
    const rhoThermo& thermo2 = phase2.thermo();

    //Info<< "heatTransfer(): get he" << endl;

    const volScalarField& he1 = thermo1.he();
    const volScalarField& he2 = thermo2.he();
    const volScalarField K1(phase1.K());
    const volScalarField K2(phase2.K());

    //Info<< "heatTransfer(): get Cpv" << endl;

    const volScalarField Cpv1("Cpv1", thermo1.Cpv());
    const volScalarField Cpv2("Cpv2", thermo2.Cpv());

    //Info<< "heatTransfer(): get H1, H2" << endl;

    volScalarField H1 = this->Kh1();
    volScalarField H2 = this->Kh2();
    //dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);
    //volScalarField HEff = (H1 * H2)/max(H1 + H2, HSmall);
    //volScalarField H1Fac(H1/max(H1 + H2, HSmall));

    //Info<< "heatTransfer(): get Tf" << endl;

    const volScalarField& Tf(this->Tf_);
    //Info<< "heatTransfer(): get hef1" << endl;
    const volScalarField hef1(thermo1.he(thermo1.p(), Tf));
    //Info<< "heatTransfer(): get hef2" << endl;
    //const volScalarField hef2(thermo2.he(thermo2.p(), Tf));
    /*
    Info<< "heatTransfer(): get hef" << endl;
    const volScalarField hef1(thermo1.he(thermo1.p(), Tf));
    const volScalarField hef2(thermo2.he(thermo2.p(), Tf));
    Info<< "heatTransfer(): get hc" << endl;
    const volScalarField hc1(thermo1.hc());
    const volScalarField hc2(thermo2.hc());

    Info<< "Cpv1: " << average(Cpv1) << endl;
    Info<< "Cpv2: " << average(Cpv2) << endl;
    
    *eqns[phase1.name()] += HEff * (thermo2.T() - thermo1.T()) 
                                + H1 * he1 / Cpv1
                                - fvm::Sp(H1/Cpv1, he1);

    *eqns[phase2.name()] += HEff * (thermo1.T() - thermo2.T()) 
                                + H2 * he2 / Cpv2
                                - fvm::Sp(H2/Cpv2, he2);
    */
    *eqns[this->phase1().name()] += H1 * (Tf - thermo1.T()) 
                                + H1 * he1 / Cpv1
                                - fvm::Sp(H1/Cpv1, he1);

    *eqns[this->phase2().name()] += H2 * (Tf - thermo2.T()) 
                                + H2 * he2 / Cpv2
                                - fvm::Sp(H2/Cpv2, he2);

    //Info<< "heatTransfer(): mass transfer" << endl;

    for (const word& member : interfaceCompositionModel2_->species())
    {
        volScalarField dmidtf
        (
            *iDmdtSu_[member]
          + *iDmdtSp_[member]*phase2.Y(member)
        );

        const volScalarField dmidtf21(negPart(dmidtf));
        const volScalarField dmidtf12(posPart(dmidtf));
/*
        volScalarField tdmdt = twoPhaseSystem::dmdt();
        const volScalarField KD(Km2()*interfaceCompositionModel2_->D(member));
        const volScalarField Yf(interfaceCompositionModel2_->Yf(member, Tf));
        tdmdt = this->phase2().rho()*KD*(Yf - this->phase2().Y(member));

        Info<< "dmidtf: " << average(dmidtf) << endl;
        Info<< "tdmdt: "  << average(tdmdt) << endl;
        Info<< "dmidtf21: " << average(dmidtf21) << endl
            << "dmidtf12: " << average(dmidtf12) << endl;
*/
        //volScalarField hefi1(hef1), hai1(hc1);

        //Info<< "Before thermo 1" << endl;
        volScalarField hefi1(hef1);
        if (isA<rhoReactionThermo>(thermo1))
        {
            //Info<< "Thermo 1 is a rhoReactionThermo" << endl;
            forAll(thermo1.p(), celli)
            {
                const basicSpecieMixture& composition1 =
                    refCast<const rhoReactionThermo>(thermo1).composition();
                hefi1[celli] =
                    composition1.HE
                    (
                        composition1.species()[member],
                        thermo1.p()[celli],
                        Tf[celli]
                    );
                /*
                hai1[celli] =
                    composition1.Ha
                    (
                        composition1.species()[member],
                        thermo1.p()[celli],
                        Tf[celli]
                    );
                */
            }
        }

        volScalarField hefi2(hef1);
        if (isA<rhoReactionThermo>(thermo2))
        {
            //Info<< "Thermo 2 is a rhoReactionThermo" << endl;
            forAll(thermo2.p(), celli)
            {
                const basicSpecieMixture& composition2 =
                    refCast<const rhoReactionThermo>(thermo2).composition();

                //Info<< "Before thermo 2 hf" << endl;

                hefi2[celli] =
                    composition2.HE
                    (
                        composition2.species()[member],
                        thermo2.p()[celli],
                        Tf[celli]
                    );

                //Info<< "Before thermo 2 ha" << endl;
                /*
                hai2[celli] =
                    composition2.Ha
                    (
                        composition2.species()[member],
                        thermo2.p()[celli],
                        Tf[celli]
                    );
                */
            }
        }
        /*
        // Create the latent heat for the transferring specie
        const volScalarField Li(hai2 - hai1);
        const volScalarField L(interfaceCompositionModel2_->L(member, Tf));

        Info<< "Li: " << average(Li) << endl
            << "L: " << average(L) << endl
            << "H1Fac: " << average(H1Fac) << endl
            << "dmidtf: " << average(dmidtf) << endl
            << "hefi1: " << average(hefi1) << endl
            << "psi.phase1: " << average(eqns[phase1.name()]->psi()) << endl;

        Info<< "enqs.phase1 source: " << average(eqns[phase1.name()]->source()) << endl;
        
        // Transfer of energy from the interface into the bulk
        *eqns[phase1.name()] += dmidtf*hefi1;
        *eqns[phase2.name()] -= dmidtf*hefi2;

        Info<< "enqs.phase1 source: " << average(eqns[phase1.name()]->source()) << endl;

        // Latent heat contribution
        

        Info<< "enqs.phase1 source: " << average(eqns[phase1.name()]->source()) << endl;
        */
        // Transfer of energy from the interface into the bulk
        //Info<< "average(hefi1): " << average(hefi1) << endl
        //    << "average(hefi2): " << average(hefi2) << endl;
/*
        *eqns[phase1.name()] += dmidtf*hefi1;
        *eqns[phase2.name()] -= dmidtf*hefi2;

        // Transfer of kinetic energy
        *eqns[phase1.name()] += dmidtf21*K2 + dmidtf12*K1;
        *eqns[phase2.name()] -= dmidtf12*K1 + dmidtf21*K2;
*/
        *eqns[phase1.name()] -= dmidtf*hefi1;
        *eqns[phase2.name()] += dmidtf*hefi2;

        // Transfer of kinetic energy
        *eqns[phase1.name()] -= dmidtf21*K2 + dmidtf12*K1;
        *eqns[phase2.name()] += dmidtf12*K1 + dmidtf21*K2;
    }
/*
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
    const volScalarField K1(this->phase1().K());
    const volScalarField K2(this->phase2().K());

    const volScalarField Cpv1("Cpv1", thermo1.Cpv());
    const volScalarField Cpv2("Cpv2", thermo2.Cpv());

    volScalarField H1 = Kh1();
    volScalarField H2 = Kh2();
    dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);
    volScalarField HEff = (H1 * H2)/max(H1 + H2, HSmall);
    volScalarField H1Fac(H1/max(H1 + H2, HSmall));

    const volScalarField& Tf(this->Tf_]);
    const volScalarField hef1(thermo1.he(thermo1.p(), Tf));
    const volScalarField hef2(thermo2.he(thermo2.p(), Tf));
    const volScalarField hc1(thermo1.hc());
    const volScalarField hc2(thermo2.hc());

    *eqns[this->phase1().name()] += HEff * (thermo2.T() - thermo1.T()) 
                                + H1 * he1 / Cpv1
                                - fvm::Sp(H1/Cpv1, he1);

    *eqns[this->phase2().name()] += HEff * (thermo1.T() - thermo2.T()) 
                                + H2 * he2 / Cpv2
                                - fvm::Sp(H2/Cpv2, he2);

    for (const word& member : interfaceCompositionModel2_->species())
    {
        tmp<volScalarField> dmidtf
        (
            *iDmdtSu_[member]
          + *iDmdtSp_[member]*this->phase2().Y(member)
        );

        const volScalarField dmidtf21(negPart(dmidtf));
        const volScalarField dmidtf12(posPart(dmidtf));

        volScalarField hefi1(hef1), hci1(hc1);

        if (isA<rhoReactionThermo>(thermo1))
        {
            const basicSpecieMixture& composition1 =
                refCast<const rhoReactionThermo>(thermo1).composition();
            hefi1 =
                composition1.HE
                (
                    composition1.species()[member],
                    thermo1.p(),
                    Tf
                );
            hci1 =
                dimensionedScalar
                (
                    dimEnergy/dimMass,
                    composition1.Hf(composition1.species()[member])
                );
        }
        volScalarField hefi2(hef2), hci2(hc2);
        if (isA<rhoReactionThermo>(thermo2))
        {
            const basicSpecieMixture& composition2 =
                refCast<const rhoReactionThermo>(thermo2).composition();
            hefi2 =
                composition2.HE
                (
                    composition2.species()[member],
                    thermo2.p(),
                    Tf
                );
            hci2 =
                dimensionedScalar
                (
                    dimEnergy/dimMass,
                    composition2.Hf(composition2.species()[member])
                );
        }

        // Create the latent heat for the transferring specie
        const volScalarField Li(hefi2 + hci2 - hefi1 - hci1);

        // Transfer of energy from the interface into the bulk
        *eqns[this->phase1().name()] += dmidtf*hefi1;
        *eqns[this->phase2().name()] -= dmidtf*hefi2;

        // Latent heat contribution
        *eqns[this->phase1().name()] += H1Fac*dmidtf*Li;
        *eqns[this->phase2().name()] += (1 - H1Fac)*dmidtf*Li;

        // Transfer of kinetic energy
        *eqns[this->phase1().name()] += dmidtf21*K2 + dmidtf12*K1;
        *eqns[this->phase2().name()] -= dmidtf12*K1 + dmidtf21*K2;
    }
*/ 
    return eqnsPtr;
}
/*
template<class BasePhaseSystem>
Foam::HashPtrTable<Foam::volScalarField> Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmidtfs() const
{
    HashPtrTable<volScalarField> dmidtfsPtr(new HashPtrTable<volScalarField>);
    HashPtrTable<volScalarField>& dmidtfs = dmidtfsPtr();

    for (const word& member : interfaceCompositionModel2_->species())
    {
        tmp<volScalarField> dmidtf
        (
            *iDmdtSu_[member]
          + *iDmdtSp_[member]*this->phase2().Y(member)
        );

        if (dmidtfs->found(member))
        {
            *dmidtfs[member] += dmidtf;
        }
        else
        {
            dmidtfs->insert(member, dmidtf.ptr());
        }
    }
}
*/

/*
template<class BasePhaseSystem>
void Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::averageTransport()
{
    const phaseModel& phase1 = this->phase1();
    const phaseModel& phase2 = this->phase2();

    PtrList<fvVectorMatrix> AEqns(phase1.nNodes());

    if (phase1.nNodes() == 1)
    {
        this->phase1().averageTransport(AEqns);
        this->phase1().correct();

        return;
    }

    // Liquid viscous stress
    volSymmTensorField taul(phase2.turbulence().devRhoReff());

    // Acceleration of liquid phase
    volVectorField DDtU2
    (
        fvc::ddt(phase2.U())
      + fvc::div(phase2.phi(), phase2.U())
      - fvc::div(phase2.phi())*phase2.U()
    );

    for (label nodei = 0; nodei < phase1.nNodes(); nodei++)
    {
        //  Build matrix to solve for velocity abscissae due to interfacial
        //  forces
        AEqns.set
        (
            nodei,
            new fvVectorMatrix
            (
                phase1.Us(nodei),
                phase1.Us(nodei).dimensions()*dimDensity*dimVol/dimTime
            )
        );
        // *****
        Info<< "Define AEqns over" << endl;
        // *****
        const volScalarField& p(mesh_.lookupObject<volScalarField>("p"));
        volScalarField alphaRhoi(phase1.alphas(nodei)*phase1.rho());

        //  Implicit drag term added to velocity abscissae equations
        volScalarField Kd(this->Kd(nodei));

        // Interfacial forces
        AEqns[nodei] +=
            // Buoyancy
            g_*alphaRhoi
          + (
              - fvc::grad(p)
              + fvc::div(taul)
            )*phase1.alphas(nodei)

            // Drag
          + Kd*phase2.U()
          - fvm::Sp(Kd, phase1.Us(nodei))

            // Virtual Mass
          + Vm(nodei)
           *(
                DDtU2
              - (
                    fvm::ddt(phase1.Us(nodei))
                  + fvm::div(phase1.phi(), phase1.Us(nodei))
                  - fvm::Sp(fvc::div(phase1.phi()), phase1.Us(nodei))
                )
            )

            // Dispersion, lift, wall lubrication, and bubble pressure
          - turbulentDispersion_->F<vector>(nodei, 0)
          - lift_->F<vector>(nodei, 0)
          - wallLubrication_->F<vector>(nodei, 0)
          + bubblePressure_->F<vector>(nodei, 0);


        //- mass transfer

        const volScalarField& Tf = this->Tf_;

        for
        (
            const word& member
            : interfaceCompositionModel2_->species()
        )
        {
            const volScalarField KD(Km2(nodei)*interfaceCompositionModel2_->D(member));

            const volScalarField Yf(interfaceCompositionModel2_->Yf(member, Tf));

            // *(*iDmdtSu_[pair])[member] = this->phase1().rho()*KD*Yf;
            // *(*iDmdtSp_[pair])[member] = - this->phase1().rho()*KD;
            // *iDmdtSu_[member] = this->phase1().rho()*KD*Yf;
            // *iDmdtSp_[member] = - this->phase1().rho()*KD;

            volScalarField dmidtf
            (
                this->phase1().rho()*KD*(Yf - this->phase2().Y(member))
            );

            const volScalarField dmidtf21(negPart(dmidtf));
            const volScalarField dmidtf12(posPart(dmidtf));

            AEqns[nodei] += dmidtf21*phase2.U() + fvm::Sp(dmidtf12, phase1.Us(nodei));
        }

    }

    phase1.averageTransport(AEqns);
    phase1.correct();

    phi_ = phase1.alphaPhi() + phase2.alphaPhi();
}
*/

/*
template<class BasePhaseSystem>
Foam::PtrList<Foam::fvVectorMatrix> 
Foam::QMOMInterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::getAEqns() const
{
    // ***
    Info<< "QMOMInterComp getEqns" << endl;
    // ***
    PtrList<fvVectorMatrix> AEqns(twoPhaseSystem::getAEqns());

    if (this->phase1().nNodes() == 1)
    {
        return AEqns;
    }

    for (label nodei = 0; nodei < this->phase1().nNodes(); nodei++)
    {
        const volScalarField& Tf = this->Tf_;
        tmp<volScalarField> dmdti = twoPhaseSystem::dmdt();

        for
        (
            const word& member
            : interfaceCompositionModel2_->species()
        )
        {
            const volScalarField KD(Km2(nodei)*interfaceCompositionModel2_->D(member));

            const volScalarField Yf(interfaceCompositionModel2_->Yf(member, Tf));

            // *(*iDmdtSu_[pair])[member] = this->phase1().rho()*KD*Yf;
            // *(*iDmdtSp_[pair])[member] = - this->phase1().rho()*KD;
            // *iDmdtSu_[member] = this->phase1().rho()*KD*Yf;
            // *iDmdtSp_[member] = - this->phase1().rho()*KD;
            dmdti.ref() += this->phase2().rho()*KD*(Yf - this->phase2().Y(member));
        }

        const volScalarField dmdti21 = posPart(dmdti);
        const volScalarField dmdti12 = negPart(dmdti);

        //- momentum
        AEqns[nodei] += dmdti21 * this->phase2().U() + fvm::Sp(dmdti12, this->phase1().Us(nodei));
    }

    return AEqns;
}
*/
// ************************************************************************* //

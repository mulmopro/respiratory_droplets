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

Class
    Foam::InterfaceCompositionPhaseChangePhaseSystem

Description
    Class to provide interfacial heat and mass transfer between a number of
    phases according to a interface composition model.

    The interface temperature is calculated such that the net rate at which the
    heat is transferred to the interface is equal to the latent heat consumed by
    the mass transfer.

SourceFiles
    InterfaceCompositionPhaseChangePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef InterfaceCompositionPhaseChangePhaseSystem_H
#define InterfaceCompositionPhaseChangePhaseSystem_H

#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class interfaceCompositionModel;
class massTransferModel;

/*---------------------------------------------------------------------------*\
                 Class InterfaceCompositionPhaseChangePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class InterfaceCompositionPhaseChangePhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected Typedefs
/*
        typedef HashTable
        <
            autoPtr<interfaceCompositionModel>,
            phasePairKey,
            phasePairKey::hash
        > interfaceCompositionModelTable;

        typedef HashTable
        <
            Pair<autoPtr<BlendedInterfacialModel<massTransferModel>>>,
            phasePairKey,
            phasePairKey::hash
        > massTransferModelTable;

        typedef HashPtrTable
        <
            HashPtrTable<volScalarField>,
            phasePairKey,
            phasePairKey::hash
        > iDmdtSuSpTable;

        typedef HashPtrTable
        <
            volScalarField,
            phasePairKey,
            phasePairKey::hash
        > iDmdtTable;
*/

        //typedef HashPtrTable
        //<
        //    HashPtrTable<volScalarField>,
        //    phasePairKey,
        //    phasePairKey::hash
        //> iDmdtSuSpTable;

    // Protected Data

        // Sub Models

            //- The number of interface correctors
            const label nInterfaceCorrectors_;

            //- used for the drift term in PBE
            volScalarField dvdtf_;

            //- Mass transfer models
            //massTransferModelTable massTransferModels_;
            autoPtr<BlendedInterfacialModel<massTransferModel>> massTransferModel2_;

            //- Interface composition models
            //interfaceCompositionModelTable interfaceCompositionModels_;
            autoPtr<interfaceCompositionModel> interfaceCompositionModel2_;

            //- The explicit part of the interfacial mass transfer rates
            //mutable iDmdtSuSpTable iDmdtSu_;
            mutable HashPtrTable<volScalarField> iDmdtSu_;

            //- The implicit part of the interfacial mass transfer rates
            //mutable iDmdtSuSpTable iDmdtSp_;
            mutable HashPtrTable<volScalarField> iDmdtSp_;

            


    // Protected Member Functions

        //- Return the interfacial mass transfer rate for a pair for a pair
        virtual tmp<volScalarField> iDmdt() const;


public:

    // Constructors

        //- Construct from fvMesh
        explicit InterfaceCompositionPhaseChangePhaseSystem(const fvMesh&, const dimensionedVector& g);


    //- Destructor
    virtual ~InterfaceCompositionPhaseChangePhaseSystem();


    // Member Functions

        //- Return the mass transfer rate for a pair
        virtual tmp<volScalarField> dmdt() const;

        //- Return the mass transfer rates for each phase
        //virtual PtrList<volScalarField> dmdts() const;

        //- Return the mass transfer matrices
        virtual autoPtr<twoPhaseSystem::massTransferTable> massTransfer() const;

        //- Add mass transfer coefficient
        //- Return mass transfer coefficient
        //tmp<volScalarField> Km1(const label) const;

        //- Return mass transfer coefficient
        //tmp<volScalarField> Km1() const;

        //- Return mass transfer coefficient
        tmp<volScalarField> Km2(const label) const;

        //- Return mass transfer coefficient
        tmp<volScalarField> Km2() const;

        //- Correct two-phase properties other than turbulence
        virtual void correct();

        //- Correct the interface temperatures
        virtual void correctInterfaceThermo();

        //- Read base phaseProperties dictionary
        virtual bool read();

        //- Return the heat transfer matrices
        //virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;
        //virtual tmp<fvScalarMatrix> heatTransfer1();

        //- Return the heat transfer matrices
        //virtual tmp<fvScalarMatrix> heatTransfer2();

        virtual autoPtr<twoPhaseSystem::heatTransferTable> heatTransfer() const;

        //virtual HashPtrTable<volScalarField> dmidtfs() const;

        //- Solve average transport of moments
        //virtual void averageTransport();
        virtual PtrList<fvVectorMatrix> getAEqns() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "InterfaceCompositionPhaseChangePhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

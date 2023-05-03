/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-26 Jeff Heylmun:    Added support of polydispersePhaseModel
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
    Foam::twoPhaseSystem

Description

SourceFiles
    twoPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseSystem_H
#define twoPhaseSystem_H

#include "IOdictionary.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "orderedPhasePair.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dragModel.H"
#include "fvMatrix.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "HashPtrTable.H"
#include "PtrListDictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class virtualMassModel;
class liftModel;
class wallLubricationModel;
class turbulentDispersionModel;
class bubblePressureModel;
// add heat transfer model
class heatTransferModel;

class blendingMethod;
template<class modelType> class BlendedInterfacialModel;

/*---------------------------------------------------------------------------*\
                      Class twoPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseSystem
:
    public IOdictionary
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Total volumetric flux
        surfaceScalarField phi_;

        //- Gravity
        const dimensionedVector g_;

        //- continuous phase
        autoPtr<phaseModel> phase2_;

        //- polydispersePhase
        autoPtr<phaseModel> phase1_;

        //- Number of nodes in the dispersed phase
        label nNodes_;

        //- Dilatation term
        volScalarField dgdt_;

        //- Optional dispersion diffusivity
        tmp<surfaceScalarField> pPrimeByA_;

        //- Unordered phase pair
        autoPtr<phasePair> pair_;

        //- Phase pair for phase 1 dispersed in phase 2
        autoPtr<orderedPhasePair> pair1In2_;

        //- Phase pair for phase 2 dispersed in phase 1
        autoPtr<orderedPhasePair> pair2In1_;

        //- Blending methods
        HashTable<autoPtr<blendingMethod>> blendingMethods_;

        //- Drag model
        autoPtr<BlendedInterfacialModel<dragModel>> drag_;

        //- Virtual mass model
        autoPtr<BlendedInterfacialModel<virtualMassModel>> virtualMass_;

        //- Lift model
        autoPtr<BlendedInterfacialModel<liftModel>> lift_;

        //- Wall lubrication model
        autoPtr<BlendedInterfacialModel<wallLubricationModel>>
            wallLubrication_;

        //- Turbulent dispersion model
        autoPtr<BlendedInterfacialModel<turbulentDispersionModel>>
            turbulentDispersion_;

        //- Bubble pressure model
        autoPtr<BlendedInterfacialModel<bubblePressureModel>>
            bubblePressure_;

        // add heat Transfer model
        //- Heat transfer model
        //autoPtr<BlendedInterfacialModel<heatTransferModel>>
        //    heatTransfer_;


    // Private member functions

        //- Return the mixture flux
        tmp<surfaceScalarField> calcPhi() const;

public:

    typedef HashPtrTable<fvScalarMatrix> heatTransferTable;

    typedef HashPtrTable<fvScalarMatrix> massTransferTable;


public:

    TypeName("twoPhaseSystem");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            twoPhaseSystem,
            dictionary,
            (
                const fvMesh& mesh, 
                const dimensionedVector& g
            ),
            (mesh, g)
        );

    // Constructors

        //- Construct from fvMesh
        twoPhaseSystem(const fvMesh&, const dimensionedVector& g);


    //- Destructor
    virtual ~twoPhaseSystem();


    // Selectors

        static autoPtr<twoPhaseSystem> New
        (
            const fvMesh& mesh,
            const dimensionedVector& g
        );


    // Member Functions

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture velocity
        tmp<volVectorField> U() const;

        //- Return the drag coefficient
        tmp<volScalarField> Kd(const label) const;

        //- Return the drag coefficient
        tmp<volScalarField> Kd() const;

        // add heat transfer coefficient
        //- Return heat transfer coefficient
        //tmp<volScalarField> Kh(const label) const;

        //- Return heat transfer coefficient
        //tmp<volScalarField> Kh() const;

        //- Return the face drag coefficient
        tmp<surfaceScalarField> Kdf(const label) const;

        //- Return the face drag coefficient
        tmp<surfaceScalarField> Kdf() const;

        //- Return the virtual mass coefficient
        tmp<volScalarField> Vm(const label) const;

        //- Return the virtual mass coefficient
        tmp<volScalarField> Vm() const;

        //- Return the face virtual mass coefficient
        tmp<surfaceScalarField> Vmf(const label) const;

        //- Return the face virtual mass coefficient
        tmp<surfaceScalarField> Vmf() const;

        //- Return the combined force (lift + wall-lubrication)
        tmp<volVectorField> F(const label) const;

        //- Return the combined force (lift + wall-lubrication)
        tmp<volVectorField> F() const;

        //- Return the combined face-force (lift + wall-lubrication)
        tmp<surfaceScalarField> Ff(const label) const;

        //- Return the combined face-force (lift + wall-lubrication)
        tmp<surfaceScalarField> Ff() const;

        //- Return the turbulent diffusivity
        //  Multiplies the phase-fraction gradient
        tmp<volScalarField> D(const label) const;

        //- Return the mean turbulent diffusivity
        //  Multiplies the phase-fraction gradient
        tmp<volScalarField> D() const;

        //- Return viscous stress matrix for phase 1
        //  Either standard or B&G based on bubble pressure
        tmp<fvVectorMatrix> divDevRhoReff1();

        //- Return viscous stress matrix for phase 2
        //  Either standard or B&G based on bubble pressure
        tmp<fvVectorMatrix> divDevRhoReff2();

        //- Solve for the two-phase-fractions
        virtual void solve();

        //- Solve average transport of moments
        virtual void averageTransport();

        //- Advect moments with deviation from mean flux
        virtual void relativeTransport();

        //- Correct two-phase properties other than turbulence
        virtual void correct();

        //- Correct two-phase turbulence
        virtual void correctTurbulence();

        //- Correct the energy transport e.g. alphat
        virtual void correctEnergyTransport();

        //- Read base phaseProperties dictionary
        virtual bool read();

        // Access

            //- Access a sub model between a phase pair
            template<class modelType>
            const modelType& lookupSubModel(const phasePair& key) const;

            //- Access a sub model between two phases
            template<class modelType>
            const modelType& lookupSubModel
            (
                const phaseModel& dispersed,
                const phaseModel& continuous
            ) const;


            //- Return the belending methods
            //HashTable<autoPtr<blendingMethod>> blendingMethods_
            const HashTable<autoPtr<blendingMethod>>& blendingMethods() const;

            //- Return the drag model for the given phase
            const dragModel& drag(const phaseModel& phase) const;

            //- Return the virtual mass model for the given phase
            const virtualMassModel& virtualMass(const phaseModel& phase) const;

            //- Return gravitational acceleration
            const dimensionedVector& g() const;

            //- Return the surface tension coefficient
            const dimensionedScalar& sigma() const;

            //const label& nNodes() const; 

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return phase model 1
            inline const phaseModel& phase1() const;

            //- Return non-const access to phase model 1
            inline phaseModel& phase1();

            //- Return phase model 2
            inline const phaseModel& phase2() const;

            //- Return non-const access to phase model 2
            inline phaseModel& phase2();

            //- Return the phase not given as an argument
            inline const phaseModel& otherPhase(const phaseModel& phase) const;

            //- Return unordered phase pair
            inline const phasePair& pair() const;

            //- Return phase 1 in phase 2 pair
            inline const orderedPhasePair& pair1In2() const;

            //- Return phase 2 in phase 1 pair
            inline const orderedPhasePair& pair2In1() const;

            //- Return the mixture flux
            inline const surfaceScalarField& phi() const;

            //- Return non-const access to the the mixture flux
            inline surfaceScalarField& phi();

            //- Return the dilatation term
            inline const volScalarField& dgdt() const;

            //- Return non-const access to the dilatation parameter
            inline volScalarField& dgdt();

            //- Return non-const access to the dispersion diffusivity
            inline tmp<surfaceScalarField>& pPrimeByA();

            //- Heat transfer
            //virtual tmp<volScalarField> Kh(const label) const = 0;

            //virtual tmp<volScalarField> Kh() const = 0;

            virtual tmp<fvScalarMatrix> heatTransfer1() = 0;

            virtual tmp<fvScalarMatrix> heatTransfer2() = 0;

            //- Mass
            //- Return the mass transfer rate for a pair
            virtual tmp<volScalarField> dmdt() const;

            //- component transport
            virtual autoPtr<massTransferTable> massTransfer() const;

            virtual autoPtr<heatTransferTable> heatTransfer() const = 0;

            //- Correct the kinematics
            virtual void correctKinematics();

            //- AEqns
            virtual PtrList<fvVectorMatrix> getAEqns() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "twoPhaseSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

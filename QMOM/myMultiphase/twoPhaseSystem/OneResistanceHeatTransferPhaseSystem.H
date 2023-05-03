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
    Foam::OneResistanceHeatTransferPhaseSystem

Description
    Class which models interfacial heat transfer between a number of phases. A
    single heat transfer model is used for each interface.

See also
    TwoResistanceHeatTransferPhaseSystem

SourceFiles
    OneResistanceHeatTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef OneResistanceHeatTransferPhaseSystem_H
#define OneResistanceHeatTransferPhaseSystem_H

#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blendingMethod;
template<class modelType> class BlendedInterfacialModel;
class heatTransferModel;

/*---------------------------------------------------------------------------*\
                    Class OneResistanceHeatTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class OneResistanceHeatTransferPhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected typedefs
        /*
        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<heatTransferModel>>,
            phasePairKey,
            phasePairKey::hash
        > heatTransferModelTable;
        */

    // Protected data

        // Sub Models

            //- Heat transfer models
            //heatTransferModelTable heatTransferModels_;
            autoPtr<BlendedInterfacialModel<heatTransferModel>>
            heatTransfer_;


public:

    // Constructors

        //- Construct from fvMesh
        OneResistanceHeatTransferPhaseSystem(const fvMesh&, const dimensionedVector& g);


    //- Destructor
    virtual ~OneResistanceHeatTransferPhaseSystem();


    // Member Functions

        // add heat transfer coefficient
        //- Return heat transfer coefficient
        tmp<volScalarField> Kh(const label) const;

        //- Return heat transfer coefficient
        tmp<volScalarField> Kh() const;

        //- Return the heat transfer matrices
        //virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;
        virtual tmp<fvScalarMatrix> heatTransfer1();

        //- Return the heat transfer matrices
        virtual tmp<fvScalarMatrix> heatTransfer2();

        //- Read base phaseProperties dictionary
        virtual bool read();

        virtual autoPtr<twoPhaseSystem::heatTransferTable> heatTransfer() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "OneResistanceHeatTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

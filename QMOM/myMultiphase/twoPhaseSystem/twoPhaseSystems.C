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

#include "addToRunTimeSelectionTable.H"

#include "twoPhaseSystem.H"
#include "OneResistanceHeatTransferPhaseSystem.H"
#include "TwoResistanceHeatTransferPhaseSystem.H"
#include "InterfaceCompositionPhaseChangePhaseSystem.H"
#include "QMOMInterfaceCompositionPhaseChangePhaseSystem.H"
#include "polyQMOMInterfaceCompositionPhaseChangePhaseSystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*
    typedef twoPhaseSystem basicTwoPhaseSystem;

    addNamedToRunTimeSelectionTable
    (
        twoPhaseSystem,
        basicTwoPhaseSystem,
        dictionary,
        basicTwoPhaseSystem
    );
    */

    typedef OneResistanceHeatTransferPhaseSystem<twoPhaseSystem> 
    oneResistanceTwoPhaseSystem;

    addNamedToRunTimeSelectionTable
    (
        twoPhaseSystem,
        oneResistanceTwoPhaseSystem,
        dictionary,
        oneResistanceTwoPhaseSystem
    );

    typedef TwoResistanceHeatTransferPhaseSystem<twoPhaseSystem> 
    twoResistanceTwoPhaseSystem;

    addNamedToRunTimeSelectionTable
    (
        twoPhaseSystem,
        twoResistanceTwoPhaseSystem,
        dictionary,
        twoResistanceTwoPhaseSystem
    );

    typedef InterfaceCompositionPhaseChangePhaseSystem
    <
        TwoResistanceHeatTransferPhaseSystem<twoPhaseSystem>
    > 
    interfaceCompositionTwoPhaseSystem;

    addNamedToRunTimeSelectionTable
    (
        twoPhaseSystem,
        interfaceCompositionTwoPhaseSystem,
        dictionary,
        interfaceCompositionTwoPhaseSystem
    );


    typedef QMOMInterfaceCompositionPhaseChangePhaseSystem
    <
        TwoResistanceHeatTransferPhaseSystem<twoPhaseSystem>
    > 
    QMOMInterfaceCompositionTwoPhaseSystem;

    addNamedToRunTimeSelectionTable
    (
        twoPhaseSystem,
        QMOMInterfaceCompositionTwoPhaseSystem,
        dictionary,
        QMOMInterfaceCompositionTwoPhaseSystem
    );

    typedef polyQMOMInterfaceCompositionPhaseChangePhaseSystem
    <
        TwoResistanceHeatTransferPhaseSystem<twoPhaseSystem>
    > 
    polyQMOMInterfaceCompositionTwoPhaseSystem;

    addNamedToRunTimeSelectionTable
    (
        twoPhaseSystem,
        polyQMOMInterfaceCompositionTwoPhaseSystem,
        dictionary,
        polyQMOMInterfaceCompositionTwoPhaseSystem
    );
}


// ************************************************************************* //

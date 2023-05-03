/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"

#include "phaseModel.H"
#include "polydispersePhaseModel.H"
#include "monodispersePhaseModel.H"
#include "multiComponentPhaseModel.H"
#include "rhoThermo.H"
#include "rhoReactionThermo.H"
#include "purePhaseModel.H"
#include "thermoPhaseModel.H"

#include "../polydispersePhaseModel/myEvapLPNormDBCoBCPolydispersePhaseModel.H"
#include "../polydispersePhaseModel/myEvapULPNormDBCoBCPolydispersePhaseModel.H"
#include "../pbeSMDLimitedDBasedM3PhaseModel/pbeSMDLimitedDBasedM3PhaseModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
    purePhaseModel
    <
        thermoPhaseModel<phaseModel, rhoThermo>
    >
    
    continuousPurePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        continuousPurePhaseModel,
        dictionary,
        continuousPurePhaseModel
    );

    typedef
    multiComponentPhaseModel
    <
        thermoPhaseModel<phaseModel, rhoReactionThermo>
    >
    
    continuousMultiComponentPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        continuousMultiComponentPhaseModel,
        dictionary,
        continuousMultiComponentPhaseModel
    );

    typedef
    polydispersePhaseModel
    <
        purePhaseModel
        <
            thermoPhaseModel<phaseModel, rhoThermo>
        >
    >
    polydispersePurePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        polydispersePurePhaseModel,
        dictionary,
        polydispersePurePhaseModel
    );
//
/*
    typedef
    myEvapLPNormDBCoBCPolydispersePhaseModel
    <
        purePhaseModel
        <
            thermoPhaseModel<phaseModel, rhoThermo>
        >
    >
    myEvapLPNormDBCoBCPolydispersePurePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        myEvapLPNormDBCoBCPolydispersePurePhaseModel,
        dictionary,
        myEvapLPNormDBCoBCPolydispersePurePhaseModel
    );

    typedef
    myEvapULPNormDBCoBCPolydispersePhaseModel
    <
        purePhaseModel
        <
            thermoPhaseModel<phaseModel, rhoThermo>
        >
    >
    myEvapULPNormDBCoBCPolydispersePurePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        myEvapULPNormDBCoBCPolydispersePurePhaseModel,
        dictionary,
        myEvapULPNormDBCoBCPolydispersePurePhaseModel
    );

//
    typedef
    monodispersePhaseModel
    <
        purePhaseModel
        <
            thermoPhaseModel<phaseModel, rhoThermo>
        >
    >
    monodispersePurePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        monodispersePurePhaseModel,
        dictionary,
        monodispersePurePhaseModel
    );

    typedef
    pbeSMDLimitedDBasedM3PhaseModel
    <
        monodispersePhaseModel
        <
            purePhaseModel
            <
                thermoPhaseModel<phaseModel, rhoThermo>
            >
        >
    >
    twoWayPBESMDLimitedDBasedM3PurePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        twoWayPBESMDLimitedDBasedM3PurePhaseModel,
        dictionary,
        twoWayPBESMDLimitedDBasedM3PurePhaseModel
    );
*/
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-04 Jeff Heylmun:    Modified drag coefficient calculation
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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

#include "TomiyamaCorrelated.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(TomiyamaCorrelated, 0);
    addToRunTimeSelectionTable(dragModel, TomiyamaCorrelated, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::TomiyamaCorrelated::TomiyamaCorrelated
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    A_("A", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::TomiyamaCorrelated::~TomiyamaCorrelated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::TomiyamaCorrelated::CdRe
(
    const label nodei,
    const label nodej
) const
{
    volScalarField Re(pair_.Re(nodei, nodej));
    volScalarField Eo(pair_.Eo(nodei, nodej));

    return
        max
        (
            A_
           *min
            (
                (1 + 0.15*pow(Re, 0.687)),
                scalar(3)
            ),
            8*Eo*Re/(3*Eo + 12)
        );
}


Foam::tmp<Foam::volScalarField>
Foam::dragModels::TomiyamaCorrelated::CdRe() const
{
    volScalarField Re(pair_.Re());
    volScalarField Eo(pair_.Eo());

    return
        max
        (
            A_
           *min
            (
                (1 + 0.15*pow(Re, 0.687)),
                scalar(3)
            ),
            8*Eo*Re/(3*Eo + 12)
);
}


// ************************************************************************* //

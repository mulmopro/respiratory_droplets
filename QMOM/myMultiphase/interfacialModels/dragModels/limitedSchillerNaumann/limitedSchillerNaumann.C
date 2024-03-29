/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "limitedSchillerNaumann.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(limitedSchillerNaumann, 0);
    addToRunTimeSelectionTable(dragModel, limitedSchillerNaumann, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::limitedSchillerNaumann::limitedSchillerNaumann
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::limitedSchillerNaumann::~limitedSchillerNaumann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::limitedSchillerNaumann::CdRe
(
    const label nodei,
    const label nodej
) const
{
    volScalarField Re(pair_.Re(nodei, nodej));

    return
        neg(Re - 1000)*24.0*(1.0 + 0.15*pow(max(Re, residualRe_), 0.687))
      + pos0(Re - 1000)*0.44*max(Re, residualRe_);
}


Foam::tmp<Foam::volScalarField> Foam::dragModels::limitedSchillerNaumann::CdRe() const
{
    volScalarField Re(pair_.Re());

    return
        neg(Re - 1000)*24.0*(1.0 + 0.15*pow(max(Re, residualRe_), 0.687))
      + pos0(Re - 1000)*0.44*max(Re, residualRe_);
}


// ************************************************************************* //

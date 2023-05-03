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

#include "limitedFrossling.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferModels
{
    defineTypeNameAndDebug(limitedFrossling, 0);
    addToRunTimeSelectionTable(massTransferModel, limitedFrossling, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferModels::limitedFrossling::limitedFrossling
(
    const dictionary& dict,
    const phasePair& pair
)
:
    massTransferModel(dict, pair),
    Le_("Le", dimless, dict),
    dmin_("dmin", dimLength, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massTransferModels::limitedFrossling::~limitedFrossling()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::limitedFrossling::K
(
    const label nodei,
    const label nodej
) const
{
    //volScalarField Sh(2 + 0.552*sqrt(pair_.Re())*cbrt(Le_*pair_.Pr()));

    //return 6*pair_.dispersed()*Sh/sqr(pair_.dispersed().d());

    volScalarField Sh(2 + 0.552*sqrt(pair_.Re(nodei, nodej))*cbrt(Le_*pair_.Pr()));
    
    return pos0(pair_.dispersed().ds(nodei)-dmin_)*6*pair_.dispersed().alphas(nodei)*Sh/sqr(pair_.dispersed().ds(nodei));
}


// ************************************************************************* //

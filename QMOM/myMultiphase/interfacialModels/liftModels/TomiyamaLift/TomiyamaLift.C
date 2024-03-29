/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenFOAM Foundation
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

#include "TomiyamaLift.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(TomiyamaLift, 0);
    addToRunTimeSelectionTable(liftModel, TomiyamaLift, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::TomiyamaLift::TomiyamaLift
(
    const dictionary& dict,
    const phasePair& pair
)
:
    liftModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::TomiyamaLift::~TomiyamaLift()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::TomiyamaLift::Cl
(
    const label nodei,
    const label nodej
) const
{
    volScalarField EoH(pair_.EoH2(nodei, nodej));

    volScalarField f
    (
        0.0010422*pow3(EoH) - 0.0159*sqr(EoH) - 0.0204*EoH + 0.474
    );

    return
        neg(EoH - scalar(4))*min(0.288*tanh(0.121*pair_.Re(nodei, nodej)), f)
      + pos0(EoH - scalar(4))*neg(EoH - scalar(10.7))*f
      + pos0(EoH - scalar(10.7))*(-0.288);
}


// ************************************************************************* //

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

#include "Frank.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(Frank, 0);
    addToRunTimeSelectionTable
    (
        wallLubricationModel,
        Frank,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::Frank::Frank
(
    const dictionary& dict,
    const phasePair& pair
)
:
    wallLubricationModel(dict, pair),
    Cwd_("Cwd", dimless, dict),
    Cwc_("Cwc", dimless, dict),
    p_(readScalar(dict.lookup("p")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModels::Frank::~Frank()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::wallLubricationModels::Frank::Fi
(
    const label nodei,
    const label nodej
) const
{
    volVectorField Ur(pair_.Ur(nodei, nodej));

    const volVectorField& n(nWall());
    const volScalarField& y(yWall());

    volScalarField Eo(pair_.Eo(nodei, nodej));
    volScalarField yTilde(y/(Cwc_*pair_.dispersed().ds(nodei)));

    return zeroGradWalls
    (
        (
            pos0(Eo - 1.0)*neg(Eo - 5.0)*exp(-0.933*Eo + 0.179)
          + pos0(Eo - 5.0)*neg(Eo - 33.0)*(0.00599*Eo - 0.0187)
          + pos0(Eo - 33.0)*0.179
        )
       *max
        (
            dimensionedScalar("zero", dimless/dimLength, 0.0),
            (1.0 - yTilde)/(Cwd_*y*pow(yTilde, p_ - 1.0))
        )
       *pair_.continuous().rho()
       *magSqr(Ur - (Ur & n)*n)
       *n
    );
}


// ************************************************************************* //

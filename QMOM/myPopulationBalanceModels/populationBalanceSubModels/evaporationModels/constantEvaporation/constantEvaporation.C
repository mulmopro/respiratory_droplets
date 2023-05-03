/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015 by Matteo Icardi and 2017 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
-------------------------------------------------------------------------------
2017-03-28 Alberto Passalacqua: Adapted to single scalar calculation.
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "constantEvaporation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace evaporationModels
{
    defineTypeNameAndDebug(constantEvaporation, 0);

    addToRunTimeSelectionTable
    (
        evaporationModel,
        constantEvaporation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::evaporationModels::constantEvaporation
::constantEvaporation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    evaporationModel(dict, mesh),
    Cg_(dict.lookupOrDefault("Cg", scalar(1.0)))
    /*
    minAbscissa_(dict.lookupOrDefault("minAbscissa", scalar(0))),
    maxAbscissa_(dict.lookupOrDefault("maxAbscissa", GREAT))
    */
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::evaporationModels::constantEvaporation
::~constantEvaporation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::evaporationModels::constantEvaporation::Kg
(
    const scalar& abscissa,
    const label celli,
    const bool lengthBased,
    const label environment
) const
{
    return Cg_*pos0(-abscissa + maxAbscissa_)
            *pos0(abscissa - minAbscissa_);
}

// ************************************************************************* //

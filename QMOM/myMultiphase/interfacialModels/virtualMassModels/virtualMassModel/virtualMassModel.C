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

#include "virtualMassModel.H"
#include "phasePair.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(virtualMassModel, 0);
    defineRunTimeSelectionTable(virtualMassModel, dictionary);
}

const Foam::dimensionSet Foam::virtualMassModel::dimK(dimDensity);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMassModel::virtualMassModel
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::virtualMassModel::~virtualMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::virtualMassModel::Ki
(
    const label nodei,
    const label nodej
) const
{
    return Cvm(nodei, nodej)*pair_.continuous().rho();
}

Foam::tmp<Foam::volScalarField> Foam::virtualMassModel::Ki() const
{
    return Cvm()*pair_.continuous().rho();
}


Foam::tmp<Foam::volScalarField> Foam::virtualMassModel::K
(
    const label nodei,
    const label nodej
) const
{
    return pair_.dispersed().alphas(nodei)*Ki(nodei, nodej);
}


Foam::tmp<Foam::volScalarField> Foam::virtualMassModel::K() const
{
    return pair_.dispersed()*Ki();
}


Foam::tmp<Foam::surfaceScalarField> Foam::virtualMassModel::Kf
(
    const label nodei,
    const label nodej
) const
{
    return
        fvc::interpolate(pair_.dispersed().alphas(nodei))
       *fvc::interpolate(Ki(nodei, nodej));
}


Foam::tmp<Foam::surfaceScalarField> Foam::virtualMassModel::Kf() const
{
    return fvc::interpolate(pair_.dispersed())*fvc::interpolate(Ki());
}


bool Foam::virtualMassModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //

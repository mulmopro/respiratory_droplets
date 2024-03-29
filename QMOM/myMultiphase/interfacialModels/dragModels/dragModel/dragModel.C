/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
2017-05-24 Jeff Heylmun:    Added return functions for acceleration
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

#include "dragModel.H"
#include "phasePair.H"
#include "swarmCorrection.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dragModel, 0);
    defineRunTimeSelectionTable(dragModel, dictionary);
}

const Foam::dimensionSet Foam::dragModel::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModel::dragModel
(
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


Foam::dragModel::dragModel
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
    pair_(pair),
    swarmCorrection_
    (
        swarmCorrection::New
        (
            dict.subDict("swarmCorrection"),
            pair
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModel::~dragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModel::Ki
(
    const label nodei,
    const label nodej
) const
{
    return
        0.75
       *CdRe(nodei, nodej)
       *swarmCorrection_->Cs(nodei, nodej)
       *pair_.continuous().rho()
       *pair_.continuous().nu()
       /sqr(pair_.dispersed().ds(nodei));
}


Foam::tmp<Foam::volScalarField> Foam::dragModel::Ki() const
{
    return
        0.75
       *CdRe()
       *swarmCorrection_->Cs()
       *pair_.continuous().rho()
       *pair_.continuous().nu()
       /sqr(pair_.dispersed().d());
}


Foam::tmp<Foam::volScalarField> Foam::dragModel::K
(
    const label nodei,
    const label nodej
) const
{
    if (pair_.continuous().nNodes() > 1)
    {
        //  Scale drag so that the sum of the drag forces from a monodisperse
        //  phase is not counted nNodes times
        return
            max
            (
                pair_.dispersed().alphas(nodei)
               *pair_.continuous().alphas(nodej)
               /max
                (
                    pair_.continuous(),
                    pair_.continuous().residualAlpha()
                ),
                pair_.dispersed().residualAlpha()/pair_.dispersed().nNodes()
            )*Ki(nodei, nodej);
    }

    return
        max
        (
            pair_.dispersed().alphas(nodei),
            pair_.dispersed().residualAlpha()
           /max
            (
                pair_.dispersed().nNodes(),
                pair_.continuous().nNodes()
            )
        )*Ki(nodei, nodej);
}


Foam::tmp<Foam::volScalarField> Foam::dragModel::K() const
{
    return
        max
        (
            pair_.dispersed(),
            pair_.dispersed().residualAlpha()
        )*Ki();
}


Foam::tmp<Foam::surfaceScalarField> Foam::dragModel::Kf
(
    const label nodei,
    const label nodej
) const
{
    return
        max
        (
            fvc::interpolate(pair_.dispersed().alphas(nodei)),
            pair_.dispersed().residualAlpha()
           /max
            (
                pair_.dispersed().nNodes(),
                pair_.continuous().nNodes()
            )
        )*fvc::interpolate(Ki(nodei, nodej));
}


Foam::tmp<Foam::surfaceScalarField> Foam::dragModel::Kf() const
{
    return
        max
        (
            fvc::interpolate(pair_.dispersed()),
            pair_.dispersed().residualAlpha()
        )*fvc::interpolate(Ki());
}


bool Foam::dragModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //

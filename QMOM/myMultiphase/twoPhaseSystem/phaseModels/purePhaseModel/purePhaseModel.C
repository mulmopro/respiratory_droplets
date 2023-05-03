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

#include "purePhaseModel.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::purePhaseModel<BasePhaseModel>::purePhaseModel
(
    const twoPhaseSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseProperties, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::purePhaseModel<BasePhaseModel>::~purePhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
bool Foam::purePhaseModel<BasePhaseModel>::pure() const
{
    return true;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::purePhaseModel<BasePhaseModel>::YiEqn(volScalarField& Yi)
{
    FatalErrorInFunction
        << "Cannot construct a species fraction equation for a pure phase"
        << exit(FatalError);

    return tmp<fvScalarMatrix>();
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::purePhaseModel<BasePhaseModel>::Y() const
{
    // Y_ has never been set, so we are returning an empty list

    return Y_;
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::purePhaseModel<BasePhaseModel>::Y(const word& name) const
{
    FatalErrorInFunction
        << "Cannot get a species fraction by name from a pure phase"
        << exit(FatalError);

    return NullObjectRef<volScalarField>();
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::purePhaseModel<BasePhaseModel>::YRef()
{
    FatalErrorInFunction
        << "Cannot access the species fractions of for a pure phase"
        << exit(FatalError);

    return Y_;
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::purePhaseModel<BasePhaseModel>::YActive() const
{
    // Y_ has never been set, so we are returning an empty list

    return Y_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::purePhaseModel<BasePhaseModel>::YActiveRef()
{
    FatalErrorInFunction
        << "Cannot access the species fractions of for a pure phase"
        << exit(FatalError);

    return Y_;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::phaseModel& Foam::phasePair::phase1() const
{
    return phase1_;
}


inline const Foam::phaseModel& Foam::phasePair::phase2() const
{
    return phase2_;
}


inline const Foam::phaseModel& Foam::phasePair::otherPhase
(
    const phaseModel& phase
) const
{
    if (&phase1_ == &phase)
    {
        return phase2_;
    }
    else if (&phase2_ == &phase)
    {
        return phase1_;
    }
    else
    {
        FatalErrorInFunction
            << "this phasePair does not contain phase " << phase.name()
            << exit(FatalError);

        return phase;
    }
}


inline const Foam::dimensionedVector& Foam::phasePair::g() const
{
    return g_;
}


inline const Foam::dimensionedScalar& Foam::phasePair::sigma() const
{
    return sigma_;
}


// ************************************************************************* //

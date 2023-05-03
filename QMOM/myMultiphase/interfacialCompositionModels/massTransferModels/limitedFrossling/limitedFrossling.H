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

Class
    Foam::massTransferModels::limitedFrossling

Description
    limitedFrossling correlation for turbulent mass transfer from the surface of a
    sphere to the surrounding fluid.

SourceFiles
    limitedFrossling.C

\*---------------------------------------------------------------------------*/

#ifndef limitedFrossling_H
#define limitedFrossling_H

#include "massTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace massTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class limitedFrossling Declaration
\*---------------------------------------------------------------------------*/

class limitedFrossling
:
    public massTransferModel
{
    // Private data

        //- Lewis number
        const dimensionedScalar Le_;
        const dimensionedScalar dmin_;


public:

    //- Runtime type information
    TypeName("limitedFrossling");


    // Constructors

        //- Construct from components
        limitedFrossling
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~limitedFrossling();


    // Member Functions

        //- The implicit mass transfer coefficient
        virtual tmp<volScalarField> K(const label nodei, const label nodej) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Class
    Foam::wallLubricationModels::noWallLubrication

Description

SourceFiles
    noWallLubrication.C

\*---------------------------------------------------------------------------*/

#ifndef noWallLubrication_H
#define noWallLubrication_H

#include "wallLubricationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace wallLubricationModels
{

/*---------------------------------------------------------------------------*\
                           Class noWallLubrication Declaration
\*---------------------------------------------------------------------------*/

class noWallLubrication
:
    public wallLubricationModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noWallLubrication
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~noWallLubrication();


    // Member Functions

        //- Return phase-intensive wall lubrication force
        tmp<volVectorField> Fi(const label, const label) const;

        //- Wall lubrication force
        tmp<volVectorField> F(const label, const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

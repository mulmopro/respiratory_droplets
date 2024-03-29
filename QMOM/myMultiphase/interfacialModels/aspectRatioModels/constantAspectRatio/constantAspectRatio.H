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
    Foam::constantAspectRatio

Description
    Constant value aspect ratio model.

SourceFiles
    constantAspectRatio.C

\*---------------------------------------------------------------------------*/

#ifndef constantAspectRatio_H
#define constantAspectRatio_H

#include "aspectRatioModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aspectRatioModels
{

/*---------------------------------------------------------------------------*\
                         Class constantAspectRatio Declaration
\*---------------------------------------------------------------------------*/

class constantAspectRatio
:
    public aspectRatioModel
{
    // Private data

        //- Constant aspect ratio value
        const dimensionedScalar E0_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary and an ordered phase pair
        constantAspectRatio
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~constantAspectRatio();


    // Member Functions

        //- Aspect ratio
        virtual tmp<volScalarField> E
        (
            const label,
            const label
        ) const;

        //- Mean aspect ratio
        virtual tmp<volScalarField> E() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aspectRatioModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::aspectRatioModels::Wellek

Description
    Aspect ratio model of Wellek et al.

    References:
    \verbatim
        "Implementation and Comparison of Correlations for interfacial Forces
         in a Gas-Liquid System within an Euler-Euler Framework"
        Otromke, M.,
        PhD Thesis, April 2013
    \endverbatim

    \verbatim
        "Shape of liquid drops moving in liquid media"
        Wellek, R.M., Agrawal, A.K., Skelland, A.H.P.,
        International Journal of Multiphase Flow
        Volume 12, Issue 5, September 1966, pp. 854-862
    \endverbatim

SourceFiles
    Wellek.C

\*---------------------------------------------------------------------------*/

#ifndef Wellek_H
#define Wellek_H

#include "aspectRatioModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aspectRatioModels
{

/*---------------------------------------------------------------------------*\
                       Class Wellek Declaration
\*---------------------------------------------------------------------------*/

class Wellek
:
    public aspectRatioModel
{
public:

    //- Runtime type information
    TypeName("Wellek");


    // Constructors

        //- Construct from a dictionary and an ordered phase pair
        Wellek
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Wellek();


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

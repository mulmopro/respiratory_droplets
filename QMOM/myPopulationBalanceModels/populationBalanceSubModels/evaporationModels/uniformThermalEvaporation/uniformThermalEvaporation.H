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
    Copyright (C) 2019 Alberto Passalacqua
-------------------------------------------------------------------------------
2017-03-28 Alberto Passalacqua: Adapted to single scalar calculation.
2020-03-12 Alberto Passalacqua: Added access to fvMesh.
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

Class
    Foam::populationBalanceSubModels::evaporationModels::uniformThermalEvaporation

Description
    Constant growth model. The value of the constant is read from dictionary
    or assumed equal to 1. The minimum and maximum value of the abscissa to
    which the model is applied can be defined by the user or it is assumed to
    be 1.0e-10 and 1.0e10.

SourceFiles
    uniformThermalEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef uniformThermalEvaporation_H
#define uniformThermalEvaporation_H

#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace evaporationModels
{

/*---------------------------------------------------------------------------*\
                    Class uniformThermalEvaporation Declaration
\*---------------------------------------------------------------------------*/

class uniformThermalEvaporation
:
    public evaporationModel
{
private:

    // Private data

        //- Minimum value of the abscissa to apply the model
        //scalar minAbscissa_;

        //- Maximum value of the abscissa to apply the model
        //scalar maxAbscissa_;

        //scalar Cg_;

public:

        //- Runtime type information
        TypeName("uniformThermal");


    // Constructors

        //- Construct from components
        uniformThermalEvaporation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~uniformThermalEvaporation();


    // Member Functions

        //- Growth model
        virtual scalar Kg
        (
            const scalar& abscissa,
            const label celli,
            const bool lengthBased,
            const label environment = 0
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace evaporationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015 by Matteo Icardi
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2020 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::evaporationModel

Description
    Abstract class for growth model.

SourceFiles
    evaporationModel.C

\*---------------------------------------------------------------------------*/

#ifndef evaporationModel_H
#define evaporationModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{

/*---------------------------------------------------------------------------*\
                    Class evaporationModel Declaration
\*---------------------------------------------------------------------------*/

class evaporationModel
{
protected:

    // Protected data

        const dictionary& dict_;

        const fvMesh& mesh_;

        //- Coefficient of aggregation kernel
        //const dimensionedScalar Cg_;

        //- Minimum value of the abscissa to apply the model
        scalar minAbscissa_;

        //- Maximum value of the abscissa to apply the model
        scalar maxAbscissa_;


public:

    //- Runtime type information
    TypeName("evaporationModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        evaporationModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        evaporationModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        evaporationModel(const evaporationModel&) = delete;


    // Selectors

        static autoPtr<evaporationModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~evaporationModel();


    //- Return the derivative of length with respect to volume
    /*
    static scalar dDdV(const scalar& d)
    {
        return
            2.0
            /max(Foam::constant::mathematical::pi*sqr(d), SMALL);
    }
    */
    // Member Functions

        //- Growth model
        virtual scalar Kg
        (
            const scalar& abscissa,
            const label celli,
            const bool lengthBased,
            const label environment = 0
        ) const = 0;

        scalar minAbscissa()
        {
            return minAbscissa_;
        }

        scalar maxAbscissa()
        {
            return maxAbscissa_;
        }

/*
        //- Return phase space convection
        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& quadrature
        );

        //- Return phase space convection
        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature
        );
*/

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const evaporationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

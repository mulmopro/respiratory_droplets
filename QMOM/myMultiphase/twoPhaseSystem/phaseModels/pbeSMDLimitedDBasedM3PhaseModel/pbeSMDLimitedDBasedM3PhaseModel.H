/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
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

Class
    Foam::pbeSMDLimitedDBasedM3PhaseModel

Description
    Class for a polydisperse phase model. Each bubble size hase a unique mean
    velocity. Size and velocity moments are stored in quadrature.

SourceFiles
    pbeSMDLimitedDBasedM3PhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef pbeSMDLimitedDBasedM3PhaseModel_H
#define pbeSMDLimitedDBasedM3PhaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "phaseModel.H"
#include "monoKineticQuadratureApproximation.H"
#include "populationBalanceModel.H"

#include "coalescenceKernel.H"
#include "breakupKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class pbeSMDLimitedDBasedM3PhaseModel
:
    public BasePhaseModel
{
    // Private data

        //- Population balance properties
        IOdictionary pbeDict_;

        //- Sauter Mean diameter
        volScalarField dsm_;

        //- Maximum diameter
        dimensionedScalar maxD_;

        //- Minimum diameter
        dimensionedScalar minD_;


protected:

    autoPtr<populationBalanceModel> populationBalance_;

public:

    // Constructors
        pbeSMDLimitedDBasedM3PhaseModel
        (
            const twoPhaseSystem& fluid,
            const dictionary& phaseProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~pbeSMDLimitedDBasedM3PhaseModel();


    // Member Functions


        //- Return the diameter for nodei
        virtual const volScalarField& ds(const label nodei) const
        {
            //return ds_[nodei];
            return dsm_;
            //return this->d_;
        }

        //- Calculate Sauter mean diameter
        //tmp<volScalarField> calcDsm();
        virtual void correctDsm();

        //- solve population balance equation
        virtual void solve();

        //- Constant access to alpha field for nodei
        virtual tmp<volScalarField> volumeFraction() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pbeSMDLimitedDBasedM3PhaseModel.C"
#endif

#endif

// ************************************************************************* //

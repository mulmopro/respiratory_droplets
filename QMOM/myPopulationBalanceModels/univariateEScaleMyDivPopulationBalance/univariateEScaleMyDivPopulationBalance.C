/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
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

\*---------------------------------------------------------------------------*/

#include "univariateEScaleMyDivPopulationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(univariateEScaleMyDivPopulationBalance, 0);
    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        univariateEScaleMyDivPopulationBalance,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::univariateEScaleMyDivPopulationBalance
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    univariateMyDivPDFTransportModel(name, dict, phi.mesh(), phi, "RPlus"),
    populationBalanceModel(name, dict, phi),
    odeType(phi.mesh(), dict),
    aggregation_(dict.lookupOrDefault("aggregation", false)),
    breakup_(dict.lookupOrDefault("breakup", false)),
    growth_(dict.lookupOrDefault("growth", false)),
    nucleation_(dict.lookupOrDefault("nucleation", false)),
    evaporation_(dict.lookupOrDefault("evaporation", false)),
    aggregationKernel_(),
    breakupKernel_(),
    growthModel_(),
    diffusionModel_
    (
        Foam::populationBalanceSubModels::diffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    ),
    nucleationModel_(),
    evaporationModel_()
{
    if (aggregation_)
    {
        aggregationKernel_ =
            Foam::populationBalanceSubModels::aggregationKernel::New
            (
                dict.subDict("aggregationKernel"),
                phi_.mesh()
            );
    }

    if (breakup_)
    {
        breakupKernel_ =
            Foam::populationBalanceSubModels::breakupKernel::New
            (
                dict.subDict("breakupKernel"),
                phi_.mesh()
            );
    }

    if (growth_)
    {
        growthModel_ =
            Foam::populationBalanceSubModels::growthModel::New
            (
                dict.subDict("growthModel"),
                phi_.mesh()
            );
    }

    if (nucleation_)
    {
        nucleationModel_ =
            Foam::populationBalanceSubModels::nucleationModel::New
            (
                dict.subDict("nucleationModel"),
                phi_.mesh()
            );
    }

    if (evaporation_)
    {
        evaporationModel_ =
            Foam::populationBalanceSubModels::evaporationModel::New
            (
                dict.subDict("evaporationModel"),
                phi_.mesh()
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::~univariateEScaleMyDivPopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::implicitMomentSource
(
    const volScalarMoment& moment
)
{
    return diffusionModel_->momentDiff(moment);
}


void
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::updateCellMomentSource(const label)
{}


Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label environment
)
{
    scalar source(0);

    if (aggregation_)
    {
        source +=
            aggregationKernel_->aggregationSource
            (
                momentOrder,
                celli,
                quadrature,
                environment
            );
    }

    if (breakup_)
    {
        source +=
            breakupKernel_->breakupSource
            (
                momentOrder,
                celli,
                quadrature
            );
    }

    if (growth_)
    {
        source +=
            growthModel_->phaseSpaceConvection
            (
                momentOrder,
                celli,
                quadrature
            );
    }

    if (nucleation_)
    {
        source += nucleationModel_->nucleationSource(momentOrder[0], celli);
    }

    return source;
}


Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::realizableCo() const
{
    return univariateMyDivPDFTransportModel::realizableCo();
}


Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::CoNum() const
{
    return 0.0;
}


bool
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::solveMomentSources() const
{
    /*
    if (evaporation_)
    {
        PtrList<volScalarNode> nodes = quadrature_.nodes();
        label sizeIndex = nodes[0].sizeIndex();
        bool lengthBased = nodes[0].lengthBased();
        //- GlobalDt
        scalar globalDt = quadrature_.moments()[0].mesh().time().deltaT().value();
        //- Minimum abscissa
        const scalar minAbscissa = evaporationModel_->minAbscissa();

        //- If QMOM
        if(!nodes[0].extended())
        {
            forAll(quadrature_.moments()[1], celli)
            {
                forAll(nodes, pNodeI)
                {
                    volScalarNode& node = nodes[pNodeI];
                    scalar& abscissa = node.primaryAbscissae()[sizeIndex][celli];

                    abscissa = abscissa - evaporationModel_->Kg(abscissa, lengthBased) * globalDt;

                    if (abscissa <= minAbscissa)
                    {
                        abscissa = minAbscissa;
                    }
                }
            }
            quadrature_.updateMoments();
            quadrature_.updateQuadrature();
            
            //return;
        }
        else
        {
            forAll(quadrature_.moments()[1], celli)
            {
                forAll(nodes, pNodeI)
                {
                    volScalarNode& node = nodes[pNodeI];
                    forAll(node.secondaryWeights()[sizeIndex], sNodei)
                    {
                        scalar& abscissa = node.secondaryAbscissae()[sizeIndex][sNodei][celli];

                        abscissa = abscissa - evaporationModel_->Kg(abscissa, lengthBased) * globalDt;

                        if (abscissa <= minAbscissa)
                        {
                            abscissa = minAbscissa;
                        }
                    }
                }
            }
            quadrature_.updateMoments();
            quadrature_.updateQuadrature();
            
            //return;
        }
    }
    */
    if (aggregation_ || breakup_ || growth_ || nucleation_)
    {
        return odeType::solveSources_;
    }

    return false;
}


bool
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::solveMomentOde() const
{
    return odeType::solveOde_;
}


void
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::explicitMomentSource()
{
    odeType::solve(quadrature_, 0);
}


void 
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::evaporation()
{
    if (evaporation_)
    {
        PtrList<volScalarNode>& nodes = quadrature_.nodes();
        label sizeIndex = nodes[0].sizeIndex();
        //bool lengthBased = nodes[0].lengthBased();
        //- GlobalDt
        scalar globalDt = quadrature_.moments()[0].mesh().time().deltaT().value();
        //- Minimum abscissa
        scalar minAbscissa = evaporationModel_->minAbscissa();

        const volScalarField& dvdtf = quadrature_.moments()[0].mesh().lookupObject<volScalarField>("dvdtf");

        const volScalarField& alpha_liquid = quadrature_.moments()[0].mesh().lookupObject<volScalarField>("alpha.liquid");

        dimensionedScalar MSmall("mSmall", inv(dimArea), SMALL);
        dimensionedScalar ASmall("ASmall", dimless, 1e-8);

        volScalarField scale
        (
            Foam::constant::mathematical::pi * quadrature_.moments()[3] / 6.0
            /Foam::max
            (
                alpha_liquid,
                ASmall
            )
        );

        //scale.min(1.0);
        Info<< "scale: "
            << "min: " << min(scale)
            << " average: " << average(scale)
            << " max: " << max(scale) << endl;

        volScalarField drift = 2 * dvdtf / (Foam::constant::mathematical::pi * max(quadrature_.moments()[1], MSmall)) * scale;

        //- If QMOM
        if(!nodes[0].extended())
        {
            forAll(quadrature_.moments()[1], celli)
            {
                forAll(nodes, pNodeI)
                {
                    volScalarNode& node = nodes[pNodeI];
                    scalar& abscissa = node.primaryAbscissae()[sizeIndex][celli];

                    if (abscissa > minAbscissa)
                    {
                        abscissa = abscissa - drift[celli] / abscissa * globalDt;
                    }

                    if (abscissa <= minAbscissa)
                    {
                        abscissa = minAbscissa;
                    }
                }
            }
            quadrature_.updateMoments();
            quadrature_.updateQuadrature();
            
            return;
        }

        Info<< "Perfoam EQMOM evaporation: " << endl;
        //Info<< "Original moments: " << quadrature_.moments() << endl; 

        forAll(quadrature_.moments()[1], celli)
        {
            forAll(nodes, pNodeI)
            {
                volScalarNode& node = nodes[pNodeI];
                forAll(node.secondaryWeights()[sizeIndex], sNodei)
                {
                    scalar& abscissa = node.secondaryAbscissae()[sizeIndex][sNodei][celli];

                    //Info<< "Node" << pNodeI << ", SNode" << sNodei << endl; 
                    //Info<< "Original abscissa: " << abscissa << endl;

                    if (abscissa > minAbscissa)
                    {
                        abscissa = abscissa - drift[celli] / abscissa * globalDt;
                    }
                    
                    //Info<< "Evaporated abscissa: " << abscissa << endl;

                    if (abscissa <= minAbscissa)
                    {
                        abscissa = minAbscissa;
                    }

                    //Info<< "Final abscissa: " << abscissa << endl;
                    //Info<< "Final node: " << node.secondaryAbscissae()[sizeIndex][sNodei][celli] << endl;
                }
            }
        }
        quadrature_.updateMoments();
        quadrature_.updateQuadrature();

        //Info<< "Final moments: " << quadrature_.moments() << endl;
        //Info<< "Final quarature: " << quadrature_.nodes() << endl;  
        
        return;
    }
    /*
    if (evaporation_)
    {
        PtrList<volScalarNode>& nodes = quadrature_.nodes();
        label sizeIndex = nodes[0].sizeIndex();
        bool lengthBased = nodes[0].lengthBased();
        //- GlobalDt
        scalar globalDt = quadrature_.moments()[0].mesh().time().deltaT().value();
        //- Minimum abscissa
        scalar minAbscissa = evaporationModel_->minAbscissa();

        //- If QMOM
        if(!nodes[0].extended())
        {
            forAll(quadrature_.moments()[1], celli)
            {
                forAll(nodes, pNodeI)
                {
                    volScalarNode& node = nodes[pNodeI];
                    scalar& abscissa = node.primaryAbscissae()[sizeIndex][celli];

                    abscissa = abscissa - evaporationModel_->Kg(abscissa, celli, lengthBased) * globalDt;

                    if (abscissa <= minAbscissa)
                    {
                        abscissa = minAbscissa;
                    }
                }
            }
            quadrature_.updateMoments();
            quadrature_.updateQuadrature();
            
            return;
        }

        Info<< "Perfoam EQMOM evaporation: " << endl;
        //Info<< "Original moments: " << quadrature_.moments() << endl; 

        forAll(quadrature_.moments()[1], celli)
        {
            forAll(nodes, pNodeI)
            {
                volScalarNode& node = nodes[pNodeI];
                forAll(node.secondaryWeights()[sizeIndex], sNodei)
                {
                    scalar& abscissa = node.secondaryAbscissae()[sizeIndex][sNodei][celli];

                    //Info<< "Node" << pNodeI << ", SNode" << sNodei << endl; 
                    //Info<< "Original abscissa: " << abscissa << endl;

                    abscissa = abscissa - evaporationModel_->Kg(abscissa, celli, lengthBased) * globalDt;

                    //Info<< "Evaporated abscissa: " << abscissa << endl;

                    if (abscissa <= minAbscissa)
                    {
                        abscissa = minAbscissa;
                    }

                    //Info<< "Final abscissa: " << abscissa << endl;
                    //Info<< "Final node: " << node.secondaryAbscissae()[sizeIndex][sNodei][celli] << endl;
                }
            }
        }
        quadrature_.updateMoments();
        quadrature_.updateQuadrature();

        //Info<< "Final moments: " << quadrature_.moments() << endl;
        //Info<< "Final quarature: " << quadrature_.nodes() << endl;  
        
        return;
    }
    */
}


void
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::solve()
{
    univariateMyDivPDFTransportModel::solve();
}


const Foam::scalarQuadratureApproximation& Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::quadrature() const
{
    return univariateMyDivPDFTransportModel::quadrature();
}


bool 
Foam::PDFTransportModels::populationBalanceModels::univariateEScaleMyDivPopulationBalance
::readIfModified()
{
    odeType::read
    (
        populationBalanceProperties_.subDict(type() + "Coeffs")
    );
    
    return true;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "JohnsonJacksonParticleSlipFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        JohnsonJacksonParticleSlipFvPatchVectorField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
JohnsonJacksonParticleSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
    specularityCoefficient_("specularityCoefficient", dimless, 0)
{}


Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
JohnsonJacksonParticleSlipFvPatchVectorField
(
    const JohnsonJacksonParticleSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(ptf, p, iF, mapper),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
JohnsonJacksonParticleSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
    specularityCoefficient_
    (
        "specularityCoefficient",
        dimless,
        dict
    ),
    internalFrictionAngle_
    (
        "internalFrictionAngle",
        dimless,
        dict
    )
{
    if
    (
        (specularityCoefficient_.value() < 0)
     || (specularityCoefficient_.value() > 1)
    )
    {
        FatalErrorInFunction
            << "The specularity coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    fvPatchVectorField::operator=
    (
        vectorField("value", dict, p.size())
    );
}


Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
JohnsonJacksonParticleSlipFvPatchVectorField
(
    const JohnsonJacksonParticleSlipFvPatchVectorField& ptf
)
:
    partialSlipFvPatchVectorField(ptf),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
JohnsonJacksonParticleSlipFvPatchVectorField
(
    const JohnsonJacksonParticleSlipFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(ptf, iF),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    partialSlipFvPatchVectorField::autoMap(m);
}


void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    partialSlipFvPatchVectorField::rmap(ptf, addr);
}


void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // lookup the fluid model and the phase
    const twoPhaseSystem& fluid = db().lookupObject<twoPhaseSystem>
    (
        "phaseProperties"
    );

    const phaseModel& phased
    (
        fluid.phase1().name() == internalField().group()
      ? fluid.phase1()
      : fluid.phase2()
    );

    // lookup all the fields on this patch
    const fvPatchScalarField& alpha
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            phased.volScalarField::name()
        )
    );

    const scalarField nu
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("nut", phased.name())
        )
    );

    word ThetaName(IOobject::groupName("Theta", phased.name()));

    const fvPatchScalarField& Theta
    (
        db().foundObject<volScalarField>(ThetaName)
      ? patch().lookupPatchField<volScalarField, scalar>(ThetaName)
      : alpha
    );

    scalarField c(alpha.size(), Zero);

    if
    (
        db().foundObject<volScalarField>
        (
            IOobject::groupName("h2Fn", phased.name())
        )
    )
    {
        const scalarField& h2Fn
        (
            patch().lookupPatchField<volScalarField, scalar>
            (
                IOobject::groupName("h2Fn", phased.name())
            )
        );

        const fvPatchScalarField& PsFric
        (
            patch().lookupPatchField<volScalarField, scalar>
            (
                IOobject::groupName("PsFric", phased.name())
            )
        );

        scalarField Vw(constant::mathematical::pi/6.0*sqrt(3.0*Theta));

        // calculate the slip value fraction
        scalarField c
        (
            (
                h2Fn*specularityCoefficient_.value()*Vw
              + PsFric*tan(internalFrictionAngle_.value())
               /max(alpha*mag(patchInternalField()), scalar(1e-4))
            )/max(nu, SMALL)
        );

    }
    else
    {
        const fvPatchScalarField& gs0
        (
            patch().lookupPatchField<volScalarField, scalar>
            (
                IOobject::groupName("gs0", phased.name())
            )
        );

        const scalarField nuFric
        (
            patch().lookupPatchField<volScalarField, scalar>
            (
                IOobject::groupName("nuFric", phased.name())
            )
        );

        // lookup the maximum allowed volume fraction
        dimensionedScalar alphaMax
        (
            "alphaMax",
            dimless,
            db()
            .lookupObject<IOdictionary>
            (
                IOobject::groupName("turbulenceProperties", phased.name())
            )
            .subDict("RAS")
            .subDict("kineticTheoryCoeffs")
        );

        scalarField c
        (
            constant::mathematical::pi
           *alpha
           *gs0
           *specularityCoefficient_.value()
           *sqrt(3.0*Theta)
           /max(6.0*(nu - nuFric)*alphaMax.value(), SMALL)
        );

        // calculate the slip value fraction
        c /= max(6.0*(nu - nuFric)*alphaMax.value(), SMALL);
    }

    this->valueFraction() = c/(c + patch().deltaCoeffs());

    partialSlipFvPatchVectorField::updateCoeffs();
}


void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("specularityCoefficient", specularityCoefficient_);
    os.writeEntry("internalFrictionAngle", internalFrictionAngle_);
    writeEntry("value", os);
}


// ************************************************************************* //

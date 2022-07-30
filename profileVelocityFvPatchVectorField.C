/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "Field.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvPatchFieldMapper.H"
#include "meshCellZonesFwd.H"
#include "profileVelocityFvPatchVectorField.H"
#include "scalar.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "wallDist.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::profileVelocityFvPatchVectorField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profileVelocityFvPatchVectorField::profileVelocityFvPatchVectorField(
    const fvPatch &p, const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchVectorField(p, iF), umax_(Zero), profile_() {}

Foam::profileVelocityFvPatchVectorField::profileVelocityFvPatchVectorField(
    const fvPatch &p, const DimensionedField<vector, volMesh> &iF,
    const dictionary &dict)
    : fixedValueFvPatchVectorField(p, iF), umax_(dict.lookup<vector>("umax")),
      profile_(Function1<scalar>::New("profile", dict)) {

    fixedValueFvPatchVectorField::evaluate();

  /*
  // Initialise with the value entry if evaluation is not possible
  fvPatchVectorField::operator=
  (
      vectorField("value", dict, p.size())
  );
  */
}

Foam::profileVelocityFvPatchVectorField::profileVelocityFvPatchVectorField(
    const profileVelocityFvPatchVectorField &ptf, const fvPatch &p,
    const DimensionedField<vector, volMesh> &iF,
    const fvPatchFieldMapper &mapper)
    : fixedValueFvPatchVectorField(ptf, p, iF, mapper), umax_(ptf.umax_),
      profile_(ptf.profile_, false) {}

Foam::profileVelocityFvPatchVectorField::profileVelocityFvPatchVectorField(
    const profileVelocityFvPatchVectorField &ptf,
    const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchVectorField(ptf, iF), umax_(ptf.umax_),
      profile_(ptf.profile_, false) {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::profileVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::profileVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);


}


void Foam::profileVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh mesh = patch().boundaryMesh().mesh();

    const volScalarField y(wallDist::New(mesh).y());

    Field<scalar> yp=patch().patchInternalField(y);
    Field<scalar> profiley=profile_->value(yp);

    profiley/=gMax(profiley);

    fixedValueFvPatchField::operator==
    (
        umax_*profiley
    );


    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::profileVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "umax", umax_);
    writeEntry(os, profile_());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        profileVelocityFvPatchVectorField
    );
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::incompressible::RASModels::kappatWallFunction4SEFoamFvPatchScalarField

Description
    Boundary condition for turbulent thermal diffusivity when using wall
    functions
    - replicates OpenFOAM v1.5 (and earlier) behaviour
    - Turbulent Prandtl number defaults to 0.85 if unspecified

SourceFiles
    kappatWallFunction4SEFoamFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef kappatWallFunction4SEFoamFvPatchScalarField_H
#define kappatWallFunction4SEFoamFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
            Class kappatWallFunction4SEFoamFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class kappatWallFunction4SEFoamFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of turbulent viscosity field
        word nutName_;

        //- Turbulent Prandtl number
        scalar Prt_;


protected:

    // Protected Member Functions

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("kappatWallFunction4SEFoam");


    // Constructors

        //- Construct from patch and internal field
        kappatWallFunction4SEFoamFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        kappatWallFunction4SEFoamFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  kappatWallFunction4SEFoamFvPatchScalarField
        //  onto a new patch
        kappatWallFunction4SEFoamFvPatchScalarField
        (
            const kappatWallFunction4SEFoamFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        kappatWallFunction4SEFoamFvPatchScalarField
        (
            const kappatWallFunction4SEFoamFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new kappatWallFunction4SEFoamFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        kappatWallFunction4SEFoamFvPatchScalarField
        (
            const kappatWallFunction4SEFoamFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new kappatWallFunction4SEFoamFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

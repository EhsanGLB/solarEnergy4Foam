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
    Foam::timeVaryingMixedRobin4SEFoam

Description

    @verbatim
        wall
        {
            type                    timeVaryingMixedRobin4SEFoam;
            Kappa                   kappa;
            absorptivity            absorptivity;
            emissivity              emissivity;
            qs
            {
                "fileName|file" "$FOAM_CASE/constant/qsVsTime";
                outOfBounds     clamp;           // (error|warn|clamp|repeat)
            }
            ho
            {
                "fileName|file" "$FOAM_CASE/constant/hovsTime";
                outOfBounds     clamp;
            }
            To
            {
                "fileName|file" "$FOAM_CASE/constant/TovsTime";
                outOfBounds     clamp;
            }
            Tsky
            {
                "fileName|file" "$FOAM_CASE/constant/TskyvsTime";
                outOfBounds     clamp;
            }
        }
    @endverbatim

SourceFiles
    timeVaryingMixedRobin4SEFoam.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMixedRobin4SEFoam_H
#define timeVaryingMixedRobin4SEFoam_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformRobinFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingMixedRobin4SEFoam
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- Name of the thermal conductivity field
        word kappaName_;

        //- absorptivity
        scalar absorptivity_;

        //- emissivity
        scalar emissivity_;

        //- Table of time vs heat flux, including the bounding treatment
        interpolationTable<scalar> timeSeriesqs_;

        //- Table of time vs convective heat transfer coefficient, including the bounding treatment
        interpolationTable<scalar> timeSeriesho_;

        //- Table of time vs temperature, including the bounding treatment
        interpolationTable<scalar> timeSeriesTo_;

        //- Table of time vs sky temperature, including the bounding treatment
        interpolationTable<scalar> timeSeriesTsky_;


public:

    //- Runtime type information
    TypeName("timeVaryingMixedRobin4SEFoam");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMixedRobin4SEFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMixedRobin4SEFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingMixedRobin4SEFoam onto
        //  a new patch
        timeVaryingMixedRobin4SEFoam
        (
            const timeVaryingMixedRobin4SEFoam&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMixedRobin4SEFoam
        (
            const timeVaryingMixedRobin4SEFoam&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new timeVaryingMixedRobin4SEFoam(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingMixedRobin4SEFoam
        (
            const timeVaryingMixedRobin4SEFoam&,
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
                new timeVaryingMixedRobin4SEFoam(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

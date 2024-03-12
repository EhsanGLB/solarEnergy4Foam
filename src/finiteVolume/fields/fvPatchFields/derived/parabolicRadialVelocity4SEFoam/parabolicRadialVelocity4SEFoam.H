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
    parabolicRadialVelocity4SEFoam

Description
    Boundary condition specifies a parabolic radial velocity inlet profile
    , given mean velocity value, R of inlet,
    and flow direction.

    @verbatim
        inlet
        {
            type                    parabolicRadialVelocity4SEFoam;
            Vm                      (12.3);
            R                       0.25;
            longDir                 (1 0 0);
        }
    @endverbatim

SourceFiles
    parabolicRadialVelocity4SEFoam.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef parabolicRadialVelocity4SEFoam_H
#define parabolicRadialVelocity4SEFoam_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class parabolicRadialVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class parabolicRadialVelocity4SEFoam
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- mean velocity magnitude depends on time
        scalarList Vm_;

        //- Radius
        scalar R_;

        //- Flow direction
        vector longDir_;


public:

    //- Runtime type information
    TypeName("parabolicRadialVelocity4SEFoam");


    // Constructors

        //- Construct from patch and internal field
        parabolicRadialVelocity4SEFoam
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicRadialVelocity4SEFoam
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given parabolicRadialVelocity4SEFoam
        //  onto a new patch
        parabolicRadialVelocity4SEFoam
        (
            const parabolicRadialVelocity4SEFoam&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicRadialVelocity4SEFoam(*this)
            );
        }

        //- Construct as copy setting internal field reference
        parabolicRadialVelocity4SEFoam
        (
            const parabolicRadialVelocity4SEFoam&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicRadialVelocity4SEFoam(*this, iF)
            );
        }


    // Member functions

        //- Calculate the qo(time)
        virtual scalar VmFunction(scalar time);

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

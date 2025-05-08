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
    Foam::parabolicCollector4SEFoam

Description
    Foam::parabolicCollector4SEFoam

    @verbatim
        hotWall
        {
            type		parabolicCollector4SEFoam;
            Kappa		kappaEff;
            range		(0 68 114 180);
            LCR1		(26.65    -0.01496    0.04643    -2.258e-3    4.175e-5    -2.588e-7);
            LCR2		(-181.1    9.099    -0.1121    4.071e-4    0.0    0.0);
            LCR3		(9.668    -0.1821    1.125e-3    -2.157e-6    0.0    0.0);
            glassRadProp	(0.01    0.02    0.97);//(Absorptivity  emissivity transmissivity)
            tubeRadProp		(0.96    0.04    0.0);//(Absorptivity  emissivity transmissivity)
            dgi			0.12;
            dgo			0.125;
            dti			0.076;
            dto			0.08;
            len			5.0;
            longDir		(1 0 0);
            radDir		(0 -1 0);
            qo			(q0 q1 q2);
            ho			(h0 h1 h2);
            To			(T0 T1 T2);
        }
    @endverbatim

SourceFiles
    parabolicCollector4SEFoam.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef parabolicCollector4SEFoams_H
#define parabolicCollector4SEFoams_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformRobinPTCFvPatch Declaration
\*---------------------------------------------------------------------------*/

class parabolicCollector4SEFoam
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word TName_;

        //- Name of the thermal conductivity field
        word kappaName_;

        //- range
        scalarList range_;

        //- LCRs
        scalarList LCR1_;
        scalarList LCR2_;
        scalarList LCR3_;

        //- Radiative properties of glass
        scalarList glassRadProp_;

        //- Radiative properties of tube
        scalarList tubeRadProp_;

        //- Inner glass diameter
        scalar dgi_;

        //- Outer glass diameter
        scalar dgo_;

        //- Inner tube diameter
        scalar dti_;

        //- Outer tube diameter
        scalar dto_;

        //- PTC length
        scalar len_;

        //- longitude direction
        vector longDir_;

        //- radiative direction
        vector radDir_;

        //- Flux depends on time
        scalarList qo_;

        //- convective heat transfer
        scalarList ho_;

        //- Out temperature
        scalarList To_;


public:

    //- Runtime type information
    TypeName("parabolicCollector4SEFoam");


    // Constructors

        //- Construct from patch and internal field
        parabolicCollector4SEFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicCollector4SEFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given parabolicCollector4SEFoam onto
        //  a new patch
        parabolicCollector4SEFoam
        (
            const parabolicCollector4SEFoam&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        parabolicCollector4SEFoam
        (
            const parabolicCollector4SEFoam&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new parabolicCollector4SEFoam(*this)
            );
        }

        //- Construct as copy setting internal field reference
        parabolicCollector4SEFoam
        (
            const parabolicCollector4SEFoam&,
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
                new parabolicCollector4SEFoam(*this, iF)
            );
        }


    // Member functions

        //- Calculate the qo(time)
        virtual scalar qoFunction(scalar time);

        //- Calculate the ho(time)
        virtual scalar hoFunction(scalar time);

        //- Calculate the To(time)
        virtual scalar ToFunction(scalar time);

        //- Calculate the LCR(teta)
        virtual scalar LCRFunction(scalar teta);

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

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

Namespace
    Foam::fluids

Description
    A namespace for various incompressible fluid implementations.

Class
    Foam::fluid

Description
    An class for incompressible fluid

SourceFiles
    fluid.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef fluid_H
#define fluid_H

#include "IOdictionary.H"
#include "volFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"

#include "material.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fluid Declaration
\*---------------------------------------------------------------------------*/

class fluid
{
protected:

    // Protected data
        dictionary fluidPropertiesDict_;
        const volVectorField& U_;
        const volScalarField& p_;
        const volScalarField& T_;

        //- The material properties
        autoPtr<material> materialPtr_;


    // Private Member Functions

        //- Disallow copy construct
        fluid(const fluid&);

        //- Disallow default bitwise assignment
        void operator=(const fluid&);


public:

    // Constructors

        //- Construct from components
        fluid
        (
            const dictionary& fluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T
        );

   // Selectors

        //- Select construct from dictionary
        static autoPtr<fluid> New
        (
            const dictionary& fluidPropertiesDict,
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField& T
        );


    // Destructor

        virtual ~fluid()
        {}


    // Member Functions

        //- Return the density [kg/m^3]
        virtual const volScalarField rho() const;

        //- Return the thermal conductivity [W/(m K)]
        virtual const volScalarField kappa() const;

        //- Return the specific heat capacity [J/(kg K)]
        virtual const volScalarField Cp() const;

        //- Return the dynamic viscosity [Pa s]
        virtual const volScalarField mu() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

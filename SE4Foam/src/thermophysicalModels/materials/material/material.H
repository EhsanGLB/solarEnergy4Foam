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
    Foam::material

Description
    The thermophysical properties of a material

SourceFiles
    material.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef material_H
#define material_H

#include "scalar.H"
#include "IOstreams.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

#include "expPolyfunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class material Declaration
\*---------------------------------------------------------------------------*/

class material
{
    // Private data

        //- Molecular weight [kg/kmol]
        scalar W_;


public:

    TypeName("material");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            material,
            ,
            (),
            ()
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            material,
            Istream,
            (Istream& is),
            (is)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            material,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );

    // Constructors

        //- Construct from components
        material(){}

        //material(Istream& is){}
        material
        (
            scalar W
        )
        :
            W_(W)
        {}

        //- Construct from Istream
        material(Istream& is)
        :
            W_(readScalar(is))
        {}

        //- Return a pointer to a new material created from Istream
        static autoPtr<material> New(Istream& is);

        //- Return a pointer to a new material created from dictionary
        static autoPtr<material> New(const dictionary& dict);


    //- Destructor
    virtual ~material()
    {}


    // Member Functions

        // Physical constants
            //- Molecular weight [kg/kmol]
            inline scalar W() const;

        // Physical property pure virtual functions

            //- material density [kg/m^3]
            virtual scalar rho(scalar p, scalar T) const = 0;

            //- material thermal conductivity [W/(m K)]
            virtual scalar kappa(scalar p, scalar T) const = 0;

            //- material heat capacity [J/(kg K)]
            virtual scalar Cp(scalar p, scalar T) const = 0;

            //- material dynamic viscosity [Pa s]
            virtual scalar mu(scalar p, scalar T) const = 0;


    // I-O

        //- Write the function coefficients
        virtual void writeData(Ostream& os) const
        {
            os  << W_;
        }

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const material& m)
        {
            m.writeData(os);
            return os;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "materialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

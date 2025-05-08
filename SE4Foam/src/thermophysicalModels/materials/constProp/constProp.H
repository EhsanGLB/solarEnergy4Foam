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
    Foam::constProp

Description
    constProp

SourceFiles
    constProp.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef constProp_H
#define constProp_H

#include "material.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace materials
{

/*---------------------------------------------------------------------------*\
                           Class constProp Declaration
\*---------------------------------------------------------------------------*/

class constProp
:
    public material
{
    // Private data

        scalar rho_;
        scalar kappa_;
        scalar Cp_;
        scalar mu_;


public:

    //- Runtime type information
    TypeName("constProp");


    // Constructors

        //- Construct from Istream
        constProp(Istream& is);


    // Member Functions

            //- material density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- material thermal conductivity  [W/(m K)]
            inline scalar kappa(scalar p, scalar T) const;

            //- material heat capacity [J/(kg K)]
            inline scalar Cp(scalar p, scalar T) const;

            //- material dynamic viscosity [Pa s]
            inline scalar mu(scalar p, scalar T) const;
};

/*
    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            material::writeData(os); os << nl;
            rho.writeData(os); os << nl;
            os  << rho_ << token::SPACE
                << kappa_ << token::SPACE
                << Cp_ << token::SPACE
                << mu_;
        }


        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const constProp& m)
        {
            m.writeData(os);
            return os;
        }
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace materials
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constPropI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

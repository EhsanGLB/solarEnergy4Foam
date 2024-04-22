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
    Foam::funcTProp

Description
    funcTProp

SourceFiles
    funcTProp.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef funcTProp_H
#define funcTProp_H

#include "material.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace materials
{

/*---------------------------------------------------------------------------*\
                           Class funcTProp Declaration
\*---------------------------------------------------------------------------*/

class funcTProp
:
    public material
{
    // Private data
        dictionary dict_;
        expPolyfunc rho_;
        expPolyfunc kappa_;
        expPolyfunc Cp_;
        expPolyfunc mu_;


public:

    //- Runtime type information
    TypeName("funcTProp");


    // Constructors

        //- Construct from Istream
        funcTProp(const dictionary& dict);


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
            rho_.writeData(os); os << nl;
            kappa_.writeData(os); os << nl;
            Cp_.writeData(os); os << nl;
            mu_.writeData(os); os << nl;
        }


        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const funcTProp& m)
        {
            m.writeData(os);
            return os;
        }
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace materials
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "funcTPropI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

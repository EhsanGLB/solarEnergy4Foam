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
    Foam::expPolyfunc

Description

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef expPolyfunc_H
#define expPolyfunc_H

#include "thermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class expPolyfunc Declaration
\*---------------------------------------------------------------------------*/

class expPolyfunc
:
    public thermophysicalFunction
{
    // Private data

        // NSRDS function 100 coefficients
        scalar a0_, a1_, a2_, a3_, a4_, a5_, a6_, b_, c_;


public:

    //- Runtime type information
    TypeName("expPolyfunc");


    // Constructors

        //- Construct from components
        expPolyfunc(scalar a0_, scalar a1_, scalar a2_, scalar a3_, scalar a4_, scalar a5_, scalar a6_, scalar b_, scalar c_)
        :
            a0_(a0_),
            a1_(a1_),
            a2_(a2_),
            a3_(a3_),
            a4_(a4_),
            a5_(a5_),
            a6_(a6_),
            b_(b_),
            c_(c_)
        {}

        //- Construct from Istream
        expPolyfunc(Istream& is)
        :
            a0_(readScalar(is)),
            a1_(readScalar(is)),
            a2_(readScalar(is)),
            a3_(readScalar(is)),
            a4_(readScalar(is)),
            a5_(readScalar(is)),
            a6_(readScalar(is)),
            b_(readScalar(is)),
            c_(readScalar(is))
        {}


    // Member Functions

        //- Evaluate the function and return the result
        scalar f(scalar p, scalar T) const
        {
            return ( a0_ + a1_ * T + a2_ * pow(T, 2) + a3_ * pow(T, 3) + a4_ * pow(T, 4) + a5_ * pow(T, 5) + a6_ * pow(T, 6) ) * exp(b_ * T) + c_;
        }


        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            os  << a0_ << token::SPACE
                << a1_ << token::SPACE
                << a2_ << token::SPACE
                << a3_ << token::SPACE
                << a4_ << token::SPACE
                << a5_ << token::SPACE
                << a6_ << token::SPACE
                << b_ << token::SPACE
                << c_;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const expPolyfunc& g)
        {
            g.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

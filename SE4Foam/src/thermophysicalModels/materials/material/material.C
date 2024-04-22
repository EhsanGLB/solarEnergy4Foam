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

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "material.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(material, 0);
    defineRunTimeSelectionTable(material,);
    defineRunTimeSelectionTable(material, Istream);
    defineRunTimeSelectionTable(material, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<material> material::New(Istream& is)
{
    if (debug)
    {
        Info<< "material::New(Istream&) : "
            << "constructing material"
            << endl;
    }
    word materialType(is);
    word coeffs(is);

    if (coeffs == "defaultCoeffs")
    {
        ConstructorTable::iterator cstrIter =
            ConstructorTablePtr_->find(materialType);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorIn("material::New(Istream&)")
                << "Unknown material type " << materialType
                << nl << nl
                << "Valid material types are:" << nl
                << ConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<material>(cstrIter()());
    }
    else if (coeffs == "coeffs")
    {
        IstreamConstructorTable::iterator cstrIter =
            IstreamConstructorTablePtr_->find(materialType);

        if (cstrIter == IstreamConstructorTablePtr_->end())
        {
            FatalErrorIn("material::New(Istream&)")
                << "Unknown material type " << materialType
                << endl << endl
                << "Valid material types are:" << nl
                << IstreamConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<material>(cstrIter()(is));
    }
    else
    {
        FatalErrorIn("material::New(Istream&)")
            << "material type " << materialType
            << ", option " << coeffs << " given"
            << ", should be coeffs or defaultCoeffs"
            << abort(FatalError);

        return autoPtr<material>(nullptr);
    }
}


autoPtr<material> material::New(const dictionary& dict)
{

    if (debug)
    {
        Info<< "material::New(dictionary&) : "
            << "constructing material"
            << endl;
    }

    Istream& matInitIs(dict.lookup("type"));
    word materialType(matInitIs);

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(materialType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
            FatalErrorIn("material::New(dictionary&)")
                << "Unknown material type " << materialType
                << nl << nl
                << "Valid material types are:" << nl
                << ConstructorTablePtr_->sortedToc()
                << abort(FatalError);
    }

    return autoPtr<material>
        (cstrIter()(dict));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

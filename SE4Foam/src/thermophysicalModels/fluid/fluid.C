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

#include "volFields.H"

#include "fluid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<fluid> fluid::New
(
    const dictionary& fluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
{
    return autoPtr<fluid>(new fluid(fluidPropertiesDict, U, p, T));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fluid::fluid
(
    const dictionary& fluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    fluidPropertiesDict_(fluidPropertiesDict),
    U_(U),
    p_(p),
    T_(T)
{
//    materialPtr_(material::New(fluidPropertiesDict))

    const dictionary* subDictPtr = fluidPropertiesDict.subDictPtr("fluidFuncTProp");

    if (subDictPtr)
    {
        materialPtr_ = material::New(subDictPtr);
    }

    else
    {
        materialPtr_ = material::New(fluidPropertiesDict.lookup("fluidConstProp"));
    }

}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField fluid::rho() const
{
    volScalarField rho_(IOobject("rhof", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("rhof", dimDensity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        rho_.internalField()[i] = materialPtr_->rho(pi_, Ti_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            rho_.boundaryField()[i][j] = materialPtr_->rho(pij_, Tij_);
	}
    }

    return rho_;
}


const volScalarField fluid::kappa() const
{
    volScalarField kappa_(IOobject("kappaf", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaf", dimThermalConductivity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        kappa_.internalField()[i] = materialPtr_->kappa(pi_, Ti_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            kappa_.boundaryField()[i][j] = materialPtr_->kappa(pij_, Tij_);
	}
    }

    return kappa_;
}


const volScalarField fluid::Cp() const
{
    volScalarField Cp_(IOobject("Cpf", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("Cpf", dimSpecificHeatCapacity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Cp_.internalField()[i] = materialPtr_->Cp(pi_, Ti_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Cp_.boundaryField()[i][j] = materialPtr_->Cp(pij_, Tij_);
	}
    }

    return Cp_;
}


const volScalarField fluid::mu() const
{
    volScalarField mu_(IOobject("muf", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muf", dimMass/dimLength/dimTime, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        mu_.internalField()[i] = materialPtr_->mu(pi_, Ti_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            mu_.boundaryField()[i][j] = materialPtr_->mu(pij_, Tij_);
	}
    }

    return mu_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

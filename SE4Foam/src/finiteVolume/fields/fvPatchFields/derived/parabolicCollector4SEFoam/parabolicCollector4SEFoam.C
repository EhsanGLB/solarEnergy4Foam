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

#include "parabolicCollector4SEFoam.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "scalar.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parabolicCollector4SEFoam::parabolicCollector4SEFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    kappaName_("undefined-Kappa"),
    range_(0),
    LCR1_(0),
    LCR2_(0),
    LCR3_(0),
    glassRadProp_(0),
    tubeRadProp_(0),
    dgi_(0.0),
    dgo_(0.0),
    dti_(0.0),
    dto_(0.0),
    len_(0.0),
    longDir_(0, 0, 0),
    radDir_(0, 0, 0),
    qo_(0.0),
    ho_(0.0),
    To_(0.0)
{}


Foam::parabolicCollector4SEFoam::parabolicCollector4SEFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    kappaName_(dict.lookup("Kappa")),
    range_(),
    LCR1_(),
    LCR2_(),
    LCR3_(),
    glassRadProp_(),
    tubeRadProp_(),
    dgi_(readScalar(dict.lookup("dgi"))),
    dgo_(readScalar(dict.lookup("dgo"))),
    dti_(readScalar(dict.lookup("dti"))),
    dto_(readScalar(dict.lookup("dto"))),
    len_(readScalar(dict.lookup("len"))),
    longDir_(dict.lookup("longDir")),
    radDir_(dict.lookup("radDir")),
    qo_(),
    ho_(),
    To_()
{
    Istream& isrange_ = dict.lookup("range");
    isrange_.format(IOstream::ASCII);
    isrange_ >> range_;

    Istream& isLCR1_ = dict.lookup("LCR1");
    isLCR1_.format(IOstream::ASCII);
    isLCR1_ >> LCR1_;

    Istream& isLCR2_ = dict.lookup("LCR2");
    isLCR2_.format(IOstream::ASCII);
    isLCR2_ >> LCR2_;

    Istream& isLCR3_ = dict.lookup("LCR3");
    isLCR3_.format(IOstream::ASCII);
    isLCR3_ >> LCR3_;

    Istream& isglassRadProp_ = dict.lookup("glassRadProp");
    isglassRadProp_.format(IOstream::ASCII);
    isglassRadProp_ >> glassRadProp_;

    Istream& istubeRadProp_= dict.lookup("tubeRadProp");
    istubeRadProp_.format(IOstream::ASCII);
    istubeRadProp_ >> tubeRadProp_;

    Istream& isqo_ = dict.lookup("qo");
    isqo_.format(IOstream::ASCII);
    isqo_ >> qo_;

    Istream& isho_ = dict.lookup("ho");
    isho_.format(IOstream::ASCII);
    isho_ >> ho_;

    Istream& isTo_ = dict.lookup("To");
    isTo_.format(IOstream::ASCII);
    isTo_ >> To_;

    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::parabolicCollector4SEFoam::parabolicCollector4SEFoam
(
    const parabolicCollector4SEFoam& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    range_(ptf.range_),
    LCR1_(ptf.LCR1_),
    LCR2_(ptf.LCR2_),
    LCR3_(ptf.LCR3_),
    glassRadProp_(ptf.glassRadProp_),
    tubeRadProp_(ptf.tubeRadProp_),
    dgi_(ptf.dgi_),
    dgo_(ptf.dgo_),
    dti_(ptf.dti_),
    dto_(ptf.dto_),
    len_(ptf.len_),
    longDir_(ptf.longDir_),
    radDir_(ptf.radDir_),
    qo_(ptf.qo_),
    ho_(ptf.ho_),
    To_(ptf.To_)
{}


Foam::parabolicCollector4SEFoam::parabolicCollector4SEFoam
(
    const parabolicCollector4SEFoam& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    range_(ptf.range_),
    LCR1_(ptf.LCR1_),
    LCR2_(ptf.LCR2_),
    LCR3_(ptf.LCR3_),
    glassRadProp_(ptf.glassRadProp_),
    tubeRadProp_(ptf.tubeRadProp_),
    dgi_(ptf.dgi_),
    dgo_(ptf.dgo_),
    dti_(ptf.dti_),
    dto_(ptf.dto_),
    len_(ptf.len_),
    longDir_(ptf.longDir_),
    radDir_(ptf.radDir_),
    qo_(ptf.qo_),
    ho_(ptf.ho_),
    To_(ptf.To_)
{}


Foam::parabolicCollector4SEFoam::parabolicCollector4SEFoam
(
    const parabolicCollector4SEFoam& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    range_(ptf.range_),
    LCR1_(ptf.LCR1_),
    LCR2_(ptf.LCR2_),
    LCR3_(ptf.LCR3_),
    glassRadProp_(ptf.glassRadProp_),
    tubeRadProp_(ptf.tubeRadProp_),
    dgi_(ptf.dgi_),
    dgo_(ptf.dgo_),
    dti_(ptf.dti_),
    dto_(ptf.dto_),
    len_(ptf.len_),
    longDir_(ptf.longDir_),
    radDir_(ptf.radDir_),
    qo_(ptf.qo_),
    ho_(ptf.ho_),
    To_(ptf.To_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::parabolicCollector4SEFoam::qoFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < qo_.size() )
    {
        f_ += qo_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}


Foam::scalar Foam::parabolicCollector4SEFoam::hoFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < ho_.size() )
    {
        f_ += ho_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}


Foam::scalar Foam::parabolicCollector4SEFoam::ToFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < To_.size() )
    {
        f_ += To_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}


Foam::scalar Foam::parabolicCollector4SEFoam::LCRFunction(scalar teta)
{
////range coefficients
    scalar teta0_(range_[0]);
    scalar teta1_(range_[1]);
    scalar teta2_(range_[2]);
    scalar teta3_(range_[3]);

    scalar f_ = 0.0;

    if ( teta0_ <= teta && teta < teta1_)
    {
        int n1_ = 0;

        while( n1_ < LCR1_.size() )
        {
            f_ += LCR1_[n1_]*pow(teta, n1_);
            n1_ += 1;
        }
    }

    else if ( teta1_ <= teta && teta < teta2_)
    {
        int n2_ = 0;

        while( n2_ < LCR2_.size() )
        {
            f_ += LCR2_[n2_]*pow(teta, n2_);
            n2_ += 1;
        }
    }

    else if ( teta2_ <= teta && teta < teta3_)
    {
        int n3_ = 0;

        while( n3_ < LCR3_.size() )
        {
            f_ += LCR3_[n3_]*pow(teta, n3_);
            n3_ += 1;
        }
    }

    return f_;
}


void Foam::parabolicCollector4SEFoam::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //- get the currernt time and time step
    scalar t_ = this->db().time().value();

////radiative properties of glass
    scalar ag_(glassRadProp_[0]);
    scalar eg_(glassRadProp_[1]);
    scalar tg_(glassRadProp_[2]);

////radiative properties of tube
    scalar at_(tubeRadProp_[0]);
    scalar et_(tubeRadProp_[1]);
    //scalar tt_(tubeRadProp_[2]);

////constant parameters
    const scalar sigmaSB_(5.67e-8);

////fixed values
    //scalar Agi_(mathematicalConstant::pi*len_*dgi_);
    scalar Ago_(mathematicalConstant::pi*len_*dgo_);
    scalar Ati_(mathematicalConstant::pi*len_*dti_);
    scalar Ato_(mathematicalConstant::pi*len_*dto_);
    scalar skyTemp_(0.0522*pow(ToFunction(t_), 1.5));

//// To calculate degree field
    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    boundBox bb_(patch().patch().localPoints(), true);
    vector ctr_ = 0.5*(bb_.max() + bb_.min());
    const vectorField& patchCent_ = patch().Cf();
    scalarField degree_(patchCent_.size(), 0.0);

    forAll(degree_, i)
    {
        vector op_ = (patchCent_[i]-ctr_)^(longDir_);
        degree_[i] = acos((op_/mag(op_)) & (radDir_/mag(radDir_))) * (180/mathematicalConstant::pi);
    }

//// To calculate LCR average and tube average temperature
    scalar LCRAvg_(0.0);
    for(int i=0; i<180; i++)
    {
        LCRAvg_ += LCRFunction(scalar(i));
    }
    LCRAvg_ /= 180;

    scalar Tt_(0.0);

    forAll(Tp, i)
    {
        Tt_ += Tp[i] * patch().magSf()[i];
    }

    Tt_ /= Ati_;

//// To calculate heat loss from tube
    scalar x1 = 10;
    scalar x2 = 1000;
    scalar x3;
    int count = 0;
    int iter = 1000;

    for (count=0; count<iter; count++)
    {
        x3 = (x1 + x2)/2;

        scalar f1 = sigmaSB_*Ato_*(pow(Tt_,4) - pow(x1,4)) / ( 1/et_ + ((1-eg_)/eg_)*(dto_/dgi_) ) + (qoFunction(t_)*LCRAvg_*ag_*Ago_) - (hoFunction(t_)*Ago_*(x1-ToFunction(t_)) + eg_*sigmaSB_*Ago_*(pow(x1,4) - pow(skyTemp_,4)));

        scalar f3 = sigmaSB_*Ato_*(pow(Tt_,4) - pow(x3,4)) / ( 1/et_ + ((1-eg_)/eg_)*(dto_/dgi_) ) + (qoFunction(t_)*LCRAvg_*ag_*Ago_) - (hoFunction(t_)*Ago_*(x3-ToFunction(t_)) + eg_*sigmaSB_*Ago_*(pow(x3,4) - pow(skyTemp_,4)));

        if(neg(f1*f3)){x2 = x3;}

        else if((mag(x1 - x2) < 1e-2)||(f3 == 1e-5)){break;}

        else{x1 = x3;}

        count++;
    }

    scalar Tg_ = x3;
Info << "Tg:" << Tg_ << endl;
    scalar Qtg_ = sigmaSB_*Ato_*(pow(Tt_,4) - pow(Tg_,4)) / ( 1/et_ + ((1-eg_)/eg_)*(dto_/dgi_) );
Info << "Qtg:" << Qtg_ << endl;
//// To calculate gradient
    const fvPatchScalarField& kappap = lookupPatchField<volScalarField, scalar>(kappaName_);

    forAll(degree_, i)
    {
        gradient()[i] = ( tg_*at_*qoFunction(t_)*LCRFunction(degree_[i]) - Qtg_/Ato_)/kappap[i];
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::parabolicCollector4SEFoam::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("Kappa") << kappaName_ << token::END_STATEMENT << nl;
    os.writeKeyword("range") << range_ << token::END_STATEMENT << nl;
    os.writeKeyword("LCR1") << LCR1_ << token::END_STATEMENT << nl;
    os.writeKeyword("LCR2") << LCR2_ << token::END_STATEMENT << nl;
    os.writeKeyword("LCR3") << LCR3_ << token::END_STATEMENT << nl;
    os.writeKeyword("glassRadProp") << glassRadProp_ << token::END_STATEMENT << nl;
    os.writeKeyword("tubeRadProp") << tubeRadProp_ << token::END_STATEMENT << nl;
    os.writeKeyword("dgi") << dgi_ << token::END_STATEMENT << nl;
    os.writeKeyword("dgo") << dgo_ << token::END_STATEMENT << nl;
    os.writeKeyword("dti") << dti_ << token::END_STATEMENT << nl;
    os.writeKeyword("dto") << dto_ << token::END_STATEMENT << nl;
    os.writeKeyword("len") << len_ << token::END_STATEMENT << nl;
    os.writeKeyword("longDir") << longDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("radDir") << radDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("qo") << qo_ << token::END_STATEMENT << nl;
    os.writeKeyword("ho") << ho_ << token::END_STATEMENT << nl;
    os.writeKeyword("To") << To_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        parabolicCollector4SEFoam
    );
}

// ************************************************************************* //

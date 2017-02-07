/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_emballement_base.cpp
// Directory : $EMBALLEMENT_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_emballement_base.h>
#include <Param.h>
#include <DoubleTab.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Matrice_Morse.h>
#include <Schema_Implicite_base.h>

Implemente_base( Source_emballement_base, "Source_emballement_base", Source_base ) ;

// XD source_emballement source_base source_emballement 1 not_set

Sortie& Source_emballement_base::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_emballement_base::readOn( Entree& is )
{
// Source_base::readOn( is );
  temps_=-1;
  Param param(que_suis_je());
  param.ajouter("Ea",&Ea_,Param::REQUIRED);  // XD_ADD_P double not_set
  param.ajouter("Asei",&Asei_,Param::REQUIRED);  // XD_ADD_P double not_set
  param.lire_avec_accolades(is);
  return is;
}
DoubleTab& Source_emballement_base::ajouter(DoubleTab& resu ) const
{
  assert(resu.dimension(0)==volumes_.size());
  assert(resu.dimension(0)==facteur_.size());
  int size=resu.dimension(0);
  // si on est en implicite on prend la valeur future
  const DoubleTab& inco=equation().inconnue().futur(is_scheme_implicite_);
  //const DoubleTab& inco=equation().inconnue().valeurs();
  for (int i=0; i<size; i++)
    resu(i)-=volumes_(i)*inco(i)*facteur_(i);
  return resu;
}

void Source_emballement_base::contribuer_a_avec(const DoubleTab& resu , Matrice_Morse& mat) const
{

  int size=resu.dimension(0);
  for (int i=0; i<size; i++)
    mat.coef(i,i)+=volumes_(i)*facteur_(i);
}
DoubleTab& Source_emballement_base::calculer(DoubleTab& resu) const
{
  resu=0;
  ajouter(resu);
  return resu;
}
void Source_emballement_base::mettre_a_jour(double temps)
{
  if (!est_egal(temps,temps_))
    {
      assert(volumes_.size()==facteur_.size());
      const DoubleTab& T =equation().probleme().get_champ("Temperature").valeurs();
      int size=T.dimension(0);
      double R=8.314;
      for (int i=0; i<size; i++)
        {

          facteur_(i)=Asei_*exp(-Ea_/R/T(i));
        }
      temps_=temps;
    }
}

void Source_emballement_base::completer()
{
  remplir_volumes();
  facteur_.resize(volumes_.size());

  is_scheme_implicite_=0;
  if (sub_type(Schema_Implicite_base,equation().probleme().schema_temps()))
    is_scheme_implicite_=1;
}
void Source_emballement_base::associer_zones(const Zone_dis& zdis,const Zone_Cl_dis& zcldis)
{
  //  Source_base::associer_zones(zdis,zcldis);
}
void Source_emballement_base::associer_pb(const Probleme_base& pb)
{

  // Source_base::associer_pb(pb ) ;
}

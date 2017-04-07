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
// File      : Puissance_emballement_base.cpp
// Directory : $EMBALLEMENT_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Puissance_emballement_base.h>
#include <Param.h>
#include <DoubleTab.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Matrice_Morse.h>
#include <Milieu_base.h>
#include <Champ_Generique_base.h>

Implemente_base( Puissance_emballement_base, "Puissance_emballement_base", Source_base ) ;
// XD puissance_emballement source_base puissance_emballement 1 -S * factor; S found by nom_champ; expression factor= $Hsei * Wc * dz  $
Sortie& Puissance_emballement_base::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Puissance_emballement_base::readOn( Entree& is )
{
// Source_base::readOn( is );
  temps_=-1;
  Param param(que_suis_je());
  param.ajouter("nom_champ",&nom_champ_,Param::REQUIRED); // XD_ADD_P chaine not_set
  param.ajouter("Hsei",&Hsei_,Param::REQUIRED); // XD_ADD_P double not_set
  param.ajouter("Wc",&Wc_,Param::REQUIRED);  // XD_ADD_P double not_set
  param.ajouter("dz",&dz_,Param::REQUIRED);  // XD_ADD_P double not_set
  param.lire_avec_accolades(is);
  return is;
}
DoubleTab& Puissance_emballement_base::ajouter(DoubleTab& resu ) const
{
  assert(resu.dimension(0)==volumes_.size());
  assert(resu.dimension(0)==facteur_.size());
  int size=resu.dimension(0);
  // const DoubleTab& inco=equation().inconnue().valeurs();
  for (int i=0; i<size; i++)
    resu(i)+=volumes_(i)*facteur_(i);
  return resu;
}

void Puissance_emballement_base::contribuer_a_avec(const DoubleTab& resu , Matrice_Morse& mat) const
{
}
DoubleTab& Puissance_emballement_base::calculer(DoubleTab& resu) const
{
  resu=0;
  ajouter(resu);
  return resu;
}
void Puissance_emballement_base::mettre_a_jour(double temps)
{
  if (!est_egal(temps,temps_))
    {
      const Champ_Generique_base& ch =equation().probleme().get_champ_post(nom_champ_);
      Champ espace_stockage;
      const Champ_base& champ_ecriture = ch.get_champ(espace_stockage);
      const DoubleTab& S = champ_ecriture.valeurs();


      int size=S.dimension(0);
      assert(volumes_.size()==facteur_.size());
      assert(size==facteur_.size());
      double f=Hsei_*Wc_*dz_/rhocp_;
      for (int i=0; i<size; i++)
        {

          facteur_(i)= -S(i)*f;
        }
      temps_=temps;
    }
}

void Puissance_emballement_base::completer()
{
  remplir_volumes();
  facteur_.resize(volumes_.size());

  double rho=equation().milieu().masse_volumique()(0,0);
  double cp=equation().milieu().capacite_calorifique()(0,0);
  rhocp_=rho*cp;
}

void Puissance_emballement_base::associer_zones(const Zone_dis& zdis,const Zone_Cl_dis& zcldis)
{
  //  Source_base::associer_zones(zdis,zcldis);
}
void Puissance_emballement_base::associer_pb(const Probleme_base& pb)
{

  // Source_base::associer_pb(pb ) ;
}

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
#include <Discretisation_base.h>
#include <Champ_Fonc_Tabule.h>
#include <Champ_Don.h>
#include <EFichier.h>

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
  Nom lnAcell,Eacell;
  avec_relaxation_=0;
  Param param(que_suis_je());
  Param& param_lineaire =param.ajouter_param("lineaire") ;// XD_ADD_P emballement_lineaire  modele lineaire
  Param& param_non_lineaire = param.ajouter_param("non_lineaire") ;// XD_ADD_P emballement_non_lineaire  modele lineaire
// 2XD emballement_lineaire objet_lecture nul 1 emballement lineaire
  param_lineaire.ajouter("Ea",&Ea_,Param::REQUIRED);  // 2XD_ADD_P double not_set
  param_lineaire.ajouter("Asei",&Asei_,Param::REQUIRED);  // 2XD_ADD_P double not_set
// 2XD emballement_non_lineaire objet_lecture nul 1 emballement non_lineaire
  param_non_lineaire.ajouter("kb",&kb_,Param::REQUIRED);  // 2XD_ADD_P double not_set
  param_non_lineaire.ajouter("lnAcell_file",&lnAcell,Param::REQUIRED);  // 2XD_ADD_P chaine not_set
  param_non_lineaire.ajouter("Eacell_file",&Eacell,Param::REQUIRED);  // 2XD_ADD_P chaine not_set
  param_non_lineaire.ajouter_flag("test_relaxation",&avec_relaxation_); // 2XD_ADD_P flag pour tester la relaxation
  param.lire_avec_accolades(is);


  const LIST(Nom)& motlus=param.get_list_mots_lus();

  if (motlus.size()!=1)
    {
      Cerr<<"Error you should read a model and only one! "<<motlus<<finl;
      exit();
    }

  if (motlus[0]==Motcle("lineaire"))
    modele_numero_=0;
  else
    {
      modele_numero_=1;
      EFichier flnAcell(lnAcell);
      flnAcell >> table_lnAcell_;

      EFichier fEacell(Eacell);
      fEacell >> table_Eacell_;
      //table_lnAcell_=ref_cast(Champ_Fonc_Tabule,lnAcell.valeur()).table();
      //table_Eacell_=ref_cast(Champ_Fonc_Tabule,Eacell.valeur()).table();
    }
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
      if (modele_numero_==0)
        {
          double R=8.314;
          for (int i=0; i<size; i++)
            {

              facteur_(i)=Asei_*exp(-Ea_/(R*T(i)));
            }

        }
      else
        {
          const DoubleTab& c =equation().inconnue().valeurs();
          for (int i=0; i<size; i++)
            {
              double ci=c(i);
              double ai=1-ci;
              double lnA=table_lnAcell_.val(ai);
              double Ea=table_Eacell_.val(ai);
              facteur_(i)=exp(lnA-Ea/(kb_*T(i)))/(dabs(c(i))+DMINFLOAT);
              if (avec_relaxation_)
                {
                  facteur_(i)*=(1 -tanh(40.*(ai-0.85)))*0.5 ;
                }
            }

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

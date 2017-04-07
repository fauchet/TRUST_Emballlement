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
// File      : Puissance_convection_base.cpp
// Directory : $EMBALLEMENT_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Puissance_convection_base.h>
#include <Param.h>
#include <DoubleTab.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Matrice_Morse.h>
#include <Milieu_base.h>
#include <Champ_Generique_base.h>
#include <Schema_Implicite_base.h>

Implemente_base( Puissance_convection_base, "Puissance_convection_base", Source_base ) ;
// XD  puissance_convection source_base puissance_convection 1 source pour emballement $ 2* [ (hconv* (Tamb-T  ) ) +  (\epsilon*Stefan*  (Tamb^4 -T^4  ) )  ]$
Sortie& Puissance_convection_base::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Puissance_convection_base::readOn( Entree& is )
{
// Source_base::readOn( is );
  implicite_=1;
  temps_=-1;
  Param param(que_suis_je());
  param.ajouter("hconv",&hconv_,Param::REQUIRED);   // XD_ADD_P double not_set
  param.ajouter("epsilon",&epsilon_,Param::REQUIRED); // XD_ADD_P double not_set
  param.ajouter("Tamb",&Tamb_,Param::REQUIRED); // XD_ADD_P double not_set
  param.ajouter("Stefan",&Stefan_,Param::REQUIRED); // XD_ADD_P double not_set
//  param.ajouter("implicite",&implicite_); // NOXD_ADD_P int not_set
  param.lire_avec_accolades(is);
  return is;
}
DoubleTab& Puissance_convection_base::ajouter(DoubleTab& resu ) const
{
  assert(resu.dimension(0)==volumes_.size());
  int size=resu.dimension(0);
  //const DoubleTab& T=equation().inconnue().valeurs();
  // si on est en implicite on prend la valeur future
  const DoubleTab& T=equation().inconnue().futur(is_scheme_implicite_);
  double Tamb4=Tamb_*Tamb_;
  Tamb4*=Tamb4;
  double f=2./rhocp_;
  for (int i=0; i<size; i++)
    {
      double ti=T(i);
      double ti3;
      if (implicite_==0)
        ti3=ti*ti*ti;
      else
        ti3=T_old_(i)*T_old_(i)*T_old_(i);

      resu(i)+=volumes_(i)*f*(hconv_*(Tamb_-ti)+ epsilon_*Stefan_*(Tamb4-ti3*ti));
    }
  return resu;
}

void Puissance_convection_base::contribuer_a_avec(const DoubleTab& resu , Matrice_Morse& mat) const
{
  if (implicite_!=0)
    {
      int size=resu.dimension(0);
      double Tamb4=Tamb_*Tamb_;
      Tamb4*=Tamb4;
      double f=2./rhocp_;
      for (int i=0; i<size; i++)
        {
          double ti3;
          ti3=T_old_(i)*T_old_(i)*T_old_(i);

          mat.coef(i,i)+=volumes_(i)*f*(hconv_*(1)+ epsilon_*Stefan_*(ti3));
        }

    }
}
DoubleTab& Puissance_convection_base::calculer(DoubleTab& resu) const
{
  resu=0;
  ajouter(resu);
  return resu;
}
void Puissance_convection_base::mettre_a_jour(double temps)
{
  if (implicite_!=0)
    if (!est_egal(temps,temps_))
      {
        T_old_=equation().inconnue().valeurs();
        temps_=temps;
      }
}

void Puissance_convection_base::completer()
{
  remplir_volumes();

  double rho=equation().milieu().masse_volumique()(0,0);
  double cp=equation().milieu().capacite_calorifique()(0,0);
  rhocp_=rho*cp;

  is_scheme_implicite_=0;
  if (sub_type(Schema_Implicite_base,equation().probleme().schema_temps()))
    is_scheme_implicite_=1;
}

void Puissance_convection_base::associer_zones(const Zone_dis& zdis,const Zone_Cl_dis& zcldis)
{
  //  Source_base::associer_zones(zdis,zcldis);
}
void Puissance_convection_base::associer_pb(const Probleme_base& pb)
{

  // Source_base::associer_pb(pb ) ;
}

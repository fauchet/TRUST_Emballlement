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
// File      : Pb_Conduction_Concentration.cpp
// Directory : $EMBALLEMENT_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Pb_Conduction_Concentration.h>
#include <Constituant.h>


Implemente_instanciable( Pb_Conduction_Concentration, "Pb_Conduction_Concentration", Probleme_base ) ;
// XD pb_conduction_concentration Pb_base pb_conduction_concentration -1 not_set
// XD   attr conduction conduction conduction 0 Energy equation (temperature diffusion convection).
// XD   attr diffusion_concentration diffusion_concentration diffusion_concentration 0 Constituent transportation equations (concentration diffusion convection).
// XD   ref constituant constituant

Sortie& Pb_Conduction_Concentration::printOn( Sortie& os ) const
{
  Probleme_base::printOn( os );
  return os;
}

Entree& Pb_Conduction_Concentration::readOn( Entree& is )
{
  Probleme_base::readOn( is );
  return is;
}
//    Renvoie le nombre d'equations du probleme.
//    Toujours egal a 1 pour un probleme de conduction standart.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: int
//    Signification: le nombre d'equations du probleme
//    Contraintes: toujours egal a 1
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
int Pb_Conduction_Concentration::nombre_d_equations() const
{
  return 2;
}

// Description:
//    Renvoie l'equation de type Conduction
//     si i = 0,
//    Provoque une erreur sinon car le probleme
//    n'a qu'une seule equation.
//    (version const)
// Precondition:
// Parametre: int i
//    Signification: l'index de l'equation a renvoyer
//    Valeurs par defaut:
//    Contraintes: i doit etre egal a zero
//    Acces:
// Retour: Equation_base&
//    Signification: l'equation de type Conduction
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
const Equation_base& Pb_Conduction_Concentration::equation(int i) const
{
  if (i==0)
    return eq_conduction_;
  if (i==1)
    return eq_concentration_;
  Cerr << "Error in Pb_Conduction_Concentration::equation"<<finl;
  throw;
}


Equation_base& Pb_Conduction_Concentration::equation(int i)
{
  if (i==0)
    return eq_conduction_;
  if (i==1)
    return eq_concentration_;
  Cerr << "Error in Pb_Conduction_Concentration::equation"<<finl;
  throw;
}
void Pb_Conduction_Concentration::associer_milieu_base(const Milieu_base& mil)
{
  if ( sub_type(Constituant,mil) )
    eq_concentration_.associer_milieu_base(mil);
  else
    {
      eq_conduction_.associer_milieu_base(mil);
    }
}


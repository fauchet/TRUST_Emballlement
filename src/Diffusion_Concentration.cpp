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
// File      : Diffusion_Concentration.cpp
// Directory : $EMBALLEMENT_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Diffusion_Concentration.h>
#include <Param.h>

Implemente_instanciable( Diffusion_Concentration, "Diffusion_Concentration", Convection_Diffusion_Concentration ) ;
// XD diffusion_concentration eqn_base diffusion_concentration 1 d/dt C = div(D grad c) +S
// XD attr convection suppress_param convection 1 not_set
Sortie& Diffusion_Concentration::printOn( Sortie& os ) const
{
  Convection_Diffusion_Concentration::printOn( os );
  return os;
}

Entree& Diffusion_Concentration::readOn( Entree& is )
{
  assert(la_concentration.non_nul());

  Convection_Diffusion_std::readOn(is);
  Nom nom;
  nom="Diffusion_";
  nom+=inconnue().le_nom();
  terme_diffusif.set_fichier(nom);
  terme_diffusif.set_description((Nom)"Diffusion mass transfer rate=Integral(alpha*grad(C)*ndS) [m"+(Nom)(dimension+bidim_axi)+".Mol.s-1]");
  return is;

}
void Diffusion_Concentration::set_param(Param& param)
{
  Convection_Diffusion_Concentration::set_param(param);
  param.supprimer_condition("condition_2");

}
int Diffusion_Concentration::nombre_d_operateurs() const
{
  return 1;
}

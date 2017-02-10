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
// File      : Pb_Conduction_Concentration.h
// Directory : $EMBALLEMENT_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Pb_Conduction_Concentration_included
#define Pb_Conduction_Concentration_included

#include <Probleme_base.h>
#include <Conduction.h>
#include <Diffusion_Concentration.h>


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Pb_Conduction_Concentration
//
// <Description of class Pb_Conduction_Concentration>
//
/////////////////////////////////////////////////////////////////////////////

class Pb_Conduction_Concentration : public Probleme_base
{

  Declare_instanciable( Pb_Conduction_Concentration ) ;

public :
  int nombre_d_equations() const;
  const Equation_base& equation(int) const ;
  Equation_base& equation(int);
  void associer_milieu_base(const Milieu_base& mil);

protected :

  Conduction  eq_conduction_;
  Diffusion_Concentration eq_concentration_;
};

#endif /* Pb_Conduction_Concentration_included */

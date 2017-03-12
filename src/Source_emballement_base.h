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
// File      : Source_emballement_base.h
// Directory : $EMBALLEMENT_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_emballement_base_included
#define Source_emballement_base_included

#include <Source_base.h>
#include <Champ_Don.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_emballement_base
//
// <Description of class Source_emballement_base>
//
/////////////////////////////////////////////////////////////////////////////

class Source_emballement_base : public Source_base
{

  Declare_base( Source_emballement_base ) ;

public :
  virtual DoubleTab& ajouter(DoubleTab& ) const;
  virtual DoubleTab& calculer(DoubleTab& ) const;
  virtual void mettre_a_jour(double temps);
  virtual void completer();
  virtual void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const ;

protected :
  virtual void associer_zones(const Zone_dis& ,const Zone_Cl_dis& ) ;
  virtual void associer_pb(const Probleme_base& ) ;
  DoubleVect volumes_,facteur_;
  double Asei_,Ea_;

  virtual void remplir_volumes() =0;
  int is_scheme_implicite_;

  int modele_numero_;

  double temps_;
  double kb_;
  Champ_Don lnAcell_,Eacell_;
};

#endif /* Source_emballement_base_included */

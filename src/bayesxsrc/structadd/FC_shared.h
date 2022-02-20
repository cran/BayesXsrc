/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2019 Christiane Belitz, Andreas Brezger,
Nadja Klein, Thomas Kneib, Stefan Lang, Nikolaus Umlauf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */



#if !defined (FCsharedINCLUDED)

#define FCsharedINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"clstring.h"
#include"FC_nonp.h"
#include"design.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_shared -----------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FC_shared  : public FC
  {

  protected:

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

//  FC_shared(void);

  // CONSTRUCTOR

  FC_shared(void);

  // COPY CONSTRUCTOR

  FC_shared(const FC_shared & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_shared & operator=(const FC_shared & m);

  // DESTRUCTOR

  ~FC_shared()
    {
    }

  void update(void);

  bool posteriormode(void);

  void reset(void);
  };


} // end: namespace MCMC

#endif



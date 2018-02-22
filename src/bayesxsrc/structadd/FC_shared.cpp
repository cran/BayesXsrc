/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2011  Christiane Belitz, Andreas Brezger,
Thomas Kneib, Stefan Lang, Nikolaus Umlauf

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



#include "FC_shared.h"


//------------------------------------------------------------------------------
//----------------- CLASS: FC_shared implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{

FC_shared::FC_shared(void)
  {
  }

FC_shared::FC_shared(const FC_shared & m)
  : FC(FC(m))
  {
  }

const FC_shared & FC_shared::operator=(const FC_shared & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  return *this;
  }


void FC_shared::update(void)
  {
  }

bool FC_shared::posteriormode(void)
  {
  return true;
  }

void FC_shared::reset(void)
  {
  }

} // end: namespace MCMC




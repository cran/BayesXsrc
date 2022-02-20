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

#include "GENERAL_OPTIONS.h"

using std::flush;

namespace MCMC
{

//------------------------------------------------------------------------------
//-------------------------- CLASS: GENERAL_OPTIONS ----------------------------
//------------------------------------------------------------------------------


GENERAL_OPTIONS::GENERAL_OPTIONS(void)
  {
  iterations = 22000;
  burnin = 2000;
  step = 20;
  nrbetween = 1000;
  nrout = 100;
  nriter = 0;
  samplesize = 0;
  logout = &cout;
  set_level1(95);
  set_level2(80);
  saveestimation = false;
  copula = false;
  rotation = 0;
  samplesel = false;
  sampleselval = 0.0;
  IWLSlineff=false;
  forceIWLS=false;
  highspeedon=false;
  }


GENERAL_OPTIONS::GENERAL_OPTIONS(
const unsigned & it,const unsigned & bu,
                         const unsigned & st, const bool & sa,
                         const bool & cop, const unsigned & rot,
                         const bool & samsel, const double & samselval, const bool & fiwls,
                         const bool & hso,
                         ostream * lo,
                         const double & l1,const double & l2)
  {
  iterations = it;
  burnin = bu;
  step = st;
  set_level1(l1);
  set_level2(l2);

  nrbetween = (iterations-burnin)/3;
  if (nrbetween < 1000)
    nrbetween = 1000;
  nrout = 1000;
  nriter = 0;
  samplesize = 0;
  logout = lo;
  saveestimation = sa;
  copula = cop;
  rotation = rot;
  samplesel = samsel;
  sampleselval = samselval;
  forceIWLS = fiwls;
  highspeedon = hso;

  (*logout) << flush;
  }


GENERAL_OPTIONS::GENERAL_OPTIONS(const GENERAL_OPTIONS & o)
  {
  iterations = o.iterations;
  burnin = o.burnin;
  step = o.step;
  level1 = o.level1;
  level2 = o.level2;
  nrbetween = o.nrbetween;
  nrout = o.nrout;
  nriter = o.nriter;
  samplesize = o.samplesize;
  logout = o.logout;
  lower1 = o.lower1;
  lower2 = o.lower2;
  upper1 = o.upper1;
  upper2 = o.upper2;
  saveestimation = o.saveestimation;
  copula = o.copula;
  rotation = o.rotation;
  samplesel = o.samplesel;
  sampleselval = o.sampleselval;
  IWLSlineff = o.IWLSlineff;
  forceIWLS = o.forceIWLS;
  highspeedon = o.highspeedon;
  }


const GENERAL_OPTIONS & GENERAL_OPTIONS::operator=(const GENERAL_OPTIONS & o)
  {
  if (this == &o)
    return *this;
  iterations = o.iterations;
  burnin = o.burnin;
  step = o.step;
  level1 = o.level1;
  level2 = o.level2;
  nrbetween = o.nrbetween;
  nrout = o.nrout;
  nriter = o.nriter;
  samplesize = o.samplesize;
  logout = o.logout;
  lower1 = o.lower1;
  lower2 = o.lower2;
  upper1 = o.upper1;
  upper2 = o.upper2;
  saveestimation = o.saveestimation;
  copula = o.copula;
  rotation = o.rotation;
  samplesel = o.samplesel;
  sampleselval = o.sampleselval;
  IWLSlineff = o.IWLSlineff;
  forceIWLS = o.forceIWLS;
  highspeedon = o.highspeedon;
  return *this;
  }


void GENERAL_OPTIONS::out(const ST::string & s,bool thick,bool italic,
                          unsigned size,int r,int g, int b)
  {
  cout << s;
  if (!(logout->fail()))
    (*logout) << s << flush;
  }


void GENERAL_OPTIONS::outoptions(void)
  {

  out("GENERAL OPTIONS:\n",true);
  out("\n");
  out("  Number of iterations:  " + ST::inttostring(iterations) + "\n");
  out("  Burn-in period:        " + ST::inttostring(burnin)+ "\n");
  out("  Thinning parameter:    " + ST::inttostring(step)+ "\n");
  if (saveestimation)
    out("  Saveestimation:        enabled\n");
  else
    out("  Saveestimation:        disabled\n");
  out("\n");
  if (copula)
    {
    out(" Copula Model specified\n");
    out("\n");
    if((rotation==0)||(rotation==90)||(rotation==180)||(rotation==270))
      out("  Copula is rotated by "+ ST::inttostring(rotation) + "\n");
    else
      out("  Invalid angle of rotation specified. Copula will not be rotated\n");
    }
  else {}

  out("\n");
  }


void GENERAL_OPTIONS::update(void)
  {

  nriter++;

  if (nriter % nrout == 0 || nriter == 1)
    {
    out("  ITERATION: " + ST::inttostring(nriter) + "\n");
    }

  if( (nriter > burnin) && ((nriter-burnin-1) % step == 0) )
      samplesize++;
  }


unsigned GENERAL_OPTIONS::compute_samplesize(void)
  {
  return 1+(iterations-burnin-1)/step;
  }


void GENERAL_OPTIONS::reset(void)
  {
  nriter = 0;
  samplesize = 0;
  }


void GENERAL_OPTIONS::set_level1(double l1)
  {
  level1 = l1;

  lower1 = (100.0-level1)/2;
  upper2 = 100.0 - lower1;
  }

void GENERAL_OPTIONS::set_level2(double l2)
  {
  level2 = l2;
  lower2 = (100.0-level2)/2;
  upper1 = 100.0 - lower2;
  }


} // end: namespace MCMC




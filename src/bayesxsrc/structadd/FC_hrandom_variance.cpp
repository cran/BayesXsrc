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



#include "FC_hrandom_variance.h"


//------------------------------------------------------------------------------
//--------- CLASS: FC_hrandom_variance implementation of member functions ------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_hrandom_variance::read_options(vector<ST::string> & op,
vector<ST::string> & vn)
  {

  int f;

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  */

  FC_nonp_variance::read_options(op,vn);

  if (op[12] == "true")
    {
    mult =true;
    f = op[9].strtodouble(lambdastart);
    f = op[10].strtodouble(a_invgamma);
    f = op[11].strtodouble(b_invgamma);
    }
  else
    {
    mult = false;
    }

  }


FC_hrandom_variance::FC_hrandom_variance(void)
  {

  }


FC_hrandom_variance::FC_hrandom_variance(MASTER_OBJ * mp,
                 GENERAL_OPTIONS * o,DISTR * lp,
                  DISTR * lpRE,
                 const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC_nonp_variance(mp,o,lp,t,fp,Dp,FCn,op,vn)
  {
  read_options(op,vn);
  likepRE = lpRE;
  }


FC_hrandom_variance::FC_hrandom_variance(const FC_hrandom_variance & m)
  : FC_nonp_variance(FC_nonp_variance(m))
  {
  likepRE = m.likepRE;
  mult = m.mult;
  }


const FC_hrandom_variance & FC_hrandom_variance::operator=(
const FC_hrandom_variance & m)
  {

  if (this==&m)
	 return *this;
  FC_nonp_variance::operator=(FC_nonp_variance(m));
  likepRE = m.likepRE;
  mult = m.mult;
  return *this;
  }


double FC_hrandom_variance::compute_quadform(void)
  {

  unsigned n;
  double sum = 0;
  double * workbeta = FCnonpp->beta.getV();
  register unsigned i;

  n = FCnonpp->beta.rows();


  double * linpredREp;
  if (likepRE->linpred_current==1)
    linpredREp = likepRE->linearpred1.getV();
  else
    linpredREp = likepRE->linearpred2.getV();

  for(i=0;i<n;i++,workbeta++,linpredREp++)
    {
    sum += pow(*workbeta-(*linpredREp),2);
    }

  return sum;

  }

  
/*
void FC_hrandom_variance::transform_beta(void)
  {
  if (mult)
    transform(0,0) = 1;
  else
    FC_nonp_variance::transform_beta();
  }
*/

void FC_hrandom_variance::update(void)
  {

  b_invgamma = masterp->level1_likep->trmult*b_invgamma_orig;

  beta(0,0) = rand_invgamma(a_invgamma+0.5*designp->rankK,
                                  b_invgamma+0.5*compute_quadform());

  beta(0,1) = likep->get_scale()/beta(0,0);

  FCnonpp->tau2 = beta(0,0);
  likepRE->sigma2=beta(0,0);

//  transform_beta();
  acceptance++;
  FC::update();

  }


bool FC_hrandom_variance::posteriormode(void)
  {
  return  FC_nonp_variance::posteriormode();
  }


} // end: namespace MCMC




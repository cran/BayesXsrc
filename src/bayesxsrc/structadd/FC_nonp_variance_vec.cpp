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



#include "FC_nonp_variance_vec.h"


//------------------------------------------------------------------------------
//------ CLASS: FC_non_variance_vec implementation of member functions ---------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_nonp_variance_vec::read_options(vector<ST::string> & op,
                                   vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  */

  FC_nonp_variance::read_options(op,vn);

  }


FC_nonp_variance_vec::FC_nonp_variance_vec(void)
  {

  }



FC_nonp_variance_vec::FC_nonp_variance_vec(MASTER_OBJ * mp, GENERAL_OPTIONS * o,
                 DISTR * lp, const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC_nonp_variance(mp,o,lp,t,fp,Dp,FCn,op,vn)
  {


  read_options(op,vn);

  datamatrix betanew(FCn->beta.rows(),1,likep->get_scale()/lambdastart);
  setbeta(betanew);

  FCnonpp->tau2 = 1;
  FCnonpp->lambda = likep->get_scale();
  Dp->compute_penalty2(beta);


  }


FC_nonp_variance_vec::FC_nonp_variance_vec(const FC_nonp_variance_vec & m)
    : FC_nonp_variance(FC_nonp_variance(m))
  {
  }


const FC_nonp_variance_vec & FC_nonp_variance_vec::operator=(
        const FC_nonp_variance_vec & m)
  {

  if (this==&m)
	 return *this;
  FC_nonp_variance::operator=(FC_nonp_variance(m));
  return *this;
  }


void FC_nonp_variance_vec::update(void)
  {

  }


bool FC_nonp_variance_vec::posteriormode(void)
  {


  b_invgamma = masterp->level1_likep->trmult*b_invgamma_orig;

  unsigned i;
  for(i=0;i<beta.rows();i++)
    beta(i,0) = likep->get_scale()/lambdastart;

  FCnonpp->lambda = likep->get_scale();
  FCnonpp->tau2 = 1;
  designp->compute_penalty2(beta);

  posteriormode_betamean();


  return true;
  }



void FC_nonp_variance_vec::outresults(ofstream & out_stata,ofstream & out_R,
                                  const ST::string & pathresults)
  {

  /*
  FC::outresults(out_stata,out_R,"");

//  optionsp->out("\n");

  ST::string l1 = ST::doubletostring(optionsp->lower1,4);
  ST::string l2 = ST::doubletostring(optionsp->lower2,4);
  ST::string u1 = ST::doubletostring(optionsp->upper1,4);
  ST::string u2 = ST::doubletostring(optionsp->upper2,4);

  ST::string nl1 = ST::doubletostring(optionsp->lower1,4);
  ST::string nl2 = ST::doubletostring(optionsp->lower2,4);
  ST::string nu1 = ST::doubletostring(optionsp->upper1,4);
  ST::string nu2 = ST::doubletostring(optionsp->upper2,4);
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  ST::string vstr;


  if (optionsp->samplesize > 1)
    {

    vstr = "    Mean:         ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betamean(0,0),6) + "\n");

    vstr = "    Std. dev.:    ";

    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(sqrt(betavar(0,0)),6) + "\n");

    vstr = "    " + l1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_lower(0,0),6) + "\n");

    vstr = "    " + l2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_lower(0,0),6) + "\n");

    vstr = "    50% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu50(0,0),6) + "\n");

    vstr = "    " + u1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_upper(0,0),6) + "\n");

    vstr = "    " + u2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_upper(0,0),6) + "\n");

    optionsp->out("\n");

    optionsp->out("    Smoothing parameter\n");

    optionsp->out("\n");

    vstr = "    Mean:         ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betamean(0,1),6) + "\n");

    vstr = "    Std. dev.:    ";

    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(sqrt(betavar(0,1)),6) + "\n");

    vstr = "    " + l1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_lower(0,1),6) + "\n");

    vstr = "    " + l2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_lower(0,1),6) + "\n");

    vstr = "    50% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu50(0,1),6) + "\n");

    vstr = "    " + u1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_upper(0,1),6) + "\n");

    vstr = "    " + u2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_upper(0,1),6) + "\n");

    optionsp->out("\n");
    }
  else
    {
    optionsp->out("    Smoothing parameter: " +
    ST::doubletostring(betamean(0,1),6) + "\n");

    optionsp->out("\n");
    }


  if (pathresults.isvalidfile() != 1)
    {

    optionsp->out("    Results for the variance component are also stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    ofstream ou(pathresults.strtochar());

    if (optionsp->samplesize > 1)
      {
      ou << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
      nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
      }
    else
      {
      ou << "pmean" << endl;
      }

    ou << betamean(0,0) << "  ";
    if (optionsp->samplesize > 1)
      {
      ou << (betavar(0,0)<0.0?0.0:sqrt(betavar(0,0))) << "  ";
      ou << betaqu_l1_lower(0,0) << "  ";
      ou << betaqu_l2_lower(0,0) << "  ";
      ou << betaqu50(0,0) << "  ";
      ou << betaqu_l2_upper(0,0) << "  ";
      ou << betaqu_l1_upper(0,0) << "  ";
      ou << betamin(0,0) << "  ";
      ou << betamax(0,0) << "  " << endl;
      }

    optionsp->out("\n");
    }

    */
  }


void FC_nonp_variance_vec::outoptions(void)
  {

  FC_nonp_variance::outoptions();

  }


void FC_nonp_variance_vec::reset(void)
  {

double t= 0;
  }


} // end: namespace MCMC




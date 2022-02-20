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

#include "FC_nonp_variance.h"

#if defined(BayesX_gsl_included)
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#endif

namespace MCMC
{

//------------------------------------------------------------------------------
//---------- CLASS: FC_nonp_variance implementation of member functions ---------
//------------------------------------------------------------------------------

void FC_nonp_variance::read_options(vector<ST::string> & op,
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

  int f;

  f = op[4].strtodouble(lambdastart);
  f = op[5].strtodouble(a_invgamma);
  f = op[6].strtodouble(b_invgamma_orig);

  if (op[38] == "true")
    {
    lambdaconst = true;
    nosamples = true;
    nosamplessave = true;
    }


  f = op[47].strtodouble(tildea);
  f = op[48].strtodouble(tildeb);
  if (op[49] == "true")
    {
    cauchy = true;
    }
  else
    cauchy = false;

  if (op[50] == "true")
    {
    wei = true;
    }
  else
    wei = false;
  f = op[51].strtodouble(scaletau2);

  if (op[58] == "iwls_tau")
    {
    proposal = 2;
    }
  else if(op[58] == "gamma")
    {
    proposal = 1;
    }
  else if(op[58] == "iwls_logtau2")
    {
    proposal = 3;
    }
  else
    {
    proposal = 0;
    }

  if (op[64]=="invgamma")
    hyperprior = invgamma;
  else if (op[64] == "scaledep")
    hyperprior = scaledep;
  else if (op[64] == "hcauchy")
    hyperprior = hcauchy;
  else if (op[64] == "hnormal")
    hyperprior = hnormal;
  else if (op[64] == "aunif")
    hyperprior = aunif;

  if(scaletau2 == 0)
    {
    if ((hyperprior == scaledep)  || (wei==true))
    scaletau2 = 0.00877812;
    else if (hyperprior == hcauchy)
    scaletau2 = 0.01034553;
    else if (hyperprior == hnormal)
    scaletau2 = 0.1457644;
    else if (hyperprior == aunif)
    scaletau2 = 0.2723532;
    }

  }





FC_nonp_variance::FC_nonp_variance(void)
  {

  }



FC_nonp_variance::FC_nonp_variance(MASTER_OBJ * mp, unsigned & enr, GENERAL_OPTIONS * o,
                 DISTR * lp, const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC(o,t,1,2,fp)
  {
  FCnonpp = FCn;
  likep = lp;
  designp = Dp;
  masterp = mp;
  equationnr = enr,
  lambdaconst = false;

  read_options(op,vn);

  datamatrix betanew(1,2);
  betanew(0,0) = likep->get_scale()/lambdastart;
  betanew(0,1) = lambdastart;
  setbeta(betanew);

  FCnonpp->tau2 = beta(0,0);
  FCnonpp->lambda = beta(0,1);

  }


FC_nonp_variance::FC_nonp_variance(const FC_nonp_variance & m)
    : FC(FC(m))
  {
  FCnonpp = m.FCnonpp;
  likep = m.likep;
  designp = m.designp;
  masterp = m.masterp;
  equationnr = m.equationnr;
  a_invgamma = m.a_invgamma;
  b_invgamma_orig = m.b_invgamma_orig;
  b_invgamma = m.b_invgamma;
  lambdastart = m.lambdastart;
  lambdaconst = m.lambdaconst;
  tildea = m.tildea;
  tildeb = m.tildeb;
  cauchy = m.cauchy;
  wei = m.wei;
  scaletau2 = m.scaletau2;
  proposal = m.proposal;
  hyperprior = m.hyperprior;
  }


const FC_nonp_variance & FC_nonp_variance::operator=(const FC_nonp_variance & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  FCnonpp = m.FCnonpp;
  likep = m.likep;
  designp = m.designp;
  masterp = m.masterp;
  equationnr = m.equationnr;
  a_invgamma = m.a_invgamma;
  b_invgamma_orig = m.b_invgamma_orig;
  b_invgamma = m.b_invgamma;
  lambdastart = m.lambdastart;
  lambdaconst = m.lambdaconst;
  tildea = m.tildea;
  tildeb = m.tildeb;
  cauchy = m.cauchy;
  wei = m.wei;
  scaletau2 = m.scaletau2;
  proposal = m.proposal;
  hyperprior = m.hyperprior;
  return *this;
  }


void FC_nonp_variance::update(void)
  {

  // TEST
  //  ofstream out("c:\\bayesx\\test\\results\\param.res");
  //  (FCnonpp->param).prettyPrint(out);
  // END: TEST

   /* double clk = (double)CLK_TCK;
    clock_t clkbegin2 = clock();

    double counti = 0;
    double nrs = 1;*/
//    double tmp = 0;
//    std::ofstream out2;
//    out2.open ("C:\\tmp\\GIG2.raw", std::ofstream::out | std::ofstream::app);
    //out2 << "GIG2"  << endl;
   // for (counti=0;counti<nrs;counti++)
  //  {
//    tmp = randnumbers::GIG2(0, 0.5, 1.5);
//    out2 << tmp << endl;
  //  }
  //  out2 << endl;
  //  clock_t clkend2 = clock();
  //  double sec2 = (clkend2-clkbegin2)/clk;

   // cout << "time GIG2: " << sec2 << endl;

   // clock_t clkbegin = clock();
   // counti = 0;
 //   std::ofstream out;
//    out.open ("C:\\tmp\\GIG.raw", std::ofstream::out | std::ofstream::app);
    //out << "GIG"  << endl;
    //for (counti=0;counti<nrs;counti++)
    //{
//    tmp = randnumbers::GIG(0, 0.5, 1.5);
 //   out << tmp << endl;
    //}
    //out << endl;

    /*clock_t clkend = clock();
    double sec = (clkend-clkbegin)/clk;

    cout << "time GIG: " << sec << endl;*/

/*    cout << "proposal: " << proposal << endl;
    cout << "tildea: " << tildea << endl;
    cout << "tildeb: " << tildeb << endl;
    cout << "scaletau2: " << scaletau2 << endl;
    cout << "cauchy: " << cauchy << endl;
    cout << "wei: " << wei << endl;*/

  b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

 /* cout << "hyperprior" << hyperprior << endl;
  cout << "wei" << wei << endl;
  cout << "cauchy" << cauchy << endl;*/

  if (lambdaconst == false)
    {

    if (cauchy == true)
      {
      double quadf = designp->penalty_compute_quadform(FCnonpp->param);
      double gamma = rand_invgamma(designp->rankK/2 +tildea ,0.5*quadf+tildeb);

      double u = log(uniform());

      double fcold = -(0.5*designp->rankK+0.5)*log(beta(0,0))-1/(2*beta(0,0))*quadf-log(1+beta(0,0)) ;
      double fcnew = -(0.5*designp->rankK+0.5)*log(gamma)-1/(2*gamma)*quadf-log(1+gamma);
      double proposalold = -(tildea+designp->rankK/2+1)*log(beta(0,0)) - (0.5*quadf+tildeb) / beta(0,0);
      double proposalnew = -(tildea+designp->rankK/2+1)*log(gamma) - (0.5*quadf+tildeb) / gamma;
      if (u <= (fcnew - fcold - proposalnew + proposalold))
        {

        beta(0,0) = gamma;
        acceptance++;
        }

      }
    else if ( (wei == true) || (hyperprior == scaledep))
      {
      double quadf = designp->penalty_compute_quadform(FCnonpp->param);
      double u = log(uniform());

      if(proposal == 0) //simply draw proposal as if tau^2 would be IG distributed with tildea and tildeb
        {
        double gamma = rand_invgamma(designp->rankK/2 +tildea ,0.5*quadf+tildeb);

        double fcold = -(0.5*designp->rankK+0.5)*log(beta(0,0))-pow(beta(0,0)/scaletau2, 0.5)-1/(2*beta(0,0))*quadf;
        double fcnew = -(0.5*designp->rankK+0.5)*log(gamma)-pow(gamma/scaletau2, 0.5)-1/(2*gamma)*quadf;
        double proposalold = -(designp->rankK/2 +tildea+1)*log(beta(0,0)) - (0.5*quadf+tildeb) / beta(0,0);
        double proposalnew = -(designp->rankK/2 +tildea+1)*log(gamma) - (0.5*quadf+tildeb) / gamma;

        if (u <= (fcnew - fcold - proposalnew + proposalold))
          {

          beta(0,0) = gamma;
          acceptance++;
          }
        }
      else if(proposal == 1) // gamma approximation of weibull distribution (such that mode is also at zero and that E_gamma = E_weibull
        {
      /*  double tmp1 = 3.0;
        double tmp2 = 5.0;
        double ew = randnumbers::gamma_exact(tmp1) *scaletau2;
        double varw = scaletau2 * scaletau2 * (randnumbers::gamma_exact(tmp2) - pow(randnumbers::gamma_exact(tmp1), 2));
        double bgamma = ew / varw;
        double agamma = bgamma * ew;*/
        double tmp1 = 3.0;
        double ew = randnumbers::gamma_exact(tmp1) *scaletau2;
        double bgamma = 1/ew; // equal expectation for agamma = 1
        double agamma = 1; //yields mode at zero
        double p = -designp->rankK/2 + agamma;
        double a = 2 * bgamma;
        double b = quadf;
        double gamma = randnumbers::GIG2(p, a, b);

        double fcold = -(0.5*designp->rankK+0.5)*log(beta(0,0))-pow(beta(0,0)/scaletau2, 0.5)-1/(2*beta(0,0))*quadf;
        double fcnew = -(0.5*designp->rankK+0.5)*log(gamma)-pow(gamma/scaletau2, 0.5)-1/(2*gamma)*quadf;
        double proposalold = (p-1)*log(beta(0,0)) - 0.5*(a*beta(0,0)+b/beta(0,0));
        double proposalnew = (p-1)*log(gamma) - 0.5*(a*gamma+b/gamma);

        if (u <= (fcnew - fcold - proposalnew + proposalold))
          {

          beta(0,0) = gamma;
          acceptance++;
          }
        }
      else if(proposal == 2)// iwls proposal for tau
        {
        double vartau = 1/ ((3 * quadf / beta(0,0) - designp->rankK) / beta(0,0));
        double mutau = sqrt(beta(0,0)) + vartau * (-designp->rankK / sqrt(beta(0,0)) + quadf / pow(beta(0, 0), 1.5) - 1/sqrt(scaletau2));

        double gamma = mutau + rand_normal() * sqrt(vartau);
        double fcold = -(0.5*designp->rankK)*log(beta(0,0))-pow(beta(0,0)/scaletau2, 0.5)-1/(2*beta(0,0))*quadf;
        double fcnew = -(designp->rankK)*log(gamma)-pow(gamma*gamma/scaletau2, 0.5)-1/(2*gamma*gamma)*quadf;

        double vartauold = 1/ ((3 * quadf / (gamma*gamma) - designp->rankK) / (gamma*gamma));
        double mutauold = gamma + vartau * (-designp->rankK / gamma + quadf / pow(gamma, 3) - 1/sqrt(scaletau2));
        double proposalold = -0.5*log(vartauold)-0.5*pow((sqrt(beta(0,0))-mutauold), 2)/vartauold;
        double proposalnew = -0.5*log(vartau)-0.5*pow((gamma-mutau), 2)/vartau;

        if (u <= (fcnew - fcold - proposalnew + proposalold))
          {
          beta(0,0) = gamma*gamma;
          acceptance++;
          }
        }
      else // iwls proposal for log(tau^2)
        {
        double vartau = 1/(0.5*quadf/beta(0,0) + 0.25*sqrt(beta(0,0))/sqrt(scaletau2));
        double mutau = log(beta(0,0)) + vartau * (1 - 0.5*(designp->rankK+1) + 0.5*quadf/beta(0,0) - 0.5*sqrt(beta(0,0))/sqrt(scaletau2));

        double gamma = mutau + rand_normal() * sqrt(vartau);
        double fcold = log(beta(0,0)) - 0.5*(designp->rankK+1)*log(beta(0,0)) - 0.5*pow(beta(0,0)/scaletau2, 0.5) - 1/(2*beta(0,0))*quadf;
        double fcnew = gamma- 0.5*(designp->rankK+1)*gamma - 0.5*pow(exp(gamma)/scaletau2, 0.5) - 1/(2*exp(gamma))*quadf;

        double vartauold = 1/(0.5*quadf/exp(gamma) + 0.25*sqrt(exp(gamma))/sqrt(scaletau2));
        double mutauold = gamma + vartauold * (1 - 0.5*(designp->rankK+1) + 0.5*quadf/exp(gamma) - 0.5*sqrt(exp(gamma))/sqrt(scaletau2));
        double proposalold = -0.5*log(vartauold)-0.5*pow((log(beta(0,0))-mutauold), 2)/vartauold;
        double proposalnew = -0.5*log(vartau)-0.5*pow((gamma-mutau), 2)/vartau;

        if (u <= (fcnew - fcold - proposalnew + proposalold))
          {
          beta(0,0) = exp(gamma);
          acceptance++;
          }
        }
      }
    else if (hyperprior == hnormal) // half normal for tau
      {
   //   cout << "hyperprior" << hyperprior << endl;
      double quadf = designp->penalty_compute_quadform(FCnonpp->param);
      double u = log(uniform());

      double vartau = 1/(0.5*quadf/beta(0,0) + 0.5*beta(0,0)/(scaletau2*scaletau2));
      double mutau = log(beta(0,0)) + vartau * (1 - 0.5*(designp->rankK+1) + 0.5*quadf/beta(0,0) - 0.5*beta(0,0)/(scaletau2*scaletau2));

      double gamma = mutau + rand_normal() * sqrt(vartau);
      double fcold = log(beta(0,0)) - 0.5*(designp->rankK+1)*log(beta(0,0)) - 0.5*beta(0,0)/(scaletau2*scaletau2) - 1/(2*beta(0,0))*quadf;
      double fcnew = gamma- 0.5*(designp->rankK+1)*gamma - 0.5*exp(gamma)/(scaletau2*scaletau2) - 1/(2*exp(gamma))*quadf;

      double vartauold = 1/(0.5*quadf/exp(gamma) + 0.5*exp(gamma)/(scaletau2*scaletau2));
      double mutauold = gamma + vartauold * (1 - 0.5*(designp->rankK+1) + 0.5*quadf/exp(gamma) - 0.5*exp(gamma)/(scaletau2*scaletau2));
      double proposalold = -0.5*log(vartauold)-0.5*pow((log(beta(0,0))-mutauold), 2)/vartauold;
      double proposalnew = -0.5*log(vartau)-0.5*pow((gamma-mutau), 2)/vartau;

      if (u <= (fcnew - fcold - proposalnew + proposalold))
        {
        beta(0,0) = exp(gamma);
        acceptance++;
        }
      }
    else if (hyperprior == hcauchy) // half cauchy for tau
      {
      double quadf = designp->penalty_compute_quadform(FCnonpp->param);
      double quadstau2=(scaletau2*scaletau2);
      double u = log(uniform());

      double vartau = 1/(0.5*quadf/beta(0,0) + (beta(0,0)/quadstau2)/pow((1+beta(0,0)/quadstau2),2));
      double mutau = log(beta(0,0)) + vartau * (1 - 0.5*(designp->rankK+1) + 0.5*quadf/beta(0,0) - (beta(0,0)/quadstau2)/(1+beta(0,0)/quadstau2));

      double gamma = mutau + rand_normal() * sqrt(vartau);
      double fcold = log(beta(0,0)) - 0.5*(designp->rankK+1)*log(beta(0,0)) - log(1+beta(0,0)/quadstau2) - 1/(2*beta(0,0))*quadf;
      double fcnew = gamma- 0.5*(designp->rankK+1)*gamma - log(1+exp(gamma)/quadstau2) - 1/(2*exp(gamma))*quadf;

      double vartauold = 1/(0.5*quadf/exp(gamma) + (exp(gamma)/quadstau2)/pow((1+exp(gamma)/quadstau2),2));
      double mutauold = gamma + vartauold * (1 - 0.5*(designp->rankK+1) + 0.5*quadf/exp(gamma) - (exp(gamma)/quadstau2)/(1+exp(gamma)/quadstau2));
      double proposalold = -0.5*log(vartauold)-0.5*pow((log(beta(0,0))-mutauold), 2)/vartauold;
      double proposalnew = -0.5*log(vartau)-0.5*pow((gamma-mutau), 2)/vartau;

      if (u <= (fcnew - fcold - proposalnew + proposalold))
        {
        beta(0,0) = exp(gamma);
        acceptance++;
        }
      }
    else if (hyperprior == aunif)  // approximation of uniform prior for tau
      {
      double tildec = 13.86294;
      double argold = sqrt(beta(0,0))*tildec/scaletau2;
      double expargold = exp(argold-tildec);
      double quadf = designp->penalty_compute_quadform(FCnonpp->param);
      double u = log(uniform());

      double vartau = 1/(0.5*quadf/beta(0,0) + ((argold*0.25*expargold)/(1+expargold))*(1+argold/(1+expargold)));
      double mutau = log(beta(0,0)) + vartau * (1 - 0.5*(designp->rankK+1) + 0.5*quadf/beta(0,0) - 0.5*argold*expargold/(1+expargold));

      double gamma = mutau + rand_normal() * sqrt(vartau);
      double argnew = sqrt(exp(gamma))*tildec/scaletau2;
      double expargnew = exp(argnew-tildec);
      double fcold = log(beta(0,0)) - 0.5*(designp->rankK+1)*log(beta(0,0)) - log(1+expargold) - 1/(2*beta(0,0))*quadf;
      double fcnew = gamma- 0.5*(designp->rankK+1)*gamma - log(1+expargnew) - 1/(2*exp(gamma))*quadf;

      double vartauold = 1/(0.5*quadf/exp(gamma) + ((argnew*0.25*expargnew)/(1+expargnew))*(1+argnew/(1+expargnew)));
      double mutauold = gamma + vartauold * (1 - 0.5*(designp->rankK+1) + 0.5*quadf/exp(gamma) - 0.5*argnew*expargnew/(1+expargnew));
      double proposalold = -0.5*log(vartauold)-0.5*pow((log(beta(0,0))-mutauold), 2)/vartauold;
      double proposalnew = -0.5*log(vartau)-0.5*pow((gamma-mutau), 2)/vartau;

      if (u <= (fcnew - fcold - proposalnew + proposalold))
        {
        beta(0,0) = exp(gamma);
        acceptance++;
        }
      }
    else
      {
      beta(0,0) = rand_invgamma(a_invgamma+0.5*designp->rankK,
                  b_invgamma+0.5*designp->penalty_compute_quadform(FCnonpp->param));

//      cout << a_invgamma << endl;
//      cout << b_invgamma << endl;
//      cout << designp->rankK << endl;
//     cout << designp->penalty_compute_quadform(FCnonpp->param) << endl;


      acceptance++;
      }

    beta(0,1) = likep->get_scale()/beta(0,0);

    FCnonpp->tau2 = beta(0,0);

    FC::update();

    }

  }


bool FC_nonp_variance::posteriormode(void)
  {
  b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

  beta(0,0) = likep->get_scale()/beta(0,1);

  FCnonpp->tau2 = beta(0,0);

  posteriormode_betamean();

  return true;
  }



void FC_nonp_variance::outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                                  const ST::string & pathresults)
  {

  FC::outresults(out_stata,out_R,out_R2BayesX,"");

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

    FC::outresults_acceptance();

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

    }
  else
    {
    optionsp->out("    Smoothing parameter: " +
    ST::doubletostring(betamean(0,1),6) + "\n");

    optionsp->out("\n");
    }

//  out_R << "term=" << title <<  ";" << endl;

  if (pathresults.isvalidfile() != 1)
    {

    optionsp->out("    Results for the variance component are also stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    ST::string paths = pathresults.substr(0,pathresults.length()-4) +
                                 "_sample.raw";

    out_R2BayesX << "pathvarsample=" << paths << endl;
//    out_R << "filetype=param; path=" << pathresults << ";" <<  endl;

    ofstream ou(pathresults.strtochar());

    if (optionsp->samplesize > 1)
      {
      ou << "pmean  pstd  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
      nu1 << "   pqu" << nu2 << endl;
      }
    else
      {
      ou << "pmean" << endl;
      }

    ou << betamean(0,0) << "  ";
    if (optionsp->samplesize > 1)
      {
      if (betavar(0,0) < 0.0000000000001)
        ou << 0 << "  ";
      else
        ou << sqrt(betavar(0,0)) << "  ";
      ou << betaqu_l1_lower(0,0) << "  ";
      ou << betaqu_l2_lower(0,0) << "  ";
      ou << betaqu50(0,0) << "  ";
      ou << betaqu_l2_upper(0,0) << "  ";
      ou << betaqu_l1_upper(0,0) << "  " << endl;
      }

    optionsp->out("\n");
    }

  }


void FC_nonp_variance::outoptions(void)
  {
  if (cauchy)
    {
    optionsp->out("  Cauchy prior\n");

    optionsp->out("  Hyperparameter tildea for proposal density: " +
                ST::doubletostring(tildea) + "\n" );

    optionsp->out("  Hyperparameter tildeb for proposal density: " +
                ST::doubletostring(tildeb) + "\n" );

    }
  else if (wei  || (hyperprior==scaledep))
    {
    optionsp->out("  Scale-dependent (weibull) prior\n");

    optionsp->out("  Scale parameter: " +
                ST::doubletostring(scaletau2) + "\n" );

    if(proposal == 0)
      {
      optionsp->out("  Inverse gamma proposal density \n" );
      optionsp->out("  Hyperparameter tildea for proposal density: " +
                ST::doubletostring(tildea) + "\n" );

      optionsp->out("  Hyperparameter tildeb for proposal density: " +
                ST::doubletostring(tildeb) + "\n" );
      }
    else if(proposal == 1)
      {
      optionsp->out("  Generalised inverse Gaussian proposal density \n" );
      }
    else if(proposal == 2)
      {
      optionsp->out("  IWLS proposal density for tau \n" );
      }
    else
      {
      optionsp->out("  IWLS proposal density for log(tau^2) \n" );
      }
    }
  else if (hyperprior==hnormal)
    {
    optionsp->out("  Gamma prior (Half normal prior for tau)\n");

    optionsp->out("  Scale parameter: " +
                ST::doubletostring(scaletau2) + "\n" );

    optionsp->out("  IWLS proposal density for log(tau^2) \n" );
    }
  else if (hyperprior==hcauchy)
    {
    optionsp->out("  Generalised beta prime prior (Half cauchy prior for tau)\n");

    optionsp->out("  Scale parameter: " +
                ST::doubletostring(scaletau2) + "\n" );

    optionsp->out("  IWLS proposal density for log(tau^2) \n" );
    }
  else if (hyperprior==aunif)
    {
    optionsp->out("  Approximation of uniform prior for tau\n");

    optionsp->out("  Scale parameter: " +
                ST::doubletostring(scaletau2) + "\n" );

    optionsp->out("  IWLS proposal density for log(tau^2) \n" );
    }
  else
    {

    optionsp->out("  Inverse gamma prior\n");

    b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

    optionsp->out("  Hyperprior a for variance parameter: " +
                ST::doubletostring(a_invgamma) + "\n" );
    optionsp->out("  Hyperprior b for variance parameter: " +
                ST::doubletostring(b_invgamma) + "\n" );
    optionsp->out("\n");
    }
  }


void FC_nonp_variance::reset(void)
  {

  datamatrix betanew(1,2);
  betanew(0,0) = likep->get_scale()/lambdastart;
  betanew(0,1) = lambdastart;
  setbeta(betanew);

  FCnonpp->tau2 = beta(0,0);
  FCnonpp->lambda = beta(0,1);
//  transform(0,0) = 1;
//  transform(1,0) = 1;

  }



//------------------------------------------------------------------------------
//--- CLASS: FC_nonp_variance_varselection implementation of member functions --
//------------------------------------------------------------------------------

void FC_nonp_variance_varselection::read_options(vector<ST::string> & op,
                                   vector<ST::string> & vn)
  {
  FC_nonp_variance::read_options(op,vn);

  datanames=vn;

  int f;

  f = op[39].strtodouble(a_omega);
  f = op[40].strtodouble(b_omega);

  f = op[41].strtodouble(r);

  f = op[53].strtodouble(v1);
  f = op[54].strtodouble(v2_orig);
  v2 = v2_orig;
//  v2 = masterp->level1_likep[equationnr]->trmult*v2_orig;

  if (op[57] == "true")
    gig = true;
  else
    gig = false;

  if (op[80]=="regcoeff")
    ssvsupdate = regcoeff;
  else if (op[80] == "sdev")
    ssvsupdate = sdev;
  else if (op[80] == "gibbs")
    ssvsupdate = gibbs;

  }


FC_nonp_variance_varselection::FC_nonp_variance_varselection(void)
  {

  }



FC_nonp_variance_varselection::FC_nonp_variance_varselection(MASTER_OBJ * mp,
                 unsigned & enr, GENERAL_OPTIONS * o,
                 DISTR * lp, bool so, const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn, double vl)
     : FC_nonp_variance(mp,enr,o,lp,t,fp,Dp,FCn,op,vn)
  {

  read_options(op,vn);

  ssvsvarlimit = vl;

  singleomega = so;

  FC_delta = FC(o,"",1,2,"");
  FC_delta.setbeta(1,1,0);
  FC_delta.setbeta(1,2,0.5);

  FC_psi2 = FC(o,"",1,1,"");
  FC_psi2.setbeta(1,1,0.5);

  omega=0.5;
  if(!singleomega)
    {
    FC_omega = FC(o,"",1,1,"");
    FC_omega.setbeta(1,1,omega);
    }

  tauold = sqrt(beta(0,0));

  if(FCnonpp->IWLS)
    {
    diff = datamatrix(likep->nrobs, 1, 0);
    Xscaled = datamatrix(likep->nrobs, 1, 0);

//    zerovec = datamatrix(likep->nrobs,1,0);
//    tauold = sqrt(likep->get_scale()/FCnonpp->lambda);
//    FCnonpp->tau2=1;
//    FCnonpp->lambda=likep->get_scale();
//    FCnonpp->designp->set_intvar(zerovec,tauold);
//    FCnonpp->designp->changingdesign=true;
    }
  }


FC_nonp_variance_varselection::FC_nonp_variance_varselection(const FC_nonp_variance_varselection & m)
    : FC_nonp_variance(FC_nonp_variance(m))
  {
  zerovec = m.zerovec;
  singleomega = m.singleomega;
  FC_delta = m.FC_delta;
  FC_psi2 = m.FC_psi2;
  FC_omega = m.FC_omega;
  a_omega = m.a_omega;
  b_omega = m.b_omega;
  omega = m.omega;
  tauold = m.tauold;
  v1 = m.v1;
  v2 = m.v2;
  v2_orig = m.v2_orig;
  r = m.r;
  X = m.X;
  Xscaled = m.Xscaled;
  diff = m.diff;
  gig = m.gig;
  r_delta = m.r_delta;
  datanames = m.datanames;
  ssvsvarlimit = m.ssvsvarlimit;
  ssvsupdate = m.ssvsupdate;
  }


const FC_nonp_variance_varselection & FC_nonp_variance_varselection::operator=(const FC_nonp_variance_varselection & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  zerovec = m.zerovec;
  singleomega = m.singleomega;
  FC_delta = m.FC_delta;
  FC_psi2 = m.FC_psi2;
  FC_omega = m.FC_omega;
  a_omega = m.a_omega;
  b_omega = m.b_omega;
  omega = m.omega;
  tauold = m.tauold;
  v1 = m.v1;
  v2 = m.v2;
  v2_orig = m.v2_orig;
  r = m.r;
  X = m.X;
  Xscaled = m.Xscaled;
  diff = m.diff;
  gig = m.gig;
  r_delta = m.r_delta;
  datanames = m.datanames;
  ssvsvarlimit = m.ssvsvarlimit;
  ssvsupdate = m.ssvsupdate;
  return *this;
  }


void FC_nonp_variance_varselection::add_linpred(datamatrix & l)
  {
  if (likep->linpred_current==1)
    likep->linearpred1.plus(l);
  else
    likep->linearpred2.plus(l);
  }


void FC_nonp_variance_varselection::update_IWLS(void)
  {
  unsigned i;

  double Sigmatauprop, Sigmataucurr;
  double mutauprop=0.0;
  double mutaucurr=0.0;

  // compute X = design matrix * regression coefficients (Z * gamma)
  FCnonpp->designp->compute_effect(X,FCnonpp->beta);

//  ofstream out1("c:/temp/X.raw");
//  X.prettyPrint(out1);
//  out1.close();
//  cout << FCnonpp->param(0,0) << endl;

  tauold = sqrt(beta(0,0));
  double tauoldinv = 1 / tauold;
  Xscaled.mult_scalar(X, tauoldinv);

  double * worklin;
  if (likep->linpred_current==1)
    worklin = likep->linearpred1.getV();
  else
    worklin = likep->linearpred2.getV();

  double logold = likep->compute_iwls(true, true);

  double * Xp = X.getV();
  double * Xsp = Xscaled.getV();
  double * responsep = likep->workingresponse.getV();
  double * weightp = likep->workingweight.getV();
  double xtx=0.0;
  for (i=0;i<X.rows();i++,Xp++,Xsp++,responsep++,worklin++,weightp++)
    {
    xtx += (*weightp) * pow(*Xsp,2);
    mutaucurr += (*weightp) * (*Xsp) * ((*responsep) - (*worklin) + (*Xp));
    }
  Sigmataucurr = 1/(xtx + 1/(FC_psi2.beta(0,0)*r_delta));

  mutaucurr *= Sigmataucurr;

  double tauprop = mutaucurr + sqrt(Sigmataucurr) * rand_normal();
  double tauratio = tauprop/tauold;

  double * diffp = diff.getV();
  Xp = X.getV();
  for(i=0; i<X.rows(); i++, Xp++, diffp++)
    {
    *diffp = (tauratio - 1) * (*Xp);
    }
  add_linpred(diff);

//  cout << "tauprop: " << tauprop << endl;
//  cout << "tauold: " << tauold << endl;
//  cout << "tauratio: " << tauratio << endl;
//  cout << "beta: " << FCnonpp->param(0,0) << endl;
//  cout << "acceptance: " << acceptance << endl;

  double logprop = likep->compute_iwls(true, true);

  Xp = X.getV();
  Xsp = Xscaled.getV();
  responsep = likep->workingresponse.getV();
  weightp = likep->workingweight.getV();
  if (likep->linpred_current==1)
    worklin = likep->linearpred1.getV();
  else
    worklin = likep->linearpred2.getV();

  xtx=0;
  for (i=0;i<X.rows();i++,Xp++,Xsp++,responsep++,worklin++,weightp++)
    {
    xtx += (*weightp) * pow(*Xsp,2);
    mutauprop += (*weightp) * (*Xsp) * ((*responsep) - (*worklin) + tauratio*(*Xp));
    }
  Sigmatauprop = 1/(xtx + 1/(FC_psi2.beta(0,0)*r_delta));
  mutauprop *= Sigmatauprop;

  logprop += -0.5*tauprop*tauprop / (FC_psi2.beta(0,0)*r_delta);
  logold += -0.5*tauold*tauold / (FC_psi2.beta(0,0)*r_delta);

  double u = log(uniform());

  double qnew = -0.5*log(Sigmataucurr) -0.5*(tauprop-mutaucurr)*(tauprop-mutaucurr)/Sigmataucurr;
  double qold = -0.5*log(Sigmatauprop) -0.5*(tauold-mutauprop)*(tauold-mutauprop)/Sigmatauprop;

  if(u <= logprop - logold - qnew + qold)
    {
    acceptance++;
    if(fabs(tauprop)<=ssvsvarlimit)
      {
      tauprop = tauold;
      tauratio = 1.0;
      }

    FCnonpp->ssvs_update(tauratio,false,false);

    beta(0,0) = tauprop*tauprop;
    beta(0,1) = 1/beta(0,0);

    FCnonpp->tau2 = beta(0,0);
    }
  else
    {
    diffp = diff.getV();
    for(i=0; i<diff.rows(); i++, diffp++)
      *diffp = -(*diffp);
    add_linpred(diff);
    tauratio=1.0;
    FCnonpp->ssvs_update(tauratio,false,true);
    }
  }



void FC_nonp_variance_varselection::update_gaussian(void)
  {
  unsigned i;
  double Sigmatau=0.0;
  double mutau=0.0;

  FCnonpp->designp->compute_effect(X,FCnonpp->beta);

  double * worklin;
  if (likep->linpred_current==1)
    worklin = likep->linearpred1.getV();
  else
    worklin = likep->linearpred2.getV();

  double * Xp = X.getV();
  double * responsep = likep->workingresponse.getV();
  double * wp = likep->workingweight.getV();
  double varinv = 1/(likep->get_scale()*beta(0,0));
  double xtx=0;
  for (i=0;i<X.rows();i++,Xp++,responsep++,worklin++,wp++)
    {
    xtx += pow(*Xp,2) * *wp;
    mutau += (*Xp) * *wp * ((*responsep) - (*worklin)+(*Xp));
    }
  Sigmatau = 1/(varinv*xtx + 1/(r_delta*FC_psi2.beta(0,0)));
  mutau *= Sigmatau/(likep->get_scale()*sqrt(beta(0,0)));

  double tau = mutau + sqrt(Sigmatau) * rand_normal();

  double tau2 = tau*tau;
  if (tau2 < 0.000000001)
    tau2 = 0.000000001;

  beta(0,0) = tau2;

  beta(0,1) = likep->get_scale()/beta(0,0);

  FCnonpp->tau2 = beta(0,0);

  acceptance++;

  }



void FC_nonp_variance_varselection::update(void)
  {
  // updating psi2
  if (FC_delta.beta(0,0) < 1)
    r_delta = r;
  else
    r_delta = 1;

  v2 = masterp->level1_likep[equationnr]->trmult*v2_orig;
//  v2 = v2_orig;

  FC_psi2.beta(0,0) = rand_invgamma(v1+0.5,v2+0.5*beta(0,0)/r_delta);
  FC_psi2.update();
  // end: updating psi2

  // updating delta
  double u = uniform();
  double L = 1/sqrt(r)*exp(- beta(0,0)/(2*FC_psi2.beta(0,0))*(1/r-1));
  double pr1 = 1/(1+ ((1-omega)/omega)*L);

  if (u <=pr1)
    {
    FC_delta.beta(0,0) = 1;
    r_delta = 1;
    }
  else
    {
    FC_delta.beta(0,0) = 0;
    r_delta = r;
    }

  FC_delta.beta(0,1) = pr1;
  FC_delta.update();
  // end: updating delta

  // updating w
  if(!singleomega)
    {
    FC_omega.beta(0,0) = randnumbers::rand_beta(a_omega+FC_delta.beta(0,0),
                                          b_omega+1-FC_delta.beta(0,0));
    omega = FC_omega.beta(0,0);

    FC_omega.update();
    }
  // end: updating w

  if(ssvsupdate==gibbs)
    {
    double p = -designp->rankK/2 + 0.5;
    double a = 1/(FC_psi2.beta(0,0)*r_delta);
    double b = designp->penalty_compute_quadform(FCnonpp->param);
    double tau2 = randnumbers::GIG2(p, a, b);
    beta(0,0) = tau2;
    beta(0,1) = 1/beta(0,0);
    tauold=sqrt(tau2);

    FCnonpp->tau2 = beta(0,0);
    double tauratio=1.0;
    FCnonpp->ssvs_update(tauratio,false,true);

    acceptance++;
    }
  else if(ssvsupdate==sdev)
    {
    if(optionsp->nriter > 100)
      {
      double qf = designp->penalty_compute_quadform(FCnonpp->param);
      double logold = -0.5*designp->rankK*log(tauold*tauold)
                      -0.5/(tauold*tauold)*qf
                      -0.5*tauold*tauold/(FC_psi2.beta(0,0)*r_delta);
      double v = -designp->rankK / tauold
                 +qf / (tauold*tauold*tauold)
                 -tauold/(FC_psi2.beta(0,0)*r_delta);
      double w = -designp->rankK / (tauold*tauold)
                 +3*qf / (tauold*tauold*tauold*tauold)
                 +1/(FC_psi2.beta(0,0)*r_delta);
      double sigmataucurr = 1/w;
      double mutaucurr = tauold + sigmataucurr*v;
      double tauprop = mutaucurr + sqrt(sigmataucurr) * rand_normal();

      double logprop = -0.5*designp->rankK*log(tauprop*tauprop)
                      -0.5/(tauprop*tauprop)*qf
                      -0.5*tauprop*tauprop/(FC_psi2.beta(0,0)*r_delta);
      v = -designp->rankK / tauprop
                 +qf / (tauprop*tauprop*tauprop)
                   -tauprop/(FC_psi2.beta(0,0)*r_delta);
      w = -designp->rankK / (tauprop*tauprop)
                 +3*qf / (tauprop*tauprop*tauprop*tauprop)
                 +1/(FC_psi2.beta(0,0)*r_delta);
      double sigmatauprop = 1/w;
      double mutauprop = tauprop + sigmatauprop*v;

      double u = log(uniform());

      double qnew = -0.5*log(sigmataucurr) -0.5*(tauprop-mutaucurr)*(tauprop-mutaucurr)/sigmataucurr;
      double qold = -0.5*log(sigmatauprop) -0.5*(tauold-mutauprop)*(tauold-mutauprop)/sigmatauprop;

      if(u <= logprop - logold - qnew + qold)
        {
        acceptance++;

        if(fabs(tauprop)<=ssvsvarlimit)
          tauprop=tauold;

        beta(0,0) = tauprop*tauprop;
        beta(0,1) = 1/beta(0,0);
        tauold=tauprop;

        FCnonpp->tau2 = beta(0,0);
        }
      else
        {
        }
      double tauratio=1.0;
      FCnonpp->ssvs_update(tauratio,false,true);
      }
    else
      {
      if (FCnonpp->IWLS)
        {
        update_IWLS();
        }
      else
        {
        update_gaussian();
        }
      }
    }
  else if(ssvsupdate==regcoeff)
    {
    if (FCnonpp->IWLS)
      {
      update_IWLS();
      }
    else
      {
      update_gaussian();
      }
    }
  FC::update();
  }

bool FC_nonp_variance_varselection::posteriormode(void)
  {
  return true;
  }

void FC_nonp_variance_varselection::outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                                  const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    ST::string pathresults_delta = pathresults.substr(0,pathresults.length()-4) + "_delta.res";
    ST::string pathresults_delta_prob = pathresults.substr(0,pathresults.length()-4) + "_delta_prob.res";
    ST::string pathresults_omega = pathresults.substr(0,pathresults.length()-4) + "_omega.res";
    ST::string pathresults_psi2 = pathresults.substr(0,pathresults.length()-4) + "_psi2.res";

    if(singleomega == false)
      {
      FC_omega.outresults(out_stata,out_R,out_R2BayesX,pathresults_omega);
      }


    FC_nonp_variance::outresults(out_stata,out_R,out_R2BayesX,pathresults);

    optionsp->out("    Inclusion probability:\n");
    optionsp->out("\n");

    FC_delta.outresults(out_stata,out_R,out_R2BayesX,pathresults_delta);
    FC_delta.outresults_help(out_stata,out_R,pathresults_delta,datanames,0);

    optionsp->out("    Results for the inclusion probabilities are also stored in file\n");
    optionsp->out("    " +  pathresults_delta + "\n");
    optionsp->out("\n");
    optionsp->out("\n");

	optionsp->out("    Rao-Blackwellised inclusion probability:\n");
    optionsp->out("\n");

    FC_delta.outresults(out_stata,out_R,out_R2BayesX,pathresults_delta_prob);
    FC_delta.outresults_help(out_stata,out_R,pathresults_delta_prob,datanames,1);

    optionsp->out("    Results for the Rao-Blackwellised inclusion probabilities are also stored in file\n");
    optionsp->out("    " +  pathresults_delta_prob + "\n");
    optionsp->out("\n");
    optionsp->out("\n");


    // deltas
/*    ofstream ou(pathresults_delta.strtochar());

    ou << "pmean_delta " << "pmean_prob" << endl;
    ou << FC_delta.betamean(0,0) << " " << FC_delta.betamean(0,1) << endl;'*/

    optionsp->out("    Variance hyperparameter psi2: \n");
    optionsp->out("\n");
    FC_psi2.outresults(out_stata,out_R,out_R2BayesX,"");
//    optionsp->out("\n");
    FC_psi2.outresults_singleparam(out_stata,out_R,pathresults_psi2);
//    optionsp->out("    psi2: " + ST::doubletostring(FC_psi2.betamean(0,0),6)  + "\n");

//    optionsp->out("\n");
    optionsp->out("    Results for psi2 are also stored in file\n");
    optionsp->out("    " +  pathresults_psi2 + "\n");
    optionsp->out("\n");
    optionsp->out("\n");

    if(singleomega == false)
      {
      optionsp->out("    Inclusion probability parameter omega:\n");
      optionsp->out("\n");
      FC_omega.outresults_singleparam(out_stata,out_R,pathresults_omega);
      optionsp->out("    Results for the inclusion probability parameter omega are also stored in file\n");
      optionsp->out("    " +  pathresults_omega + "\n");
      optionsp->out("\n");
      optionsp->out("\n");
      }

    }

}


void FC_nonp_variance_varselection::get_samples(
   const ST::string & filename,ofstream & outg) const
  {
  FC_nonp_variance::get_samples(filename,outg);

  ST::string filename_delta = filename.substr(0,filename.length()-11) + "_delta_sample.raw";
  FC_delta.get_samples(filename_delta,outg);

  if(singleomega == false)
  {
    ST::string filename_omega = filename.substr(0,filename.length()-11) + "_omega_sample.raw";
    FC_omega.get_samples(filename_omega,outg);

  }
  ST::string filename_psi2 = filename.substr(0,filename.length()-11) + "_psi2_sample.raw";
  FC_psi2.get_samples(filename_psi2,outg);
  }

void FC_nonp_variance_varselection::outoptions(void)
  {
//  FC_nonp_variance::outoptions();
  optionsp->out("  IG prior for psi2 \n");

  optionsp->out("  Hyperparameter v1 for prior: " +
              ST::doubletostring(v1) + "\n" );

  optionsp->out("  Hyperparameter v2 for prior: " +
              ST::doubletostring(v2) + "\n" );

  optionsp->out("  Parameter r for spike: " +
              ST::doubletostring(r) + "\n" );

  if(singleomega==false)
    {
    optionsp->out("\n");
    optionsp->out("  Prior inclusion probability omega\n");
    optionsp->out("  Parameter a_omega: " +
                ST::doubletostring(a_omega) + "\n" );
    optionsp->out("  Parameter b_omega: " +
                ST::doubletostring(b_omega) + "\n\n" );

    }
  }


void FC_nonp_variance_varselection::reset(void)
  {

  FC_nonp_variance::reset();

  }

//-------------------------------------------------------------------------------------------------------------
//------------------------------------  FC_varselection_omega -------------------------------------------------
//-------------------------------------------------------------------------------------------------------------



FC_varselection_omega::FC_varselection_omega(void)
  {

  }


FC_varselection_omega::FC_varselection_omega(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
           const ST::string & t)
     : FC(o,t,1,1,"")

  {
  likep = lp;
  masterp = mp;

  setbeta(1,1,0.5);
  a_omega = 1;
  b_omega = 1;
  }


FC_varselection_omega::FC_varselection_omega(const FC_varselection_omega & m)
    : FC(FC(m))
   {
   FC_tau2s = m.FC_tau2s;
   a_omega = m.a_omega;
   b_omega = m.b_omega;
   likep = m.likep;
   masterp = m.masterp;
   }


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_varselection_omega & FC_varselection_omega::operator=(const FC_varselection_omega & m)
   {

   if (this==&m)
   	 return *this;
   FC::operator=(FC(m));
   FC_tau2s = m.FC_tau2s;
   a_omega = m.a_omega;
   b_omega = m.b_omega;
   likep = m.likep;
   masterp = m.masterp;
   return *this;
   }


  void FC_varselection_omega::update(void)
    {
        //update omega when having one single omega for each equation.
        double sumdelta=0;
        unsigned nfunc = FC_tau2s.size();
        unsigned j;
        for (j=0;j<nfunc;j++)
        {
            sumdelta +=  FC_tau2s[j]->FC_delta.beta(0,0);
        }

        //draw new omega
        beta(0,0) = randnumbers::rand_beta(a_omega+sumdelta,
                                          b_omega+nfunc-sumdelta);


        for (j=0;j<nfunc;j++)
        {
        FC_tau2s[j]->omega=beta(0,0);
        }


        FC::update();
    }

//bool FC_varselection_omega::posteriormode(void)
//  {
//
//  }
//
//
  void FC_varselection_omega::outoptions(void)
    {
    optionsp->out("  Prior inclusion probability omega\n");
    optionsp->out("  Parameter a_omega: " +
                ST::doubletostring(a_omega) + "\n" );

    optionsp->out("  Parameter b_omega: " +
                ST::doubletostring(b_omega) + "\n" );
    }


  void FC_varselection_omega::outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                  const ST::string & pathresults)
   {
      ST::string pathresults_omega = pathresults.substr(0,pathresults.length()-4) + "_omega.res";
//        optionsp->out("    Prior inclusion probability parameter omega: " + ST::doubletostring(FC::betamean(0,0),6)  + "\n");
//        optionsp->out("\n");
        FC::outresults(out_stata,out_R,out_R2BayesX,pathresults_omega);
        FC::outresults_singleparam(out_stata,out_R,pathresults_omega);

        optionsp->out("    Results for the prior inclusion probability parameter omega are also stored in file\n");
        optionsp->out("    " +  pathresults_omega + "\n");
        optionsp->out("\n");
   }

  void FC_varselection_omega::reset(void)
    {
        FC::reset();
    }


//  void FC_varselection_omega::read_options(vector<ST::string> & op,vector<ST::string> & vn)
//    {
//
//    }
//

  void FC_varselection_omega::get_samples(const ST::string & filename,ofstream & outg) const
   {
        ST::string filename_omega = filename.substr(0,filename.length()-11) + "_omega_sample.raw";
        FC::get_samples(filename_omega,outg);

   }

//------------------------------------------------------------------------------
//------------------------ CLASS: FC_tensor_omega ------------------------------
//------------------------------------------------------------------------------

  // DEFAULT CONSTRUCTOR

  FC_tensor_omega::FC_tensor_omega(void){}

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // nro  : number of different values for omega
  // fp   : file path for storing sampled parameters

  FC_tensor_omega::FC_tensor_omega(DESIGN_userdefined_tensor * d,
           FC_nonp * f,
           GENERAL_OPTIONS * o,const ST::string & t,
           const unsigned & nro, const ST::string & fp)
           : FC(o,t,1,1,fp)
    {
    nromega=d->nromega;
    unsigned i;
    omegas = datamatrix(nromega,1,0.0);
    for(i=0; i<nromega; i++)
      omegas(i,0) = d->omegas[i];
    omegaindex = (unsigned)((nromega-1)/2);
    dut = d;
    fcn = f;
    logdets = dut->logdets;
    }
  // COPY CONSTRUCTOR

  FC_tensor_omega::FC_tensor_omega(const FC_tensor_omega & m)
         : FC(FC(m))
    {
    omegas = m.omegas;
    nromega = m.nromega;
    omegaindex = m.omegaindex;
    dut = m.dut;
    logdets = m.logdets;
    fcn = m.fcn;
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_tensor_omega & FC_tensor_omega::operator=(const FC_tensor_omega & m)
    {
    if (this==&m)
    	 return *this;
    FC::operator=(FC(m));
    omegas = m.omegas;
    nromega = m.nromega;
    omegaindex = m.omegaindex;
    dut = m.dut;
    logdets = m.logdets;
    fcn = m.fcn;
    return *this;
    }

  // FUNCTION: get_samples
  // TASK: stores the sampled parameters in ASCII format

  void FC_tensor_omega::get_samples(const ST::string & filename,ofstream & outg) const
    {
    if(nromega>1)
      {
    ST::string filename_omega = filename.substr(0,filename.length()-11) + "_omega_sample.raw";
        //FC::get_samples(filename_omega,outg);
    if ((nosamples == false) && (nosamplessave == false))
    {

    unsigned i,j,k;

    unsigned nrpar = beta.rows()*beta.cols();

    ofstream out(filename_omega.strtochar());
    assert(!out.fail());

    out << "intnr " << " ";
    if (beta.cols() > 1)
      {
      for (j=0;j<beta.rows();j++)
        for(k=0;k<beta.cols();k++)
          out << "b_" << (j+1) << "_" << (k+1) << " ";
      }
    else
      {
      for (j=0;j<nrpar;j++)
        out << "b_" << (j+1) << " ";
      }

    out << endl;

    double * sampled_betap = sampled_beta.getV();
    for(i=0;i<optionsp->samplesize;i++)
      {
      out << (i+1) << " ";
      for (j=0;j<nrpar;j++,sampled_betap++)
        out << ST::doubletostring(omegas(*sampled_betap,0),4) << " ";

      out << endl;
      }

    out.close();

    optionsp->out(filename_omega + "\n");

    outg << "_d.infile using " << filename_omega << endl;
    ST::string pathps = filename_omega.substr(0,filename_omega.length()-4) + ".ps";
    outg << "_g.plotsample , outfile=" <<  pathps.strtochar() <<  " using _d" << endl;
    outg << endl;

    } // end: if (nosamples == false)
      }
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void FC_tensor_omega::update(void)
    {
    if(nromega>1)
      {
    datamatrix logprobs(nromega,1,0.0);
    datamatrix probs(nromega,1,0.0);
    double probsum=0;
    double maxlogprob;
    unsigned i;
    for(i=0; i<nromega; i++)
      {
      dut->omegaindex=i;
      logprobs(i,0) = 0.5*logdets(i,0)-0.5/(fcn->tau2)*(dut->penalty_compute_quadform(fcn->param));

/*cout << "i: " << i << "\n";
cout << "tau2: " << fcn->tau2 << "\n";
cout << "bKb: " << dut->penalty_compute_quadform(fcn->param) << "\n";
cout << "0.5*logdets(i,0): " << 0.5*logdets(i,0) << "\n";
cout << "-0.5/(fcn->tau2)*(dut->penalty_compute_quadform(fcn->param)): " << -0.5/(fcn->tau2)*(dut->penalty_compute_quadform(fcn->param)) << "\n\n";*/

      if(i==0)
        maxlogprob = logprobs(i,0);
      else
        {
        if(logprobs(i,0) > maxlogprob)
          maxlogprob = logprobs(i,0);
        }
      }
    for(i=0; i<nromega; i++)
      {
      probs(i,0) = exp(logprobs(i,0)-maxlogprob);
      probsum += probs(i,0);
      }
    for(i=0; i<nromega; i++)
      {
      probs(i,0) /= probsum;
      }
    datamatrix cumprobs(nromega,1,0.0);
    cumprobs(0,0) = probs(0,0);
    for(i=1; i<nromega; i++)
      cumprobs(i,0) = cumprobs(i-1,0) + probs(i,0);

/*    ofstream out1("c:/temp/userdefined_tensor/probs.raw");
    probs.prettyPrint(out1);
    out1.close();*/

    double u = uniform();
    omegaindex=0;
    while(cumprobs(omegaindex,0)<=u)
      omegaindex++;

    dut->omegaindex = omegaindex;

    beta(0,0) = omegaindex;
    acceptance++;
    FC::update();
      }
    }

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void FC_tensor_omega::outoptions(void)
   {

   }
  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void FC_tensor_omega::outresults(ofstream & out_stata, ofstream & out_R, ofstream & out_R2BayesX,
               const ST::string & pathresults)
    {
    if(nromega>1)
      {
    FC::outresults(out_stata, out_R, out_R2BayesX, "");
    datamatrix omegafreq(nromega,1,0.0);
    unsigned i;
    for(i=0; i<sampled_beta.rows(); i++)
       omegafreq((int)sampled_beta(i,0),0) += 1.0;

    optionsp->out("\n");
    optionsp->out("Frequencies of anisotropy values:\n");
    optionsp->out("\n");
    for(i=0; i<nromega; i++)
      optionsp->out("  " + ST::doubletostring(omegas(i,0),4) + ": " +
      ST::doubletostring(omegafreq(i,0)/((double)sampled_beta.rows()),3) + "\n");
    optionsp->out("\n");

    ST::string pathresults_omega = pathresults.substr(0,pathresults.length()-4) + "_omega.res";
    ofstream outomegas(pathresults_omega.strtochar());
    outomegas << "omega   freq" << endl;
    for(i=0; i<nromega; i++)
      {
      outomegas << ST::doubletostring(omegas(i,0),4) << "   ";
      outomegas << ST::doubletostring(omegafreq(i,0)/((double)sampled_beta.rows()),3);
      outomegas << endl;
      }
    outomegas.close();

    optionsp->out("\n");
    optionsp->out("    Results for the frequencies of anisotropy values are also stored in file\n");
    optionsp->out("    " +  pathresults_omega + "\n");
    optionsp->out("\n");
    optionsp->out("\n");
      }
    }

  void FC_tensor_omega::read_options(vector<ST::string> & op,vector<ST::string> & vn)
    {

    }

} // end: namespace MCMC




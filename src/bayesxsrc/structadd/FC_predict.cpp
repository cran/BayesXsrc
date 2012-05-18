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



#include "FC_predict.h"
#include "clstring.h"

//------------------------------------------------------------------------------
//---------------- CLASS: FC_predict implementation of member functions --------
//------------------------------------------------------------------------------


namespace MCMC
{

void FC_predict::compute_autocorr_all(const ST::string & path, unsigned lag,
                                    ofstream & outg) const
  {

  }


void FC_predict::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {

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
  13      samplemult
  14      constraints
  */



  }

FC_predict::FC_predict(void)
  {
  }


FC_predict::FC_predict(GENERAL_OPTIONS * o,DISTR * lp,const ST::string & t,
     const ST::string & fp, const ST::string & fpd, datamatrix & dm,
      vector<ST::string> & dn)
  : FC(o,t,1,1,fp)
  {
//  nosamples = true;
  MSE = noMSE;
  MSEparam = 0.5;
  likep = lp;
  designmatrix= dm;
  varnames = dn;
  if (likep->maindistribution == true)
    {
    setbeta(lp->nrobs,2,0);
    }
  else
    {
    setbeta(lp->nrobs,2,0);
    }

  if (likep->maindistribution == true)
    {
    FC_deviance = FC(o,"",2,1,fpd);
    }

  int d = betamean.rows();

  }


FC_predict::FC_predict(const FC_predict & m)
  : FC(FC(m))
  {
  MSE = m.MSE;
  MSEparam = m.MSEparam;
  likep = m.likep;
  designmatrix = m.designmatrix;
  varnames = m.varnames;
  FC_deviance = m.FC_deviance;
  deviance = m.deviance;
  deviancesat = m.deviancesat;
  }


const FC_predict & FC_predict::operator=(const FC_predict & m)
  {
  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  MSE = m.MSE;
  MSEparam = m.MSEparam;
  likep = m.likep;
  designmatrix = m.designmatrix;
  varnames = m.varnames;
  FC_deviance = m.FC_deviance;
  deviance = m.deviance;
  deviancesat = m.deviancesat;
  return *this;
  }


void  FC_predict::update(void)
  {

  likep->FCpredict_betamean = &betamean;

  get_predictor();

  acceptance++;

  FC::update();

  if (likep->maindistribution == true)
    {
    FC_deviance.beta(0,0) = deviance;
    FC_deviance.beta(1,0) = deviancesat;
    FC_deviance.acceptance++;
    FC_deviance.update();
    }

  }


void FC_predict::get_predictor(void)
  {

  unsigned i;
  double * betap = beta.getV();

  double * worklinp;
  if (likep->linpred_current==1)
    worklinp = likep->linearpred1.getV();
  else
    worklinp = likep->linearpred2.getV();

  deviance=0;
  deviancesat=0;
  double deviancehelp;
  double deviancesathelp;

  double * workresponse = likep->response_untransformed.getV();
  double * workweight = likep->weight.getV();
  double muhelp;
  double scalehelp=likep->get_scale();
  double weightone = 1;

  double predlik;

  for(i=0;i<likep->nrobs;i++,worklinp++,workresponse++,workweight++,betap++)
    {

    likep->compute_mu(worklinp,&muhelp);

    if (likep->maindistribution == true)
      {
      likep->compute_deviance(workresponse,workweight,&muhelp,&deviancehelp,
      &deviancesathelp,&scalehelp);

      deviance+=deviancehelp;
      deviancesat+=deviancesathelp;

//      if (*workweight==0)
//        predlik = likep->loglikelihood(workresponse,worklinp,&weightone);
//      else
//        predlik = likep->loglikelihood(workresponse,worklinp,workweight);

      }

    *betap = *worklinp;
    betap++;
    *betap = muhelp;
//    if (likep->maindistribution == true)
//      {
//      betap++;
//      *betap = predlik;
//      }

    }

  }


bool FC_predict::posteriormode(void)
  {

  get_predictor();

  posteriormode_betamean();

  return true;
  }


void FC_predict::outoptions(void)
  {

  }



void FC_predict::outresults_DIC(const ST::string & pathresults)
  {

  ST::string pathresultsdic = pathresults.substr(0,pathresults.length()-4) + "_DIC.res";
  ofstream out(pathresultsdic.strtochar());

  double deviance2=0;
  double deviance2_sat=0;


  double scalehelp = likep->get_scalemean();


  double devhelp;
  double devhelp_sat;

  double mu_meanlin;
  double * workmeanlin = betamean.getV();
  double * workresponse = likep->response_untransformed.getV();
  double * workweight = likep->weight.getV();

  unsigned i;
  for (i=0;i<likep->nrobs;i++,workmeanlin+=beta.cols(),workresponse++,workweight++)
    {

    likep->compute_mu(workmeanlin,&mu_meanlin);

    likep->compute_deviance(workresponse,workweight,&mu_meanlin,
    &devhelp,&devhelp_sat,&scalehelp);

    deviance2 += devhelp;
    deviance2_sat += devhelp_sat;
    }


  double devhelpm = FC_deviance.betamean(0,0);

  unsigned d;
  if (devhelpm > 1000000000)
    d = 14;
  else if (devhelpm > 1000000)
    d = 11;
  else
    d = 8;

  out << "deviance   pd dic" << endl;

  optionsp->out("  ESTIMATION RESULTS FOR THE DIC: \n",true);
  optionsp->out("\n");

  optionsp->out("    DIC based on the unstandardized deviance\n");
  optionsp->out("\n");

  optionsp->out("    Deviance(bar_mu):           " +
  ST::doubletostring(deviance2,d) + "\n");
  out << deviance2 << "   ";

  optionsp->out("    pD:                         " +
  ST::doubletostring(devhelpm-deviance2,d) + "\n");
  out << (devhelpm-deviance2) << "   ";

  optionsp->out("    DIC:                        " +
  ST::doubletostring(2*devhelpm-deviance2,d) + "\n");
  optionsp->out("\n");
  out << (2*devhelpm-deviance2) << "   " << endl;

  optionsp->out("    DIC based on the saturated deviance\n");
  optionsp->out("\n");


  double devhelpm_sat = FC_deviance.betamean(1,0);

  optionsp->out("    Deviance(bar_mu):           " +
  ST::doubletostring(deviance2_sat,d) + "\n");

  optionsp->out("    pD:                         " +
  ST::doubletostring(devhelpm_sat-deviance2_sat,d) + "\n");

  optionsp->out("    DIC:                        " +
  ST::doubletostring(2*devhelpm_sat-deviance2_sat,d) + "\n");

  optionsp->out("\n");

  }


void FC_predict::outresults_deviance(void)
    {

    ST::string l1 = ST::doubletostring(optionsp->lower1,4);
    ST::string l2 = ST::doubletostring(optionsp->lower2,4);
    ST::string u1 = ST::doubletostring(optionsp->upper1,4);
    ST::string u2 = ST::doubletostring(optionsp->upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');


    ST::string meanstr = "    Mean:          ";
    unsigned l_meanstr = meanstr.length();

    ST::string stdstr =  "    Std. Dev:      ";
    unsigned l_stdstr = stdstr.length();

    ST::string l1str = "    " + l1 + "% Quantile: ";
    unsigned l_l1str = l1str.length();

    ST::string l2str = "    " + l2 + "% Quantile: ";
    unsigned l_l2str = l2str.length();

    ST::string medianstr = "    50% Quantile: ";
    unsigned l_medianstr = medianstr.length();

    ST::string u1str = "    " + u1 + "% Quantile: ";
    unsigned l_u1str = u1str.length();

    ST::string u2str = "    " + u2 + "% Quantile: ";
    unsigned l_u2str = u2str.length();


    optionsp->out("  ESTIMATION RESULT FOR THE DEVIANCE: \n",true);
    optionsp->out("\n");

    optionsp->out("    Unstandardized Deviance (-2*Loglikelihood(y|mu))\n");
    optionsp->out("\n");

    double devhelpm = FC_deviance.betamean(0,0);
    double devhelp;

    unsigned d;
    if (devhelpm > 1000000000)
      d = 14;
    else if (devhelpm > 1000000)
      d = 11;
    else
      d = 8;

    optionsp->out(meanstr + ST::string(' ',20-l_meanstr) +
    ST::doubletostring(devhelpm,d) + "\n");


    devhelp = sqrt(FC_deviance.betavar(0,0));
    optionsp->out(stdstr + ST::string(' ',20-l_stdstr) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l1_lower(0,0);
    optionsp->out(l1str +  ST::string(' ',20-l_l1str) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l2_lower(0,0);
    optionsp->out(l2str +  ST::string(' ',20-l_l2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu50(0,0);
    optionsp->out(medianstr +  ST::string(' ',20-l_medianstr) +
    ST::doubletostring(devhelp,d) +  "\n");


    devhelp = FC_deviance.betaqu_l2_upper(0,0);
    optionsp->out(u1str +  ST::string(' ',20-l_u1str) +
    ST::doubletostring(devhelp,d) +  "\n");


    devhelp = FC_deviance.betaqu_l1_upper(0,0);
    optionsp->out(u2str +  ST::string(' ',20-l_u2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    optionsp->out("\n");


    optionsp->out("  Saturated Deviance (-2*Loglikelihood(y|mu) + 2*Loglikelihood(y|mu=y))\n");
    optionsp->out("\n");

    devhelpm = FC_deviance.betamean(1,0);

    optionsp->out(meanstr + ST::string(' ',20-l_meanstr) +
    ST::doubletostring(devhelpm,d) + "\n");

    devhelp = sqrt(FC_deviance.betavar(1,0));
    optionsp->out(stdstr + ST::string(' ',20-l_stdstr) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l1_lower(1,0);;
    optionsp->out(l1str +  ST::string(' ',20-l_l1str) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l2_lower(1,0);
    optionsp->out(l2str +  ST::string(' ',20-l_l2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu50(1,0);
    optionsp->out(medianstr +  ST::string(' ',20-l_medianstr) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l2_upper(1,0);
    optionsp->out(u1str +  ST::string(' ',20-l_u1str) +
    ST::doubletostring(devhelp,d) +  "\n");

    devhelp = FC_deviance.betaqu_l1_upper(1,0);
    optionsp->out(u2str +  ST::string(' ',20-l_u2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    optionsp->out("\n");

    optionsp->out("\n");

    }


ST::string FC_predict::getloss(void)
  {
  if (MSE==quadraticMSE)
    return "quadratic";
  else if (MSE==checkMSE)
    return "check";
  else
    return " ";
  }

void FC_predict::compute_MSE(const ST::string & pathresults)
  {

  unsigned nrobs = designmatrix.rows();
  unsigned nrzeroweights;
  double meanmse;
  double meanmse_zeroweight;

  likep->compute_MSE_all(betamean, meanmse, meanmse_zeroweight,nrzeroweights,
                         MSE,MSEparam);

  ST::string h;
  optionsp->out("  EMPIRICAL MSE: \n",true);
  optionsp->out("\n");
  ST::string t = getloss();

  optionsp->out("    Loss function: " + t +  "\n");
  h = ST::doubletostring(meanmse_zeroweight+meanmse,10);
  optionsp->out("    sum MSE (all observations):            " +
  h +   "\n");
  optionsp->out("    sum MSE (zero weight observations):    " +
  ST::doubletostring(meanmse_zeroweight,10) +  "\n");
  optionsp->out("    sum MSE (nonzero weight observations): " +
  ST::doubletostring(meanmse,10) +  "\n");
  optionsp->out("\n");

  ST::string pathmse = pathresults.substr(0,pathresults.length()-4) + "_MSE.res";
  ofstream out(pathmse.strtochar());
  out << "nrobs  nrzeroweights  sum_MSE  sum_MSE_zeroweights   sum_MSE_nonzeroweights" << endl;
  out << nrobs  << "  "
      << nrzeroweights << "  "
      <<  (meanmse_zeroweight+meanmse)
      << "  " << meanmse_zeroweight
      << "  " << meanmse << endl;

  optionsp->out("    Results for the MSE are also stored in file\n");
  optionsp->out("    " + pathmse + "\n");
  optionsp->out("\n");

  }

/*
  for(i=0;i<nrobs;i++,responsep++,weightp++,linpredp+=2)
    if (*weightp==0)
      {
      meanmse_zeroweight += likep->compute_MSE(responsep,weightp,linpredp,MSE,MSEparam);
      nrzeroweights++;
      }
    else
      meanmse += likep->compute_MSE(responsep,weightp,linpredp,MSE,MSEparam);
*/

//void FC_predict::get_samples(const ST::string & filename,ofstream & outg) const
//  {
//  double f=0;
//  }

void FC_predict::outresults(ofstream & out_stata, ofstream & out_R,
                            const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    FC::outresults(out_stata,out_R,"");

    if (likep->maindistribution == true)
      {
      FC_deviance.outresults(out_stata,out_R,"");
      }

    optionsp->out("  PREDICTED VALUES: \n",true);
    optionsp->out("\n");

    optionsp->out("    Results for the predictor, mean are stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    if ((likep->maindistribution == true) && (MSE != noMSE))
      {
      compute_MSE(pathresults);
      }

    ofstream outres(pathresults.strtochar());

    optionsp->out("\n");

    unsigned i,j;

    ST::string l1 = ST::doubletostring(optionsp->lower1,4);
    ST::string l2 = ST::doubletostring(optionsp->lower2,4);
    ST::string u1 = ST::doubletostring(optionsp->upper1,4);
    ST::string u2 = ST::doubletostring(optionsp->upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');

    outres << "intnr" << "   ";

    for (i=0;i<varnames.size();i++)
      outres << varnames[i] << "   ";

    outres << "pmean_pred   ";

    if (optionsp->samplesize > 1)
      {
      outres << "pqu"  << l1  << "_pred   ";
      outres << "pqu"  << l2  << "_pred   ";
      outres << "pmed_pred   ";
      outres << "pqu"  << u1  << "_pred   ";
      outres << "pqu"  << u2  << "_pred   ";
      }

    outres << "pmean_mu   ";

    if (optionsp->samplesize > 1)
      {
      outres << "pqu"  << l1  << "_mu   ";
      outres << "pqu"  << l2  << "_mu   ";
      outres << "pmed_mu   ";
      outres << "pqu"  << u1  << "_mu   ";
      outres << "pqu"  << u2  << "_mu   ";
      }

/*
    if (likep->maindistribution == true)
      {

      outres << "pmean_predl   ";

      if (optionsp->samplesize > 1)
        {
        outres << "pqu"  << l1  << "_predl   ";
        outres << "pqu"  << l2  << "_predl   ";
        outres << "pmed_predl   ";
        outres << "pqu"  << u1  << "_predl   ";
        outres << "pqu"  << u2  << "_predl   ";
        }

      }
*/

    outres << endl;

    double * workmean = betamean.getV();
    double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
    double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
    double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
    double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
    double * workbetaqu50 = betaqu50.getV();


    for(i=0;i<designmatrix.rows();i++,
          workmean++,
          workbetaqu_l1_lower_p++,
          workbetaqu_l2_lower_p++,
          workbetaqu50++,
          workbetaqu_l1_upper_p++,
          workbetaqu_l2_upper_p++)
      {

      outres << (i+1) << "   ";

      for (j=0;j<designmatrix.cols();j++)
        outres << designmatrix(i,j) << "   ";

      outres << *workmean << "   ";

      if (optionsp->samplesize > 1)
        {
        outres << *workbetaqu_l1_lower_p << "   ";
        outres << *workbetaqu_l2_lower_p << "   ";
        outres << *workbetaqu50 << "   ";
        outres << *workbetaqu_l2_upper_p << "   ";
        outres << *workbetaqu_l1_upper_p << "   ";
        }

      // mu
      workmean++;
      workbetaqu_l1_lower_p++;
      workbetaqu_l2_lower_p++;
      workbetaqu50++;
      workbetaqu_l1_upper_p++;
      workbetaqu_l2_upper_p++;

      outres << *workmean << "   ";

      if (optionsp->samplesize > 1)
        {
        outres << *workbetaqu_l1_lower_p << "   ";
        outres << *workbetaqu_l2_lower_p << "   ";
        outres << *workbetaqu50 << "   ";
        outres << *workbetaqu_l2_upper_p << "   ";
        outres << *workbetaqu_l1_upper_p << "   ";
        }
      // end mu

/*
      // predictive likelihood
      if (likep->maindistribution == true)
        {
        workmean++;
        workbetaqu_l1_lower_p++;
        workbetaqu_l2_lower_p++;
        workbetaqu50++;
        workbetaqu_l1_upper_p++;
        workbetaqu_l2_upper_p++;

        outres << *workmean << "   ";

        if (optionsp->samplesize > 1)
          {
          outres << *workbetaqu_l1_lower_p << "   ";
          outres << *workbetaqu_l2_lower_p << "   ";
          outres << *workbetaqu50 << "   ";
          outres << *workbetaqu_l2_upper_p << "   ";
          outres << *workbetaqu_l1_upper_p << "   ";
          }
        }
      // end predictive likelihood
*/
      outres << endl;

     }

    if (likep->maindistribution == true)
      {
      outresults_deviance();
      outresults_DIC(pathresults);
      }

    }   // end if (pathresults.isvalidfile() != 1)

  }


void FC_predict::reset(void)
  {

  }



} // end: namespace MCMC




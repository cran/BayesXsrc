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

#include "distr_gamlss.h"

#if defined(BayesX_gsl_included)
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#endif

namespace MCMC
{


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_negbin_delta --------------------------
//------------------------------------------------------------------------------

DISTR_negbin_delta::DISTR_negbin_delta(GENERAL_OPTIONS * o,
                                       const datamatrix & r,
                                       double & ss, int & strmax,
                                       int & sts,
                                       bool & sl,
                                       const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {

  predictor_name = "delta";
  outpredictor = true;
  outexpectation = false;


  family = "Negative_Binomial - delta";

  double responsemax = response.max(0);

  stopsum = ss;
  stoprmax = strmax;
  if (stoprmax < responsemax)
    stoprmax = responsemax;
  nrbetween = sts;

  slow=sl;

  E_dig_y_delta_m = datamatrix(nrobs,1,0);
  E_trig_y_delta_m = datamatrix(nrobs,1,0);

  linpredminlimit=-10;
  linpredmaxlimit=10;
  }


DISTR_negbin_delta::DISTR_negbin_delta(const DISTR_negbin_delta & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  E_dig_y_delta = nd.E_dig_y_delta;
  E_trig_y_delta = nd.E_trig_y_delta;
  delta = nd.delta;
  log_delta_div_delta_plus_mu = nd.log_delta_div_delta_plus_mu;
  lngamma_delta = nd.lngamma_delta;
  delta_plus_mu = nd.lngamma_delta;

  stopsum = nd.stopsum;
  stoprmax = nd.stoprmax;
  nrbetween = nd.nrbetween;

  slow = nd.slow;

  E_dig_y_delta_m = nd.E_dig_y_delta_m;
  E_trig_y_delta_m = nd.E_trig_y_delta_m;
  Ep = nd.Ep;
  Ep_trig = nd.Ep_trig;
  }


const DISTR_negbin_delta & DISTR_negbin_delta::operator=(
                            const DISTR_negbin_delta & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  E_dig_y_delta = nd.E_dig_y_delta;
  E_trig_y_delta = nd.E_trig_y_delta;
  delta = nd.delta;
  log_delta_div_delta_plus_mu = nd.log_delta_div_delta_plus_mu;
  lngamma_delta = nd.lngamma_delta;
  delta_plus_mu = nd.lngamma_delta;

  stopsum = nd.stopsum;
  stoprmax = nd.stoprmax;
  nrbetween = nd.nrbetween;

  slow = nd.slow;

  E_dig_y_delta_m = nd.E_dig_y_delta_m;
  E_trig_y_delta_m = nd.E_trig_y_delta_m;
  Ep = nd.Ep;
  Ep_trig = nd.Ep_trig;
  return *this;
  }



double DISTR_negbin_delta::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_negbin_delta::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[0]);
  }

double DISTR_negbin_delta::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double delta;
  double l;

  if (*linpred <= linpredminlimit)
    delta = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    delta = exp(linpredmaxlimit);
  else
    delta = exp(*linpred);

  double log_mu_plus_delta = log((*worktransformlin[0]) + delta);

  if (*response==0)
    l = -delta*(log_mu_plus_delta - log(delta));
  else
    {
    double resp_plus_delta = (*response) + delta;

    l = randnumbers::lngamma_exact(resp_plus_delta) -
        randnumbers::lngamma_exact(delta) -
        (resp_plus_delta)* log_mu_plus_delta  +
        delta*log(delta);
    }

  modify_worklin();

  return l;

  }


void DISTR_negbin_delta::compute_expectation(void)
  {

  int k=1;
  double k_delta;
  double kplus1;
  double psum;

  double L = exp(delta*log_delta_div_delta_plus_mu);
  E_dig_y_delta = randnumbers::digamma_exact(delta)*L;
  E_trig_y_delta = randnumbers::trigamma_exact(delta)*L;

  psum = L;

  while ((psum < stopsum) && (k <=stoprmax))
    {
    k_delta = k + delta;
    kplus1 = k + 1;

    L = exp(randnumbers::lngamma_exact(k_delta) -
            randnumbers::lngamma_exact(kplus1) -
            lngamma_delta + delta*log_delta_div_delta_plus_mu +
            k* log(*worktransformlin[0]/delta_plus_mu)
           );

    psum += L;

    E_dig_y_delta += randnumbers::digamma_exact(k_delta)*L;

    E_trig_y_delta += randnumbers::trigamma_exact(k_delta)*L;

    k++;
    }

  E_dig_y_delta -=  randnumbers::digamma_exact(delta);

  E_trig_y_delta -= randnumbers::trigamma_exact(delta);

  E_dig_y_delta *=  delta;

  E_trig_y_delta *= delta*delta;

  *Ep = E_dig_y_delta;
  *Ep_trig = E_trig_y_delta;

  }

double DISTR_negbin_delta::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {
  if (*weight == 0)
    {
    if (counter==0)
      {
      set_worklin();
      Ep = E_dig_y_delta_m.getV();
      Ep_trig = E_trig_y_delta_m.getV();
      }

    *workingweight = 0;

    modify_worklin();

    return 0;
    }
  else
    {

    double l=0;
    compute_iwls_wweightschange_weightsone(response,linpred, workingweight,
                                           workingresponse, l,
                                           like);
     *workingweight *= (*weight);

     return (*weight)*l;

    }

  }

void DISTR_negbin_delta::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    Ep = E_dig_y_delta_m.getV();
    Ep_trig = E_trig_y_delta_m.getV();
    }

  if (*linpred <= linpredminlimit)
    delta = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    delta = exp(linpredmaxlimit);
  else
    delta = exp(*linpred);

  delta_plus_mu = delta + (*worktransformlin[0]);

  log_delta_div_delta_plus_mu = log(delta/delta_plus_mu);

  lngamma_delta = randnumbers::lngamma_exact(delta);

  double delta_plus_response = delta + (*response);

  double nu = delta*(randnumbers::digamma_exact(delta_plus_response) -
                    randnumbers::digamma_exact(delta) +
                     log_delta_div_delta_plus_mu +
                    (*worktransformlin[0]-(*response))/delta_plus_mu);

  if ((optionsp->nriter < 1) ||
      slow ||
      (optionsp->nriter % nrbetween == 0)
      )
    compute_expectation();
  else
    {
    E_dig_y_delta = (*Ep);
    E_trig_y_delta = (*Ep_trig);
    }

  *workingweight = -delta*(log_delta_div_delta_plus_mu
                   +(*worktransformlin[0])/delta_plus_mu)
                   -E_dig_y_delta-E_trig_y_delta;

  if (*workingweight <= 0)
    *workingweight = 0.0001;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    double resp_plus_delta = (*response) + delta;

    like += randnumbers::lngamma_exact(resp_plus_delta) -
            lngamma_delta -
            resp_plus_delta*log(delta_plus_mu) +
            delta*log(delta);

    }

  modify_worklin();
  Ep++;
  Ep_trig++;

  }


void DISTR_negbin_delta::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (delta): exponential\n");
  optionsp->out("  Stop criteria for approximating expected values\n");
  optionsp->out("  in working weights of delta equation:\n");
  optionsp->out("    cumulative probability:"  + ST::doubletostring(stopsum) +  "\n");
  optionsp->out("    Maximum values:"  + ST::inttostring(stoprmax) +  "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_negbin_delta::update_end(void)
  {
  DISTR_gamlss::update_end();
  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_negbin_mu -----------------------------
//------------------------------------------------------------------------------

void DISTR_negbin_mu::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {
        if (*workresp != int(*workresp))
          {
          errors=true;
          errormessages.push_back("ERROR: response must be integer values\n");
          }

        if (*workresp < 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


DISTR_negbin_mu::DISTR_negbin_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {

  predictor_name = "mu";
  outpredictor = true;
  outexpectation = true;

  family = "Negative_Binomial - mu";

  linpredminlimit=-10;
  linpredmaxlimit=15;

  check_errors();
  }


DISTR_negbin_mu::DISTR_negbin_mu(const DISTR_negbin_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  }


const DISTR_negbin_mu & DISTR_negbin_mu::operator=(
                            const DISTR_negbin_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_negbin_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_delta
   // *linpred[1] = eta_mu

   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double delta = exp(*linpred[0]);
     double mu = exp(*linpred[1]);
     double delta_plus_mu = delta + mu;
     double resp_plus_one = (*response[1]) + 1;

     double l;

     if (*response[1]==0)
       l = - randnumbers::lngamma_exact(resp_plus_one)
           + delta*log(delta/delta_plus_mu);
     else
       {
       double delta_plus_response = delta+(*response[1]);

       l = randnumbers::lngamma_exact(delta_plus_response) -
           randnumbers::lngamma_exact(resp_plus_one) -
           randnumbers::lngamma_exact(delta) +
           delta*log(delta/delta_plus_mu)+
           (*response[1])*log(mu/delta_plus_mu);

       }
    *deviance = -2*l;
    }

  }


double DISTR_negbin_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_negbin_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]);
  }


double DISTR_negbin_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
  //  double p =  (*param[1])/((*param[0])+(*param[1]));
  //  double r = (*param[0]);
  //  double u = gsl_ran_negative_binomial_pdf(*response[1], p, r);
   // return u;
      return 0;
    }

double DISTR_negbin_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
  //  double p =  (*param[1])/((*param[0])+(*param[1]));
  //  double r = (*param[0]);
  //  double u = gsl_cdf_negative_binomial_P(*response[1], p, r);
  //  return u;
      return 0;
    }


double DISTR_negbin_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of delta equation
  // *worktransformlin[0] = exp(eta_delta);

  if (counter==0)
    {
    set_worklin();
    }

  double mu;
  if (*linpred <= linpredminlimit)
    mu = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    mu = exp(linpredmaxlimit);
  else
    mu = exp(*linpred);

  double l;
  if (*response==0)
     l = - (*worktransformlin[0])*log((*worktransformlin[0])+mu);
  else
     l = - ((*worktransformlin[0]) + (*response))*
           log((*worktransformlin[0])+mu) +(*response)*log(mu);

  modify_worklin();

  return l;

  }


void DISTR_negbin_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of delta equation
  // *worktransformlin[0] = exp(eta_delta);

  if (counter==0)
    set_worklin();

   double mu;
  if (*linpred <= linpredminlimit)
    mu = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    mu = exp(linpredmaxlimit);
  else
    mu = exp(*linpred);

  double delta_plus_mu = (*worktransformlin[0]) + mu;

  double nu = (*worktransformlin[0])*((*response)-mu)/delta_plus_mu;

  *workingweight = (*worktransformlin[0])*mu/delta_plus_mu;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    if (*response==0)
      like -= (*worktransformlin[0])*log(delta_plus_mu);
    else
    like += -((*worktransformlin[0])+(*response))*log(delta_plus_mu) +
             (*response)*log(mu);

    }

  modify_worklin();

  }


void DISTR_negbin_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = exp(*linpred[predstart_mumult+1]);
  }


void DISTR_negbin_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_negbin_mu::update_end(void)
  {
  DISTR_gamlss::update_end();
  }


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_zip_cloglog_pi --------------------------
//------------------------------------------------------------------------------


DISTR_zip_cloglog_pi::DISTR_zip_cloglog_pi(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {

  predictor_name = "pi";
  outpredictor = true;
  outexpectation = false;

  family = "Zero_Inflated_Poisson - pi";

  helpmat1 = datamatrix(nrobs, 1, 1.0-exp(-1.0));

  linpredminlimit=-10;
  linpredmaxlimit=10;

  }


DISTR_zip_cloglog_pi::DISTR_zip_cloglog_pi(const DISTR_zip_cloglog_pi & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_zip_cloglog_pi & DISTR_zip_cloglog_pi::operator=(
                            const DISTR_zip_cloglog_pi & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_zip_cloglog_pi::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_zip_cloglog_pi::compute_param_mult(vector<double *>  linpred,double * param)
  {
  double el = exp(*linpred[0]);
  *param = exp(-exp(el));
  }

double DISTR_zip_cloglog_pi::loglikelihood_weightsone(double * response,
                                                      double * linpred)
  {

  if (counter==0)
    set_worklin();

  double explinpi = exp(*linpred);
  double oneminuspi = 1 - exp(-explinpi);
  double pi = 1-oneminuspi;
  double expminuslambda = exp(-(*worktransformlin[0]));
  double denompart = pi+oneminuspi*expminuslambda;

  modify_worklin();

  if (*response == 0)
    return log(denompart);
  else
    return log(oneminuspi);

  }


void DISTR_zip_cloglog_pi::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    set_worklin();

  double explinpi = exp(*linpred);
  double oneminuspi = 1 - exp(-explinpi);
  double pi = 1-oneminuspi;
  double expminuslambda = exp(-(*worktransformlin[0]));
  double denompart = pi+oneminuspi*expminuslambda;
  double denom =  denompart*oneminuspi;
  double explinpi_pi = explinpi*pi;

  double nu = explinpi_pi/oneminuspi;
  if (*response == 0)
    nu -= explinpi_pi/denom;

  *workingweight =  pow(explinpi,2)*pow(pi,2)*(1-expminuslambda)/denom;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    if (*response == 0)
      like += log(denompart);
    else
      like += log(oneminuspi);

    }

  modify_worklin();

  }


void DISTR_zip_cloglog_pi::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (pi): complementary log log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_zip_cloglog_pi::update_end(void)
  {

  // helpmat1 stores 1-pi

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();


  double * ppi = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,ppi++,worklin++)
    {
    *ppi = 1-exp(-exp(*worklin));
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_zip_cloglog_mu --------------------------
//------------------------------------------------------------------------------

void DISTR_zip_cloglog_mu::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {
        if (*workresp != int(*workresp))
          {
          errors=true;
          errormessages.push_back("ERROR: response must be integer values\n");
          }

        if (*workresp < 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


DISTR_zip_cloglog_mu::DISTR_zip_cloglog_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {

  predictor_name = "lambda";
  outpredictor = true;
  outexpectation = true;

  family = "Zero_Inflated_Poisson - lambda";


  linpredminlimit=-10;
  linpredmaxlimit=15;

  check_errors();
  }


DISTR_zip_cloglog_mu::DISTR_zip_cloglog_mu(const DISTR_zip_cloglog_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  }


const DISTR_zip_cloglog_mu & DISTR_zip_cloglog_mu::operator=(
                            const DISTR_zip_cloglog_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_zip_cloglog_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

  double l;
  double explinpi = exp(*linpred[0]);
  double pi = exp(-explinpi);
  double lambda = exp(*linpred[1]);

  if (*response[1]==0)
    {
    l=  log(pi+(1-pi)*exp(-lambda));
    }
  else // response > 0
    {
    double help1 = *response[1]+1;
    l= log(1-pi) + (*response[1])*(*linpred[1])- lambda
       - randnumbers::lngamma_exact(help1);
    }

  *deviance = -2*l;

  }


double DISTR_zip_cloglog_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_zip_cloglog_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]);
  }

double DISTR_zip_cloglog_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_zip_cloglog_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_zip_cloglog_mu::loglikelihood_weightsone(double * response,
                                                      double * linpred)
  {

  // *worklin[0] = linear predictor of pi equation
  // *worktransformlin[0] = 1-pi = 1-exp(-exp(eta_pi));

  if (counter==0)
    {
    set_worklin();
    }

  double lambda;
  double expminuslambda;

  if (*linpred <= linpredminlimit)
    lambda  = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    lambda  = exp(linpredmaxlimit);
  else
    lambda = exp(*linpred);

  expminuslambda = exp(-lambda);

  double denom = 1-(*worktransformlin[0])+(*worktransformlin[0])*expminuslambda;

  modify_worklin();

  if (*response==0)
    return log(denom);
  else
    return (*response)*(*linpred)-lambda;

  }


void DISTR_zip_cloglog_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of pi equation
  // *worktransformlin[0] = 1-pi = 1-exp(-exp(eta_pi));

  if (counter==0)
    {
    set_worklin();
    }

  double lambda;
  double expminuslambda;

  if (*linpred <= linpredminlimit)
    lambda  = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    lambda  = exp(linpredmaxlimit);
  else
    lambda = exp(*linpred);

  expminuslambda = exp(-lambda);

  double pi = 1-(*worktransformlin[0]);
  double denom = pi+(*worktransformlin[0])*expminuslambda;

  double nu = (*response) - lambda;
  if (*response == 0)
    nu += pi*lambda/denom;

  *workingweight = (lambda* (*worktransformlin[0])*(denom-expminuslambda*lambda*pi))/denom;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    if (*response==0)
      like += log(denom);
    else // response > 0
      like += (*response)*(*linpred)-lambda;

    }

  modify_worklin();

  }


void DISTR_zip_cloglog_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double el = exp(-exp(*linpred[predstart_mumult]));

  *mu = (1-el)*exp(*linpred[predstart_mumult+1]);
  }


void DISTR_zip_cloglog_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (lambda): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_zip_cloglog_mu::update_end(void)
  {
  DISTR_gamlss::update_end();
  }


//------------------------------------------------------------------------------
//----------------------------- CLASS DISTR_gamlss -----------------------------
//------------------------------------------------------------------------------


DISTR_gamlss::DISTR_gamlss(GENERAL_OPTIONS * o, const datamatrix & r,
                           unsigned nrdistr,
                           const datamatrix & w)
  : DISTR(o,r,w)

  {

  gamlss = true;

  predictor_name = "mu";

  predict_mult = true;

  if (check_weightsone() && !optionsp->copula  && !optionsp->samplesel)
    {
    wtype = wweightschange_weightsone;
    weightsone = true;
    }
  else
    {
    wtype = wweightschange_weightsneqone;
    weightsone = false;
    }


  helpmat1 = datamatrix(nrobs,1,1);

  worklin = vector<double*>(nrdistr);
  worktransformlin = vector<double*>(nrdistr);

  updateIWLS = true;

  }


const DISTR_gamlss & DISTR_gamlss::operator=(
const DISTR_gamlss & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  worklin = nd.worklin;
  worktransformlin = nd.worktransformlin;
  return *this;
  }


DISTR_gamlss::DISTR_gamlss(const DISTR_gamlss & nd)
   : DISTR(DISTR(nd))
  {
  worklin = nd.worklin;
  worktransformlin = nd.worktransformlin;
  }


void DISTR_gamlss::outoptions(void)
  {
  DISTR::outoptions();
  }


void DISTR_gamlss::set_worklin(void)
  {

  unsigned i;
  for (i=0;i<worklin.size();i++)
    {

    if (distrp[i]->linpred_current==1)
      worklin[i] = distrp[i]->linearpred1.getV();
    else
      worklin[i] = distrp[i]->linearpred2.getV();

    worktransformlin[i] = distrp[i]->helpmat1.getV();
    }

  }


void DISTR_gamlss::modify_worklin(void)
  {

  if (counter<nrobs-1)
    {
    counter++;
    unsigned i;
    for (i=0;i<worklin.size();i++)
      {
      worklin[i]++;
      worktransformlin[i]++;
      }
    }
  else
    {
    counter=0;
    }

  }


double DISTR_gamlss::get_intercept_start(void)
  {
  return 0;
  }


double DISTR_gamlss::loglikelihood(double * response, double * linpred,
                                         double * weight)
  {
  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (*weight == 0 && !optionsp->copula)
    {
    if (counter==0)
      {
      set_worklin();
      }

    modify_worklin();

    return 0;
    }
  else
    {
    double ll = loglikelihood_weightsone(response, linpred);
    if (*weight == 0.0)
      {
      return 0.0;
      }
    return (*weight) * ll;
    }
  }


double DISTR_gamlss::loglikelihood_weightsone(double * response,
                                                    double * linpred)
  {
  return 0;
  }


void DISTR_gamlss::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  }


void DISTR_gamlss::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[0]);
  }

void DISTR_gamlss::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
  //  *deviance = -2*l;
  }


double DISTR_gamlss::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {
  if (*weight == 0 && !optionsp->copula)
    {
    if (counter==0)
      {
      set_worklin();
      }

    *workingweight = 0;

    modify_worklin();

    return 0.0;
    }
  else
    {
    double l = 0.0;
    compute_iwls_wweightschange_weightsone(response,linpred, workingweight,
                                           workingresponse, l, like);
    if (*weight == 0.0)
      {
      *workingweight = 0.0;
      *workingresponse = 1.0;
      l = 0.0;
      return l;
      }
    *workingweight *= (*weight);

    return (*weight)*l;
    }
  }


void DISTR_gamlss::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    set_worklin();

//  modify_worklin();

  }


void DISTR_gamlss::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_gamlss::update_end(void)
  {

  // helpmat1 stores exp(linpred)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    if (*worklin <= linpredminlimit)
      *pmu  = exp(linpredminlimit);
//    else if (*worklin >= linpredmaxlimit)
//      *pmu  = exp(linpredmaxlimit);
    else
      *pmu = exp(*worklin);

    }

  }

//------------------------------------------------------------------------------
//------------------------- CLASS DISTR_negbinzip_mu ---------------------------
//------------------------------------------------------------------------------

void DISTR_negbinzip_mu::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {
        if (*workresp != int(*workresp))
          {
          errors=true;
          errormessages.push_back("ERROR: response must be integer values\n");
          }

        if (*workresp < 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


DISTR_negbinzip_mu::DISTR_negbinzip_mu(GENERAL_OPTIONS * o, const datamatrix & r,
                                 const datamatrix & w)
  : DISTR(o,r,w)

  {

  predictor_name = "mu";
  outpredictor = true;
  outexpectation = true;
  predict_mult = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;

  helpmat1 = datamatrix(nrobs,1,1);

  family = "Zero_Inflated_Negative_Binomial - mu";

  updateIWLS = true;

  linpredminlimit=-10;
  linpredmaxlimit=15;
  }


const DISTR_negbinzip_mu & DISTR_negbinzip_mu::operator=(
const DISTR_negbinzip_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  worklinpi = nd.worklinpi;
  workexplinpi = nd.workexplinpi;
  workonempi = nd.workonempi;
  worklindelta = nd.worklindelta;
  workexplindelta = nd.workexplindelta;
  distrpi = nd.distrpi;
  distrdelta = nd.distrdelta;
  counter = nd.counter;
  return *this;
  }


DISTR_negbinzip_mu::DISTR_negbinzip_mu(const DISTR_negbinzip_mu & nd)
   : DISTR(DISTR(nd))
  {
  worklinpi = nd.worklinpi;
  workonempi = nd.workonempi;
  workexplinpi = nd.workexplinpi;
  worklindelta = nd.worklindelta;
  workexplindelta = nd.workexplindelta;
  distrpi = nd.distrpi;
  distrdelta = nd.distrdelta;
  counter = nd.counter;
  }


void DISTR_negbinzip_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_negbinzip_mu::set_worklinpidelta(void)
  {

  if (distrpi->linpred_current==1)
    worklinpi = distrpi->linearpred1.getV();
  else
    worklinpi = distrpi->linearpred2.getV();

  workonempi = distrpi->helpmat1.getV();
  workexplinpi = distrpi->helpmat2.getV();

  if (distrdelta->linpred_current==1)
    worklindelta = distrdelta->linearpred1.getV();
  else
    worklindelta = distrdelta->linearpred2.getV();

  workexplindelta = distrdelta->helpmat1.getV();

  }


void DISTR_negbinzip_mu::modify_worklinpidelta(void)
  {

  if (counter<nrobs-1)
    {
    counter++;
    worklinpi++;
    workonempi++;
    workexplinpi++;
    worklindelta++;
    workexplindelta++;
    }
  else
    {
    counter=0;
    }

  }


double DISTR_negbinzip_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_negbinzip_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[2]));
  }

double DISTR_negbinzip_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_negbinzip_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
//    double p =  (*param[2])/((*param[0])+(*param[2]));
//    double r = (*param[0]);
//    double kplusone = 1 + (*response[2]);

    return 0;
//    return ( (*param[1])+(1-(*param[1]))*(1-randnumbers::incomplete_beta(kplusone,r,p)) );
    }

double DISTR_negbinzip_mu::loglikelihood(double * response, double * linpred,
                                         double * weight)
  {
  if (*weight == 0)
    {
    if (counter==0)
      set_worklinpidelta();

    modify_worklinpidelta();

    return 0;
    }
  else
    return (*weight)*loglikelihood_weightsone(response,linpred);
  }


double DISTR_negbinzip_mu::loglikelihood_weightsone(double * response,
                                                    double * linpred)
  {

  if (counter==0)
    set_worklinpidelta();

  double mu;
  if (*linpred <= linpredminlimit)
    mu = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    mu = exp(linpredmaxlimit);
  else
    mu = exp(*linpred);

  double deltaplusmu = (*workexplindelta)+mu;
  double pot = pow((*workexplindelta)/deltaplusmu,(*workexplindelta));

  double l;

  if (*response==0)
    {
    l = log((*workexplinpi)+pot);
    }
  else // response > 0
    {
    double deltaplusy = (*response)+(*workexplindelta);
    l = (*response)*(*linpred) - deltaplusy * log((*workexplindelta)+mu);
    }

  modify_worklinpidelta();

  return l;

  }


void DISTR_negbinzip_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

  double lambda;
  if (*linpred[predstart_mumult+2] <= linpredminlimit)
    lambda = exp(linpredminlimit);
  else
    lambda = exp(*linpred[predstart_mumult+2]);


  double explinpi;
  double min = distrpi->linpredminlimit;
  if (*linpred[predstart_mumult+1] <= min)
    explinpi = exp(min);
  else
    explinpi = exp(*linpred[predstart_mumult+1]);

  *mu = 1/(1+explinpi)*lambda;

  }


void DISTR_negbinzip_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

  double mu;
  if (*linpred[2] <= linpredminlimit)
    mu = exp(linpredminlimit);
  else
    mu = exp(*linpred[2]);


  double explinpi;
  double min = distrpi->linpredminlimit;
  if (*linpred[1] <= min)
    explinpi = exp(min);
  else
    explinpi = exp(*linpred[1]);


  double explindelta;
  min = distrdelta->linpredminlimit;
  if (*linpred[0] <= min)
    explindelta = exp(min);
  else
    explindelta = exp(*linpred[0]);


  double l= -log(1+explinpi);

  if (*response[2]==0)
    {
    double pot = pow(explindelta/(explindelta+mu),explindelta);
    l += log(explinpi+pot);
    }
  else // response > 0
    {
    double help1 = *response[2]+explindelta;
    double help2 = *response[2]+1;
    l +=    randnumbers::lngamma_exact(help1)
          - randnumbers::lngamma_exact(help2)
          - randnumbers::lngamma_exact(explindelta)
          + explindelta*(*linpred[0])
          + (*response[2])*(*linpred[2])
          - (explindelta+(*response[2]))*log(explindelta+mu);
    }

  *deviance = -2*l;
  }


double DISTR_negbinzip_mu::compute_iwls(double * response, double * linpred,
                                     double * weight, double * workingweight,
                                     double * workingresponse, const bool & like)
  {
  if (*weight == 0)
    {
    if (counter==0)
      set_worklinpidelta();

    *workingweight = 0;

    modify_worklinpidelta();

    return 0;
    }
  else
    {

    double l=0;
    compute_iwls_wweightschange_weightsone(response,linpred, workingweight,
                                           workingresponse, l,
                                           like);
     *workingweight *= (*weight);

     return (*weight)*l;

    }

  }


void DISTR_negbinzip_mu::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    set_worklinpidelta();

  double mu;
  if (*linpred <= linpredminlimit)
    mu = exp(linpredminlimit);
  else
    mu = exp(*linpred);

  double pi = 1-(*workonempi);

  double deltaplusmu = (*workexplindelta)+mu;

  double pot = pow((*workexplindelta)/deltaplusmu,(*workexplindelta));
  double denom = pi+(*workonempi)*pot;
  double denomfull = denom*deltaplusmu;
  double denomfull2 = denomfull*deltaplusmu;

  double nu = ((*response)*(*workexplindelta) - (*workexplindelta)*mu)/
               deltaplusmu;
  if (*response == 0)
    nu += (pi*(*workexplindelta)*mu ) /denomfull;

  *workingweight = (*workexplindelta)*mu*(*workonempi)/deltaplusmu -
                    (pi* (*workonempi) * pow((*workexplindelta),2) * pow(mu,2) * pot)  /denomfull2;

  *workingresponse = *linpred + nu/(*workingweight);


  if (compute_like)
    {

    if (*response==0)
      {
      like += log((*workexplinpi)+pot);
      }
    else // response > 0
      {
      double deltaplusy = (*response)+(*workexplindelta);
      like += (*response)*(*linpred) - deltaplusy * log((*workexplindelta)+mu);
      }

    }

  modify_worklinpidelta();

  }


void DISTR_negbinzip_mu::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_negbinzip_mu::update_end(void)
  {

  // helpmat1 stores mu, i.e. exp(linpred)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {

    if (*worklin <= linpredminlimit)
      *pmu = exp(linpredminlimit);
//    else if (*worklin >= linpredmaxlimit)
//      *pmu = exp(linpredmaxlimit);
    else
      *pmu = exp(*worklin);

    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS DISTR_negbinzip_pi ---------------------------
//------------------------------------------------------------------------------


DISTR_negbinzip_pi::DISTR_negbinzip_pi(GENERAL_OPTIONS * o, const datamatrix & r,
                                 const datamatrix & w)
  : DISTR(o,r,w)

  {

  predictor_name = "pi";
  outpredictor = true;
  outexpectation = false;

  predict_mult = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;

  // helpmat1 stores 1-pi
  // helpmat2 stores exp(eta_pi)
  helpmat1=datamatrix(nrobs,1,0.5);
  helpmat2=datamatrix(nrobs,1,1.0);

  family = "Zero_Inflated_Negative_Binomial - pi";
  updateIWLS = true;

  linpredminlimit=-10;
  linpredmaxlimit= 10;

  }


const DISTR_negbinzip_pi & DISTR_negbinzip_pi::operator=(
const DISTR_negbinzip_pi & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  worklinmu = nd.worklinmu;
  workexplinmu = nd.workexplinmu;
  worklindelta = nd.worklindelta;
  workexplindelta = nd.workexplindelta;
  distrmu = nd.distrmu;
  distrdelta = nd.distrdelta;
  counter = nd.counter;
  return *this;
  }


DISTR_negbinzip_pi::DISTR_negbinzip_pi(const DISTR_negbinzip_pi & nd)
   : DISTR(DISTR(nd))
  {
  worklinmu = nd.worklinmu;
  workexplinmu = nd.workexplinmu;
  worklindelta = nd.worklindelta;
  workexplindelta = nd.workexplindelta;
  distrmu = nd.distrmu;
  distrdelta = nd.distrdelta;
  counter = nd.counter;
  }


void DISTR_negbinzip_pi::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (pi): logistic distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_negbinzip_pi::set_worklinmudelta(void)
  {
  if (distrmu->linpred_current==1)
    worklinmu = distrmu->linearpred1.getV();
  else
    worklinmu = distrmu->linearpred2.getV();

  workexplinmu = distrmu->helpmat1.getV();


  if (distrdelta->linpred_current==1)
    worklindelta = distrdelta->linearpred1.getV();
  else
    worklindelta = distrdelta->linearpred2.getV();

  workexplindelta = distrdelta->helpmat1.getV();
  }


void DISTR_negbinzip_pi::modify_worklinmudelta(void)
  {
  if (counter<nrobs-1)
    {
    counter++;
    worklinmu++;
    workexplinmu++;
    worklindelta++;
    workexplindelta++;
    }
  else
    {
    counter=0;
    }
  }


double DISTR_negbinzip_pi::get_intercept_start(void)
  {
  /*
  unsigned i;
  double * responsep = response.getV();
  double m = 0;
  for (i=0;i<nrobs;i++,responsep++)
    {
    if (*responsep==0)
      m += 1;
    }

  m /= nrobs;

  return log(m/(1-m));
  */
  return 0;
  }

  void DISTR_negbinzip_pi::compute_param_mult(vector<double *>  linpred,double * param)
  {
  double el = exp(*linpred[1]);
  *param = el/(1+el);
  }


double DISTR_negbinzip_pi::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {
  if (*weight == 0)
    {
    if (counter==0)
      set_worklinmudelta();

    modify_worklinmudelta();

    return 0;
    }
  else
    return (*weight)*loglikelihood_weightsone(response,linpred);
  }


double DISTR_negbinzip_pi::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {

  if (counter==0)
    set_worklinmudelta();

  double explinpredpi;
  if (*linpred <= linpredminlimit)
    explinpredpi = exp(linpredminlimit);
  else
    explinpredpi = exp(*linpred);

  double l = -log(1+explinpredpi);

  if (*response==0)
    {
    double deltaplusmu = (*workexplindelta)+(*workexplinmu);
    double pot = pow((*workexplindelta)/deltaplusmu,(*workexplindelta));

    l += log(explinpredpi+pot);
    }

  modify_worklinmudelta();

  return l;

  }

double DISTR_negbinzip_pi::compute_iwls(double * response, double * linpred,
                                     double * weight, double * workingweight,
                                     double * workingresponse, const bool & like)
  {
  if (*weight == 0)
    {
    if (counter==0)
      set_worklinmudelta();

    *workingweight = 0;

    modify_worklinmudelta();

    return 0;
    }
  else
    {

    double l=0;
    compute_iwls_wweightschange_weightsone(response,linpred, workingweight,
                                           workingresponse, l,
                                           like);
     *workingweight *= (*weight);

     return (*weight)*l;

    }

  }


void DISTR_negbinzip_pi::compute_iwls_wweightschange_weightsone(
                                         double * response,
                                         double * linpred,
                                         double * workingweight,
                                         double * workingresponse,
                                         double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    set_worklinmudelta();

  double explinpredpi;
  if (*linpred <= linpredminlimit)
    explinpredpi = exp(linpredminlimit);
  else
    explinpredpi = exp(*linpred);

  double oneminuspi = 0.001+0.998/(1+explinpredpi);
  double pi = 1-oneminuspi;

  double deltaplusmu = (*workexplindelta)+(*workexplinmu);

  double pot = pow((*workexplindelta)/deltaplusmu,(*workexplindelta));
  double denom = pi+oneminuspi*pot;

  double nu = -pi;
  if (*response == 0)
    nu +=  pi/denom;

  *workingweight = (pow(pi,2)*oneminuspi*(1-pot)) /denom;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    like -= log(1+explinpredpi);

    if (*response==0)
      like += log(explinpredpi+pot);
    }

  modify_worklinmudelta();

  }


void DISTR_negbinzip_pi::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_negbinzip_pi::update_end(void)
  {

  // helpmat1 stores 1-pi
  // helpmat2 stores exp(eta_pi)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * ppi = helpmat1.getV();
  double * workexplin = helpmat2.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,ppi++,worklin++,workexplin++)
    {

    if (*worklin <= linpredminlimit)
      *workexplin = exp(linpredminlimit);
//    else if (*worklin >= linpredmaxlimit)
//      *workexplin = exp(linpredmaxlimit);
    else
      *workexplin = exp(*worklin);

    *ppi = 0.001+0.998/(1+(*workexplin));
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS DISTR_negbinzip_delta ------------------------
//------------------------------------------------------------------------------

DISTR_negbinzip_delta::DISTR_negbinzip_delta(GENERAL_OPTIONS * o,
                                             const datamatrix & r,
                                             double & ss, int & strmax,
                                             int & sts, bool & sl,
                                             const datamatrix & w)
  : DISTR(o,r,w)

  {

  predictor_name = "delta";
  outpredictor = true;
  outexpectation = false;

  responsemax = response.max(0);

  predict_mult = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;
  helpmat1=datamatrix(nrobs,1,1);

  family = "Zero_Inflated_Negative_Binomial - delta";
  updateIWLS = true;

  stoprmax = strmax;
  if (stoprmax < responsemax)
    stoprmax = responsemax;

  stopsum = ss;
  nrbetween = sts;

  slow = sl;

  E_dig_y_delta_m = datamatrix(nrobs,1,0);
  E_trig_y_delta_m = datamatrix(nrobs,1,0);

  linpredminlimit=-10;
  linpredmaxlimit=10;
  }


const DISTR_negbinzip_delta & DISTR_negbinzip_delta::operator=(
const DISTR_negbinzip_delta & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));

  responsemax = nd.responsemax;

  worklinmu = nd.worklinmu;
  workexplinmu = nd.workexplinmu;
  worklinpi = nd.worklinpi;
  workonempi = nd.workonempi;
  workexplinpi = nd.workexplinpi;;

  distrmu = nd.distrmu;
  distrpi = nd.distrpi;
  counter = nd.counter;

  stopsum = nd.stopsum;
  stoprmax = nd.stoprmax;
  fraclimit = nd.fraclimit;
  slow = nd.slow;
  nrbetween = nd.nrbetween;

  E_dig_y_delta_m = nd.E_dig_y_delta_m;
  E_trig_y_delta_m = nd.E_trig_y_delta_m;
  Ep = nd.Ep;
  Ep_trig = nd.Ep_trig;

  return *this;
  }


DISTR_negbinzip_delta::DISTR_negbinzip_delta(const DISTR_negbinzip_delta & nd)
   : DISTR(DISTR(nd))
  {

  responsemax = nd.responsemax;

  worklinmu = nd.worklinmu;
  workexplinmu = nd.workexplinmu;
  worklinpi = nd.worklinpi;
  workonempi = nd.workonempi;
  workexplinpi = nd.workexplinpi;;

  distrmu = nd.distrmu;
  distrpi = nd.distrpi;
  counter = nd.counter;

  stopsum = nd.stopsum;
  stoprmax = nd.stoprmax;
  fraclimit = nd.fraclimit;
  slow = nd.slow;
  nrbetween = nd.nrbetween;


  E_dig_y_delta_m = nd.E_dig_y_delta_m;
  E_trig_y_delta_m = nd.E_trig_y_delta_m;
  Ep = nd.Ep;
  Ep_trig = nd.Ep_trig;

  }


void DISTR_negbinzip_delta::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (delta): exponential\n");
  if (slow)
    optionsp->out("  Estimation type: slow\n");
  else
    optionsp->out("  Estimation type: fast\n");
  optionsp->out("  Stop criteria for approximating expected values\n");
  optionsp->out("  in working weights of delta equation:\n");
  optionsp->out("    Maximum values          : "  + ST::inttostring(stoprmax) +  "\n");
  optionsp->out("    Cumulative probabilities: "  + ST::doubletostring(stopsum) +  "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_negbinzip_delta::set_worklinmupi(void)
  {

  if (distrmu->linpred_current==1)
    worklinmu = distrmu->linearpred1.getV();
  else
    worklinmu = distrmu->linearpred2.getV();

  workexplinmu = distrmu->helpmat1.getV();


  if (distrpi->linpred_current==1)
    worklinpi = distrpi->linearpred1.getV();
  else
    worklinpi = distrpi->linearpred2.getV();

  workonempi = distrpi->helpmat1.getV();
  workexplinpi = distrpi->helpmat2.getV();

  }


void DISTR_negbinzip_delta::modify_worklinmupi(void)
  {

  if (counter<nrobs-1)
    {
    counter++;
    worklinmu++;
    workexplinmu++;
    worklinpi++;
    workonempi++;
    workexplinpi++;
    }
  else
    {
    counter=0;
    }

  }

double DISTR_negbinzip_delta::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_negbinzip_delta::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[0]));
  }

double DISTR_negbinzip_delta::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {

  if (*weight == 0)
    {
    if (counter==0)
      set_worklinmupi();

    modify_worklinmupi();

    return 0;
    }
  else
    return (*weight)*loglikelihood_weightsone(response,linpred);

  }


double DISTR_negbinzip_delta::loglikelihood_weightsone(double * response,
                                                       double * linpred)
  {

  if (counter==0)
    set_worklinmupi();

  double delta;
  if (*linpred <= linpredminlimit)
    delta = exp(linpredminlimit);
  else
    delta = exp(*linpred);

  double y_plus_delta = (*response)+delta;
  double delta_plus_mu = delta + (*workexplinmu);
  double delta_div_delta_plus_mu = delta/delta_plus_mu;
  double pot = pow(delta_div_delta_plus_mu,delta);
  double delta_linpred = delta*(*linpred);
  double log_delta_plus_mu = log(delta_plus_mu);


  double l;
  if (*response==0)
    l = log((*workexplinpi)+ pot);
  else
    l = randnumbers::lngamma_exact(y_plus_delta) -
        randnumbers::lngamma_exact(delta) + delta_linpred-
             (y_plus_delta)*log_delta_plus_mu;


  modify_worklinmupi();

  return l;

  }


double DISTR_negbinzip_delta::compute_iwls(double * response, double * linpred,
                                     double * weight, double * workingweight,
                                     double * workingresponse, const bool & like)
  {
  if (*weight == 0)
    {
    if (counter==0)
      set_worklinmupi();

    *workingweight = 0;

    modify_worklinmupi();

    return 0;
    }
  else
    {

    double l=0;
    compute_iwls_wweightschange_weightsone(response,linpred, workingweight,
                                           workingresponse, l,
                                           like);
     *workingweight *= (*weight);

     return (*weight)*l;

    }

  }



void DISTR_negbinzip_delta::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

    if (counter==0)
      {
      set_worklinmupi();
      Ep = E_dig_y_delta_m.getV();
      Ep_trig = E_trig_y_delta_m.getV();
      }

  double delta;
  if (*linpred <= linpredminlimit)
    delta = exp(linpredminlimit);
  else
    delta = exp(*linpred);


    double delta2 = delta*delta;

    double pi = 1-(*workonempi);
    double log_one_plus_explinpi = log(1+(*workexplinpi));

    double y_plus_delta = (*response)+delta;
    double digamma_delta = randnumbers::digamma_exact(delta);
    double trigamma_delta = randnumbers::trigamma_exact(delta);
    double digamma_y_plus_delta = randnumbers::digamma_exact(y_plus_delta);
    double delta_plus_mu = delta + (*workexplinmu);
    double delta_div_delta_plus_mu = delta/delta_plus_mu;
    double log_delta_div_delta_plus_mu = log(delta_div_delta_plus_mu);
    double mu_minus_y = (*workexplinmu) - (*response);
    double mu_div_delta_plus_mu = (*workexplinmu)/delta_plus_mu;
    double pot = pow(delta_div_delta_plus_mu,delta);
    double denom = pi+(*workonempi)*pot;
    double log_plus_term = log_delta_div_delta_plus_mu + mu_div_delta_plus_mu;
    double log_plus_term2 = pow(log_plus_term,2);
    double log_explinpi_plus_pot = log((*workexplinpi)+ pot);

    double nu = delta*(digamma_y_plus_delta -
                       digamma_delta +
                       log_delta_div_delta_plus_mu+
                       mu_minus_y/delta_plus_mu);

    if (*response == 0)
      nu -= (delta*pi*log_plus_term)/
            denom;
    // -------------------------------------------------------------------------

    double E_digamma_y_delta;
    double E_trigamma_y_delta;
    double lngamma_delta = randnumbers::lngamma_exact(delta);
    double log_delta_plus_mu = log(delta_plus_mu);
    double delta_linpred = delta*(*linpred);

    if ((optionsp->nriter < 1) ||
        slow ||
        (optionsp->nriter % nrbetween == 0)
        )
      {

      double sum =0;
      double L = exp(-log_one_plus_explinpi+log_explinpi_plus_pot);
      sum+= L;
      E_digamma_y_delta  = digamma_delta*L;
      E_trigamma_y_delta = trigamma_delta*L;
      int k = 1;
      double k_plus_delta;
      double k_plus_one;

      while ((sum < stopsum) && (k<=stoprmax))
        {
        k_plus_delta = k+delta;
        k_plus_one = k+1;

        L = exp(-log_one_plus_explinpi +
                 randnumbers::lngamma_exact(k_plus_delta) -
                 randnumbers::lngamma_exact(k_plus_one) -
                 lngamma_delta + delta_linpred+k*(*worklinmu)-
                 (delta+k)*log_delta_plus_mu);

        sum+= L;
        E_digamma_y_delta  +=  randnumbers::digamma_exact(k_plus_delta)*L;
        E_trigamma_y_delta +=  randnumbers::trigamma_exact(k_plus_delta)*L;
        k++;
        }

      E_digamma_y_delta  -= digamma_delta;
      E_trigamma_y_delta -= trigamma_delta;

      E_digamma_y_delta  *= delta;
      E_trigamma_y_delta *= delta2;


      *Ep = E_digamma_y_delta;
      *Ep_trig = E_trigamma_y_delta;
      }
    else
      {
      E_digamma_y_delta  = (*Ep);
      E_trigamma_y_delta = (*Ep_trig);
      }

    //--------------------------------------------------------------------------

  *workingweight = -delta*(*workonempi)*
                   (log_delta_div_delta_plus_mu+mu_div_delta_plus_mu) -
                    ((*workonempi)*pi*delta2*pot*log_plus_term2)/denom -
                    E_digamma_y_delta - E_trigamma_y_delta;

  if (*workingweight <=0)
    *workingweight = 0.001;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    if (*response==0)
      {

      like += log_explinpi_plus_pot;

      }
    else // response > 0
      {

      like += randnumbers::lngamma_exact(y_plus_delta) -
              lngamma_delta + delta_linpred-
             (y_plus_delta)*log_delta_plus_mu;

      }

    }

  modify_worklinmupi();
  Ep++;
  Ep_trig++;

  }




void DISTR_negbinzip_delta::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_negbinzip_delta::update_end(void)
  {

  // helpmat1 stores exp(linpred)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * l = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,worklin++,l++)
    {

    if (*worklin <= linpredminlimit)
      *l = exp(linpredminlimit);
    else
      *l = exp(*worklin);

    }

  }


//------------------------------------------------------------------------------
//---------------------- CLASS DISTRIBUTION_ziplambda --------------------------
//------------------------------------------------------------------------------

void DISTR_ziplambda::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {
        if (*workresp != int(*workresp))
          {
          errors=true;
          errormessages.push_back("ERROR: response must be integer values\n");
          }

        if (*workresp < 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


DISTR_ziplambda::DISTR_ziplambda(GENERAL_OPTIONS * o, const datamatrix & r,
                                 const datamatrix & w)
  : DISTR(o,r,w)

  {

  predictor_name = "lambda";
  outpredictor = true;
  outexpectation = true;
  predict_mult = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;

  family = "ZIP";
  updateIWLS = true;

  linpredminlimit=-10;
  linpredmaxlimit=15;

  check_errors();
  }


const DISTR_ziplambda & DISTR_ziplambda::operator=(const DISTR_ziplambda & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  distrpi = nd.distrpi;
  counter = nd.counter;
  worklinpi = nd.worklinpi;
  workonempi = nd.workonempi;
  workexplinpi = nd.workexplinpi;
  return *this;
  }


DISTR_ziplambda::DISTR_ziplambda(const DISTR_ziplambda & nd)
   : DISTR(DISTR(nd))
  {
  distrpi = nd.distrpi;
  counter = nd.counter;
  worklinpi = nd.worklinpi;
  workonempi = nd.workonempi;
  workexplinpi = nd.workexplinpi;
  }


void DISTR_ziplambda::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (lambda): exponential\n");
  optionsp->out("  Response function (pi): logistic distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_ziplambda::set_worklinpi(void)
  {
  if (distrpi->linpred_current==1)
    worklinpi = distrpi->linearpred1.getV();
  else
    worklinpi = distrpi->linearpred2.getV();

  workexplinpi = distrpi->helpmat1.getV();
  workonempi = distrpi->helpmat2.getV();

  }


void DISTR_ziplambda::modify_worklinpi(void)
  {
  if (counter<nrobs-1)
    {
    counter++;
    worklinpi++;
    workexplinpi++;
    workonempi++;
    }
  else
    {
    counter=0;
    }
  }


double DISTR_ziplambda::get_intercept_start(void)
  {
  return 0; //log(response.mean(0));
  }

double DISTR_ziplambda::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    if((*response[1]) == 0)
    {
    //    double p = (*param[0])+(1-(*param[0]))*gsl_ran_poisson_pdf((*response[1]), (*param[1])) ;
    //    return p;
        return 0;
    }
    else
    {
        //double p = (1-(*param[0]))*gsl_ran_poisson_pdf((*response[1]), (*param[1]));
        //return p;
        return 0;
    }

    }

double DISTR_ziplambda::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
  /*  unsigned int y = int (*response[1]);
    double cl = (*param[0]) + (1-(*param[0]))*gsl_cdf_poisson_P(((*response[1])-1), *param[1]);
    if ((*response[1]) == 0)
    {
        cl = 0;
    }
    double cr = (*param[0]) + (1-(*param[0]))*gsl_cdf_poisson_P((*response[1]), *param[1]);
    double u = randnumbers::uniform_ab(cl, cr);
    return u;*/
    return 0;
    }

void DISTR_ziplambda::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = exp((*linpred[1]));
  *param = arg;
  }

double DISTR_ziplambda::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {

  if (*weight == 0)
    {
    if (counter==0)
      set_worklinpi();

    modify_worklinpi();

    return 0;
    }
  else
    return (*weight)*loglikelihood_weightsone(response,linpred);

  }


double DISTR_ziplambda::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {

  if (counter==0)
    {
    set_worklinpi();
    }

  double lambda;
  double expminuslambda;

  if (*linpred <= linpredminlimit)
    lambda = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    lambda = exp(linpredmaxlimit);
  else
    lambda = exp(*linpred);

  expminuslambda = exp(-lambda);

  double l;

  if (*response==0)
    l = -log(1+(*workexplinpi)) + log((*workexplinpi)+expminuslambda);
  else // response > 0
    l = -log(1+(*workexplinpi)) + (*response)*(*linpred)-lambda;

  modify_worklinpi();

  return l;

  }


void DISTR_ziplambda::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

  double lambda;
  if (*linpred[predstart_mumult+1] <= linpredminlimit)
    lambda = exp(linpredminlimit);
  else
    lambda = exp(*linpred[predstart_mumult+1]);


  double explinpi;
  if (*linpred[predstart_mumult] <= distrpi->linpredminlimit)
    explinpi = exp(distrpi->linpredminlimit);
  else
    explinpi = exp(*linpred[predstart_mumult]);

  *mu = 1/(1+explinpi)*lambda;

  }


void DISTR_ziplambda::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   if (*weight[1] == 0)
     *deviance=0;
   else
     {

     double l;

     double lambda;
     double expminuslambda;
     double explinpi;

     if (*linpred[1] <= linpredminlimit)
       lambda = exp(linpredminlimit);
//     else if (*linpred[1] >= linpredmaxlimit)
//       lambda = exp(linpredmaxlimit);
     else
       lambda = exp(*linpred[1]);

     expminuslambda = exp(-lambda);

     if (*linpred[0] <= distrpi->linpredminlimit)
       explinpi = exp(distrpi->linpredminlimit);
//     else if (*linpred[0] >= distrpi->linpredmaxlimit)
//       explinpi = exp(distrpi->linpredmaxlimit);
     else
       explinpi = exp(*linpred[0]);

      if (*response[1]==0)
         {
         l= -log(1+ explinpi) + log(explinpi+ expminuslambda);
         }
       else // response > 0
         {
         double help1 = *response[1]+1;
         l= -log(1+ explinpi) + (*response[1])*(*linpred[1])- lambda
            - randnumbers::lngamma_exact(help1);
         }

     *deviance = -2*l;
     }

  }


double DISTR_ziplambda::compute_iwls(double * response, double * linpred,
                                     double * weight, double * workingweight,
                                     double * workingresponse, const bool & like)
  {
  if (*weight == 0)
    {
    if (counter==0)
      set_worklinpi();

    *workingweight = 0;

    modify_worklinpi();

    return 0;
    }
  else
    {

    double l=0;
    compute_iwls_wweightschange_weightsone(response,linpred, workingweight,
                                           workingresponse, l,
                                           like);
     *workingweight *= (*weight);

     return (*weight)*l;

    }

  }


void DISTR_ziplambda::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    {
    set_worklinpi();
    }

  double lambda;
  double expminuslambda;

  if (*linpred <= linpredminlimit)
    lambda = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    lambda = exp(linpredmaxlimit);
  else
    lambda = exp(*linpred);

  expminuslambda = exp(-lambda);

  double pi = 1-(*workonempi);
  double denom = pi+(*workonempi)*expminuslambda;

  double nu = (*response) - lambda;
  if (*response == 0)
    nu += pi*lambda/denom;

  *workingweight = (lambda* (*workonempi)*(denom-expminuslambda*lambda*pi))/denom;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    if (*response==0)
      {
      like += -log(1+ (*workexplinpi)) + log((*workexplinpi)+expminuslambda);
      }
    else // response > 0
      {
      like += -log(1+(*workexplinpi)) + (*response)*(*linpred)-lambda;
      }

    }

  modify_worklinpi();

  }


void DISTR_ziplambda::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_ziplambda::update_end(void)
  {

  // helpmat1 stores exp(-lambda)
  // helpmat2 stores lambda

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  if (helpmat1.rows() == 1)
    {
    helpmat1 = datamatrix(nrobs,1,0);
    helpmat2 = datamatrix(nrobs,1,0);
    }

  double * ph = helpmat1.getV();
  double * l = helpmat2.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,ph++,worklin++,l++)
    {
    if (*worklin <= linpredminlimit)
      *l = exp(linpredminlimit);
//    else if (*worklin >= linpredmaxlimit)
//      *l = exp(linpredmaxlimit);
    else
      *l = exp(*worklin);

    *ph = exp(-(*l));
    }

  }


//------------------------------------------------------------------------------
//------------------------ CLASS DISTRIBUTION_zippi ----------------------------
//------------------------------------------------------------------------------


DISTR_zippi::DISTR_zippi(GENERAL_OPTIONS * o, const datamatrix & r,
                                 const datamatrix & w)
  : DISTR(o,r,w)

  {

  predictor_name = "pi";
  outpredictor = true;
  outexpectation = false;


  maindistribution=false;
  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  counter = 0;
  helpmat1=datamatrix(nrobs,1,1);
  helpmat2=datamatrix(nrobs,1,0.5);


  family = "ZIP_pi";
  updateIWLS = true;

  linpredminlimit=-10;
  linpredmaxlimit= 10;
  }


const DISTR_zippi & DISTR_zippi::operator=(const DISTR_zippi & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  distrlambda = nd.distrlambda;
  counter = nd.counter;
  worklinlambda = nd.worklinlambda;
  workexpmlambda = nd.workexpmlambda;
  return *this;
  }


DISTR_zippi::DISTR_zippi(const DISTR_zippi & nd)
   : DISTR(DISTR(nd))
  {
  distrlambda = nd.distrlambda;
  counter = nd.counter;
  worklinlambda = nd.worklinlambda;
  workexpmlambda = nd.workexpmlambda;
  }


void DISTR_zippi::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (pi): logistic distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_zippi::get_intercept_start(void)
  {
  /*
  unsigned i;
  double * responsep = response.getV();
  double m = 0;
  for (i=0;i<nrobs;i++,responsep++)
    {
    if (*responsep==0)
      m += 1;
    }

  m /= nrobs;
  */

  return 0; // log(m/(1-m));
  }

void DISTR_zippi::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = exp((*linpred[0]));
  *param = arg/(1+arg);
  }


void DISTR_zippi::set_worklinlambda(void)
  {
  if (distrlambda->linpred_current==1)
    worklinlambda = distrlambda->linearpred1.getV();
  else
    worklinlambda = distrlambda->linearpred2.getV();

  workexpmlambda = distrlambda->helpmat1.getV();
  worklambda = distrlambda->helpmat2.getV();
  }


void DISTR_zippi::modify_worklinlambda(void)
  {
  if (counter<nrobs-1)
    {
    counter++;
    worklinlambda++;
    worklambda++;
    workexpmlambda++;
    }
  else
    {
    counter=0;
    }
  }


double DISTR_zippi::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {
  if (*weight == 0)
    {
    if (counter==0)
      set_worklinlambda();

    modify_worklinlambda();

    return 0;
    }
  else
    return (*weight)*loglikelihood_weightsone(response,linpred);
  }


double DISTR_zippi::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {

  if (counter==0)
    set_worklinlambda();

  double exptildeeta;
  if (*linpred <= linpredminlimit)
    exptildeeta = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    exptildeeta = exp(linpredmaxlimit);
  else
    exptildeeta = exp(*linpred);

  double l;

  if (*response==0)
    {
    l = -log(1+exptildeeta) + log(exptildeeta+(*workexpmlambda));
    }
  else // response > 0
    {
    l = -log(1+exptildeeta) + (*response)*(*worklinlambda)-(*worklambda);
    }

  modify_worklinlambda();

  return l;

  }


double DISTR_zippi::compute_iwls(double * response, double * linpred,
                                     double * weight, double * workingweight,
                                     double * workingresponse, const bool & like)
  {
  if (*weight == 0)
    {
    if (counter==0)
      set_worklinlambda();

    *workingweight = 0;

    modify_worklinlambda();

    return 0;
    }
  else
    {

    double l=0;
    compute_iwls_wweightschange_weightsone(response,linpred, workingweight,
                                           workingresponse, l,
                                           like);
     *workingweight *= (*weight);

     return (*weight)*l;

    }

  }


void DISTR_zippi::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  if (counter==0)
    {
    set_worklinlambda();
    }

  double exptildeeta;
  if (*linpred <= linpredminlimit)
    exptildeeta = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    exptildeeta = exp(linpredmaxlimit);
  else
    exptildeeta = exp(*linpred);

  double oneminuspi = 0.001+0.998/(1+exptildeeta);
  double pi = 1-oneminuspi;
  double denom = pi+oneminuspi* (*workexpmlambda);

  double nu = - pi;
  if (*response == 0)
    nu += pi/denom;

  *workingweight = (pi*pi*(1-(*workexpmlambda))*oneminuspi)/denom;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    if (*response==0)
      like += -log(1+exptildeeta) + log(exptildeeta+ (*workexpmlambda));
    else // response > 0
      like += -log(1+exptildeeta) + (*response)*(*worklinlambda)-(*worklambda);
    }

  modify_worklinlambda();

  }


void DISTR_zippi::posteriormode_end(void)
  {
  update_end();
  }


void DISTR_zippi::update_end(void)
  {

  // helpmat1 stores exp(tildeeta)
  // helpmat2 stores 1-pi

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  if (helpmat1.rows() == 1)
    {
    helpmat1 = datamatrix(nrobs,1,0);
    helpmat2 = datamatrix(nrobs,1,0);
    }

  double * ete = helpmat1.getV();
  double * wpi = helpmat2.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,ete++,worklin++,wpi++)
    {

    if (*worklin <= linpredminlimit)
      *ete = exp(linpredminlimit);
 //   else if (*worklin >= linpredmaxlimit)
 //     *ete = exp(linpredmaxlimit);
    else
      *ete = exp(*worklin);

    *wpi = 0.001+0.998/(1+(*ete));
    }

  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_normal2_sigma -------------------------
//------------------------------------------------------------------------------


DISTR_cnormal_sigma::DISTR_cnormal_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Censored Normal Distribution - sigma";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
  linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_cnormal_sigma::DISTR_cnormal_sigma(const DISTR_cnormal_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_cnormal_sigma & DISTR_cnormal_sigma::operator=(
                            const DISTR_cnormal_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_cnormal_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_cnormal_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[0]));
  }

double DISTR_cnormal_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma = exp(*linpred);
  double sigma_2 = pow(sigma,2);

  double l;

  if (*response <= 0)
    {
    l = log(randnumbers::Phi2(-*worklin[0]/sigma));
    }
  else
    {
    l = -0.5*log(2*M_PI)-pow(*response-(*worklin[0]),2)/(2*sigma_2) - log(sigma);
    }


  modify_worklin();

  return l;

  }

void DISTR_cnormal_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu)/(1+exp(eta_mu));

  if (counter==0)
    {
    set_worklin();
    }

    double sigma = exp(*linpred);
    double sigma_2 = pow(sigma,2);
    double p = pow(*response - *worklin[0],2);

    double nu;
    double ddist, pdist, mills;


     if (*response <= 0)
       {
       ddist = 1/sqrt(2*M_PI) *exp(-0.5*pow(*worklin[0],2)/sigma_2);
       pdist = randnumbers::Phi2(-(*worklin[0])/sigma);
       mills = ddist/pdist;

       nu = (*worklin[0])*mills/sigma;

       *workingweight = nu - pow(*worklin[0],3) * mills/sigma_2 + pow(*worklin[0],2) * pow(mills,2)/sigma_2;

/*
    cout << "response: " << *response << endl;
    cout << "linpred: " << *worklin[0] << endl;
    cout << "p: " << p << endl;
    cout << "sigma: " << sigma << endl;
    cout << "sigma_2: " << sigma_2 << endl;
    cout << "ddist: " << ddist << endl;
    cout << "pdist: " << pdist << endl;
    cout << "mills: " << mills << endl;
    cout << "mu_sigma: " << mu_sigma << endl;
    cout << "nu: " << nu << endl;
    cout << "workingweight: " << *workingweight << endl;
*/
       }
     else
       {
       nu = p/sigma_2-1;
       *workingweight = 2/sigma_2*p;

/*
    cout << "response: " << *response << endl;
    cout << "linpred: " << *worklin[0] << endl;
    cout << "p: " << p << endl;
    cout << "sigma: " << sigma << endl;
    cout << "sigma_2: " << sigma_2 << endl;
    cout << "nu: " << nu << endl;
    cout << "workingweight: " << *workingweight << endl;
*/
       }

    *workingresponse = *linpred + nu/(*workingweight);

/*
    cout << "workingresponse: " << *workingresponse << endl << endl;

    if(isnan(*workingresponse))
      cout << "arg (workingresponse, sigma)";
    if(isnan(*workingweight))
      cout << "arg (workingweight, sigma)";
*/

    if (compute_like)
      {

      if (*response <= 0)
        {
        like += log(randnumbers::Phi2(-(*worklin[0])/sigma));
        }
      else
        {
        like += -0.5*log(2*M_PI)-p/(2*sigma_2)- log(sigma);
        }

      }

  modify_worklin();

  }


void DISTR_cnormal_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_cnormal_sigma::update_end(void)
  {

  // helpmat1 stores sigma^2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = pow(exp(*worklin),2);
    }
  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_cnormal_mu --------------------------
//------------------------------------------------------------------------------
void DISTR_cnormal_mu::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {
      if(*workresp < 0)
        {
        errors=true;
        errormessages.push_back("ERROR: negative responses encountered\n");
        }
      if (*workweight > 0)
        {
        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;
      }
    }
  }


DISTR_cnormal_mu::DISTR_cnormal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Censored Normal Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  check_errors();
 // linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_cnormal_mu::DISTR_cnormal_mu(const DISTR_cnormal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  }


const DISTR_cnormal_mu & DISTR_cnormal_mu::operator=(
                            const DISTR_cnormal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


// FEHLT NOCH - NICHT KORREKT
void DISTR_cnormal_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma2
   // *linpred[1] = eta_mu

   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double sigma_2 = pow(exp(*linpred[0]),2);
     double mu = (*linpred[1]);

     double l;

       l = -0.5*log(2*PI)-0.5*log(sigma_2)-pow((((*response[0]))-mu),2)/(2*sigma_2);


    *deviance = -2*l;
    }

  }


double DISTR_cnormal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_cnormal_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[1]);
  }

 double DISTR_cnormal_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }


// FEHLT NOCH - NICHT KORREKT
double DISTR_cnormal_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = ((*response[1])-(*param[1]))/(*param[0]) ;

    return (randnumbers::Phi2(arg));
    }


double DISTR_cnormal_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double sigma = sqrt(*worktransformlin[0]);

  double l;


  if (*response <= 0)
    {
    l = log(randnumbers::Phi2(-mu/sigma));
    }
  else
    {
    l = -0.5*log(2*M_PI)-pow(*response-mu,2)/(2*(*worktransformlin[0])) - log(sigma);
    }


  modify_worklin();

  return l;

  }


void DISTR_cnormal_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

  double sigma = sqrt(*worktransformlin[0]);
  double p = *response - (*linpred);
  double nu;
  double ddist, pdist, mills;

/*
  cout << *worklin[0] << endl;
  cout << *worktransformlin[0] << endl;
  cout << distrp[0]->helpmat1(0,0) << endl;
  cout << sigma << endl;
*/

   if (*response <= 0)
     {
     ddist = 1/(sqrt(2*M_PI))*exp(-0.5*pow(*linpred,2)/(*worktransformlin[0]));
     pdist = randnumbers::Phi2(-(*linpred)/sigma);
     mills = ddist/pdist;
     nu = -mills/sigma;

     *workingweight =  -(*linpred)/(pow(sigma,3.0)) * mills + pow(mills,2)/(*worktransformlin[0]);
     }
   else
     {
     nu = p/(*worktransformlin[0]);
     *workingweight = 1/(*worktransformlin[0]);
     }


    *workingresponse = *linpred + nu/(*workingweight);

/*
    if(isnan(*workingresponse))
      {
      cout << counter << endl;
      cout << "arg (workingresponse, mu)" << endl;
      cout << "sigma: " << sigma << endl;
      cout << "mu: " << *linpred << endl;
      cout << "pdist:" << pdist << endl;
      cout << "ddist:" << ddist << endl;
      cout << "mills:" << mills << endl;
      }
    if(isnan(*workingweight))
      {
      cout << counter << endl;
      cout << "arg (workingweight, mu)" << endl;
      cout << "sigma: " << sigma << endl;
      cout << "mu: " << *linpred << endl;
      cout << "pdist:" << pdist << endl;
      cout << "ddist:" << ddist << endl;
      cout << "mills:" << mills << endl;
      }
*/

    if (compute_like)
      {


      if (*response <= 0)
        {
        like += log(randnumbers::Phi2(-(*linpred)/sigma));
        }
      else
        {
        like += -0.5*log(2*M_PI)-pow(p,2)/(2*(*worktransformlin[0])) - log(sigma);
        }

      }


  modify_worklin();

  }


void DISTR_cnormal_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+1]));
  }


void DISTR_cnormal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }

void DISTR_cnormal_mu::update_end(void)
  {
  // helpmat1 stores (eta_mu)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin);
//    double t = 0;
    }

  }
} // end: namespace MCMC

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




#include "distr_categorical.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//------------------ CLASS: DISTRIBUTION_logit_fruehwirth ----------------------
//------------------------------------------------------------------------------

DISTR_logit_fruehwirth::DISTR_logit_fruehwirth(const int h, GENERAL_OPTIONS * o,
                                             const datamatrix r,
 																				     const datamatrix & w)
 	: DISTR_binomial(o, r, w), H(h)
	{


  family = "Binomial_l1";
  updateIWLS = false;


  SQ = datamatrix(6,5,0);
  SQ(0,1) = 1/1.2131;
  SQ(1,1) = 1/2.9955;
  SQ(2,1) = 1/7.5458;

  SQ(0,4) = 1/0.68159;
  SQ(1,4) = 1/1.2419;
  SQ(2,4) = 1/2.2388;
  SQ(3,4) = 1/4.0724;
  SQ(4,4) = 1/7.4371;
  SQ(5,4) = 1/13.772;


  weights_mixed = datamatrix(6,5,0);
  weights_mixed(0,1) = 0.2522;
  weights_mixed(1,1) = 0.58523;
  weights_mixed(2,1) = 0.16257;


  weights_mixed(0,4) = 0.018446;
  weights_mixed(1,4) = 0.17268;
  weights_mixed(2,4) = 0.37393;
  weights_mixed(3,4) = 0.31697;
  weights_mixed(4,4) = 0.1089;
  weights_mixed(5,4) = 0.0090745;

	}

DISTR_logit_fruehwirth::DISTR_logit_fruehwirth(const DISTR_logit_fruehwirth & nd)
	: DISTR_binomial(DISTR_binomial(nd)) , H(nd.H), SQ(nd.SQ), weights_mixed(nd.weights_mixed)
	{
	}

const DISTR_logit_fruehwirth & DISTR_logit_fruehwirth::operator=(const DISTR_logit_fruehwirth & nd)
	{
	if (this==&nd)
  	return *this;
  DISTR_binomial::operator=(DISTR_binomial(nd));
  H = nd.H;
  SQ = nd.SQ;
  weights_mixed = nd.weights_mixed;
  return *this;
	}

/*
double DISTR_logit_fruehwirth::compute_MSE()
  {

  }

void DISTR_logit_fruehwirth::compute_mu()
  {
  // datamatrix test(3,3,0);
  // ofstream out("c:\\bayesx\\temp\\r.res");
  // test.prettyPrint(out);
  }
*/




void DISTR_logit_fruehwirth::outoptions()
  {
  DISTR::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("  Number of mixture components: " + ST::inttostring(H) + "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_logit_fruehwirth::update(void)
  {

  double * workresp;
  double * workwresp;
  double * weightwork;
  double * wweightwork;

  workresp = response.getV();
  workwresp = workingresponse.getV();
  weightwork = weight.getV();
  wweightwork = workingweight.getV();

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double lambda;
  double U;
  datamatrix weights_aux(H,1);


  for(int i=0;i<nrobs;i++,worklin++,workresp++,weightwork++,workwresp++,wweightwork++)
    {
    lambda= exp(*worklin);
    U = uniform();
    *workwresp = log(lambda*U+*workresp)-log(1-U+lambda*(1-*workresp));

    //weights_mixed
    for(int j=0; j < H; j++)
    	{
      weights_aux(j,0) = weights_mixed(j,H-2)*sqrt(SQ(j,H-2)) * exp(-1/2*  pow((*workresp - *worklin), 2)*SQ(j,H-2) );
      }

    //distribution function
    for(int j=1; j <H; j++)
    	{
       weights_aux(j,0) = weights_aux(j-1,0) + weights_aux(j,0);
      }

    U = uniform();
    U = U*weights_aux(H-1,0);	//scale to [0, max]


    int iaux = 0;
    while (U > weights_aux(iaux,0))
      {
      iaux++;
      }

    *wweightwork =  SQ(iaux,H-2);

    }
  }


bool DISTR_logit_fruehwirth::posteriormode(void)
  {
  return DISTR_binomial::posteriormode();
  }

//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_binomial --------------------------
//------------------------------------------------------------------------------


DISTR_binomial::DISTR_binomial(GENERAL_OPTIONS * o, const datamatrix & r,
                               const datamatrix & w)
  : DISTR(o,r,w)

  {

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  family = "Binomial";
  updateIWLS = true;
  }


const DISTR_binomial & DISTR_binomial::operator=(
                                      const DISTR_binomial & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  return *this;
  }


DISTR_binomial::DISTR_binomial(const DISTR_binomial & nd)
   : DISTR(DISTR(nd))
  {
  }


void DISTR_binomial::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_binomial::loglikelihood(double * response, double * linpred,
                                     double * weight) const
  {

  if (*linpred >= 10)
    return *weight *(*response * *linpred - *linpred);
  else
    return *weight *(*response * *linpred - log(1+exp(*linpred)));
  }


double DISTR_binomial::loglikelihood_weightsone(
                                  double * response, double * linpred) const
  {

  if (*linpred >= 10)
    return *response * (*linpred) - *linpred;
  else
    return *response * (*linpred) - log(1+exp(*linpred));
  }


void DISTR_binomial::compute_mu(const double * linpred,double * mu)
  {
  double el = exp(*linpred);
  *mu = el/(1+el);
  }


void DISTR_binomial::compute_deviance(const double * response,
                   const double * weight,const double * mu,double * deviance,
                   double * deviancesat, double * scale) const
  {

  if (*response==0)
    {
    *deviance = -2* *weight * log(1-*mu);
    *deviancesat = *deviance;
    }
  else if (*response == 1)
    {
    *deviance = -2* *weight*log(*mu);
    *deviancesat = *deviance;
    }
  else
    {
    *deviance = -2* *weight*( *response*log(*mu)+(1-*response)*log(1-*mu) );
    *deviancesat = *deviance +
    2* *weight*( *response*log(*response)+(1-*response)*log(1-*response) );
    }

  }


double DISTR_binomial::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  double el = exp(*linpred);
  double mu = el/(1+el);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;
  double v = mu*(1-mu);

  *workingweight = *weight * v;

  *workingresponse = *linpred + (*response - mu)/v;

  if (like)
    {
    if (*linpred >= 10)
      return *weight *(*response * *linpred - *linpred);
    else
      return *weight *(*response * *linpred - log(1+el));
    }
  else
    {
    return 0;
    }

  }


void DISTR_binomial::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {


  double el = exp(*linpred);
  double mu = el/(1+el);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;
  double v = mu*(1-mu);

  *workingweight =  v;

  *workingresponse = *linpred + (*response - mu)/v;

  if (compute_like)
    {
    if (*linpred >= 10)
      like += *response * (*linpred) - *linpred;
    else
      like += *response * (*linpred) - log(1+el);
    }

  }


void DISTR_binomial::compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  }

void DISTR_binomial::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  }


void DISTR_binomial::sample_responses(unsigned i,datamatrix & sr)
  {
  double * linpredp;

  if (linpred_current==1)
    linpredp = linearpred1.getV();
  else
    linpredp = linearpred2.getV();

  double * rp = sr.getV()+i;
  double mu;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    {
    compute_mu(linpredp,&mu);

    *rp = randnumbers::rand_binom(1,mu);

    }

  }


void DISTR_binomial::sample_responses_cv(unsigned i,datamatrix & linpred,
                                         datamatrix & sr)
  {

  double * linpredp;

  linpredp = linpred.getV();

  double * rp = sr.getV()+i;
  double mu;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    {
    compute_mu(linpredp,&mu);

    *rp = randnumbers::rand_binom(1,mu);
    }

  }



//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_binomialprobit --------------------
//------------------------------------------------------------------------------


DISTR_binomialprobit::DISTR_binomialprobit(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR(o,r,w)

  {

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  family = "Binomial";
  updateIWLS = false;
  }


const DISTR_binomialprobit & DISTR_binomialprobit::operator=(
                                      const DISTR_binomialprobit & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  return *this;
  }


DISTR_binomialprobit::DISTR_binomialprobit(const DISTR_binomialprobit & nd)
   : DISTR(DISTR(nd))
  {
  }


void DISTR_binomialprobit::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: standard normal (probit link)\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_binomialprobit::loglikelihood(double * response, double * linpred,
                                     double * weight) const
  {

  if (*weight!=0)
    {
    double mu = randnumbers::Phi2(*linpred);
    if (*response > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
    return 0;

  }



double DISTR_binomialprobit::loglikelihood_weightsone(
                                  double * response, double * linpred) const
  {

  double mu = randnumbers::Phi2(*linpred);
  if (*response > 0)
    return log(mu);
  else
    return log(1-mu);

  }


void DISTR_binomialprobit::compute_mu(const double * linpred,double * mu)
  {
  *mu = randnumbers::Phi2(*linpred);
  }


void DISTR_binomialprobit::compute_deviance(const double * response,
                   const double * weight,const double * mu,double * deviance,
                   double * deviancesat, double * scale) const
  {

  if (*weight !=  0)
    {

    if (*response<=0)
      {
      *deviance = -2*log(1-*mu);
      *deviancesat = *deviance;
      }
    else if (*response > 0)
      {
      *deviance = -2*log(*mu);
      *deviancesat = *deviance;
      }

    }
  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }

  }


double DISTR_binomialprobit::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  double  mu = randnumbers::Phi2(*linpred);

  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  double g = 1/pow(h,2);

  *workingweight = *weight / (mu*(1-mu) * g);


  *workingresponse = *linpred + (*response - mu)/h;

  if (like)
    {

    if (*response > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
    {
    return 0;
    }


  }



void DISTR_binomialprobit::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  double  mu = randnumbers::Phi2(*linpred);
  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  double g = 1/pow(h,2);

  *workingweight = 1.0 / (mu*(1-mu) * g);

  *workingresponse = *linpred + (*response - mu)/h;

  if (compute_like)
    {

    if (*response > 0)
      like+= log(mu);
    else
      like+= log(1-mu);
    }

  }


void DISTR_binomialprobit::compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  }

void DISTR_binomialprobit::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  }


void DISTR_binomialprobit::update(void)
  {

  double * worklin;
  double * workresp;
  double * workwresp;
  double * weightwork;
  double * workingweightwork;

  register unsigned i;


  if (optionsp->nriter==1)
    {

    weightwork = weight.getV();
    workingweightwork = workingweight.getV();

    for (i=0;i<nrobs;i++,weightwork++,workingweightwork++)
      {
      *workingweightwork = *weightwork;
      }

    }
  else
    {
    if (check_weightsone())
      wtype = wweightsnochange_one;
    else
      wtype = wweightsnochange_constant;
    }



  workresp = response.getV();
  workwresp = workingresponse.getV();
  weightwork = weight.getV();

  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();


  for(i=0;i<nrobs;i++,worklin++,workresp++,weightwork++,workwresp++)
    {

    if (*weightwork != 0)
      {
      if (*workresp > 0)
        *workwresp = trunc_normal2(0,20,*worklin,1);
      else
        *workwresp = trunc_normal2(-20,0,*worklin,1);
      }

    }

  DISTR::update();

  }



//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_poisson ---------------------------
//------------------------------------------------------------------------------


DISTR_poisson::DISTR_poisson(GENERAL_OPTIONS * o, const datamatrix & r,
                               const datamatrix & w)
  : DISTR(o,r,w)

  {

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  family = "Poisson";
  updateIWLS = true;
  }


const DISTR_poisson & DISTR_poisson::operator=(
                                      const DISTR_poisson & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  return *this;
  }


DISTR_poisson::DISTR_poisson(const DISTR_poisson & nd)
   : DISTR(DISTR(nd))
  {
  }


void DISTR_poisson::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_poisson::loglikelihood(double * response, double * linpred,
                                     double * weight) const
  {

  return *weight * (*response * *linpred - exp(*linpred));

  }


double DISTR_poisson::loglikelihood_weightsone(
                                  double * response, double * linpred) const
  {
  return  *response * (*linpred) - exp(*linpred);
  }


void DISTR_poisson::compute_mu(const double * linpred,double * mu)
  {
  *mu = exp(*linpred);
  }


void DISTR_poisson::compute_deviance(const double * response,
                   const double * weight,const double * mu,double * deviance,
                   double * deviancesat, double * scale) const
  {

  if (*response==0)
    {
    *deviance = 2* *weight * *mu;
    *deviancesat = *deviance;
    }
  else
    {
    *deviance = -2* *weight*(*response*log(*mu)-*mu);
    *deviancesat = *deviance+2 *
                     *weight*(*response * log(*response) - *response);
    }

  }


double DISTR_poisson::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  double mu = exp(*linpred);

  *workingweight = *weight * mu;

  *workingresponse = *linpred + (*response - mu)/mu;

   if (like)
     {
     return *weight * (*response * (*linpred) - mu);
     }
   else
     return 0;

  }


void DISTR_poisson::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  *workingweight = exp(*linpred);

  *workingresponse = *linpred + (*response - *workingweight)/(*workingweight);

   if (compute_like)
     {
     like+=  *response * (*linpred) - (*workingweight);
     }

  }


void DISTR_poisson::compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  }

void DISTR_poisson::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  }


void DISTR_poisson::sample_responses(unsigned i,datamatrix & sr)
  {
  double * linpredp;

  if (linpred_current==1)
    linpredp = linearpred1.getV();
  else
    linpredp = linearpred2.getV();

  double * rp = sr.getV()+i;
  double mu;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    {
    compute_mu(linpredp,&mu);

    *rp = randnumbers::rand_pois(mu);

    }

  }


void DISTR_poisson::sample_responses_cv(unsigned i,datamatrix & linpred,
                                         datamatrix & sr)
  {

  double * linpredp;

  linpredp = linpred.getV();

  double * rp = sr.getV()+i;
  double mu;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    {
    compute_mu(linpredp,&mu);

    *rp = randnumbers::rand_pois(mu);
    }

  }



} // end: namespace MCMC




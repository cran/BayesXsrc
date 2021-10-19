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

#include "distr_gamlss_nadja.h"
#include "Random.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_betainf1_tau ------------------------
//------------------------------------------------------------------------------


DISTR_betainf1_tau::DISTR_betainf1_tau(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,0,w)
  {
  family = "Beta One Inflated Distribution - tau";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "tau";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_betainf1_tau::DISTR_betainf1_tau(const DISTR_betainf1_tau & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf1_tau & DISTR_betainf1_tau::operator=(
                            const DISTR_betainf1_tau & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_betainf1_tau::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_tau
   // *linpred[1] = eta_nu
   // *linpred[2] = eta_sigma2
   // *linpred[3] = eta_mu

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double taup = exp(*linpred[2]);
     double sigma_2 = exp(*linpred[1])/(1+exp(*linpred[1]));
     double mu = exp(*linpred[0])/(1+exp(*linpred[0]));
     double help = (1-sigma_2)/sigma_2;
     double one_minus_mu_help = (1-mu)*help;
     double mu_help = mu*help;
     double one_taup = (1+taup);

     double l;

     if ((*response[2])==1)
     {
         l = log(taup) - log(one_taup);
     }
      else
       l = (mu_help-1)*log(*response[1]) +
			(one_minus_mu_help-1)*log(1-(*response[1]))-
			randnumbers::lngamma_exact(mu_help)-
			randnumbers::lngamma_exact(one_minus_mu_help)+
			randnumbers::lngamma_exact(help)- log(one_taup);


    *deviance = -2*l;
    }

  }


double DISTR_betainf1_tau::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_betainf1_tau::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[2]);
  }

 double DISTR_betainf1_tau::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_betainf1_tau::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
//    double a =  (*param[2])*(*param[3]);
//    double b = (*param[2])*(1-(*param[3]));
 //   double frac = 1 + (*param[0]) + (*param[1]);
    return 0;
//
//    return ( ((*param[0])+(*param[1]))/frac+((1-(*param[0])-(*param[1]))/frac)*randnumbers::incomplete_beta(a,b,(*response[3])) );
    }


double DISTR_betainf1_tau::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);

  if (counter==0)
    {
    set_worklin();
    }

  double taup = exp((*linpred));

  double l;

    if ((*response)==1)
     {
         l = log(taup) - log(1+taup);
     }
     else
     {
        l =  - log(1+taup);
     }

  modify_worklin();

  return l;

  }


void DISTR_betainf1_tau::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double taup = exp(*linpred);

    double hilfs = 1+ taup;

    double nu = -taup/hilfs;

    if ((*response)==1)
    {
        nu += 1;
    }

    *workingweight = taup/(pow(hilfs,2));

    *workingresponse = *linpred + nu/(*workingweight);

        if (compute_like)
      {

        if ((*response)==0)
        {
            like += log(taup) - log(1+taup);
        }
        else
        {
            like -=   log(1+taup);
        }
        }


  modify_worklin();

  }


void DISTR_betainf1_tau::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double exp_lin_taup = exp(*linpred[predstart_mumult+2]);
  double exp_lin_mu = exp(*linpred[predstart_mumult]);
  double hilfs = (1+exp_lin_taup);
  *mu = (1-(exp_lin_taup)/hilfs)*(exp_lin_mu/(1+exp_lin_mu))+ exp_lin_taup/hilfs;
  }


void DISTR_betainf1_tau::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (nu): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf1_tau::update_end(void)
  {


  // helpmat1 stores exp(eta_nu)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();


  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }



//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_betainf0_nu -------------------------
//------------------------------------------------------------------------------


DISTR_betainf0_nu::DISTR_betainf0_nu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,0,w)
  {
  family = "Beta Zero Inflated Distribution - nu";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "nu";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_betainf0_nu::DISTR_betainf0_nu(const DISTR_betainf0_nu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf0_nu & DISTR_betainf0_nu::operator=(
                            const DISTR_betainf0_nu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_betainf0_nu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_tau
   // *linpred[1] = eta_nu
   // *linpred[2] = eta_sigma2
   // *linpred[3] = eta_mu

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double nup = exp(*linpred[2]);
     double sigma_2 = exp(*linpred[1])/(1+exp(*linpred[1]));
     double mu = exp(*linpred[0])/(1+exp(*linpred[0]));
     double help = (1-sigma_2)/sigma_2;
     double one_minus_mu_help = (1-mu)*help;
     double mu_help = mu*help;
     double one_nup = (1+nup);

     double l;

     if ((*response[2])==0)
     {
         l = log(nup) - log(one_nup);
     }
      else
       l = (mu_help-1)*log(*response[1]) +
			(one_minus_mu_help-1)*log(1-(*response[1]))-
			randnumbers::lngamma_exact(mu_help)-
			randnumbers::lngamma_exact(one_minus_mu_help)+
			randnumbers::lngamma_exact(help)- log(one_nup);


    *deviance = -2*l;
    }

  }


double DISTR_betainf0_nu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_betainf0_nu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[2]);
  }
 double DISTR_betainf0_nu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_betainf0_nu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
//    double a =  (*param[2])*(*param[3]);
//    double b = (*param[2])*(1-(*param[3]));
 //   double frac = 1 + (*param[0]) + (*param[1]);
    return 0;
//
//    return ( ((*param[0])+(*param[1]))/frac+((1-(*param[0])-(*param[1]))/frac)*randnumbers::incomplete_beta(a,b,(*response[3])) );
    }


double DISTR_betainf0_nu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);

  if (counter==0)
    {
    set_worklin();
    }

  double nup = exp((*linpred));

  double l;

    if ((*response)==0)
     {
         l = log(nup) - log(1+nup);
     }
     else
     {
        l =  - log(1+nup);
     }

  modify_worklin();

  return l;

  }


void DISTR_betainf0_nu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double nup = exp(*linpred);

    double hilfs = 1+ nup;

    double nu = -nup/hilfs;

    if ((*response)==0)
    {
        nu += 1;
    }

    *workingweight = nup/(pow(hilfs,2));

    *workingresponse = *linpred + nu/(*workingweight);

        if (compute_like)
      {

        if ((*response)==0)
        {
            like += log(nup) - log(1+nup);
        }
        else
        {
            like -=   log(1+nup);
        }
        }


  modify_worklin();

  }


void DISTR_betainf0_nu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double exp_lin_nup = exp(*linpred[predstart_mumult+2]);
  double exp_lin_mu = exp(*linpred[predstart_mumult]);
  double hilfs = (1+exp_lin_nup);
  *mu = (1-(exp_lin_nup)/hilfs)*(exp_lin_mu/(1+exp_lin_mu));
  }


void DISTR_betainf0_nu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (nu): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf0_nu::update_end(void)
  {


  // helpmat1 stores exp(eta_nu)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();


  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_t_df ----------------------------------
//------------------------------------------------------------------------------


DISTR_t_df::DISTR_t_df(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "t-Distribution - Degrees of Freedom";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "df";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_t_df::DISTR_t_df(const DISTR_t_df & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_t_df & DISTR_t_df::operator=(
                            const DISTR_t_df & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_t_df::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_t_df::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[0]);
  }
double DISTR_t_df::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = exp(eta_sigma2);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

    double degf = exp(*linpred);
    double arg = (degf+1)/2;
    double arg1 = (degf)/2;
    double l;

     l = (randnumbers::lngamma_exact(arg))-(randnumbers::lngamma_exact(arg1))-
            0.5*log(degf) - (arg)*log(1+pow((*response)-(*worklin[1]),2)/((*worktransformlin[0])*degf));


  modify_worklin();

  return l;

  }

void DISTR_t_df::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = exp(eta_sigma2);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);


  if (counter==0)
    {
    set_worklin();
    }

    double degf = exp((*linpred));
    double denom = (*worktransformlin[0])*degf+pow((*response)-(*worklin[1]),2);
    double frac = pow((*response)-(*worktransformlin[1]),2)/((*worktransformlin[0])*degf);
    double arg = (degf+1)/2;
    double arg1 = degf/2;
    double nu = -0.5 + arg*pow(((*response)-(*worklin[1])),2)/denom + 0.5*degf*(randnumbers::digamma_exact(arg)-randnumbers::digamma_exact(arg1)-log(1+frac));

  //  *workingweight = pow(nu,2);
    *workingweight =  - 0.25*pow(degf,2)*(randnumbers::trigamma_exact(arg)-randnumbers::trigamma_exact(arg1)) - degf/((degf+1)) + degf/(2*(degf+3));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += (randnumbers::lngamma_exact(arg))-(randnumbers::lngamma_exact(arg1))-
            0.5*log(degf) - (arg)*log(1+frac);

      }

  modify_worklin();

  }


void DISTR_t_df::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (df): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_t_df::update_end(void)
  {

  // helpmat1 stores tau

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_t_sigma2 ------------------------------
//------------------------------------------------------------------------------


DISTR_t_sigma2::DISTR_t_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "t-Distribution - sigma2";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_t_sigma2::DISTR_t_sigma2(const DISTR_t_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_t_sigma2 & DISTR_t_sigma2::operator=(
                            const DISTR_t_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_t_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

  void DISTR_t_sigma2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  // this was for Alex paper
  //*param = pow(exp(*linpred[1]), 0.5);
  *param = exp(*linpred[1]);
  }


double DISTR_t_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of n equation
  // *worktransformlin[0] = exp(eta_n);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = exp(*linpred);

  double l;

     l =  - 0.5*log(sigma_2) - (((*worktransformlin[0])+1)/(2))*log(1+pow((*response)-(*worklin[1]),2)/(sigma_2*(*worktransformlin[0])));


  modify_worklin();

  return l;

  }

void DISTR_t_sigma2::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = exp(eta_tau);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

    double sig = exp(*linpred);
    double denom = sig*(*worktransformlin[0])+pow(((*response)-(*worklin[1])),2);

    double nu = -0.5 + (((*worktransformlin[0])+1)/(2)*pow(((*response)-(*worklin[1])),2))/denom;

  //  *workingweight = pow(nu,2);

 //  *workingweight = 0.5*sig*(*worktransformlin[0])*((*worktransformlin[0])+1)*pow((*response)-(*worktransformlin[1]),2)/pow(denom,2);

    *workingweight = (*worktransformlin[0])/(2*((*worktransformlin[0])+3));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - 0.5*log(sig) - (((*worktransformlin[0])+1)/(2))*log(1+pow((*response)-(*worklin[1]),2)/(sig*(*worktransformlin[0])));

      }

  modify_worklin();

  }


void DISTR_t_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_t_sigma2::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_t_mu --------------------------------
//------------------------------------------------------------------------------


void DISTR_t_mu::check_errors(void)
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


DISTR_t_mu::DISTR_t_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "t Distribution - mu";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
   // linpredminlimit=-10;
  //linpredmaxlimit=15;
  check_errors();
  }


DISTR_t_mu::DISTR_t_mu(const DISTR_t_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_t_mu & DISTR_t_mu::operator=(
                            const DISTR_t_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_t_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_df
   // *linpred[1] = eta_sigma2
   // *linpred[2] = eta_mu

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
	 double degf = exp(*linpred[0]);
     double sigma_2 = exp(*linpred[1]);
     double mu = (*linpred[2]);
     double hilfs = (degf+1)/2;
     double hilfs2 = degf/2;

     double l;

       l = randnumbers::lngamma_exact((hilfs)) -log(sqrt(PI)) -randnumbers::lngamma_exact(hilfs2) -0.5*log(degf) - 0.5*log(sigma_2) -
            hilfs*log(1+(pow((*response[2])-mu,2))/(sigma_2*degf));


    *deviance = -2*l;
    }

  }


double DISTR_t_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

  void DISTR_t_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2]);
  }

 double DISTR_t_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_t_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
 //   double x = ((*response[2])-(*param[2]))/pow((*param[1]),0.5);
 //   double u = gsl_cdf_tdist_P(x, (*param[0]));
   // double x = (*param[0])/(pow((((*response[2])-(*param[2]))/pow((*param[1]))),2)+(*param[0]));
    return 0;
//    return ( 1- 0.5*randnumbers::incomplete_beta(a,b,x));
    }

double DISTR_t_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of degf equation
  // *worktransformlin[0] = degf;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double hilfs =((*worktransformlin[0])+1)/2;
  double l;

     l =  - hilfs*log(1+(pow((*response)-mu,2))/((*worktransformlin[1])*(*worktransformlin[0])));

  modify_worklin();

  return l;

  }


void DISTR_t_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of degf equation
  // *worktransformlin[0] = degf;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);
    double hilfs = ((*worktransformlin[0])+1)/2;
//    double denom1 = (*worktransformlin[1])*(*worktransformlin[0])-pow((*response)-mu,2);
    double denom2 = (*worktransformlin[1])*(*worktransformlin[0])+pow((*response)-mu,2);

    double nu = ((*worktransformlin[0])+1)*((*response)-mu)/denom2;

   // *workingweight = pow(nu,2);
    *workingweight = -1/(*worktransformlin[1]) + 2*((*worktransformlin[0])+2)/((*worktransformlin[1])*((*worktransformlin[0])+3));
    //*workingweight = (((*worktransformlin[0])+1)*denom1)/(denom2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=   - hilfs*log(1+(pow((*response)-mu,2))/((*worktransformlin[1])*(*worktransformlin[0])));

      }


  modify_worklin();

  }


void DISTR_t_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

   if(exp((*linpred[predstart_mumult]))>1)
   {
       *mu = (*linpred[predstart_mumult+2]);
   } else
   {
       *mu = 0;
   }



  }


void DISTR_t_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_t_mu::update_end(void)
  {


  // helpmat1 stores exp(eta_a)

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



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_invgaussian_sigma2 --------------------
//------------------------------------------------------------------------------


DISTR_invgaussian_sigma2::DISTR_invgaussian_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Inverse Gaussian Distribution - sigma2";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_invgaussian_sigma2::DISTR_invgaussian_sigma2(const DISTR_invgaussian_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_invgaussian_sigma2 & DISTR_invgaussian_sigma2::operator=(
                            const DISTR_invgaussian_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_invgaussian_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_invgaussian_sigma2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[0]));
  }

double DISTR_invgaussian_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = exp((*linpred));

  double l;

     l = -0.5*log(sigma_2)-pow((((*response))-(*worktransformlin[0])),2)/(2*(*response)*pow((*worktransformlin[0]),2)*sigma_2);


  modify_worklin();

  return l;

  }

void DISTR_invgaussian_sigma2::compute_iwls_wweightschange_weightsone(
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
    }

    double sigma_2 = exp(*linpred);


    double nu = -0.5 + (pow(((*response)-(*worktransformlin[0])),2))/(2*(*response)*(pow((*worktransformlin[0]),2))*sigma_2);



    *workingweight = 0.5;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((((*response))-(*worktransformlin[0])),2)/(2*(*response)*pow((*worktransformlin[0]),2)*sigma_2);

      }

  modify_worklin();

  }


void DISTR_invgaussian_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_invgaussian_sigma2::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_invgaussian_mu ----------------------
//------------------------------------------------------------------------------
void DISTR_invgaussian_mu::check_errors(void)
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

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
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


DISTR_invgaussian_mu::DISTR_invgaussian_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Inverse Gaussian Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_invgaussian_mu::DISTR_invgaussian_mu(const DISTR_invgaussian_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_invgaussian_mu & DISTR_invgaussian_mu::operator=(
                            const DISTR_invgaussian_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_invgaussian_mu::compute_deviance_mult(vector<double *> response,
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
     double sigma_2 = exp(*linpred[0]);
     double mu = exp(*linpred[1]);

     double l;

       l = -0.5*log(2*PI)-0.5*log(sigma_2)-1.5*log((*response[1]))-(1/(2*(*response[1])))*pow((((*response[0]))-mu),2)/(pow(mu,2)*sigma_2);


    *deviance = -2*l;
    }

  }


double DISTR_invgaussian_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_invgaussian_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]);
  }
 double DISTR_invgaussian_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_invgaussian_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = pow(1/((*param[0])*(*response[1])),0.5);
    double arg0 = (*response[1])/((*param[1]));
    double arg1 = arg*(arg0-1);
    double arg2 = -arg*(arg0+1);
    double arg3 = 2/((*param[0])*(*param[0]));

    return (randnumbers::Phi2(arg1)+exp(arg3)*randnumbers::Phi2(arg2));
    }


double DISTR_invgaussian_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);

  double l;

     l = -pow((((*response))-mu),2)/(2*(*response)*pow(mu,2)*(*worktransformlin[0]));

  modify_worklin();

  return l;

  }


void DISTR_invgaussian_mu::compute_iwls_wweightschange_weightsone(
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

    double mu = exp(*linpred);

    double nu = ((*response)-mu)/(pow(mu,2)*(*worktransformlin[0]));

    *workingweight = 1/(mu*(*worktransformlin[0]));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((((*response))-mu),2)/(2*(*response)*pow(mu,2)*(*worktransformlin[0]));

      }


  modify_worklin();

  }


void DISTR_invgaussian_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

  *mu = exp((*linpred[predstart_mumult+1]));

  }


void DISTR_invgaussian_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_invgaussian_mu::update_end(void)
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
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_betainf_mu --------------------------
//------------------------------------------------------------------------------


DISTR_betainf_mu::DISTR_betainf_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta Zero One Inflated Distribution - mu";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "mu";
    linpredminlimit=-10;
  linpredmaxlimit=10;
  }


DISTR_betainf_mu::DISTR_betainf_mu(const DISTR_betainf_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf_mu & DISTR_betainf_mu::operator=(
                            const DISTR_betainf_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_betainf_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_betainf_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  double arg = exp(*linpred[0]);
  *param = arg/(1+arg);
  }

double DISTR_betainf_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = (1-sigma2)/sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double exp_lin = exp((*linpred));
  double mu = exp_lin/(1+exp_lin);
  double mu_worktrans = mu*(*worktransformlin[0]);
  double one_minus_mu_worktrans = (1-mu)*(*worktransformlin[0]);

  double l;

     l = mu*(*worktransformlin[0])*log(*response) +
		 (1-mu)*(*worktransformlin[0])*log(1-(*response))-
		 randnumbers::lngamma_exact(mu_worktrans)-
		 randnumbers::lngamma_exact(one_minus_mu_worktrans);

  modify_worklin();

  return l;

  }


void DISTR_betainf_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = (1-sigma2)/sigma2;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double exp_lin = exp(*linpred);

    double mu = exp_lin/(1+exp_lin);

    double mu_worktrans = mu*(*worktransformlin[0]);

    double one_minus_mu_worktrans = (1-mu)*(*worktransformlin[0]);

    double nu = mu*one_minus_mu_worktrans*log((*response))-mu*one_minus_mu_worktrans*log((1-(*response)))+
				mu*one_minus_mu_worktrans*(randnumbers::digamma_exact(one_minus_mu_worktrans)-randnumbers::digamma_exact(mu_worktrans));

    *workingweight = pow((*worktransformlin[0]),2)*pow(mu,2)*pow((1-mu),2)*(randnumbers::trigamma_exact(one_minus_mu_worktrans)+randnumbers::trigamma_exact(mu_worktrans));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += mu_worktrans*log((*response))+one_minus_mu_worktrans*log(1-(*response))-
		        randnumbers::lngamma_exact(mu_worktrans)-randnumbers::lngamma_exact((one_minus_mu_worktrans));

      }


  modify_worklin();

  }



void DISTR_betainf_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf_mu::update_end(void)
  {


  // helpmat1 stores exp(eta_mu)/(1+exp(eta_mu))

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin)/(1+exp(*worklin));
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_betainf_sigma2 ------------------------
//------------------------------------------------------------------------------


DISTR_betainf_sigma2::DISTR_betainf_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta Zero One Inflated Distribution - sigma2";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=10;
  }


DISTR_betainf_sigma2::DISTR_betainf_sigma2(const DISTR_betainf_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf_sigma2 & DISTR_betainf_sigma2::operator=(
                            const DISTR_betainf_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_betainf_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_betainf_sigma2::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = exp(*linpred[1]);
  *param = arg/(1+arg);
  }

double DISTR_betainf_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu)/(1+exp(eta_mu));

  if (counter==0)
    {
    set_worklin();
    }

  double exp_lin = exp((*linpred));
  double sigma_2 = exp_lin/(1+exp_lin);
  double help = (1-sigma_2)/sigma_2;
  double mu_help = (*worktransformlin[0])*help;
  double one_minus_mu_help = (1-(*worktransformlin[0]))*help;

  double l;

     l = mu_help*log(*response) +
         one_minus_mu_help*log(1-(*response)) -
         randnumbers::lngamma_exact(mu_help) - randnumbers::lngamma_exact(one_minus_mu_help) +
         randnumbers::lngamma_exact(help);


  modify_worklin();

  return l;

  }

void DISTR_betainf_sigma2::compute_iwls_wweightschange_weightsone(
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

    double exp_lin = exp(*linpred);

    double sigma_2 = exp_lin/(1+exp_lin);

    double help = (1-sigma_2)/sigma_2;
    double mu_help = (*worktransformlin[0])*help;
    double one_minus_mu_help = (1-(*worktransformlin[0]))*help;

    double nu = -help*( -(*worktransformlin[0])*randnumbers::digamma_exact(mu_help)
                       - (1-(*worktransformlin[0]))*randnumbers::digamma_exact(one_minus_mu_help) +
                       randnumbers::digamma_exact(help) + (*worktransformlin[0])*log(*response)
                       + (1-(*worktransformlin[0]))*log(1-(*response)) );



    *workingweight =  pow(help,2)*( pow((1-(*worktransformlin[0])),2)*randnumbers::trigamma_exact(one_minus_mu_help)+
                                   pow(((*worktransformlin[0])),2)*randnumbers::trigamma_exact(mu_help) -
                                   randnumbers::trigamma_exact(help) );

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += mu_help*log(*response) +
         one_minus_mu_help*log(1-(*response)) -
         randnumbers::lngamma_exact(mu_help) - randnumbers::lngamma_exact(one_minus_mu_help) +
         randnumbers::lngamma_exact(help);

      }

  modify_worklin();

  }


void DISTR_betainf_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf_sigma2::update_end(void)
  {

  // helpmat1 stores (1-sigma2)/sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double exp_lin;
  double sigma_2;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    exp_lin = exp(*worklin);
    sigma_2 = exp_lin/(1+exp_lin);
    *pmu = (1-sigma_2)/sigma_2;
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_betainf_nu --------------------------
//------------------------------------------------------------------------------


DISTR_betainf_nu::DISTR_betainf_nu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta Zero One Inflated Distribution - nu";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "nu";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_betainf_nu::DISTR_betainf_nu(const DISTR_betainf_nu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf_nu & DISTR_betainf_nu::operator=(
                            const DISTR_betainf_nu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_betainf_nu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_betainf_nu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[2]);
  }

double DISTR_betainf_nu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = exp(eta_tau);

  if (counter==0)
    {
    set_worklin();
    }

  double nup = exp((*linpred));

  double l;

    if ((*response)==0)
     {
         l = log(nup) - log(1+(*worktransformlin[0])+nup);
     }
     else
     {
        l =  - log(1+(*worktransformlin[0])+nup);
     }

  modify_worklin();

  return l;

  }


void DISTR_betainf_nu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = exp(eta_tau);

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double nup = exp(*linpred);

    double hilfs = 1+ (*worktransformlin[0]) + nup;

    double nu = -nup/hilfs;

    if ((*response)==0)
    {
        nu += 1;
    }

    *workingweight = nup*(1+(*worktransformlin[0]))/(pow(hilfs,2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        if ((*response)==0)
        {
            like += log(nup) - log(1+(*worktransformlin[0])+nup);
        }
        else
        {
            like +=  - log(1+(*worktransformlin[0])+nup);
        }
      }


  modify_worklin();

  }



void DISTR_betainf_nu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (nu): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf_nu::update_end(void)
  {


  // helpmat1 stores exp(eta_tau)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();


  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_betainf_tau -------------------------
//------------------------------------------------------------------------------


DISTR_betainf_tau::DISTR_betainf_tau(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta Zero One Inflated Distribution - tau";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "tau";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_betainf_tau::DISTR_betainf_tau(const DISTR_betainf_tau & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_betainf_tau & DISTR_betainf_tau::operator=(
                            const DISTR_betainf_tau & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_betainf_tau::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_tau
   // *linpred[1] = eta_nu
   // *linpred[2] = eta_sigma2
   // *linpred[3] = eta_mu

   if (*weight[3] == 0)
     *deviance=0;
   else
     {
     double tau =  exp(*linpred[3]);
     double nup = exp(*linpred[2]);
     double sigma_2 = exp(*linpred[1])/(1+exp(*linpred[1]));
     double mu = exp(*linpred[0])/(1+exp(*linpred[0]));
     double help = (1-sigma_2)/sigma_2;
     double one_minus_mu_help = (1-mu)*help;
     double mu_help = mu*help;
     double one_nup_tau = (1+nup+tau);

     double l;

     if ((*response[3])==0)
     {
         l = log(nup) - log(one_nup_tau);
     }
     else if ((*response[3])==1)
     {
        l = log(tau) - log(one_nup_tau);
     }
     else
       l = (mu_help-1)*log(*response[3]) +
			(one_minus_mu_help-1)*log(1-(*response[3]))-
			randnumbers::lngamma_exact(mu_help)-
			randnumbers::lngamma_exact(one_minus_mu_help)+
			randnumbers::lngamma_exact(help)- log(one_nup_tau);


    *deviance = -2*l;
    }

  }


double DISTR_betainf_tau::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_betainf_tau::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[3]);
  }
 double DISTR_betainf_tau::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_betainf_tau::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
//    double a =  (*param[2])*(*param[3]);
//    double b = (*param[2])*(1-(*param[3]));
//    double frac = 1 + (*param[0]) + (*param[1]);
    return 0;
//
//    return ( ((*param[0])+(*param[1]))/frac+((1-(*param[0])-(*param[1]))/frac)*randnumbers::incomplete_beta(a,b,(*response[3])) );
    }


double DISTR_betainf_tau::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);

  if (counter==0)
    {
    set_worklin();
    }

  double tau = exp((*linpred));

  double l;

    if ((*response)==1)
     {
         l = log(tau) - log(1+(*worktransformlin[0])+tau);
     }
     else
     {
        l =  - log(1+(*worktransformlin[0])+tau);
     }

  modify_worklin();

  return l;

  }


void DISTR_betainf_tau::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double tau = exp(*linpred);

    double hilfs = 1+ (*worktransformlin[0]) + tau;

    double nu = -tau/hilfs;

    if ((*response)==1)
    {
        nu += 1;
    }

    *workingweight = tau*(1+(*worktransformlin[0]))/(pow(hilfs,2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        if ((*response)==1)
        {
            like += log(tau) - log(1+(*worktransformlin[0])+tau);
        }
        else
        {
            like +=  - log(1+(*worktransformlin[0])+tau);
        }
      }


  modify_worklin();

  }


void DISTR_betainf_tau::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double exp_lin_tau = exp(*linpred[predstart_mumult+3]);
  double exp_lin_nup = exp(*linpred[predstart_mumult+2]);
  double exp_lin_mu = exp(*linpred[predstart_mumult]);
  double hilfs = (1+exp_lin_tau+exp_lin_nup);
  *mu = (1-(exp_lin_tau+exp_lin_nup)/hilfs)*(exp_lin_mu/(1+exp_lin_mu))+exp_lin_tau/(hilfs);
  }


void DISTR_betainf_tau::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (tau): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_betainf_tau::update_end(void)
  {


  // helpmat1 stores exp(eta_nu)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();


  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_pareto_p ------------------------------
//------------------------------------------------------------------------------


DISTR_pareto_p::DISTR_pareto_p(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Pareto Distribution - p";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "p";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_pareto_p::DISTR_pareto_p(const DISTR_pareto_p & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_pareto_p & DISTR_pareto_p::operator=(
                            const DISTR_pareto_p & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_pareto_p::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_pareto_p::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[0]);
  }

double DISTR_pareto_p::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of b equation
  // *worktransformlin[0] = exp(eta_b);

  if (counter==0)
    {
    set_worklin();
    }

  double p = exp((*linpred));

  double l;

     l = log(p) + p*log(*worktransformlin[0]) - (p)*log((*response)+(*worktransformlin[0]));


  modify_worklin();

  return l;

  }

void DISTR_pareto_p::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of b equation
  // *worktransformlin[0] = exp(eta_b);

  if (counter==0)
    {
    set_worklin();
    }

    double p = exp((*linpred));

    double nu = 1  + p*log(*worktransformlin[0]) - (p)*log((*response)+(*worktransformlin[0]));



    *workingweight = 1;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  log(p) + p*log(*worktransformlin[0]) - (p)*log((*response)+(*worktransformlin[0]));

      }

  modify_worklin();

  }


void DISTR_pareto_p::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (p): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_pareto_p::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_pareto_b ----------------------------
//------------------------------------------------------------------------------

void DISTR_pareto_b::check_errors(void)
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

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
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

DISTR_pareto_b::DISTR_pareto_b(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Pareto Distribution - b";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "b";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_pareto_b::DISTR_pareto_b(const DISTR_pareto_b & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_pareto_b & DISTR_pareto_b::operator=(
                            const DISTR_pareto_b & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_pareto_b::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_p
   // *linpred[1] = eta_b

   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double p = exp(*linpred[0]);
     double b = exp(*linpred[1]);

     double l;

       l =   log(p) + p*log(b) - (p+1)*log((*response[1])+b);


    *deviance = -2*l;
    }

  }


double DISTR_pareto_b::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_pareto_b::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]);
  }

 double DISTR_pareto_b::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_pareto_b::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double b = *param[1];
    double p = *param[0];
    double res = 1- pow(b,p)*pow((*response[1])+b,-p);
    return ( res );
    }

double DISTR_pareto_b::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = p;

  if (counter==0)
    {
    set_worklin();
    }

  double b = exp(*linpred);

  double l;

     l =  (*worktransformlin[0])*log(b) - ((*worktransformlin[0])+1)*log((*response)+b) ;

  modify_worklin();


  return l;

  }


void DISTR_pareto_b::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = p;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double b = exp(*linpred);

    double nu = (*worktransformlin[0]) - (((*worktransformlin[0])+1)*b)/((*response)+b);

    *workingweight = ((*worktransformlin[0]))/((*worktransformlin[0])+2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += (*worktransformlin[0])*log(b) - ((*worktransformlin[0])+1)*log((*response)+b) ;

      }


  modify_worklin();

  }


void DISTR_pareto_b::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

   double p = exp((*linpred[predstart_mumult]));
   double b = exp((*linpred[predstart_mumult+1]));
  *mu = b/(p-1);

  }


void DISTR_pareto_b::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (b): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_pareto_b::update_end(void)
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
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_dagum_p -------------------------------
//------------------------------------------------------------------------------

DISTR_dagum_p::DISTR_dagum_p(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Dagum Distribution - p";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "p";
  linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_dagum_p::DISTR_dagum_p(const DISTR_dagum_p & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_dagum_p & DISTR_dagum_p::operator=(
                            const DISTR_dagum_p & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_dagum_p::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

 //  double test = *linpred;
// compute cdf (might work more efficiently)
  double res,a,b,p;
  p = exp(*linpredp);
  a = *worktransformlin[1];//exp(*linpred[0]);
  b = *worktransformlin[0];
  res = pow(1+pow(resp/b,-a),-p);

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_dagum_p::cdf(const double & resp, const double & linpred)
  {
  double res,a,b,p;
  p = exp(linpred);
  a = *worktransformlin[1];//exp(*linpred[0]);
  b = *worktransformlin[0];
  res = pow(1+pow(resp/b,-a),-p);
  return res;
  }


double DISTR_dagum_p::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_dagum_p::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = exp(*linpred[copulaoffset+0]);
  }

double DISTR_dagum_p::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of b equation
  // *worktransformlin[0] = exp(eta_b);
  // *worklin[1] = linear predictor of a equation
  // *worktransformlin[1] = exp(eta_a);

  if (counter==0)
    {
    set_worklin();
    }

  double p = exp((*linpred));
  double hilfs = pow((*response)/(*worktransformlin[0]),(*worktransformlin[1]));

  double l;

     l = log(p) + (*worktransformlin[1])*p*log((*response)) - (*worktransformlin[1])*p*log((*worktransformlin[0]))
        -p*log(1+hilfs);

  if(optionsp->copula)
    {
    //implement loglik for copula models, i.e. add part logc
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();

  return l;

  }

void DISTR_dagum_p::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of b equation
  // *worktransformlin[0] = exp(eta_b);
  // *worklin[1] = linear predictor of a equation
  // *worktransformlin[1] = exp(eta_a);


  if (counter==0)
    {
    set_worklin();
    }

    double p = exp((*linpred));
    double hilfs = pow((*response)/(*worktransformlin[0]),(*worktransformlin[1]));

    double nu = 1 + (*worktransformlin[1])*p*log((*response)) - (*worktransformlin[1])*p*log((*worktransformlin[0]))
                -p*log(1+hilfs);

    *workingweight = 1;

    if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
//    double lyb = log(*response/(*worktransformlin[0]));
    double ybpma = pow(*response/(*worktransformlin[0]),-(*worktransformlin[1]));
    double ybpmap = pow(1+ybpma,-p);
    double dF = -p*log(1+ybpma)*ybpmap;
    double ddF = p*log(1+ybpma)*ybpmap*(p*log(1+ybpma)-1);
  /*  if(copularotate)
      {
      dF = -dF;
      ddF = -ddF;
      }*/
    nu += logcandderivs[1]*dF;

   // *workingweight = (*worktransformlin[0])*(*worktransformlin[0])*hilfs1-logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
   *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;

    if (*workingweight <=0)
      *workingweight = 0.001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += log(p) + (*worktransformlin[1])*p*log((*response)) - (*worktransformlin[1])*p*log((*worktransformlin[0]))
        -p*log(1+hilfs);

      }

  modify_worklin();

  }


void DISTR_dagum_p::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (p): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_dagum_p::update_end(void)
  {

  // helpmat1 stores tau

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_dagum_b -------------------------------
//------------------------------------------------------------------------------


DISTR_dagum_b::DISTR_dagum_b(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Dagum Distribution - b";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "b";
  linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_dagum_b::DISTR_dagum_b(const DISTR_dagum_b & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_dagum_b & DISTR_dagum_b::operator=(
                            const DISTR_dagum_b & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_dagum_b::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

 //  double test = *linpred;
// compute cdf (might work more efficiently)
  double res,a,b,p;
  b = exp(*linpredp);
  a = *worktransformlin[1];//exp(*linpred[0]);
  p = *worktransformlin[0];
  res = pow(1+pow(resp/b,-a),-p);

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_dagum_b::cdf(const double & resp, const double & linpred)
  {
  double res,a,b,p;
  b = exp(linpred);
  a = *worktransformlin[1];//exp(*linpred[0]);
  p = *worktransformlin[0];
  res = pow(1+pow(resp/b,-a),-p);
  return res;
  }


double DISTR_dagum_b::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_dagum_b::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = exp(*linpred[copulaoffset+1]);
  }

double DISTR_dagum_b::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = exp(eta_p);
  // *worklin[1] = linear predictor of a equation
  // *worktransformlin[1] = exp(eta_a);

  if (counter==0)
    {
    set_worklin();
    }

  double b = exp((*linpred));
  double hilfs = pow((*response)/b,(*worktransformlin[1]));
  double l;

     l = - (*worktransformlin[1])*(*worktransformlin[0])*log(b) - ((*worktransformlin[0])+1)*log(1+hilfs);

  if(optionsp->copula)
    {
    //implement loglik for copula models, i.e. add part logc
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();

  return l;

  }

void DISTR_dagum_b::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = exp(eta_p);
  // *worklin[1] = linear predictor of a equation
  // *worktransformlin[1] = exp(eta_a);

  if (counter==0)
    {
    set_worklin();
    }

    double b = exp((*linpred));
    double hilfs = pow((*response)/b,(*worktransformlin[1]));

    double nu = (*worktransformlin[1]) - (((*worktransformlin[0])+1)*(*worktransformlin[1]))/(1+hilfs) ;


    *workingweight = (pow((*worktransformlin[1]),2)*(*worktransformlin[0]))/((*worktransformlin[0])+2);
   // *workingweight = (((*worktransformlin[0])+1)*pow((*worktransformlin[1]),2)*hilfs)/pow((1+hilfs),2);

   if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
//    double lyb = log(*response/b);
    double ybpa = pow(*response/b,(*worktransformlin[1]));
    double ybpmap = pow(1+pow(*response/b,-(*worktransformlin[1])),-(*worktransformlin[0])-1);
    double ybpmap2 = pow(1+pow(*response/b,-(*worktransformlin[1])),-(*worktransformlin[0]));
    double dF = -(*worktransformlin[1])*(*worktransformlin[0])*ybpmap/ybpa;
    double ddF = (*worktransformlin[1])*(*worktransformlin[1])*(*worktransformlin[0])*((*worktransformlin[0])-ybpa)*ybpmap2/((ybpa+1)*(ybpa+1));
  /*  if(copularotate)
      {
      dF = -dF;
      ddF = -ddF;
      }*/
    nu += logcandderivs[1]*dF;

   // *workingweight = (*worktransformlin[0])*(*worktransformlin[0])*hilfs1-logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
   *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;

    if (*workingweight <=0)
      *workingweight = 0.001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - (*worktransformlin[1])*(*worktransformlin[0])*log(b) - ((*worktransformlin[0])+1)*log(1+hilfs);

      }

  modify_worklin();

  }


void DISTR_dagum_b::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (b): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_dagum_b::update_end(void)
  {

  // helpmat1 stores b

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_dagum_a -----------------------------
//------------------------------------------------------------------------------


void DISTR_dagum_a::check_errors(void)
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

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
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


DISTR_dagum_a::DISTR_dagum_a(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Dagum Distribution - a";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "a";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_dagum_a::DISTR_dagum_a(const DISTR_dagum_a & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_dagum_a & DISTR_dagum_a::operator=(
                            const DISTR_dagum_a & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_dagum_a::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

 //  double test = *linpred;
// compute cdf (might work more efficiently)
  double res,a,b,p;
  a = exp(*linpredp);
  b = *worktransformlin[1];//exp(*linpred[0]);
  p = *worktransformlin[0];
  res = pow(1+pow(resp/b,-a),-p);

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_dagum_a::cdf(const double & resp, const double & linpred)
  {
  double res,a,b,p;
  a = exp(linpred);
  b = *worktransformlin[1];//exp(*linpred[0]);
  p = *worktransformlin[0];
  res = pow(1+pow(resp/b,-a),-p);
  return res;
  }

double DISTR_dagum_a::cdf(const double & resp, vector<double *>  linpred)
  {
  //linpred has always (also with copula) size=3
  double res,a,b,p;
  a = exp(*linpred[2]);
  b = exp(*linpred[1]);
  p = exp(*linpred[0]);
  res = pow(1+pow(resp/b,-a),-p);
  return res;
  }

void DISTR_dagum_a::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_p
   // *linpred[1] = eta_b
   // *linpred[2] = eta_a

  //linpred,weight,response has always (also with copula) size=3
   if (*weight[2] == 0)
     *deviance=0;
   else
     {
	 double p = exp(*linpred[0]);
     double b = exp(*linpred[1]);
     double a = exp(*linpred[2]);
     double hilfs = pow((*response[2])/b,a);

     double l;

       l = log(a) + log(p) +(a*p-1)*log((*response[2])) - a*p*log(b) - (p+1)*log(1+hilfs);


    *deviance = -2*l;
    }

  }


double DISTR_dagum_a::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_dagum_a::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = exp(*linpred[copulaoffset+2]);
  }

 double DISTR_dagum_a::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_dagum_a::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    //weight, response and param has size>1 if copula model is specified!!
    return ( pow((1+pow((*response[copulaoffset+1])/(*param[copulaoffset+1]),-(*param[copulaoffset+2]))),-(*param[copulaoffset+0])) );
    }

double DISTR_dagum_a::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = p;
  // *worklin[1] = linear predictor of b equation
  // *worktransformlin[1] = b;

  if (counter==0)
    {
    set_worklin();
    }

  double a = exp(*linpred);
  double hilfs = pow((*response)/(*worktransformlin[1]),a);
  double l;

     l =  log(a) +(a*(*worktransformlin[0]))*log((*response)) - a*(*worktransformlin[0])*log((*worktransformlin[1]))
            - ((*worktransformlin[0])+1)*log(1+hilfs);

  if(optionsp->copula)
    {
    //implement loglik for copula models, i.e. add part logc
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }
  modify_worklin();

  return l;

  }


void DISTR_dagum_a::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of p equation
  // *worktransformlin[0] = p;
  // *worklin[1] = linear predictor of b equation
  // *worktransformlin[1] = b;


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double a = exp((*linpred));

    double hilfs = pow((*response)/(*worktransformlin[1]),a);

    double nu = 1 + a*(*worktransformlin[0])*log((*response)/(*worktransformlin[1]))
                - (((*worktransformlin[0])+1)*a*hilfs*log((*response)/(*worktransformlin[1])))/(1+hilfs);

    *workingweight = 1 + (((*worktransformlin[0])+1)*pow(a,2)*hilfs*pow(log((*response)/(*worktransformlin[1])),2))/pow((1+hilfs),2);

    if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
    double ybpa = pow((*response/(*worktransformlin[1])),a);
    double lyb = log(*response/(*worktransformlin[1]));
    double ybpma = pow(*response/(*worktransformlin[1]),-a);
    double ybpmap = pow(1+ybpma,-(*worktransformlin[0])-1);
    double ybpmap2 = pow(1+ybpma,-(*worktransformlin[0]));
    double dF = (*worktransformlin[0])*lyb*a*ybpma*ybpmap;
    double ddF = a*(*worktransformlin[0])*lyb*ybpmap2*(a*lyb*((*worktransformlin[0])-ybpa)+ybpa+1)/pow(1+ybpa,2);
  /*  if(copularotate)
      {
      dF = -dF;
      ddF = -ddF;
      }*/
    nu += logcandderivs[1]*dF;

   // *workingweight = (*worktransformlin[0])*(*worktransformlin[0])*hilfs1-logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
   *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;

    if (*workingweight <=0)
      *workingweight = 0.001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  log(a) +(a*(*worktransformlin[0]))*log((*response)) - a*(*worktransformlin[0])*log((*worktransformlin[1]))
            - ((*worktransformlin[0])+1)*log(1+hilfs);

      }


  modify_worklin();

  }


void DISTR_dagum_a::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  //weight, response and param has size>1 if copula model is specified!!
  double exp_lin_a = exp((*linpred[copulaoffset+predstart_mumult+2]));
  double exp_lin_b = exp((*linpred[copulaoffset+predstart_mumult+1]));
  double exp_lin_p = exp((*linpred[copulaoffset+predstart_mumult]));
  double help1 = -1/exp_lin_a;
  double help2 = -help1 + exp_lin_p;

  *mu = 0;

  if (exp_lin_a>1)
  {
      *mu = -(exp_lin_b/exp_lin_a)*(randnumbers::gamma_exact(help1)*randnumbers::gamma_exact(help2))/(randnumbers::gamma_exact(exp_lin_p));
  }

  }


void DISTR_dagum_a::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (a): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_dagum_a::update_end(void)
  {


  // helpmat1 stores exp(eta_a)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_weibull_alpha -------------------------
//------------------------------------------------------------------------------


DISTR_weibull_alpha::DISTR_weibull_alpha(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Weibull Distribution - alpha";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "alpha";
    linpredminlimit=-7;
  linpredmaxlimit=7;
  }


DISTR_weibull_alpha::DISTR_weibull_alpha(const DISTR_weibull_alpha & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_weibull_alpha & DISTR_weibull_alpha::operator=(
                            const DISTR_weibull_alpha & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_weibull_alpha::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
      if (linpred_current==1)
        linpredp = linearpred1.getV();
      else
        linpredp = linearpred2.getV();
    }
// compute cdf (might work more efficiently)
  double res,lambda,alpha;
  alpha = exp(*linpredp);
  lambda = *worktransformlin[0];
  res = 1 - exp(-pow(resp*lambda,alpha));

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_weibull_alpha::cdf(const double & resp, const double & linpred)
  {
  double res,lambda,alpha;
  alpha = exp(linpred);
  lambda = *worktransformlin[0];//exp(*linpred[0]);
  res = 1 - exp(-pow(resp*lambda,alpha));
  return res;
  }

double DISTR_weibull_alpha::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_weibull_alpha::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = exp(*linpred[copulaoffset+0]);
  }

double DISTR_weibull_alpha::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double alpha = exp((*linpred));

  double l;

  l = (alpha-1)*log(*response) - pow((*response)*(*worktransformlin[0]),alpha) +alpha*log((*worktransformlin[0])) + log(alpha);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();

  return l;

  }

void DISTR_weibull_alpha::compute_iwls_wweightschange_weightsone(
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
    }

    double alpha = exp((*linpred));
    double hilfs1 = pow((*response)*(*worktransformlin[0]),alpha);

    double nu = 1 + alpha*log((*response)*(*worktransformlin[0]))*(1-hilfs1);

  //  *workingweight = 1 + pow(alpha,2)*pow((log((*response)/(*worktransformlin[0]))),2)*hilfs1;
    *workingweight =  1.823681;

    if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
    double dF = alpha*hilfs1*exp(-hilfs1)*log((*worktransformlin[0])*(*response));
    double ddF = -dF*dF/exp(-hilfs1)+dF
                +hilfs1*log((*worktransformlin[0])*(*response))*log((*worktransformlin[0])*(*response))*exp(-hilfs1)*alpha*alpha;
  /*  if(copularotate)
      {
      dF = -dF;
      ddF = -ddF;
      }*/
    nu += logcandderivs[1]*dF;

  /*  *workingweight = -alpha*log((*response)*(*worktransformlin[0])) * (1-hilfs1-hilfs1*alpha*log((*response)*(*worktransformlin[0])))
                        -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;*/
    *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;

    if (*workingweight <=0)
      *workingweight = 0.0001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {
        like +=  (alpha-1)*log(*response) - hilfs1 +alpha*log((*worktransformlin[0])) + log(alpha);

      }

  modify_worklin();

  }


void DISTR_weibull_alpha::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (alpha): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_weibull_alpha::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_weibull_lambda ----------------------
//------------------------------------------------------------------------------

void DISTR_weibull_lambda::check_errors(void)
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


DISTR_weibull_lambda::DISTR_weibull_lambda(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Weibull Distribution - lambda";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "lambda";
  linpredminlimit=-7;
  linpredmaxlimit=7;
  check_errors();
  }


DISTR_weibull_lambda::DISTR_weibull_lambda(const DISTR_weibull_lambda & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_weibull_lambda & DISTR_weibull_lambda::operator=(
                            const DISTR_weibull_lambda & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_weibull_lambda::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

 //  double test = *linpred;
// compute cdf (might work more efficiently)
  double res,lambda,alpha;
  lambda = exp(*linpredp);
  alpha = *worktransformlin[0];//exp(*linpred[0]);
  res = 1 - exp(-pow(resp*lambda,alpha));

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_weibull_lambda::cdf(const double & resp, const double & linpred)
  {
  double res,lambda,alpha;
  lambda = exp(linpred);
  alpha = *worktransformlin[0];//exp(*linpred[0]);
  res = 1 - exp(-pow(resp*lambda,alpha));
  return res;
  }

double DISTR_weibull_lambda::cdf(const double & resp, vector<double *>  linpred)
  {
  //linpred has only size=2!
  double res,lambda,alpha;
  lambda = exp(*linpred[1]);
  alpha = exp(*linpred[0]);//exp(*linpred[0]);
  res = 1 - exp(-pow(resp*lambda,alpha));
  return res;
  }


/*double DISTR_weibull_lambda::logpdf(const double & resp)
  {
  if(counter==0)
    {
    set_worklin();

    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

 //  double test = *linpred;
// compute cdf (might work more efficiently)
  double res,lambda,alpha;
  lambda = exp(*linpredp);
  alpha = *worktransformlin[0];//exp(*linpred[0]);
  res = (alpha-1)*log(resp) - pow(resp*lambda,alpha) +alpha*log(lambda) + log(alpha);

  modify_worklin();
  linpredp++;
  return res;
  }*/

void DISTR_weibull_lambda::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma
   // *linpred[1] = eta_mu

  // weight, linpred, response has only size=2!
   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double alpha = exp(*linpred[0]);
     double lambda = exp(*linpred[1]);

     double l;

      l =   (alpha-1)*log(*response[1]) - pow((*response[1])*lambda,alpha) +alpha*log(lambda) + log(alpha) ;


    *deviance = -2*l;
    }

  }


double DISTR_weibull_lambda::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_weibull_lambda::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = exp(*linpred[copulaoffset+1]);
  }

 double DISTR_weibull_lambda::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_weibull_lambda::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    //weight, response and param has size>1 if copula model is specified!!
    return ( 1 - exp(-pow((*response[copulaoffset+1])*(*param[copulaoffset+1]),(*param[copulaoffset+0]))) );
    }

double DISTR_weibull_lambda::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  if (counter==0)
    {
    set_worklin();
    }

  double lambda = exp(*linpred);

  double l;
  l = - pow((*response)*lambda,(*worktransformlin[0])) + (*worktransformlin[0])*log(lambda) ;

  if(optionsp->copula)
    {
    //implement loglik for copula models, i.e. add part logc
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();


  return l;

  }


void DISTR_weibull_lambda::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  if (counter==0)
    {
    set_worklin();
    }
  double lambda = exp(*linpred);

  double hilfs1 = pow((*response)*lambda,(*worktransformlin[0]));
 //   double hilfs2 = (*worktransformlin[0])+1;

  double nu = (*worktransformlin[0])*(1-hilfs1);

  //  *workingweight = (*worktransformlin[0])*randnumbers::gamma_exact(hilfs2)*((*worktransformlin[0])-1)+(*worktransformlin[0]);

  *workingweight = pow((*worktransformlin[0]),2);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
    double dF = hilfs1*exp(-hilfs1)*(*worktransformlin[0]);
    double ddF = dF*(*worktransformlin[0])-hilfs1*hilfs1*exp(-hilfs1)*(*worktransformlin[0])*(*worktransformlin[0]);
  /*  if(copularotate)
      {
      dF = -dF;
      ddF = -ddF;
      }*/
    nu += logcandderivs[1]*dF;

   // *workingweight = (*worktransformlin[0])*(*worktransformlin[0])*hilfs1-logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
   *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
    if (*workingweight <=0)
      *workingweight = 0.001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {
      like += - hilfs1 + (*worktransformlin[0])*log(lambda);
      }

    modify_worklin();

  }


void DISTR_weibull_lambda::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  //weight, response and param has size>1 if copula model is specified!!
   double hilfs = 1+1/exp((*linpred[copulaoffset+predstart_mumult]));
  *mu = (1/exp((*linpred[copulaoffset+predstart_mumult+1])))*randnumbers::gamma_exact(hilfs);

  }


void DISTR_weibull_lambda::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (lambda): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_weibull_lambda::update_end(void)
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
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_weibull2_alpha -------------------------
//------------------------------------------------------------------------------


DISTR_weibull2_alpha::DISTR_weibull2_alpha(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Weibull2 Distribution - alpha";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "alpha";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_weibull2_alpha::DISTR_weibull2_alpha(const DISTR_weibull2_alpha & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_weibull2_alpha & DISTR_weibull2_alpha::operator=(
                            const DISTR_weibull2_alpha & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_weibull2_alpha::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
      if (linpred_current==1)
        linpredp = linearpred1.getV();
      else
        linpredp = linearpred2.getV();
    }
 //  double test = *linpred;
// compute cdf (might work more efficiently)
  double res,lambda,alpha;
  alpha = exp(*linpredp);
  lambda = *worktransformlin[0];//exp(*linpred[0]);
  res = 1 - exp(-pow(resp/lambda,alpha));

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_weibull2_alpha::cdf(const double & resp, const double & linpred)
  {
  double res,lambda,alpha;
  alpha = exp(linpred);
  lambda = *worktransformlin[0];//exp(*linpred[0]);
  res = 1 - exp(-pow(resp/lambda,alpha));
  return res;
  }

double DISTR_weibull2_alpha::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_weibull2_alpha::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = exp(*linpred[copulaoffset+0]);
  }

double DISTR_weibull2_alpha::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  if (counter==0)
    {
    set_worklin();
    }

  double alpha = exp((*linpred));

  double l;

  l = (alpha-1)*log(*response) - pow((*response)/(*worktransformlin[0]),alpha) - alpha*log((*worktransformlin[0])) + log(alpha);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();

  return l;

  }

void DISTR_weibull2_alpha::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  if (counter==0)
    {
    set_worklin();
    }

    double alpha = exp((*linpred));
    double hilfs1 = pow((*response)/(*worktransformlin[0]),alpha);

    double nu = 1 + alpha*log((*response)/(*worktransformlin[0]))*(1-hilfs1);

  //  *workingweight = 1 + pow(alpha,2)*pow((log((*response)/(*worktransformlin[0]))),2)*hilfs1;
    *workingweight =  1.823681;

    if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
    double dF = alpha*hilfs1*exp(-hilfs1)*log((*response)/(*worktransformlin[0]));
    double ddF = -dF*dF/exp(-hilfs1)+dF
                +hilfs1*log((*response)/(*worktransformlin[0]))*log((*response)/(*worktransformlin[0]))*exp(-hilfs1)*alpha*alpha;
  /*  if(copularotate)
      {
      dF = -dF;
      ddF = -ddF;
      }*/
    nu += logcandderivs[1]*dF;

  /*  *workingweight = -alpha*log((*response)*(*worktransformlin[0])) * (1-hilfs1-hilfs1*alpha*log((*response)*(*worktransformlin[0])))
                        -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;*/
    *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;

    if (*workingweight <=0)
      *workingweight = 0.0001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {
        like +=  (alpha-1)*log(*response) - hilfs1 -alpha*log((*worktransformlin[0])) + log(alpha);

      }

  modify_worklin();

  }


void DISTR_weibull2_alpha::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (alpha): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_weibull2_alpha::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_weibull2_lambda ----------------------
//------------------------------------------------------------------------------

void DISTR_weibull2_lambda::check_errors(void)
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


DISTR_weibull2_lambda::DISTR_weibull2_lambda(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Weibull2 Distribution - lambda";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "lambda";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_weibull2_lambda::DISTR_weibull2_lambda(const DISTR_weibull2_lambda & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_weibull2_lambda & DISTR_weibull2_lambda::operator=(
                            const DISTR_weibull2_lambda & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_weibull2_lambda::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

 //  double test = *linpred;
// compute cdf (might work more efficiently)
  double res,lambda,alpha;
  lambda = exp(*linpredp);
  alpha = *worktransformlin[0];//exp(*linpred[0]);
  res = 1 - exp(-pow(resp/lambda,alpha));

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_weibull2_lambda::cdf(const double & resp, const double & linpred)
  {
  double res,lambda,alpha;
  lambda = exp(linpred);
  alpha = *worktransformlin[0];//exp(*linpred[0]);
  res = 1 - exp(-pow(resp/lambda,alpha));
  return res;
  }

double DISTR_weibull2_lambda::cdf(const double & resp, vector<double *>  linpred)
  {
  //linpred has always size=2!
  double res,lambda,alpha;
  lambda = exp(*linpred[1]);
  alpha = exp(*linpred[0]);//exp(*linpred[0]);
  res = 1 - exp(-pow(resp/lambda,alpha));
  return res;
  }


void DISTR_weibull2_lambda::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma
   // *linpred[1] = eta_mu

    //linpred, reponse, weight has always size=2!
   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double alpha = exp(*linpred[0]);
     double lambda = exp(*linpred[1]);

     double l;

      l =   (alpha-1)*log(*response[1]) - pow((*response[1])/lambda,alpha) -alpha*log(lambda) + log(alpha) ;


    *deviance = -2*l;
    }

  }


double DISTR_weibull2_lambda::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_weibull2_lambda::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]);
  }

 double DISTR_weibull2_lambda::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_weibull2_lambda::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    //weight, response and param has size>1 if copula model is specified!!
    return ( 1 - exp(-pow((*response[copulaoffset+1])/(*param[copulaoffset+1]),(*param[copulaoffset+0]))) );
    }

double DISTR_weibull2_lambda::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
 // *worklin[0] = linear predictor of alpha equation
  // *worktransformlin[0] = alpha;

  if (counter==0)
    {
    set_worklin();
    }

  double lambda = exp(*linpred);

  double l;
  l = - pow((*response)/lambda,(*worktransformlin[0])) - (*worktransformlin[0])*log(lambda) ;

  if(optionsp->copula)
    {
    //implement loglik for copula models, i.e. add part logc
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();


  return l;

  }


void DISTR_weibull2_lambda::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

 // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }
  double lambda = exp(*linpred);

  double hilfs1 = pow((*response)/lambda,(*worktransformlin[0]));
 //   double hilfs2 = (*worktransformlin[0])+1;

  double nu = -(*worktransformlin[0])*(1-hilfs1);

  //  *workingweight = (*worktransformlin[0])*randnumbers::gamma_exact(hilfs2)*((*worktransformlin[0])-1)+(*worktransformlin[0]);

  *workingweight = pow((*worktransformlin[0]),2);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
    double dF = -hilfs1*exp(-hilfs1)*(*worktransformlin[0]);
    double ddF = dF*(*worktransformlin[0])-hilfs1*hilfs1*exp(-hilfs1)*(*worktransformlin[0])*(*worktransformlin[0]);
  /*  if(copularotate)
      {
      dF = -dF;
      ddF = -ddF;
      }*/
    nu += logcandderivs[1]*dF;

   // *workingweight = (*worktransformlin[0])*(*worktransformlin[0])*hilfs1-logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
   *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
    if (*workingweight <=0)
      *workingweight = 0.001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {
      like += - hilfs1 - (*worktransformlin[0])*log(lambda);
      }

    modify_worklin();

  }


void DISTR_weibull2_lambda::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  //weight, response and param has size>1 if copula model is specified!!
   double hilfs = 1+1/exp((*linpred[copulaoffset+predstart_mumult]));
  *mu = (1/exp((*linpred[copulaoffset+predstart_mumult+1])))*randnumbers::gamma_exact(hilfs);

  }


void DISTR_weibull2_lambda::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (lambda): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_weibull2_lambda::update_end(void)
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
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gumbel_sigma -------------------------
//------------------------------------------------------------------------------


DISTR_gumbel_sigma::DISTR_gumbel_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "gumbel Distribution - sigma";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
  linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_gumbel_sigma::DISTR_gumbel_sigma(const DISTR_gumbel_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gumbel_sigma & DISTR_gumbel_sigma::operator=(
                            const DISTR_gumbel_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_gumbel_sigma::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
      if (linpred_current==1)
        linpredp = linearpred1.getV();
      else
        linpredp = linearpred2.getV();
    }
// compute cdf (might work more efficiently)
  double res,mu,sigma;
  sigma = exp(*linpredp);
  mu = *worktransformlin[0];
//  res = 1 - exp(-exp(-(resp-mu)/sigma));
  res = exp(-exp(-(resp-mu)/sigma));

/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF < 0.001" << endl;
    }*/

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_gumbel_sigma::cdf(const double & resp, const double & linpred)
  {
  double res,mu,sigma;
  sigma = exp(linpred);
  mu = *worktransformlin[0];
//  res = 1 - exp(-exp(-(resp-mu)/sigma));
  res = exp(-exp(-(resp-mu)/sigma));

/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF < 0.001" << endl;
    }*/

  return res;
  }

double DISTR_gumbel_sigma::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_gumbel_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = exp(*linpred[copulaoffset+1]);
  }

double DISTR_gumbel_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }

  double sigma = exp((*linpred));
  double hilfs = (*response-(*worktransformlin[0]))/sigma;

  double l;

//  l = log(sigma)-hilfs-exp(-hilfs);
  l = -log(sigma)-hilfs-exp(-hilfs);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();

  return l;
  }

void DISTR_gumbel_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma = exp((*linpred));
  double hilfs = (*response-(*worktransformlin[0]))/sigma;

  double nu = -1+hilfs-exp(-hilfs)*hilfs;

//  *workingweight =  hilfs*(1-exp(-hilfs))+hilfs*hilfs*exp(-hilfs);
  *workingweight =  hilfs - hilfs*exp(-hilfs) + hilfs*hilfs*exp(-hilfs);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
//    double dF = exp(-exp(-hilfs))*exp(-hilfs)*hilfs;
    double dF = -exp(-exp(-hilfs))*exp(-hilfs)*hilfs;
    double ddF = -dF*(1+exp(-hilfs)*hilfs-hilfs);
//    double ddF = -dF*(1+exp(-hilfs)*hilfs+hilfs);

    nu += logcandderivs[1]*dF;

    *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;

    if (*workingweight <=0)
      *workingweight = 0.0001;
    }

/*  if(isnan(*workingweight) || *workingweight > 100)
    {
    *workingweight = 1.0;
    cout << "counter: " << counter << endl;
    cout << "Gumbel sigma *workingweight NAN" << endl;
    }*/

  *workingresponse = *linpred + nu/(*workingweight);

/*  if(isnan(*workingresponse) || *workingresponse > 100 || *workingresponse < -100)
    {
    *workingresponse = 0.0;
    cout << "counter: " << counter << endl;
    cout << "Gumbel sigma *workingresponse NAN" << endl;
    }*/

  if (compute_like)
    {
//    like +=  log(sigma)-hilfs-exp(-hilfs);
    like += -log(sigma)-hilfs-exp(-hilfs);
    }

  modify_worklin();

  }


void DISTR_gumbel_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gumbel_sigma::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_gumbel_mu ---------------------------
//------------------------------------------------------------------------------

void DISTR_gumbel_mu::check_errors(void)
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

/*        if (*workresp < 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative response values encountered\n");
          }*/


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


DISTR_gumbel_mu::DISTR_gumbel_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "gumbel Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  check_errors();
  }


DISTR_gumbel_mu::DISTR_gumbel_mu(const DISTR_gumbel_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gumbel_mu & DISTR_gumbel_mu::operator=(
                            const DISTR_gumbel_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_gumbel_mu::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

 // compute cdf (might work more efficiently)
  double res,mu,sigma;
  mu = (*linpredp);
  sigma = *worktransformlin[0];
//  res = 1 - exp(-exp(-(resp-mu)/sigma));
  res = exp(-exp(-(resp-mu)/sigma));

/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF < 0.001" << endl;
    }*/

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_gumbel_mu::cdf(const double & resp, const double & linpred)
  {
  double res,mu,sigma;
  mu = (linpred);
  sigma = *worktransformlin[0];
//  res = 1 - exp(-exp(-(resp-mu)/sigma));
  res = exp(-exp(-(resp-mu)/sigma));

/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF < 0.001" << endl;
    }*/

  return res;
  }

double DISTR_gumbel_mu::cdf(const double & resp, vector<double *>  linpred)
  {
  //linpred has only size=2!
  double res,mu,sigma;
  mu = (*linpred[0]);
  sigma = exp(*linpred[1]);
//  res = 1 - exp(-exp(-(resp-mu)/sigma));
  res = exp(-exp(-(resp-mu)/sigma));

/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF < 0.001" << endl;
    }*/

  return res;
  }

void DISTR_gumbel_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma
   // *linpred[1] = eta_mu

  // weight, linpred, response has only size=2!
   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double sigma = exp(*linpred[1]);
     double mu = (*linpred[0]);

     double l;

//     l =   -(*response[1]-mu)/sigma - exp(-(*response[1]-mu)/sigma);
     l =  -log(sigma) -(*response[1]-mu)/sigma - exp(-(*response[1]-mu)/sigma);

    *deviance = -2*l;
    }

  }


double DISTR_gumbel_mu::get_intercept_start(void)
  {
  return 0.0;
  }

void DISTR_gumbel_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = (*linpred[copulaoffset+0]);
  }

 double DISTR_gumbel_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_gumbel_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    //weight, response and param has size>1 if copula model is specified!!
    return 0;
    }

double DISTR_gumbel_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;
//  l = -(*response-mu)/(*worktransformlin[0])-exp(-(*response-mu)/(*worktransformlin[0]));
  l = -log(*worktransformlin[0]) - (*response-mu)/(*worktransformlin[0]) - exp(-(*response-mu)/(*worktransformlin[0]));

  if(optionsp->copula)
    {
    //implement loglik for copula models, i.e. add part logc
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();
  return l;

  }


void DISTR_gumbel_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  if (counter==0)
    {
    set_worklin();
    }
  double mu = (*linpred);

  double hilfs =(*response-mu)/(*worktransformlin[0]);

  double nu = 1/(*worktransformlin[0])-exp(-hilfs)/(*worktransformlin[0]);

//  *workingweight = exp(-hilfs)/pow((*worktransformlin[0]),2);
  *workingweight = exp(-hilfs)/((*worktransformlin[0])*(*worktransformlin[0]));

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
//    double dF = exp(-exp(-hilfs))*exp(-hilfs)/(*worktransformlin[0]);
    double dF = -exp(-exp(-hilfs))*exp(-hilfs)/(*worktransformlin[0]);
    double ddF = dF*(1-exp(-hilfs))/(*worktransformlin[0]);
    nu += logcandderivs[1]*dF;

   *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
    if (*workingweight <=0)
      *workingweight = 0.001;
    }

/*  if(isnan(*workingweight) || *workingweight > 100)
    {
    *workingweight = 1.0;
    cout << "counter: " << counter << endl;
    cout << "Gumbel mu *workingweight NAN" << endl;
    }*/

  *workingresponse = *linpred + nu/(*workingweight);

/*  if(isnan(*workingresponse) || *workingresponse > 100 || *workingresponse < -100)
    {
    *workingresponse = 0.0;
    cout << "counter: " << counter << endl;
    cout << "Gumbel mu *workingresponse NAN" << endl;
    }*/

  if (compute_like)
    {
//      like += -(*response-mu)/(*workingweight)-exp(-(*response-mu)/(*workingweight));
    like += -log(*worktransformlin[0]) - hilfs - exp(-hilfs);
    }

    modify_worklin();
  }


void DISTR_gumbel_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = (*linpred[copulaoffset+predstart_mumult+1])+exp(*linpred[copulaoffset+predstart_mumult+2])* 0.5772156649;
  }


void DISTR_gumbel_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gumbel_mu::update_end(void)
  {

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu =(*worklin);
    }

  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gumbel2_sigma2 -------------------------
//------------------------------------------------------------------------------


DISTR_gumbel2_sigma2::DISTR_gumbel2_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "gumbel2 Distribution - sigma2";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
  linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_gumbel2_sigma2::DISTR_gumbel2_sigma2(const DISTR_gumbel2_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gumbel2_sigma2 & DISTR_gumbel2_sigma2::operator=(
                            const DISTR_gumbel2_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_gumbel2_sigma2::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
      if (linpred_current==1)
        linpredp = linearpred1.getV();
      else
        linpredp = linearpred2.getV();
    }
// compute cdf (might work more efficiently)
  double res,mu,sigma2;
  sigma2 = exp((*linpredp));
  mu = *worktransformlin[0];
//  res = 1 - exp(-exp(-(resp-mu)/sigma));
//  double hilfs=exp(log(resp-mu)-0.5*log(sigma2));
  double hilfs = (resp-mu)/sqrt(sigma2);
  res = 1-exp(-exp(hilfs));

/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF < 0.001" << endl;
    }*/

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_gumbel2_sigma2::cdf(const double & resp, const double & linpred)
  {
  double res,mu,sigma2;
  sigma2 = exp(linpred);
  mu = *worktransformlin[0];
//  res = 1 - exp(-exp(-(resp-mu)/sigma));
//  double hilfs=exp(log(resp-mu)-0.5*log(sigma2));
  double hilfs = (resp-mu)/sqrt(sigma2);
    res = 1-exp(-exp(hilfs));
/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF < 0.001" << endl;
    }*/

  return res;
  }

double DISTR_gumbel2_sigma2::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_gumbel2_sigma2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = exp(*linpred[copulaoffset+0]);
  }

double DISTR_gumbel2_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }

  double sigma = exp(0.5*(*linpred));
//  double sigma2 = exp((*linpred));
  double hilfs = (*response-(*worktransformlin[0]))/sigma;

  double l;

//  l = log(sigma)-hilfs-exp(-hilfs);
  l = -0.5*(*linpred)+hilfs-exp(hilfs);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();

  return l;
  }

void DISTR_gumbel2_sigma2::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma = exp(0.5*(*linpred));
//  double sigma2 = exp(*linpred);
  double hilfs = (*response-(*worktransformlin[0]))/sigma;
  double ehilfs = exp(hilfs);
  double nu = -0.5 - 0.5*hilfs + 0.5*ehilfs*hilfs;

//  *workingweight =  hilfs*(1-exp(-hilfs))+hilfs*hilfs*exp(-hilfs);
  *workingweight =  -0.25*hilfs + 0.25*hilfs*ehilfs + 0.25*hilfs*hilfs*ehilfs;

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
//    double dF = exp(-exp(-hilfs))*exp(-hilfs)*hilfs;
    double dF = -exp(-ehilfs+hilfs)*hilfs/2;
    double ddF = -dF*(1-ehilfs*hilfs+hilfs)/2;
//    double ddF = -dF*(1+exp(-hilfs)*hilfs+hilfs);

    nu += logcandderivs[1]*dF;

    *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;

    if (*workingweight <=0)
      *workingweight = 0.0001;
    }

/*  if(isnan(*workingweight) || *workingweight > 100)
    {
    *workingweight = 1.0;
   d cout << "counter: " << counter << endl;
    cout << "Gumbel sigma *workingweight NAN" << endl;
    }*/

  *workingresponse = *linpred + nu/(*workingweight);


  // cout << "workingweight sigma: " << *workingweight << endl;

  if (compute_like)
    {
//    like +=  log(sigma)-hilfs-exp(-hilfs);
    like += -0.5*(*linpred)+hilfs-ehilfs;
    }

  modify_worklin();

  }


void DISTR_gumbel2_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gumbel2_sigma2::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp((*worklin));
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_gumbel2_mu ---------------------------
//------------------------------------------------------------------------------

void DISTR_gumbel2_mu::check_errors(void)
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

/*        if (*workresp < 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative response values encountered\n");
          }*/


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


DISTR_gumbel2_mu::DISTR_gumbel2_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "gumbel2 Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  check_errors();
  }


DISTR_gumbel2_mu::DISTR_gumbel2_mu(const DISTR_gumbel2_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gumbel2_mu & DISTR_gumbel2_mu::operator=(
                            const DISTR_gumbel2_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_gumbel2_mu::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

 // compute cdf (might work more efficiently)
  double res,mu,sigma;
  mu = (*linpredp);
//  sigma2 = *worktransformlin[0];
  sigma=exp(0.5*log(*worktransformlin[0]));
//  double hilfs=exp(log(resp-mu)-0.5*log(sigma2));
  double hilfs = (resp-mu)/sigma;
  res = 1-exp(-exp(hilfs));

/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF < 0.001" << endl;
    }*/

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_gumbel2_mu::cdf(const double & resp, const double & linpred)
  {
  double res,mu,sigma;
  mu = (linpred);
//  sigma2 = *worktransformlin[0];
  sigma=exp(0.5*log(*worktransformlin[0]));
//  double hilfs=exp(log(resp-mu)-0.5*log(sigma2));
  double hilfs = (resp-mu)/sigma;
  res = 1-exp(-exp(hilfs));


/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF < 0.001" << endl;
    }*/

  return res;
  }

double DISTR_gumbel2_mu::cdf(const double & resp, vector<double *>  linpred)
  {
  //linpred has only size=2!
  double res,mu,sigma;
  mu = (*linpred[1]);
  sigma = exp(0.5*(*linpred[0]));
//  double hilfs=exp(log(resp-mu)-0.5*log(sigma2));
  double hilfs = (resp-mu)/(sigma);
  res = 1-exp(-exp(hilfs));
//  res = 1 - exp(-exp(-(resp-mu)/sigma));


/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Gumbel CDF < 0.001" << endl;
    }*/

  return res;
  }

void DISTR_gumbel2_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma
   // *linpred[1] = eta_mu

  // weight, linpred, response has only size=2!
   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double sigma = exp(0.5*(*linpred[0]));
     double mu = (*linpred[1]);
     double hilfs =(*response[1]-mu)/sigma;
     double l;

//     l =   -(*response[1]-mu)/sigma - exp(-(*response[1]-mu)/sigma);
     l = hilfs - exp(hilfs);

    *deviance = -2*l;
    }

  }


double DISTR_gumbel2_mu::get_intercept_start(void)
  {
  return 0.0;
  }

void DISTR_gumbel2_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = (*linpred[copulaoffset+1]);
  }

 double DISTR_gumbel2_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_gumbel2_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    //weight, response and param has size>1 if copula model is specified!!
    return 0;
    }

double DISTR_gumbel2_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;
//  l = -(*response-mu)/(*worktransformlin[0])-exp(-(*response-mu)/(*worktransformlin[0]));
//  double hilfs=exp(log(*response-mu)-0.5*log(*worktransformlin[0]));
  double hilfs = (*response-mu)/sqrt(*worktransformlin[0]);

  l =  hilfs - exp(hilfs);

  if(optionsp->copula)
    {
    //implement loglik for copula models, i.e. add part logc
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();
  return l;

  }


void DISTR_gumbel2_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  if (counter==0)
    {
    set_worklin();
    }
  double mu = (*linpred);
  double sigma=sqrt(*worktransformlin[0]);

//  double hilfs =exp(log(*response-mu)-0.5*log(*worktransformlin[0]));
  double hilfs = (*response-mu)/sigma;
  double ehilfs = exp(hilfs);
  double nu = -1/sigma+ehilfs/sigma;
//  *workingweight = exp(-hilfs)/pow((*worktransformlin[0]),2);
  *workingweight = ehilfs/((*worktransformlin[0]));

 /* cout << "mu: " << mu << endl;
  cout << "hilfs: " << hilfs << endl;
  cout << "worktransformlin[0]: " << *worktransformlin[0] << endl;
  cout << "ehilfs: " << ehilfs << endl;
  cout << "sigma: " << sigma << endl;
  cout << "nu: " << nu << endl;*/

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
//    cout << "F mu: " << F << endl;
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
//    double dF = exp(-exp(-hilfs))*exp(-hilfs)/(*worktransformlin[0]);
//cout << "in logcandderivs mu: " << nu << endl;
//      cout << "logcandderivs[1] mu: " << logcandderivs[2] << endl;
//      cout << "logcandderivs[2] mu: " << logcandderivs[2] << endl;
//    double dF = -exp(-ehilfs)*ehilfs/sigma;
//    double ddF = -dF*(1-ehilfs)/sigma;
//    nu += logcandderivs[1]*dF;
//
//   *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
//    if (*workingweight <=0)
//      *workingweight = 0.001;
//
//    cout << "dF mu: " << dF << endl;
//    cout << "ddF mu: " << ddF << endl;
//    cout << "nu mu: " << nu << endl;
//    cout << "workingweight mu: " << *workingweight << endl;
    }


/*  if(isnan(*workingweight) || *workingweight > 100)
    {
    *workingweight = 1.0;
    cout << "counter: " << counter << endl;
    cout << "Gumbel mu *workingweight NAN" << endl;
    }*/

  *workingresponse = *linpred + nu/(*workingweight);

/*  if(isnan(*workingresponse) || *workingresponse > 100 || *workingresponse < -100)
    {
    *workingresponse = 0.0;
    cout << "counter: " << counter << endl;
    cout << "Gumbel mu *workingresponse NAN" << endl;
    }*/

  if (compute_like)
    {
//      like += -(*response-mu)/(*workingweight)-exp(-(*response-mu)/(*workingweight));
    like += hilfs - ehilfs;
    }
  modify_worklin();
  }


void DISTR_gumbel2_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = (*linpred[copulaoffset+predstart_mumult+1])+exp(*linpred[copulaoffset+predstart_mumult+0])* 0.5772156649;
  }


void DISTR_gumbel2_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gumbel2_mu::update_end(void)
  {

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu =(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gengamma_tau --------------------------
//------------------------------------------------------------------------------


DISTR_gengamma_tau::DISTR_gengamma_tau(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Generalized Gamma Distribution - tau";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "tau";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_gengamma_tau::DISTR_gengamma_tau(const DISTR_gengamma_tau & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gengamma_tau & DISTR_gengamma_tau::operator=(
                            const DISTR_gengamma_tau & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_gengamma_tau::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gengamma_tau::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[0]);
  }

double DISTR_gengamma_tau::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = exp(eta_sigma);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double tau = exp((*linpred));

  double l;

     l = log(tau) + (tau*(*worktransformlin[0]))*log((*response)) -pow((((*worktransformlin[0])/(*worktransformlin[1]))*(*response)),tau) +
     	 tau*(*worktransformlin[0])*log((*worktransformlin[0])/(*worktransformlin[1]));


  modify_worklin();

  return l;

  }

void DISTR_gengamma_tau::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = exp(eta_sigma);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);


  if (counter==0)
    {
    set_worklin();
    }

    double tau = exp((*linpred));


    double nu = 1 + (*worktransformlin[0])*tau*log((*response)) -
				pow(((*worktransformlin[0])/(*worktransformlin[1]))*(*response),tau)*tau*log(((*worktransformlin[0])/(*worktransformlin[1]))*(*response)) +
				tau*(*worktransformlin[0])*log((*worktransformlin[0])/(*worktransformlin[1]));

	double exp_linsigma_plus1 = ((*worktransformlin[0])+1);

    *workingweight = (*worktransformlin[0])*((randnumbers::trigamma_exact(exp_linsigma_plus1))+pow((randnumbers::digamma_exact(exp_linsigma_plus1)),2))+1;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += log(tau) + (tau*(*worktransformlin[0]))*log((*response)) -pow((((*worktransformlin[0])/(*worktransformlin[1]))*(*response)),tau) +
     	 tau*(*worktransformlin[0])*log((*worktransformlin[0])/(*worktransformlin[1]));

      }

  modify_worklin();

  }


void DISTR_gengamma_tau::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (tau): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gengamma_tau::update_end(void)
  {

  // helpmat1 stores tau

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gengamma_sigma ------------------------
//------------------------------------------------------------------------------


DISTR_gengamma_sigma::DISTR_gengamma_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Generalized Gamma Distribution - sigma";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_gengamma_sigma::DISTR_gengamma_sigma(const DISTR_gengamma_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gengamma_sigma & DISTR_gengamma_sigma::operator=(
                            const DISTR_gengamma_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_gengamma_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gengamma_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]);
  }

double DISTR_gengamma_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = exp(eta_tau);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sig = exp((*linpred));

  double l;

     l = (sig*(*worktransformlin[0])-1)*log((*response)) -pow(((sig/(*worktransformlin[1]))*(*response)),(*worktransformlin[0])) +
     	 sig*(*worktransformlin[0])*log(sig/(*worktransformlin[1])) -randnumbers::lngamma_exact(sig);


  modify_worklin();

  return l;

  }

void DISTR_gengamma_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = exp(eta_tau);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

    double sig = exp((*linpred));


    double nu = sig*(*worktransformlin[0])*log((*response)) - (*worktransformlin[0])*pow((sig/(*worktransformlin[1]))*(*response),(*worktransformlin[0])) +
				sig*(*worktransformlin[0])*log(sig/((*worktransformlin[1]))) + sig*(*worktransformlin[0]) - sig*(randnumbers::digamma_exact(sig));



    *workingweight = sig*(sig*randnumbers::trigamma_exact(sig) - 2*(*worktransformlin[0]) + pow((*worktransformlin[0]),2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  (sig*(*worktransformlin[0]))*log((*response)) -pow(((sig/(*worktransformlin[1]))*(*response)),(*worktransformlin[0])) +
     	 sig*(*worktransformlin[0])*log(sig/(*worktransformlin[1])) -randnumbers::lngamma_exact(sig);

      }

  modify_worklin();

  }


void DISTR_gengamma_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gengamma_sigma::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_gengamma_mu -------------------------
//------------------------------------------------------------------------------

void DISTR_gengamma_mu::check_errors(void)
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

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
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

DISTR_gengamma_mu::DISTR_gengamma_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Generalized Gamma Distribution - mu";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_gengamma_mu::DISTR_gengamma_mu(const DISTR_gengamma_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gengamma_mu & DISTR_gengamma_mu::operator=(
                            const DISTR_gengamma_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_gengamma_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_tau
   // *linpred[1] = eta_sigma
   // *linpred[2] = eta_mu

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
	 double tau = exp(*linpred[0]);
     double sig = exp(*linpred[1]);
     double mu = exp(*linpred[2]);

     double l;

       l = log(tau) + (sig*tau-1)*log((*response[2])) -pow(((sig/mu)*(*response[2])),tau) + sig*tau*log(sig/mu) -randnumbers::lngamma_exact(sig);


    *deviance = -2*l;
    }

  }


double DISTR_gengamma_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gengamma_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[2]);
  }
 double DISTR_gengamma_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_gengamma_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {

    return ( 0 );
    }


double DISTR_gengamma_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = tau;
  // *worklin[1] = linear predictor of sigma equation
  // *worktransformlin[1] = sigma;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);

  double l;

     l =  -pow((((*worktransformlin[1])/mu)*(*response)),(*worktransformlin[0])) - (*worktransformlin[1])*(*worktransformlin[0])*log(mu) ;
  modify_worklin();

  return l;

  }


void DISTR_gengamma_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of tau equation
  // *worktransformlin[0] = tau;
  // *worklin[1] = linear predictor of sigma equation
  // *worktransformlin[1] = sigma;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = exp((*linpred));

    double exponent = (*worktransformlin[0]);

    double nu = (*worktransformlin[0])*pow(((*worktransformlin[1])/mu)*(*response),exponent) - ((*worktransformlin[1]))*((*worktransformlin[0]));

    *workingweight = (*worktransformlin[1])*pow((*worktransformlin[0]),2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -pow((((*worktransformlin[1])/mu)*(*response)),(*worktransformlin[0])) - (*worktransformlin[1])*(*worktransformlin[0])*log(mu);

      }


  modify_worklin();

  }


void DISTR_gengamma_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

  double exp_linmu = exp((*linpred[predstart_mumult+2]));
  double exp_linsigma = exp((*linpred[predstart_mumult+1]));
  double exp_lintau = exp((*linpred[predstart_mumult]));
  double help = exp_linsigma+1/exp_lintau;
  *mu = (randnumbers::lngamma_exact(help))*exp_linmu/(exp_linsigma*(randnumbers::lngamma_exact(exp_linsigma)));

  }


void DISTR_gengamma_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gengamma_mu::update_end(void)
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
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gamma_sigma ---------------------------
//------------------------------------------------------------------------------


DISTR_gamma_sigma::DISTR_gamma_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Gamma Distribution - sigma";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
  linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_gamma_sigma::DISTR_gamma_sigma(const DISTR_gamma_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gamma_sigma & DISTR_gamma_sigma::operator=(
                            const DISTR_gamma_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_gamma_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

double DISTR_gamma_sigma::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
      if (linpred_current==1)
        linpredp = linearpred1.getV();
      else
        linpredp = linearpred2.getV();
    }
  double const sigma = exp(*linpredp);
  double const mu = *worktransformlin[0];
  double const res = randnumbers::gamma_cdf(resp, mu, sigma);;

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_gamma_sigma::cdf(const double & resp, const double & linpred)
  {
  double const sigma = exp(linpred);
  double const mu = *worktransformlin[0];
  double const res = randnumbers::gamma_cdf(resp, mu, sigma);;
  return res;
  }

void DISTR_gamma_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[copulaoffset+0]));
  }

double DISTR_gamma_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sig = exp((*linpred));

  double l;
  l = sig*log(sig) - sig*log((*worktransformlin[0])) - randnumbers::lngamma_exact(sig) + (sig-1)*log(*response) - (sig/(*worktransformlin[0]))*(*response);
  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();

  return l;
  }

void DISTR_gamma_sigma::compute_iwls_wweightschange_weightsone(
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
    }

  double sig = exp((*linpred));
  double nu = sig*log(sig) +  sig - sig*log((*worktransformlin[0])) - sig*(randnumbers::digamma_exact(sig)) +
		sig*log((*response)) - (sig/(*worktransformlin[0]))*(*response);

  *workingweight = sig*(sig*randnumbers::trigamma_exact(sig) - 1);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute dF/deta, d^2 F/deta ^2
    // TODO better implementation
    // numerical derivatives
    double eps = 1e-06;
    double dF = (cdf(*response, *linpred + 0.5*eps) - cdf(*response, *linpred - 0.5*eps))/eps;
    double ddF = (cdf(*response, *linpred + eps) + cdf(*response, *linpred - eps) - 2.0 * F)/(eps*eps);

    nu += logcandderivs[1]*dF;

    *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;

    if (*workingweight <=0)
      *workingweight = 0.0001;
    }

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    like += sig*log(sig) - sig*log((*worktransformlin[0])) - randnumbers::lngamma_exact(sig) + (sig-1)*log(*response) - (sig/(*worktransformlin[0]))*(*response);
    }
  modify_worklin();
  }


void DISTR_gamma_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gamma_sigma::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_gamma_mu ----------------------------
//------------------------------------------------------------------------------

void DISTR_gamma_mu::check_errors(void)
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

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
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

DISTR_gamma_mu::DISTR_gamma_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Gamma Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_gamma_mu::DISTR_gamma_mu(const DISTR_gamma_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_gamma_mu & DISTR_gamma_mu::operator=(
                            const DISTR_gamma_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_gamma_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma
   // *linpred[1] = eta_mu

   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double sig = exp(*linpred[0]);
     double mu = exp(*linpred[1]);

     double l;

       l = sig*log(sig) - sig*log(mu) - randnumbers::lngamma_exact(sig) + (sig-1)*log(*response[1]) - (sig/mu)*(*response[1]);


    *deviance = -2*l;
    }

  }


double DISTR_gamma_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

double DISTR_gamma_mu::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

  double const mu = exp(*linpredp);
  double const sigma = *worktransformlin[0];
  double const res = randnumbers::gamma_cdf(resp, mu, sigma);;

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_gamma_mu::cdf(const double & resp, const double & linpred)
  {
  double const mu = exp(*linpredp);
  double const sigma = *worktransformlin[0];
  double const res = randnumbers::gamma_cdf(resp, mu, sigma);;
  return res;
  }

double DISTR_gamma_mu::cdf(const double & resp, vector<double *>  linpred)
  {
  double const mu = exp(*linpred[0]);
  double const sigma = exp(*linpred[1]);
  double const res = randnumbers::gamma_cdf(resp, mu, sigma);;
  return res;
  }

void DISTR_gamma_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //this was for Alex paper
  //*param = exp((*linpred[1])) / exp((*linpred[0]));
  *param = (*linpred[copulaoffset+0]);
  }

 double DISTR_gamma_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_gamma_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {

    return ( 0 );
    }

double DISTR_gamma_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);

  double l;
  l =  - (*worktransformlin[0])*log(mu) - ((*worktransformlin[0])/mu)*(*response);
  if(optionsp->copula)
    {
    //implement loglik for copula models, i.e. add part logc
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();

  return l;

  }


void DISTR_gamma_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);
  double nu = -(*worktransformlin[0]) + ((*worktransformlin[0])/mu)*(*response);
  *workingweight = (*worktransformlin[0]);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute dF/deta, d^2 F/deta ^2
    // TODO better implementation
    // numerical derivatives
    double eps = 1e-06;
    double dF = (cdf(*response, *linpred + 0.5*eps) - cdf(*response, *linpred - 0.5*eps))/eps;
    double ddF = (cdf(*response, *linpred + eps) + cdf(*response, *linpred - eps) - 2.0 * F)/(eps*eps);
    nu += logcandderivs[1]*dF;

   *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
    if (*workingweight <=0)
      *workingweight = 0.001;
    }

  *workingresponse = *linpred + nu/(*workingweight);
  if (compute_like)
    {
    like +=  - (*worktransformlin[0])*log(mu) - ((*worktransformlin[0])/mu)*(*response);
    }
  modify_worklin();
  }


void DISTR_gamma_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

  *mu = exp((*linpred[copulaoffset+predstart_mumult+1]));

  }


void DISTR_gamma_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gamma_mu::update_end(void)
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
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_lognormal2_sigma ----------------------
//------------------------------------------------------------------------------


DISTR_lognormal2_sigma::DISTR_lognormal2_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Lognormal Distribution - sigma";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_lognormal2_sigma::DISTR_lognormal2_sigma(const DISTR_lognormal2_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_lognormal2_sigma & DISTR_lognormal2_sigma::operator=(
                            const DISTR_lognormal2_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_lognormal2_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_lognormal2_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[0]));
  }

double DISTR_lognormal2_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = pow(exp((*linpred)),2);

  double l;

     l = -0.5*log(sigma_2)-pow((log((*response))-(*worklin[0])),2)/(2*sigma_2);


  modify_worklin();

  return l;

  }

void DISTR_lognormal2_sigma::compute_iwls_wweightschange_weightsone(
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

    double sigma_2 = pow(exp(*linpred),2);


    double nu = -1 + (pow((log(*response)-(*worklin[0])),2))/(sigma_2);



    *workingweight = 2;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((log((*response))-(*worklin[0])),2)/(2*sigma_2);

      }

  modify_worklin();

  }


void DISTR_lognormal2_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_lognormal2_sigma::update_end(void)
  {

  // helpmat1 stores sigma2

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
//--------------------------- CLASS: DISTR_lognormal2_mu -----------------------
//------------------------------------------------------------------------------
void DISTR_lognormal2_mu::check_errors(void)
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

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
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


DISTR_lognormal2_mu::DISTR_lognormal2_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Lognormal Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  updateIWLS =false;
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
 check_errors();
  }


DISTR_lognormal2_mu::DISTR_lognormal2_mu(const DISTR_lognormal2_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_lognormal2_mu & DISTR_lognormal2_mu::operator=(
                            const DISTR_lognormal2_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_lognormal2_mu::compute_deviance_mult(vector<double *> response,
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

       l = -0.5*log(2*PI)-0.5*log(sigma_2)-log((*response[0]))-pow((log((*response[0]))-mu),2)/(2*sigma_2);


    *deviance = -2*l;
    }

  }


double DISTR_lognormal2_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_lognormal2_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[1]);
  }

 double DISTR_lognormal2_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_lognormal2_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = (log(*response[1])-(*param[1]))/(*param[0]) ;

    return (randnumbers::Phi2(arg));
    }

double DISTR_lognormal2_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;

     l = -pow((log((*response))-mu),2)/(2*(*worktransformlin[0]));

  modify_worklin();

  return l;

  }


void DISTR_lognormal2_mu::compute_iwls_wweightschange_weightsone(
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

    double mu = (*linpred);

    double nu = (log(*response)-mu)/(*worktransformlin[0]);

//    if(updateIWLS)
      *workingweight = 1/(*worktransformlin[0]);
//    else
//      *workingweight = (*weightp)/(*worktransformlin[0]);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((log((*response))-mu),2)/(2*(*worktransformlin[0]));

      }


  modify_worklin();

  }


void DISTR_lognormal2_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double s2 = pow(exp(*linpred[predstart_mumult]),2);
  *mu = exp((*linpred[predstart_mumult+1])+s2/2);
  }


void DISTR_lognormal2_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }

/*void DISTR_lognormal2_mu::update(void)
  {

   unsigned i;

  double help;

  double * worktransformlinp;
  double * workweight;

  worktransformlinp = distrp[0]->helpmat1.getV();
  workweight = workingweight.getV();

  if(updateIWLS)
    {
    for (i=0;i<nrobs;i++,worktransformlinp++,workweight++)
      {
      *workweight = 1/(*worktransformlinp);
      }
    }
  else
    {
    weightp = weight.getV();
    for (i=0;i<nrobs;i++,worktransformlinp++,workweight++,weightp++)
      {
      *workweight = (*weightp)/(*worktransformlinp);
      }
    }

  }*/

void DISTR_lognormal2_mu::update_end(void)
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

void DISTR_lognormal2_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  weightp = weight.getV();

  }

void DISTR_lognormal2_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    weightp++;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_lognormal_sigma2 ----------------------
//------------------------------------------------------------------------------


DISTR_lognormal_sigma2::DISTR_lognormal_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Lognormal Distribution - sigma2";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_lognormal_sigma2::DISTR_lognormal_sigma2(const DISTR_lognormal_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_lognormal_sigma2 & DISTR_lognormal_sigma2::operator=(
                            const DISTR_lognormal_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_lognormal_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_lognormal_sigma2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //standarddeviation
  //*param = pow(exp((*linpred[0])),0.5);
  *param = exp((*linpred[0]));
  }


double DISTR_lognormal_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = exp((*linpred));

  double l;

     l = -0.5*log(sigma_2)-pow((log((*response))-(*worklin[0])),2)/(2*sigma_2);


  modify_worklin();

  return l;

  }

void DISTR_lognormal_sigma2::compute_iwls_wweightschange_weightsone(
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

    double sigma_2 = exp(*linpred);


    double nu = -0.5 + (pow((log(*response)-(*worklin[0])),2))/(2*sigma_2);



    *workingweight = 0.5;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((log((*response))-(*worklin[0])),2)/(2*sigma_2);

      }

  modify_worklin();

  }


void DISTR_lognormal_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_lognormal_sigma2::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_lognormal_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_lognormal_mu::check_errors(void)
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

        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
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


DISTR_lognormal_mu::DISTR_lognormal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Lognormal Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  if(o->forceIWLS)
    {
    updateIWLS = true;
    }
  else
    {
    updateIWLS = false;
    }
  check_errors();
  }


DISTR_lognormal_mu::DISTR_lognormal_mu(const DISTR_lognormal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_lognormal_mu & DISTR_lognormal_mu::operator=(
                            const DISTR_lognormal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_lognormal_mu::compute_deviance_mult(vector<double *> response,
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
     double sigma_2 = exp(*linpred[0]);
     double mu = (*linpred[1]);

     double l;

       l = -0.5*log(2*PI)-0.5*log(sigma_2)-log((*response[0]))-pow((log((*response[0]))-mu),2)/(2*sigma_2);


    *deviance = -2*l;
    }

  }


double DISTR_lognormal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }



void DISTR_lognormal_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[1]);
  }

 double DISTR_lognormal_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_lognormal_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = (log(*response[1])-(*param[1]))/((*param[0])) ;

    return (randnumbers::Phi2(arg));
    }

double DISTR_lognormal_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;

     l = -pow((log((*response))-mu),2)/(2*(*worktransformlin[0]));

  modify_worklin();

  return l;

  }


void DISTR_lognormal_mu::compute_iwls_wweightschange_weightsone(
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

    double mu = (*linpred);

    double nu = (log(*response)-mu)/(*worktransformlin[0]);

//    if(updateIWLS)
      *workingweight = 1/(*worktransformlin[0]);
//    else
//      *workingweight = (*weightp)/(*worktransformlin[0]);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {
      like += -pow((log((*response))-mu),2)/(2*(*worktransformlin[0]));
      }


  modify_worklin();

  }


void DISTR_lognormal_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double s = exp(*linpred[predstart_mumult]);
  *mu = exp((*linpred[predstart_mumult+1])+s/2);
  }


void DISTR_lognormal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }

/*void DISTR_lognormal_mu::update(void)
  {

   unsigned i;

  double help;

  double * worktransformlinp;
  double * workweight;

  worktransformlinp = distrp[0]->helpmat1.getV();
  workweight = workingweight.getV();

  if(updateIWLS)
    {
    for (i=0;i<nrobs;i++,worktransformlinp++,workweight++)
    {
      *workweight = 1/(*worktransformlinp);
      }
    }
  else
    {
    weightp = weight.getV();
    for (i=0;i<nrobs;i++,worktransformlinp++,workweight++,weightp++)
      {
      *workweight = (*weightp)/(*worktransformlinp);
      }
    }
  }*/

void DISTR_lognormal_mu::update_end(void)
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

void DISTR_lognormal_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  weightp = weight.getV();

  }



void DISTR_lognormal_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    weightp++;
    }

  }

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_truncnormal2_sigma --------------------
//------------------------------------------------------------------------------


DISTR_truncnormal2_sigma::DISTR_truncnormal2_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Truncated normal2 Distribution - sigma";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_truncnormal2_sigma::DISTR_truncnormal2_sigma(const DISTR_truncnormal2_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_truncnormal2_sigma & DISTR_truncnormal2_sigma::operator=(
                            const DISTR_truncnormal2_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_truncnormal2_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_truncnormal2_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[0]));
  }

double DISTR_truncnormal2_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sig = exp((*linpred));

  double l;

     l = -log(sig)-pow((((*response))-(*worklin[0])),2)/(2*pow(sig,2)) - log(randnumbers::Phi2((*worklin[0])/sig));


  modify_worklin();

  return l;

  }

void DISTR_truncnormal2_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = mu

  if (counter==0)
    {
    set_worklin();
    }

    double sig = exp(*linpred);

    double dens = (1/(pow(2*PI,0.5)))*exp(-0.5*pow(((*worklin[0])/sig),2));

    double hilfs = ((*worklin[0])/sig)*(dens/randnumbers::Phi2((*worklin[0])/sig));

    double nu = -1 + (pow(((*response)-(*worklin[0])),2))/pow(sig,2) + hilfs;


   // *workingweight = 2;
    *workingweight = 2 - hilfs*( 1 + pow(((*worklin[0])/sig),2) + hilfs );

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -log(sig)-pow((((*response))-(*worklin[0])),2)/(2*pow(sig,2)) - log(randnumbers::Phi2((*worklin[0])/sig));

      }

  modify_worklin();

  }


void DISTR_truncnormal2_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_truncnormal2_sigma::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_truncnormal2_mu ---------------------
//------------------------------------------------------------------------------
void DISTR_truncnormal2_mu::check_errors(void)
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


DISTR_truncnormal2_mu::DISTR_truncnormal2_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Truncated normal2 Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  check_errors();
//    linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_truncnormal2_mu::DISTR_truncnormal2_mu(const DISTR_truncnormal2_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_truncnormal2_mu & DISTR_truncnormal2_mu::operator=(
                            const DISTR_truncnormal2_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_truncnormal2_mu::compute_deviance_mult(vector<double *> response,
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
     double sig = exp(*linpred[0]);
     double mu = (*linpred[1]);

     double l;

       l = -0.5*log(2*PI)-log(sig)-pow((((*response[0]))-mu),2)/(2*pow(sig,2)) - log(randnumbers::Phi2(mu/sig));


    *deviance = -2*l;
    }

  }


double DISTR_truncnormal2_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_truncnormal2_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[1]);
  }

 double DISTR_truncnormal2_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_truncnormal2_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = ((*response[1])-(*param[1]))/(*param[0]) ;
    double res = (randnumbers::Phi2(arg)-randnumbers::Phi2(-(*param[1])/(*param[0])))/(randnumbers::Phi2(-(*param[1])/(*param[0])));
    return (res);
    }

double DISTR_truncnormal2_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;

     l = -pow((((*response))-mu),2)/(2*pow((*worktransformlin[0]),2)) - log(randnumbers::Phi2(mu/(*worktransformlin[0])));

  modify_worklin();

  return l;

  }


void DISTR_truncnormal2_mu::compute_iwls_wweightschange_weightsone(
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

    double mu = (*linpred);

    double dens = (1/(pow(2*PI,0.5)))*exp(-0.5*pow((mu/(*worktransformlin[0])),2));

    double hilfs = (1/(*worktransformlin[0]))*(dens/randnumbers::Phi2(mu/(*worktransformlin[0])));

    double nu = ((*response)-mu)/pow((*worktransformlin[0]),2) - hilfs;

    *workingweight = 1/pow((*worktransformlin[0]),2) - (mu/pow((*worktransformlin[0]),2))*hilfs - pow(hilfs,2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((((*response))-mu),2)/(2*pow((*worktransformlin[0]),2)) - log(randnumbers::Phi2(mu/(*worktransformlin[0])));

      }


  modify_worklin();

  }


void DISTR_truncnormal2_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
   double arg = (*linpred[predstart_mumult+1])/exp(*linpred[predstart_mumult]);
   double hilfs = (exp(*linpred[predstart_mumult])*(1/(pow(2*PI,0.5)))*exp(-0.5*pow(arg,2)))/randnumbers::Phi2(arg);
  *mu = ((*linpred[predstart_mumult+1])) + hilfs;
  }


void DISTR_truncnormal2_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_truncnormal2_mu::update_end(void)
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


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_normal2_sigma -------------------------
//------------------------------------------------------------------------------


DISTR_normal2_sigma::DISTR_normal2_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Normal Distribution - sigma";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_normal2_sigma::DISTR_normal2_sigma(const DISTR_normal2_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_normal2_sigma & DISTR_normal2_sigma::operator=(
                            const DISTR_normal2_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_normal2_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_normal2_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[0]));
  }

double DISTR_normal2_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = pow(exp((*linpred)),2);

  double l;

     l = -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*sigma_2);


  modify_worklin();

  return l;

  }

void DISTR_normal2_sigma::compute_iwls_wweightschange_weightsone(
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

    double sigma_2 = pow(exp(*linpred),2);


    double nu = -1 + (pow(((*response)-(*worklin[0])),2))/(sigma_2);



    *workingweight = 2;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*sigma_2);

      }

  modify_worklin();

  }


void DISTR_normal2_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_normal2_sigma::update_end(void)
  {

  // helpmat1 stores sigma

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
//--------------------------- CLASS: DISTR_normal2_mu --------------------------
//------------------------------------------------------------------------------
void DISTR_normal2_mu::check_errors(void)
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


DISTR_normal2_mu::DISTR_normal2_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Normal Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  if(o->copula || o->forceIWLS)
    {
    updateIWLS = true;
    }
  else
    {
    updateIWLS = false;
    }
  check_errors();
//    linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_normal2_mu::DISTR_normal2_mu(const DISTR_normal2_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_normal2_mu & DISTR_normal2_mu::operator=(
                            const DISTR_normal2_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_normal2_mu::compute_deviance_mult(vector<double *> response,
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


double DISTR_normal2_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_normal2_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[1]);
  }

 double DISTR_normal2_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_normal2_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    double arg = ((*response[1])-(*param[1]))/(*param[0]) ;

    return (randnumbers::Phi2(arg));
    }

double DISTR_normal2_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;

     l = -pow((((*response))-mu),2)/(2*(*worktransformlin[0]));

  modify_worklin();

  return l;

  }


void DISTR_normal2_mu::compute_iwls_wweightschange_weightsone(
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

    double mu = (*linpred);

    double nu = ((*response)-mu)/(*worktransformlin[0]);

//    if(updateIWLS)
      *workingweight = 1/(*worktransformlin[0]);
//    else
//      *workingweight = (*weightp)/(*worktransformlin[0]);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((((*response))-mu),2)/(2*(*worktransformlin[0]));

      }


  modify_worklin();

  }


void DISTR_normal2_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+1]));
  }


void DISTR_normal2_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


/*void DISTR_normal2_mu::update(void)
  {

   unsigned i;

  double help;

  double * worktransformlinp;
  double * workweight;

  worktransformlinp = distrp[0]->helpmat1.getV();
  workweight = workingweight.getV();

  if(updateIWLS)
    {
    for (i=0;i<nrobs;i++,worktransformlinp++,workweight++)
    {
      *workweight = 1/(*worktransformlinp);
      }
    }
  else
    {
    weightp = weight.getV();
    for (i=0;i<nrobs;i++,worktransformlinp++,workweight++,weightp++)
      {
      *workweight = (*weightp)/(*worktransformlinp);
      }
    }

  }*/

void DISTR_normal2_mu::update_end(void)
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

void DISTR_normal2_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  weightp = weight.getV();

  }

void DISTR_normal2_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    weightp++;
    }

  }



//--------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_normal_sigma2 ---------------------------
//--------------------------------------------------------------------------------


DISTR_normal_sigma2::DISTR_normal_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Normal Distribution - sigma2";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_normal_sigma2::DISTR_normal_sigma2(const DISTR_normal_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_normal_sigma2 & DISTR_normal_sigma2::operator=(
                            const DISTR_normal_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_normal_sigma2::cdf(const double & resp, const double & linpred)
  {
  double res,mu,sigma2;
  sigma2 = exp(linpred);
  mu = *worklin[0];
  double z = (resp-mu)/sqrt(sigma2);
  res = randnumbers::Phi2(z);
  return res;
  }

double DISTR_normal_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_normal_sigma2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //standarddeviations
  //*param = pow(exp((*linpred[0])),0.5);
  //linpred has always size=2!
  *param = exp((*linpred[0]));
  }

double DISTR_normal_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = exp((*linpred));

  double l;

     l = -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*sigma_2);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();

  return l;

  }

void DISTR_normal_sigma2::compute_iwls_wweightschange_weightsone(
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

    double sigma_2 = exp(*linpred);


    double nu = -0.5 + (pow(((*response)-(*worklin[0])),2))/(2*sigma_2);



    *workingweight = 0.5;

    if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
    double z = ((*response)-(*worklin[0]))/sqrt(sigma_2);
    double dF = -0.5*z*exp(-0.5*z*z)*0.3989423;
    double ddF = -0.25*z*exp(-0.5*z*z)*(1-z*z)*0.3989423;
   /* if(copularotate)
      {
      dF = -dF;
      ddF = -ddF;
      }*/
    nu += logcandderivs[1]*dF;

  //  *workingweight = 0.5*z*z-logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
   *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
    if (*workingweight <=0)
      *workingweight = 0.0001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*sigma_2);

      }

  modify_worklin();

  }


void DISTR_normal_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_normal_sigma2::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_normal_mu ---------------------------
//------------------------------------------------------------------------------
void DISTR_normal_mu::check_errors(void)
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


DISTR_normal_mu::DISTR_normal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Normal Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  check_errors();
  if(o->copula || o->forceIWLS)
    {
    updateIWLS = true;
    }
  else
    {
    updateIWLS = true;
    }
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_normal_mu::DISTR_normal_mu(const DISTR_normal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_normal_mu & DISTR_normal_mu::operator=(
                            const DISTR_normal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

double DISTR_normal_mu::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

  double res,mu,sigma2;
  mu = (*linpredp);
  sigma2 = *worktransformlin[0];
  double z = (resp-mu)/sqrt(sigma2);
  res = randnumbers::Phi2(z);

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_normal_mu::cdf(const double & resp, const double & linpred)
  {
  double res,mu,sigma2;
  mu = linpred;
  sigma2 = *worktransformlin[0];
  double z = (resp-mu)/sqrt(sigma2);
  res = randnumbers::Phi2(z);
  return res;
  }

double DISTR_normal_mu::cdf(const double & resp, vector<double *>  linpred)
  {
    //linpred has always size=2
  double res,mu,sigma2;
  mu = (*linpred[1]);
  sigma2 = exp(*linpred[0]);
  double z = (resp-mu)/sqrt(sigma2);
  res = randnumbers::Phi2(z);
  return res;
  }


/*double DISTR_normal_mu::logpdf(const double & resp)
  {
  if(counter==0)
    {
    set_worklin();

    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

 //  double test = *linpred;
// compute cdf (might work more efficiently)
  double res,mu,sigma2;
  mu = (*linpredp);
  sigma2 = *worktransformlin[0];
  double z = (resp-mu)/sqrt(sigma2);
  res = -0.9189385-0.5*log(sigma2)-0.5*z*z;

  modify_worklin();
  linpredp++;
  return res;
  }*/


void DISTR_normal_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_sigma2
   // *linpred[1] = eta_mu
  //linpred has always size=2!
   if (*weight[1] == 0)
     *deviance=0;
   else
     {
     double sigma_2 = exp(*linpred[0]);
     double mu = (*linpred[1]);

     double l;

       l = -0.5*log(2*PI)-0.5*log(sigma_2)-pow((((*response[0]))-mu),2)/(2*sigma_2);


    *deviance = -2*l;
    }

  }


double DISTR_normal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_normal_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = (*linpred[copulaoffset+1]);
  }

 double DISTR_normal_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_normal_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    //weight, response and param has size>1 if copula model is specified!!
    double arg = ((*response[copulaoffset+1])-(*param[copulaoffset+1]))/(sqrt(*param[copulaoffset+0])) ;

    return (randnumbers::Phi2(arg));
    }

//double DISTR_normal_mu::compute_quantile_residual_mult(vector<double *> response,
//                                             vector<double *> param,
//                                             vector<double *> weight,
//                                             vector<datamatrix *> aux)
//  {
//  double u_est;
//  u_est = cdf_mult(response,param,weight,aux);
//  double res_est = randnumbers::invPhi2(u_est);
//  return res_est;
//  }

double DISTR_normal_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;

     l = -pow((((*response))-mu),2)/(2*(*worktransformlin[0]));

  if(optionsp->copula)
    {
    //implement loglik for copula models, i.e. add part logc
    double F = cdf(*response,*linpred);
    l += (distrcopulap[0]->logc(F,copulapos,false))[0];
    }

  modify_worklin();

  return l;

  }


void DISTR_normal_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;


  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double nu = ((*response)-mu)/(*worktransformlin[0]);

//  if(updateIWLS)
    *workingweight = 1/(*worktransformlin[0]);
//  else
//    *workingweight = (*weightp)/(*worktransformlin[0]);

  if(optionsp->copula)
    {
    double F = cdf(*response,*linpred);
    vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);
    if (compute_like)
      {
      like += logcandderivs[0];
      }
    // compute and implement dF/deta, d^2 F/deta ^2
    double z = ((*response)-mu)/sqrt((*worktransformlin[0]));
    double dF = -0.3989423*exp(-0.5*z*z)/(sqrt(*worktransformlin[0]));
    double ddF = -0.3989423*z*exp(-0.5*z*z)/((*worktransformlin[0]));
    /*if(copularotate)
      {
      dF = -dF;
      ddF = -ddF;
      }*/
    nu += logcandderivs[1]*dF;

   // *workingweight = (*worktransformlin[0])*(*worktransformlin[0])*hilfs1-logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
   *workingweight += -logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;
    if (*workingweight <=0)
      *workingweight = 0.001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((((*response))-mu),2)/(2*(*worktransformlin[0]));

      }


  modify_worklin();

  }


void DISTR_normal_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *mu = ((*linpred[copulaoffset+predstart_mumult+1]));
  }


void DISTR_normal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


/*void DISTR_normal_mu::update(void)
  {

   unsigned i;

  double help;

  double * worktransformlinp;
  double * workweight;

  worktransformlinp = distrp[0]->helpmat1.getV();
  workweight = workingweight.getV();

  if(updateIWLS)
    {
    for (i=0;i<nrobs;i++,worktransformlinp++,workweight++)
      {
      *workweight = 1/(*worktransformlinp);
      }
    }
  else
    {
    weightp = weight.getV();
    for (i=0;i<nrobs;i++,worktransformlinp++,workweight++,weightp++)
      {
      *workweight = (*weightp)/(*worktransformlinp);
      }
    }
  }*/

void DISTR_normal_mu::update_end(void)
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

void DISTR_normal_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  weightp = weight.getV();

  }



void DISTR_normal_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    weightp++;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_beta_sigma2 ---------------------------
//------------------------------------------------------------------------------


DISTR_beta_sigma2::DISTR_beta_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta Distribution - sigma2";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=10;
  maindistribution=false;
  check_errors();
  }


DISTR_beta_sigma2::DISTR_beta_sigma2(const DISTR_beta_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_beta_sigma2 & DISTR_beta_sigma2::operator=(
                            const DISTR_beta_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_beta_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_beta_sigma2::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = exp((*linpred[0]));
  *param = arg/(1+arg);
  }

double DISTR_beta_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu)/(1+exp(eta_mu));

  if (counter==0)
    {
    set_worklin();
    }

  double exp_lin = exp((*linpred));
  double sigma_2 = exp_lin/(1+exp_lin);
  double help = (1-sigma_2)/sigma_2;
  double mu_help = (*worktransformlin[0])*help;
  double one_minus_mu_help = (1-(*worktransformlin[0]))*help;

  double l;

     l = mu_help*log(*response) +
         one_minus_mu_help*log(1-(*response)) -
         randnumbers::lngamma_exact(mu_help) - randnumbers::lngamma_exact(one_minus_mu_help) +
         randnumbers::lngamma_exact(help);


  modify_worklin();

  return l;

  }

void DISTR_beta_sigma2::compute_iwls_wweightschange_weightsone(
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

    double exp_lin = exp(*linpred);

    double sigma_2 = exp_lin/(1+exp_lin);

    double help = (1-sigma_2)/sigma_2;
    double mu_help = (*worktransformlin[0])*help;
    double one_minus_mu_help = (1-(*worktransformlin[0]))*help;

    double nu = -help*( -(*worktransformlin[0])*randnumbers::digamma_exact(mu_help)
                       - (1-(*worktransformlin[0]))*randnumbers::digamma_exact(one_minus_mu_help) +
                       randnumbers::digamma_exact(help) + (*worktransformlin[0])*log(*response)
                       + (1-(*worktransformlin[0]))*log(1-(*response)) );



    *workingweight =  pow(help,2)*( pow((1-(*worktransformlin[0])),2)*randnumbers::trigamma_exact(one_minus_mu_help)+
                                   pow(((*worktransformlin[0])),2)*randnumbers::trigamma_exact(mu_help) -
                                   randnumbers::trigamma_exact(help) );

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += mu_help*log(*response) +
         one_minus_mu_help*log(1-(*response)) -
         randnumbers::lngamma_exact(mu_help) - randnumbers::lngamma_exact(one_minus_mu_help) +
         randnumbers::lngamma_exact(help);

      }

  modify_worklin();

  }


void DISTR_beta_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_beta_sigma2::update_end(void)
  {

  // helpmat1 stores (1-sigma2)/sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double exp_lin;
  double sigma_2;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    exp_lin = exp(*worklin);
    sigma_2 = exp_lin/(1+exp_lin);
    *pmu = (1-sigma_2)/sigma_2;
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_beta_mu -----------------------------
//------------------------------------------------------------------------------

void DISTR_beta_mu::check_errors(void)
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
        //cout << "workweight: " << *workweight << endl;
        if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
          }

        if (*workresp >= 1)
          {
          errors=true;
          errormessages.push_back("ERROR: response values greater or equal to one encountered\n");
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


DISTR_beta_mu::DISTR_beta_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "Beta Distribution - mu";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  linpredminlimit=-10;
  linpredmaxlimit=10;

  check_errors();
  }


DISTR_beta_mu::DISTR_beta_mu(const DISTR_beta_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_beta_mu & DISTR_beta_mu::operator=(
                            const DISTR_beta_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_beta_mu::compute_deviance_mult(vector<double *> response,
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
     double sigma_2 = exp(*linpred[0])/(1+exp(*linpred[0]));
     double mu = exp(*linpred[1])/(1+exp(*linpred[1]));
     double help = (1-sigma_2)/sigma_2;
     double one_minus_mu_help = (1-mu)*help;
     double mu_help = mu*help;

     double l;

       l = (mu_help-1)*log(*response[1]) +
			(one_minus_mu_help-1)*log(1-(*response[1]))-
			randnumbers::lngamma_exact(mu_help)-
			randnumbers::lngamma_exact(one_minus_mu_help)+
			randnumbers::lngamma_exact(help);


    *deviance = -2*l;

 /*   double a = sigma_2*mu;
    double b = sigma_2*(1-mu);
    std::ofstream out;
  // helpmat1.prettyPrint(out);
    for (int i=0; i<(a+b); i++) {
    out.open ("C:\\Urs\\nchoosek.raw", std::ofstream::out | std::ofstream::app);
    out << a ;
    out << " " ;
    out << b ;
    out << " " ;
    out << i ;
    out << " " ;
    out << randnumbers::n_choose_k(a+b-1,i) ;
    out << " " ;
    out << randnumbers::incomplete_beta(a,b,(*response[1])) << endl;
    out.close();
    }*/

    }

  }


double DISTR_beta_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

 double DISTR_beta_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

void DISTR_beta_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = exp((*linpred[1]));
  *param = arg/(1+arg);
  }

double DISTR_beta_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
//    double a =  (*param[0])*(*param[1]);
//    double b = (*param[0])*(1-(*param[1]));

     return 0;
//   return ( randnumbers::incomplete_beta(a,b,(*response[1])) );
    }

double DISTR_beta_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = (1-sigma2)/sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double exp_lin = exp((*linpred));
  double mu = exp_lin/(1+exp_lin);
  double mu_worktrans = mu*(*worktransformlin[0]);
  double one_minus_mu_worktrans = (1-mu)*(*worktransformlin[0]);

  double l;

     l = mu*(*worktransformlin[0])*log(*response) +
		 (1-mu)*(*worktransformlin[0])*log(1-(*response))-
		 randnumbers::lngamma_exact(mu_worktrans)-
		 randnumbers::lngamma_exact(one_minus_mu_worktrans);

  modify_worklin();

  return l;

  }


void DISTR_beta_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = (1-sigma2)/sigma2;

  if (counter==0)
    {
    set_worklin();
    }

    double exp_lin = exp(*linpred);

    double mu = exp_lin/(1+exp_lin);

    double mu_worktrans = mu*(*worktransformlin[0]);

    double one_minus_mu_worktrans = (1-mu)*(*worktransformlin[0]);

    double nu = mu*one_minus_mu_worktrans*log((*response))-mu*one_minus_mu_worktrans*log((1-(*response)))+
				mu*one_minus_mu_worktrans*(randnumbers::digamma_exact(one_minus_mu_worktrans)-randnumbers::digamma_exact(mu_worktrans));

    *workingweight = pow((*worktransformlin[0]),2)*pow(mu,2)*pow((1-mu),2)*(randnumbers::trigamma_exact(one_minus_mu_worktrans)+randnumbers::trigamma_exact(mu_worktrans));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += mu_worktrans*log((*response))+one_minus_mu_worktrans*log(1-(*response))-
		        randnumbers::lngamma_exact(mu_worktrans)-randnumbers::lngamma_exact((one_minus_mu_worktrans));

      }


  modify_worklin();

  }


void DISTR_beta_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double exp_lin = exp(*linpred[predstart_mumult+1]);
  *mu = exp_lin/(1+exp_lin);
  }


void DISTR_beta_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_beta_mu::update_end(void)
  {


  // helpmat1 stores exp(eta_mu)/(1+exp(eta_mu))

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double exp_lin;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    exp_lin = exp(*worklin);
    *pmu = exp_lin/(1+exp_lin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_cloglog -----------------------------
//------------------------------------------------------------------------------

void DISTR_cloglog::check_errors(void)
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

        if (((*workresp)!= 0.0) && ((*workresp)!= 1.0) )
          {
          errors=true;
          errormessages.push_back("ERROR: response has to be equal to zero or one\n");
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


DISTR_cloglog::DISTR_cloglog(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,0,w)
  {
  family = "Binomial Distribution - cloglog";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "pi";
  linpredminlimit=-7.5;
  linpredmaxlimit=2.2;
  check_errors();
  }


DISTR_cloglog::DISTR_cloglog(const DISTR_cloglog & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_cloglog & DISTR_cloglog::operator=(
                            const DISTR_cloglog & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_cloglog::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response


   if (*weight[0] == 0)
     *deviance=0;
   else
     {
     double oneminuspi = exp(-exp((*linpred[0])));
     double pi = 1 - oneminuspi;

     double l;

     if ((*response[0])==0) {
        l = log(oneminuspi);
     } else {
        l = log(pi);
     }


    *deviance = -2*l;
    }

  }


double DISTR_cloglog::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_cloglog::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = exp((*linpred[0]));
  *param = 1-exp(-arg);
  }

double DISTR_cloglog::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
  {
  double Fy = 0;
  if (*response[0]>0)
    Fy = 1;
  else
    Fy = 1-(*param[0]);
  return Fy;
  }

double DISTR_cloglog::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_cloglog::compute_quantile_residual_mult(vector<double *> response,
                                         vector<double *> param,
                                         vector<double *> weight,
                                          vector<datamatrix *> aux)
    {
    double u_est;
    if(*response[0]<=0) {
        u_est = randnumbers::uniform_ab(0,1-(*param[0]));
    }
    else {
        u_est = randnumbers::uniform_ab(1-(*param[0]), 1);
    }
    double res_est = randnumbers::invPhi2(u_est);
    return res_est;
    }


double DISTR_cloglog::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = (1-sigma2)/sigma2;

  if (counter==0)
    {
    set_worklin();
    }

     double oneminuspi = exp(-exp((*linpred)));
     double pi = 1 - oneminuspi;

  double l;

   if ((*response)==0) {
        l = log(oneminuspi);
     } else {
        l = log(pi);
     }

  modify_worklin();

  return l;

  }


void DISTR_cloglog::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

     double oneminuspi = exp(-exp((*linpred)));
     double pi = 1 - oneminuspi;


    double nu = exp((*linpred))*oneminuspi/pi;

    if ((*response)==0) {
        nu -= exp((*linpred))/pi;
    }

    *workingweight = pow(exp((*linpred)),2)*oneminuspi/pi;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

   if (*response==0)
      like += log(oneminuspi);
    else // response == 1
      like += log(pi);

      }


  modify_worklin();

  }


void DISTR_cloglog::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double exp_lin = exp(*linpred[predstart_mumult]);
  *mu = 1-exp(-exp_lin);
  }


void DISTR_cloglog::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (pi):cloglog\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_cloglog::update_end(void)
  {


  // helpmat1 stores 1-exp(-exp_lin)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double exp_lin;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    exp_lin = exp(*worklin);
    *pmu = 1-exp(-exp_lin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//----------------- CLASS: DISTR_binomialprobit_copula -------------------------
//------------------------------------------------------------------------------

void DISTR_binomialprobit_copula::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = responseorig.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {

        if (((*workresp)!= 0.0) && ((*workresp)!= 1.0) )
          {
          errors=true;
          errormessages.push_back("ERROR: response has to be equal to zero or one\n");
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


DISTR_binomialprobit_copula::DISTR_binomialprobit_copula(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,0,w)
  {
  family = "Binomial Distribution - probit";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  responseorig = r;
  linpredminlimit=-7;
  linpredmaxlimit=7;
  posteriormodemode=true;
  updateIWLS = true;
  check_errors();

  }


DISTR_binomialprobit_copula::DISTR_binomialprobit_copula(const DISTR_binomialprobit_copula & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  responseorig = nd.responseorig;
  responsecopmat = nd.responsecopmat;
  posteriormodemode = nd.posteriormodemode;
  }


const DISTR_binomialprobit_copula & DISTR_binomialprobit_copula::operator=(
                            const DISTR_binomialprobit_copula & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  responseorig = nd.responseorig;
  responsecopmat = nd.responsecopmat;
  posteriormodemode = nd.posteriormodemode;
  return *this;
  }

/*bool DISTR_binomialprobit_copula::posteriormode(void)
  {
  update();
  return true;
  }*/


void DISTR_binomialprobit_copula::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
  //weight, linpred and response have only size=1!!!
   if (*weight[0] == 0)
     *deviance=0;
   else
     {
     double pi = randnumbers::Phi2(*linpred[0]);

     double l;

     if ((*response[0])<=0) {
        l = log(1-pi);
     } else {
        l = log(pi);
     }
    *deviance = -2*l;
    }

  }

double DISTR_binomialprobit_copula::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_binomialprobit_copula::compute_param_mult(vector<double *>  linpred,double * param)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *param = randnumbers::Phi2(*linpred[copulaoffset + 0]);
  }

double DISTR_binomialprobit_copula::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
  {
  //weight, response and param has size>1 if copula model is specified!!
  double Fy = 0;
  if (*response[copulaoffset + 0]>0)
    Fy = 1;
  else
    Fy = 1-(*param[copulaoffset + 0]);

  return Fy;
  }

double DISTR_binomialprobit_copula::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    //weight, response and param has size>1 if copula model is specified!!
    return 0;
    }

double DISTR_binomialprobit_copula::compute_quantile_residual_mult(vector<double *> response,
                                         vector<double *> param,
                                         vector<double *> weight,
                                          vector<datamatrix *> aux)
    {
    //weight, response and param has size>1 if copula model is specified!!
    double u_est;
    if(*response[copulaoffset + 0]<=0) {
        u_est = randnumbers::uniform_ab(0,1-(*param[copulaoffset + 0]));
    }
    else {
        u_est = randnumbers::uniform_ab(1-(*param[copulaoffset + 0]), 1);
    }
    double res_est = randnumbers::invPhi2(u_est);
    return res_est;
    }


double DISTR_binomialprobit_copula::cdf(const double & resp, const bool & ifcop)
  {
  if(counter==0)
    {
    if(ifcop)
      {
      set_worklin();
      }
    if (linpred_current==1)
      linpredp = linearpred1.getV();
    else
      linpredp = linearpred2.getV();
    }

  double res;
  res = randnumbers::Phi2(resp-*linpredp);

/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Probit CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Probit CDF < 0.001" << endl;
    }*/

  if(ifcop)
    {
    modify_worklin();
    }
  linpredp++;
  return res;
  }

double DISTR_binomialprobit_copula::cdf(const double & resp, const double & linpred)
  {
  double res = randnumbers::Phi2(resp-linpred);

/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Probit CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Probit CDF < 0.001" << endl;
    }*/

  return res;
  }

double DISTR_binomialprobit_copula::cdf(const double & resp, vector<double *> linpred)
  {
  //linpred has only size=1!
  double res = randnumbers::Phi2(resp-*linpred[0]);

/*  if(res>0.999)
    {
    res = 0.999;
    cout << "counter: " << counter << endl;
    cout << "Probit CDF > 0.999" << endl;
    }
  if(res<0.001)
    {
    res = 0.001;
    cout << "counter: " << counter << endl;
    cout << "Probit CDF < 0.001" << endl;
    }*/

  return res;
  }

double DISTR_binomialprobit_copula::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }

  double F = cdf(*response,*linpred);

  vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,false);

  if(optionsp->samplesel)
    {
    if(*response<=0)
      {
      logcandderivs[0] = 0;
      }
    }

  double l =-0.5*(*response-*linpred)*(*response-*linpred) + logcandderivs[0];

  modify_worklin();

  return l;
  }

void DISTR_binomialprobit_copula::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  if (counter==0)
    {
    set_worklin();
    }

  double F = cdf(*response,*linpred);

  vector<double> logcandderivs = distrcopulap[0]->logc(F,copulapos,true);

  if(optionsp->samplesel)
    {
    if(*response<=0)
      {
      for(unsigned j=0; j<logcandderivs.size(); j++)
        {
        logcandderivs[j] = 0;
        }
      }
    }

   //   cout << "distrcopulap.size():" << distrcopulap.size() << endl;
  // compute and implement dF/deta, d^2 F/deta ^2
  double z = (*response-*linpred);
  double dF = -0.3989423*exp(-0.5*z*z);
  double ddF = z*dF;

  double nu = *response - *linpred + logcandderivs[1]*dF;
  *workingweight = 1-logcandderivs[2]*dF*dF-logcandderivs[1]*ddF;

  if (*workingweight <=0)
    *workingweight = 0.0001;

/*  if(isnan(*workingweight) || *workingweight > 100)
    {
    *workingweight = 1.0;
    cout << "counter: " << counter << endl;
    cout << "Probit *workingweight NAN" << endl;
    }*/

  *workingresponse = *linpred + nu/(*workingweight);

/*  if(isnan(*workingresponse) || *workingresponse > 100 || *workingresponse < -100)
    {
    *workingresponse = 0.0;
    cout << "counter: " << counter << endl;
    cout << "Probit *workingresponse NAN" << endl;
    }*/

  if (compute_like)
    {
    like += -0.5*(*response-*linpred)*(*response-*linpred) + logcandderivs[0];
    }
  modify_worklin();
  }


void DISTR_binomialprobit_copula::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  //weight, response and param has size>1 if copula model is specified!!
  *mu = randnumbers::Phi2(*linpred[copulaoffset + predstart_mumult]);
  }


void DISTR_binomialprobit_copula::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): probit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_binomialprobit_copula::update_end(void)
  {
  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = randnumbers::Phi2(*worklin);
    }
  }

void DISTR_binomialprobit_copula::update(void)
  {
  posteriormodemode = false;

  double * workresporig = responseorig.getV();
  double * workresp = response.getV();
  double * weightwork = weight.getV();
  double * weightworkcopula = (distrcopulap[0]->weight).getV();

  double * worklin_current;
  if (linpred_current==1)
    worklin_current = linearpred1.getV();
  else
    worklin_current = linearpred2.getV();

  double * responsecop = responsecopmat->getV();

  unsigned i;
  for(i=0;i<nrobs;i++,worklin_current++,workresp++,weightwork++,responsecop++,workresporig++,weightworkcopula++)
    {
    double x = randnumbers::uniform();
    *workresp = distrcopulap[0]->condfc(x, *worklin_current, *workresporig, copulapos);
    *responsecop = *workresp;
    }
//  ofstream out1("c://temp//samplesel3//response_generic.raw");
//  response.prettyPrint(out1);
//  out1.close();
  }

//--------------------------------------------------------------------------------
//----------------- CLASS: DISTR_claytoncopula2_normal_sigma2 --------------------
//--------------------------------------------------------------------------------


DISTR_claytoncopula2_normal_sigma2::DISTR_claytoncopula2_normal_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Clayton Copula - sigma2";
  pos =p;
  outpredictor = true;
  outexpectation = true;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_claytoncopula2_normal_sigma2::DISTR_claytoncopula2_normal_sigma2(const DISTR_claytoncopula2_normal_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_claytoncopula2_normal_sigma2 & DISTR_claytoncopula2_normal_sigma2::operator=(
                            const DISTR_claytoncopula2_normal_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_claytoncopula2_normal_sigma2::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_claytoncopula2_normal_sigma2::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_claytoncopula2_normal_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_claytoncopula2_normal_sigma2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = pow(exp((*linpred[2*pos])),0.5);
  }

double DISTR_claytoncopula2_normal_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = exp((*linpred));
  double arg1 = ((*response) - (*worktransformlin[2])) / pow(sigma_2, 0.5);
  double arg2 = ((*response2p) - (*worktransformlin[3])) / pow((*worktransformlin[1]), 0.5);
  double prop1 = randnumbers::Phi2(arg1);
  double prop2 = randnumbers::Phi2(arg2);

  double l;

     l = -0.5*log(sigma_2)-pow((((*response))-(*worklin[2])),2)/(2*sigma_2)- (1 + (*worktransformlin[0])) * log(prop1)
          - (2 + 1 / (*worktransformlin[0])) * log(pow(prop1, -(*worktransformlin[0])) + pow(prop2, -(*worktransformlin[0])) -1);


  modify_worklin();

  return l;

  }

void DISTR_claytoncopula2_normal_sigma2::compute_iwls_wweightschange_weightsone(
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

    double sigma_2 = exp(*linpred);
    double arg1 = ((*response) - (*worktransformlin[2])) / pow(sigma_2, 0.5);
    double arg2 = ((*response2p) - (*worktransformlin[3])) / pow((*worktransformlin[1]), 0.5);
    double prop1 = randnumbers::Phi2(arg1);
    double prop2 = randnumbers::Phi2(arg2);
    double hilfs = pow(prop1, -(*worktransformlin[0])) + pow(prop2, -(*worktransformlin[0])) - 1;
    double hilfs2 = (2 * (*worktransformlin[0]) + 1);

    double d1 = - 0.398942280401433 * exp(- 0.5 * pow(arg1, 2)) * arg1 * 0.5 ;
    double d2 = - 0.5 * d1 * (1 - pow(arg1, 2));


    double nu = -0.5 + (pow(((*response)-(*worklin[2])),2))/(2*sigma_2) - (1 + (*worktransformlin[0])) * d1 / prop1
                + hilfs2 * pow(prop1, -(*worktransformlin[0])-1) * d1 / hilfs;



    *workingweight = 0.5 - (1 + (*worktransformlin[0])) * pow(d1 / prop1, 2)
                    + (1 + (*worktransformlin[0])) * d2 / prop1 + hilfs2 * ((*worktransformlin[0]) + 1) * pow(prop1, (-(*worktransformlin[0])-2)) * pow(d1, 2) / hilfs
                    - hilfs2 * pow(prop1, -(*worktransformlin[0])-1) * d2 / hilfs
                    - (*worktransformlin[0]) * hilfs2 * pow(pow(prop1, -(*worktransformlin[0])-1) * d1 / hilfs, 2);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((((*response))-(*worklin[2])),2)/(2*sigma_2)- (1 + (*worktransformlin[0])) * log(prop1)
          - (2 + 1 / (*worktransformlin[0])) * log(hilfs);

      }

  modify_worklin();

  }

void DISTR_claytoncopula2_normal_sigma2::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double arg = - 1 / exp((*linpred[predstart_mumult+4]));
  *mu = pow(2, arg);
  }


void DISTR_claytoncopula2_normal_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_claytoncopula2_normal_sigma2::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_claytoncopula2_normal_mu ----------------
//------------------------------------------------------------------------------


DISTR_claytoncopula2_normal_mu::DISTR_claytoncopula2_normal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  pos=p;
  family = "Clayton Copula - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_claytoncopula2_normal_mu::DISTR_claytoncopula2_normal_mu(const DISTR_claytoncopula2_normal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
    pos = nd.pos;
    response2 = nd.response2;
    response2p = nd.response2p;
  }


const DISTR_claytoncopula2_normal_mu & DISTR_claytoncopula2_normal_mu::operator=(
                            const DISTR_claytoncopula2_normal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_claytoncopula2_normal_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_claytoncopula2_normal_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_claytoncopula2_normal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_claytoncopula2_normal_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2*pos + 1]);
  }

 double DISTR_claytoncopula2_normal_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_claytoncopula2_normal_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    return 0;
    }

//double DISTR_claytoncopula2_normal_mu::compute_quantile_residual_mult(vector<double *> response,
//                                             vector<double *> param,
//                                             vector<double *> weight,
//                                             vector<datamatrix *> aux)
//  {
//  double u_est;
//  u_est = cdf_mult(response,param,weight,aux);
//  double res_est = randnumbers::invPhi2(u_est);
//  return res_est;
//  }

double DISTR_claytoncopula2_normal_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;
  double prop1 = randnumbers::Phi2(((*response) - mu) / pow((*worktransformlin[2]), 0.5));
  double prop2 = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[3]), 0.5));

  l = -pow((((*response))-mu),2)/(2*(*worktransformlin[2])) - (1 + (*worktransformlin[0])) * log(prop1)
          - (2 + 1 / (*worktransformlin[0])) * log(pow(prop1, -(*worktransformlin[0])) + pow(prop2, -(*worktransformlin[0])) -1);

  modify_worklin();

  return l;

  }


void DISTR_claytoncopula2_normal_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;


  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);
    double arg1 = ((*response) - mu) / pow((*worktransformlin[2]), 0.5);
    double arg2 = ((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[3]), 0.5);
    double prop1 = randnumbers::Phi2(arg1);
    double prop2 = randnumbers::Phi2(arg2);
    double hilfs = pow(prop1, -(*worktransformlin[0])) + pow(prop2, -(*worktransformlin[0])) - 1;
    double hilfs2 = (2 * (*worktransformlin[0]) + 1);

    double d1 = - 0.398942280401433 * exp(- 0.5 * pow(arg1, 2)) / pow((*worktransformlin[2]), 0.5);
    double d2 = d1 * arg1 / pow((*worktransformlin[2]), 0.5);

    double nu = ((*response)-mu)/(*worktransformlin[2]) - (1 + (*worktransformlin[0])) * d1 / prop1
                + hilfs2 * pow(prop1, -(*worktransformlin[0])-1) * d1 / hilfs;

    *workingweight = 1/(*worktransformlin[2]) - (1 + (*worktransformlin[0])) * pow(d1 / prop1, 2)
                    + (1 + (*worktransformlin[0])) * d2 / prop1 + hilfs2 * ((*worktransformlin[0]) + 1) * pow(prop1, (-(*worktransformlin[0])-2)) * pow(d1, 2) / hilfs
                    - hilfs2 * pow(prop1, -(*worktransformlin[0])-1) * d2 / hilfs
                    - (*worktransformlin[0]) * hilfs2 * pow(pow(prop1, -(*worktransformlin[0])-1) * d1 / hilfs, 2);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((((*response))-mu),2)/(2*(*worktransformlin[2])) - (1 + (*worktransformlin[0])) * log(prop1)
                - (2 + 1 / (*worktransformlin[0])) * log(hilfs);

      }


  modify_worklin();

  }


void DISTR_claytoncopula2_normal_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+2 * pos+1]));
  }


void DISTR_claytoncopula2_normal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


/*void DISTR_claytoncopula2_normal_mu::update(void)
  {

   unsigned i;

  double help;

  double * worktransformlinp;
  double * workweight;

  worktransformlinp = distrp[0]->helpmat1.getV();
  workweight = workingweight.getV();

  for (i=0;i<nrobs;i++,worktransformlinp++,workweight++)
    {
        *workweight = 1/(*worktransformlinp);
    }

  }*/


void DISTR_claytoncopula2_normal_mu::update_end(void)
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


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_claytoncopula2_rho --------------------
//------------------------------------------------------------------------------

void DISTR_claytoncopula2_rho::check_errors(void)
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

DISTR_claytoncopula2_rho::DISTR_claytoncopula2_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Clayton Copula - rho";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_claytoncopula2_rho::DISTR_claytoncopula2_rho(const DISTR_claytoncopula2_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_claytoncopula2_rho & DISTR_claytoncopula2_rho::operator=(
                            const DISTR_claytoncopula2_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_claytoncopula2_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_claytoncopula2_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = exp(*linpred[4]);
  *param = arg;
  }

void DISTR_claytoncopula2_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_claytoncopula2_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_claytoncopula2_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {


   if (*weight[4] == 0)
     *deviance=0;
   else
     {
     double rho = exp((*linpred[4]));
     double u = randnumbers::Phi2(((*response[3]) - (*linpred[3])) / pow(exp(*linpred[2]), 0.5));
     double v = randnumbers::Phi2(((*response[0]) - (*linpred[1])) / pow(exp(*linpred[0]), 0.5));
     double urho = pow(u, -rho);
     double vrho = pow(v, -rho);
     double arg = urho + vrho - 1;
     double l;

       l = log(rho + 1) - (1 + rho) * (log(u) + log(v)) - (2 + 1 / rho) * log(arg)
            +log(randnumbers::phi(((*response[3]) - (*linpred[3])) / pow(exp(*linpred[2]), 0.5))) +
                log(randnumbers::phi(((*response[0]) - (*linpred[1])) / pow(exp(*linpred[0]), 0.5)));


    *deviance = -2*l;
    }

  }

double DISTR_claytoncopula2_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {



  if (counter==0)
    {
    set_worklin();
    }
    double rho = exp((*linpred));
    double u = randnumbers::Phi2(((*response) - (*worktransformlin[3])) / pow((*worktransformlin[2]), 0.5));
    double v = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[0]), 0.5));
    double logu = log(u);
    double logv = log(v);
    double urho = pow(u, -rho);
    double vrho = pow(v, -rho);
    double arg = urho + vrho - 1;
    double l;

    l = log(rho + 1) - (1 + rho) * (logu + logv) - (2 + 1 / rho) * log(arg);

  modify_worklin();

  return l;

  }

void DISTR_claytoncopula2_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }

    double rho = exp((*linpred));
    double u = randnumbers::Phi2(((*response) - (*worktransformlin[3])) / pow((*worktransformlin[2]), 0.5));
    double v = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[0]), 0.5));
    double logu = log(u);
    double logv = log(v);
    double urho = pow(u, -rho);
    double vrho = pow(v, -rho);
    double arg = urho + vrho - 1;

    double nu = rho / (rho + 1) - rho * (logu + logv) + log(arg) / rho + (2 * rho + 1) * (logu * urho + logv * vrho) / arg;

    *workingweight = -rho / pow(rho + 1, 2) + rho * (logu + logv) + log(arg) / rho + (1 - 2 * rho) * (logu * urho + logv * vrho) / arg
                        - ((pow(rho, 2) * (2 + 1 / rho)) / (arg)) * (pow((logu * urho + logv * vrho), 2) / arg - pow(logu, 2) * urho - pow(logv, 2) * vrho );

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

          like += log(rho + 1) - (1 + rho) * (logu + logv) - (2 + 1 / rho) * log(arg);
      }


  modify_worklin();

  }

void DISTR_claytoncopula2_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double arg = exp((*linpred[predstart_mumult+4]));
  *mu = arg / (arg + 2);
  }


void DISTR_claytoncopula2_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_claytoncopula2_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_claytoncopula_rho ---------------------
//------------------------------------------------------------------------------

void DISTR_claytoncopula_rho::check_errors(void)
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

            if((*workresp) > 1 )
            {
                errors=true;
                errormessages.push_back("ERROR: cdfs of marginals take values inbetween zero and one!\n");
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

DISTR_claytoncopula_rho::DISTR_claytoncopula_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Clayton Copula - rho";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();

  }


DISTR_claytoncopula_rho::DISTR_claytoncopula_rho(const DISTR_claytoncopula_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_claytoncopula_rho & DISTR_claytoncopula_rho::operator=(
                            const DISTR_claytoncopula_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_claytoncopula_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_claytoncopula_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[2]));
  }

void DISTR_claytoncopula_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_claytoncopula_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_claytoncopula_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {


   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double rho = exp((*linpred[2]));
     double logu = log((*response[1]));
     double logv = log((*response[0]));
     double urho = pow((*response[1]), -rho);
     double vrho = pow((*response[0]), -rho);
     double arg = urho + vrho - 1;
     double l;

       l = log(rho + 1) - (1 + rho) * (logu + logv) - (2 + 1 / rho) * log(arg);


    *deviance = -2*l;
    }

  }

double DISTR_claytoncopula_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {



  if (counter==0)
    {
    set_worklin();
    }
    double rho = exp((*linpred));
    double logu = log((*response));
    double logv = log((*response2p));
    double urho = pow((*response), -rho);
    double vrho = pow((*response2p), -rho);
    double arg = urho + vrho - 1;
    double l;

    l = log(rho + 1) - (1 + rho) * (logu + logv) - (2 + 1 / rho) * log(arg);


  modify_worklin();

  return l;

  }

void DISTR_claytoncopula_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }

    double rho = exp((*linpred));
    double logu = log((*response));
    double logv = log((*response2p));
    double urho = pow((*response), -rho);
    double vrho = pow((*response2p), -rho);
    double arg = urho + vrho - 1;

    double nu = rho / (rho + 1) - rho * (logu + logv) + log(arg) / rho + (2 * rho + 1) * (logu * urho + logv * vrho) / arg;

    *workingweight = -rho / pow(rho + 1, 2) + rho * (logu + logv) + log(arg) / rho + (1 - 2 * rho) * (logu * urho + logv * vrho) / arg
                        - ((pow(rho, 2) * (2 + 1 / rho)) / (arg)) * (pow((logu * urho + logv * vrho), 2) / arg - pow(logu, 2) * urho - pow(logv, 2) * vrho );

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

          like += log(rho + 1) - (1 + rho) * (logu + logv) - (2 + 1 / rho) * log(arg);
      }

  modify_worklin();

  }

void DISTR_claytoncopula_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double arg = exp((*linpred[predstart_mumult+2]));
  *mu = arg / (arg + 2);
  }


void DISTR_claytoncopula_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_claytoncopula_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }

//--------------------------------------------------------------------------------
//----------------- CLASS: DISTR_gumbelcopula2_normal_sigma2 --------------------
//--------------------------------------------------------------------------------


DISTR_gumbelcopula2_normal_sigma2::DISTR_gumbelcopula2_normal_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "gumbel Copula - sigma2";
  pos =p;
  outpredictor = true;
  outexpectation = true;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_gumbelcopula2_normal_sigma2::DISTR_gumbelcopula2_normal_sigma2(const DISTR_gumbelcopula2_normal_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_gumbelcopula2_normal_sigma2 & DISTR_gumbelcopula2_normal_sigma2::operator=(
                            const DISTR_gumbelcopula2_normal_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_gumbelcopula2_normal_sigma2::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gumbelcopula2_normal_sigma2::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_gumbelcopula2_normal_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gumbelcopula2_normal_sigma2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[2*pos]));
  }

double DISTR_gumbelcopula2_normal_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

    double sigma_2 = exp((*linpred));
    double arg1 = ((*response) - (*worktransformlin[2])) / pow(sigma_2, 0.5);
    double arg2 = ((*response2p) - (*worktransformlin[3])) / pow((*worktransformlin[1]), 0.5);
    double u = randnumbers::Phi2(arg1);
    double v = randnumbers::Phi2(arg2);

    double logu = log(u);
    double logv = log(v);
    double logurho = pow(-logu, (*worktransformlin[0]));
    double logvrho = pow(-logv, (*worktransformlin[0]));
    double arg = logurho + logvrho;
    double lcopula = (-pow(arg, (1 / (*worktransformlin[0]))));
    double l;

    l = -0.5*log(sigma_2)-pow((((*response))-(*worklin[2])),2)/(2*sigma_2) + lcopula + ((*worktransformlin[0]) -1) * (log(-logu) + log(-logv)) - logu - logv +
                (2 / (*worktransformlin[0]) - 2) * log(logurho + logvrho) + log(1 + ((*worktransformlin[0]) - 1) * pow((logurho + logvrho), (-1 / (*worktransformlin[0]))));

  modify_worklin();

  return l;

  }

void DISTR_gumbelcopula2_normal_sigma2::compute_iwls_wweightschange_weightsone(
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

    double sigma_2 = exp(*linpred);
    double arg1 = ((*response) - (*worktransformlin[2])) / pow(sigma_2, 0.5);
    double arg2 = ((*response2p) - (*worktransformlin[3])) / pow((*worktransformlin[1]), 0.5);
    double u = randnumbers::Phi2(arg1);
    double v = randnumbers::Phi2(arg2);
    double logu = log(u);
    double logv = log(v);
    double logurho = pow(-logu, (*worktransformlin[0]));
    double logvrho = pow(-logv, (*worktransformlin[0]));
    double arg = logurho + logvrho;

//    double dd1 = - 0.398942280401433 * exp(- 0.5 * pow(arg1, 2)) * arg1 * 0.5 ;
//    double dd2 = - 0.5 * dd1 * (1 - pow(arg1, 2));

//    double darg = -(*worktransformlin[0]) * dd1 * pow(-logu, ((*worktransformlin[0]) - 1)) / u;
//    double ddarg = - darg * ((dd1 / u) * (((*worktransformlin[0]) - 1) * pow(-logu, -1) + 1) - dd2 / dd1);
//    double zer = 1 + ((*worktransformlin[0]) - 1) * pow(arg, (- 1 / (*worktransformlin[0])));
//    double ner = (1 / (*worktransformlin[0])) * ((*worktransformlin[0]) - 1) * pow(arg, (-1 / (*worktransformlin[0]) - 1));
//    double rest1 = - pow(arg, (1 / (*worktransformlin[0]) - 1)) / (*worktransformlin[0])
//                    + 2 * (1 / (*worktransformlin[0]) - 1) / arg - ner / zer;
//    double rest2 = ((*worktransformlin[0]) - 1) / logu - 1;
//    double rest3 = - (1 / (*worktransformlin[0]) - 1) * pow(arg, (1 / (*worktransformlin[0]) - 2)) / (*worktransformlin[0])
//                    - 2 * (1 / (*worktransformlin[0]) - 1) / pow(arg, 2)
//                    - (1 / (*worktransformlin[0])) * (-1 / (*worktransformlin[0]) - 1) * ((*worktransformlin[0]) - 1) * pow(arg, (-1 / (*worktransformlin[0]) - 2)) / zer
 //                   - pow(((ner) / (zer)), 2);

    double nu = -0.5 + (pow(((*response)-(*worklin[2])),2))/(2*sigma_2) ;
    //+ darg * rest1 + dd1 * rest2 / u;

    *workingweight = 0.5 ;
    //- ddarg * rest1 - pow(darg, 2) * rest3 + pow((dd1 / u), 2) * (rest2 + ((*worktransformlin[0])- 1) / pow(logu, 2)) - dd2 * rest2 / u;

    *workingresponse = *linpred + nu/(*workingweight);
 /*             std::ofstream out;
  // helpmat1.prettyPrint(out);
   out.open ("C:\\tmp\\3.raw", std::ofstream::out | std::ofstream::app);
    out << *linpred ;
    out << " " ;
    out << *workingweight ;
    out << " " ;
    out << nu  ;
    out << " " ;
    out << u << endl;
    out.close();*/
    if (compute_like)

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((((*response))-(*worklin[2])),2)/(2*sigma_2)
                    -pow(arg, (1 / (*worktransformlin[0]))) + ((*worktransformlin[0]) -1) * (log(-logu) + log(-logv)) - logu - logv +
                (2 / (*worktransformlin[0]) - 2) * log(arg) + log(1 + ((*worktransformlin[0]) - 1) * pow(arg, (-1 / (*worktransformlin[0]))));

      }

  modify_worklin();

  }

void DISTR_gumbelcopula2_normal_sigma2::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_gumbelcopula2_normal_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gumbelcopula2_normal_sigma2::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gumbelcopula2_normal_mu ----------------
//------------------------------------------------------------------------------


DISTR_gumbelcopula2_normal_mu::DISTR_gumbelcopula2_normal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  pos=p;
  family = "gumbel Copula - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_gumbelcopula2_normal_mu::DISTR_gumbelcopula2_normal_mu(const DISTR_gumbelcopula2_normal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
    pos = nd.pos;
    response2 = nd.response2;
    response2p = nd.response2p;
  }


const DISTR_gumbelcopula2_normal_mu & DISTR_gumbelcopula2_normal_mu::operator=(
                            const DISTR_gumbelcopula2_normal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_gumbelcopula2_normal_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gumbelcopula2_normal_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_gumbelcopula2_normal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gumbelcopula2_normal_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2*pos + 1]);
  }

 double DISTR_gumbelcopula2_normal_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_gumbelcopula2_normal_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    return 0;
    }

//double DISTR_gumbelcopula2_normal_mu::compute_quantile_residual_mult(vector<double *> response,
//                                             vector<double *> param,
//                                             vector<double *> weight,
//                                             vector<datamatrix *> aux)
//  {
//  double u_est;
//  u_est = cdf_mult(response,param,weight,aux);
//  double res_est = randnumbers::invPhi2(u_est);
//  return res_est;
//  }

double DISTR_gumbelcopula2_normal_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);

    double u = randnumbers::Phi2(((*response) - mu) / pow((*worktransformlin[2]), 0.5));
    double v = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[3]), 0.5));

    double logu = log(u);
    double logv = log(v);
    double logurho = pow(-logu, (*worktransformlin[0]));
    double logvrho = pow(-logv, (*worktransformlin[0]));
    double arg = logurho + logvrho;
    double lcopula = (-pow(arg, (1 / (*worktransformlin[0]))));
    double l;

    l = -pow((((*response))-mu),2)/(2*(*worktransformlin[2])) + lcopula + ((*worktransformlin[0]) -1) * (log(-logu)) - logu +
                (2 / (*worktransformlin[0]) - 2) * log(arg) + log(1 + ((*worktransformlin[0]) - 1) * pow((arg), (-1 / (*worktransformlin[0]))));

    modify_worklin();

    return l;

  }


void DISTR_gumbelcopula2_normal_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;


  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);
    double arg1 = ((*response) - mu) / pow((*worktransformlin[2]), 0.5);
    double arg2 = ((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[3]), 0.5);
    double u = randnumbers::Phi2(arg1);
    double v = randnumbers::Phi2(arg2);
    double logu = log(u);
    double logv = log(v);
    double logurho = pow(-logu, (*worktransformlin[0]));
    double logvrho = pow(-logv, (*worktransformlin[0]));
    double arg = logurho + logvrho;

//    double dd1 = - 0.398942280401433 * exp(- 0.5 * pow(arg1, 2)) / pow((*worktransformlin[2]), 0.5);
//    double dd2 = dd1 * arg1 / pow((*worktransformlin[2]), 0.5);
//    double darg = -(*worktransformlin[0]) * dd1 * pow(-logu, ((*worktransformlin[0]) - 1)) / u;
//    double ddarg = - darg * ((dd1 / u) * (((*worktransformlin[0]) - 1) * pow(-logu, -1) + 1) - dd2 / dd1);
//    double zer = 1 + ((*worktransformlin[0]) - 1) * pow(arg, (- 1 / (*worktransformlin[0])));
//    double ner = (1 / (*worktransformlin[0])) * ((*worktransformlin[0]) - 1) * pow(arg, (-1 / (*worktransformlin[0]) - 1));
//    double rest1 = - pow(arg, (1 / (*worktransformlin[0]) - 1)) / (*worktransformlin[0])
//                    + 2 * (1 / (*worktransformlin[0]) - 1) / arg - ner / zer;
//    double rest2 = ((*worktransformlin[0]) - 1) / logu - 1;
//    double rest3 = - (1 / (*worktransformlin[0]) - 1) * pow(arg, (1 / (*worktransformlin[0]) - 2)) / (*worktransformlin[0])
//                    - 2 * (1 / (*worktransformlin[0]) - 1) / pow(arg, 2)
//                    - (1 / (*worktransformlin[0])) * (-1 / (*worktransformlin[0]) - 1) * ((*worktransformlin[0]) - 1) * pow(arg, (-1 / (*worktransformlin[0]) - 2)) / zer
//                    - pow(((ner) / (zer)), 2);

    double nu = ((*response)-mu)/(*worktransformlin[2]);
    //  + darg * rest1 + dd1 * rest2 / u;

    *workingweight = 1/(*worktransformlin[2]);
    // - ddarg * rest1 - pow(darg, 2) * rest3
      //                + pow((dd1 / u), 2) * (rest2 + ((*worktransformlin[0])- 1) / pow(logu, 2)) - dd2 * rest2 / u;

    *workingresponse = *linpred + nu/(*workingweight);

 /*             std::ofstream out;
  // helpmat1.prettyPrint(out);
   out.open ("C:\\tmp\\2.raw", std::ofstream::out | std::ofstream::app);
    out << *linpred ;
    out << " " ;
    out << *workingweight ;
    out << " " ;
    out << nu  ;
    out << " " ;
    out << u << endl;
    out.close();*/
    if (compute_like)
      {

        like += -pow((((*response))-mu),2)/(2*(*worktransformlin[2])) -pow(arg, (1 / (*worktransformlin[0]))) + ((*worktransformlin[0]) -1) * (log(-logu) + log(-logv)) - logu - logv +
                (2 / (*worktransformlin[0]) - 2) * log(arg) + log(1 + ((*worktransformlin[0]) - 1) * pow((arg), (-1 / (*worktransformlin[0]))));
      }


  modify_worklin();

  }


void DISTR_gumbelcopula2_normal_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+2 * pos+1]));
  }


void DISTR_gumbelcopula2_normal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


/*void DISTR_gumbelcopula2_normal_mu::update(void)
  {

   unsigned i;

  double help;

  double * worktransformlinp;
  double * workweight;

  worktransformlinp = distrp[0]->helpmat1.getV();
  workweight = workingweight.getV();

  for (i=0;i<nrobs;i++,worktransformlinp++,workweight++)
    {
        *workweight = 1/(*worktransformlinp);
    }

  }*/


void DISTR_gumbelcopula2_normal_mu::update_end(void)
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
//    double t = *pmu;
    *pmu = (*worklin);
//    double t = 0;
    }

  }


//--------------------------------------------------------------------------------
//----------------- CLASS: DISTR_gumbelcopula2_normal_sigma2_2 --------------------
//--------------------------------------------------------------------------------


DISTR_gumbelcopula2_normal_sigma2_2::DISTR_gumbelcopula2_normal_sigma2_2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "gumbel Copula - sigma2";
  pos =p;
  outpredictor = true;
  outexpectation = true;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_gumbelcopula2_normal_sigma2_2::DISTR_gumbelcopula2_normal_sigma2_2(const DISTR_gumbelcopula2_normal_sigma2_2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_gumbelcopula2_normal_sigma2_2 & DISTR_gumbelcopula2_normal_sigma2_2::operator=(
                            const DISTR_gumbelcopula2_normal_sigma2_2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_gumbelcopula2_normal_sigma2_2::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gumbelcopula2_normal_sigma2_2::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_gumbelcopula2_normal_sigma2_2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gumbelcopula2_normal_sigma2_2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[2*pos]));
  }

double DISTR_gumbelcopula2_normal_sigma2_2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

    double sigma_2 = exp((*linpred));

  double l;

  l = -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*sigma_2);

  modify_worklin();

  return l;

  }

void DISTR_gumbelcopula2_normal_sigma2_2::compute_iwls_wweightschange_weightsone(
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

  double sigma_2 = exp(*linpred);

  double nu = -0.5 + (pow(((*response)-(*worklin[0])),2))/(2*sigma_2);

  *workingweight = 0.5;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    like +=  -0.5*log(sigma_2)-pow((((*response))-(*worklin[0])),2)/(2*sigma_2);
    }

  modify_worklin();

  }

void DISTR_gumbelcopula2_normal_sigma2_2::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_gumbelcopula2_normal_sigma2_2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gumbelcopula2_normal_sigma2_2::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gumbelcopula2_normal_mu_2 ----------------
//------------------------------------------------------------------------------


DISTR_gumbelcopula2_normal_mu_2::DISTR_gumbelcopula2_normal_mu_2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  pos=p;
  family = "gumbel Copula - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";

  if(o->forceIWLS)
    {
    updateIWLS = true;
    }
  else
    {
    updateIWLS = false;
    }
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_gumbelcopula2_normal_mu_2::DISTR_gumbelcopula2_normal_mu_2(const DISTR_gumbelcopula2_normal_mu_2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
    pos = nd.pos;
    response2 = nd.response2;
    response2p = nd.response2p;
  }


const DISTR_gumbelcopula2_normal_mu_2 & DISTR_gumbelcopula2_normal_mu_2::operator=(
                            const DISTR_gumbelcopula2_normal_mu_2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_gumbelcopula2_normal_mu_2::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gumbelcopula2_normal_mu_2::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_gumbelcopula2_normal_mu_2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gumbelcopula2_normal_mu_2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2*pos + 1]);
  }

 double DISTR_gumbelcopula2_normal_mu_2::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_gumbelcopula2_normal_mu_2::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    return 0;
    }

//double DISTR_gumbelcopula2_normal_mu_2::compute_quantile_residual_mult(vector<double *> response,
//                                             vector<double *> param,
//                                             vector<double *> weight,
//                                             vector<datamatrix *> aux)
//  {
//  double u_est;
//  u_est = cdf_mult(response,param,weight,aux);
//  double res_est = randnumbers::invPhi2(u_est);
//  return res_est;
//  }

double DISTR_gumbelcopula2_normal_mu_2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;
  l = -pow((((*response))-mu),2)/(2*(*worktransformlin[0]));

  modify_worklin();

  return l;

  }


void DISTR_gumbelcopula2_normal_mu_2::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;


  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double nu = ((*response)-mu)/(*worktransformlin[0]);

  *workingweight = 1/(*worktransformlin[0]);

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    like += -pow((((*response))-mu),2)/(2*(*worktransformlin[0]));
    }

  modify_worklin();

  }


void DISTR_gumbelcopula2_normal_mu_2::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+2 * pos+1]));
  }


void DISTR_gumbelcopula2_normal_mu_2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gumbelcopula2_normal_mu_2::update(void)
  {

   unsigned i;

//  double help;

  double * worktransformlinp;
  double * workweight;

  worktransformlinp = distrp[0]->helpmat1.getV();
  workweight = workingweight.getV();

  for (i=0;i<nrobs;i++,worktransformlinp++,workweight++)
    {
        *workweight = 1/(*worktransformlinp);
    }

  }


void DISTR_gumbelcopula2_normal_mu_2::update_end(void)
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
    }

  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gumbelcopula2_rho ---------------------
//------------------------------------------------------------------------------

void DISTR_gumbelcopula2_rho::check_errors(void)
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

DISTR_gumbelcopula2_rho::DISTR_gumbelcopula2_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Gumbel Copula - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "true";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();

  }


DISTR_gumbelcopula2_rho::DISTR_gumbelcopula2_rho(const DISTR_gumbelcopula2_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_gumbelcopula2_rho & DISTR_gumbelcopula2_rho::operator=(
                            const DISTR_gumbelcopula2_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_gumbelcopula2_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gumbelcopula2_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = exp(*linpred[4]);
  *param = arg + 1;
  }

void DISTR_gumbelcopula2_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gumbelcopula2_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_gumbelcopula2_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[2] = *response[4] = first component of two dimensional reponse
   // *linpred[0] = eta_rho
   // *linpred[1] = eta_sigma_2
   // *linpred[2] = eta_mu_2
   // *linpred[3] = eta_sigma_1
   // *linpred[4] = eta_mu_1

   if (*weight[4] == 0)
     *deviance=0;
   else
     {
     double rho = exp((*linpred[4])) + 1;
     double u = randnumbers::Phi2(((*response[3]) - (*linpred[3])) / pow(exp(*linpred[2]), 0.5));
     double v = randnumbers::Phi2(((*response[0]) - (*linpred[1])) / pow(exp(*linpred[0]), 0.5));
     double logu = log(u);
     double logv = log(v);
     double logurho = pow(-logu, rho);
     double logvrho = pow(-logv, rho);
     double arg = logurho + logvrho;
     double copula = exp(-pow(arg, (1 / rho)));
     double l;

       l = log(copula) + (rho -1) * (log(-logu) + log(-logv)) - logu - logv +
                (2 / rho - 2) * log(logurho + logvrho) + log(1 + (rho - 1) * pow((logurho + logvrho), (-1 / rho)))
                 + log(randnumbers::phi(((*response[3]) - (*linpred[3])) / pow(exp(*linpred[2]), 0.5))) +
                log(randnumbers::phi(((*response[0]) - (*linpred[1])) / pow(exp(*linpred[0]), 0.5)));


    *deviance = -2*l;
    }

  }

double DISTR_gumbelcopula2_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma_2 equation
  // *worktransformlin[0] = sigma_2;
  // *worklin[1] = linear predictor of sigma_1 equation
  // *worktransformlin[1] = sigma_1;
  // *worklin[2] = linear predictor of mu_2 equation
  // *worktransformlin[2] = mu_2;
  // *worklin[3] = linear predictor of mu_1 equation
  // *worktransformlin[3] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }
    double rho = exp((*linpred)) + 1;
    double u = randnumbers::Phi2(((*response) - (*worktransformlin[3])) / pow((*worktransformlin[2]), 0.5));
    double v = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[0]), 0.5));
    double logu = log(u);
    double logv = log(v);
    double logurho = pow(-logu, rho);
    double logvrho = pow(-logv, rho);
    double arg = logurho + logvrho;
    double lcopula = (-pow(arg, (1 / rho)));
    double l;

   /*  l = log(copula) + (rho -1) * (log(-logu) + log(-logv)) - logu - logv +
                log(pow(arg, (2 / rho - 2))) + (rho -1 ) * pow(arg, (1 / rho - 2)); */

    l = lcopula + (rho -1) * (log(-logu) + log(-logv)) - logu - logv +
                (2 / rho - 2) * log(logurho + logvrho) + log(1 + (rho - 1) * pow((logurho + logvrho), (-1 / rho)));


  modify_worklin();

  return l;

  }

void DISTR_gumbelcopula2_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma_2 equation
  // *worktransformlin[0] = sigma_2;
  // *worklin[1] = linear predictor of sigma_1 equation
  // *worktransformlin[1] = sigma_1;
  // *worklin[2] = linear predictor of mu_2 equation
  // *worktransformlin[2] = mu_2;
  // *worklin[3] = linear predictor of mu_1 equation
  // *worktransformlin[3] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

    double rho = exp((*linpred)) + 1;
//    double t1 = (*worktransformlin[0]);
//    double t2 = (*worktransformlin[1]);
//    double t3 = (*worktransformlin[2]);
//    double t4 = (*worktransformlin[3]);
//    double t5 = (*response);
//    double t6 = (*response2p);
    double u = randnumbers::Phi2(((*response) - (*worktransformlin[3])) / pow((*worktransformlin[2]), 0.5));
    double v = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[0]), 0.5));
    double logu = log(u);
    double logv = log(v);
    double logurho = pow(-logu, rho);
    double logvrho = pow(-logv, rho);
    double arg = logurho + logvrho;
    double onerho = (1 / rho);
    double hilfs = (rho - 1) / pow(rho, 2);
    double hilfs2 = (2 / rho - 2);
    double lcopula = (-pow(arg, onerho));
    double grho = logurho + logvrho;
    double hilfs4 = pow(grho, onerho);
    double hilfs5 = pow(grho, (-onerho));
    double abl = logurho * (rho - 1) * log(-logu) + logvrho * (rho - 1) * log(-logv);
    double hilfs3 = abl / (grho * rho) - log(grho) * hilfs;
    double abl2 = abl + pow(rho - 1, 2) * (logurho * pow(log(-logu), 2) + logvrho * pow(log(-logv), 2));
    double zaehler = (rho - 1) * hilfs5 * (1 + (-hilfs3));
    double nenner = 1 + (rho - 1) * hilfs5;

    double nu = -hilfs4 * (hilfs3) +
                (rho - 1) * (log(-logu) + log(-logv)) + hilfs2 * abl / grho - 2 * hilfs * log(grho) + zaehler / nenner;

    *workingweight = hilfs4 * (pow(hilfs3, 2) + abl2 * onerho / grho - pow(abl / grho, 2) * onerho - 2 * abl * hilfs / grho - hilfs * log(grho) + 2 * hilfs * onerho * log(grho))
                      - (rho - 1) * (log(-logu) + log(-logv)) + 4 * hilfs * abl /grho - hilfs2 * (abl2 / grho - pow(abl / grho, 2)) +
                      2 * hilfs * log(grho) - 4 * hilfs * (rho - 1) * onerho * log(grho) -
                      (rho - 1) * hilfs5 * (1 - hilfs3) * (1 - hilfs3)/ nenner -
                      (rho - 1) * hilfs5 * (-1) * (abl2 * onerho / grho - pow(abl / grho, 2) * onerho - 2 * abl * hilfs / grho - hilfs * log(grho) + 2 * hilfs * onerho * log(grho)) / nenner
                      + (zaehler * (nenner - 1) + (rho - 1) * hilfs5 * (-hilfs3)) / pow(nenner, 2);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

      //  like += lcopula + (rho -1) * (log(-logu) + log(-logv)) - logu - logv +
       //         log(arg2 + (rho -1 ) * arg3);
          like += lcopula + (rho -1) * (log(-logu) + log(-logv)) - logu - logv +
                (2 / rho - 2) * log(logurho + logvrho) + log(1 + (rho - 1) * pow((logurho + logvrho), (-1 / rho)));
      }

  modify_worklin();

  }

void DISTR_gumbelcopula2_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 1 - 1 / exp((*linpred[predstart_mumult+4]));
  }


void DISTR_gumbelcopula2_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): shifted exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gumbelcopula2_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin) + 1;
    }

  }



//--------------------------------------------------------------------------------
//----------------- CLASS: DISTR_gaussiancopula_dagum_p --------------------
//--------------------------------------------------------------------------------


DISTR_gaussiancopula_dagum_p::DISTR_gaussiancopula_dagum_p(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,6,w)
  {
  family = "Gaussian Copula - Dagum p";
  pos =p;
  outpredictor = true;
  outexpectation = true;
  predictor_name = "p";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_gaussiancopula_dagum_p::DISTR_gaussiancopula_dagum_p(const DISTR_gaussiancopula_dagum_p & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_gaussiancopula_dagum_p & DISTR_gaussiancopula_dagum_p::operator=(
                            const DISTR_gaussiancopula_dagum_p & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_gaussiancopula_dagum_p::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gaussiancopula_dagum_p::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_gaussiancopula_dagum_p::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_dagum_p::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[3*pos]);
  }

double DISTR_gaussiancopula_dagum_p::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

    double p = exp((*linpred));
    double bcurrent = (*worktransformlin[2]);
    double acurrent = (*worktransformlin[4]);
    double respdivb = (*response) / bcurrent;
    double hilfs = pow(respdivb,acurrent);
    double hilfs2 = pow(respdivb, -acurrent);


    double u = pow((1 + hilfs2), -p);
    double v = pow((1 + pow(((*response2p) / (*worktransformlin[3])), -(*worktransformlin[5]))), -(*worktransformlin[1]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = randnumbers::invPhi2(v);
    double orho = 1 - pow((*worktransformlin[0]), 2);
    double l;

    l = (*worktransformlin[0]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[0]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho +
                log(p) + acurrent*p*log((*response)) - acurrent*p*log(bcurrent)-p*log(1+hilfs);

  modify_worklin();

  return l;

  }

void DISTR_gaussiancopula_dagum_p::compute_iwls_wweightschange_weightsone(
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

    double p = exp((*linpred));
    double bcurrent = (*worktransformlin[2]);
    double acurrent = (*worktransformlin[4]);
    double respdivb = (*response) / bcurrent;
    double hilfs = pow(respdivb,acurrent);
    double hilfs2 = pow(respdivb, -acurrent);


    double u = pow((1 + hilfs2), -p);
    double v = pow((1 + pow(((*response2p) / (*worktransformlin[3])), -(*worktransformlin[5]))), -(*worktransformlin[1]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = randnumbers::invPhi2(v);
    double orho = 1 - pow((*worktransformlin[0]), 2);

    double dphiinvu = pow(2*PI, 0.5) / exp(-0.5*pow(phinvu, 2));
    double ddphiinvu = phinvu * 2 * PI / pow(exp(-0.5*pow(phinvu, 2)), 2);


    double d1 =  -p * pow(1 +  hilfs2, -p) * log(1 +  hilfs2);
//    double d2 = d1 - +p * p * pow(1 +  hilfs2, -p) * pow(log(1 +  hilfs2), 2);

    double nu = 1 + acurrent*p*log((*response)) - acurrent*p*log(bcurrent)
                -p*log(1+hilfs) +
                ((*worktransformlin[0]) * d1 * dphiinvu / orho) * (phinvv - (*worktransformlin[0]) * phinvu);

    *workingweight = 1  -
                    ((*worktransformlin[0]) * d1 * d1 * ddphiinvu / orho) * (phinvv - (*worktransformlin[0]) * phinvu) +
                    (*worktransformlin[0]) * (*worktransformlin[0]) * pow(d1 * dphiinvu, 2) / orho;

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  (*worktransformlin[0]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[0]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho +
                log(p) + acurrent*p*log((*response)) - acurrent*p*log(bcurrent)-p*log(1+hilfs);

      }


  modify_worklin();

  }

void DISTR_gaussiancopula_dagum_p::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_gaussiancopula_dagum_p::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (p): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_dagum_p::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }




//--------------------------------------------------------------------------------
//----------------- CLASS: DISTR_gaussiancopula_dagum_b --------------------
//--------------------------------------------------------------------------------


DISTR_gaussiancopula_dagum_b::DISTR_gaussiancopula_dagum_b(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,6,w)
  {
  family = "Gaussian Copula - Dagum b";
  pos =p;
  outpredictor = true;
  outexpectation = true;
  predictor_name = "b";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_gaussiancopula_dagum_b::DISTR_gaussiancopula_dagum_b(const DISTR_gaussiancopula_dagum_b & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_gaussiancopula_dagum_b & DISTR_gaussiancopula_dagum_b::operator=(
                            const DISTR_gaussiancopula_dagum_b & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_gaussiancopula_dagum_b::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gaussiancopula_dagum_b::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_gaussiancopula_dagum_b::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_dagum_b::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[3*pos + 1]));
  }

double DISTR_gaussiancopula_dagum_b::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

    double b = exp((*linpred));
    double respdivb = (*response) / b;
    double acurrent = (*worktransformlin[4]);
    double hilfs = pow(respdivb,acurrent);
//    double hilfs2 = pow(respdivb, -acurrent);
    double pcurrent = (*worktransformlin[2]);

    double u = pow((1 + pow((respdivb), -acurrent)), -pcurrent);
    double v = pow((1 + pow(((*response2p) / (*worktransformlin[1])), -(*worktransformlin[5]))), -(*worktransformlin[3]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = randnumbers::invPhi2(v);
    double orho = 1 - pow((*worktransformlin[0]), 2);

    double l;

    l = (*worktransformlin[0]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[0]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho
                - acurrent*pcurrent*log(b) - (pcurrent+1)*log(1+hilfs);

  modify_worklin();

  return l;

  }

void DISTR_gaussiancopula_dagum_b::compute_iwls_wweightschange_weightsone(
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

    double b = exp((*linpred));
    double respdivb = (*response) / b;
    double acurrent = (*worktransformlin[4]);
    double hilfs = pow(respdivb,acurrent);
    double hilfs2 = pow(respdivb, -acurrent);
    double pcurrent = (*worktransformlin[2]);

    double u = pow((1 + pow((respdivb), -acurrent)), -pcurrent);
    double v = pow((1 + pow(((*response2p) / (*worktransformlin[1])), -(*worktransformlin[5]))), -(*worktransformlin[3]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = randnumbers::invPhi2(v);
    double orho = 1 - pow((*worktransformlin[0]), 2);

    double dphiinvu = pow(2*PI, 0.5) / exp(-0.5*pow(phinvu, 2));
    double ddphiinvu = phinvu * 2 * PI / pow(exp(-0.5*pow(phinvu, 2)), 2);


    double d1 =  -acurrent * pcurrent * pow(1 +  hilfs2, -pcurrent-1) * hilfs2;
//    double d2 = -acurrent * acurrent * pow(1+ hilfs2, -pcurrent-1) * pcurrent * hilfs2 +
//                acurrent * acurrent * pow(1+ hilfs2, -pcurrent-2) * (pcurrent+1) * pow(hilfs2, 2);

    double nu = acurrent - ((pcurrent+1)*acurrent)/(1+hilfs)+
                ((*worktransformlin[0]) * d1 * dphiinvu / orho) * (phinvv - (*worktransformlin[0]) * phinvu);

    *workingweight = (pow(acurrent,2)*pcurrent)/(pcurrent+2)-
                    ((*worktransformlin[0]) * d1 * d1 * ddphiinvu / orho) * (phinvv - (*worktransformlin[0]) * phinvu)
                    + (*worktransformlin[0]) * (*worktransformlin[0]) * pow(d1 * dphiinvu, 2) / orho;

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  (*worktransformlin[0]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[0]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho
                - acurrent*pcurrent*log(b) - (pcurrent+1)*log(1+hilfs);

      }



  modify_worklin();

  }

void DISTR_gaussiancopula_dagum_b::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_gaussiancopula_dagum_b::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (b): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_dagum_b::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_gaussiancopula_dagum_a ------------------
//------------------------------------------------------------------------------


DISTR_gaussiancopula_dagum_a::DISTR_gaussiancopula_dagum_a(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,6,w)
  {
  pos=p;
  family = "Gaussian Copula - Dagum a";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "a";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_gaussiancopula_dagum_a::DISTR_gaussiancopula_dagum_a(const DISTR_gaussiancopula_dagum_a & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
    pos = nd.pos;
    response2 = nd.response2;
    response2p = nd.response2p;
  }


const DISTR_gaussiancopula_dagum_a & DISTR_gaussiancopula_dagum_a::operator=(
                            const DISTR_gaussiancopula_dagum_a & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_gaussiancopula_dagum_a::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gaussiancopula_dagum_a::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_gaussiancopula_dagum_a::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_dagum_a::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[3 * pos + 2]);
  }

 double DISTR_gaussiancopula_dagum_a::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_gaussiancopula_dagum_a::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    return 0;
    }

//double DISTR_gaussiancopula_dagum_a::compute_quantile_residual_mult(vector<double *> response,
//                                             vector<double *> param,
//                                             vector<double *> weight,
//                                             vector<datamatrix *> aux)
//  {
//  double u_est;
//  u_est = cdf_mult(response,param,weight,aux);
//  double res_est = randnumbers::invPhi2(u_est);
//  return res_est;
//  }

double DISTR_gaussiancopula_dagum_a::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

    double a = exp((*linpred));
    double respdivb = (*response) / (*worktransformlin[4]);
    double hilfs = pow(respdivb,a);
    double pcurrent = (*worktransformlin[2]);

    double u = pow((1 + pow((respdivb), -a)), -pcurrent);
    double v = pow((1 + pow(((*response2p) / (*worktransformlin[5])), -(*worktransformlin[1]))), -(*worktransformlin[3]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = randnumbers::invPhi2(v);
    double orho = 1 - pow((*worktransformlin[0]), 2);
    double l;

    l = (*worktransformlin[0]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[0]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho +
                log(a) +(a*pcurrent)*log((*response)) - a*pcurrent*log((*worktransformlin[4]))
                - (pcurrent+1)*log(1+hilfs);


    modify_worklin();

    return l;

  }


void DISTR_gaussiancopula_dagum_a::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;


  if (counter==0)
    {
    set_worklin();
    }

    double a = exp((*linpred));
    double respdivb = (*response) / (*worktransformlin[4]);
    double hilfs = pow(respdivb,a);
    double pcurrent = (*worktransformlin[2]);

    double u = pow((1 + pow((respdivb), -a)), -pcurrent);
    double v = pow((1 + pow(((*response2p) / (*worktransformlin[5])), -(*worktransformlin[1]))), -(*worktransformlin[3]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = randnumbers::invPhi2(v);
    double orho = 1 - pow((*worktransformlin[0]), 2);

    double dphiinvu = pow(2*PI, 0.5) / exp(-0.5*pow(phinvu, 2));
    double ddphiinvu = phinvu * 2 * PI / pow(exp(-0.5*pow(phinvu, 2)), 2);


    double d1 =  pcurrent * pow(1 +  pow(respdivb, -a), -pcurrent-1) * a * pow(respdivb, -a) * log(respdivb);
//    double d2 = d1 - pcurrent * pow(1+pow(respdivb, -a), -pcurrent-1) * a * a * log(respdivb) * log(respdivb) * pow(respdivb, -a) +
//                pcurrent * (pcurrent + 1) *pow(1+pow(respdivb, -a), -pcurrent-2) * pow((pow(respdivb, -a) * log(respdivb) * a), 2);

    double nu = 1 + a*pcurrent*log(respdivb)
                - ((pcurrent+1)*a*hilfs*log(respdivb))/(1+hilfs) +
                ((*worktransformlin[0]) * d1 * dphiinvu / orho) * (phinvv - (*worktransformlin[0]) * phinvu);

    *workingweight = 1 + ((pcurrent+1)*pow(a,2)*hilfs*pow(log(respdivb),2))/pow((1+hilfs),2) -
                    ((*worktransformlin[0]) * d1 * d1 * ddphiinvu / orho) * (phinvv - (*worktransformlin[0]) * phinvu) + (*worktransformlin[0]) * (*worktransformlin[0]) * pow(d1 * dphiinvu, 2) / orho;

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  (*worktransformlin[0]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[0]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho +
                log(a) +(a*pcurrent)*log((*response)) - a*pcurrent*log((*worktransformlin[4]))
                - (pcurrent+1)*log(1+hilfs);

      }


  modify_worklin();

  }


void DISTR_gaussiancopula_dagum_a::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+2 * pos+1]));
  }


void DISTR_gaussiancopula_dagum_a::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (a): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


/*void DISTR_gaussiancopula_dagum_a::update(void)
  {

   unsigned i;

  double help;

  double * worktransformlinp;
  double * workweight;

  worktransformlinp = distrp[0]->helpmat1.getV();
  workweight = workingweight.getV();

  for (i=0;i<nrobs;i++,worktransformlinp++,workweight++)
    {
        *workweight = 1/(*worktransformlinp);
    }

  }*/


void DISTR_gaussiancopula_dagum_a::update_end(void)
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
//    double t = *pmu;
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gaussiancopula_dagum_rho ----------------
//------------------------------------------------------------------------------

void DISTR_gaussiancopula_dagum_rho::check_errors(void)
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
            if (*workresp <= 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative/zero response values encountered\n");
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

DISTR_gaussiancopula_dagum_rho::DISTR_gaussiancopula_dagum_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,6,w)
  {
  family = "Gaussian Copula - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
    linpredminlimit=-100;
  linpredmaxlimit=100;
  check_errors();
  }


DISTR_gaussiancopula_dagum_rho::DISTR_gaussiancopula_dagum_rho(const DISTR_gaussiancopula_dagum_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_gaussiancopula_dagum_rho & DISTR_gaussiancopula_dagum_rho::operator=(
                            const DISTR_gaussiancopula_dagum_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_gaussiancopula_dagum_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_dagum_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[6]) / pow(1+ pow((*linpred[6]), 2), 0.5);
  }

void DISTR_gaussiancopula_dagum_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gaussiancopula_dagum_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_gaussiancopula_dagum_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {


   if (*weight[6] == 0)
     *deviance=0;
   else
     {
     double rho = ((*linpred[6]))/pow((1 + pow((*linpred[6]), 2)), 0.5);
     double a1 = exp((*linpred[5]));
     double b1 = exp((*linpred[4]));
     double p1 = exp((*linpred[3]));
     double a2 = exp((*linpred[2]));
     double b2 = exp((*linpred[1]));
     double p2 = exp((*linpred[0]));
     double u = pow((1 + pow(((*response[5]) / b1), -a1)), -p1);
     double v = pow((1 + pow(((*response[1]) / b2), -a2)), -p2);
     double lfy1 = log(a1 * p1 * pow(((*response[5]) / b1), a1 * p1) / ((*response[5]) * pow((pow(((*response[5]) / b1), a1) + 1), p1 + 1)));
     double lfy2 = log(a2 * p2 * pow(((*response[1]) / b2), a2 * p2) / ((*response[1]) * pow((pow(((*response[1]) / b2), a2) + 1), p2 + 1)));
     double phinvu = randnumbers::invPhi2(u);
     double phinvv = randnumbers::invPhi2(v);
     double orho = 1 - pow(rho, 2);
     double l;

       l = - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho
                 + lfy1 + lfy2;


    *deviance = -2*l;
    }

  }

double DISTR_gaussiancopula_dagum_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {



  if (counter==0)
    {
    set_worklin();
    }
    double u = pow((1 + pow(((*response) / (*worktransformlin[4])), -(*worktransformlin[5]))), -(*worktransformlin[3]));
    double v = pow((1 + pow(((*response2p) / (*worktransformlin[1])), -(*worktransformlin[2]))), -(*worktransformlin[0]));
    double hilfs = pow(1 + pow((*linpred), 2), 0.5);
    double rho = (*linpred) / hilfs;
    if (*linpred <= -100)
        rho  = -0.99995;
    else if (*linpred >= 100)
        rho  = 0.99995;

    double orho = 1 - pow(rho, 2);
    double phinvu = randnumbers::invPhi2(u);
    double phinvv = randnumbers::invPhi2(v);
    double l;

    l = - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;


  modify_worklin();

  return l;

  }

void DISTR_gaussiancopula_dagum_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }

    double u = pow((1 + pow(((*response) / (*worktransformlin[4])), -(*worktransformlin[5]))), -(*worktransformlin[3]));
    double v = pow((1 + pow(((*response2p) / (*worktransformlin[1])), -(*worktransformlin[2]))), -(*worktransformlin[0]));
    double hilfs = pow(1 + pow((*linpred), 2), 0.5);
    double rho = (*linpred) / hilfs;
    if (*linpred <= -100)
        rho  = -0.99995;
    else if (*linpred >= 100)
        rho  = 0.99995;

    double orho = 1 - pow(rho, 2);
    double phinvu = randnumbers::invPhi2(u);
    double phinvv = randnumbers::invPhi2(v);

    double nu = rho * pow(orho, 0.5) + (hilfs + rho * (*linpred)) * (phinvu * phinvv) - (*linpred) * (pow(phinvu, 2) + pow(phinvv, 2));

  //  *workingweight = (orho * pow(rho, 2) - pow(orho, 2)) + (pow(phinvu, 2) + pow(phinvv, 2))
   //                 - (2 * rho + rho * orho) * (phinvu * phinvv);
    *workingweight = 1 - pow(rho, 4);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

          like += - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;
      }

  modify_worklin();

  }

void DISTR_gaussiancopula_dagum_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 2 * std::asin((*linpred[predstart_mumult+6]) / (pow(1 + pow((*linpred[predstart_mumult+6]), 2), 0.5))) / PI ;
  }


void DISTR_gaussiancopula_dagum_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): \n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_dagum_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin) / pow(1 + pow((*worklin), 2), 0.5);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gumbelcopula_rho ----------------------
//------------------------------------------------------------------------------

void DISTR_gumbelcopula_rho::check_errors(void)
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

            if((*workresp) > 1 )
            {
                errors=true;
                errormessages.push_back("ERROR: cdfs of marginals take values inbetween zero and one!\n");
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

DISTR_gumbelcopula_rho::DISTR_gumbelcopula_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Gumbel Copula - rho";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();

  }


DISTR_gumbelcopula_rho::DISTR_gumbelcopula_rho(const DISTR_gumbelcopula_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_gumbelcopula_rho & DISTR_gumbelcopula_rho::operator=(
                            const DISTR_gumbelcopula_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_gumbelcopula_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gumbelcopula_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = exp(*linpred[2]);
  *param = arg + 1;
  }

void DISTR_gumbelcopula_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gumbelcopula_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_gumbelcopula_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double rho = exp((*linpred[2])) + 1;
     double logu = log((*response[1]));
     double logv = log((*response[0]));
     double logurho = pow(-logu, rho);
     double logvrho = pow(-logv, rho);
     double arg = logurho + logvrho;
     double copula = exp(-pow(arg, (1 / rho)));
     double l;

       l = log(copula) + (rho -1) * (log(-logu) + log(-logv)) - logu - logv +
                (2 / rho - 2) * log(logurho + logvrho) + log(1 + (rho - 1) * pow((logurho + logvrho), (-1 / rho)));


    *deviance = -2*l;
    }

  }

double DISTR_gumbelcopula_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma_2 equation
  // *worktransformlin[0] = sigma_2;
  // *worklin[1] = linear predictor of sigma_1 equation
  // *worktransformlin[1] = sigma_1;
  // *worklin[2] = linear predictor of mu_2 equation
  // *worktransformlin[2] = mu_2;
  // *worklin[3] = linear predictor of mu_1 equation
  // *worktransformlin[3] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }
    double rho = exp((*linpred)) + 1;
    double logu = log((*response));
    double logv = log((*response2p));
    double logurho = pow(-logu, rho);
    double logvrho = pow(-logv, rho);
    double arg = logurho + logvrho;
    double copula = exp(-pow(arg, (1 / rho)));
    double l;

   /*  l = log(copula) + (rho -1) * (log(-logu) + log(-logv)) - logu - logv +
                log(pow(arg, (2 / rho - 2))) + (rho -1 ) * pow(arg, (1 / rho - 2)); */

    l = log(copula) + (rho -1) * (log(-logu) + log(-logv)) - logu - logv +
                (2 / rho - 2) * log(logurho + logvrho) + log(1 + (rho - 1) * pow((logurho + logvrho), (-1 / rho)));


  modify_worklin();

  return l;

  }

void DISTR_gumbelcopula_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }
  /*  double rho = exp((*linpred)) + 1;
    double logu = log((*response));
    double logv = log((*response2p));
    double logurho = pow(-logu, rho);
    double logvrho = pow(-logv, rho);
    double arg = logurho + logvrho;
    double onerho = (1 / rho);
    double lcopula = (-pow(arg, onerho));
    double hilfs = (1 / rho - 2);
    double hilfs2 = (2 / rho - 2);
    double arg2 = pow(arg, hilfs2);
    double arg3 = pow(arg, hilfs);
    double divi = (rho - 1) / rho;
    double abl = logurho * (rho - 1) * log(-logu) + logvrho * (rho - 1) * log(-logv);
    double abl2 = abl2 = abl + pow(rho - 1, 2) * (logurho * pow(log(-logu), 2) + logvrho * pow(log(-logv), 2));
    double nenner = arg2 + (rho -1) *  arg3;
    double zaehler1 = arg2 * (abl * (2 / rho - 2) / arg + log(arg) * (-2 * (rho - 1) / pow(rho, 2)));
    double zaehler2 = (rho - 1) * arg3 + (rho - 1) * arg3 * (abl * hilfs / arg + log(arg) * (-(rho - 1)/pow(rho, 2)));
    double zaehler = zaehler1 + zaehler2;

    double nu = -lcopula * (log(arg) * (-(rho - 1) / pow(rho, 2)) + (1 / rho) * abl / arg) +
                (rho - 1) * (log(-logu) + log(-logv)) +
                (zaehler1 + zaehler2) / (arg2 + (rho -1) * arg3);

    *workingweight =  lcopula * pow((abl * onerho / arg - (rho - 1) * log(arg) / pow(rho, 2)), 2) +
                      lcopula * (abl2 * onerho / arg - pow(abl / arg, 2) * onerho - 2 * abl * (rho - 1) / (arg * pow(rho, 2)) + log(arg) * (rho -1) * (-1 + 2 * (rho - 1) / rho) / pow(rho, 2)) -
                      (rho - 1) * (log(-logu) + log(-logv)) -
                      pow(arg2 * (hilfs2 * abl / arg + log(arg) * (-2 * (rho - 1) / pow(rho, 2))), 2) / nenner -
                      - (arg2 * (-4 * divi * abl / (arg * rho) + hilfs2 * abl2 / arg - hilfs2 * pow(abl / arg, 2) - 2 * divi * log(arg) / rho + 3 * pow(divi, 2) * log(arg) / rho)) / nenner
                      - zaehler2 / nenner
                      - ((rho - 1) * arg3 * (abl * hilfs / arg - divi * log(arg) / rho) * (1 + abl * hilfs / arg - divi * log(arg) / rho)) / nenner
                      - (arg3 * (abl2 * hilfs / arg - hilfs * pow(abl / arg, 2) - 2 * divi * abl / (arg * rho) - log(arg) * divi / rho + 2 *  log(arg) * pow(divi, 2) / rho)) / nenner
                      + pow(zaehler / nenner, 2); */
    double rho = exp((*linpred)) + 1;
    double logu = log((*response));
    double logv = log((*response2p));
    double logurho = pow(-logu, rho);
    double logvrho = pow(-logv, rho);
    double arg = logurho + logvrho;
    double onerho = (1 / rho);
    double hilfs = (rho - 1) / pow(rho, 2);
    double hilfs2 = (2 / rho - 2);
    double lcopula = (-pow(arg, onerho));
    double grho = logurho + logvrho;
    double hilfs4 = pow(grho, onerho);
    double hilfs5 = pow(grho, (-onerho));
    double abl = logurho * (rho - 1) * log(-logu) + logvrho * (rho - 1) * log(-logv);
    double hilfs3 = abl / (grho * rho) - log(grho) * hilfs;
    double abl2 = abl + pow(rho - 1, 2) * (logurho * pow(log(-logu), 2) + logvrho * pow(log(-logv), 2));
    double zaehler = (rho - 1) * hilfs5 * (1 + (-hilfs3));
    double nenner = 1 + (rho - 1) * hilfs5;

    double nu = -hilfs4 * (hilfs3) +
                (rho - 1) * (log(-logu) + log(-logv)) + hilfs2 * abl / grho - 2 * hilfs * log(grho) + zaehler / nenner;

    *workingweight = hilfs4 * (pow(hilfs3, 2) + abl2 * onerho / grho - pow(abl / grho, 2) * onerho - 2 * abl * hilfs / grho - hilfs * log(grho) + 2 * hilfs * onerho * log(grho))
                      - (rho - 1) * (log(-logu) + log(-logv)) + 3 * hilfs * abl /grho - hilfs2 * (abl2 / grho - pow(abl / grho, 2)) +
                      2 * hilfs * log(grho) - 2 * hilfs * onerho * log(grho) -
                      (rho - 1) * hilfs5 * (1 - hilfs3) * (1 - hilfs3)/ nenner -
                      (rho - 1) * hilfs5 * (-1) * (abl2 * onerho / grho - pow(abl / grho, 2) * onerho - 2 * abl * hilfs / grho - hilfs * log(grho) + 2 * hilfs * onerho * log(grho)) / nenner
                      + (zaehler * (nenner - 1) + (rho - 1) * hilfs5 * (-hilfs3)) / pow(nenner, 2);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

      //  like += lcopula + (rho -1) * (log(-logu) + log(-logv)) - logu - logv +
       //         log(arg2 + (rho -1 ) * arg3);
          like += lcopula + (rho -1) * (log(-logu) + log(-logv)) - logu - logv +
                (2 / rho - 2) * log(logurho + logvrho) + log(1 + (rho - 1) * pow((logurho + logvrho), (-1 / rho)));
      }

  modify_worklin();

  }

void DISTR_gumbelcopula_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 1 - 1 / exp((*linpred[predstart_mumult+2]));
  }


void DISTR_gumbelcopula_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): shifted exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gumbelcopula_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin) + 1;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gaussiancopula_rho --------------------
//------------------------------------------------------------------------------
void DISTR_gaussiancopula_rho::check_errors(void)
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

            if((*workresp) > 1 )
            {
                errors=true;
                errormessages.push_back("ERROR: cdfs of marginals take values inbetween zero and one!\n");
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


DISTR_gaussiancopula_rho::DISTR_gaussiancopula_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Gaussian Copula - rho";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
  linpredminlimit=-100;
  linpredmaxlimit=100;
  check_errors();
  }


DISTR_gaussiancopula_rho::DISTR_gaussiancopula_rho(const DISTR_gaussiancopula_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_gaussiancopula_rho & DISTR_gaussiancopula_rho::operator=(
                            const DISTR_gaussiancopula_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_gaussiancopula_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2]) / pow(1+ pow((*linpred[2]), 2), 0.5);
  }

void DISTR_gaussiancopula_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gaussiancopula_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_gaussiancopula_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {



   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double rho = (*linpred[2]) / pow(1 + pow((*linpred[2]), 2), 0.5);
     if (*linpred[0] <= -100)
        rho  = -0.99995;
     else if (*linpred[0] >= 100)
        rho  = 0.99995;

     double orho = 1 - pow(rho, 2);
     double phinvu = randnumbers::invPhi2((*response[1]));
     double phinvv = randnumbers::invPhi2((*response[0]));
     double l;

      l = - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;


    *deviance = -2*l;
    }

  }

double DISTR_gaussiancopula_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {


  if (counter==0)
    {
    set_worklin();
    }
    double rho = (*linpred) / pow(1 + pow((*linpred), 2), 0.5);
    if (*linpred <= -100)
        rho  = -0.99995;
    else if (*linpred >= 100)
        rho  = 0.99995;

    double orho = 1 - pow(rho, 2);
    double phinvu = randnumbers::invPhi2((*response));
    double phinvv = randnumbers::invPhi2((*response2p));
    double l;


    l = - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;


  modify_worklin();

  return l;

  }

void DISTR_gaussiancopula_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }
    double hilfs = pow(1 + pow((*linpred), 2), 0.5);
    double rho = (*linpred) / hilfs;
    if (*linpred <= -100)
        rho  = -0.99995;
    else if (*linpred >= 100)
        rho  = 0.99995;

    double orho = 1 - pow(rho, 2);
    double phinvu = randnumbers::invPhi2((*response));
    double phinvv = randnumbers::invPhi2((*response2p));

    double nu = rho * pow(orho, 0.5) + (hilfs + rho * (*linpred)) * (phinvu * phinvv) - (*linpred) * (pow(phinvu, 2) + pow(phinvv, 2));

  //  *workingweight = (orho * pow(rho, 2) - pow(orho, 2)) + (pow(phinvu, 2) + pow(phinvv, 2))
   //                 - (2 * rho + rho * orho) * (phinvu * phinvv);
    *workingweight = 1 - pow(rho, 4);

  //  if((*workingweight) <= 0)
   //     *workingweight = 0.001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

          like += - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;
      }

  modify_worklin();

  }

void DISTR_gaussiancopula_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 2 * std::asin((*linpred[predstart_mumult+2]) / (pow(1 + pow((*linpred[predstart_mumult+2]), 2), 0.5))) / PI ;
  }


void DISTR_gaussiancopula_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): \n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_rho::update_end(void)
  {

  // helpmat1 stores rho

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin) / pow(1 + pow((*worklin), 2), 0.5);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_gaussiancopula_rhofz ------------------
//------------------------------------------------------------------------------
void DISTR_gaussiancopula_rhofz::check_errors(void)
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

            if((*workresp) > 1 )
            {
                errors=true;
                errormessages.push_back("ERROR: cdfs of marginals take values inbetween zero and one!\n");
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


DISTR_gaussiancopula_rhofz::DISTR_gaussiancopula_rhofz(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Gaussian Copula - rho";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
    linpredminlimit=-20.6;
  linpredmaxlimit=20.6;
  check_errors();
  }


DISTR_gaussiancopula_rhofz::DISTR_gaussiancopula_rhofz(const DISTR_gaussiancopula_rhofz & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_gaussiancopula_rhofz & DISTR_gaussiancopula_rhofz::operator=(
                            const DISTR_gaussiancopula_rhofz & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_gaussiancopula_rhofz::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_rhofz::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (exp(2 * (*linpred[2])) - 1) / (exp(2 * (*linpred[2])) + 1);
  }

void DISTR_gaussiancopula_rhofz::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gaussiancopula_rhofz::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_gaussiancopula_rhofz::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {



   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double rho;

 /* if (*linpred <= -20.6)
    rho  = -0.99995;
  else if (*linpred >= 20.6)
   rho  = 0.99995;
  else*/
    rho = (exp(2 * (*linpred[2])) - 1) / (exp(2 * (*linpred[2])) + 1);

     double orho = 1 - pow(rho, 2);
     double phinvu = randnumbers::invPhi2((*response[1]));
     double phinvv = randnumbers::invPhi2((*response[0]));
     double l;

      l = - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;


    *deviance = -2*l;
    }

  }

double DISTR_gaussiancopula_rhofz::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {



  if (counter==0)
    {
    set_worklin();
    }
    double rho;

 /* if (*linpred <= -20.6)
    rho  = -0.99995;
  else if (*linpred >= 20.6)
   rho  = 0.99995;
  else*/
    rho = (exp(2 * (*linpred)) - 1) / (exp(2 * (*linpred)) + 1);

    double orho = 1 - pow(rho, 2);
    double phinvu = randnumbers::invPhi2((*response));
    double phinvv = randnumbers::invPhi2((*response2p));
    double l;


    l = - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;


  modify_worklin();

  return l;

  }

void DISTR_gaussiancopula_rhofz::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }
    double rho;
   /* if (*linpred <= -20.6)
    rho  = -0.99995;
  else if (*linpred >= 20.6)
   rho  = 0.99995;
  else*/
    rho = (exp(2 * (*linpred)) - 1) / (exp(2 * (*linpred)) + 1);


  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;

    double orho = 1 - pow(rho, 2);
    double phinvu = randnumbers::invPhi2((*response));
    double phinvv = randnumbers::invPhi2((*response2p));

    double nu = rho - (rho / orho) * (pow(phinvu, 2) + pow(phinvv, 2)) + 2 * (1 / oneminusrho2 - 0.5) * (phinvu * phinvv) ;

  //  *workingweight = (orho * pow(rho, 2) - pow(orho, 2)) + (pow(phinvu, 2) + pow(phinvv, 2))
   //                 - (2 * rho + rho * orho) * (phinvu * phinvv);
    *workingweight = (rho2) + 1;

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

          like += - 0.5 * log(orho) + rho * phinvu * phinvv / orho - 0.5 * pow(rho, 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho;
      }

  modify_worklin();

  }

void DISTR_gaussiancopula_rhofz::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 2 * std::asin((exp( 2 * (*linpred[predstart_mumult+2])) - 1) / (exp( 2 * (*linpred[predstart_mumult+2])) - 1)) / PI ;
  }


void DISTR_gaussiancopula_rhofz::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): Fisher z-transformation\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_rhofz::update_end(void)
  {

  // helpmat1 stores rho

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (exp(2 * (*worklin)) - 1) / (exp(2 * (*worklin)) + 1);
    }

  }

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_tcopula_df -------------------------------
//------------------------------------------------------------------------------


DISTR_tcopula_df::DISTR_tcopula_df(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "t copula - Degrees of Freedom";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "df";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_tcopula_df::DISTR_tcopula_df(const DISTR_tcopula_df & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_tcopula_df & DISTR_tcopula_df::operator=(
                            const DISTR_tcopula_df & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_tcopula_df::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_tcopula_df::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[1]));
  }

void DISTR_tcopula_df::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_tcopula_df::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }


double DISTR_tcopula_df::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }

  double degf = exp((*linpred));

  double oneminusrho2 = 1- pow((*worktransformlin[0]),2);
  double np2d2 = (degf+2)/2;
  double nd2 = degf/2;
  double X_1 = (*response)-(*worklin[4]);
  double X_2 = (*response2p)-(*worklin[3]);
  double l;


     l = randnumbers::lngamma_exact(np2d2)-randnumbers::lngamma_exact(nd2)-log(degf)
     -np2d2*log( 1 + (1/(degf*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[0])*X_1*X_2+pow(X_2,2) ) );


  modify_worklin();

  return l;

  }

void DISTR_tcopula_df::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double degf = exp((*linpred));

  double oneminusrho2 = 1- pow((*worktransformlin[0]),2);
  double np2d2 = (degf+2)/2;
  double nd2 = degf/2;
  double X_1 = (*response);
  double X_2 = (*response2p);
  double nenner_C = 1 + (1/(degf*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[0])*X_1*X_2+pow(X_2,2) );

    double nu = nd2*( randnumbers::digamma_exact(np2d2)-randnumbers::digamma_exact(nd2)-log(nenner_C) ) - 1 + np2d2*(nenner_C-1)/nenner_C;

    *workingweight =  - pow(nd2,2)*( randnumbers::trigamma_exact(np2d2) - randnumbers::trigamma_exact(nd2) ) - degf/(degf+2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  randnumbers::lngamma_exact(np2d2)-randnumbers::lngamma_exact(nd2)-log(degf)
     -np2d2*log( nenner_C );

      }

  modify_worklin();

  }


void DISTR_tcopula_df::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (df): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_tcopula_df::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_tcopula_rho ------------------------------
//------------------------------------------------------------------------------


DISTR_tcopula_rho::DISTR_tcopula_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "t copula - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
    linpredminlimit=-100;
  linpredmaxlimit=100;

  }


DISTR_tcopula_rho::DISTR_tcopula_rho(const DISTR_tcopula_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_tcopula_rho & DISTR_tcopula_rho::operator=(
                            const DISTR_tcopula_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_tcopula_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_tcopula_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = (*linpred[0]);
  *param = arg/pow(1+pow(arg,2),0.5);
  }

void DISTR_tcopula_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_tcopula_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_tcopula_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

 // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }


  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double X_1 = (*response);
  double X_2 = (*response2p);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*rho*X_1*X_2+pow(X_2,2) );
  double l;


     l = -0.5*log(oneminusrho2)-(((*worktransformlin[0])+2)/2)*log( nenner_C );



  modify_worklin();

  return l;

  }

void DISTR_tcopula_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }

  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double X_1 = (*response);
  double X_2 = (*response2p);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*rho*X_1*X_2+pow(X_2,2) );

  double nu = oneminusrho2*(*linpred) - (((*worktransformlin[0])+2)/((*worktransformlin[0])*nenner_C))*( (*linpred)*(pow(X_1,2)+pow(X_2,2)) -
                                                                                                      (pow(1/oneminusrho2,0.5)+rho*(*linpred))*X_1*X_2 );



    *workingweight = 1-pow(rho2,2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(oneminusrho2)-(((*worktransformlin[0])+2)/2)*log( nenner_C );

      }

  modify_worklin();

  }

void DISTR_tcopula_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 2 * std::asin((*linpred[predstart_mumult+1])) / PI;
  }

void DISTR_tcopula_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho): fisher z-transformation\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_tcopula_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin)/pow((1+pow((*worklin),2)),0.5);
    }

  }


//--------------------------------------------------------------------------------
//----------------- CLASS: DISTR_frankcopula2_normal_sigma2 --------------------
//--------------------------------------------------------------------------------


DISTR_frankcopula2_normal_sigma2::DISTR_frankcopula2_normal_sigma2(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "frank Copula - sigma2";
  pos =p;
  outpredictor = true;
  outexpectation = true;
  predictor_name = "sigma2";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_frankcopula2_normal_sigma2::DISTR_frankcopula2_normal_sigma2(const DISTR_frankcopula2_normal_sigma2 & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_frankcopula2_normal_sigma2 & DISTR_frankcopula2_normal_sigma2::operator=(
                            const DISTR_frankcopula2_normal_sigma2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_frankcopula2_normal_sigma2::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_frankcopula2_normal_sigma2::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_frankcopula2_normal_sigma2::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_frankcopula2_normal_sigma2::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = pow(exp((*linpred[pos])),0.5);
  }

double DISTR_frankcopula2_normal_sigma2::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = (eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

  double sigma_2 = exp((*linpred));
  double arg1 = ((*response) - (*worktransformlin[2])) / pow(sigma_2, 0.5);
  double arg2 = ((*response2p) - (*worktransformlin[3])) / pow((*worktransformlin[1]), 0.5);
  double prop1 = randnumbers::Phi2(arg1);
  double prop2 = randnumbers::Phi2(arg2);
  double hilfs = 1 - exp(-(*worktransformlin[0]));

  double l;

     l = -0.5*log(sigma_2)-pow((((*response))-(*worklin[2])),2)/(2*sigma_2)+log((*worktransformlin[0]) * hilfs
        * exp(- (*worktransformlin[0]) * (prop1 + prop2)) / pow(hilfs - (exp(-(*worktransformlin[0]) * prop1) - 1) * (exp(-(*worktransformlin[0]) * prop2) - 1), 2));

  modify_worklin();

  return l;

  }

void DISTR_frankcopula2_normal_sigma2::compute_iwls_wweightschange_weightsone(
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

    double sigma_2 = exp(*linpred);
    double arg1 = ((*response) - (*worktransformlin[2])) / pow(sigma_2, 0.5);
    double arg2 = ((*response2p) - (*worktransformlin[3])) / pow((*worktransformlin[1]), 0.5);
    double prop1 = randnumbers::Phi2(arg1);
    double prop2 = randnumbers::Phi2(arg2);
    double hilfs = 1 - exp(-(*worktransformlin[0]));
    double hilfs2 = exp(-(*worktransformlin[0]) * prop1) * (exp(-(*worktransformlin[0]) * prop2) - 1);
    double nenner = hilfs - (exp(-(*worktransformlin[0]) * prop1) - 1) * (exp(-(*worktransformlin[0]) * prop2) - 1);
    double rest = 1 + 2 * hilfs2 / nenner;

    double d1 = - 0.398942280401433 * exp(- 0.5 * pow(arg1, 2)) * arg1 * 0.5 ;
    double d2 = - 0.5 * d1 * (1 - pow(arg1, 2));


    double nu = -0.5 + (pow(((*response)-(*worklin[2])),2))/(2*sigma_2) - (*worktransformlin[0]) * d1 * rest;

    *workingweight = 0.5 + (*worktransformlin[0]) * d2 * (rest - 2 * (*worktransformlin[0]) * hilfs2 * (1 + hilfs2 / nenner) / nenner);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sigma_2)-pow((((*response))-(*worklin[2])),2)/(2*sigma_2)+log((*worktransformlin[0]) * hilfs
        * exp(- (*worktransformlin[0]) * (prop1 + prop2)) / pow(nenner, 2));

      }

  modify_worklin();

  }

void DISTR_frankcopula2_normal_sigma2::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_frankcopula2_normal_sigma2::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma2): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_frankcopula2_normal_sigma2::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS: DISTR_frankcopula2_normal_mu ----------------
//------------------------------------------------------------------------------


DISTR_frankcopula2_normal_mu::DISTR_frankcopula2_normal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  pos=p;
  family = "frank Copula - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_frankcopula2_normal_mu::DISTR_frankcopula2_normal_mu(const DISTR_frankcopula2_normal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
    pos = nd.pos;
    response2 = nd.response2;
    response2p = nd.response2p;
  }


const DISTR_frankcopula2_normal_mu & DISTR_frankcopula2_normal_mu::operator=(
                            const DISTR_frankcopula2_normal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }

void DISTR_frankcopula2_normal_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_frankcopula2_normal_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_frankcopula2_normal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_frankcopula2_normal_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[pos + 2]);
  }

 double DISTR_frankcopula2_normal_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_frankcopula2_normal_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    return 0;
    }

//double DISTR_frankcopula2_normal_mu::compute_quantile_residual_mult(vector<double *> response,
//                                             vector<double *> param,
//                                             vector<double *> weight,
//                                             vector<datamatrix *> aux)
//  {
//  double u_est;
//  u_est = cdf_mult(response,param,weight,aux);
//  double res_est = randnumbers::invPhi2(u_est);
//  return res_est;
//  }

double DISTR_frankcopula2_normal_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double l;
  double prop1 = randnumbers::Phi2(((*response) - mu) / pow((*worktransformlin[2]), 0.5));
  double prop2 = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[3]), 0.5));
  double hilfs = 1 - exp(-(*worktransformlin[0]));

  l = -pow((((*response))-mu),2)/(2*(*worktransformlin[2])) +log((*worktransformlin[0]) * hilfs
        * exp(- (*worktransformlin[0]) * (prop1 + prop2)) / pow(hilfs - (exp(-(*worktransformlin[0]) * prop1) - 1) * (exp(-(*worktransformlin[0]) * prop2) - 1), 2));

  modify_worklin();

  return l;

  }


void DISTR_frankcopula2_normal_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;


  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);
    double arg1 = ((*response) - mu) / pow((*worktransformlin[2]), 0.5);
    double arg2 = ((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[3]), 0.5);
    double prop1 = randnumbers::Phi2(arg1);
    double prop2 = randnumbers::Phi2(arg2);
    double hilfs = 1 - exp(-(*worktransformlin[0]));
    double hilfs2 = exp(-(*worktransformlin[0]) * prop1) * (exp(-(*worktransformlin[0]) * prop2) - 1);
    double nenner = hilfs - (exp(-(*worktransformlin[0]) * prop1) - 1) * (exp(-(*worktransformlin[0]) * prop2) - 1);
    double rest = 1 + 2 * hilfs2 / nenner;

    double d1 = - 0.398942280401433 * exp(- 0.5 * pow(arg1, 2)) / pow((*worktransformlin[2]), 0.5);
    double d2 = d1 * arg1 / pow((*worktransformlin[2]), 0.5);

    double nu = ((*response)-mu)/(*worktransformlin[2]) - (*worktransformlin[0]) * d1 * rest;

    *workingweight = 1/(*worktransformlin[2])
                     + (*worktransformlin[0]) * d2 * (rest - 2 * (*worktransformlin[0]) * hilfs2 * (1 + hilfs2 / nenner) / nenner);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((((*response))-mu),2)/(2*(*worktransformlin[2])) + log((*worktransformlin[0]) * hilfs
        * exp(- (*worktransformlin[0]) * (prop1 + prop2)) / pow(nenner, 2));

      }


  modify_worklin();

  }


void DISTR_frankcopula2_normal_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+2 * pos+1]));
  }


void DISTR_frankcopula2_normal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


/*void DISTR_frankcopula2_normal_mu::update(void)
  {

   unsigned i;

  double help;

  double * worktransformlinp;
  double * workweight;

  worktransformlinp = distrp[0]->helpmat1.getV();
  workweight = workingweight.getV();

  for (i=0;i<nrobs;i++,worktransformlinp++,workweight++)
    {
        *workweight = 1/(*worktransformlinp);
    }

  }*/


void DISTR_frankcopula2_normal_mu::update_end(void)
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


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_frankcopula2_rho --------------------
//------------------------------------------------------------------------------

void DISTR_frankcopula2_rho::check_errors(void)
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

DISTR_frankcopula2_rho::DISTR_frankcopula2_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Frank Copula - rho";

 datamatrix d(nrobs,1,0.0001);

  if (linpred_current==1)
    linearpred1.plus(d);
  else
    linearpred2.plus(d);
  check_errors();
  }


DISTR_frankcopula2_rho::DISTR_frankcopula2_rho(const DISTR_frankcopula2_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_frankcopula2_rho & DISTR_frankcopula2_rho::operator=(
                            const DISTR_frankcopula2_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_frankcopula2_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_frankcopula2_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[4]);
  }

void DISTR_frankcopula2_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_frankcopula2_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_frankcopula2_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {


   if (*weight[4] == 0)
     *deviance=0;
   else
     {
     double u = randnumbers::Phi2(((*response[3]) - (*linpred[3])) / pow(exp(*linpred[2]), 0.5));
     double v = randnumbers::Phi2(((*response[0]) - (*linpred[1])) / pow(exp(*linpred[0]), 0.5));

     double e1 = exp(-(*linpred[4]));
     double e1m1 = 1- e1;
     double e2 = exp(-(*linpred[4]) * u);
     double e3 = exp(-(*linpred[4]) * v);
     double e2m1 = e2 - 1;
     double e3m1 = e3 - 1;
     double l;

       l = log((*linpred[4]) * e1m1 * exp(- (*linpred[4]) * (u + v)) /  pow((e1m1 - e2m1 * e3m1), 2))
            +log(randnumbers::phi(((*response[3]) - (*linpred[3])) / pow(exp(*linpred[2]), 0.5))) +
                log(randnumbers::phi(((*response[0]) - (*linpred[1])) / pow(exp(*linpred[0]), 0.5)));



    *deviance = -2*l;
    }

  }

double DISTR_frankcopula2_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {



  if (counter==0)
    {
    set_worklin();
    }

    double u = randnumbers::Phi2(((*response) - (*worktransformlin[3])) / pow((*worktransformlin[2]), 0.5));
    double v = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[0]), 0.5));

    double e1m1 = 1 - exp(-(*linpred));
    double e2 = exp(-(*linpred) * u);
    double e3 = exp(-(*linpred) * v);
    double e2m1 = e2 - 1;
    double e3m1 = e3 - 1;
    double l;


    l =  log((*linpred) * e1m1 * exp(- (*linpred) * (u + v)) /  pow((e1m1 - e2m1 * e3m1), 2));

  modify_worklin();

  return l;

  }

void DISTR_frankcopula2_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }

    double u = randnumbers::Phi2(((*response) - (*worktransformlin[3])) / pow((*worktransformlin[2]), 0.5));
    double v = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[0]), 0.5));

    double e1 = exp(-(*linpred));
    double e1m1 = 1 - e1;
    double e2 = exp(-(*linpred) * u);
    double e3 = exp(-(*linpred) * v);
    double e2m1 = e2 - 1;
    double e3m1 = e3 - 1;

    double nu =  1 / (*linpred) + e1 / e1m1 - (u + v)
                - 2 * (e1 + (u + v) * e2 * e3 - u * e2 - v * e3) / ((e1m1 - e2m1 * e3m1));


    *workingweight =  1 / pow((*linpred), 2) + e1 / pow(e1m1, 2)
                    - 2 * (e1 + pow((u + v), 2) * e2 * e3 - pow(u, 2) * e2 - pow(v, 2) * e3) / ((e1m1 - e2m1 * e3m1))
                    - 2 * pow(((e1 + (u + v) * e2 * e3 - u * e2 - v * e3) / ((e1m1 - e2m1 * e3m1))), 2);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

          like += log((*linpred) * e1m1 * exp(- (*linpred) * (u + v)) /  pow((e1m1 - e2m1 * e3m1), 2));
      }


  modify_worklin();

  }

void DISTR_frankcopula2_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_frankcopula2_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_frankcopula2_rho::update_end(void)
  {

  // helpmat1 stores rho2

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
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_frankcopula2_exp_rho --------------------
//------------------------------------------------------------------------------

void DISTR_frankcopula2_exp_rho::check_errors(void)
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

DISTR_frankcopula2_exp_rho::DISTR_frankcopula2_exp_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Frank Copula - rho";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "true";
   linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_frankcopula2_exp_rho::DISTR_frankcopula2_exp_rho(const DISTR_frankcopula2_exp_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_frankcopula2_exp_rho & DISTR_frankcopula2_exp_rho::operator=(
                            const DISTR_frankcopula2_exp_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_frankcopula2_exp_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_frankcopula2_exp_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[4]);
  }

void DISTR_frankcopula2_exp_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_frankcopula2_exp_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_frankcopula2_exp_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {


   if (*weight[4] == 0)
     *deviance=0;
   else
     {
     double u = randnumbers::Phi2(((*response[3]) - (*linpred[3])) / pow(exp(*linpred[2]), 0.5));
     double v = randnumbers::Phi2(((*response[0]) - (*linpred[1])) / pow(exp(*linpred[0]), 0.5));

     double e1 = exp(-exp(*linpred[4]));
     double e1m1 = 1- e1;
     double e2 = exp(-exp(*linpred[4]) * u);
     double e3 = exp(-exp(*linpred[4]) * v);
     double e2m1 = e2 - 1;
     double e3m1 = e3 - 1;
     double l;

       l = log(exp(*linpred[4]) * e1m1 * exp(- exp(*linpred[4]) * (u + v)) /  pow((e1m1 - e2m1 * e3m1), 2))
            +log(randnumbers::phi(((*response[3]) - (*linpred[3])) / pow(exp(*linpred[2]), 0.5))) +
                log(randnumbers::phi(((*response[0]) - (*linpred[1])) / pow(exp(*linpred[0]), 0.5)));



    *deviance = -2*l;
    }

  }

double DISTR_frankcopula2_exp_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {



  if (counter==0)
    {
    set_worklin();
    }

    double u = randnumbers::Phi2(((*response) - (*worktransformlin[3])) / pow((*worktransformlin[2]), 0.5));
    double v = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[0]), 0.5));
    double thet = exp((*linpred));
    double e1m1 = 1 - exp(-thet);
    double e2 = exp(-thet * u);
    double e3 = exp(-thet * v);
    double e2m1 = e2 - 1;
    double e3m1 = e3 - 1;
    double l;


    l =  log(thet * e1m1 * exp(- thet * (u + v)) /  pow((e1m1 - e2m1 * e3m1), 2));

  modify_worklin();

  return l;

  }

void DISTR_frankcopula2_exp_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }

    double u = randnumbers::Phi2(((*response) - (*worktransformlin[3])) / pow((*worktransformlin[2]), 0.5));
    double v = randnumbers::Phi2(((*response2p) - (*worktransformlin[1])) / pow((*worktransformlin[0]), 0.5));
    double thet = exp((*linpred));
    double e1 = exp(-thet);
    double e1m1 = 1 - e1;
    double e2 = exp(-thet * u);
    double e3 = exp(-thet * v);
    double e2m1 = e2 - 1;
    double e3m1 = e3 - 1;


    double zaehler = thet * (e1 + (u + v) * e2 * e3 - u * e2 - v * e3) ;
    double stars =  (e1 + pow((u + v), 2) * e2 * e3 - pow(u, 2) * e2 - pow(v, 2) * e3);
    double nenner = e1m1 - e2m1 * e3m1;

    double nu =  1 + thet * e1 / e1m1 - thet * (u + v) - 2 * zaehler / nenner;

    *workingweight =  -nu + 1 + thet * e1 / e1m1 * (thet + e1 / e1m1) - 2 * thet * thet * stars / nenner - 2 * pow(zaehler / nenner, 2);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

          like += log(thet * e1m1 * exp(- thet * (u + v)) /  pow((nenner), 2));
      }


  modify_worklin();

  }

void DISTR_frankcopula2_exp_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_frankcopula2_exp_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_frankcopula2_exp_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_frankcopula_rho -----------------------
//------------------------------------------------------------------------------
void DISTR_frankcopula_rho::check_errors(void)
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

            if((*workresp) > 1 )
            {
                errors=true;
                errormessages.push_back("ERROR: cdfs of marginals take values inbetween zero and one!\n");
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


DISTR_frankcopula_rho::DISTR_frankcopula_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Frank Copula - rho";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
 datamatrix d(nrobs,1,0.0001);

  if (linpred_current==1)
    linearpred1.plus(d);
  else
    linearpred2.plus(d);

  }


DISTR_frankcopula_rho::DISTR_frankcopula_rho(const DISTR_frankcopula_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_frankcopula_rho & DISTR_frankcopula_rho::operator=(
                            const DISTR_frankcopula_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_frankcopula_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_frankcopula_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2]);
  }

void DISTR_frankcopula_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_frankcopula_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_frankcopula_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {



   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double e1 = exp(-(*linpred[2]));
     double e1m1 = 1- e1;
     double e2 = exp(-(*linpred[2]) * (*response[1]));
     double e3 = exp(-(*linpred[2]) * (*response[0]));
     double e2m1 = e2 - 1;
     double e3m1 = e3 - 1;
     double l;

      l = log((*linpred[2]) * e1m1 * exp(- (*linpred[2]) * ((*response[1]) + (*response[0]))) /  pow((e1m1 - e2m1 * e3m1), 2));


    *deviance = -2*l;
    }

  }

double DISTR_frankcopula_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {



  if (counter==0)
    {
    set_worklin();
    }
    double e1m1 = 1 - exp(-(*linpred));
    double e2 = exp(-(*linpred) * (*response));
    double e3 = exp(-(*linpred) * (*response2p));
    double e2m1 = e2 - 1;
    double e3m1 = e3 - 1;
    double l;


    l =  log((*linpred) * e1m1 * exp(- (*linpred) * ((*response) + (*response2p))) /  pow((e1m1 - e2m1 * e3m1), 2));


  modify_worklin();

  return l;

  }

void DISTR_frankcopula_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }
    double e1 = exp(-(*linpred));
    double e1m1 = 1 - e1;
    double e2 = exp(-(*linpred) * (*response));
    double e3 = exp(-(*linpred) * (*response2p));
    double e2m1 = e2 - 1;
    double e3m1 = e3 - 1;

    double nu =  1 / (*linpred) + e1 / e1m1 - ((*response) + (*response2p))
                - 2 * (e1 + ((*response) + (*response2p)) * e2 * e3 - (*response) * e2 - (*response2p) * e3) / ((e1m1 - e2m1 * e3m1));


    *workingweight =  1 / pow((*linpred), 2) + e1 / pow(e1m1, 2)
                    - 2 * (e1 + pow(((*response) + (*response2p)), 2) * e2 * e3 - pow((*response), 2) * e2 - pow((*response2p), 2) * e3) / ((e1m1 - e2m1 * e3m1))
                    - 2 * pow(((e1 + ((*response) + (*response2p)) * e2 * e3 - (*response) * e2 - (*response2p) * e3) / ((e1m1 - e2m1 * e3m1))), 2);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

          like +=  log((*linpred) * e1m1 * exp(- (*linpred) * ((*response) + (*response2p))) /  pow((e1m1 - e2m1 * e3m1), 2));
      }

  modify_worklin();

  }

void DISTR_frankcopula_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_frankcopula_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_frankcopula_rho::update_end(void)
  {

  // helpmat1 stores rho

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
    }

  }

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_frankcopula_exp_rho -----------------------
//------------------------------------------------------------------------------
void DISTR_frankcopula_exp_rho::check_errors(void)
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

            if((*workresp) > 1 )
            {
                errors=true;
                errormessages.push_back("ERROR: cdfs of marginals take values inbetween zero and one!\n");
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


DISTR_frankcopula_exp_rho::DISTR_frankcopula_exp_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Frank Copula - rho";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "rho";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_frankcopula_exp_rho::DISTR_frankcopula_exp_rho(const DISTR_frankcopula_exp_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_frankcopula_exp_rho & DISTR_frankcopula_exp_rho::operator=(
                            const DISTR_frankcopula_exp_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_frankcopula_exp_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_frankcopula_exp_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[2]);
  }

void DISTR_frankcopula_exp_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_frankcopula_exp_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

void DISTR_frankcopula_exp_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {



   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double thet = exp((*linpred[2]));
     double e1 = exp(-thet);
     double e1m1 = 1- e1;
     double e2 = exp(-thet * (*response[1]));
     double e3 = exp(-thet * (*response[0]));
     double e2m1 = e2 - 1;
     double e3m1 = e3 - 1;
     double l;

      l = log(thet * e1m1 * exp(- thet * ((*response[1]) + (*response[0]))) /  pow((e1m1 - e2m1 * e3m1), 2));


    *deviance = -2*l;
    }

  }

double DISTR_frankcopula_exp_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {



  if (counter==0)
    {
    set_worklin();
    }
    double thet = exp((*linpred));
    double e1m1 = 1 - exp(-thet);
    double e2 = exp(-thet * (*response));
    double e3 = exp(-thet * (*response2p));
    double e2m1 = e2 - 1;
    double e3m1 = e3 - 1;
    double l;


    l =  log(thet * e1m1 * exp(- thet * ((*response) + (*response2p))) /  pow((e1m1 - e2m1 * e3m1), 2));


  modify_worklin();

  return l;

  }

void DISTR_frankcopula_exp_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  if (counter==0)
    {
    set_worklin();
    }
    double thet = exp((*linpred));
    double e1 = exp(-thet);
    double e1m1 = 1 - e1;
    double e2 = exp(-thet * (*response));
    double e3 = exp(-thet * (*response2p));
    double e2m1 = e2 - 1;
    double e3m1 = e3 - 1;


    double zaehler = thet * (e1 + ((*response) + (*response2p)) * e2 * e3 - (*response) * e2 - (*response2p) * e3) ;
    double stars =  (e1 + pow(((*response) + (*response2p)), 2) * e2 * e3 - pow((*response), 2) * e2 - pow((*response2p), 2) * e3);
    double nenner = e1m1 - e2m1 * e3m1;

    double nu =  1 + thet * e1 / e1m1 - thet * ((*response) + (*response2p)) - 2 * zaehler / nenner;

    *workingweight =  -nu + 1 + thet * thet * e1 / e1m1 * (1 + e1 / e1m1) - 2 * thet * thet * stars / nenner - 2 * pow(zaehler / nenner, 2);

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

          like += log(thet * e1m1 * exp(- thet * ((*response) + (*response2p))) /  pow((nenner), 2));
      }

  modify_worklin();

  }

void DISTR_frankcopula_exp_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_frankcopula_exp_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (rho): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_frankcopula_exp_rho::update_end(void)
  {

  // helpmat1 stores rho

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_copula ------------------------------
//------------------------------------------------------------------------------


DISTR_copula::DISTR_copula(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  pos =p;
  family = "Copula";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "u";
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  }


DISTR_copula::DISTR_copula(const DISTR_copula & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_copula & DISTR_copula::operator=(
                            const DISTR_copula & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }





double DISTR_copula::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_copula::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = 0;
  }

void DISTR_copula::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_copula::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_copula::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {


  if (counter==0)
    {
    set_worklin();
    }

  double l;


     l = 0;

  modify_worklin();

  return l;

  }


void DISTR_copula::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {



  if (counter==0)
    {
    set_worklin();
    }


    double nu = 0;

    *workingweight = 1;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += 0;

      }


  modify_worklin();

  }




void DISTR_copula::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_copula::update_end(void)
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


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_dirichlet ---------------------------
//------------------------------------------------------------------------------

void DISTR_dirichlet::check_errors(void)
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

        if ((*workresp>1) | (*workresp<0) )
          {
          errors=true;
          errormessages.push_back("ERROR: response has to be between zero and one\n");
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


DISTR_dirichlet::DISTR_dirichlet(GENERAL_OPTIONS * o,
                                           const datamatrix & r, int & nrc, unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,nrc-1,w)
  {
  family = "Dirichlet Distribution";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "alpha";
  linpredminlimit=-10;
  linpredmaxlimit=15;
  nrcat = nrc;
  pos = p;
  check_errors();
  }


DISTR_dirichlet::DISTR_dirichlet(const DISTR_dirichlet & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
    nrcat = nd.nrcat;
    pos = nd.pos;
  }


const DISTR_dirichlet & DISTR_dirichlet::operator=(
                            const DISTR_dirichlet & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  nrcat = nd.nrcat;
  pos = nd.pos;
  return *this;
  }


void DISTR_dirichlet::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response


   if (*weight[0] == 0)
     *deviance=0;
   else
     {
     double sum_alpha = 0;
     double sum_log_gamma = 0;
     double sum_rest = 0;

     unsigned i;
     for(i=0;i<nrcat;i++) {
        double hilfs = exp(*linpred[i]);
        sum_log_gamma += randnumbers::lngamma_exact(hilfs);
        sum_alpha += hilfs;
        sum_rest += (hilfs-1)*log(*response[i]);
     }

     double l = -sum_log_gamma + randnumbers::lngamma_exact(sum_alpha) + sum_rest;


    *deviance = -2*l;
    }

  }


 double DISTR_dirichlet::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

  void DISTR_dirichlet::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[pos]);
  }

 double DISTR_dirichlet::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

 double DISTR_dirichlet::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_dirichlet::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha_i equation

  if (counter==0)
    {
    set_worklin();
    }
     double alpha_current = exp(*linpred);
     double sum_alpha = alpha_current;

     unsigned i;
     for(i=0;i<(nrcat-1);i++) {
        sum_alpha += (*worktransformlin[i]);
     }

     double l = -randnumbers::lngamma_exact(alpha_current) + randnumbers::lngamma_exact(sum_alpha) + alpha_current*log(*response);

  modify_worklin();

  return l;

  }


void DISTR_dirichlet::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

     double alpha_current = exp(*linpred);
     double sum_alpha = alpha_current;

     unsigned i;
     for(i=0;i<(nrcat-1);i++) {
        sum_alpha += (*worktransformlin[i]);
     }

    double nu = alpha_current*( -randnumbers::digamma_exact(alpha_current) +randnumbers::digamma_exact(sum_alpha) + log(*response) );

    *workingweight = pow(alpha_current,2)*( randnumbers::trigamma_exact(alpha_current) - randnumbers::trigamma_exact(sum_alpha) );

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

      like += -randnumbers::lngamma_exact(alpha_current) + randnumbers::lngamma_exact(sum_alpha) + alpha_current*log(*response);

      }


  modify_worklin();

  }


void DISTR_dirichlet::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
     double alpha_current = exp(*linpred[predstart_mumult+pos]);
     double sum_alpha = 0;

     unsigned i;
     for(i=0;i<nrcat;i++) {
        double hilfs = exp(*linpred[predstart_mumult+i]);
        sum_alpha += hilfs;
        }

     *mu = alpha_current/sum_alpha;
  }


void DISTR_dirichlet::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (alpha): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_dirichlet::update_end(void)
  {


  // helpmat1 stores 1-exp(-exp_lin)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  // double exp_lin;
  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_bivt_df -------------------------------
//------------------------------------------------------------------------------


DISTR_bivt_df::DISTR_bivt_df(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,5,w)
  {
  family = "Bivariate t-Distribution - Degrees of Freedom";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "df";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_bivt_df::DISTR_bivt_df(const DISTR_bivt_df & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivt_df & DISTR_bivt_df::operator=(
                            const DISTR_bivt_df & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivt_df::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivt_df::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[0]));
  }

void DISTR_bivt_df::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivt_df::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }


double DISTR_bivt_df::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }

  double degf = exp((*linpred));

  double oneminusrho2 = 1- pow((*worktransformlin[0]),2);
  double np2d2 = (degf+2)/2;
  double nd2 = degf/2;
  double X_1 = ((*response)-(*worklin[4]))/(*worktransformlin[2]);
  double X_2 = ((*response2p)-(*worklin[3]))/(*worktransformlin[1]);
  double l;


     l = randnumbers::lngamma_exact(np2d2)-randnumbers::lngamma_exact(nd2)-log(degf)
     -np2d2*log( 1 + (1/(degf*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[0])*X_1*X_2+pow(X_2,2) ) );


  modify_worklin();

  return l;

  }

void DISTR_bivt_df::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double degf = exp((*linpred));

  double oneminusrho2 = 1- pow((*worktransformlin[0]),2);
  double np2d2 = (degf+2)/2;
  double nd2 = degf/2;
  double X_1 = ((*response)-(*worklin[4]))/(*worktransformlin[2]);
  double X_2 = ((*response2p)-(*worklin[3]))/(*worktransformlin[1]);
  double nenner_C = 1 + (1/(degf*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[0])*X_1*X_2+pow(X_2,2) );

    double nu = nd2*( randnumbers::digamma_exact(np2d2)-randnumbers::digamma_exact(nd2)-log(nenner_C) ) - 1 + np2d2*(nenner_C-1)/nenner_C;

    *workingweight =  - pow(nd2,2)*( randnumbers::trigamma_exact(np2d2) - randnumbers::trigamma_exact(nd2) ) -2*degf/(degf+2)-degf*(degf+2)/(2*(degf+4))+degf/2; //- degf/(degf+2);

    *workingresponse = *linpred + nu/(*workingweight);

  //  cout << "df equation y1: " << *response << endl;
  //  cout << "df equation y2: " << *response2p << endl;


    if (compute_like)
      {

        like +=  randnumbers::lngamma_exact(np2d2)-randnumbers::lngamma_exact(nd2)-log(degf)
     -np2d2*log( nenner_C );

      }

  modify_worklin();

  }


void DISTR_bivt_df::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (df): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivt_df::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_bivt_rho ------------------------------
//------------------------------------------------------------------------------


DISTR_bivt_rho::DISTR_bivt_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,5,w)
  {
  family = "Bivariate t-Distribution - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
    linpredminlimit=-100;
  linpredmaxlimit=100;

  }


DISTR_bivt_rho::DISTR_bivt_rho(const DISTR_bivt_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivt_rho & DISTR_bivt_rho::operator=(
                            const DISTR_bivt_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivt_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivt_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = (*linpred[1]);
  *param = arg/pow(1+pow(arg,2),0.5);
  }

void DISTR_bivt_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivt_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_bivt_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

 // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }


  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double X_1 = ((*response)-(*worklin[4]))/(*worktransformlin[2]);
  double X_2 = ((*response2p)-(*worklin[3]))/(*worktransformlin[1]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*rho*X_1*X_2+pow(X_2,2) );
  double l;


     l = -0.5*log(oneminusrho2)-(((*worktransformlin[0])+2)/2)*log( nenner_C );



  modify_worklin();

  return l;

  }

void DISTR_bivt_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

 // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of sigma_2 equation
  // *worktransformlin[1] = sigma_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double X_1 = ((*response)-(*worklin[4]))/(*worktransformlin[2]);
  double X_2 = ((*response2p)-(*worklin[3]))/(*worktransformlin[1]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*rho*X_1*X_2+pow(X_2,2) );

  double nu = oneminusrho2*(*linpred) - (((*worktransformlin[0])+2)/((*worktransformlin[0])*nenner_C))*( (*linpred)*(pow(X_1,2)+pow(X_2,2)) -
                                                                                                      (pow(1/oneminusrho2,0.5)+rho*(*linpred))*X_1*X_2 );



    *workingweight = 1-pow(rho2,2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(oneminusrho2)-(((*worktransformlin[0])+2)/2)*log( nenner_C );

      }

  modify_worklin();

  }


void DISTR_bivt_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho): fisher z-transformation\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivt_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin)/pow((1+pow((*worklin),2)),0.5);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_bivt_sigma ----------------------------
//------------------------------------------------------------------------------


DISTR_bivt_sigma::DISTR_bivt_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,5,w)
  {
  family = "Bivariate t-Distribution - sigma";

  pos = p;
  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_bivt_sigma::DISTR_bivt_sigma(const DISTR_bivt_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivt_sigma & DISTR_bivt_sigma::operator=(
                            const DISTR_bivt_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivt_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivt_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[pos + 2]));
  }

void DISTR_bivt_sigma::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivt_sigma::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }


double DISTR_bivt_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of rho equation
  // *worktransformlin[1] = rho;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }

  double sig_current = exp((*linpred));

  double oneminusrho2 = 1- pow((*worktransformlin[1]),2);
  double X_1 = ((*response)-(*worklin[3]))/(sig_current);
  double X_2 = ((*response2p)-(*worklin[4]))/(*worktransformlin[2]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[1])*X_1*X_2+pow(X_2,2) );
  double l;


     l = -log(sig_current)
     -(((*worktransformlin[0])+2)/2)*log( nenner_C );


  modify_worklin();

  return l;

  }

void DISTR_bivt_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of rho equation
  // *worktransformlin[1] = rho;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of mu_2 equation
  // *worktransformlin[3] = mu_2;
  // *worklin[4] = linear predictor of mu_1 equation
  // *worktransformlin[4] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double sig_current = exp((*linpred));

  double oneminusrho2 = 1- pow((*worktransformlin[1]),2);
  double X_1 = ((*response)-(*worklin[3]))/(sig_current);
  double X_2 = ((*response2p)-(*worklin[4]))/(*worktransformlin[2]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[1])*X_1*X_2+pow(X_2,2) );

    double nu = -1 - (((*worktransformlin[0])+2)/(nenner_C*(*worktransformlin[0])*oneminusrho2))*((*worktransformlin[1])*X_1*X_2-pow(X_1,2));



    *workingweight = 1+1/oneminusrho2 ;
//- 2/(*worktransformlin[0]) - 2*pow((*worktransformlin[1]),2)/((*worktransformlin[0])*oneminusrho2)
    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -log(sig_current)-(((*worktransformlin[0])+2)/2)*log( nenner_C );

      }

  modify_worklin();

  }


void DISTR_bivt_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivt_sigma::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_bivt_mu -----------------------------
//------------------------------------------------------------------------------
void DISTR_bivt_mu::check_errors(void)
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


DISTR_bivt_mu::DISTR_bivt_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,5,w)
  {
  pos =p;
  family = "Bivariate t-Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  check_errors();
  }


DISTR_bivt_mu::DISTR_bivt_mu(const DISTR_bivt_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivt_mu & DISTR_bivt_mu::operator=(
                            const DISTR_bivt_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


void DISTR_bivt_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[2] = *response[4] = first component of two dimensional reponse
   // *linpred[0] = eta_df
   // *linpred[1] = eta_rho
   // *linpred[2] = eta_sigma_2
   // *linpred[3] = eta_mu_2
   // *linpred[4] = eta_sigma_1
   // *linpred[5] = eta_mu_1

   if (*weight[5] == 0)
     *deviance=0;
   else
     {
     double rho = (*linpred[1])/pow(1+pow((*linpred[1]),2),0.5);
     double degf = exp(*linpred[0]);
     double np2d2 = (degf+2)/2;
     double nd2 = degf/2;
     double sigma_2 = exp(*linpred[2]);
     double mu_2 = (*linpred[4]);
     double sigma_1 = exp(*linpred[3]);
     double mu_1 = (*linpred[5]);
     double hilfs1 = 1-pow(rho,2);
     double X_1 = ((*response[5])-mu_1)/sigma_1;
     double X_2 = ((*response[4])-mu_2)/sigma_2;
     double nenner_C = 1+(1/(degf*hilfs1))*(pow(X_1,2)-2*rho*X_1*X_2+pow(X_2,2));
     double l;

       l = randnumbers::lngamma_exact(np2d2)-randnumbers::lngamma_exact(nd2)-log(degf)-log(PI)-log(sigma_1)-log(sigma_2)-0.5*log(hilfs1)-
           np2d2*log(nenner_C);


    *deviance = -2*l;
    }

  }


double DISTR_bivt_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivt_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[pos + 4]);
  }

 double DISTR_bivt_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_bivt_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
  //  double arg = ((*response[5])-(*param[5]))/pow((*param[3]),0.5) ;
  //  double u = gsl_cdf_tdist_P(arg, (*param[0]));
   // return (u);
   return 0;
    }

void DISTR_bivt_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivt_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_bivt_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of rho equation
  // *worktransformlin[1] = rho;
  // *worklin[2] = linear predictor of mu_1 equation
  // *worktransformlin[2] = mu_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;
  // *worklin[4] = linear predictor of sigma_1 equation
  // *worktransformlin[4] = sigma_1;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double oneminusrho2 = 1- pow((*worktransformlin[1]),2);
  double X_1 = ((*response)-mu)/((*worktransformlin[3]));
  double X_2 = ((*response2p)-(*worklin[2]))/(*worktransformlin[4]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[1])*X_1*X_2+pow(X_2,2) );
  double l;


     l = -(((*worktransformlin[0])+2)/2)*log( nenner_C );

  modify_worklin();

  return l;

  }


void DISTR_bivt_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of df equation
  // *worktransformlin[0] = df;
  // *worklin[1] = linear predictor of rho equation
  // *worktransformlin[1] = rho;
  // *worklin[2] = linear predictor of mu_1 equation
  // *worktransformlin[2] = mu_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;
  // *worklin[4] = linear predictor of sigma_1 equation
  // *worktransformlin[4] = sigma_1;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);

  double oneminusrho2 = 1- pow((*worktransformlin[1]),2);
  double X_1 = ((*response)-mu)/((*worktransformlin[3]));
  double X_2 = ((*response2p)-(*worklin[2]))/(*worktransformlin[4]);
  double nenner_C = 1 + (1/((*worktransformlin[0])*oneminusrho2))*( pow(X_1,2)-2*(*worktransformlin[1])*X_1*X_2+pow(X_2,2) );

    double nu =(((*worktransformlin[0])+2)/((*worktransformlin[0])*oneminusrho2*(*worktransformlin[3])*nenner_C))*(X_1-(*worktransformlin[1])*X_2);

    *workingweight = (1)/(oneminusrho2*pow((*worktransformlin[3]),2)) ;
    //- ((*worktransformlin[0])+2)/(pow((*worktransformlin[0]),2)*oneminusrho2*pow((*worktransformlin[3]),2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -(((*worktransformlin[0])+2)/2)*log( nenner_C );

      }


  modify_worklin();

  }


void DISTR_bivt_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
    if(exp((*linpred[predstart_mumult]))>2)
    {
        *mu = ((*linpred[predstart_mumult+4+pos]));
    } else
    {
        *mu = 0;
    }

  }


void DISTR_bivt_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }



void DISTR_bivt_mu::update_end(void)
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

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_bivnormal_rhofz ----------------------
//------------------------------------------------------------------------------


DISTR_bivnormal_rhofz::DISTR_bivnormal_rhofz(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Bivariate Normal Distribution - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
    linpredminlimit=-20.6;
  linpredmaxlimit=20.6;

  }


DISTR_bivnormal_rhofz::DISTR_bivnormal_rhofz(const DISTR_bivnormal_rhofz & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivnormal_rhofz & DISTR_bivnormal_rhofz::operator=(
                            const DISTR_bivnormal_rhofz & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivnormal_rhofz::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivnormal_rhofz::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = exp(2 * (*linpred[0]));
  *param = (arg - 1) / (arg + 1);
  }

void DISTR_bivnormal_rhofz::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivnormal_rhofz::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_bivnormal_rhofz::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma_2 equation
  // *worktransformlin[0] = sigma_2;
  // *worklin[1] = linear predictor of sigma_1 equation
  // *worktransformlin[1] = sigma_1;
  // *worklin[2] = linear predictor of mu_2 equation
  // *worktransformlin[2] = mu_2;
  // *worklin[3] = linear predictor of mu_1 equation
  // *worktransformlin[3] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }
  double rho;

 /* if (*linpred <= -20.6)
    rho  = -0.99995;
  else if (*linpred >= 20.6)
   rho  = 0.99995;
  else*/
    rho = (exp(2 * (*linpred)) - 1) / (exp(2 * (*linpred)) + 1);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response))-(*worklin[3])),2)/pow((*worktransformlin[1]),2) -
                                 2*rho*(((*response)-(*worklin[3]))/((*worktransformlin[1])))*(((*response2p)-(*worklin[2]))/((*worktransformlin[0])))
                                +  pow((((*response2p))-(*worklin[2])),2)/pow((*worktransformlin[0]),2) );


  modify_worklin();

  return l;

  }

void DISTR_bivnormal_rhofz::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma_2 equation
  // *worktransformlin[0] = sigma_2;
  // *worklin[1] = linear predictor of sigma_1 equation
  // *worktransformlin[1] = sigma_1;
  // *worklin[2] = linear predictor of mu_2 equation
  // *worktransformlin[2] = mu_2;
  // *worklin[3] = linear predictor of mu_1 equation
  // *worktransformlin[3] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double rho;
   /* if (*linpred <= -20.6)
    rho  = -0.99995;
  else if (*linpred >= 20.6)
   rho  = 0.99995;
  else*/
    rho = (exp(2 * (*linpred)) - 1) / (exp(2 * (*linpred)) + 1);


  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;


  double nu = rho - (rho / oneminusrho2) *( pow((((*response))-(*worktransformlin[3])),2)/pow((*worktransformlin[1]),2)
                                                      +  pow((((*response2p))-(*worktransformlin[2])),2)/pow((*worktransformlin[0]),2) )
                + 2 * (1 / oneminusrho2 - 0.5)*( (((*response)-(*worktransformlin[3]))/((*worktransformlin[1])))*(((*response2p)-(*worktransformlin[2]))/((*worktransformlin[0]))) );



    *workingweight =  (rho2) + 1;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response))-(*worklin[3])),2)/pow((*worktransformlin[1]),2) -
                                 2*rho*(((*response)-(*worklin[3]))/((*worktransformlin[1])))*(((*response2p)-(*worklin[2]))/((*worktransformlin[0])))
                                +  pow((((*response2p))-(*worklin[2])),2)/pow((*worktransformlin[0]),2) );

      }

  modify_worklin();

  }


void DISTR_bivnormal_rhofz::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho): fishers z-transformation\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivnormal_rhofz::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
 /*  if (*worklin <= -10.6)
     *pmu  = -0.99995;
   else if (*worklin >= 10.6)
     *pmu  = 0.99995;
   else*/
    *pmu = (exp(2 * (*worklin)) - 1) / (exp(2 *(*worklin)) + 1);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_bivnormal_mufz ------------------------
//------------------------------------------------------------------------------
void DISTR_bivnormal_mufz::check_errors(void)
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


DISTR_bivnormal_mufz::DISTR_bivnormal_mufz(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  pos =p;
  family = "Bivariate Normal Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  if(o->forceIWLS)
    {
    updateIWLS = true;
    }
  else
    {
    updateIWLS = false;
    }
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  check_errors();
  }


DISTR_bivnormal_mufz::DISTR_bivnormal_mufz(const DISTR_bivnormal_mufz & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivnormal_mufz & DISTR_bivnormal_mufz::operator=(
                            const DISTR_bivnormal_mufz & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


void DISTR_bivnormal_mufz::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[2] = *response[4] = first component of two dimensional reponse
   // *linpred[0] = eta_rho
   // *linpred[1] = eta_sigma_2
   // *linpred[2] = eta_mu_2
   // *linpred[3] = eta_sigma_1
   // *linpred[4] = eta_mu_1

   if (*weight[4] == 0)
     *deviance=0;
   else
     {
     double rho = (exp(2 * (*linpred[0])) - 1)/(exp(2 * (*linpred[0])) + 1);
     double sigma_2 = exp(*linpred[1]);
     double mu_2 = (*linpred[3]);
     double sigma_1 = exp(*linpred[2]);
     double mu_1 = (*linpred[4]);
     double hilfs1 = 1-pow(rho,2);
     double l;

       l = -log(2*PI)-log(sigma_1)-log(sigma_2)-0.5*log(hilfs1)-(1/(2*hilfs1))*( pow((((*response[4]))-mu_1),2)/pow(sigma_1,2) -
                                                                                2*rho*(((*response[4])-mu_1)/(sigma_1))*(((*response[3])-mu_2)/(sigma_2))
                                                                                + pow((((*response[3]))-mu_2),2)/pow(sigma_2,2) );


    *deviance = -2*l;
    }

  }


double DISTR_bivnormal_mufz::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivnormal_mufz::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[3 + pos]);
  }

void DISTR_bivnormal_mufz::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();
  weightp = weight.getV();

  }



void DISTR_bivnormal_mufz::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    weightp++;
    }

  }



double DISTR_bivnormal_mufz::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double rho2 = pow((*worktransformlin[0]),2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -(1/(2*oneminusrho2))*( pow((((*response))-mu),2)/pow((*worktransformlin[2]),2) -
                                 2*(*worktransformlin[0])*(((*response)-mu)/((*worktransformlin[2])))*(((*response2p)-(*worktransformlin[1]))/((*worktransformlin[3]))) );

  modify_worklin();

  return l;

  }


void DISTR_bivnormal_mufz::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);
    double rho2 = pow((*worktransformlin[0]),2);
    double oneminusrho2 = 1- rho2;


    double nu = (1/(oneminusrho2))*( (((*response))-mu)/pow((*worktransformlin[2]),2) -
                                 ((*worktransformlin[0])/(*worktransformlin[2]))*(((*response2p)-(*worktransformlin[1]))/((*worktransformlin[3]))) );

//    if(updateIWLS)
      *workingweight = 1/(oneminusrho2*pow((*worktransformlin[2]),2));
//    else
//      *workingweight = (*weightp)/(oneminusrho2*pow((*worktransformlin[2]),2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {
      like += -(1/(2*oneminusrho2))*( pow((((*response))-mu),2)/pow((*worktransformlin[2]),2) -
                                 2*(*worktransformlin[0])*(((*response)-mu)/((*worktransformlin[2])))*(((*response2p)-(*worktransformlin[1]))/((*worktransformlin[3]))) );
      }


  modify_worklin();

  }


void DISTR_bivnormal_mufz::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+3+pos]));
  }


void DISTR_bivnormal_mufz::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }

/*void DISTR_bivnormal_mufz::update(void)
  {

   unsigned i;

  double help;

  double * worktransformlinr;
  double * worktransformlins;
  double * workweight;

  worktransformlinr = distrp[0]->helpmat1.getV();
  worktransformlins = distrp[2]->helpmat1.getV();
  workweight = workingweight.getV();

  if(updateIWLS)
    {
    for (i=0;i<nrobs;i++,worktransformlinr++,worktransformlins++,workweight++)
      {
      *workweight = 1/((1 - pow((*worktransformlinr),2))*pow((*worktransformlins),2));
      }
    }
  else
    {
    weightp = weight.getV();
    for (i=0;i<nrobs;i++,worktransformlinr++,worktransformlins++,workweight++,weightp++)
      {
      *workweight = (*weightp)/((1 - pow((*worktransformlinr),2))*pow((*worktransformlins),2));
      }
    }
  }*/

void DISTR_bivnormal_mufz::update_end(void)
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



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_bivnormal_rho -------------------------
//------------------------------------------------------------------------------


DISTR_bivnormal_rho::DISTR_bivnormal_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Bivariate Normal Distribution - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
    linpredminlimit=-100;
  linpredmaxlimit=100;

  }


DISTR_bivnormal_rho::DISTR_bivnormal_rho(const DISTR_bivnormal_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivnormal_rho & DISTR_bivnormal_rho::operator=(
                            const DISTR_bivnormal_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivnormal_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivnormal_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = (*linpred[0]);
  *param = arg/pow(1+pow(arg,2),0.5);
  }

void DISTR_bivnormal_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivnormal_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_bivnormal_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma_2 equation
  // *worktransformlin[0] = sigma_2;
  // *worklin[1] = linear predictor of sigma_1 equation
  // *worktransformlin[1] = sigma_1;
  // *worklin[2] = linear predictor of mu_2 equation
  // *worktransformlin[2] = mu_2;
  // *worklin[3] = linear predictor of mu_1 equation
  // *worktransformlin[3] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }
  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response))-(*worklin[3])),2)/pow((*worktransformlin[1]),2) -
                                 2*rho*(((*response)-(*worklin[3]))/((*worktransformlin[1])))*(((*response2p)-(*worklin[2]))/((*worktransformlin[0])))
                                +  pow((((*response2p))-(*worklin[2])),2)/pow((*worktransformlin[0]),2) );


  modify_worklin();

  return l;

  }

void DISTR_bivnormal_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma_2 equation
  // *worktransformlin[0] = sigma_2;
  // *worklin[1] = linear predictor of sigma_1 equation
  // *worktransformlin[1] = sigma_1;
  // *worklin[2] = linear predictor of mu_2 equation
  // *worktransformlin[2] = mu_2;
  // *worklin[3] = linear predictor of mu_1 equation
  // *worktransformlin[3] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double rho;
  double hilfs;

  if (*linpred <= -100) {
    rho  = -0.99995;
    hilfs = 100.05;
  }
  else if (*linpred >= 100) {
    rho  = 0.99995;
    hilfs = 100.05;
  }
  else {
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);
    hilfs = pow((1+pow((*linpred),2)),0.5);
  }


  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;


    double nu = oneminusrho2*(*linpred) - (*linpred)*( pow((((*response))-(*worktransformlin[3])),2)/pow((*worktransformlin[1]),2)
                                                      +  pow((((*response2p))-(*worktransformlin[2])),2)/pow((*worktransformlin[0]),2) )
                +(hilfs+rho*(*linpred))*( (((*response)-(*worktransformlin[3]))/((*worktransformlin[1])))*(((*response2p)-(*worktransformlin[2]))/((*worktransformlin[0]))) );



    *workingweight = 1-pow(rho2,2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response))-(*worklin[3])),2)/pow((*worktransformlin[1]),2) -
                                 2*rho*(((*response)-(*worklin[3]))/((*worktransformlin[1])))*(((*response2p)-(*worklin[2]))/((*worktransformlin[0])))
                                +  pow((((*response2p))-(*worklin[2])),2)/pow((*worktransformlin[0]),2) );

      }

  modify_worklin();

  }


void DISTR_bivnormal_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho): rhogit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivnormal_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin)/pow((1+pow((*worklin),2)),0.5);
    }

  }

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_bivnormal_sigma ----------------------
//------------------------------------------------------------------------------


DISTR_bivnormal_sigma::DISTR_bivnormal_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Bivariate Normal Distribution - sigma";

  pos = p;
  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_bivnormal_sigma::DISTR_bivnormal_sigma(const DISTR_bivnormal_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivnormal_sigma & DISTR_bivnormal_sigma::operator=(
                            const DISTR_bivnormal_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivnormal_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivnormal_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[pos + 1]));
  }
void DISTR_bivnormal_sigma::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivnormal_sigma::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }


double DISTR_bivnormal_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;


  if (counter==0)
    {
    set_worklin();
    }

  double sigma = exp((*linpred));

  // double rho2 = pow((*worktransformlin[0]),2);
  double oneminusrho2 = 1- pow((*worktransformlin[0]),2);
  double l;

//hier ist jetzt das problem, dass ich die aktuelle response gleichung brauche und die von dem mu2 oder sigma2

     l = -log(sigma) -(1/(2*oneminusrho2))*( pow((((*response))-(*worktransformlin[2])),2)/pow(sigma,2) -
                                 2*(*worktransformlin[0])*(((*response)-(*worktransformlin[2]))/(sigma))*(((*response2p)-(*worktransformlin[3]))/((*worktransformlin[1]))) );


  modify_worklin();

  return l;

  }

void DISTR_bivnormal_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;

  if (counter==0)
    {
    set_worklin();
    }

    double sigma = exp((*linpred));

    double rho2 = pow((*worktransformlin[0]),2);
    double oneminusrho2 = 1- rho2;


    double nu = -1 + (1/oneminusrho2)*(pow(((*response)-(*worklin[2])),2))/pow(sigma,2)
                - ((*worktransformlin[0])/oneminusrho2)*(((*response)-(*worktransformlin[2]))/(sigma))*(((*response2p)-(*worklin[3]))/((*worktransformlin[1])));



    *workingweight = 1+1/oneminusrho2;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -log(sigma) -(1/(2*oneminusrho2))*( pow((((*response))-(*worktransformlin[2])),2)/pow(sigma,2) -
                                 2*(*worktransformlin[0])*(((*response)-(*worktransformlin[2]))/(sigma))*(((*response2p)-(*worklin[3]))/((*worktransformlin[1]))) );

      }

  modify_worklin();

  }


void DISTR_bivnormal_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivnormal_sigma::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_bivnormal_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_bivnormal_mu::check_errors(void)
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


DISTR_bivnormal_mu::DISTR_bivnormal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  pos =p;
  family = "Bivariate Normal Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  if(o->forceIWLS)
    {
    updateIWLS = true;
    }
  else
    {
    updateIWLS = false;
    }
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  check_errors();
  }


DISTR_bivnormal_mu::DISTR_bivnormal_mu(const DISTR_bivnormal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivnormal_mu & DISTR_bivnormal_mu::operator=(
                            const DISTR_bivnormal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


void DISTR_bivnormal_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[2] = *response[4] = first component of two dimensional reponse
   // *linpred[0] = eta_rho
   // *linpred[1] = eta_sigma_2
   // *linpred[2] = eta_mu_2
   // *linpred[3] = eta_sigma_1
   // *linpred[4] = eta_mu_1

   if (*weight[4] == 0)
     *deviance=0;
   else
     {
     double rho = (*linpred[0])/pow(1+pow((*linpred[0]),2),0.5);
     double sigma_2 = exp(*linpred[1]);
     double mu_2 = (*linpred[3]);
     double sigma_1 = exp(*linpred[2]);
     double mu_1 = (*linpred[4]);
     double hilfs1 = 1-pow(rho,2);
     double l;

       l = -log(2*PI)-log(sigma_1)-log(sigma_2)-0.5*log(hilfs1)-(1/(2*hilfs1))*( pow((((*response[4]))-mu_1),2)/pow(sigma_1,2) -
                                                                                2*rho*(((*response[4])-mu_1)/(sigma_1))*(((*response[3])-mu_2)/(sigma_2))
                                                                                + pow((((*response[3]))-mu_2),2)/pow(sigma_2,2) );


    *deviance = -2*l;
    }

  }


double DISTR_bivnormal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivnormal_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[3 + pos]);
  }

 double DISTR_bivnormal_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_bivnormal_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
 //   double arg = ((*response[4])-(*param[4]))/pow((*param[2]),0.5) ;
 //   double u = gsl_cdf_ugaussian_P(arg);
 //   return u;
      return 0;
    }

void DISTR_bivnormal_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();
  weightp = weight.getV();

  }



void DISTR_bivnormal_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    weightp++;
    }

  }



double DISTR_bivnormal_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double rho2 = pow((*worktransformlin[0]),2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -(1/(2*oneminusrho2))*( pow((((*response))-mu),2)/pow((*worktransformlin[2]),2) -
                                 2*(*worktransformlin[0])*(((*response)-mu)/((*worktransformlin[2])))*(((*response2p)-(*worktransformlin[1]))/((*worktransformlin[3]))) );

  modify_worklin();

  return l;

  }


void DISTR_bivnormal_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);
    double rho2 = pow((*worktransformlin[0]),2);
   double oneminusrho2 = 1- rho2;


    double nu = (1/(oneminusrho2))*( (((*response))-mu)/pow((*worktransformlin[2]),2) -
                                 ((*worktransformlin[0])/(*worktransformlin[2]))*(((*response2p)-(*worktransformlin[1]))/((*worktransformlin[3]))) );

//    if(updateIWLS)
      *workingweight = 1/(oneminusrho2*pow((*worktransformlin[2]),2));
//    else
//      *workingweight = (*weightp)/(oneminusrho2*pow((*worktransformlin[2]),2));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -(1/(2*oneminusrho2))*( pow((((*response))-mu),2)/pow((*worktransformlin[2]),2) -
                                 2*(*worktransformlin[0])*(((*response)-mu)/((*worktransformlin[2])))*(((*response2p)-(*worktransformlin[1]))/((*worktransformlin[3]))) );

      }


  modify_worklin();

  }


void DISTR_bivnormal_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+3+pos]));
  }


void DISTR_bivnormal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }

/*void DISTR_bivnormal_mu::update(void)
  {

   unsigned i;

  double help;

  double * worktransformlinr;
  double * worktransformlins;
  double * workweight;

  worktransformlinr = distrp[0]->helpmat1.getV();
  worktransformlins = distrp[2]->helpmat1.getV();
  workweight = workingweight.getV();

  if(updateIWLS)
    {
    for (i=0;i<nrobs;i++,worktransformlinr++,worktransformlins++,workweight++)
      {
      *workweight = 1/((1 - pow((*worktransformlinr),2))*pow((*worktransformlins),2));
      }
    }
  else
    {
    weightp = weight.getV();
    for (i=0;i<nrobs;i++,worktransformlinr++,worktransformlins++,workweight++,weightp++)
      {
      *workweight = (*weightp)/((1 - pow((*worktransformlinr),2))*pow((*worktransformlins),2));
      }
    }
  }*/

void DISTR_bivnormal_mu::update_end(void)
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


//------------------------------------------------------------------------------
//------------------------- CLASS: #_rho -------------------------
//------------------------------------------------------------------------------


DISTR_bivprobit_rho::DISTR_bivprobit_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Bivariate Probit Distribution with Latent Variable Observation Model - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
    linpredminlimit=-100;
  linpredmaxlimit=100;
 // responseorig = response;

  }


DISTR_bivprobit_rho::DISTR_bivprobit_rho(const DISTR_bivprobit_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
 // responseorig = nd.responseorig;
//  response2 = nd.response2;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  response1p = nd.response1p;
  workingresponse1p = nd.workingresponse1p;
  }


const DISTR_bivprobit_rho & DISTR_bivprobit_rho::operator=(
                            const DISTR_bivprobit_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
 // responseorig = nd.responseorig;
 // response2 = nd.response2;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
    response1p = nd.response1p;
  workingresponse1p = nd.workingresponse1p;
  return *this;
  }


double DISTR_bivprobit_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivprobit_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = (*linpred[0]);
  *param = arg/pow(1+pow(arg,2),0.5);
  }

void DISTR_bivprobit_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

   response2p = workingresponse2p->getV();
   response1p = workingresponse1p->getV();
  }



void DISTR_bivprobit_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    response1p++;
    }

  }



double DISTR_bivprobit_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu_2 equation
  // *worktransformlin[0] = mu_2;
  // *worklin[1] = linear predictor of mu_1 equation
  // *worktransformlin[1] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }
  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response1p))-(*worktransformlin[1])),2) -
                                 2*rho*(((*response1p)-(*worktransformlin[1])))*(((*response2p)-(*worktransformlin[0])))
                                +  pow((((*response2p))-(*worktransformlin[0])),2) );


  modify_worklin();

  return l;

  }

void DISTR_bivprobit_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu_2 equation
  // *worktransformlin[0] = mu_2;
  // *worklin[1] = linear predictor of mu_1 equation
  // *worktransformlin[1] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double rho;
  double hilfs;

  if (*linpred <= -100) {
    rho  = -0.99995;
    hilfs = 100.05;
  }
  else if (*linpred >= 100) {
    rho  = 0.99995;
    hilfs = 100.05;
  }
  else {
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);
    hilfs = pow((1+pow((*linpred),2)),0.5);
  }

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;


    double nu = oneminusrho2*(*linpred) - (*linpred)*( pow((((*response1p))-(*worktransformlin[1])),2)
                                                      +  pow((((*response2p))-(*worktransformlin[0])),2) )
                +(hilfs+rho*(*linpred))*( (((*response1p)-(*worktransformlin[1])))*(((*response2p)-(*worktransformlin[0]))) );


     *workingweight = 1-pow(rho2,2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response1p))-(*worktransformlin[1])),2) -
                                 2*rho*(((*response1p)-(*worktransformlin[1])))*(((*response2p)-(*worktransformlin[0])))
                                +  pow((((*response2p))-(*worktransformlin[0])),2) );

      }


  modify_worklin();



  }


void DISTR_bivprobit_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho): fisher z-transformation\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivprobit_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin)/pow((1+pow((*worklin),2)),0.5);
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_bivprobit_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_bivprobit_mu::check_errors(void)
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
            if(((*workresp)!=0) && ((*workresp)!=1.0)) {
                errors=true;
                errormessages.push_back("ERROR: response has to be zero or one\n");
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


DISTR_bivprobit_mu::DISTR_bivprobit_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  pos =p;
  family = "Bivariate Probit Distribution with Latent Variable Observation Model - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  responseorig = response;
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  check_errors();
  }


DISTR_bivprobit_mu::DISTR_bivprobit_mu(const DISTR_bivprobit_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  responseorig = nd.responseorig;
//  response2 = nd.response2;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  }


const DISTR_bivprobit_mu & DISTR_bivprobit_mu::operator=(
                            const DISTR_bivprobit_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  responseorig = nd.responseorig;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  return *this;
  }


void DISTR_bivprobit_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[2] = *response[4] = first component of two dimensional reponse
   // *linpred[0] = eta_rho
   // *linpred[1] = eta_mu_2
   // *linpred[2] = eta_mu_1

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     double rho = (*linpred[0])/pow(1+pow((*linpred[0]),2),0.5);
     double mu_2 = (*linpred[1]);
     double mu_1 = (*linpred[2]);
     double hilfs1 = 1-pow(rho,2);
     double l;

       l = -log(2*PI)-0.5*log(hilfs1)-(1/(2*hilfs1))*( pow((((*response[2]))-mu_1),2) -
                                                                                2*rho*(((*response[2])-mu_1))*(((*response[1])-mu_2))
                                                                                + pow((((*response[1]))-mu_2),2) );


    *deviance = -2*l;
    }

  }


double DISTR_bivprobit_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivprobit_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[1 + pos]);
  }

void DISTR_bivprobit_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = workingresponse2p->getV();
  weightp = weight.getV();

  }



void DISTR_bivprobit_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    weightp++;
    }

  }



void DISTR_bivprobit_mu::update(void)
  {


  double * workresp = response.getV();
  double * workresporig =responseorig.getV();
  double * weightwork = weight.getV();

  double * worklin_current;
  if (linpred_current==1)
    worklin_current = linearpred1.getV();
  else
    worklin_current = linearpred2.getV();

  set_worklin();

  unsigned i;
  for(i=0;i<nrobs;i++,worklin_current++,workresp++,weightwork++,
           response2p++,workresporig++,worktransformlin[0]++,worklin[1]++)
    {

    if (*weightwork != 0)
      {
      if (*workresporig > 0)
        *workresp = trunc_normal2(0,20,*worklin_current+(*worktransformlin[0])*((*response2p)-(*worklin[1])),pow(1-pow(*worktransformlin[0],2),0.5));
      else
        *workresp = trunc_normal2(-20,0,*worklin_current+(*worktransformlin[0])*((*response2p)-(*worklin[1])),pow(1-pow(*worktransformlin[0],2),0.5));
      }

//          std::ofstream out;
//  // helpmat1.prettyPrint(out);
//    out.open ("C:\\tmp\\bivprobit.raw", std::ofstream::out | std::ofstream::app);
//    out << *workresp ;
//    out << " " ;
//    out << *workresporig ;
//    out << " " ;
//    out << *worklin[1] ;
//    out << " " ;
//    out << *worktransformlin[0] ;
//    out << " " ;
//    out << *worktransformlin[1] ;
//    out << " " ;
//    out << *response2p << endl;

    }




  }



double DISTR_bivprobit_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;


  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double rho2 = pow((*worktransformlin[0]),2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -(1/(2*oneminusrho2))*( pow((((*response))-mu),2) -
                                 2*(*worktransformlin[0])*(((*response)-mu))*(((*response2p)-(*worktransformlin[1]))) );

  modify_worklin();

  return l;

  }


void DISTR_bivprobit_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = (*linpred);
    double rho2 = pow((*worktransformlin[0]),2);
    double oneminusrho2 = 1- rho2;
 // cout << "mu equation y1: " << *response << endl;
 // cout << "mu equation y2: " << *response2p << endl;

    double nu = (1/(oneminusrho2))*( (((*response))-mu) -
                                 ((*worktransformlin[0]))*(((*response2p)-(*worktransformlin[1]))) );

//    if(updateIWLS)
      *workingweight = 1/(oneminusrho2);
//    else
//      *workingweight = (*weightp)/(oneminusrho2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -(1/(2*oneminusrho2))*( pow((((*response))-mu),2) -
                                 2*(*worktransformlin[0])*(((*response)-mu))*(((*response2p)-(*worktransformlin[1]))) );

      }

  modify_worklin();

  }


void DISTR_bivprobit_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+1+pos]));
  }


void DISTR_bivprobit_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivprobit_mu::update_end(void)
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

//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_bivprobit2_rho -------------------------
//------------------------------------------------------------------------------


DISTR_bivprobit2_rho::DISTR_bivprobit2_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Bivariate Probit Distribution with Binary Observation Model - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
    linpredminlimit=-100;
  linpredmaxlimit=100;
 // responseorig = response;

  }


DISTR_bivprobit2_rho::DISTR_bivprobit2_rho(const DISTR_bivprobit2_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
 // responseorig = nd.responseorig;
//  response2 = nd.response2;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  response1p = nd.response1p;
  workingresponse1p = nd.workingresponse1p;
  }


const DISTR_bivprobit2_rho & DISTR_bivprobit2_rho::operator=(
                            const DISTR_bivprobit2_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
 // responseorig = nd.responseorig;
 // response2 = nd.response2;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
    response1p = nd.response1p;
  workingresponse1p = nd.workingresponse1p;
  return *this;
  }


double DISTR_bivprobit2_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivprobit2_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = (*linpred[0]);
  *param = arg/pow(1+pow(arg,2),0.5);
  }

void DISTR_bivprobit2_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

   response2p = workingresponse2p->getV();
   response1p = workingresponse1p->getV();
  }



void DISTR_bivprobit2_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    response1p++;
    }

  }



double DISTR_bivprobit2_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu_2 equation
  // *worktransformlin[0] = mu_2;
  // *worklin[1] = linear predictor of mu_1 equation
  // *worktransformlin[1] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }

  vector<double> lower(3);
  vector<double> upper(3);
  lower[1] = -DBL_MAX;
  lower[2] = -DBL_MAX;
  upper[1] = (*worklin[1]);
  upper[2] = (*worklin[0]);
  double r;
  if (*linpred <= -100)
    r  = -0.99995;
  else if (*linpred >= 100)
    r  = 0.99995;
  else
    r = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double p = 0;
  double p1 = randnumbers::Phi2(upper[1]);
  double p2 = randnumbers::Phi2(upper[2]);
  double p11 = randnumbers::pbivn(lower[1], upper[1], lower[2], upper[2], r);
  double p10 = p1-p11;//randnumbers::pbivn(lower[1], upper[1], lower[2], -upper[2], -r);
  double p01 = p2-p11;//randnumbers::pbivn(lower[1], -upper[1], lower[2], upper[2], -r);
  double p00 = 1-p11-p10-p01;//randnumbers::pbivn(lower[1], upper[1], lower[2], upper[2], r);
  if((*response1p)>0)
    {
    if((*response2p)<=0)
      {
      p = p10;
      }
    else
      {
      p = p11;
      }
    }
  else
    {
    if((*response2p)<=0)
      {
      p = p00;
      }
    else
      {
      p = p01;
      }
    }
  double l = log(p);

  modify_worklin();

  return l;

  }

void DISTR_bivprobit2_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu_2 equation
  // *worktransformlin[0] = mu_2;
  // *worklin[1] = linear predictor of mu_1 equation
  // *worktransformlin[1] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  vector<double> lower(3);
  vector<double> upper(3);
  lower[1] = -DBL_MAX;
  lower[2] = -DBL_MAX;
  upper[1] = (*worklin[1]);
  upper[2] = (*worklin[0]);
  double r;
  if (*linpred <= -100)
    r  = -0.99995;
  else if (*linpred >= 100)
    r  = 0.99995;
  else
    r = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double dr = 1/pow((1+(*linpred)*(*linpred)), 1.5);
//  double ddr = -3*(*linpred)/pow((1+(*linpred)*(*linpred)), 2.5);
  double nucont;
  double dens = randnumbers::dbivn((*worklin[1]),(*worklin[0]),r);
  double p = 0;
  double p1 = randnumbers::Phi2(upper[1]);
  double p2 = randnumbers::Phi2(upper[2]);
  double p11 = randnumbers::pbivn(lower[1], upper[1], lower[2], upper[2], r);
  double p10 = p1-p11;//randnumbers::pbivn(lower[1], upper[1], lower[2], -upper[2], -r);
  double p01 = p2-p11;//randnumbers::pbivn(lower[1], -upper[1], lower[2], upper[2], -r);
  double p00 = 1-p11-p10-p01;//randnumbers::pbivn(lower[1], upper[1], lower[2], upper[2], r);
  if((*response1p)>0)
    {
    if((*response2p)<=0)
      {
      p = p10;
      nucont = -dens*dr/p10;
      }
    else
      {
      p = p11;
      nucont = dens*dr/p11;
      }
    }
  else
    {
    if((*response2p)<=0)
      {
      p = p00;
      nucont = dens*dr/p00;
      }
    else
      {
      p = p01;
      nucont = -dens*dr/p01;
      }
    }

  double nu = nucont;

  *workingweight = dens*dens*dr*dr*(1/p11 + 1/p10 + 1/p01 + 1/p00);

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    like +=  log(p);
    }

  modify_worklin();

  }


void DISTR_bivprobit2_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho): fisher z-transformation\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivprobit2_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin)/pow((1+pow((*worklin),2)),0.5);
    }

  }

//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_bivprobit2_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_bivprobit2_mu::check_errors(void)
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
            if(((*workresp)!=0) && ((*workresp)!=1.0)) {
                errors=true;
                errormessages.push_back("ERROR: response has to be zero or one\n");
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


DISTR_bivprobit2_mu::DISTR_bivprobit2_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  pos =p;
  family = "Bivariate Probit Distribution with Binary Observation Model - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  responseorig = response;
  check_errors();
  }


DISTR_bivprobit2_mu::DISTR_bivprobit2_mu(const DISTR_bivprobit2_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  responseorig = nd.responseorig;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  }


const DISTR_bivprobit2_mu & DISTR_bivprobit2_mu::operator=(
                            const DISTR_bivprobit2_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  responseorig = nd.responseorig;
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  return *this;
  }


void DISTR_bivprobit2_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[2] = *response[4] = first component of two dimensional reponse
   // *linpred[0] = eta_rho
   // *linpred[1] = eta_mu_2
   // *linpred[2] = eta_mu_1

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
     // double rho = (*linpred[0])/pow(1+pow((*linpred[0]),2),0.5);
     double mu_2 = (*linpred[1]);
     double mu_1 = (*linpred[2]);

     vector<double> lower(3);
    vector<double> upper(3);
    lower[1] = -DBL_MAX;
    lower[2] = -DBL_MAX;
    upper[1] = mu_1;
    upper[2] = mu_2;
    double r = (*worktransformlin[0]);
    double p = 0;
    double p1 = randnumbers::Phi2(upper[1]);
    double p2 = randnumbers::Phi2(upper[2]);
    double p11 = randnumbers::pbivn(lower[1], upper[1], lower[2], upper[2], r);
    double p10 = p1-p11;//randnumbers::pbivn(lower[1], upper[1], lower[2], -upper[2], -r);
    double p01 = p2-p11;//randnumbers::pbivn(lower[1], -upper[1], lower[2], upper[2], -r);
    double p00 = 1-p11-p10-p01;//randnumbers::pbivn(lower[1], upper[1], lower[2], upper[2], r);
    if((*response[2])>0)
      {
      if((*response[1])<=0)
        {
        p = p10;
        }
      else
        {
        p = p11;
        }
      }
    else
      {
      if((*response[1])<=0)
        {
        p = p00;
        }
      else
        {
        p = p01;
        }
      }

     double l;

     l = log(p);

     *deviance = -2*l;
     }

  }


double DISTR_bivprobit2_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivprobit2_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[1 + pos]);
  }

void DISTR_bivprobit2_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

    response2p = workingresponse2p->getV();

  }



void DISTR_bivprobit2_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    response2p++;

  }

/*void DISTR_bivprobit2_mu::update(void)
  {


  double * workresp = response.getV();
  double * workresporig =responseorig.getV();
  double * weightwork = weight.getV();

  double * worklin_current;
  if (linpred_current==1)
    worklin_current = linearpred1.getV();
  else
    worklin_current = linearpred2.getV();

  set_worklin();

  unsigned i;
  for(i=0;i<nrobs;i++,worklin_current++,workresp++,weightwork++,
           response2p++,workresporig++,worktransformlin[0]++,worklin[1]++)
    {

    if (*weightwork != 0)
      {
      if (*workresporig > 0)
        *workresp = trunc_normal2(0,20,*worklin_current+(*worktransformlin[0])*((*response2p)-(*worklin[1])),pow(1-pow(*worktransformlin[0],2),0.5));
      else
        *workresp = trunc_normal2(-20,0,*worklin_current+(*worktransformlin[0])*((*response2p)-(*worklin[1])),pow(1-pow(*worktransformlin[0],2),0.5));
      }

    }
  }*/

double DISTR_bivprobit2_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  if (counter==0)
    {
    set_worklin();
    }

  vector<double> lower(3);
  vector<double> upper(3);
  lower[1] = -DBL_MAX;
  lower[2] = -DBL_MAX;
  upper[1] = (*linpred);
  upper[2] = (*worklin[1]);
  double r = (*worktransformlin[0]);
  double p = 0;
  double p1 = randnumbers::Phi2(upper[1]);
  double p2 = randnumbers::Phi2(upper[2]);
  double p11 = randnumbers::pbivn(lower[1], upper[1], lower[2], upper[2], r);
  double p10 = p1-p11;//randnumbers::pbivn(lower[1], upper[1], lower[2], -upper[2], -r);
  double p01 = p2-p11;//randnumbers::pbivn(lower[1], -upper[1], lower[2], upper[2], -r);
  double p00 = 1-p11-p10-p01;//randnumbers::pbivn(lower[1], upper[1], lower[2], upper[2], r);
  if((*response)>0)
    {
    if((*response2p)<=0)
      {
      p = p10;
      }
    else
      {
      p = p11;
      }
    }
  else
    {
    if((*response2p)<=0)
      {
      p = p00;
      }
    else
      {
      p = p01;
      }
    }
  double l = log(p);

  modify_worklin();

  return l;

  }


void DISTR_bivprobit2_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;

  if (counter==0)
    {
    set_worklin();
    }

  vector<double> lower(3);
  vector<double> upper(3);
  lower[1] = -DBL_MAX;
  lower[2] = -DBL_MAX;
  upper[1] = (*linpred);
  upper[2] = (*worklin[1]);
  double r = (*worktransformlin[0]);
  double p;
  double arg = ((*worklin[1])-r*(*linpred))/pow((1-r*r),0.5);
  double nucont;
  double prop = randnumbers::Phi2(arg);
  double dens = randnumbers::phi(upper[1]);
  double p1 = randnumbers::Phi2(upper[1]);
  double p2 = randnumbers::Phi2(upper[2]);
  double p11 = randnumbers::pbivn(lower[1], upper[1], lower[2], upper[2], r);
  double p10 = p1-p11;//randnumbers::pbivn(lower[1], upper[1], lower[2], -upper[2], -r);
  double p01 = p2-p11;//randnumbers::pbivn(lower[1], -upper[1], lower[2], upper[2], -r);
  double p00 = 1-p11-p10-p01;//randnumbers::pbivn(lower[1], upper[1], lower[2], upper[2], r);
  if((*response)>0)
    {
    if((*response2p)<=0)
      {
      p = p10;
      nucont = dens*(1-prop)/p10;
      }
    else
      {
      p = p11;
      nucont = dens*prop/p11;
      }
    }
  else
    {
    if((*response2p)<=0)
      {
      p = p00;
      nucont = -dens*(1-prop)/p00;
      }
    else
      {
      p = p01;
      nucont = -dens*prop/p01;
      }
    }


  double nu = nucont;

  *workingweight = -dens*dens * ( prop*prop*(-1/p11 - 1/p01) + (1-prop)*(1-prop)*(-1/p10 - 1/p00) );

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    like += log(p);
    }

  modify_worklin();

  }


void DISTR_bivprobit2_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+1+pos]));
  }


void DISTR_bivprobit2_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivprobit2_mu::update_end(void)
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


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_bivlogit_or -------------------------
//------------------------------------------------------------------------------


DISTR_bivlogit_or::DISTR_bivlogit_or(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Bivariate Logit Distribution - Odds Ratio";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "oddsratio";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  datamatrix d(nrobs,1,0.00001);

  if (linpred_current==1)
    linearpred1.plus(d);
  else
    linearpred2.plus(d);


  }


DISTR_bivlogit_or::DISTR_bivlogit_or(const DISTR_bivlogit_or & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivlogit_or & DISTR_bivlogit_or::operator=(
                            const DISTR_bivlogit_or & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_bivlogit_or::get_intercept_start(void)
  {
  return log(response.mean(0));
  }

void DISTR_bivlogit_or::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[0]));
  }

void DISTR_bivlogit_or::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivlogit_or::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_bivlogit_or::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu_2 equation
  // *worktransformlin[0] = mu_2;
  // *worklin[1] = linear predictor of mu_1 equation
  // *worktransformlin[1] = mu_1;


  if (counter==0)
    {
    set_worklin();
    }
  double odds = exp((*linpred));
  double p1 = (*worktransformlin[1]);
  double psiminone = odds - 1;
  double hilfs1 = 1 + (p1 + (*worktransformlin[0]))*psiminone;
  double hilfs2 = -4*odds*psiminone*p1*(*worktransformlin[0]);
  double p11 = 0.5*pow(psiminone, -1)*( hilfs1 - pow((pow(hilfs1,2) + hilfs2), 0.5));
  if(odds == 1)
  {
      p11 = p1*(*worktransformlin[0]);
  }

  double l;

  if(((*response) == 0) && ((*response2p) == 0)) {
    l = log(1 + p11 - p1 - (*worktransformlin[0]));
  }
  else if(((*response) == 0) && ((*response2p) == 1)) {
    l = log( (*worktransformlin[0]) - p11 );
  }
  else if(((*response) == 1) && ((*response2p) == 0)) {
    l = log( p1 - p11 );
  }
  else {
    l = log( p11 );
  }


  modify_worklin();

  return l;

  }

void DISTR_bivlogit_or::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma_2 equation
  // *worktransformlin[0] = sigma_2;
  // *worklin[1] = linear predictor of sigma_1 equation
  // *worktransformlin[1] = sigma_1;
  // *worklin[2] = linear predictor of mu_2 equation
  // *worktransformlin[2] = mu_2;
  // *worklin[3] = linear predictor of mu_1 equation
  // *worktransformlin[3] = mu_1;

  if (counter==0)
    {
    set_worklin();
    }

  double odds = exp((*linpred));
  double p1 = (*worktransformlin[1]);
  double p2 = (*worktransformlin[0]);
  double psiminone = odds - 1;
  double hilfs1 = 1 + (p1 + p2)*psiminone;
  double hilfs2 = -4*odds*psiminone*p1*p2;
  double p11 = 0.5*pow(psiminone, -1)*( hilfs1 - pow((pow(hilfs1,2) + hilfs2), 0.5));
  if(odds == 1)
  {
      p11 = p1*p2;
  }

  double dp11 = -p11*odds/psiminone + 0.5*pow(psiminone, -1)*( odds*(p1+p2) - ( hilfs1*odds*(p1+p2) - 2*(2*pow(odds,2)-odds)*p2*p1 )/( pow((pow(hilfs1, 2) + hilfs2), 0.5) ) );

  if(odds == 1)
  {
      dp11 = 0;
  }

  double nu;

  if(((*response) == 0) && ((*response2p) == 0)) {
    nu = ( dp11 )/( 1+ p11 - p1- p2 );
  }
  else if(((*response) == 0) && ((*response2p) == 1)) {
    nu = ( -dp11 )/( p2 - p11 );
  }
  else if(((*response) == 1) && ((*response2p) == 0)) {
    nu = ( - dp11 )/( p1 - p11 );
  }
  else {
    nu = ( dp11 )/( p11 );
  }

  *workingweight = pow(dp11, 2)/((1 + p11 -p1 -p2)*( p2 - p11 )) + pow(dp11, 2)/(p11*( p1 - p11 ));

  *workingresponse = *linpred + nu/(*workingweight);

  if(odds == 1)
  {
      *workingresponse = *linpred;
  }

    if (compute_like)
      {

        if(((*response) == 0) && ((*response2p) == 0)) {
            like += log(1 + p11 - p1 - p2);
        }
        else if(((*response) == 0) && ((*response2p) == 1)) {
            like += log( p2 - p11 );
        }
        else if(((*response) == 1) && ((*response2p) == 0)) {
            like += log( p1 - p11 );
        }
        else {
            like += log( p11 );
        }

      }

  modify_worklin();

  }


void DISTR_bivlogit_or::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (odds ratio): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivlogit_or::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
      *pmu = exp((*worklin));
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_bivlogit_mu ------------------------
//------------------------------------------------------------------------------
void DISTR_bivlogit_mu::check_errors(void)
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

            if (((*workresp)!= 0.0) && ((*workresp)!= 1.0) )
            {
                errors=true;
                errormessages.push_back("ERROR: response has to be equal to zero or one\n");
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


DISTR_bivlogit_mu::DISTR_bivlogit_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           unsigned & p,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  pos =p;
  family = "Bivariate Logit Distribution - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_bivlogit_mu::DISTR_bivlogit_mu(const DISTR_bivlogit_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  }


const DISTR_bivlogit_mu & DISTR_bivlogit_mu::operator=(
                            const DISTR_bivlogit_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  pos = nd.pos;
  response2 = nd.response2;
  response2p = nd.response2p;
  return *this;
  }


void DISTR_bivlogit_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[2] = *response[4] = first component of two dimensional reponse
   // *linpred[0] = eta_or
   // *linpred[1] = eta_mu_2

   if (*weight[2] == 0)
     *deviance=0;
   else
     {

    double odds = exp((*linpred[0]));
    double p1 = exp( (*linpred[2]) );
    p1 = p1 / (1+p1);
    double p2 = exp( (*linpred[1]) );
    p2 = p2 / (1+p2);
    double psiminone = odds - 1;
    double hilfs1 = 1 + (p1 +p2)*psiminone;
    double hilfs2 = -4*odds*psiminone*p1*p2;
    double p11 = 0.5*pow(psiminone, -1)*( hilfs1 - pow((pow(hilfs1,2) + hilfs2), 0.5));

     double l = 0;

    if(((*response[2]) == 0) && ((*response[1]) == 0)) {
        l += log(1 + p11 - p1 - p2);
    }
    else if(((*response[2]) == 0) && ((*response[1]) == 1)) {
        l += log( p2 - p11 );
    }
    else if(((*response[2]) == 1) && ((*response[1]) == 0)) {
        l += log( p1 - p11 );
    }
    else {
        l += log( p11 );
    }


    *deviance = -2*l;
    }

  }


double DISTR_bivlogit_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_bivlogit_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  double el = exp(*linpred[1 + pos]);
  *param = el/(1+el);
  }

void DISTR_bivlogit_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_bivlogit_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_bivlogit_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of or equation
  // *worktransformlin[0] = or;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;


  if (counter==0)
    {
    set_worklin();
    }
  double el = exp((*linpred));
  double p1 = el/(1+el);
  double psiminone = (*worktransformlin[0]) - 1;
  double hilfs1 = 1 + (p1 + (*worktransformlin[1]))*psiminone;
  double hilfs2 = -4*(*worktransformlin[0])*psiminone*p1*(*worktransformlin[1]);
  double p11 = 0.5*pow(psiminone, -1)*( hilfs1 - pow((pow(hilfs1,2) + hilfs2), 0.5));
  if((*worktransformlin[0]) == 1)
  {
      p11 = p1*(*worktransformlin[1]);
  }

  double l;

  if(((*response) == 0) && ((*response2p) == 0)) {
    l = log(1 + p11 - p1 - (*worktransformlin[1]));
  }
  else if(((*response) == 0) && ((*response2p) == 1)) {
    l = log( (*worktransformlin[1]) - p11 );
  }
  else if(((*response) == 1) && ((*response2p) == 0)) {
    l = log( p1 - p11 );
  }
  else {
    l = log( p11 );
  }


  modify_worklin();

  return l;

  }


void DISTR_bivlogit_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of rho equation
  // *worktransformlin[0] = rho;
  // *worklin[1] = linear predictor of mu_2 equation
  // *worktransformlin[1] = mu_2;
  // *worklin[2] = linear predictor of sigma_1 equation
  // *worktransformlin[2] = sigma_1;
  // *worklin[3] = linear predictor of sigma_2 equation
  // *worktransformlin[3] = sigma_2;

  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

  double el = exp((*linpred));
  double p1 = el/(1+el);
  double psiminone = (*worktransformlin[0]) - 1;
  double hilfs1 = 1 + (p1 + (*worktransformlin[1]))*psiminone;
  double hilfs2 = -4*(*worktransformlin[0])*psiminone*p1*(*worktransformlin[1]);
  double p11 = 0.5*pow(psiminone, -1)*( hilfs1 - pow((pow(hilfs1,2) + hilfs2), 0.5));
  if((*worktransformlin[0]) == 1)
  {
      p11 = p1*(*worktransformlin[1]);
  }
  double p2 = (*worktransformlin[1]);
  double p1oneminusp1 = p1*(1-p1);
  double dp11 = 0.5*p1oneminusp1*( 1 - ( hilfs1 - 2*(*worktransformlin[0])*p2 )/( pow((pow(hilfs1, 2) + hilfs2), 0.5) ) );

  if((*worktransformlin[0]) == 1)
  {
      dp11 = p1oneminusp1*p2;
  }

  double nu;

  if(((*response) == 0) && ((*response2p) == 0)) {
    nu = ( dp11 - p1oneminusp1 )/( 1+ p11 - p1- p2 );
  }
  else if(((*response) == 0) && ((*response2p) == 1)) {
    nu = ( -dp11 )/( p2 - p11 );
  }
  else if(((*response) == 1) && ((*response2p) == 0)) {
    nu = ( p1oneminusp1 - dp11 )/( p1 - p11 );
  }
  else {
    nu = ( dp11 )/( p11 );
  }

  *workingweight = pow(dp11 - p1oneminusp1, 2)/(1 + p11 -p1 -p2) + pow(dp11, 2)/(p11*( p2 - p11 )) + pow(p1oneminusp1 - dp11, 2)/( p1 - p11 );

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
  {

    if(((*response) == 0) && ((*response2p) == 0)) {
        like += log(1 + p11 - p1 - p2);
    }
    else if(((*response) == 0) && ((*response2p) == 1)) {
        like += log( p2 - p11 );
    }
    else if(((*response) == 1) && ((*response2p) == 0)) {
        like += log( p1 - p11 );
    }
    else {
        like += log( p11 );
    }

  }


  modify_worklin();

  }


void DISTR_bivlogit_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double p1 = exp( (*linpred[predstart_mumult+2]) );
  p1 = p1 / (1+p1);
  double p2 = exp( (*linpred[predstart_mumult+1]) );
  p2 = p2 / (1+p2);

  double psi = exp((*linpred[predstart_mumult]));

  double a = 1 + ( p1 + p2 ) * ( psi - 1 );
  double b = - 4 * psi * ( psi - 1 ) * p1 * p2;

  if(psi == 1) {
    *mu = ( p1 * p2 );
  } else {
    *mu = 0.5 * pow(( psi - 1), -1) * ( a - pow(( pow(a, 2) + b ), 0.5) );
  }

  }


void DISTR_bivlogit_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_bivlogit_mu::update_end(void)
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
     double el = exp((*worklin));
    *pmu = el/(1+el);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_BCCG_nu -------------------------------
//------------------------------------------------------------------------------


DISTR_BCCG_nu::DISTR_BCCG_nu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "BCCG Distribution - nu";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "nu";
    linpredminlimit=-100;
  linpredmaxlimit=150;
  datamatrix d(nrobs,1,0.0001);

  if (linpred_current==1)
    linearpred1.plus(d);
  else
    linearpred2.plus(d);

  }


DISTR_BCCG_nu::DISTR_BCCG_nu(const DISTR_BCCG_nu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_BCCG_nu & DISTR_BCCG_nu::operator=(
                            const DISTR_BCCG_nu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }

  void DISTR_BCCG_nu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[0]);
  }

double DISTR_BCCG_nu::get_intercept_start(void)
  {
  return 0;
  }


double DISTR_BCCG_nu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = exp(eta_sigma);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

	 double nup = (*linpred);
	 double hilfs2 = (*response)/(*worktransformlin[1]);
     double hilfs = pow(hilfs2, nup);


     double l;

       l =   (nup)*log((*response)) - nup*log((*worktransformlin[1]))
            - ((1)/(2*pow((*worktransformlin[0])*nup, 2)))*pow((hilfs-1) ,2);// - log(randnumbers::Phi2(1 / (abs(nup) * (*worktransformlin[0]))));


  modify_worklin();

  return l;

  }

void DISTR_BCCG_nu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = exp(eta_sigma);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);


  if (counter==0)
    {
    set_worklin();
    }

	 double nup = (*linpred);
     double hilfs2 = (*response)/(*worktransformlin[1]);
     double hilfs = pow(hilfs2, nup);
     double nu2sig2 =  pow(nup * (*worktransformlin[0]), 2);
     // double arg = abs(nup) * (*worktransformlin[0]);

     double nu = log((*response)) - log((*worktransformlin[1])) + (1/(pow(nup,3)*pow((*worktransformlin[0]),2)))*pow(hilfs-1,2) -
                    (1/(nu2sig2))*(hilfs-1)*hilfs*log(hilfs2);// + randnumbers::phi(1 / arg) * randnumbers::sgn(nup) / (pow(nup, 2) * randnumbers::Phi2(1 / arg));


    *workingweight = 3 * pow(hilfs - 1, 2) / (pow(nup,2) * nu2sig2) - 4 * (hilfs - 1) * hilfs * log(hilfs2) / (nup * nu2sig2) +
                     2 * pow(hilfs, 2) * pow(log(hilfs2), 2) / nu2sig2 -
                     hilfs * pow(log(hilfs2), 2) / nu2sig2;// + pow(randnumbers::phi(1 / arg) * randnumbers::sgn(nup) / (pow(nup, 2) * randnumbers::Phi2(1 / arg)), 2) +
                   //  randnumbers::phi(1 / arg) / (arg * pow(nup, 2) * randnumbers::Phi2(1 / arg)) +
                   //  randnumbers::phi(1 / arg) * randnumbers::sgn(nup) / (pow(nup, 3) * randnumbers::Phi2(1 / arg));

    if((*workingweight<=0))
        *workingweight = 0.0001;


       *workingresponse = *linpred + nu/(*workingweight);


    if (compute_like)
      {

        like +=  (nup)*log((*response)) - nup*log((*worktransformlin[1]))
                -((1)/(2*pow((*worktransformlin[0])*nup, 2)))*pow((hilfs-1) ,2);// - log(randnumbers::Phi2(1 / (abs(nup) * (*worktransformlin[0]))));
      }

 /*     std::ofstream out;
  // helpmat1.prettyPrint(out);
    out.open ("C:\\tmp\\res.raw", std::ofstream::out | std::ofstream::app);
    out << nu ;
    out << " " ;
    out << *workingresponse ;
    out << " " ;
    out << *workingweight ;
    out << " " ;
    out << *worktransformlin[0] ;
    out << " " ;
    out << *worktransformlin[1] ;
    out << " " ;
    out << *linpred ;
    out << " " ;
    out << nup  << endl;
*/

  modify_worklin();

  }


void DISTR_BCCG_nu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (nu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_BCCG_nu::update_end(void)
  {

  // helpmat1 stores tau

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
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_BCCG_sigma ----------------------------
//------------------------------------------------------------------------------


DISTR_BCCG_sigma::DISTR_BCCG_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "BCCG Distribution - sigma";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_BCCG_sigma::DISTR_BCCG_sigma(const DISTR_BCCG_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_BCCG_sigma & DISTR_BCCG_sigma::operator=(
                            const DISTR_BCCG_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_BCCG_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_BCCG_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[1]));
  }
double DISTR_BCCG_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

     double sig = exp(*linpred);
     double hilfs2 = (*response)/(*worktransformlin[1]);
     double hilfs = pow(hilfs2, (*worktransformlin[0]));

     double l;

       l =  - log(sig) - ((1)/(2*pow(sig*(*worktransformlin[0]), 2)))*pow((hilfs-1) ,2);// - log(randnumbers::Phi2(1 / (abs((*worktransformlin[0])) * sig)));


  modify_worklin();

  return l;

  }

void DISTR_BCCG_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = exp(eta_nu);
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }

     double sig = exp(*linpred);
     double hilfs2 = (*response)/(*worktransformlin[1]);
     double hilfs = pow(hilfs2, (*worktransformlin[0]));
     // double arg = sig * abs((*worktransformlin[0]));


    double nu = -1 +  ((1)/(pow(sig*(*worktransformlin[0]), 2)))*pow((hilfs-1) ,2);// + randnumbers::phi(1 / arg) / (arg * randnumbers::Phi2(1 / arg));

  //  *workingweight = ((2)/(pow(sig*(*worktransformlin[0]), 2)))*pow((hilfs-1) ,2);

    *workingweight = 2;// + pow(randnumbers::phi(1 / arg) / (arg * randnumbers::Phi2(1 / arg)), 2) + 4 * abs(*worktransformlin[0]) * sig / pow(2 * PI, 0.5)
                      // + randnumbers::phi(1 / arg) / (arg * randnumbers::Phi2(1 / arg)) - randnumbers::phi(1 / arg) / (pow(arg, 3) * randnumbers::Phi2(1 / arg));

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=   - log(sig) - ((1)/(2*pow(sig*(*worktransformlin[0]), 2)))*pow((hilfs-1) ,2);// - log(randnumbers::Phi2(1 / (abs((*worktransformlin[0])) * sig)));

      }

  modify_worklin();

  }


void DISTR_BCCG_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_BCCG_sigma::update_end(void)
  {

  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_BCCG_mu -----------------------------
//------------------------------------------------------------------------------


void DISTR_BCCG_mu::check_errors(void)
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


DISTR_BCCG_mu::DISTR_BCCG_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
    family = "BCCG Distribution - mu";
    outpredictor = true;
    outexpectation = true;
    predictor_name = "mu";
    linpredminlimit=-10;
    linpredmaxlimit=15;
	check_errors();
  }


DISTR_BCCG_mu::DISTR_BCCG_mu(const DISTR_BCCG_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_BCCG_mu & DISTR_BCCG_mu::operator=(
                            const DISTR_BCCG_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_BCCG_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_nu
   // *linpred[1] = eta_sigma
   // *linpred[2] = eta_mu

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
	 double nup = (*linpred[0]);
     double sig = exp(*linpred[1]);
     double mu = exp(*linpred[2]);
     double hilfs = pow((*response[2])/mu, nup);
     // double hilfs2 = (*response[2])/mu;

     double l;

       l =  -0.5*log(2*(PI)) - log(sig) + (nup-1)*log((*response[2])) - nup*log(mu)
            - ((1)/(2*pow(sig*nup, 2)))*pow((hilfs-1) ,2);// - log(randnumbers::Phi2(1 / (abs(nup) * sig)));



    *deviance = -2*l;
    }

  }


double DISTR_BCCG_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_BCCG_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[2]));
  }

 double DISTR_BCCG_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_BCCG_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
 //   double x = ((*response[2])-(*param[2]))/pow((*param[1]),0.5);
 //   double u = gsl_cdf_tdist_P(x, (*param[0]));
 //   return u;
    return 0;
    }

double DISTR_BCCG_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = nu;
  // *worklin[1] = linear predictor of sigma equation
  // *worktransformlin[1] = sigma;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = exp(*linpred);
  double hilfs2 = (*response)/mu;
  double hilfs = pow(hilfs2, (*worktransformlin[0]));

  double l;

       l = - (*worktransformlin[0])*log(mu) -  ((1)/(2*pow((*worktransformlin[1])*(*worktransformlin[0]), 2)))*pow((hilfs-1) ,2);


  modify_worklin();

  return l;

  }


void DISTR_BCCG_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of nu equation
  // *worktransformlin[0] = nu;
  // *worklin[1] = linear predictor of sigma equation
  // *worktransformlin[1] = sigma;


  // ofstream out("d:\\_sicher\\papzip\\results\\helpmat1.raw");
  // helpmat1.prettyPrint(out);
  // for (i=0;i<helpmat1.rows();i++)
  //   out << helpmat1(i,0) << endl;

  if (counter==0)
    {
    set_worklin();
    }

    double mu = exp(*linpred);
    double hilfs2 = (*response)/mu;
    double hilfs = pow(hilfs2, (*worktransformlin[0]));

    double nu = -(*worktransformlin[0]) +  ((1)/((*worktransformlin[0])*pow((*worktransformlin[1]),2)))*(hilfs-1)*hilfs;


  /*  if ((*worktransformlin[0])==0)
    {
        *workingweight = ((1)/(pow((*worktransformlin[1]),2)));
    }
    else
    {
        *workingweight = ((1)/(pow((*worktransformlin[1]),2)))*(2*pow(hilfs,2)-hilfs);
    } */

    *workingweight = 2 * pow((*worktransformlin[0]) ,2) + 2/pow((*worktransformlin[1]),2);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - (*worktransformlin[0])*log(mu) -  ((1)/(2*pow((*worktransformlin[1])*(*worktransformlin[0]), 2)))*pow((hilfs-1) ,2);


      }


  modify_worklin();

  }


void DISTR_BCCG_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

       *mu = 0;


  }


void DISTR_BCCG_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_BCCG_mu::update_end(void)
  {


  // helpmat1 stores exp(eta_a)

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa0_sigma_v --------------------------
//------------------------------------------------------------------------------


DISTR_sfa0_sigma_v::DISTR_sfa0_sigma_v(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Stochastic Frontier Analysis - sigma_v";
    outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma_v";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_sfa0_sigma_v::DISTR_sfa0_sigma_v(const DISTR_sfa0_sigma_v & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa0_sigma_v & DISTR_sfa0_sigma_v::operator=(
                            const DISTR_sfa0_sigma_v & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa0_sigma_v::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa0_sigma_v::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[0]));
  }

double DISTR_sfa0_sigma_v::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma_u equation
  // *worktransformlin[0] = exp(eta_sigma_u);
  // *worklin[1] = linear predictor of mu_y equation
  // *worktransformlin[1] =(eta_mu_y);

  if (counter==0)
    {
    set_worklin();
    }

  double sigv = exp((*linpred));
    double sig2 = pow(sigv, 2) + pow((*worktransformlin[0]), 2);
    double epsi = (*response) - (*worktransformlin[1]);
    double proparg2 = randnumbers::Phi2(-epsi*(*worktransformlin[0])/(pow(sig2,0.5)*sigv));

  double l;

     l = - 0.5*log(sig2) - 0.5*pow(epsi, 2)/sig2 + log(proparg2);


  modify_worklin();

  return l;

  }

void DISTR_sfa0_sigma_v::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma_u equation
  // *worktransformlin[0] = exp(eta_sigma_u);
  // *worklin[1] = linear predictor of mu_y equation
  // *worktransformlin[1] =(eta_mu_y);


  if (counter==0)
    {
    set_worklin();
    }

    double sigv = exp((*linpred));
    double sig2 = pow(sigv, 2) + pow((*worktransformlin[0]), 2);
    double epsi = (*response) - (*worktransformlin[1]);
    double arg2 = -epsi*(*worktransformlin[0])/(pow(sig2,0.5)*sigv);
    double darg2 = -arg2 + epsi*sigv*(*worktransformlin[0])/(pow(sig2,1.5));
    double ddarg2 = arg2 - 3*epsi*pow(sigv,3)*(*worktransformlin[0]) /(pow(sig2,2.5)) ;
    double proparg2 = randnumbers::Phi2(arg2);
    double rat = randnumbers::phi(arg2)/proparg2;

    double nu = -pow(sigv, 2)/sig2 + pow(sigv*epsi, 2)/pow(sig2, 2) + rat*darg2;

    *workingweight = 2*pow(sigv*(*worktransformlin[0]), 2)/pow(sig2, 2) + 2*pow(sigv*epsi, 2)*(pow(sigv, 2) - pow((*worktransformlin[0]), 2))/pow(sig2, 3) -
                     rat*(-arg2)*pow(darg2, 2) - rat*ddarg2 + pow(rat*darg2, 2);

    if((*workingweight)<=0)
    {
        *workingweight = 0.0001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += - 0.5*log(sig2) - 0.5*pow(epsi, 2)/sig2 + log(proparg2);

      }

  modify_worklin();

  }


void DISTR_sfa0_sigma_v::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma_v): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa0_sigma_v::update_end(void)
  {

  // helpmat1 stores tau

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa0_sigma_u --------------------------
//------------------------------------------------------------------------------


DISTR_sfa0_sigma_u::DISTR_sfa0_sigma_u(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Stochastic Frontier Analysis - sigma_u";
    outpredictor = true;
  outexpectation = true;
  predictor_name = "sigma_u";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  }


DISTR_sfa0_sigma_u::DISTR_sfa0_sigma_u(const DISTR_sfa0_sigma_u & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa0_sigma_u & DISTR_sfa0_sigma_u::operator=(
                            const DISTR_sfa0_sigma_u & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa0_sigma_u::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa0_sigma_u::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[1]));
  }

double DISTR_sfa0_sigma_u::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma_v equation
  // *worktransformlin[0] = exp(eta_sigma_v);
  // *worklin[1] = linear predictor of mu_y equation
  // *worktransformlin[1] =(eta_mu_y);

  if (counter==0)
    {
    set_worklin();
    }

  double sigu = exp((*linpred));
  double sig2 = pow(sigu, 2) + pow((*worktransformlin[0]), 2);
  double epsi = (*response) - (*worktransformlin[1]);
  double proparg2 = randnumbers::Phi2(-epsi*sigu/(pow(sig2,0.5)*(*worktransformlin[0])));

  double l;

     l = -0.5*log(sig2) - 0.5*pow(epsi,2)/sig2 + log(proparg2) ;


  modify_worklin();

  return l;

  }

void DISTR_sfa0_sigma_u::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma_v equation
  // *worktransformlin[0] = exp(eta_sigma_v);
  // *worklin[1] = linear predictor of mu_y equation
  // *worktransformlin[1] =(eta_mu_y);

  if (counter==0)
    {
    set_worklin();
    }

    double sigu = exp((*linpred));
    double sig2 = pow(sigu, 2) + pow((*worktransformlin[0]), 2);
    double epsi = (*response) - (*worktransformlin[1]);
    double arg2 = -epsi*sigu/(pow(sig2,0.5)*(*worktransformlin[0]));
    double darg2 = arg2 + epsi*pow(sigu,3)/(pow(sig2,1.5)*(*worktransformlin[0]));
    double ddarg2 = darg2 + 3*epsi*pow(sigu,3)*(*worktransformlin[0])/(pow(sig2,2.5)) ;
    double proparg2 = randnumbers::Phi2(arg2);
    double rat = randnumbers::phi(arg2)/proparg2;

    double nu = -pow(sigu,2)/sig2 + pow(epsi*sigu,2)/pow(sig2, 2) + rat*darg2;



    *workingweight = 2*pow(sigu*(*worktransformlin[0]) ,2)/pow(sig2, 2) + 2*pow(epsi*sigu, 2)*(pow(sigu, 2) - pow((*worktransformlin[0]),2))/pow(sig2, 3) -
                     rat*(-arg2)*pow(darg2, 2) - rat*ddarg2 + pow(rat*darg2, 2);

    if((*workingweight)<=0)
    {
        *workingweight = 0.0001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*log(sig2) - 0.5*pow(epsi,2)/sig2 + log(proparg2) ;

      }

  modify_worklin();

  }

void DISTR_sfa0_sigma_u::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double muu = 0;
  double sigu2 = pow(exp((*linpred[predstart_mumult+1])), 2);
  double sigv2 = pow(exp((*linpred[predstart_mumult])), 2);
  double sig2 = sigu2 + sigv2;
  double ga = sigu2/sig2;
  double sigast = pow(ga * (1 - ga)*sig2, 0.5);
  double epsi = (*response[2]) - (*linpred[predstart_mumult+2]);
  double muast = (1 - ga) * muu - ga * epsi;

  *mu = muast + sigast * randnumbers::phi(muast / sigast) / randnumbers::Phi2(muast / sigast);

  }


void DISTR_sfa0_sigma_u::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma_u): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa0_sigma_u::update_end(void)
  {
  // helpmat1 stores sigma

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_sfa0_mu_y ---------------------------
//------------------------------------------------------------------------------

void DISTR_sfa0_mu_y::check_errors(void)
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

DISTR_sfa0_mu_y::DISTR_sfa0_mu_y(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w)
  {
  family = "Stochastic Frontier Analysis - mu_y";
    outpredictor = true;
  outexpectation = true;
//  outcondE = false;
  predictor_name = "mu_y";
//   linpredminlimit=-10;
//  linpredmaxlimit=15;
  check_errors();
  }


DISTR_sfa0_mu_y::DISTR_sfa0_mu_y(const DISTR_sfa0_mu_y & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa0_mu_y & DISTR_sfa0_mu_y::operator=(
                            const DISTR_sfa0_mu_y & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_sfa0_mu_y::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_tau
   // *linpred[1] = eta_sigma
   // *linpred[2] = eta_mu

   if (*weight[2] == 0)
     *deviance=0;
   else
     {
	 double sigv = exp(*linpred[0]);
     double sigu = exp(*linpred[1]);
     double epsi = (*response[2]) - (*linpred[2]);
     double hilfs = pow(sigu, 2) + pow(sigv, 2);


     double l;

       l =  -0.5*log(2.0*PI) + log(2.0) - 0.5*log(hilfs) - 0.5*pow(epsi,2)/hilfs + log(randnumbers::Phi2(-epsi*sigu/(pow(hilfs,0.5)*sigv)));


    *deviance = -2*l;
    }

  }


double DISTR_sfa0_mu_y::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa0_mu_y::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2]);
  }

 double DISTR_sfa0_mu_y::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_sfa0_mu_y::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {

    return ( 0 );
    }


double DISTR_sfa0_mu_y::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma_v equation
  // *worktransformlin[0] = sigma_v;
  // *worklin[1] = linear predictor of sigma_u equation
  // *worktransformlin[1] = sigma_u;

  if (counter==0)
    {
    set_worklin();
    }

  double epsi = (*response) - (*linpred);
  double sig2 = pow((*worktransformlin[0]),2) + pow((*worktransformlin[1]),2);
  double proparg2 = randnumbers::Phi2(-(epsi*(*worktransformlin[1]))/(pow(sig2,0.5)*(*worktransformlin[0])));
  double l;

     l =  -0.5*pow(epsi,2)/sig2 + log(proparg2);

  modify_worklin();

  return l;

  }


void DISTR_sfa0_mu_y::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma_v equation
  // *worktransformlin[0] = sigma_v;
  // *worklin[1] = linear predictor of sigma_u equation
  // *worktransformlin[1] = sigma_u;


  if (counter==0)
    {
    set_worklin();
    }

    double epsi = (*response) - (*linpred);

    double sig2 = pow((*worktransformlin[0]),2) + pow((*worktransformlin[1]),2);
    double arg2 = (-(epsi*(*worktransformlin[1]))/(pow(sig2,0.5)*(*worktransformlin[0])));
    double densarg2 = randnumbers::phi(arg2);
    double proparg2 = randnumbers::Phi2(arg2);
    double rat = densarg2/proparg2;
    double darg2 = (*worktransformlin[1])/(pow(sig2,0.5)*(*worktransformlin[0]));

    double nu = epsi/sig2 + rat*darg2;

    *workingweight = 1/sig2 - rat*(-arg2)*pow(darg2,2) + pow(rat*darg2, 2);

    if((*workingweight)<=0)
    {
        *workingweight = 0.0001;
    }

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  -0.5*pow(epsi,2)/sig2 + log(proparg2);

      }


  modify_worklin();

  }


void DISTR_sfa0_mu_y::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
    double muu = 0;
  double sigu2 = pow(exp((*linpred[predstart_mumult+1])), 2);
  double sigv2 = pow(exp((*linpred[predstart_mumult])), 2);
  double sig2 = sigu2 + sigv2;
  double ga = sigu2/sig2;
  double sigast = pow(ga * (1 - ga)*sig2, 0.5);
  double epsi = (*response[2]) - (*linpred[predstart_mumult+2]);
  double muast = (1 - ga) * muu - ga * epsi;

  *mu = randnumbers::Phi2(muast / sigast - sigast) * exp(-muast + 0.5 * pow(sigast, 2)) / randnumbers::Phi2(muast / sigast);

  }



void DISTR_sfa0_mu_y::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu_y): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa0_mu_y::update_end(void)
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

//--------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa_alpha -------------------------------
//--------------------------------------------------------------------------------


DISTR_sfa_alpha::DISTR_sfa_alpha(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Stochastic Frontier Analysis - alpha";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "alpha";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_sfa_alpha::DISTR_sfa_alpha(const DISTR_sfa_alpha & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa_alpha & DISTR_sfa_alpha::operator=(
                            const DISTR_sfa_alpha & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa_alpha::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa_alpha::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[0]));
  }


double DISTR_sfa_alpha::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma_v equation
  // *worklin[1] = linear predictor of sigma_u* equation
  // *worklin[2] = linear predictor of mu_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = sigma_v;
  // *worktransformlin[1] = sigma_u*;
  // *worktransformlin[2] = mu_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }
    double alpha = exp((*linpred));
    double muu = (*worktransformlin[2])*alpha;
    double sigu = (*worktransformlin[1])*alpha;
    double hilfs = pow(sigu, 2) + pow((*worktransformlin[0]), 2);
    double epsi = (*response) - (*worktransformlin[3]);

    double l;

    l =  - 0.5*log(hilfs) - pow(( epsi + muu ),2)/(2*hilfs) +
            log(randnumbers::Phi2( ((sigu*(*worktransformlin[0]))/pow(hilfs, 0.5))*(  - epsi/pow((*worktransformlin[0]),2) + muu/pow(sigu, 2) )) );


  modify_worklin();

  return l;

  }

void DISTR_sfa_alpha::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma_v equation
  // *worklin[1] = linear predictor of sigma_u* equation
  // *worklin[2] = linear predictor of mu_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = sigma_v;
  // *worktransformlin[1] = sigma_u*;
  // *worktransformlin[2] = mu_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }
    double alp = exp(*linpred);
    double muu = (*worktransformlin[2])*alp;
    double sigu = alp*(*worktransformlin[1]);
    double sigv = (*worktransformlin[0]);
    double epsi = (*response) - (*worktransformlin[3]);
    double hilfs = pow(sigu,2) + pow(sigv,2);
    double hilfs2 = sigu*sigv;
  //  double darg2 = -pow(sigu,2)/pow(hilfs, 1.5)*( - (epsi*sigu)/(sigv) + muu*sigv/sigu ) +
  //                  (1/pow(hilfs,0.5))*( - (epsi*sigu)/(sigv)  );
    double arg2 =  (hilfs2/pow(hilfs,0.5))*(  - (epsi)/pow(sigv,2) + muu/pow(sigu,2) );
    double darg2 = - (*worktransformlin[1]) * sigv * (pow(alp, 2) * (*worktransformlin[2]) + alp * epsi) / pow(hilfs, 1.5);
    double ddarg2 = 3 * pow(alp, 2) * pow((*worktransformlin[1]), 3) * sigv * (pow(alp, 2) * (*worktransformlin[2]) + epsi) / pow(hilfs, 2.5)
                    - 2 * pow(alp, 2) * (*worktransformlin[1]) * (*worktransformlin[2]) * sigv / pow(hilfs, 1.5);
  //  double ddarg2 = ((3*pow(sigu,2))/pow(hilfs,2.5))*( - (epsi*pow(sigu,3))/sigv + muu*sigu*sigv ) -
   //                 (1/pow(hilfs, 1.5))*( - (4*epsi*pow(sigu,3))/sigv ) - 2*muu*sigu*sigv/pow(hilfs, 1.5) - (epsi*sigu/sigv)/pow(hilfs, 0.5);

    double densarg2 = randnumbers::phi(arg2);

    double proparg2 = randnumbers::Phi2(arg2);

    double nu = - ( pow(sigu,2) + (epsi+muu)*muu )/hilfs + pow(((epsi+muu)*sigu/hilfs), 2) + (densarg2*darg2)/(randnumbers::Phi2(arg2));


    *workingweight = (2*pow(sigu,2)*pow(sigv,2))/pow(hilfs,2) + (2*(pow(epsi+muu,2))*(pow(sigu,4) - pow(sigu,2)*pow(sigv,2)))/pow(hilfs,3)  -
                        (densarg2*(pow(darg2,2)*(-arg2)+ddarg2))/proparg2 + pow((densarg2*darg2)/proparg2,2) + ( muu*(epsi+muu) + pow(muu, 2) )/hilfs
                        -( 4*(epsi+muu)*muu*pow(sigu, 2) )/pow(hilfs, 2);
 /*   *workingweight = (2*pow(sigu,2)*pow(sigv,2))/pow(hilfs,2) + (2*(pow(epsi+muu,2))*(pow(sigu,4) - pow(sigu,2)*pow(sigv,2)))/pow(hilfs,3)  -
                        (densarg2*(pow(darg2,2)*(-arg2)+ddarg2))/proparg2 + pow((densarg2*darg2)/proparg2,2) + (  pow(muu, 2) )/hilfs ;
                        -( 2*(epsi+muu)*muu*pow(sigu, 2) )/pow(hilfs, 2) + 2 * muu * (epsi + muu) * pow(sigv, 2) / pow(hilfs, 2); */

    if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - 0.5*log(hilfs) - pow(epsi+muu,2)/(2*hilfs) + log(proparg2);


      }

  modify_worklin();

  }


void DISTR_sfa_alpha::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (alpha): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa_alpha::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//--------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa_sigma_v ---------------------------
//--------------------------------------------------------------------------------


DISTR_sfa_sigma_v::DISTR_sfa_sigma_v(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Stochastic Frontier Analysis - sigma_v";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma_v";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_sfa_sigma_v::DISTR_sfa_sigma_v(const DISTR_sfa_sigma_v & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa_sigma_v & DISTR_sfa_sigma_v::operator=(
                            const DISTR_sfa_sigma_v & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa_sigma_v::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa_sigma_v::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[1]));
  }

double DISTR_sfa_sigma_v::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_u* equation
  // *worklin[2] = linear predictor of mu_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_u*;
  // *worktransformlin[2] = mu_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = (*worktransformlin[2])*(*worktransformlin[0]);
    double sigu = (*worktransformlin[1])*(*worktransformlin[0]);
    double sigu2 = pow(sigu,2);
    double sigv = exp((*linpred));
    double hilfs = sigu2 + pow(sigv, 2);
	double epsi = (*response) - (*worktransformlin[3]);

    double l;

    l = - 0.5*log(hilfs) - pow(( epsi + muu),2)/(2*hilfs) +
            log(randnumbers::Phi2( ((sigu*sigv)/pow(hilfs, 0.5))*(  - (epsi)/pow(sigv,2) + muu/sigu2 )) );



  modify_worklin();

  return l;

  }

void DISTR_sfa_sigma_v::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of mu_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = mu_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = (*worktransformlin[2])*(*worktransformlin[0]);
    double sigu = (*worktransformlin[1])*(*worktransformlin[0]);
    double sigu2 = pow(sigu,2);
    double sigv = exp((*linpred));
    double epsi = (*response) - (*worktransformlin[3]);
    double hilfs = sigu2 + pow(sigv,2);
    double hilfs2 = sigu*sigv;
    double darg2 = (-1/pow(hilfs, 1.5))*( - epsi*sigv*sigu + muu*pow(sigv,3)/sigu ) +
                    (1/pow(hilfs,0.5))*( (epsi*sigu)/(sigv) + muu*sigv/sigu );
    double arg2 =  (sigu*sigv/pow(hilfs,0.5))*(  - (epsi)/pow(sigv,2) + muu/sigu2 );
    double ddarg2 = ((3*pow(sigv,2))/pow(hilfs,2.5))*( muu*pow(sigv,3)/sigu - (epsi*hilfs2))  + arg2 -
                    4*muu*pow(sigv,3)/(sigu*pow(hilfs, 1.5));

    double densarg2 = randnumbers::phi(arg2);

    double proparg2 = randnumbers::Phi2(arg2);

    double nu = - pow(sigv,2)/hilfs + (pow((epsi+muu),2)*pow(sigv,2))/(pow(hilfs,2)) +
                    (densarg2*darg2)/(proparg2);


    *workingweight = (2*sigu2*pow(sigv,2))/pow(hilfs,2) + (2*(pow(epsi+muu,2))*(pow(sigv,4) - sigu2*pow(sigv,2)))/pow(hilfs,3)  -
                        (densarg2*(pow(darg2,2)*(-arg2)+ddarg2))/proparg2 + pow((densarg2*darg2)/proparg2,2);

    if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - 0.5*log(hilfs) - pow(epsi+muu,2)/(2*hilfs) + log(proparg2);

      }
  modify_worklin();

  }


void DISTR_sfa_sigma_v::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma_v): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa_sigma_v::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }



//--------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa_sigma_u -----------------------------
//--------------------------------------------------------------------------------


DISTR_sfa_sigma_u::DISTR_sfa_sigma_u(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Stochastic Frontier Analysis - sigma_u";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "sigma_u";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_sfa_sigma_u::DISTR_sfa_sigma_u(const DISTR_sfa_sigma_u & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa_sigma_u & DISTR_sfa_sigma_u::operator=(
                            const DISTR_sfa_sigma_u & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa_sigma_u::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));

  }

void DISTR_sfa_sigma_u::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[2]));
  }


double DISTR_sfa_sigma_u::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of mu_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = mu_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = (*worktransformlin[2])*(*worktransformlin[0]);
    double sigu = exp((*linpred))*(*worktransformlin[0]);
    double hilfs = pow(sigu, 2) + pow((*worktransformlin[1]), 2);
	double epsi = (*response) - (*worktransformlin[3]);

    double l;

    l =  - 0.5*log(hilfs) - pow((epsi + muu),2)/(2*hilfs) - log(randnumbers::Phi2(muu/sigu)) +
            log(randnumbers::Phi2( ((sigu*(*worktransformlin[1]))/pow(hilfs, 0.5))*(  - (epsi)/pow((*worktransformlin[1]),2) + muu/pow(sigu,2) )) );


  modify_worklin();

  return l;

  }

void DISTR_sfa_sigma_u::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of mu_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = mu_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = (*worktransformlin[2])*(*worktransformlin[0]);
    double sigu = exp((*linpred))*(*worktransformlin[0]);
    double epsi = (*response) - (*worktransformlin[3]);
    double hilfs = pow(sigu,2) + pow((*worktransformlin[1]),2);
	double hilfs2 = - (epsi*sigu)/(*worktransformlin[1]) + ((*worktransformlin[1])*muu)/sigu;
    double arg2 =  hilfs2/pow(hilfs,0.5);
    double darg2 = -(pow(sigu,2)/pow(hilfs, 1.5))*( hilfs2 )
					+ (- (epsi*sigu)/(*worktransformlin[1]) - (*worktransformlin[1])*muu/sigu )/pow(hilfs,0.5);

	double ddarg2 = ((3*pow(sigu,3))/pow(hilfs,2.5))*(muu*(*worktransformlin[1]) - (epsi*pow(sigu,2))/(*worktransformlin[1])) -
                    (1/pow(hilfs, 1.5))*( - (4*epsi*pow(sigu,3))/(*worktransformlin[1])) + arg2;

	double arg1 = muu/sigu;
	double darg1 = -arg1;

    double densarg1 = randnumbers::phi(arg1);
	double proparg1 = randnumbers::Phi2(arg1);
	double densarg2 = randnumbers::phi(arg2);
	double proparg2 = randnumbers::Phi2(arg2);

    double nu = - pow(sigu,2)/hilfs + (pow((epsi+muu),2)*pow(sigu,2))/(pow(hilfs,2)) +
                    (densarg2*darg2)/proparg2 - (densarg1*darg1)/proparg1;


    *workingweight = (2*pow(sigu,2)*pow((*worktransformlin[1]),2))/pow(hilfs,2) + (2*(pow(epsi+muu,2))*(pow(sigu,4) - pow(sigu,2)*pow((*worktransformlin[1]),2)))/pow(hilfs,3)  -
                        (densarg2*(pow(darg2,2)*(-arg2)+ddarg2))/proparg2 + pow((densarg2*darg2)/proparg2,2) +
						densarg1*(-arg1)*pow(darg1,2)/proparg1 + densarg1*arg1/proparg1 - pow(densarg1*darg1/proparg1,2);

    if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - 0.5*log(hilfs) - pow(epsi+muu,2)/(2*hilfs) - log(proparg1) + log(proparg2 );


      }

  modify_worklin();

  }

void DISTR_sfa_sigma_u::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double muu = exp((*linpred[predstart_mumult])) * exp(*linpred[predstart_mumult+3]);
  double sigu2 = pow(exp((*linpred[predstart_mumult+2])) * exp((*linpred[predstart_mumult])), 2);
  double sigv2 = pow(exp((*linpred[predstart_mumult+1])), 2);
  double sig2 = sigu2 + sigv2;
  double ga = sigu2/sig2;
  double sigast = pow(ga * (1 - ga)*sig2, 0.5);
  double epsi = (*response[4]) - (*linpred[predstart_mumult+4]);
  double muast = (1 - ga) * muu - ga * epsi;

  *mu = muast + sigast * randnumbers::phi(muast / sigast) / randnumbers::Phi2(muast / sigast);
  }

void DISTR_sfa_sigma_u::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma_u): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa_sigma_u::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//--------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa_mu_u --------------------------------
//--------------------------------------------------------------------------------


DISTR_sfa_mu_u::DISTR_sfa_mu_u(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Stochastic Frontier Analysis - mu_u";

  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu_u";
   linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_sfa_mu_u::DISTR_sfa_mu_u(const DISTR_sfa_mu_u & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa_mu_u & DISTR_sfa_mu_u::operator=(
                            const DISTR_sfa_mu_u & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa_mu_u::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa_mu_u::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[3]));
  }

double DISTR_sfa_mu_u::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }
  //  double muu = ((*linpred));
    double muu = exp((*linpred))*(*worktransformlin[0]);
	double epsi = (*response) -(*worktransformlin[3]);
	double sigu = (*worktransformlin[2])*(*worktransformlin[0]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[1]), 2);
    double hilfs = sigu2 + sigv2;

    double l;

    l = -pow((epsi + muu),2)/(2*hilfs) - log(randnumbers::Phi2(muu/sigu)) +
            log(randnumbers::Phi2( (((*worktransformlin[1])*sigu)/pow(hilfs, 0.5))*( (muu)/(sigu2) - (epsi)/sigv2 )) );



  modify_worklin();

  return l;

  }

void DISTR_sfa_mu_u::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }
  //  double muu = (*linpred);
    double muu = exp((*linpred))*(*worktransformlin[0]);
	double epsi = (*response) -(*worktransformlin[3]);
	double sigu = (*worktransformlin[2])*(*worktransformlin[0]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[1]), 2);
    double hilfs = sigu2 + sigv2;

	double arg2 =  (sigu*(*worktransformlin[1])/pow(hilfs,0.5))*( ( muu)/sigu2 - (epsi)/sigv2);

    double darg2 = (muu*(*worktransformlin[1]))/(sigu*pow(hilfs,0.5));
  //  double darg2 = ((*worktransformlin[1]))/((*worktransformlin[2])*pow(hilfs,0.5));
    double arg1 = muu/sigu;
    // double darg1 = arg1;
	double densarg1 = randnumbers::phi(arg1);
	double proparg1 = randnumbers::Phi2(arg1);
	double densarg2 = randnumbers::phi(arg2);
	double proparg2 = randnumbers::Phi2(arg2);

    double nu = - (muu*(epsi + muu))/(hilfs) - (densarg1*arg1)/proparg1 + (densarg2*darg2)/proparg2;
  //  double nu = - ((*worktransformlin[0])*(epsi + muu))/(hilfs) - (densarg1*darg1)/proparg1 + (densarg2*darg2)/proparg2;

 //   *workingweight =  pow((*worktransformlin[0]), 2)/hilfs + (densarg1*((-arg1)*pow(darg1,2)))/proparg1 -
 //                       pow((densarg1*darg1)/proparg1,2)
 //                   - (densarg2*((-arg2)*pow(darg2,2)))/proparg2 + pow((densarg2*darg2)/proparg2,2);
    *workingweight = -nu + pow(muu,2)/hilfs + (densarg1*((-arg1)*pow(arg1,2)))/proparg1 -
                        pow((densarg1*arg1)/proparg1,2)
                    - (densarg2*((-arg2)*pow(darg2,2)))/proparg2 + pow((densarg2*darg2)/proparg2,2);
	if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((epsi + muu),2)/(2*hilfs) - log(proparg1) + log(proparg2);


      }

  modify_worklin();

  }

  void DISTR_sfa_mu_u::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {


  double muu = exp((*linpred[predstart_mumult])) * exp(*linpred[predstart_mumult+3]);
  double sigu2 = pow(exp((*linpred[predstart_mumult+2])) * exp((*linpred[predstart_mumult])), 2);
  double sigv2 = pow(exp((*linpred[predstart_mumult+1])), 2);
  double sig2 = sigu2 + sigv2;
  double ga = sigu2/sig2;
  double sigast = pow(ga * (1 - ga)*sig2, 0.5);
  double epsi = (*response[4]) - (*linpred[predstart_mumult+4]);
  double muast = (1 - ga) * muu - ga * epsi;

  *mu = randnumbers::Phi2(muast / sigast - sigast) * exp(-muast + 0.5 * pow(sigast, 2)) / randnumbers::Phi2(muast / sigast);
  }



void DISTR_sfa_mu_u::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu_u): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa_mu_u::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_sfa_mu_y ----------------------------
//------------------------------------------------------------------------------
void DISTR_sfa_mu_y::check_errors(void)
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


DISTR_sfa_mu_y::DISTR_sfa_mu_y(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Stochastic Frontier Analysis -  mu_y";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu_y";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  check_errors();
  }


DISTR_sfa_mu_y::DISTR_sfa_mu_y(const DISTR_sfa_mu_y & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa_mu_y & DISTR_sfa_mu_y::operator=(
                            const DISTR_sfa_mu_y & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_sfa_mu_y::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response[2] = *response[3] = *response[4]
   // *linpred[0] = eta_alpha
   // *linpred[1] = eta_sigma_v
   // *linpred[2] = eta_sigma_u*
   // *linpred[3] = eta_mu_u*
   // *linpred[4] = eta_mu_y

   if (*weight[4] == 0)
     *deviance=0;
   else
     {
     double alpha = exp(*linpred[0]);
     double sigv = exp(*linpred[1]);
     double sigu = exp(*linpred[2])*alpha;
     double muu = exp(*linpred[3])*alpha;
     double epsi = (*response[4]) - (*linpred[4]);
     double sigu2 = pow(sigu, 2);
     double sigv2 = pow(sigv, 2);
     double hilfs1 = sigu2 + sigv2;
     double l;

       l = -0.5*log(2*PI) - 0.5*log(hilfs1) - pow((epsi+muu),2)/(2*hilfs1) -
            log(randnumbers::Phi2(muu/sigu)) + log(randnumbers::Phi2(((sigu*sigv)/pow(hilfs1,0.5))*( - epsi/sigv2 + muu/sigu2)));


    *deviance = -2*l;
    }

  }


double DISTR_sfa_mu_y::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa_mu_y::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[4]);
  }

 double DISTR_sfa_mu_y::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_sfa_mu_y::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
 //   double arg = ((*response[1])-(*param[1]))/((*param[0])) ;

 //   return (randnumbers::Phi2(arg));
      return 0;
    }


double DISTR_sfa_mu_y::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_u* equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_u*;

  if (counter==0)
    {
    set_worklin();
    }

  double epsi = (*response) - (*linpred);
  double sigu = (*worktransformlin[2])*(*worktransformlin[0]);
  double hilfs = pow(sigu, 2) + pow((*worktransformlin[1]), 2);
  double muu = (*worktransformlin[3])*(*worktransformlin[0]);

  double l;

     l = -pow((epsi + muu ),2)/(2*hilfs) +
            log(randnumbers::Phi2( ((sigu*(*worktransformlin[1]))/pow(hilfs, 0.5))*(  - epsi/pow((*worktransformlin[1]),2) + muu/pow(sigu, 2) )) );

  modify_worklin();

  return l;

  }


void DISTR_sfa_mu_y::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_u* equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_u*;


  if (counter==0)
    {
    set_worklin();
    }

    double epsi = (*response) - (*linpred);
    double sigu = (*worktransformlin[0])*(*worktransformlin[2]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[1]), 2);
    double hilfs = sigu2 + sigv2;
    double muu = (*worktransformlin[0])*(*worktransformlin[3]);
    double darg2 = sigu/((*worktransformlin[1])*pow(hilfs,0.5));
    double arg2 =  (((*worktransformlin[1])*sigu)/(pow(hilfs,0.5)))*(  - (epsi)/sigv2 + muu/sigu2 );

    double densarg2 = randnumbers::phi(arg2);

    double proparg2 = randnumbers::Phi2(arg2);

    double nu = (epsi+muu)/(hilfs) + (densarg2*darg2)/proparg2;

    *workingweight = 1/hilfs - (densarg2*(-arg2)*pow(darg2,2))/(proparg2) + pow(((densarg2*darg2)/proparg2),2);

    if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow(epsi+muu,2)/(2*hilfs) + log(proparg2);

      }


  modify_worklin();

  }


void DISTR_sfa_mu_y::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
 // double muu = exp((*linpred[predstart_mumult])) * exp(*linpred[predstart_mumult+3]);
  double sigu2 = pow(exp((*linpred[predstart_mumult+2])) * exp((*linpred[predstart_mumult])), 2);
  double sigv2 = pow(exp((*linpred[predstart_mumult+1])), 2);
  double sig2 = sigu2 + sigv2;

  *mu = sigu2/sig2;
//  double ga = sigu2/sig2;
 // double sigast = pow(ga * (1 - ga)*sig2, 0.5);
 // double epsi = (*response[4]) - (*linpred[predstart_mumult+4]);
 // double muast = (1 - ga) * muu - ga * epsi;

 // *mu = randnumbers::Phi2(muast / sigast - sigast) * exp(-muast + 0.5 * pow(sigast, 2)) / randnumbers::Phi2(muast / sigast);
  }


void DISTR_sfa_mu_y::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu_y): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa_mu_y::update_end(void)
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


//--------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa_mu_u_id -----------------------------
//--------------------------------------------------------------------------------


DISTR_sfa_mu_u_id::DISTR_sfa_mu_u_id(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Stochastic Frontier Analysis - mu_u";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "mu_u";
 //  linpredminlimit=-10;
 // linpredmaxlimit=15;

  }


DISTR_sfa_mu_u_id::DISTR_sfa_mu_u_id(const DISTR_sfa_mu_u_id & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa_mu_u_id & DISTR_sfa_mu_u_id::operator=(
                            const DISTR_sfa_mu_u_id & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa_mu_u_id::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa_mu_u_id::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = ((*linpred[3]));
  }

double DISTR_sfa_mu_u_id::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = ((*linpred))*(*worktransformlin[0]);
	double epsi = (*response) -(*worktransformlin[3]);
	double sigu = (*worktransformlin[2])*(*worktransformlin[0]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[1]), 2);
    double hilfs = sigu2 + sigv2;

    double l;

    l = -pow((epsi + muu),2)/(2*hilfs) - log(randnumbers::Phi2(muu/sigu)) +
            log(randnumbers::Phi2( (((*worktransformlin[1])*sigu)/pow(hilfs, 0.5))*( (muu)/(sigu2) - (epsi)/sigv2 )) );



  modify_worklin();

  return l;

  }

void DISTR_sfa_mu_u_id::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = ((*linpred))*(*worktransformlin[0]);
	double epsi = (*response) -(*worktransformlin[3]);
	double sigu = (*worktransformlin[2])*(*worktransformlin[0]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[1]), 2);
    double hilfs = sigu2 + sigv2;


	double arg2 =  (sigu*(*worktransformlin[1])/pow(hilfs,0.5))*( ( muu)/sigu2 - (epsi)/sigv2);

  //  double darg2 = (muu*(*worktransformlin[1]))/(sigu*pow(hilfs,0.5));
    double darg2 = ((*worktransformlin[1]))/((*worktransformlin[2])*pow(hilfs,0.5));
    double arg1 = muu/sigu;
    double darg1 = 1/(*worktransformlin[2]);
	double densarg1 = randnumbers::phi(arg1);
	double proparg1 = randnumbers::Phi2(arg1);
	double densarg2 = randnumbers::phi(arg2);
	double proparg2 = randnumbers::Phi2(arg2);

  //  double nu = - (muu*(epsi + muu))/(hilfs) - (densarg1*arg1)/proparg1 + (densarg2*darg2)/proparg2;
    double nu = - ((*worktransformlin[0])*(epsi + muu))/(hilfs) - (densarg1*darg1)/proparg1 + (densarg2*darg2)/proparg2;

    *workingweight =  pow((*worktransformlin[0]), 2)/hilfs + (densarg1*((-arg1)*pow(darg1,2)))/proparg1 -
                        pow((densarg1*darg1)/proparg1,2)
                    - (densarg2*((-arg2)*pow(darg2,2)))/proparg2 + pow((densarg2*darg2)/proparg2,2);
 /*   *workingweight = -nu + pow(muu,2)/hilfs + (densarg1*((-arg1)*pow(arg1,2)))/proparg1 -
                        pow((densarg1*arg1)/proparg1,2)
                    - (densarg2*((-arg2)*pow(darg2,2)))/proparg2 + pow((densarg2*darg2)/proparg2,2);*/
	if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((epsi + muu),2)/(2*hilfs) - log(proparg1) + log(proparg2);


      }

  modify_worklin();

  }


void DISTR_sfa_mu_u_id::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu_u): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa_mu_u_id::update_end(void)
  {

  // helpmat1 stores sigma2

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
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_sfa_mu_y_id -------------------------
//------------------------------------------------------------------------------
void DISTR_sfa_mu_y_id::check_errors(void)
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


DISTR_sfa_mu_y_id::DISTR_sfa_mu_y_id(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Stochastic Frontier Analysis -  mu_y";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu_y";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  check_errors();
  }


DISTR_sfa_mu_y_id::DISTR_sfa_mu_y_id(const DISTR_sfa_mu_y_id & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa_mu_y_id & DISTR_sfa_mu_y_id::operator=(
                            const DISTR_sfa_mu_y_id & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_sfa_mu_y_id::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response[2] = *response[3] = *response[4]
   // *linpred[0] = eta_alpha
   // *linpred[1] = eta_sigma_v
   // *linpred[2] = eta_sigma_u*
   // *linpred[3] = eta_mu_u*
   // *linpred[4] = eta_mu_y

   if (*weight[4] == 0)
     *deviance=0;
   else
     {
     double alpha = exp(*linpred[0]);
     double sigv = exp(*linpred[1]);
     double sigu = exp(*linpred[2])*alpha;
     double muu = (*linpred[3])*alpha;
     double epsi = (*response[4]) - (*linpred[4]);
     double sigu2 = pow(sigu, 2);
     double sigv2 = pow(sigv, 2);
     double hilfs1 = sigu2 + sigv2;
     double l;

       l = -0.5*log(2*PI) - 0.5*log(hilfs1) - pow((epsi+muu),2)/(2*hilfs1) -
            log(randnumbers::Phi2(muu/sigu)) + log(randnumbers::Phi2(((sigu*sigv)/pow(hilfs1,0.5))*( - epsi/sigv2 + muu/sigu2)));


    *deviance = -2*l;
    }

  }


double DISTR_sfa_mu_y_id::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa_mu_y_id::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[4]);
  }

 double DISTR_sfa_mu_y_id::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_sfa_mu_y_id::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
 //   double arg = ((*response[1])-(*param[1]))/((*param[0])) ;

 //   return (randnumbers::Phi2(arg));
      return 0;
    }


double DISTR_sfa_mu_y_id::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_u* equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_u*;

  if (counter==0)
    {
    set_worklin();
    }

  double epsi = (*response) - (*linpred);
  double sigu = (*worktransformlin[2])*(*worktransformlin[0]);
  double hilfs = pow(sigu, 2) + pow((*worktransformlin[1]), 2);
  double muu = (*worktransformlin[3])*(*worktransformlin[0]);

  double l;

     l = -pow((epsi + muu ),2)/(2*hilfs) +
            log(randnumbers::Phi2( ((sigu*(*worktransformlin[1]))/pow(hilfs, 0.5))*(  - epsi/pow((*worktransformlin[1]),2) + muu/pow(sigu, 2) )) );

  modify_worklin();

  return l;

  }


void DISTR_sfa_mu_y_id::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_u* equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_u*;


  if (counter==0)
    {
    set_worklin();
    }

    double epsi = (*response) - (*linpred);
    double sigu = (*worktransformlin[0])*(*worktransformlin[2]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[1]), 2);
    double hilfs = sigu2 + sigv2;
    double muu = (*worktransformlin[0])*(*worktransformlin[3]);
    double darg2 = sigu/((*worktransformlin[1])*pow(hilfs,0.5));
    double arg2 =  (((*worktransformlin[1])*sigu)/(pow(hilfs,0.5)))*(  - (epsi)/sigv2 + muu/sigu2 );

    double densarg2 = randnumbers::phi(arg2);

    double proparg2 = randnumbers::Phi2(arg2);

    double nu = (epsi+muu)/(hilfs) + (densarg2*darg2)/proparg2;

    *workingweight = 1/hilfs - (densarg2*(-arg2)*pow(darg2,2))/(proparg2) + pow(((densarg2*darg2)/proparg2),2);

    if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow(epsi+muu,2)/(2*hilfs) + log(proparg2);

      }


  modify_worklin();

  }


void DISTR_sfa_mu_y_id::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double muu = ((*linpred[predstart_mumult])) * exp(*linpred[predstart_mumult+3]);
  double sigu2 = pow(exp((*linpred[predstart_mumult+2])) * exp((*linpred[predstart_mumult])), 2);
  double sigv2 = pow(exp((*linpred[predstart_mumult+1])), 2);
  double sig2 = sigu2 + sigv2;
  double ga = sigu2/sig2;
  double sigast = pow(ga * (1 - ga)*sig2, 0.5);
  double epsi = (*response[4]) - (*linpred[predstart_mumult+4]);
  double muast = (1 - ga) * muu - ga * epsi;

  *mu = randnumbers::Phi2(muast / sigast - sigast) * exp(-muast + 0.5 * pow(sigast, 2)) / randnumbers::Phi2(muast / sigast);
  }


void DISTR_sfa_mu_y_id::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu_y): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa_mu_y_id::update_end(void)
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


//--------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa2_sigma_v ---------------------------
//--------------------------------------------------------------------------------


DISTR_sfa2_sigma_v::DISTR_sfa2_sigma_v(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,3,w)
  {
  family = "Stochastic Frontier Analysis - sigma_v";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma_v";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_sfa2_sigma_v::DISTR_sfa2_sigma_v(const DISTR_sfa2_sigma_v & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa2_sigma_v & DISTR_sfa2_sigma_v::operator=(
                            const DISTR_sfa2_sigma_v & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa2_sigma_v::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa2_sigma_v::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[0]));
  }


double DISTR_sfa2_sigma_v::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_u* equation
  // *worklin[2] = linear predictor of mu_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_u*;
  // *worktransformlin[2] = mu_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = (*worktransformlin[1]);
    double sigu = (*worktransformlin[0]);
    double sigu2 = pow(sigu,2);
    double sigv = exp((*linpred));
    double hilfs = sigu2 + pow(sigv, 2);
	double epsi = (*response) - (*worktransformlin[2]);

    double l;

    l = - 0.5*log(hilfs) - pow(( epsi + muu),2)/(2*hilfs) +
            log(randnumbers::Phi2( ((sigu*sigv)/pow(hilfs, 0.5))*(  - (epsi)/pow(sigv,2) + muu/sigu2 )) );



  modify_worklin();

  return l;

  }

void DISTR_sfa2_sigma_v::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of mu_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = mu_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = (*worktransformlin[1]);
    double sigu = (*worktransformlin[0]);
    double sigu2 = pow(sigu,2);
    double sigv = exp((*linpred));
    double epsi = (*response) - (*worktransformlin[2]);
    double hilfs = sigu2 + pow(sigv,2);
    double hilfs2 = sigu*sigv;
    double darg2 = (-1/pow(hilfs, 1.5))*( - epsi*sigv*sigu + muu*pow(sigv,3)/sigu ) +
                    (1/pow(hilfs,0.5))*( (epsi*sigu)/(sigv) + muu*sigv/sigu );
    double arg2 =  (sigu*sigv/pow(hilfs,0.5))*(  - (epsi)/pow(sigv,2) + muu/sigu2 );
    double ddarg2 = ((3*pow(sigv,2))/pow(hilfs,2.5))*( muu*pow(sigv,3)/sigu - (epsi*hilfs2))  + arg2 -
                    4*muu*pow(sigv,3)/(sigu*pow(hilfs, 1.5));

    double densarg2 = randnumbers::phi(arg2);

    double proparg2 = randnumbers::Phi2(arg2);

    double nu = - pow(sigv,2)/hilfs + (pow((epsi+muu),2)*pow(sigv,2))/(pow(hilfs,2)) +
                    (densarg2*darg2)/(proparg2);


    *workingweight = (2*sigu2*pow(sigv,2))/pow(hilfs,2) + (2*(pow(epsi+muu,2))*(pow(sigv,4) - sigu2*pow(sigv,2)))/pow(hilfs,3)  -
                        (densarg2*(pow(darg2,2)*(-arg2)+ddarg2))/proparg2 + pow((densarg2*darg2)/proparg2,2);

    if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - 0.5*log(hilfs) - pow(epsi+muu,2)/(2*hilfs) + log(proparg2);

      }
  modify_worklin();

  }


void DISTR_sfa2_sigma_v::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma_v): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa2_sigma_v::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }



//--------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa2_sigma_u -----------------------------
//--------------------------------------------------------------------------------


DISTR_sfa2_sigma_u::DISTR_sfa2_sigma_u(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,3,w)
  {
  family = "Stochastic Frontier Analysis - sigma_u";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma_u";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_sfa2_sigma_u::DISTR_sfa2_sigma_u(const DISTR_sfa2_sigma_u & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa2_sigma_u & DISTR_sfa2_sigma_u::operator=(
                            const DISTR_sfa2_sigma_u & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa2_sigma_u::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));

  }

void DISTR_sfa2_sigma_u::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[1]));
  }


double DISTR_sfa2_sigma_u::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of mu_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = mu_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = (*worktransformlin[1]);
    double sigu = exp((*linpred));
    double hilfs = pow(sigu, 2) + pow((*worktransformlin[0]), 2);
	double epsi = (*response) - (*worktransformlin[2]);

    double l;

    l =  - 0.5*log(hilfs) - pow((epsi + muu),2)/(2*hilfs) - log(randnumbers::Phi2(muu/sigu)) +
            log(randnumbers::Phi2( ((sigu*(*worktransformlin[0]))/pow(hilfs, 0.5))*(  - (epsi)/pow((*worktransformlin[0]),2) + muu/pow(sigu,2) )) );


  modify_worklin();

  return l;

  }

void DISTR_sfa2_sigma_u::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of mu_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = mu_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = (*worktransformlin[1]);
    double sigu = exp((*linpred));
    double epsi = (*response) - (*worktransformlin[2]);
    double hilfs = pow(sigu,2) + pow((*worktransformlin[0]),2);
	double hilfs2 = - (epsi*sigu)/(*worktransformlin[0]) + ((*worktransformlin[0])*muu)/sigu;
    double arg2 =  hilfs2/pow(hilfs,0.5);
    double darg2 = -(pow(sigu,2)/pow(hilfs, 1.5))*( hilfs2 )
					+ (- (epsi*sigu)/(*worktransformlin[0]) - (*worktransformlin[0])*muu/sigu )/pow(hilfs,0.5);

	double ddarg2 = ((3*pow(sigu,3))/pow(hilfs,2.5))*(muu*(*worktransformlin[0]) - (epsi*pow(sigu,2))/(*worktransformlin[0])) -
                    (1/pow(hilfs, 1.5))*( - (4*epsi*pow(sigu,3))/(*worktransformlin[0])) + arg2;

	double arg1 = muu/sigu;
	double darg1 = -arg1;

    double densarg1 = randnumbers::phi(arg1);
	double proparg1 = randnumbers::Phi2(arg1);
	double densarg2 = randnumbers::phi(arg2);
	double proparg2 = randnumbers::Phi2(arg2);

    double nu = - pow(sigu,2)/hilfs + (pow((epsi+muu),2)*pow(sigu,2))/(pow(hilfs,2)) +
                    (densarg2*darg2)/proparg2 - (densarg1*darg1)/proparg1;


    *workingweight = (2*pow(sigu,2)*pow((*worktransformlin[0]),2))/pow(hilfs,2) + (2*(pow(epsi+muu,2))*(pow(sigu,4) - pow(sigu,2)*pow((*worktransformlin[0]),2)))/pow(hilfs,3)  -
                        (densarg2*(pow(darg2,2)*(-arg2)+ddarg2))/proparg2 + pow((densarg2*darg2)/proparg2,2)
                         +densarg1*(-arg1)*pow(darg1,2)/proparg1 + densarg1*arg1/proparg1 - pow(densarg1*darg1/proparg1,2);

    if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  - 0.5*log(hilfs) - pow(epsi+muu,2)/(2*hilfs) - log(proparg1) + log(proparg2 );


      }

  modify_worklin();

  }


void DISTR_sfa2_sigma_u::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma_u): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa2_sigma_u::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }



//--------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa2_mu_u -------------------------------
//--------------------------------------------------------------------------------


DISTR_sfa2_mu_u::DISTR_sfa2_mu_u(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,3,w)
  {
  family = "Stochastic Frontier Analysis - mu_u";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "mu_u";
   linpredminlimit=-15;
  linpredmaxlimit=15;

  }


DISTR_sfa2_mu_u::DISTR_sfa2_mu_u(const DISTR_sfa2_mu_u & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa2_mu_u & DISTR_sfa2_mu_u::operator=(
                            const DISTR_sfa2_mu_u & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa2_mu_u::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa2_mu_u::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2]);
  }

double DISTR_sfa2_mu_u::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = exp(*linpred);
	double epsi = (*response) -(*worktransformlin[2]);
	double sigu = (*worktransformlin[1]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[0]), 2);
    double hilfs = sigu2 + sigv2;

    double l;

    l = -pow((epsi + muu),2)/(2*hilfs) - log(randnumbers::Phi2(muu/sigu)) +
            log(randnumbers::Phi2( (((*worktransformlin[0])*sigu)/pow(hilfs, 0.5))*( (muu)/(sigu2) - (epsi)/sigv2 )) );



  modify_worklin();

  return l;

  }

void DISTR_sfa2_mu_u::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }
    double muu = exp(*linpred);
	double epsi = (*response) -(*worktransformlin[2]);
	double sigu = (*worktransformlin[1]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[0]), 2);
    double hilfs = sigu2 + sigv2;

	double arg2 =  (sigu*(*worktransformlin[1])/pow(hilfs,0.5))*( ( muu)/sigu2 - (epsi)/sigv2);

    double darg2 = (muu*(*worktransformlin[1]))/(sigu*pow(hilfs,0.5));
  //  double darg2 = ((*worktransformlin[1]))/((*worktransformlin[2])*pow(hilfs,0.5));
    double arg1 = muu/sigu;
    // double darg1 = 1/(*worktransformlin[2]);
	double densarg1 = randnumbers::phi(arg1);
	double proparg1 = randnumbers::Phi2(arg1);
	double densarg2 = randnumbers::phi(arg2);
	double proparg2 = randnumbers::Phi2(arg2);

    double nu = - (muu*(epsi + muu))/(hilfs) - (densarg1*arg1)/proparg1 + (densarg2*darg2)/proparg2;
  //  double nu = - ((*worktransformlin[0])*(epsi + muu))/(hilfs) - (densarg1*darg1)/proparg1 + (densarg2*darg2)/proparg2;

 //   *workingweight =  pow((*worktransformlin[0]), 2)/hilfs + (densarg1*((-arg1)*pow(darg1,2)))/proparg1 -
 //                       pow((densarg1*darg1)/proparg1,2)
 //                   - (densarg2*((-arg2)*pow(darg2,2)))/proparg2 + pow((densarg2*darg2)/proparg2,2);
    *workingweight = -nu + pow(muu,2)/hilfs + (densarg1*((-arg1)*pow(arg1,2)))/proparg1 -
                        pow((densarg1*arg1)/proparg1,2)
                    - (densarg2*((-arg2)*pow(darg2,2)))/proparg2 + pow((densarg2*darg2)/proparg2,2);
	if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((epsi + muu),2)/(2*hilfs) - log(proparg1) + log(proparg2);


      }

  modify_worklin();

  }


void DISTR_sfa2_mu_u::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu_u): log\n");
  optionsp->out("\n");
  optionsp->out("\n") ;
  }


void DISTR_sfa2_mu_u::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_sfa2_mu_y ----------------------------
//------------------------------------------------------------------------------
void DISTR_sfa2_mu_y::check_errors(void)
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


DISTR_sfa2_mu_y::DISTR_sfa2_mu_y(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,3,w)
  {
  family = "Stochastic Frontier Analysis -  mu_y";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu_y";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  check_errors();
  }


DISTR_sfa2_mu_y::DISTR_sfa2_mu_y(const DISTR_sfa2_mu_y & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa2_mu_y & DISTR_sfa2_mu_y::operator=(
                            const DISTR_sfa2_mu_y & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_sfa2_mu_y::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response[2] = *response[3] = *response[4]
   // *linpred[0] = eta_alpha
   // *linpred[1] = eta_sigma_v
   // *linpred[2] = eta_sigma_u*
   // *linpred[3] = eta_mu_u*
   // *linpred[4] = eta_mu_y

   if (*weight[3] == 0)
     *deviance=0;
   else
     {
     double sigv = exp(*linpred[0]);
     double sigu = exp(*linpred[1]);
     double muu = exp(*linpred[2]);
     double epsi = (*response[3]) - (*linpred[3]);
     double sigu2 = pow(sigu, 2);
     double sigv2 = pow(sigv, 2);
     double hilfs1 = sigu2 + sigv2;
     double l;

       l = -0.5*log(2*PI) - 0.5*log(hilfs1) - pow((epsi+muu),2)/(2*hilfs1) -
            log(randnumbers::Phi2(muu/sigu)) + log(randnumbers::Phi2(((sigu*sigv)/pow(hilfs1,0.5))*( - epsi/sigv2 + muu/sigu2)));


    *deviance = -2*l;
    }

  }


double DISTR_sfa2_mu_y::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa2_mu_y::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[3]);
  }

 double DISTR_sfa2_mu_y::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_sfa2_mu_y::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
 //   double arg = ((*response[1])-(*param[1]))/((*param[0])) ;

 //   return (randnumbers::Phi2(arg));
      return 0;
    }


double DISTR_sfa2_mu_y::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_u* equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_u*;

  if (counter==0)
    {
    set_worklin();
    }

  double epsi = (*response) - (*linpred);
  double sigu = (*worktransformlin[1]);
  double hilfs = pow(sigu, 2) + pow((*worktransformlin[0]), 2);
  double muu = (*worktransformlin[2]);

  double l;

     l = -pow((epsi + muu ),2)/(2*hilfs) +
            log(randnumbers::Phi2( ((sigu*(*worktransformlin[0]))/pow(hilfs, 0.5))*(  - epsi/pow((*worktransformlin[0]),2) + muu/pow(sigu, 2) )) );

  modify_worklin();

  return l;

  }


void DISTR_sfa2_mu_y::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_u* equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_u*;


  if (counter==0)
    {
    set_worklin();
    }

    double epsi = (*response) - (*linpred);
    double sigu = (*worktransformlin[1]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[0]), 2);
    double hilfs = sigu2 + sigv2;
    double muu = (*worktransformlin[2]);
    double darg2 = sigu/((*worktransformlin[0])*pow(hilfs,0.5));
    double arg2 =  (((*worktransformlin[0])*sigu)/(pow(hilfs,0.5)))*(  - (epsi)/sigv2 + muu/sigu2 );

    double densarg2 = randnumbers::phi(arg2);

    double proparg2 = randnumbers::Phi2(arg2);

    double nu = (epsi+muu)/(hilfs) + (densarg2*darg2)/proparg2;

    *workingweight = 1/hilfs - (densarg2*(-arg2)*pow(darg2,2))/(proparg2) + pow(((densarg2*darg2)/proparg2),2);

    if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow(epsi+muu,2)/(2*hilfs) + log(proparg2);

      }


  modify_worklin();

  }


void DISTR_sfa2_mu_y::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double muu = exp(*linpred[predstart_mumult+2]);
  double sigu2 = pow(exp((*linpred[predstart_mumult+1])), 2);
  double sigv2 = pow(exp((*linpred[predstart_mumult+0])), 2);
  double sig2 = sigu2 + sigv2;
  double ga = sigu2/sig2;
  double sigast = pow(ga * (1 - ga)*sig2, 0.5);
  double epsi = (*response[3]) - (*linpred[predstart_mumult+3]);
  double muast = (1 - ga) * muu - ga * epsi;

  *mu = randnumbers::Phi2(muast / sigast - sigast) * exp(-muast + 0.5 * pow(sigast, 2)) / randnumbers::Phi2(muast / sigast);

  }


void DISTR_sfa2_mu_y::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu_y): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa2_mu_y::update_end(void)
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


//--------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_sfa2_mu_u_id --------------------------------
//--------------------------------------------------------------------------------


DISTR_sfa2_mu_u_id::DISTR_sfa2_mu_u_id(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,3,w)
  {
  family = "Stochastic Frontier Analysis - mu_u";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "mu_u";
   linpredminlimit=-15;
  linpredmaxlimit=15;

  }


DISTR_sfa2_mu_u_id::DISTR_sfa2_mu_u_id(const DISTR_sfa2_mu_u_id & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa2_mu_u_id & DISTR_sfa2_mu_u_id::operator=(
                            const DISTR_sfa2_mu_u_id & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sfa2_mu_u_id::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa2_mu_u_id::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2]);
  }

double DISTR_sfa2_mu_u_id::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }

    double muu = (*linpred);
	double epsi = (*response) -(*worktransformlin[2]);
	double sigu = (*worktransformlin[1]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[0]), 2);
    double hilfs = sigu2 + sigv2;

    double l;

    l = -pow((epsi + muu),2)/(2*hilfs) - log(randnumbers::Phi2(muu/sigu)) +
            log(randnumbers::Phi2( (((*worktransformlin[0])*sigu)/pow(hilfs, 0.5))*( (muu)/(sigu2) - (epsi)/sigv2 )) );



  modify_worklin();

  return l;

  }

void DISTR_sfa2_mu_u_id::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_y equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_y;

  if (counter==0)
    {
    set_worklin();
    }
    double muu = (*linpred);
	double epsi = (*response) -(*worktransformlin[2]);
	double sigu = (*worktransformlin[1]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[0]), 2);
    double hilfs = sigu2 + sigv2;

	double arg2 =  (sigu*(*worktransformlin[0])/pow(hilfs,0.5))*( ( muu)/sigu2 - (epsi)/sigv2);

    double darg2 = ((*worktransformlin[0]))/(sigu*pow(hilfs,0.5));
    double arg1 = muu/sigu;
    double darg1 = 1/sigu;
	double densarg1 = randnumbers::phi(arg1);
	double proparg1 = randnumbers::Phi2(arg1);
	double densarg2 = randnumbers::phi(arg2);
	double proparg2 = randnumbers::Phi2(arg2);

    double nu = - ((epsi + muu))/(hilfs) - (densarg1*darg1)/proparg1 + (densarg2*darg2)/proparg2;

    *workingweight =  1/hilfs - (densarg2*((-arg2)*pow(darg2,2)))/proparg2 + pow((densarg2*darg2)/proparg2,2);
                    //- pow((densarg1*darg1)/proparg1,2) + (densarg1*((-arg1)*pow(darg1,2)))/proparg1 ;

	if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow((epsi + muu),2)/(2*hilfs) - log(proparg1) + log(proparg2);


      }

  modify_worklin();

  }


void DISTR_sfa2_mu_u_id::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (mu_u): identity\n");
  optionsp->out("\n");
  optionsp->out("\n") ;
  }


void DISTR_sfa2_mu_u_id::update_end(void)
  {

  // helpmat1 stores sigma2

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
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_sfa2_mu_y_id ------------------------
//------------------------------------------------------------------------------
void DISTR_sfa2_mu_y_id::check_errors(void)
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


DISTR_sfa2_mu_y_id::DISTR_sfa2_mu_y_id(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,3,w)
  {
  family = "Stochastic Frontier Analysis -  mu_y";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu_y";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  check_errors();
  }


DISTR_sfa2_mu_y_id::DISTR_sfa2_mu_y_id(const DISTR_sfa2_mu_y_id & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_sfa2_mu_y_id & DISTR_sfa2_mu_y_id::operator=(
                            const DISTR_sfa2_mu_y_id & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_sfa2_mu_y_id::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response[2] = *response[3] = *response[4]
   // *linpred[0] = eta_alpha
   // *linpred[1] = eta_sigma_v
   // *linpred[2] = eta_sigma_u*
   // *linpred[3] = eta_mu_u*
   // *linpred[4] = eta_mu_y

   if (*weight[3] == 0)
     *deviance=0;
   else
     {
     double sigv = exp(*linpred[0]);
     double sigu = exp(*linpred[1]);
     double muu = (*linpred[2]);
     double epsi = (*response[3]) - (*linpred[3]);
     double sigu2 = pow(sigu, 2);
     double sigv2 = pow(sigv, 2);
     double hilfs1 = sigu2 + sigv2;
     double l;

       l = -0.5*log(2*PI) - 0.5*log(hilfs1) - pow((epsi+muu),2)/(2*hilfs1) -
            log(randnumbers::Phi2(muu/sigu)) + log(randnumbers::Phi2(((sigu*sigv)/pow(hilfs1,0.5))*( - epsi/sigv2 + muu/sigu2)));


    *deviance = -2*l;
    }

  }


double DISTR_sfa2_mu_y_id::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_sfa2_mu_y_id::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[3]);
  }

 double DISTR_sfa2_mu_y_id::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_sfa2_mu_y_id::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
 //   double arg = ((*response[1])-(*param[1]))/((*param[0])) ;

 //   return (randnumbers::Phi2(arg));
      return 0;
    }


double DISTR_sfa2_mu_y_id::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_u* equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_u*;

  if (counter==0)
    {
    set_worklin();
    }

  double epsi = (*response) - (*linpred);
  double sigu = (*worktransformlin[1]);
  double hilfs = pow(sigu, 2) + pow((*worktransformlin[0]), 2);
  double muu = (*worktransformlin[2]);

  double l;

     l = -pow((epsi + muu ),2)/(2*hilfs) +
            log(randnumbers::Phi2( ((sigu*(*worktransformlin[0]))/pow(hilfs, 0.5))*(  - epsi/pow((*worktransformlin[0]),2) + muu/pow(sigu, 2) )) );

  modify_worklin();

  return l;

  }


void DISTR_sfa2_mu_y_id::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of alpha equation
  // *worklin[1] = linear predictor of sigma_v equation
  // *worklin[2] = linear predictor of sigma_u* equation
  // *worklin[3] = linear predictor of mu_u* equation
  // *worktransformlin[0] = alpha;
  // *worktransformlin[1] = sigma_v;
  // *worktransformlin[2] = sigma_u*;
  // *worktransformlin[3] = mu_u*;


  if (counter==0)
    {
    set_worklin();
    }

    double epsi = (*response) - (*linpred);
    double sigu = (*worktransformlin[1]);
    double sigu2 = pow(sigu, 2);
    double sigv2 = pow((*worktransformlin[0]), 2);
    double hilfs = sigu2 + sigv2;
    double muu = (*worktransformlin[2]);
    double darg2 = sigu/((*worktransformlin[0])*pow(hilfs,0.5));
    double arg2 =  (((*worktransformlin[0])*sigu)/(pow(hilfs,0.5)))*(  - (epsi)/sigv2 + muu/sigu2 );

    double densarg2 = randnumbers::phi(arg2);

    double proparg2 = randnumbers::Phi2(arg2);

    double nu = (epsi+muu)/(hilfs) + (densarg2*darg2)/proparg2;

    *workingweight = 1/hilfs - (densarg2*(-arg2)*pow(darg2,2))/(proparg2) + pow(((densarg2*darg2)/proparg2),2);

    if ((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += -pow(epsi+muu,2)/(2*hilfs) + log(proparg2);

      }


  modify_worklin();

  }


void DISTR_sfa2_mu_y_id::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double muu = (*linpred[predstart_mumult+2]);
  double sigu2 = pow(exp((*linpred[predstart_mumult+1])), 2);
  double sigv2 = pow(exp((*linpred[predstart_mumult+0])), 2);
  double sig2 = sigu2 + sigv2;
  double ga = sigu2/sig2;
  double sigast = pow(ga * (1 - ga)*sig2, 0.5);
  double epsi = (*response[3]) - (*linpred[predstart_mumult+3]);
  double muast = (1 - ga) * muu - ga * epsi;

  *mu = randnumbers::Phi2(muast / sigast - sigast) * exp(-muast + 0.5 * pow(sigast, 2)) / randnumbers::Phi2(muast / sigast);

  }


void DISTR_sfa2_mu_y_id::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu_y): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sfa2_mu_y_id::update_end(void)
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


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_hurdle_pi ----------------------------
//------------------------------------------------------------------------------


DISTR_hurdle_pi::DISTR_hurdle_pi(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,0,w)
  {
  family = "hurdle - pi";
  outpredictor = true;
  outexpectation = false;
  predictor_name = "pi";
    linpredminlimit=-10;
  linpredmaxlimit=10;

  }


DISTR_hurdle_pi::DISTR_hurdle_pi(const DISTR_hurdle_pi & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_hurdle_pi & DISTR_hurdle_pi::operator=(
                            const DISTR_hurdle_pi & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_hurdle_pi::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_hurdle_pi::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[0]))/(1+exp((*linpred[0])));
  }

double DISTR_hurdle_pi::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = exp(eta_mu);

  if (counter==0)
    {
    set_worklin();
    }
    double explinpi;
   if (*linpred <= linpredminlimit)
    explinpi = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    exptildeeta = exp(linpredmaxlimit);
   else
    explinpi = exp(*linpred);
    double pi = explinpi/(1+explinpi);

    double l;
    if((*response) == 0)
    {
        l = log(pi);
    } else
    {
        l = log(1-pi);
    }


  modify_worklin();

  return l;

  }

void DISTR_hurdle_pi::compute_iwls_wweightschange_weightsone(
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
    }

  double explinpi;
  if (*linpred <= linpredminlimit)
    explinpi = exp(linpredminlimit);
//  else if (*linpred >= linpredmaxlimit)
//    exptildeeta = exp(linpredmaxlimit);
  else
    explinpi = exp(*linpred);
    double pi = explinpi/(1+explinpi);


    double nu = -pi;
    if((*response)==0) {
        nu += 1;
    }


    *workingweight = pi*(1-pi);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {
        if((*response)==0)
        {
            like += log(pi);
        }
        else{
            like += log(1-pi);
        }

      }

  modify_worklin();

  }


void DISTR_hurdle_pi::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (pi): logit\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_hurdle_pi::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin)/(1+exp(*worklin));
    }

  }


//------------------------------------------------------------------------------
//--------------------------- CLASS: DISTR_hurdle_lambda ----------------------
//------------------------------------------------------------------------------
void DISTR_hurdle_lambda::check_errors(void)
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


DISTR_hurdle_lambda::DISTR_hurdle_lambda(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {
  family = "hurdle poisson - lambda";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "lamba";
    linpredminlimit=-10;
  linpredmaxlimit=15;
  check_errors();
  }


DISTR_hurdle_lambda::DISTR_hurdle_lambda(const DISTR_hurdle_lambda & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }


const DISTR_hurdle_lambda & DISTR_hurdle_lambda::operator=(
                            const DISTR_hurdle_lambda & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_hurdle_lambda::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_pi
   // *linpred[1] = eta_lamba

   if (*weight[0] == 0)
     *deviance=0;
   else
     {
    double lambda;
    double expminuslambda;
    double l;
     double explinpi = exp(*linpred[0]);
     // double p =  explinpi / (1 + explinpi);
     if (*linpred[1] <= linpredminlimit)
       lambda = exp(linpredminlimit);
//     else if (*linpred[1] >= linpredmaxlimit)
//       lambda = exp(linpredmaxlimit);
     else
       lambda = exp(*linpred[1]);

     expminuslambda = exp(-lambda);

      if (*response[0]==0)
         {
         l= -log(1+ explinpi) + (*linpred[0]);
         }
       else // response > 0
         {
         double help1 = (*response[1])+1;
         l= -log(1+ explinpi) + (*response[1])*(*linpred[1])- lambda
            - randnumbers::lngamma_exact(help1)-log(1-exp(-lambda));
         }


    *deviance = -2*l;
    }

  }


double DISTR_hurdle_lambda::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_hurdle_lambda::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]);
  }
 double DISTR_hurdle_lambda::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_hurdle_lambda::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {

    return 0;
    }


double DISTR_hurdle_lambda::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

  double lambda = exp(*linpred);

  // double expminuslambda = exp(-lambda);

  double l = -log(1-exp(-lambda)) + (*response)*(*linpred)-lambda;


  modify_worklin();

  return l;

  }


void DISTR_hurdle_lambda::compute_iwls_wweightschange_weightsone(
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

    double w = (*linpred);
   // double w0 = exp(*linpred);

    double lambda = exp(w);
    double expminuslambda = exp(-lambda);

    double nu = (*response) - lambda  - lambda*expminuslambda/(1-expminuslambda);


  *workingweight = -lambda*(-1+expminuslambda*(1+lambda))/pow((1-expminuslambda), 2); //((1-(*worktransformlin[0]))*

  *workingresponse = *linpred + nu/(*workingweight);




    if (compute_like)
      {

            like += -log(1-expminuslambda) + (*response)*(*linpred)-lambda;


      }
  modify_worklin();

  }


void DISTR_hurdle_lambda::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
   double lambda = exp((*linpred[predstart_mumult+1]));
   double p = exp((*linpred[predstart_mumult])) / (1 + exp((*linpred[predstart_mumult])));
  *mu = (1-p)*lambda/(1-exp(-lambda));

  }


void DISTR_hurdle_lambda::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (lambda): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_hurdle_lambda::update_end(void)
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
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_hurdle_delta --------------------------
//------------------------------------------------------------------------------

DISTR_hurdle_delta::DISTR_hurdle_delta(GENERAL_OPTIONS * o,
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


  family = "hurdle negative binomial - delta";

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
  linpredmaxlimit=15;
  }


DISTR_hurdle_delta::DISTR_hurdle_delta(const DISTR_hurdle_delta & nd)
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


const DISTR_hurdle_delta & DISTR_hurdle_delta::operator=(
                            const DISTR_hurdle_delta & nd)
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



double DISTR_hurdle_delta::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_hurdle_delta::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]);
  }

double DISTR_hurdle_delta::loglikelihood_weightsone(double * response,
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

    double resp_plus_delta = (*response) + delta;

    l = randnumbers::lngamma_exact(resp_plus_delta) -
        randnumbers::lngamma_exact(delta) -
        (*response)* log_mu_plus_delta
        -log(pow((delta+(*worktransformlin[0]))/delta, delta)-1);


  modify_worklin();

  return l;

  }


void DISTR_hurdle_delta::compute_expectation(void)
  {

  int k=1;
  double k_delta;
  double kplus1;
  double psum;

  double L = 0;
  E_dig_y_delta = 0;
  E_trig_y_delta = 0;

  psum = L;

  while ((psum < stopsum) && (k <=stoprmax))
    {
    k_delta = k + delta;
    kplus1 = k + 1;

    L = exp(randnumbers::lngamma_exact(k_delta) -
            randnumbers::lngamma_exact(kplus1) -
            lngamma_delta -log(pow(delta_plus_mu/delta, delta)-1) +
            k* log((*worktransformlin[0])/delta_plus_mu)
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


void DISTR_hurdle_delta::compute_iwls_wweightschange_weightsone(
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

  lngamma_delta = randnumbers::lngamma_exact(delta);

  double delta_plus_response = delta + (*response);
  double hilfs = pow(delta/delta_plus_mu, delta);

  double nu = delta*(randnumbers::digamma_exact(delta_plus_response) -
                    randnumbers::digamma_exact(delta))
                    -delta*(*response)/(delta_plus_mu)-
                    delta*(log(delta_plus_mu/delta)-(*worktransformlin[0])/delta_plus_mu)/(1-hilfs);

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

  *workingweight = delta*(log(delta_plus_mu/delta)-(*worktransformlin[0])/delta_plus_mu)/(1-hilfs) - pow(((*worktransformlin[0])/delta_plus_mu), 2)*delta/(1-hilfs)
                   -pow((delta*hilfs*(log(delta/delta_plus_mu)+(*worktransformlin[0])/delta_plus_mu)/(1-hilfs)),2)
                   -E_dig_y_delta-E_trig_y_delta + delta*pow((*worktransformlin[0]),2)/(pow(delta+(*worktransformlin[0]),2)*(1-hilfs));

  if (*workingweight <= 0)
    *workingweight = 0.0001;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    double resp_plus_delta = (*response) + delta;
    double log_mu_plus_delta = log((*worktransformlin[0]) + delta);
    like += randnumbers::lngamma_exact(resp_plus_delta) -
        randnumbers::lngamma_exact(delta) -
        (*response)* log_mu_plus_delta
        -log(pow((delta+(*worktransformlin[0]))/delta, delta)-1);

    }

  modify_worklin();
  Ep++;
  Ep_trig++;

  }


void DISTR_hurdle_delta::outoptions(void)
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


void DISTR_hurdle_delta::update_end(void)
  {
  DISTR_gamlss::update_end();
  }


//------------------------------------------------------------------------------
//------------------------- CLASS: DISTR_hurdle_mu -----------------------------
//------------------------------------------------------------------------------

void DISTR_hurdle_mu::check_errors(void)
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


DISTR_hurdle_mu::DISTR_hurdle_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,1,w)
  {

  predictor_name = "mu";
  outpredictor = true;
  outexpectation = true;

  family = "Hurdle negative binomial - mu";

  linpredminlimit=-10;
  linpredmaxlimit=15;

  check_errors();
  }


DISTR_hurdle_mu::DISTR_hurdle_mu(const DISTR_hurdle_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  }


const DISTR_hurdle_mu & DISTR_hurdle_mu::operator=(
                            const DISTR_hurdle_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


void DISTR_hurdle_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[0] = *response[1] = response
   // *linpred[0] = eta_delta
   // *linpred[1] = eta_pi
   // *linpred[2] = eta_mu

   if (*weight[0] == 0)
     *deviance=0;
   else
     {
     double explinpi = exp(*linpred[0]);
     double delta = exp(*linpred[1]);
     double mu = exp(*linpred[2]);
     double delta_plus_mu = delta + mu;
     double resp_plus_one = (*response[0]) + 1;

     double l;

     if (*response[0]==0)
         {
         l= -log(1+ explinpi) + (*linpred[0]);
         }
       else // response > 0
         {
       double delta_plus_response = delta+(*response[0]);

       l = randnumbers::lngamma_exact(delta_plus_response) -
           randnumbers::lngamma_exact(resp_plus_one) -
           randnumbers::lngamma_exact(delta) +
           (*response[0])*log(mu/delta_plus_mu)-log(1+ explinpi) - log(pow(delta_plus_mu/delta, delta)-1);

       }
    *deviance = -2*l;
    }

  }


double DISTR_hurdle_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_hurdle_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[2]);
  }


double DISTR_hurdle_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
 //   double p =  (*param[1])/((*param[0])+(*param[1]));
 //   double r = (*param[0]);
 //   double u = gsl_ran_negative_binomial_pdf(*response[1], p, r);
 //   return u;
      return 0;
    }

double DISTR_hurdle_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
  //  double p =  (*param[1])/((*param[0])+(*param[1]));
  //  double r = (*param[0]);
  //  double u = gsl_cdf_negative_binomial_P(*response[1], p, r);
  //  return u;
      return 0.0;
    }


double DISTR_hurdle_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of delta equation
  // *worktransformlin[0] = exp(eta_delta);
  // *worklin[1] = linear predictor of pi equation
  // *worktransformlin[1] = pi;

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
     l = -(*response)*log((*worktransformlin[0])+mu) +(*response)*log(mu)
           - log(pow((mu + (*worktransformlin[0])) / (*worktransformlin[0]), (*worktransformlin[0])) - 1);

  modify_worklin();

  return l;

  }


void DISTR_hurdle_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

   // *worklin[0] = linear predictor of delta equation
  // *worktransformlin[0] = exp(eta_delta);
  // *worklin[1] = linear predictor of pi equation
  // *worktransformlin[1] = pi;

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
    double hilfs = pow((*worktransformlin[0])/delta_plus_mu,(*worktransformlin[0]));

    double nu = (*worktransformlin[0])*((*response))/delta_plus_mu
                    - mu*(*worktransformlin[0])/ (delta_plus_mu*(1-hilfs));

  *workingweight = (*worktransformlin[0])*mu*(mu+(*worktransformlin[0]))/(pow(delta_plus_mu, 2)*(1-hilfs))+
                    -pow(mu*(*worktransformlin[0]), 2)*hilfs/(pow(delta_plus_mu*(1-hilfs),2));

  if (*workingweight <= 0)
    *workingweight = 0.0001;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    like += -(*response)*log((*worktransformlin[0])+mu) +(*response)*log(mu)
           - log(pow((mu + (*worktransformlin[0])) / (*worktransformlin[0]), (*worktransformlin[0])) - 1);

    }

  modify_worklin();

  }


void DISTR_hurdle_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double delta = exp(*linpred[predstart_mumult+1]);
  double p = exp(*linpred[predstart_mumult])/(1+exp(*linpred[predstart_mumult]));
  double m = exp(*linpred[predstart_mumult+2]);
  *mu = (1-p)*m/(1-pow(delta/(m+delta), delta));
  }


void DISTR_hurdle_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_hurdle_mu::update_end(void)
  {
  DISTR_gamlss::update_end();
  }


//--------SNDP------------------------------

// bereits ersetzt
//  double ll_sn_dp(double x, double xi, double om, double al){
//    double l = 0;
//    l = -0.5*log(2*PI) - log(om) - 0.5*pow(x-xi, 2)/om/om+log(randnumbers::Phi2(al*(x-xi)/om));
//    return l;
//  }

 double sn_a0(double al){
    double a = al * 0.96;
    return 2/PI/sqrt(1 + 8*a*a/PI/PI);
  }

  double sn_a1(double al){
    double a = al;
    double b = 2/PI;
    return - 0.975 * pow(sqrt(b),3)*a/pow(sqrt(1 + 8*a*a / (PI*PI)),3) / sqrt(1 + a*a/(1 + 8*a*a/(PI*PI)));
  }

  double sn_a2(double al){
    double a = al * 0.96;
    return 2/PI/pow(sqrt(1 + 8*a*a/PI/PI), 3);
  }


// Constructor:
DISTR_sndp_alpha::DISTR_sndp_alpha(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w) // <-----
  {
  family = "Skew Normal Distribution, Direct Parametrization - alpha";
  outpredictor = true;
  outexpectation = false; // only true if main parameter
  predictor_name = "alpha";

  // linpredminlimit=-10; // necessary for log-link
  // linpredmaxlimit=15; // necessary for log-link

  // updateIWLS = false; // only false and necessary if gibbs sampling?
  }


// Copy Constructor:
DISTR_sndp_alpha::DISTR_sndp_alpha(const DISTR_sndp_alpha & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }

// Operator:
const DISTR_sndp_alpha & DISTR_sndp_alpha::operator=(
                            const DISTR_sndp_alpha & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }



double DISTR_sndp_alpha::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_sndp_alpha::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[0]); // <--- c.f. linpred info
  }



double DISTR_sndp_alpha::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }

  double x = (*response);
  double xi = *worklin[1];
  double om = exp(*worklin[0]);
  double al = (*linpred);

  double l;

  l = - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // calculate log-likelihood (possibly simplified)

  modify_worklin();

  return l;

  }


void DISTR_sndp_alpha::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  if (counter==0)
    {
    set_worklin();
    }

    double x = (*response);
    double xi = (*worklin[1]);
    double om = exp(*worklin[0]);
    double al = (*linpred);
    double z = (x - xi)/om;

    double nu = randnumbers::phi(al*z) / randnumbers::Phi2(al*z) * z; // v

//    *workingweight = zeta2(al*z)*z*z; //
    *workingweight = sn_a2(al); // w

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // loglikelihood

      }


  modify_worklin();

  }


void DISTR_sndp_alpha::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("Link function (alpha): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sndp_alpha::update_end(void)
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


// Constructor:
DISTR_sndp_omega::DISTR_sndp_omega(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w) // <-----
  {
  family = "Skew Normal Distribution, Direct Parametrization - omega";
  outpredictor = true;
  outexpectation = false; // only true if main parameter
  predictor_name = "omega";

   linpredminlimit=-10; // necessary for log-link
   linpredmaxlimit=15; // necessary for log-link

  // updateIWLS = false; // only false and necessary if gibbs sampling?
  }


// Copy Constructor:
DISTR_sndp_omega::DISTR_sndp_omega(const DISTR_sndp_omega & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }

// Operator:
const DISTR_sndp_omega & DISTR_sndp_omega::operator=(
                            const DISTR_sndp_omega & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


double DISTR_sndp_omega::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_sndp_omega::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]); // <--- c.f. linpred info
  }


double DISTR_sndp_omega::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }

  double x = (*response);
  double xi = (*worklin[1]);
  double om = exp(*linpred);
  double al = (*worklin[0]);

  double l;

  l = - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // calculate log-likelihood (possibly simplified)

  modify_worklin();

  return l;

  }


void DISTR_sndp_omega::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  if (counter==0)
    {
    set_worklin();
    }

    double x = (*response);
    double xi = (*worklin[1]);
    double om = exp(*linpred);
    double al = (*worklin[0]);
    double z = (x - xi)/om;

    double nu = (z*z - 1 - al*randnumbers::phi(al*z) / randnumbers::Phi2(al*z)*z)/om/om; // v

    *workingweight = (2 + al*al*sn_a2(al)); // w /pow(om,4) ???

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // loglikelihood

      }


  modify_worklin();

  }



void DISTR_sndp_omega::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("Link function (omega): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sndp_omega::update_end(void)
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
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }

// ONLY FOR THE MAIN PARAMETER:
void DISTR_sndp_xi::check_errors(void)
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

// Constructor:
DISTR_sndp_xi::DISTR_sndp_xi(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w) // <-----
  {
  family = "Skew Normal Distribution, Direct Parametrization - xi";
  outpredictor = true;
  outexpectation = true; // only true if main parameter
  predictor_name = "xi";

  // linpredminlimit=-10; // necessary for log-link
  // linpredmaxlimit=15; // necessary for log-link

  // updateIWLS = false; // only false and necessary if gibbs sampling?
  }


// Copy Constructor:
DISTR_sndp_xi::DISTR_sndp_xi(const DISTR_sndp_xi & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }

// Operator:
const DISTR_sndp_xi & DISTR_sndp_xi::operator=(
                            const DISTR_sndp_xi & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


// ONLY FOR THE MAIN PARAMETER:
void DISTR_sndp_xi::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
   if (*weight[2] == 0) // <---
     *deviance=0;
   else
     {

     double x = (*response[0]);
     double xi = (*linpred[2]);
     double om = exp(*linpred[1]);
     double al = (*linpred[0]);
     double l;

     l = -0.5*log(2*PI) - log(om) - 0.5*pow(x-xi, 2)/om/om + log(2*randnumbers::Phi2(al*(x-xi)/om)); // Compute log-likelihood (complete)

    *deviance = -2*l;
    }
  }


double DISTR_sndp_xi::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_sndp_xi::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2]); // <--- c.f. linpred info
  }

// ONLY FOR THE MAIN PARAMETER:
double DISTR_sndp_xi::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    return 0; // calculate cdf for observation using param!
    }

// ONLY FOR THE MAIN PARAMETER:
double DISTR_sndp_xi::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0; // usually set to be zero
    }




double DISTR_sndp_xi::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }
  double x = (*response);
  double xi = (*linpred);
  double om = exp(*worklin[1]);
  double al = (*worklin[0]);

  double l;

  l = - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // calculate log-likelihood (possibly simplified)

  modify_worklin();

  return l;

  }


void DISTR_sndp_xi::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  if (counter==0)
    {
    set_worklin();
    }

    double x = (*response);
    double xi = (*linpred);
    double om = exp(*worklin[1]);
    double al = (*worklin[0]);
    double z = (x - xi)/om;

    double nu = (z - al*randnumbers::phi(al*z) / randnumbers::Phi2(al*z))/om; // v

    *workingweight = (1 + al * al * sn_a0(al))/om/om; // w

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // loglikelihood

      }


  modify_worklin();

  }


// ONLY FOR THE MAIN PARAMETER:
void DISTR_sndp_xi::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+2])); // <--- depends on the number of parameters
  }


void DISTR_sndp_xi::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("Link function (xi): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sndp_xi::update_end(void)
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


DISTR_sncp_gamma::DISTR_sncp_gamma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w) // <-----
  {
  family = "Skew Normal Distribution, Centered Parametrization - Gamma";
  outpredictor = true;
  outexpectation = false; // only true if main parameter
  predictor_name = "gamma";

  linpredminlimit=-10; // necessary for log-link
  linpredmaxlimit=10; // necessary for log-link

  // updateIWLS = false; // only false and necessary if gibbs sampling?
  }


// Copy Constructor:
DISTR_sncp_gamma::DISTR_sncp_gamma(const DISTR_sncp_gamma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }

// Operator:
const DISTR_sncp_gamma & DISTR_sncp_gamma::operator=(
                            const DISTR_sncp_gamma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }




double DISTR_sncp_gamma::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_sncp_gamma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = 0.9952716*((*linpred[0])/(sqrt(1 + (*linpred[0])*(*linpred[0])))); // <--- c.f. linpred info
  }




double DISTR_sncp_gamma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }

  double x = (*response);
  double m = *worklin[1];
  double s = exp(*worklin[0]);
  double g = 0.9952716*((*linpred)/(sqrt(1 + (*linpred)*(*linpred))));
  double R = cbrt(2*g/(4-PI));
  double b2 = 2.0/PI;
  double al = R/sqrt(b2-(1-b2)*R*R);
  double d = al/sqrt(1 + al*al);
  double om = s/sqrt(1 - b2*d*d);
  double xi = m - sqrt(b2)*om*d;

  double l = 0;

  l = - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // calculate log-likelihood (possibly simplified)

  modify_worklin();

  return l;

  }


void DISTR_sncp_gamma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  if (counter==0)
    {
    set_worklin();
    }

  double x = (*response);
  double m = (*worklin[1]);
  double s = exp(*worklin[0]);
  double g = 0.9952716*((*linpred)/(sqrt(1 + (*linpred)*(*linpred))));
  double g1 = g;
  if(abs(g) < 0.001){
    g = 0.001;
  }
  double R = cbrt(2*g/(4-PI));
  double R1 = cbrt(2*g1/(4-PI));
  double b2 = 2.0/PI;
  double b = sqrt(b2);
  double al = R/sqrt(b2-(1-b2)*R*R);
  double al1 = R1/sqrt(b2-(1-b2)*R1*R1);
  double d = al/sqrt(1 + al*al);
  double d2 = d*d;
  double om = s/sqrt(1 - b2*d2);
  double xi = m - b*om*d;
  double z = (x - xi) / om;
  double zt1 = randnumbers::phi(al*z) / randnumbers::Phi2(al*z);
  double a0 = sn_a0(al);
  double a1 = sn_a1(al);
  double a2 = sn_a2(al);

  double nu = 0;

      nu = (((al * b * zt1 - z * b + b2 * d / (1 - b2 * d2) * (al * (b * d - z) * zt1 - z * b * d - 1 + z * z))
                 / pow(1 + al * al, 1.5) + z * zt1)
                 * (1 / sqrt(b2 - (1 - b2) * R * R) + (1 - b2) * R * R / pow(sqrt(b2 - (1 - b2) * R * R),3)) * cbrt(2 / (4 - PI) / pow(g,2)) / 3
                 * 0.9952716 / pow(sqrt(1 + (*linpred) * (*linpred)), 3));//v_g(x, m, s, g); // v
      *workingweight = ((a2 + 2 / pow(sqrt(1 + al * al), 3) * (b * al * a1 + (b * d * al * a1 - al * a2) * b2 * d / (1 - b2 * d2))
                       + (b2 * (al * al * a0 - 1) - 2 * b2 * b * d / (1 - b2 * d2) * (b * d2 * d + al * al *a1 + b * d * (1 - al * al * a0))
                          + b2 * b2 * d2 / pow(1 - b2 * d2, 2) * (2 - b2 * d2 * (1 + 2 * d2) + al * al * (a2 - 2 * b * d * a1 + b2 * d2 * a0)))
                      / pow(1 + al * al, 3))
                      * pow((1 / sqrt(b2 - (1 - b2) * R * R) + (1 - b2) * R * R / pow(sqrt(b2 - (1 - b2) * R * R),3)) * cbrt(2 / (4 - PI) / pow(g,2)) / 3
                      * 0.9952716 / pow(sqrt(1 + (*linpred) * (*linpred)), 3), 2)); //w_g(x, m, s, g); // w


    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += - log(om) - 0.5*pow(x-xi, 2)/om/om + log(2*randnumbers::Phi2(al1*(x-xi)/om)); // loglikelihood

      }


  modify_worklin();

  }



void DISTR_sncp_gamma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("Link function (gamma): custom\n"); // HIER custom ersetzen
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sncp_gamma::update_end(void)
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
    *pmu = 0.9952716*((*worklin)/(sqrt(1 + (*worklin)*(*worklin))));
//    double t = 0;
    }

  }


  // Constructor:
DISTR_sncp_sigma::DISTR_sncp_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w) // <-----
  {
  family = "Skew Normal Distribution, Centered Parametrization - sigma";
  outpredictor = true;
  outexpectation = false; // only true if main parameter
  predictor_name = "sigma";

   linpredminlimit=-10; // necessary for log-link
   linpredmaxlimit=15; // necessary for log-link

  // updateIWLS = false; // only false and necessary if gibbs sampling?
  }


// Copy Constructor:
DISTR_sncp_sigma::DISTR_sncp_sigma(const DISTR_sncp_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }

// Operator:
const DISTR_sncp_sigma & DISTR_sncp_sigma::operator=(
                            const DISTR_sncp_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }



double DISTR_sncp_sigma::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_sncp_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]); // <--- c.f. linpred info
  }



double DISTR_sncp_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }

  double x = (*response);
  double m = (*worklin[1]);
  double s = exp(*linpred);
  double g = 0.9952716*((*worklin[0])/(sqrt(1 + (*worklin[0])*(*worklin[0]))));
  double R = cbrt(2*g/(4-PI));
  double b2 = 2.0/PI;
  double al = R/sqrt(b2-(1-b2)*R*R);
  double d = al/sqrt(1 + al*al);
  double om = s/sqrt(1 - b2*d*d);
  double xi = m - sqrt(b2)*om*d;

  double l = 0;

  l = - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // calculate log-likelihood (possibly simplified)

  modify_worklin();

  return l;

  }


void DISTR_sncp_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  if (counter==0)
    {
    set_worklin();
    }

    double x = (*response);
    double m = (*worklin[1]);
    double s = exp(*linpred);
    double g = 0.9952716*((*worklin[0])/(sqrt(1 + (*worklin[0])*(*worklin[0]))));
    double R = cbrt(2*g/(4-PI));
    double b2 = 2.0/PI;
    double b = sqrt(b2);
    double al = R/sqrt(b2-(1-b2)*R*R);
    double d = al/sqrt(1 + al*al);
    double d2 = d*d;
    double om = s/sqrt(1 - b2*d2);
    double xi = m - b*om*d;
    double z = (x - xi) / om;

    double nu = (z * z - b * d * z - 1 - al * randnumbers::phi(al*z) / randnumbers::Phi2(al*z) * (z - b * d));//v_s(x, m, s, g);

    *workingweight = b2 * d2 * (1 + al * al * sn_a0(al) - 2 * (1 + 2 * al * al) / (1 + al * al)) - 2 * b * d * al * al * sn_a1(al) + 2 + al * al * sn_a2(al);
    //(b2 * d * (d * (1 + al * al * sn_a0(al) - 2 * (1 - al * al) / (1 + al * al)) - 2 * al * al * sn_a1(al)) + 2 + al * al * sn_a2(al));//w_s(x, m, s, g);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // loglikelihood

      }


  modify_worklin();

  }




void DISTR_sncp_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sncp_sigma::update_end(void)
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
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }


// ONLY FOR THE MAIN PARAMETER:
void DISTR_sncp_mu::check_errors(void)
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

// Constructor:
DISTR_sncp_mu::DISTR_sncp_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,2,w) // <-----
  {
  family = "Skew Normal Distribution, Centered Parametrization - mu";
  outpredictor = true;
  outexpectation = true; // only true if main parameter
  predictor_name = "mu";

  // linpredminlimit=-10; // necessary for log-link
  // linpredmaxlimit=15; // necessary for log-link

  // updateIWLS = false; // only false and necessary if gibbs sampling?
  }


// Copy Constructor:
DISTR_sncp_mu::DISTR_sncp_mu(const DISTR_sncp_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {

  }

// Operator:
const DISTR_sncp_mu & DISTR_sncp_mu::operator=(
                            const DISTR_sncp_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  return *this;
  }


// ONLY FOR THE MAIN PARAMETER:
void DISTR_sncp_mu::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {
   if (*weight[2] == 0) // <---
     *deviance=0;
   else
     {
     double x = (*response[0]);
     double m = (*linpred[2]);
     double s = exp(*linpred[1]);
     double g = 0.9952716*((*linpred[0])/(sqrt(1 + (*linpred[0])*(*linpred[0]))));
     double R = cbrt(2*g/(4-PI));
     double b2 = 2.0/PI;
     double al = R/sqrt(b2-(1-b2)*R*R);
     double d = al/sqrt(1 + al*al);
     double om = s/sqrt(1 - b2*d*d);
     double xi = m - sqrt(b2)*om*d;

     double l;

     l = -0.5*log(2*PI) - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // Compute log-likelihood (complete)

    *deviance = -2*l;
    }
  }


double DISTR_sncp_mu::get_intercept_start(void)
  {
  return 0;
  }

void DISTR_sncp_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2]); // <--- c.f. linpred info
  }

// ONLY FOR THE MAIN PARAMETER:
double DISTR_sncp_mu::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    return 0; // calculate cdf for observation using param!
    }

// ONLY FOR THE MAIN PARAMETER:
double DISTR_sncp_mu::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0; // usually set to be zero
    }




double DISTR_sncp_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {
  if (counter==0)
    {
    set_worklin();
    }
  double x = (*response);
  double m = (*linpred);
  double s = exp(*worklin[1]);
  double g = 0.9952716*((*worklin[0])/(sqrt(1 + (*worklin[0])*(*worklin[0]))));
  double R = cbrt(2*g/(4-PI));
  double b2 = 2.0/PI;
  double al = R/sqrt(b2-(1-b2)*R*R);
  double d = al/sqrt(1 + al*al);
  double om = s/sqrt(1 - b2*d*d);
  double xi = m - sqrt(b2)*om*d;

  double l;

  l = - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // calculate log-likelihood (possibly simplified)

  modify_worklin();

  return l;

  }


void DISTR_sncp_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  if (counter==0)
    {
    set_worklin();
    }
  double x = (*response);
  double m = (*linpred);
  double s = exp(*worklin[1]);
  double g = 0.9952716*((*worklin[0])/(sqrt(1 + (*worklin[0])*(*worklin[0]))));
  double R = cbrt(2*g/(4-PI));
  double b2 = 2.0/PI;
  double al = R/sqrt(b2-(1-b2)*R*R);
  double d = al/sqrt(1 + al*al);
  double om = s/sqrt(1 - b2*d*d);
  double xi = m - sqrt(b2)*om*d;
  double z = (x - xi) / om;

    double nu = z / om - al / om * randnumbers::phi(al*z) / randnumbers::Phi2(al*z);//v_m(x, m, s, g);

    *workingweight = (1 + al * al * sn_a0(al)) / om / om;//w_m(x, m, s, g);

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like += - log(om) - 0.5*pow(x-xi, 2)/om/om+log(2*randnumbers::Phi2(al*(x-xi)/om)); // loglikelihood

      }


  modify_worklin();

  }


// ONLY FOR THE MAIN PARAMETER:
void DISTR_sncp_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+2])); // <--- depends on the number of parameters
  }


void DISTR_sncp_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("Link function (mu): identitiy\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_sncp_mu::update_end(void)
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

//------------------------------------------------------------------------------
//--------- CLASS: DISTR_gaussiancopula_binary_dagum_latent -------------------
//------------------------------------------------------------------------------
void DISTR_gaussiancopula_binary_dagum_latent::check_errors(void)
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
           // cout << "y2: " << *workresp << endl;
            if(((*workresp)!=0.0) && ((*workresp)!=1.0)) {
                errors=true;
                errormessages.push_back("ERROR: response of latent margin has to be zero or one\n");
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


DISTR_gaussiancopula_binary_dagum_latent::DISTR_gaussiancopula_binary_dagum_latent(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Gaussiancopula latent marginal - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  responseorig = response;
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  check_errors();
  }


DISTR_gaussiancopula_binary_dagum_latent::DISTR_gaussiancopula_binary_dagum_latent(const DISTR_gaussiancopula_binary_dagum_latent & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  responseorig = nd.responseorig;
  response2p = nd.response2p;
  response2 = nd.response2;
  }


const DISTR_gaussiancopula_binary_dagum_latent & DISTR_gaussiancopula_binary_dagum_latent::operator=(
                            const DISTR_gaussiancopula_binary_dagum_latent & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  responseorig = nd.responseorig;
  response2p = nd.response2p;
  response2 = nd.response2;
  return *this;
  }

double DISTR_gaussiancopula_binary_dagum_latent::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_binary_dagum_latent::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = randnumbers::Phi2(*linpred[0]);
  }

void DISTR_gaussiancopula_binary_dagum_latent::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gaussiancopula_binary_dagum_latent::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    response2p++;

  }



void DISTR_gaussiancopula_binary_dagum_latent::update(void)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = mu;
  // *worklin[2] = linear predictor of rho equation
  // *worktransformlin[2] = rho;

  double * workresp = response.getV();
  double * workresporig =responseorig.getV();
  double * weightwork = weight.getV();

  double * worklin_current;
  if (linpred_current==1)
    worklin_current = linearpred1.getV();
  else
    worklin_current = linearpred2.getV();

  set_worklin();

  unsigned i;
  for(i=0;i<nrobs;i++,worklin_current++,workresp++,weightwork++,
           response2p++,workresporig++,worktransformlin[3]++,worktransformlin[2]++,worktransformlin[1]++,worktransformlin[0]++)
    {

    if (*weightwork != 0)
      {
      if(optionsp->samplesel && *workresporig == 0)
        {
        // double x = randnumbers::uniform();
        *workresp = trunc_normal2(-20,0,*worklin_current, 1);
        }
      else
        {
        double Fdagum = pow( ( 1 + pow( *response2p / *worktransformlin[1], - *worktransformlin[2] ) ), - *worktransformlin[0] );
        double help2 = sqrt(1-(*worktransformlin[3])*(*worktransformlin[3]));
        double help3 = randnumbers::invPhi2(Fdagum);

        double help = randnumbers::Phi2( (- *worklin_current - (*worktransformlin[3])*help3) / (help2));
        double x = randnumbers::uniform();
        double xstar;
        if(*workresporig==0)
          xstar = x*help;
        else
          xstar = x*(1-help) + help;


        *workresp = randnumbers::invPhi2(xstar)*help2 + (*worktransformlin[3])*help3 + *worklin_current;
        }
      }
    }
//  ofstream out1("c://temp//samplesel3//response.raw");
//  response.prettyPrint(out1);
//  out1.close();
  }

double DISTR_gaussiancopula_binary_dagum_latent::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = mu;
  // *worklin[2] = linear predictor of rho equation
  // *worktransformlin[2] = rho;
  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double l;



  l = -0.5 * pow(((*response)-mu),2) ;

  if((optionsp->samplesel) && (*response > 0))
    {
    double rho2 = pow((*worktransformlin[3]),2);
    double oneminusrho2 = 1- rho2;
    double v = pow( ( 1 + pow( *response2p / *worktransformlin[1], - *worktransformlin[2] ) ), - *worktransformlin[0] );
 //   double u = randnumbers::Phi2(*response-mu);
    double phiinvu = *response-mu;

    l += ((*worktransformlin[3])/(oneminusrho2))*phiinvu*( randnumbers::invPhi2(v) - 0.5*(*worktransformlin[3])*phiinvu );
    }
  modify_worklin();

  return l;

  }


void DISTR_gaussiancopula_binary_dagum_latent::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = mu;
  // *worklin[2] = linear predictor of rho equation
  // *worktransformlin[2] = rho;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double rho2;
  double oneminusrho2 = 0.0;
  double v;
  double u;
  double phiinvu = 0.0;
  double phiinvv = 0.0;
  double nu;
  double dphiinvu;
  double ddphiinvu;
  double du;
  double ddu;


  nu = *response-mu ;
  *workingweight = 1;

  if((!optionsp->samplesel) || (*response > 0))
    {
    rho2 = pow((*worktransformlin[3]),2);
    oneminusrho2 = 1- rho2;
    v = pow( ( 1 + pow( *response2p / *worktransformlin[1], - *worktransformlin[2] ) ), - *worktransformlin[0] );
    u = randnumbers::Phi2(*response-mu);
    phiinvu = *response-mu;
    phiinvv = randnumbers::invPhi2(v);
    dphiinvu = pow(2*PI,0.5)*exp(0.5*pow(phiinvu,2));
    ddphiinvu = 2*PI*phiinvu/(exp(-0.5*pow(phiinvu,2))*exp(-0.5*pow(phiinvu,2)));
    du = -exp(-0.5*(*response-mu)*(*response-mu))/pow(2*PI,0.5);
    ddu = (*response-mu)*du;

    nu += ((*worktransformlin[3])/oneminusrho2)*dphiinvu*du*( phiinvv - (*worktransformlin[3])*phiinvu );
    *workingweight += -((*worktransformlin[3])/oneminusrho2)*( phiinvv - (*worktransformlin[3])*phiinvu )*(ddphiinvu*du*du + dphiinvu*ddu) +
                            ((*worktransformlin[3])/oneminusrho2)*(*worktransformlin[3])*dphiinvu*dphiinvu*du*du;
    }

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
    like += -0.5 * pow(((*response)-mu),2) ;

    if((!optionsp->samplesel) || (*response > 0))
      {
      like += ((*worktransformlin[3])/(oneminusrho2))*phiinvu*( phiinvv - 0.5*(*worktransformlin[3])*phiinvu );
      }

    }
  modify_worklin();

  }


void DISTR_gaussiancopula_binary_dagum_latent::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = randnumbers::Phi2((*linpred[predstart_mumult]));
  }


void DISTR_gaussiancopula_binary_dagum_latent::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_binary_dagum_latent::update_end(void)
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
    }

  }


//--------------------------------------------------------------------------------
//----------------- CLASS: DISTR_gaussiancopula_binary_dagum_p -------------------
//--------------------------------------------------------------------------------


DISTR_gaussiancopula_binary_dagum_p::DISTR_gaussiancopula_binary_dagum_p(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Gaussian Copula - Dagum p";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "p";
    linpredminlimit=-10;
  linpredmaxlimit=15;

  }

DISTR_gaussiancopula_binary_dagum_p::DISTR_gaussiancopula_binary_dagum_p(const DISTR_gaussiancopula_binary_dagum_p & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  }


const DISTR_gaussiancopula_binary_dagum_p & DISTR_gaussiancopula_binary_dagum_p::operator=(
                            const DISTR_gaussiancopula_binary_dagum_p & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  return *this;
  }

void DISTR_gaussiancopula_binary_dagum_p::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = workingresponse2p->getV();

  }



void DISTR_gaussiancopula_binary_dagum_p::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_gaussiancopula_binary_dagum_p::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_binary_dagum_p::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[1]);
  }

double DISTR_gaussiancopula_binary_dagum_p::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[2] = linear predictor of latent equation
  // *worktransformlin[0] = b
  // *worktransformlin[1] = a
  // *worktransformlin[3]=rho

  if (counter==0)
    {
    set_worklin();
    }

    double p = exp((*linpred));
    double bcurrent = (*worktransformlin[0]);
    double acurrent = (*worktransformlin[1]);
    double respdivb = (*response) / bcurrent;
    double hilfs = pow(respdivb,acurrent);
    double hilfs2 = pow(respdivb, -acurrent);


    double u = pow((1 + hilfs2), -p);
//    double v = randnumbers::Phi2(*response2p-(*worklin[2]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = *response2p-(*worklin[2]);
    double orho = 1 - pow((*worktransformlin[3]), 2);
    double l;

    l = (*worktransformlin[3]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[3]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho +
                log(p) + acurrent*p*log((*response)) - acurrent*p*log(bcurrent)-p*log(1+hilfs);

  modify_worklin();

  return l;

  }

void DISTR_gaussiancopula_binary_dagum_p::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

    // *worklin[2] = linear predictor of latent equation
  // *worktransformlin[0] = b
  // *worktransformlin[1] = a
  // *worktransformlin[3]=rho

  if (counter==0)
    {
    set_worklin();
    }

    double p = exp((*linpred));
    double bcurrent = (*worktransformlin[0]);
    double acurrent = (*worktransformlin[1]);
    double respdivb = (*response) / bcurrent;
    double hilfs = pow(respdivb,acurrent);
    double hilfs2 = pow(respdivb, -acurrent);


    double u = pow((1 + hilfs2), -p);
    // double v = randnumbers::Phi2(*response2p-(*worklin[2]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = *response2p-(*worklin[2]);
    double orho = 1 - pow((*worktransformlin[3]), 2);

    double dphiinvu = pow(2*PI, 0.5) / exp(-0.5*pow(phinvu, 2));
    double ddphiinvu = phinvu * 2 * PI / pow(exp(-0.5*pow(phinvu, 2)), 2);

    double ybpma = hilfs2;//pow(*response/(*worktransformlin[0]),-(*worktransformlin[1]));
    double ybpmap = u;//pow(1+ybpma,-p);
    double dF = -p*log(1+ybpma)*ybpmap;
    double ddF = p*log(1+ybpma)*ybpmap*(p*log(1+ybpma)-1);

//    double d1 =  -p * pow(1 +  hilfs2, -p) * log(1 +  hilfs2);
//    double d2 = d1 + p * p * pow(1 +  hilfs2, -p) * pow(log(1 +  hilfs2), 2);

    double nu = 1 + acurrent*p*log((*response)) - acurrent*p*log(bcurrent)
                -p*log(1+hilfs) +
                ((*worktransformlin[3]) * dF * dphiinvu / orho) * (phinvv - (*worktransformlin[3]) * phinvu);

    *workingweight = 1  -
                    ((*worktransformlin[3]) * (dF * dF * ddphiinvu + dphiinvu*ddF) / orho) * (phinvv - (*worktransformlin[3]) * phinvu) +
                    (*worktransformlin[3]) * (*worktransformlin[3]) * pow(dF * dphiinvu, 2) / orho;

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  (*worktransformlin[3]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[3]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho +
                log(p) + acurrent*p*log((*response)) - acurrent*p*log(bcurrent)-p*log(1+hilfs);

      }


  modify_worklin();

  }

void DISTR_gaussiancopula_binary_dagum_p::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_gaussiancopula_binary_dagum_p::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (p): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_binary_dagum_p::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }




//--------------------------------------------------------------------------------
//----------------- CLASS: DISTR_gaussiancopula_binary_dagum_b -------------------
//--------------------------------------------------------------------------------


DISTR_gaussiancopula_binary_dagum_b::DISTR_gaussiancopula_binary_dagum_b(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Gaussian Copula - Dagum b";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "b";
  linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_gaussiancopula_binary_dagum_b::DISTR_gaussiancopula_binary_dagum_b(const DISTR_gaussiancopula_binary_dagum_b & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  }


const DISTR_gaussiancopula_binary_dagum_b & DISTR_gaussiancopula_binary_dagum_b::operator=(
                            const DISTR_gaussiancopula_binary_dagum_b & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  return *this;
  }

void DISTR_gaussiancopula_binary_dagum_b::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = workingresponse2p->getV();

  }



void DISTR_gaussiancopula_binary_dagum_b::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_gaussiancopula_binary_dagum_b::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_binary_dagum_b::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp((*linpred[2]));
  }

double DISTR_gaussiancopula_binary_dagum_b::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

    // *worklin[2] = linear predictor of latent equation
  // *worktransformlin[0] = p
  // *worktransformlin[1] = a
  // *worktransformlin[3]=rho

  if (counter==0)
    {
    set_worklin();
    }

    double b = exp((*linpred));
    double respdivb = (*response) / b;
    double acurrent = (*worktransformlin[1]);
    double hilfs = pow(respdivb,acurrent);
    // double hilfs2 = pow(respdivb, -acurrent);
    double pcurrent = (*worktransformlin[0]);

    double u = pow((1 + pow((respdivb), -acurrent)), -pcurrent);
    // double v = randnumbers::Phi2(*response2p-(*worklin[2]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = *response2p-(*worklin[2]);
    double orho = 1 - pow((*worktransformlin[3]), 2);

    double l;

    l = (*worktransformlin[3]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[3]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho
                - acurrent*pcurrent*log(b) - (pcurrent+1)*log(1+hilfs);

  modify_worklin();

  return l;

  }

void DISTR_gaussiancopula_binary_dagum_b::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

    // *worklin[2] = linear predictor of latent equation
  // *worktransformlin[0] = p
  // *worktransformlin[1] = a
  // *worktransformlin[3]=rho;

  if (counter==0)
    {
    set_worklin();
    }

  /*  cout << worktransformlin.size() << endl;
    for(unsigned i=0; i< worktransformlin.size(); i++)
      cout << *worktransformlin[i] << endl;
    cout << *linpred << endl;*/
    double b = exp((*linpred));
   // cout << b << endl;
   // cout << *response << endl;
    double respdivb = (*response) / b;
  //  cout << respdivb << endl;
  //  cout << *worktransformlin[4] << endl;
    double acurrent = (*worktransformlin[1]);
    double hilfs = pow(respdivb,acurrent);
    double hilfs2 = pow(respdivb, -acurrent);
    double pcurrent = (*worktransformlin[0]);

    double u = pow((1 + pow((respdivb), -acurrent)), -pcurrent);
    // double v = randnumbers::Phi2(*response2p-(*worklin[2]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = *response2p-(*worklin[2]);
    double orho = 1 - pow((*worktransformlin[3]), 2);

    double dphiinvu = pow(2*PI, 0.5) / exp(-0.5*pow(phinvu, 2));
    double ddphiinvu = phinvu * 2 * PI / pow(exp(-0.5*pow(phinvu, 2)), 2);


    // double d1 =  -acurrent * pcurrent * pow(1 +  hilfs2, -pcurrent-1) * hilfs2;
    double d2 = -acurrent * acurrent * pow(1+ hilfs2, -pcurrent-1) * pcurrent * hilfs2 +
                acurrent * acurrent * pow(1+ hilfs2, -pcurrent-2) * (pcurrent+1) * pcurrent * pow(hilfs2, 2);

    double ybpa = pow(*response/b,(*worktransformlin[1]));
    double ybpmap = pow(1+pow(*response/b,-(*worktransformlin[1])),-(*worktransformlin[0])-1);
    // double ybpmap2 = pow(1+pow(*response/b,-(*worktransformlin[1])),-(*worktransformlin[0]));
    double dF = -(*worktransformlin[1])*(*worktransformlin[0])*ybpmap/ybpa;
    // double ddF = (*worktransformlin[1])*(*worktransformlin[1])*(*worktransformlin[0])*((*worktransformlin[0])-ybpa)*ybpmap2/((ybpa+1)*(ybpa+1));

    double nu = acurrent - ((pcurrent+1)*acurrent)/(1+hilfs)+
                ((*worktransformlin[3]) * dF * dphiinvu / orho) * (phinvv - (*worktransformlin[3]) * phinvu);

    *workingweight = (pow(acurrent,2)*pcurrent)/(pcurrent+2)-
                    ((*worktransformlin[3]) * (dF * dF * ddphiinvu + dphiinvu*d2) / orho) * (phinvv - (*worktransformlin[3]) * phinvu)
                    + (*worktransformlin[3]) * (*worktransformlin[3]) * pow(dF * dphiinvu, 2) / orho;

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  (*worktransformlin[3]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[3]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho
                - acurrent*pcurrent*log(b) - (pcurrent+1)*log(1+hilfs);

      }



  modify_worklin();

  }

void DISTR_gaussiancopula_binary_dagum_b::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 0;
  }


void DISTR_gaussiancopula_binary_dagum_b::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (b): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_binary_dagum_b::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }


//------------------------------------------------------------------------------
//---------------- CLASS: DISTR_gaussiancopula_binary_dagum_a ------------------
//------------------------------------------------------------------------------


DISTR_gaussiancopula_binary_dagum_a::DISTR_gaussiancopula_binary_dagum_a(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Gaussian Copula - Dagum a";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "a";
 //   linpredminlimit=-10;
 // linpredmaxlimit=15;
  }


DISTR_gaussiancopula_binary_dagum_a::DISTR_gaussiancopula_binary_dagum_a(const DISTR_gaussiancopula_binary_dagum_a & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
    response2p = nd.response2p;
    workingresponse2p = nd.workingresponse2p;
  }


const DISTR_gaussiancopula_binary_dagum_a & DISTR_gaussiancopula_binary_dagum_a::operator=(
                            const DISTR_gaussiancopula_binary_dagum_a & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  return *this;
  }

void DISTR_gaussiancopula_binary_dagum_a::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = workingresponse2p->getV();

  }



void DISTR_gaussiancopula_binary_dagum_a::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }

double DISTR_gaussiancopula_binary_dagum_a::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_binary_dagum_a::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = exp(*linpred[3]);
  }

 double DISTR_gaussiancopula_binary_dagum_a::pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)
    {
    return 0;
    }

double DISTR_gaussiancopula_binary_dagum_a::cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    return 0;
    }

//double DISTR_gaussiancopula_binary_dagum_a::compute_quantile_residual_mult(vector<double *> response,
//                                             vector<double *> param,
//                                             vector<double *> weight,
//                                             vector<datamatrix *> aux)
//  {
//  double u_est;
//  u_est = cdf_mult(response,param,weight,aux);
//  double res_est = randnumbers::invPhi2(u_est);
//  return res_est;
//  }

double DISTR_gaussiancopula_binary_dagum_a::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;

  if (counter==0)
    {
    set_worklin();
    }

    double a = exp((*linpred));
    double respdivb = (*response) / (*worktransformlin[1]);
    double hilfs = pow(respdivb,a);
    double pcurrent = (*worktransformlin[0]);

    double u = pow((1 + pow((respdivb), -a)), -pcurrent);
    // double v = randnumbers::Phi2(*response2p-(*worklin[2]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = *response2p-(*worklin[2]);
    double orho = 1 - pow((*worktransformlin[3]), 2);
    double l;

    l = (*worktransformlin[3]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[3]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho +
                log(a) +(a*pcurrent-1)*log((*response)) - a*pcurrent*log((*worktransformlin[1]))
                - (pcurrent+1)*log(1+hilfs);


    modify_worklin();

    return l;

  }


void DISTR_gaussiancopula_binary_dagum_a::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma2 equation
  // *worktransformlin[0] = sigma2;


  if (counter==0)
    {
    set_worklin();
    }

    double a = exp((*linpred));
    double respdivb = (*response) / (*worktransformlin[1]);
    double hilfs = pow(respdivb,a);
    double pcurrent = (*worktransformlin[0]);

    double u = pow((1 + pow((respdivb), -a)), -pcurrent);
    // double v = randnumbers::Phi2(*response2p-(*worklin[2]));

    double phinvu = randnumbers::invPhi2(u);
    double phinvv = *response2p-(*worklin[2]);
    double orho = 1 - pow((*worktransformlin[3]), 2);

    double dphiinvu = pow(2*PI, 0.5) / exp(-0.5*pow(phinvu, 2));
    double ddphiinvu = phinvu * 2 * PI / pow(exp(-0.5*pow(phinvu, 2)), 2);


    double d1 =  pcurrent * pow(1 +  pow(respdivb, -a), -pcurrent-1) * a * pow(respdivb, -a) * log(respdivb);
    double d2 = d1 - pcurrent * pow(1+pow(respdivb, -a), -pcurrent-1) * a * a * log(respdivb) * log(respdivb) * pow(respdivb, -a) +
                pcurrent * (pcurrent + 1) *pow(1+pow(respdivb, -a), -pcurrent-2) * pow((pow(respdivb, -a) * log(respdivb) * a), 2);

    // double ybpa = pow((*response/(*worktransformlin[1])),a);
    double lyb = log(*response/(*worktransformlin[1]));
    double ybpma = pow(*response/(*worktransformlin[1]),-a);
    double ybpmap = pow(1+ybpma,-(*worktransformlin[0])-1);
    // double ybpmap2 = pow(1+ybpma,-(*worktransformlin[0]));
    double dF = (*worktransformlin[0])*lyb*a*ybpma*ybpmap;
    // double ddF = a*(*worktransformlin[0])*lyb*ybpmap2*(a*lyb*((*worktransformlin[0])-ybpa)+ybpa+1)/pow(1+ybpa,2);

    double nu = 1 + a*pcurrent*log(respdivb)
                - ((pcurrent+1)*a*hilfs*log(respdivb))/(1+hilfs) +
                ((*worktransformlin[3]) * dF * dphiinvu / orho) * (phinvv - (*worktransformlin[3]) * phinvu);

    *workingweight = 1 + ((pcurrent+1)*pow(a,2)*hilfs*pow(log(respdivb),2))/pow((1+hilfs),2) -
                    ((*worktransformlin[3]) * (dF * dF * ddphiinvu + dphiinvu*d2) / orho) * (phinvv - (*worktransformlin[3]) * phinvu)
                    + (*worktransformlin[3]) * (*worktransformlin[3]) * pow(dF * dphiinvu, 2) / orho;

    if((*workingweight) <= 0)
        *workingweight = 0.0001;

    *workingresponse = *linpred + nu/(*workingweight);

    if (compute_like)
      {

        like +=  (*worktransformlin[3]) * phinvu * phinvv / orho - 0.5 * pow((*worktransformlin[3]), 2) * (pow(phinvu, 2) + pow(phinvv, 2)) / orho +
                log(a) +(a*pcurrent-1)*log((*response)) - a*pcurrent*log((*worktransformlin[1]))
                - (pcurrent+1)*log(1+hilfs);

      }


  modify_worklin();

  }


void DISTR_gaussiancopula_binary_dagum_a::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+3]));
  }


void DISTR_gaussiancopula_binary_dagum_a::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (a): exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_binary_dagum_a::update_end(void)
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
//    double t = *pmu;
    *pmu = exp(*worklin);
//    double t = 0;
    }

  }

//------------------------------------------------------------------------------
//---------- CLASS: DISTR_gaussiancopula_binary_dagum_rho ----------------------
//------------------------------------------------------------------------------

void DISTR_gaussiancopula_binary_dagum_rho::check_errors(void)
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

DISTR_gaussiancopula_binary_dagum_rho::DISTR_gaussiancopula_binary_dagum_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,4,w)
  {
  family = "Gaussiancopula - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
  linpredminlimit=-100;
  linpredmaxlimit=100;

  }


DISTR_gaussiancopula_binary_dagum_rho::DISTR_gaussiancopula_binary_dagum_rho(const DISTR_gaussiancopula_binary_dagum_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  }


const DISTR_gaussiancopula_binary_dagum_rho & DISTR_gaussiancopula_binary_dagum_rho::operator=(
                            const DISTR_gaussiancopula_binary_dagum_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  return *this;
  }

void DISTR_gaussiancopula_binary_dagum_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[1] = *response[2] = *response[3] = first and continuous component of two dimensional reponse
   // *response[0] = workingresponse of binary component
   // *linpred[0] = eta_mu_latent
   // *linpred[1] = eta_sigma
   // *linpred[2] = eta_mu
   // *linpred[3] = eta_rho
   // *linpred[4] = eta_mu_1

   if (*weight[0]==0)
     *deviance=0;
   else
     {
     double rho = (*linpred[4])/pow(1+pow((*linpred[4]),2),0.5);
     double p = exp(*linpred[1]);
     double b = exp(*linpred[2]);
     double a = exp(*linpred[3]);
     double mu = (*linpred[0]);
     // double hilfs1 = 1-pow(rho,2);
     double u = pow(1+pow(*response[1]/b,-a),-p);
     double phiinvu = randnumbers::invPhi2(u);
     double phiinv = *response[0]-mu;
     double orho = 1-rho*rho;
     double l;

     if(*weight[4]==0)
       l = -0.5*log(2*PI)-(*response[0]-mu)*(*response[0]-mu);
     else
       {
       l = log(a)+log(p)-a*p*log(b)+(a*p-1)*log(*response[1])-(p+1)*log(1+pow(*response[1]/b,a))
           -0.5*log(2*PI)-0.5*log(orho) + (rho*phiinvu*phiinv-0.5*rho*rho*(phiinv*phiinv+phiinvu*phiinvu))/orho;
       }


    *deviance = -2*l;
    }

  }

double DISTR_gaussiancopula_binary_dagum_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_binary_dagum_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = (*linpred[4]);
  *param = arg/pow(1+pow(arg,2),0.5);
  }

void DISTR_gaussiancopula_binary_dagum_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

   response2p = workingresponse2p->getV();
  }



void DISTR_gaussiancopula_binary_dagum_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_gaussiancopula_binary_dagum_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  if (counter==0)
    {
    set_worklin();
    }

  double rho;
  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);


  double p = *worktransformlin[1];
  double b = *worktransformlin[2];
  double a = *worktransformlin[3];
  double mu = *worklin[0];
//  double hilfs1 = 1-pow(rho,2);
  double u = pow(1+pow(*response/b,-a),-p);
  double phiinvu = randnumbers::invPhi2(u);
  double phiinv = *response2p-mu;
  double orho = 1-rho*rho;
  double l;

  l = -0.5*log(orho) + (rho*phiinvu*phiinv-0.5*rho*rho*(phiinv*phiinv+phiinvu*phiinvu))/orho;//log(a)+log(p)-a*p*log(b)+(a*p-1)*log(*response)-(p+1)*log(1+pow(*response/b,a))
        //   -0.5*log(2*PI)-(*response2p-mu)*(*response2p-mu)



  modify_worklin();

  return l;

  }

void DISTR_gaussiancopula_binary_dagum_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of latent equation
  // *worktransformlin[0] = latent;
  // *worklin[1] = linear predictor of sigma equation
  // *worktransformlin[1] = sigma;
  // *worklin[2] = linear predictor of mu equation
  // *worktransformlin[2] = mu;
  if (counter==0)
    {
    set_worklin();
    }

  double rho;
  double hilfs;

  if (*linpred <= -100) {
    rho  = -0.99995;
    hilfs = 100.05;
  }
  else if (*linpred >= 100) {
    rho  = 0.99995;
    hilfs = 100.05;
  }
  else {
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);
    hilfs = pow((1+pow((*linpred),2)),0.5);
  }

  double p = *worktransformlin[1];
  double b = *worktransformlin[2];
  double a = *worktransformlin[3];
  double mu = *worklin[0];
//  double hilfs1 = 1-pow(rho,2);
  double u = pow(1+pow(*response/b,-a),-p);
  double phiinvu = randnumbers::invPhi2(u);
  double phiinv = (*response2p)-mu;
  double orho = 1-rho*rho;

  double nu = rho * pow(orho, 0.5) + (hilfs + rho * (*linpred)) * (phiinvu * phiinv) - (*linpred) * (pow(phiinvu, 2) + pow(phiinv, 2));

  *workingweight = 1 - pow(rho, 4);

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    like += -0.5*log(orho) + (rho*phiinvu*phiinv-0.5*rho*rho*(phiinv*phiinv+phiinvu*phiinvu))/orho;// log(a)+log(p)-a*p*log(b)+(a*p-1)*log(*response)-(p+1)*log(1+pow(*response/b,a))
          // -0.5*log(2*PI)-(*response2p-mu)*(*response2p-mu)


    }


  modify_worklin();



  }

void DISTR_gaussiancopula_binary_dagum_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 2 * std::asin((*linpred[predstart_mumult+4]) / (pow(1 + pow((*linpred[predstart_mumult+4]), 2), 0.5))) / PI ;
  }

void DISTR_gaussiancopula_binary_dagum_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho): fisher z-transformation\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_binary_dagum_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin)/pow((1+pow((*worklin),2)),0.5);
    }

  }





//------------------------------------------------------------------------------
//--------- CLASS: DISTR_gaussiancopula_binary_normal_latent -------------------
//------------------------------------------------------------------------------
void DISTR_gaussiancopula_binary_normal_latent::check_errors(void)
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
           // cout << "y2: " << *workresp << endl;
            if(((*workresp)!=0.0) && ((*workresp)!=1.0)) {
                errors=true;
                errormessages.push_back("ERROR: response of latent margin has to be zero or one\n");
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


DISTR_gaussiancopula_binary_normal_latent::DISTR_gaussiancopula_binary_normal_latent(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,3,w)
  {
  family = "Gaussiancopula latent marginal - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  responseorig = response;
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  check_errors();
  }


DISTR_gaussiancopula_binary_normal_latent::DISTR_gaussiancopula_binary_normal_latent(const DISTR_gaussiancopula_binary_normal_latent & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  responseorig = nd.responseorig;
  response2p = nd.response2p;
  response2 = nd.response2;
  }


const DISTR_gaussiancopula_binary_normal_latent & DISTR_gaussiancopula_binary_normal_latent::operator=(
                            const DISTR_gaussiancopula_binary_normal_latent & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  responseorig = nd.responseorig;
  response2p = nd.response2p;
  response2 = nd.response2;
  return *this;
  }

double DISTR_gaussiancopula_binary_normal_latent::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_binary_normal_latent::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[0]);
  }

void DISTR_gaussiancopula_binary_normal_latent::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = response2.getV();

  }



void DISTR_gaussiancopula_binary_normal_latent::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    response2p++;

  }



void DISTR_gaussiancopula_binary_normal_latent::update(void)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = mu;
  // *worklin[2] = linear predictor of rho equation
  // *worktransformlin[2] = rho;

  double * workresp = response.getV();
  double * workresporig =responseorig.getV();
  double * weightwork = weight.getV();

  double * worklin_current;
  if (linpred_current==1)
    worklin_current = linearpred1.getV();
  else
    worklin_current = linearpred2.getV();

  set_worklin();

  unsigned i;
  for(i=0;i<nrobs;i++,worklin_current++,workresp++,weightwork++,
           response2p++,workresporig++,worktransformlin[2]++,worklin[1]++,worktransformlin[0]++)
    {

    if (*weightwork != 0)
      {
      if (*workresporig > 0)
        *workresp = trunc_normal2(0,20,*worklin_current+(*worktransformlin[2])*((*response2p)-(*worklin[1]))/(*worktransformlin[0]),pow(1-pow(*worktransformlin[2],2),0.5));
      else
        *workresp = trunc_normal2(-20,0,*worklin_current+(*worktransformlin[2])*((*response2p)-(*worklin[1]))/(*worktransformlin[0]),pow(1-pow(*worktransformlin[2],2),0.5));
      }
    }
  }

double DISTR_gaussiancopula_binary_normal_latent::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = mu;
  // *worklin[2] = linear predictor of rho equation
  // *worktransformlin[2] = rho;
  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double rho2 = pow((*worktransformlin[2]),2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -(1/(2*oneminusrho2))*( pow((((*response))-mu),2) -
                                 2*(*worktransformlin[2])*(((*response)-mu))*(((*response2p)-(*worktransformlin[1]))/(*worktransformlin[0])) );

  modify_worklin();

  return l;

  }


void DISTR_gaussiancopula_binary_normal_latent::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;
  // *worklin[1] = linear predictor of mu equation
  // *worktransformlin[1] = mu;
  // *worklin[2] = linear predictor of rho equation
  // *worktransformlin[2] = rho;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double rho2 = pow((*worktransformlin[2]),2);
  double oneminusrho2 = 1- rho2;


  double nu = (1/(oneminusrho2))*( (((*response))-mu) -
                                 ((*worktransformlin[2]))*(((*response2p)-(*worktransformlin[1]))/(*worktransformlin[0])) );

  *workingweight = 1/(oneminusrho2);

//  cout << "latent equation y1: " << *response << endl;
//  cout << "latent equation y2: " << *response2p << endl;
//  cout << "latent equation y2orig: " << *workresporig << endl;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    like += -(1/(2*oneminusrho2))*( pow((((*response))-mu),2) -
                                 2*(*worktransformlin[2])*(((*response)-mu))*(((*response2p)-(*worktransformlin[1]))/(*worktransformlin[0])) );

    }

  modify_worklin();

  }


void DISTR_gaussiancopula_binary_normal_latent::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult]));
  }


void DISTR_gaussiancopula_binary_normal_latent::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_binary_normal_latent::update_end(void)
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
    }

  }

//------------------------------------------------------------------------------
//---------- CLASS: DISTR_gaussiancopula_binary_normal_sigma ------------------
//------------------------------------------------------------------------------


DISTR_gaussiancopula_binary_normal_sigma::DISTR_gaussiancopula_binary_normal_sigma(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,3,w)
  {
  family = "Gaussiancopula Normal marginal - sigma";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "sigma";
  linpredminlimit=-10;
  linpredmaxlimit=15;

  }


DISTR_gaussiancopula_binary_normal_sigma::DISTR_gaussiancopula_binary_normal_sigma(const DISTR_gaussiancopula_binary_normal_sigma & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  workingresponse2p = nd.workingresponse2p;
  response2p = nd.response2p;
  }


const DISTR_gaussiancopula_binary_normal_sigma & DISTR_gaussiancopula_binary_normal_sigma::operator=(
                            const DISTR_gaussiancopula_binary_normal_sigma & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  workingresponse2p = nd.workingresponse2p;
  response2p = nd.response2p;
  return *this;
  }


double DISTR_gaussiancopula_binary_normal_sigma::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_binary_normal_sigma::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param =  exp((*linpred[1]));
  }
void DISTR_gaussiancopula_binary_normal_sigma::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = workingresponse2p->getV();

  }



void DISTR_gaussiancopula_binary_normal_sigma::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }


double DISTR_gaussiancopula_binary_normal_sigma::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = mu;
  // *worklin[1] = linear predictor of latent equation
  // *worktransformlin[1] = mu_latent;
  // *worklin[2] = linear predictor of rho equation
  // *worktransformlin[2] = rho;

  if (counter==0)
    {
    set_worklin();
    }

  double sigma = exp((*linpred));

  // double rho2 = pow((*worktransformlin[2]),2);
  double oneminusrho2 = 1- pow((*worktransformlin[2]),2);
  double l;

  l = -log(sigma) -(1/(2*oneminusrho2))*( pow((((*response))-(*worktransformlin[0])),2)/pow(sigma,2) -
                                 2*(*worktransformlin[2])*(((*response)-(*worktransformlin[0]))/(sigma))*(((*response2p)-(*worktransformlin[1]))/(1)) );
  modify_worklin();

  return l;

  }

void DISTR_gaussiancopula_binary_normal_sigma::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of mu equation
  // *worktransformlin[0] = mu;
  // *worklin[1] = linear predictor of latent equation
  // *worktransformlin[1] = mu_latent;
  // *worklin[2] = linear predictor of rho equation
  // *worktransformlin[2] = rho;

  if (counter==0)
    {
    set_worklin();
    }

  double sigma = exp((*linpred));

  double rho2 = pow((*worktransformlin[2]),2);
  double oneminusrho2 = 1- rho2;


  double nu = -1 + (1/oneminusrho2)*(pow(((*response)-(*worklin[0])),2))/pow(sigma,2)
                - ((*worktransformlin[2])/oneminusrho2)*(((*response)-(*worktransformlin[0]))/(sigma))*(((*response2p)-(*worklin[1]))/(1));



  *workingweight = 1+1/oneminusrho2;

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {
	like +=  -log(sigma) -(1/(2*oneminusrho2))*( pow((((*response))-(*worktransformlin[0])),2)/pow(sigma,2) -
                                 2*(*worktransformlin[2])*(((*response)-(*worktransformlin[0]))/(sigma))*(((*response2p)-(*worktransformlin[1]))/(1)) );
    }

  modify_worklin();

  }


void DISTR_gaussiancopula_binary_normal_sigma::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (sigma): log\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_binary_normal_sigma::update_end(void)
  {

  // helpmat1 stores sigma2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = exp(*worklin);
    }

  }

//------------------------------------------------------------------------------
//-------------- CLASS: DISTR_gaussiancopula_binary_normal_mu ------------------
//------------------------------------------------------------------------------

DISTR_gaussiancopula_binary_normal_mu::DISTR_gaussiancopula_binary_normal_mu(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,3,w)
  {
  family = "Gaussiancopula normal marginal - mu";
  outpredictor = true;
  outexpectation = true;
  predictor_name = "mu";
  if(o->forceIWLS)
    {
    updateIWLS = true;
    }
  else
    {
    updateIWLS = false;
    }
//    linpredminlimit=-10;
//  linpredmaxlimit=15;
  }


DISTR_gaussiancopula_binary_normal_mu::DISTR_gaussiancopula_binary_normal_mu(const DISTR_gaussiancopula_binary_normal_mu & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  }


const DISTR_gaussiancopula_binary_normal_mu & DISTR_gaussiancopula_binary_normal_mu::operator=(
                            const DISTR_gaussiancopula_binary_normal_mu & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  return *this;
  }

double DISTR_gaussiancopula_binary_normal_mu::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_binary_normal_mu::compute_param_mult(vector<double *>  linpred,double * param)
  {
  *param = (*linpred[2]);
  }


void DISTR_gaussiancopula_binary_normal_mu::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

  response2p = workingresponse2p->getV();

  }



void DISTR_gaussiancopula_binary_normal_mu::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_gaussiancopula_binary_normal_mu::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;
  // *worklin[1] = linear predictor of latent equation
  // *worktransformlin[1] = mu_latent;
  // *worklin[2] = linear predictor of rho equation
  // *worktransformlin[2] = rho;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double rho2 = pow((*worktransformlin[2]),2);
  double oneminusrho2 = 1- rho2;
  double l;


  l = -(1/(2*oneminusrho2))*( pow((((*response))-mu),2)/pow((*worktransformlin[0]),2) -
                                 2*(*worktransformlin[2])*(((*response)-mu)/((*worktransformlin[0])))*(((*response2p)-(*worktransformlin[1]))/(1)) );

  modify_worklin();

  return l;

  }


void DISTR_gaussiancopula_binary_normal_mu::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  // *worklin[0] = linear predictor of sigma equation
  // *worktransformlin[0] = sigma;
  // *worklin[1] = linear predictor of latent equation
  // *worktransformlin[1] = mu_latent;
  // *worklin[2] = linear predictor of rho equation
  // *worktransformlin[2] = rho;

  if (counter==0)
    {
    set_worklin();
    }

  double mu = (*linpred);
  double rho2 = pow((*worktransformlin[2]),2);
  double oneminusrho2 = 1- rho2;

  double nu = (1/(oneminusrho2))*( (((*response))-mu)/pow((*worktransformlin[0]),2) -
                                 ((*worktransformlin[2])/(*worktransformlin[0]))*(((*response2p)-(*worktransformlin[1]))/(1)) );

  *workingweight = 1/(oneminusrho2*pow((*worktransformlin[0]),2));

  *workingresponse = *linpred + nu/(*workingweight);

//  cout << "mu equation y1: " << *response << endl;
//  cout << "mu equation y2: " << *response2p << endl;

  if (compute_like)
    {
    like += -(1/(2*oneminusrho2))*( pow((((*response))-mu),2)/pow((*worktransformlin[0]),2) -
                                 2*(*worktransformlin[2])*(((*response)-mu)/((*worktransformlin[0])))*(((*response2p)-(*worktransformlin[1]))/(1)) );
    }


  modify_worklin();

  }


void DISTR_gaussiancopula_binary_normal_mu::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = ((*linpred[predstart_mumult+2]));
  }


void DISTR_gaussiancopula_binary_normal_mu::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function (mu): identity\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }

void DISTR_gaussiancopula_binary_normal_mu::update(void)
  {

   unsigned i;

  // double help;

  double * worktransformlinr;
  double * worktransformlins;
  double * workweight;

//  cout << "mu equation y2: " << *response2p << endl;

  worktransformlinr = distrp[2]->helpmat1.getV();
  worktransformlins = distrp[0]->helpmat1.getV();
  workweight = workingweight.getV();

  for (i=0;i<nrobs;i++,worktransformlinr++,worktransformlins++,workweight++)
    {
        *workweight = 1/((1 - pow((*worktransformlinr),2))*pow((*worktransformlins),2));
    }

  }

void DISTR_gaussiancopula_binary_normal_mu::update_end(void)
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
    }

  }

//------------------------------------------------------------------------------
//---------- CLASS: DISTR_gaussiancopula_binary_normal_rho ---------------------
//------------------------------------------------------------------------------
void DISTR_gaussiancopula_binary_normal_rho::check_errors(void)
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

DISTR_gaussiancopula_binary_normal_rho::DISTR_gaussiancopula_binary_normal_rho(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_gamlss(o,r,3,w)
  {
  family = "Gaussiancopula - rho";

  outpredictor = true;
  outexpectation = false;
  predictor_name = "rho";
  linpredminlimit=-100;
  linpredmaxlimit=100;

  }


DISTR_gaussiancopula_binary_normal_rho::DISTR_gaussiancopula_binary_normal_rho(const DISTR_gaussiancopula_binary_normal_rho & nd)
   : DISTR_gamlss(DISTR_gamlss(nd))
  {
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  }


const DISTR_gaussiancopula_binary_normal_rho & DISTR_gaussiancopula_binary_normal_rho::operator=(
                            const DISTR_gaussiancopula_binary_normal_rho & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gamlss::operator=(DISTR_gamlss(nd));
  response2p = nd.response2p;
  workingresponse2p = nd.workingresponse2p;
  return *this;
  }

void DISTR_gaussiancopula_binary_normal_rho::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

   // *response[1] = *response[2] = *response[3] = first and continuous component of two dimensional reponse
   // *response[0] = workingresponse of binary component
   // *linpred[0] = eta_mu_latent
   // *linpred[1] = eta_sigma
   // *linpred[2] = eta_mu
   // *linpred[3] = eta_rho
   // *linpred[4] = eta_mu_1

   if (*weight[3] == 0)
     *deviance=0;
   else
     {
     double rho = (*linpred[3])/pow(1+pow((*linpred[3]),2),0.5);
     double sigma = exp(*linpred[1]);
     double mu = (*linpred[2]);
     double mu_latent = (*linpred[0]);
     double hilfs1 = 1-pow(rho,2);
     double l;

       l = -log(2*PI)-log(sigma)-0.5*log(hilfs1)-(1/(2*hilfs1))*( pow((((*response[2]))-mu),2)/pow(sigma,2) -
                                                                                2*rho*(((*response[2])-mu)/(sigma))*(((*response[0])-mu_latent))
                                                                                + pow((((*response[0]))-mu_latent),2) );

   // cout << "deviance y1: " << *response[2] << endl;
   // cout << "deviance y2: " << *response[0] << endl;
   // cout << "deviance y2: " << workingresponse2p << endl;
    *deviance = -2*l;
    }

  }

double DISTR_gaussiancopula_binary_normal_rho::get_intercept_start(void)
  {
  return 0; // log(response.mean(0));
  }

void DISTR_gaussiancopula_binary_normal_rho::compute_param_mult(vector<double *>  linpred,double * param)
  {
   double arg = (*linpred[3]);
  *param = arg/pow(1+pow(arg,2),0.5);
  }

void DISTR_gaussiancopula_binary_normal_rho::set_worklin(void)
  {

  DISTR_gamlss::set_worklin();

   response2p = workingresponse2p->getV();
  }



void DISTR_gaussiancopula_binary_normal_rho::modify_worklin(void)
  {

  DISTR_gamlss::modify_worklin();

  if (counter<nrobs)
    {
    response2p++;
    }

  }



double DISTR_gaussiancopula_binary_normal_rho::loglikelihood_weightsone(double * response,
                                                 double * linpred)
  {

  // *worklin[0] = linear predictor of latent equation
  // *worktransformlin[0] = latent;
  // *worklin[1] = linear predictor of sigma equation
  // *worktransformlin[1] = sigma;
  // *worklin[2] = linear predictor of mu equation
  // *worktransformlin[2] = mu;
  if (counter==0)
    {
    set_worklin();
    }
  double rho;

  if (*linpred <= -100)
    rho  = -0.99995;
  else if (*linpred >= 100)
    rho  = 0.99995;
  else
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);

  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;
  double l;


     l = -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response))-(*worktransformlin[2])),2)/pow((*worktransformlin[1]),2) -
                                 2*rho*(((*response)-(*worktransformlin[2]))/(*worktransformlin[1]))*(((*response2p)-(*worktransformlin[0])))
                                +  pow((((*response2p))-(*worktransformlin[0])),2) );


  modify_worklin();

  return l;

  }

void DISTR_gaussiancopula_binary_normal_rho::compute_iwls_wweightschange_weightsone(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  /// *worklin[0] = linear predictor of latent equation
  // *worktransformlin[0] = latent;
  // *worklin[1] = linear predictor of sigma equation
  // *worktransformlin[1] = sigma;
  // *worklin[2] = linear predictor of mu equation
  // *worktransformlin[2] = mu;
  if (counter==0)
    {
    set_worklin();
    }

  double rho;
  double hilfs;

  if (*linpred <= -100) {
    rho  = -0.99995;
    hilfs = 100.05;
  }
  else if (*linpred >= 100) {
    rho  = 0.99995;
    hilfs = 100.05;
  }
  else {
    rho = (*linpred)/pow((1+pow((*linpred),2)),0.5);
    hilfs = pow((1+pow((*linpred),2)),0.5);
  }


  double rho2 = pow(rho,2);
  double oneminusrho2 = 1- rho2;


  double nu = oneminusrho2*(*linpred) - (*linpred)*( pow((((*response))-(*worktransformlin[2])),2)/pow((*worktransformlin[1]),2)
                                                      +  pow((((*response2p))-(*worktransformlin[0])),2) )
                +(hilfs+rho*(*linpred))*( (((*response)-(*worktransformlin[2]))/(*worktransformlin[1]))*(((*response2p)-(*worktransformlin[0]))) );

  // cout << "rho equation y1: " << *response << endl;
  // cout << "rho equation y2: " << *response2p << endl;

  *workingweight = 1-pow(rho2,2);

  *workingresponse = *linpred + nu/(*workingweight);

  if (compute_like)
    {

    like +=  -0.5*log(oneminusrho2) -(1/(2*oneminusrho2))*( pow((((*response))-(*worktransformlin[2])),2)/pow((*worktransformlin[1]),2) -
                                 2*rho*(((*response)-(*worktransformlin[2]))/(*worktransformlin[1]))*(((*response2p)-(*worktransformlin[0])))
                                +  pow((((*response2p))-(*worktransformlin[0])),2) );

    }


  modify_worklin();



  }

void DISTR_gaussiancopula_binary_normal_rho::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = 2 * std::asin((*linpred[predstart_mumult+3]) / (pow(1 + pow((*linpred[predstart_mumult+3]), 2), 0.5))) / PI ;
  }

void DISTR_gaussiancopula_binary_normal_rho::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Link function (rho): fisher z-transformation\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussiancopula_binary_normal_rho::update_end(void)
  {

  // helpmat1 stores rho2

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * pmu = helpmat1.getV();

  unsigned i;
  for (i=0;i<nrobs;i++,pmu++,worklin++)
    {
    *pmu = (*worklin)/pow((1+pow((*worklin),2)),0.5);
    }

  }


} // end: namespace MCMC




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



#if !defined (DISTR_gamlss_copula_INCLUDED)
#define DISTR_gamlss_copula_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"distr_gamlss.h"
#include"distr_gamlss_nadja.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//-------------------------- CLASS: DISTR_copula_basis --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_copula_basis : public DISTR_gamlss
  {

  protected:

  void set_worklin(void);
  void modify_worklin(void);

  public:

  datamatrix response1;
  datamatrix response2;
  double * response1p;
  double * response2p;
  double * weightp;

  double * linpredp;

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_copula_basis(void) : DISTR_gamlss()
    {
    }

   // CONSTRUCTOR

  DISTR_copula_basis(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_copula_basis(const DISTR_copula_basis & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_copula_basis & operator=(const DISTR_copula_basis & nd);

   // DESTRUCTOR

  ~DISTR_copula_basis() {}

  double get_intercept_start(void);

  void compute_param_mult(vector<double *>  linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  //vector<double> derivative(double & F, int & copulapos);

  virtual vector<double> derivative(double & F1, double & F2, double * linpred);

  vector<double> logc(double & F, int & copulapos, const bool & deriv);

  virtual double logc(double & F1, double & F2, double * linpred);

  virtual double condfc(double & x, double & linpred_F, double & y, double & F2, double * linpred);

  double condfc(double & x, double & linpred_F, double & y, int & copulapos);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//-------------------------- CLASS: DISTR_gausscopula --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gausscopula : public DISTR_copula_basis
  {

  protected:

  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_gausscopula(void) : DISTR_copula_basis()
    {
    }

   // CONSTRUCTOR

  DISTR_gausscopula(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gausscopula(const DISTR_gausscopula & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gausscopula & operator=(const DISTR_gausscopula & nd);

   // DESTRUCTOR

  ~DISTR_gausscopula() {}

  double get_intercept_start(void);

  void compute_param_mult(vector<double *>  linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  //vector<double> derivative(double & F, int & copulapos);

  vector<double> derivative(double & F1, double & F2, double * linpred);

  double logc(double & F1, double & F2, double * linpred);

  double condfc(double & x, double & linpred_F, double & y, double & F2, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//-------------------------- CLASS: DISTR_gausscopula2 --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gausscopula2 : public DISTR_copula_basis
  {

  protected:

  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_gausscopula2(void) : DISTR_copula_basis()
    {
    }

   // CONSTRUCTOR

  DISTR_gausscopula2(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gausscopula2(const DISTR_gausscopula2 & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gausscopula2 & operator=(const DISTR_gausscopula2 & nd);

   // DESTRUCTOR

  ~DISTR_gausscopula2() {}

  double get_intercept_start(void);

  void compute_param_mult(vector<double *>  linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  //vector<double> derivative(double & F, int & copulapos);

  vector<double> derivative(double & F1, double & F2, double * linpred);

  double logc(double & F1, double & F2, double * linpred);

  double condfc(double & x, double & linpred_F, double & y, double & F2, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };


//------------------------------------------------------------------------------
//-------------------------- CLASS: DISTR_clayton_copula --------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_clayton_copula : public DISTR_copula_basis
  {

  protected:

  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_clayton_copula(void) : DISTR_copula_basis()
    {
    }

   // CONSTRUCTOR

  DISTR_clayton_copula(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_clayton_copula(const DISTR_clayton_copula & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_clayton_copula & operator=(const DISTR_clayton_copula & nd);

   // DESTRUCTOR

  ~DISTR_clayton_copula() {}

  double get_intercept_start(void);

  void compute_param_mult(vector<double *>  linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  vector<double> derivative(double & F1, double & F2, double * linpred);

  double logc(double & F1, double & F2, double * linpred);

  double condfc(double & x, double & linpred_F, double & y, double & F2, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);
  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };

//------------------------------------------------------------------------------
//-------------------------- CLASS: DISTR_gumbel_copula ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gumbel_copula : public DISTR_copula_basis
  {

  protected:

  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

  DISTR_gumbel_copula(void) : DISTR_copula_basis()
    {
    }

   // CONSTRUCTOR

  DISTR_gumbel_copula(GENERAL_OPTIONS * o, const datamatrix & r,
                       const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

  DISTR_gumbel_copula(const DISTR_gumbel_copula & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_gumbel_copula & operator=(const DISTR_gumbel_copula & nd);

   // DESTRUCTOR

  ~DISTR_gumbel_copula() {}

  double get_intercept_start(void);

  void compute_param_mult(vector<double *>  linpred,double * param);

  double loglikelihood_weightsone(double * response, double * linpred);

  //vector<double> derivative(double & F, int & copulapos);

  vector<double> derivative(double & F1, double & F2, double * linpred);

  double logc(double & F1, double & F2, double * linpred);

  double condfc(double & x, double & linpred_F, double & y, double & F2, double * linpred);

  void compute_iwls_wweightschange_weightsone(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  void update_end(void);

  };



} // end: namespace MCMC


#endif

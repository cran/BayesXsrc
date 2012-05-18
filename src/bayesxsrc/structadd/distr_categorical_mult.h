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



#if !defined (DISTRcategorical_mult_INCLUDED)
#define DISTRcategorical_mult_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"distr.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_multinomprobit ----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_multinomprobit : public DISTR
  {

  protected:

  bool master;                                 // master equation that updates
                                               // utilities yes/no

  vector<DISTR*> othercat;
  unsigned nrcat;                              // total number of categories
                                               // including master
  unsigned nrothercat;                         // number of servant categories                                                

  double maxutility(vector<datamatrix*>
                    responsep,
                    const unsigned & i,
                    const unsigned & cat);


  public:

  datamatrix responsecat;

  void create_responsecat(void);

  void assign_othercat(DISTR* o);


   // DEFAULT CONSTRUCTOR

   DISTR_multinomprobit(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_multinomprobit(GENERAL_OPTIONS * o, bool mast,
                        const datamatrix & r,
                        const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_multinomprobit(const DISTR_multinomprobit & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_multinomprobit & operator=(const DISTR_multinomprobit & nd);

   // DESTRUCTOR

   ~DISTR_multinomprobit() {}

  void compute_mu(const double * linpred,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * deviancesat, double * scale) const;

  double loglikelihood(double * response, double * linpred,
                       double * weight) const;

  double loglikelihood_weightsone(double * response, double * linpred) const;

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void outoptions(void);

  void update(void);

  };


} // end: namespace MCMC


#endif

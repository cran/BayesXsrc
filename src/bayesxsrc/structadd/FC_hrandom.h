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



#if !defined (FChrandomINCLUDED)

#define FChrandomINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC_nonp.h"
#include"design.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_hrandom --------------------------------
//------------------------------------------------------------------------------


enum hrandomtype {additive,mult,multexp};

class __EXPORT_TYPE FC_hrandom  : public FC_nonp
  {

  protected:

  datamatrix pmodemat;
  datamatrix varmat;

  datamatrix beta_prior;                            // sampled re's from prior

  datamatrix likelihoodc,likelihoodn;

  bool simplerandom;
  datamatrix simplerandom_linpred;

  DISTR * likep_RE;

  FC FCrcoeff;

  hrandomtype rtype;

  datamatrix response_o;
  datamatrix linpred_o;

  void set_rcoeff(void);

  void update_linpred_multexp(void);
  void update_response_multexp(void);

  bool posteriormode_multexp(void);
  bool posteriormode_additive(void);

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_hrandom(void);

  // CONSTRUCTOR

  FC_hrandom(MASTER_OBJ * mp, unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp, DISTR * lp_RE,
             const ST::string & t,
           const ST::string & fp, const ST::string & fp2, DESIGN * dp,
           vector<ST::string> & op,vector<ST::string> & vn);

  FC_hrandom(MASTER_OBJ * mp, unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
             const ST::string & t,
           const ST::string & fp, const ST::string & fp2, DESIGN * dp,
           vector<ST::string> & op,vector<ST::string> & vn);


  // COPY CONSTRUCTOR

  FC_hrandom(const FC_hrandom & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_hrandom & operator=(const FC_hrandom & m);

  // DESTRUCTOR

  ~FC_hrandom()
    {
    }


  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  void update_IWLS(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

//  void transform_beta(void);

  void compute_autocorr_all(const ST::string & path, unsigned lag,
                            ofstream & outg) const;

  void get_samples(const ST::string & filename,ofstream & outg) const;

    // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                 const ST::string & pathresults);

  void outgraphs(ofstream & out_stata, ofstream & out_R, ofstream & out_R2BayesX,
                         const ST::string & path);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  //----------------------------------------------------------------------------
  //----------------------- For Cross Validation stuff -------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: sample_for_cv
  // TASK: samples for cv-score

  void sample_for_cv(datamatrix & pred);

  void compute_effect_cv(datamatrix & effect);

  //----------------------------------------------------------------------------
  //-------------------- End: For Cross Validation stuff -----------------------
  //----------------------------------------------------------------------------

  };

// -----------------------------------------------------------------//
// class: FC_hrandom_distributional (for non-normal random effects)
// -----------------------------------------------------------------//

class __EXPORT_TYPE FC_hrandom_distributional  : public FC_hrandom
  {

  protected:

  datamatrix offset_RE;
  datamatrix offsetold_RE;

  // FCs für Hyperparameter

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_hrandom_distributional(void);

  // CONSTRUCTOR

  FC_hrandom_distributional(MASTER_OBJ * mp, unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
             DISTR * lp_RE, const ST::string & t,
           const ST::string & fp, const ST::string & fp2, DESIGN * dp,
           vector<ST::string> & op,vector<ST::string> & vn);

  FC_hrandom_distributional(MASTER_OBJ * mp, unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
             const ST::string & t,
           const ST::string & fp, const ST::string & fp2, DESIGN * dp,
           vector<ST::string> & op,vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_hrandom_distributional(const FC_hrandom_distributional & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_hrandom_distributional & operator=(const FC_hrandom_distributional & m);

  // DESTRUCTOR

  ~FC_hrandom_distributional()
    {
    }

  // FUNCTION: update

  void update(void);

  void compute_autocorr_all(const ST::string & path, unsigned lag,
                            ofstream & outg) const;

  void get_samples(const ST::string & filename,ofstream & outg) const;

    // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                 const ST::string & pathresults);

  void outgraphs(ofstream & out_stata, ofstream & out_R, ofstream & out_R2BayesX,
                         const ST::string & path);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };

} // end: namespace MCMC

#endif



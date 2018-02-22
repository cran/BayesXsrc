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



#if !defined (FCNONPVARIANCEINCLUDED)

#define FCNONPVARIANCEINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC_nonp.h"
#include"design.h"
#include"design_userdefined.h"
#include<cmath>

namespace MCMC
{


/*
vector<ST::string> hyperpriors;
  hyperpriors.push_back("invgamma"); // standard inverse gamma
  hyperpriors.push_back("scaledep"); // scale-dependent as proposed by Klein&Kneib, corresponds to weibull distribution
  hyperpriors.push_back("hcauchy"); // generalised beta prime based on half-cauchy for tau
  hyperpriors.push_back("hnormal");  // gamma prior based on half normal for tau
  hyperpriors.push_back("aunif");      // based on approximation for uniform prior for tau
  hyperprior = stroption("hyperprior",hyperpriors,"invgamma");
*/
enum hyperpr {invgamma,scaledep,hcauchy,hnormal,
              aunif};

enum ssvsupdates {regcoeff, sdev, gibbs};

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_nonp_variance --------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_nonp_variance  : public FC
  {

  protected:

  FC_nonp * FCnonpp;                         // Pointer to corresponding
                                             // FC_nonp object
  DISTR * likep;                             // Pointer to DISTR obejct
  DESIGN * designp;                          // Pointer to design object

  MASTER_OBJ * masterp;
  unsigned equationnr;

  double a_invgamma;
  double b_invgamma_orig;
  double b_invgamma;
  double lambdastart;

  bool lambdaconst;

  double tildea;
  double tildeb;
  bool cauchy;
  bool wei;
  double scaletau2;

  int proposal;

  hyperpr hyperprior;

  ssvsupdates ssvsupdate;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp_variance(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp_variance(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
           const ST::string & t,
           const ST::string & fp,DESIGN * dp,FC_nonp * FCn,
           vector<ST::string> & op,vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_nonp_variance(const FC_nonp_variance & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp_variance & operator=(const FC_nonp_variance & m);

  // DESTRUCTOR

  ~FC_nonp_variance()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                  const ST::string & pathresults);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  // virtual void transform_beta(void);

  };






class __EXPORT_TYPE FC_nonp_variance_varselection  : public FC_nonp_variance
  {

  protected:

  double ssvsvarlimit;

  datamatrix zerovec;
  vector<ST::string> datanames;

  bool singleomega;

  FC FC_psi2;
  FC FC_omega;

  double a_omega;
  double b_omega;

  double v1;
  double v2;
  double v2_orig;

  double r_delta;

  bool gig;

  double r;

  datamatrix X;
  datamatrix Xscaled;
  datamatrix diff;

  double tauold;

  public:

  double omega;
  FC FC_delta;

  //----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp_variance_varselection(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp_variance_varselection(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
          bool so,
           const ST::string & t,
           const ST::string & fp,DESIGN * dp,FC_nonp * FCn,
           vector<ST::string> & op,vector<ST::string> & vn, double vl);

  // COPY CONSTRUCTOR

  FC_nonp_variance_varselection(const FC_nonp_variance_varselection & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp_variance_varselection & operator=(const FC_nonp_variance_varselection & m);

  // DESTRUCTOR

  ~FC_nonp_variance_varselection()
    {
    }

  void add_linpred(datamatrix & l);

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  void update_IWLS(void);
  void update_gaussian(void);



  bool posteriormode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                  const ST::string & pathresults);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void get_samples(const ST::string & filename,ofstream & outg) const;


  // virtual void transform_beta(void);

  };




//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_varselection_omega ---------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FC_varselection_omega  : public FC
  {

  protected:

  double a_omega;
  double b_omega;

  DISTR * likep;                             // Pointer to DISTR obejct

  MASTER_OBJ * masterp;


  public:

  vector<FC_nonp_variance_varselection*> FC_tau2s;


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_varselection_omega(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_varselection_omega(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
           const ST::string & t);

  // COPY CONSTRUCTOR

  FC_varselection_omega(const FC_varselection_omega & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_varselection_omega & operator=(const FC_varselection_omega & m);

  // DESTRUCTOR

  ~FC_varselection_omega()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

//  bool posteriormode(void);
//
//  // FUNCTION: outoptions
//  // TASK: writes estimation options (hyperparameters, etc.) to outputstream
//
  void outoptions(void);
//
//  // FUNCTION: outresults
//  // TASK: writes estimation results to logout or into a file
//
  void outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                  const ST::string & pathresults);
//
//  // FUNCTION: reset
//  // TASK: resets all parameters
//
  void reset(void);
//
//  void read_options(vector<ST::string> & op,vector<ST::string> & vn);
//
  void get_samples(const ST::string & filename,ofstream & outg) const;
//
//
  // virtual void transform_beta(void);

  };

//------------------------------------------------------------------------------
//------------------------ CLASS: FC_tensor_omega ------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_tensor_omega : public FC
  {

  protected:

  datamatrix omegas;
  unsigned nromega;
  unsigned omegaindex;

  DESIGN_userdefined_tensor * dut;
  FC_nonp * fcn;

  datamatrix logdets;

  public:

  // DEFAULT CONSTRUCTOR

  FC_tensor_omega(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // nro  : number of different values for omega
  // fp   : file path for storing sampled parameters

  FC_tensor_omega(DESIGN_userdefined_tensor * d, FC_nonp * f,
           GENERAL_OPTIONS * o,const ST::string & t,
           const unsigned & nro, const ST::string & fp);

  // COPY CONSTRUCTOR

  FC_tensor_omega(const FC_tensor_omega & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_tensor_omega & operator=(const FC_tensor_omega & m);

  // DESTRUCTOR

  ~FC_tensor_omega()
    {
    }

  // FUNCTION: get_samples
  // TASK: stores the sampled parameters in ASCII format

  void get_samples(const ST::string & filename,ofstream & outg) const;

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata, ofstream & out_R, ofstream & out_R2BayesX,
               const ST::string & pathresults);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);
  };

} // end: namespace MCMC

#endif



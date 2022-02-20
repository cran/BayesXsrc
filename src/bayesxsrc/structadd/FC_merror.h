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



#if !defined (FCmerrorINCLUDED)

#define FCmerrorINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"clstring.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"FC_nonp.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_merror ---------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_merror  : public FC
  {

  protected:

  datamatrix xobs;  // observed covariate values with measurement error (n x M)
  datamatrix xmean; // mean of the observed covariare values (n x 1)
  double merror;       // number of replicates (M)

  double minbin;    // smallest value of the binning grid
  double maxbin;    // largest value of the binning grid
  double deltabin;  // length of the binning intervals
  double binning;   // number of binning intervals

  statmatrix<unsigned> indexold;
  statmatrix<unsigned> indexprop;

  statmatrix<int> countmat;

  datamatrix mevar; // measurement error variances (n x 1)
  vector<datamatrix> mecovinv; // vector of inverse covariance matrices
                              // length n, individual elements of dim (M x M)
  datamatrix mesd;  // measurement error standard deviations or trace of the
                    // measurement error covariance matrices (n x 1)


  FC_nonp * FCp;    // pointer to the P-spline full conditional

  FC FC_mu_x;
  double m_mu_x;
  double s_mu_x;
  double s2_mu_x;

  FC FC_tau2_x;
  double a_tau2_x;
  double b_tau2_x;

  double mepropscale;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_merror(void);

  // CONSTRUCTOR

  FC_merror(GENERAL_OPTIONS * o, const ST::string & t,
            const ST::string & fp, vector<ST::string> & op,
            vector<ST::string> & vn, datamatrix & xo, datamatrix & mv,
            datamatrix & xd,
            FC_nonp * fcn);

  // COPY CONSTRUCTOR

  FC_merror(const FC_merror & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_merror & operator=(const FC_merror & m);

  // DESTRUCTOR

  ~FC_merror()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  void outresults(ofstream & out_stata,ofstream & out_R, ofstream & out_R2BayesX,
                 const ST::string & pathresults);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void get_samples(const ST::string & filename, ofstream & outg) const;

  void outoptions(void);

  };

} // end: namespace MCMC

#endif



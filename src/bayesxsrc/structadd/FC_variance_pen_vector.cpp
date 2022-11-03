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



#include"FC_variance_pen_vector.h"



//------------------------------------------------------------------------------
//-- CLASS: FC_variance_pen_vector implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{

using randnumbers::rand_inv_gaussian;
using randnumbers::rand_gamma;


void FC_variance_pen_vector::add_variable(datamatrix & x,vector<ST::string> & op,
                                          vector<ST::string> & vn)
  {

  int f;
  double ta;
  double shrink;
  bool sfix;
  bool adsh;
  double a,b;
  double shw;

  // read options
  f = op[33].strtodouble(ta);
  f = op[29].strtodouble(shrink);
  if (op[30]=="true")
    sfix = true;
  else
    sfix = false;
  if (op[32]=="true")
    adsh = true;
  else
    adsh = false;
  f = op[5].strtodouble(a);
  f = op[6].strtodouble(b);
  f = op[31].strtodouble(shw);


  // options from each term
  tau2.push_back(ta);
  shrinkageweight.push_back(shw);

  // options from first term or all terms. Depends on adaptiveshrinkage
  shrinkagefix.push_back(sfix);
  adaptiveshrinkage.push_back(adsh);

  // shrinkagefix and adaptiveshrinkage are set from first term
  is_fix = shrinkagefix[0];
  is_adaptive = adaptiveshrinkage[0];

  if(is_ridge==1)                        // for Ridge
    shrinkagestart.push_back(1/(2*ta));
  if(is_ridge==0)                        // for Lasso
    shrinkagestart.push_back(shrink);

  a_shrinkagegamma.push_back(a);
  b_shrinkagegamma.push_back(b);

  if(is_adaptive==false)
    {
    shrinkagestart[nrpen] = shrinkagestart[0];
    a_shrinkagegamma[nrpen] = a_shrinkagegamma[0];
    b_shrinkagegamma[nrpen] = b_shrinkagegamma[0];
    }

  // optionsp->out("\n");
  // optionsp->out("OPTIONS FROM TERMS\n");
  // optionsp->out("penalty = " + ST::doubletostring(shelp(0,0)) + "\n");
  // optionsp->out("nrpen = " + ST::doubletostring(nrpen) + "\n");
  // optionsp->out("tau2 = " + ST::doubletostring(tau2[nrpen]) + "\n");
  // optionsp->out("weight = " + ST::doubletostring(shrinkageweight[nrpen]) + "\n");
  // optionsp->out("shrinkage = " + ST::doubletostring(shrinkagestart[nrpen]) + "\n");
  // optionsp->out("a_shrinkage = " + ST::doubletostring(a_shrinkagegamma[nrpen]) + "\n");
  // optionsp->out("b_shrinkage = " + ST::doubletostring(b_shrinkagegamma[nrpen]) + "\n");
  // optionsp->out("\n");

  nrpen++;
  Cp->tau2 = datamatrix(nrpen,1,0);
  Cp->tau2oldinv = datamatrix(nrpen,1,0);
  for(int i=0;i<nrpen;i++)
    {
    Cp->tau2(i,0) = tau2[i];
    }
  }


//______________________________________________________________________________
//
// CONSTRUCTOR with Parameters
//
// mp                : pointer to MASTER_object
// o                 : pointer to GENERAL_options object
// p                 : pointer to FULLCOND_lin_pen object
// d                 : pointer to DISTRIBUTION object
// ti                : reference to title of the full conditional(for example "fixed effects")
// fp                : reference to file path for storing sampled parameters
// isr               : variable,indicates if L2- or L1- penalty is uesd
// fr                : reference to filename for storing results
// shrinkage_start   : Starting value for the shrinkageparameter
// a_shrinkage_gamma : Hyperparameter for gammaprior of shrinkageparameter
// b_shrinkagegamma  : Hyperparameter for gammaprior of shrinkageparameter
// shrinkage_fix     : Should the shrinkageparameter be fixed at the value "shrinkage_start"
// ct                : blocks of regression coefficients
// c                 : reference to responsecategory (important in the case of multivariate response)
//______________________________________________________________________________
//==============================================================================
// KOMMENTAR Susanne:
//-------------------
// Derzeit ist die scale-dependent Variante implementiert, d.h. die
// Shrinkageprioris haengen in Gaussfall vom Skalenparameter sigma2 ab.
// Dies hat zur Folge, das zum update von sigma2 der Wert sum(beta^2/tau^2),
// der in lassosum und ridgesum berechnet wird, benoetigt wird.
// Fuer den Fall das kein adaptives shrinkage gewaehlt wird, werden die
// ersten Eintraege aus den Vektoroptionen fuer das Verfahren verwendet.
// Einentsprechender Hinweis sollte in das Handbuch.
//==============================================================================

FC_variance_pen_vector::FC_variance_pen_vector(MASTER_OBJ * mp,
                        GENERAL_OPTIONS * o, FC_linear_pen * p,
                        DISTR * d,const ST::string & ti,
                        const ST::string & fp, shrinktype sh)
                        :FC(o,ti,1,1,fp)
  {

  masterp = mp;

  shtype = sh;
  if (shtype==ridge)
    is_ridge = true;
  else
    is_ridge = false;

  update_sigma2 = true;

  Cp = p;
  shelp = datamatrix(3,1,0);
  if (is_ridge)
    shelp(0,0) = 2;
  else
    shelp(0,0) = 1;

  distrp = d;

  priorassumptions.push_back("\\\\");

  // vector<ST::string> vnames = FC_shrinkage.get_datanames();
  //Initialisieren der Betamatrizen fuer die Varianzan + Uebergabe der Startwerte
  //is_fix = shrinkagefix[0];
  //is_adaptive = adaptiveshrinkage[0];

  nrpen = 0;
  pensum = 0;

  }
//______________________________________________________________________________
//
// COPY CONSTRUCTOR
//______________________________________________________________________________

FC_variance_pen_vector::FC_variance_pen_vector(const FC_variance_pen_vector & t)
    : FC(FC(t))
  {
  masterp = t.masterp;
  tau2 = t.tau2;
  shelp = t.shelp;
  update_sigma2 = t.update_sigma2;
  Cp = t.Cp;
  distrp = t.distrp;
  FC_shrinkage = t.FC_shrinkage;
  shrinkagefix = t.shrinkagefix;
  adaptiveshrinkage = t.adaptiveshrinkage;
  a_shrinkagegamma = t.a_shrinkagegamma;
  b_shrinkagegamma = t.b_shrinkagegamma;
  shrinkagestart = t.shrinkagestart;
  shrinkageweight = t.shrinkageweight;
  nrpen = t.nrpen;
  pensum = t.pensum;
  is_ridge = t.is_ridge;
  is_fix = t.is_fix;
  is_adaptive = t.is_adaptive;
  }


//______________________________________________________________________________
//
// OVERLOADED ASSIGNMENT OPERATOR
//______________________________________________________________________________

const FC_variance_pen_vector & FC_variance_pen_vector::operator=(
                                        const FC_variance_pen_vector & t)
  {
  if (this == &t)
    return *this;
  FC::operator=(FC(t));
  masterp = t.masterp;
  tau2 = t.tau2;
  shelp = t.shelp;
  update_sigma2 = t.update_sigma2;
  Cp = t.Cp;
  distrp = t.distrp;
  FC_shrinkage = t.FC_shrinkage;
  shrinkagefix = t.shrinkagefix;
  adaptiveshrinkage = t.adaptiveshrinkage;
  a_shrinkagegamma = t.a_shrinkagegamma;
  b_shrinkagegamma = t.b_shrinkagegamma;
  shrinkagestart = t.shrinkagestart;
  shrinkageweight = t.shrinkageweight;
  nrpen = t.nrpen;
  pensum = t.pensum;
  is_ridge = t.is_ridge;
  is_fix = t.is_fix;
  is_adaptive = t.is_adaptive;
  return *this;
  }


bool FC_variance_pen_vector::posteriormode(void)
  {

  shelp(1,0) = double(tau2.size());

  datamatrix tau2inv(tau2.size(),1,0);
  unsigned i;
  for (i=0;i<tau2inv.rows();i++)
    tau2inv(i,0) = 1/tau2[i];

  setbeta(tau2inv);

  return true;
  }


//______________________________________________________________________________
//
// FUNCTION: update
// TASK: - stores sampled parameters in file 'samplepath'
//         storing order: first row, second row, ...
//______________________________________________________________________________

void FC_variance_pen_vector::update(void)
  {
  acceptance++;
  int i;
  double helpshrinkage;


  // reset Variables for summs
//  int nrpen = beta.rows();
  double rand_invgaussian = 0;
  double sumvariances = 0;
  double sumregcoeff = 0;

  pensum = 0;

  if (optionsp->nriter == 1)
    {

    FC_shrinkage = FC(optionsp,"",shrinkagefix.size(),1,samplepath + ".shrinkage");
//    Cp->tau2 = datamatrix(nrpen,1,0);
  // get current value of first shrinkagearameter
    double * shrinkagep = FC_shrinkage.beta.getV();
    for(i=0;i<nrpen;i++,shrinkagep++)
      {
      *shrinkagep = shrinkagestart[i];
      }
    }

  // get current value of first shrinkagearameter
  double * shrinkagep = FC_shrinkage.beta.getV();

  // get current value of first regressionparameter
  double * workbeta = Cp->beta.getV();

  // getcurrent value of sqrt(scale) parameter
  double sigma = sqrt(distrp->get_scale());

//TEMP:BEGIN--------------------------------------------------------------------
  // ofstream outpenreg("c:/bayesx/test/outpenreg.txt", ios::out|ios::app);
  // outpenreg << "beta ";
  // for(unsigned i=0; i<nrpen; i++,workbeta++) {outpenreg << *workbeta  << " ";}
  // outpenreg << "\n";
  // outpenreg << "tau2 ";
  // for(unsigned i=0; i<nrpen; i++) {outpenreg << beta(i,0)  << " ";}
  // outpenreg << "\n";
  // outpenreg << "lambda ";
  // for(unsigned i=0; i<nrpen; i++,shrinkagep++) {outpenreg << *shrinkagep  << " ";}
  // outpenreg << "\n";
  // outpenreg << "sigma "<< sigma <<"\n";

  // shrinkagep = FC_shrinkage.beta.getV();
  // workbeta = Cp->beta.getV();
  // sigma = sqrt(distrp->get_scale());
//TEMP:END----------------------------------------------------------------------



  // Gibbs-Update of varianceparameters tau^2
  //-----------------------------------------
  if (is_ridge == 0 && is_adaptive == false)   // Lasso L1-Penalty
   {
   for(i=0; i<nrpen; i++, workbeta++, shrinkagep++)
    {
    if (*workbeta>0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(sigma*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
      beta(i,0) = 1.0/rand_invgaussian;
      }
    if (*workbeta<0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(-1.0*sigma*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
      beta(i,0) = 1.0/(rand_invgaussian);
      }
    if (*workbeta==0 || *shrinkagep<=0)
      {
      beta(i,0) = 1E-6;
      }
    pensum += ((*workbeta)*(*workbeta))/beta(i,0);                          // sum(beta^2/tau^2)
    sumvariances = sumvariances + beta(i,0);  // sum(tau^2/weights^2) of current variances
   }
  }

  if (is_ridge == 0 && is_adaptive == true)   // Lasso L1-Penalty adaptive
   {
   for(i=0; i<nrpen; i++, workbeta++, shrinkagep++)
    {
    if (*workbeta>0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(sigma*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
      beta(i,0) = 1.0/rand_invgaussian;
      }
    if (*workbeta<0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(-1.0*sigma*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
      beta(i,0) = 1.0/(rand_invgaussian);
      }
    if (*workbeta==0 || *shrinkagep<=0)
      {
      beta(i,0) = 1E-6;
      }
    pensum +=  ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
    //sumvariances = sumvariances + beta(i,0);                    // sum(tau^2) of current variances
   }
  }

  if (is_ridge == 1 && is_adaptive == false)   // Ridge L2-penalty
    {
    for(i=0; i<nrpen; i++, workbeta++, shrinkagep++)
      {
      beta(i,0) = 1/(2*(*shrinkagep));
      pensum +=  ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
      sumregcoeff = sumregcoeff + (*workbeta)*(*workbeta);
      }
    }

  if (is_ridge == 1 && is_adaptive == true)   // Ridge L2-penalty adaptive
    {
    for(i=0; i<nrpen; i++, workbeta++, shrinkagep++)
      {
      beta(i,0) = 1/(2*(*shrinkagep));
      pensum += ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
      //sumregcoeff = sumregcoeff + (*workbeta)*(*workbeta);
      }
    }

  // Gibbs-Update of Shrinkageparameter with Gammadistribution
  //----------------------------------------------------------
  if(is_fix==false)
    {
    shrinkagep = FC_shrinkage.beta.getV();
    workbeta = Cp->beta.getV();
    if(is_ridge == 0 && is_adaptive == false)            // Lasso L1-penalty
      {
      helpshrinkage = sqrt(rand_gamma(nrpen + a_shrinkagegamma[0], b_shrinkagegamma[0] + 0.5*sumvariances));
      for(i=0; i<nrpen; i++, shrinkagep++)
        {
        *shrinkagep = helpshrinkage;
        }
      }
    if(is_ridge == 0 && is_adaptive == true)            // Lasso L1-penalty adaptive
      {
      for(i=0; i<nrpen; i++, shrinkagep++)
        {
        *shrinkagep = sqrt(rand_gamma(1 + a_shrinkagegamma[i], b_shrinkagegamma[i] + 0.5*beta(i,0)));
        }
      }


    if(is_ridge == 1 && is_adaptive == false)            // Ridge L2-penalty
      {
      helpshrinkage = rand_gamma(0.5*nrpen + a_shrinkagegamma[0], b_shrinkagegamma[0] + sumregcoeff/(sigma*sigma));
      for(i=0; i<nrpen; i++, shrinkagep++)
        {
        *shrinkagep = helpshrinkage;
        }
      }
    if(is_ridge == 1 && is_adaptive == true)            // Ridge L2-penalty adaptive
      {
      for(i=0; i<nrpen; i++, shrinkagep++, workbeta++)
        {
        *shrinkagep = rand_gamma(0.5*1 + a_shrinkagegamma[i], b_shrinkagegamma[i] + ((*workbeta)*(*workbeta))/(sigma*sigma));
        }
      }
    }


  shelp(2,0) = pensum;
  distrp->update_scale_hyperparameters(shelp);

  Cp->tau2.assign(beta);

  // Update Shrinkageparameter
  FC_shrinkage.acceptance++;
  FC_shrinkage.update();

//  acceptance++;
  FC::update();


  }


//______________________________________________________________________________
//
// FUNCTION: get_samples
// TASK: - write samples to files
//______________________________________________________________________________

void FC_variance_pen_vector::get_samples(const ST::string & filename,ofstream & outg) const
  {

  FC::get_samples(filename,outg);

  ST::string pathhelp = filename.substr(0,filename.length()-14)+"shrinkage_sample.raw";

  FC_shrinkage.get_samples(pathhelp,outg);

  }

//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results for varianceparameters to output window and files
//______________________________________________________________________________

void FC_variance_pen_vector::outresults(ofstream & out_stata, ofstream & out_R, ofstream & out_R2BayesX,
                             const ST::string & pathresults)
  {

//  int i;
  vector<ST::string> vnames;
  vnames = Cp->datanames;

  FC::outresults(out_stata,out_R,out_R2BayesX,"");
  FC::outresults_help(out_stata,out_R,pathresults,vnames);

  optionsp->out("\n");
  optionsp->out("    Results for variances are also stored in file\n");
  optionsp->out("    " + pathresults + "\n");
  optionsp->out("\n");


  // outresults shrinkage
  ST::string shrinkage_pathresults = pathresults.substr(0,pathresults.length()-7) + "shrinkage.res";
  outresults_shrinkage(shrinkage_pathresults);

  }

//______________________________________________________________________________
//
// FUNCTION: outresults_shrinkage
// TASK: - write results for shrinkageparameter to output window and files
//______________________________________________________________________________

//void FC_variance_pen_vector::outresults_shrinkage(void)
void FC_variance_pen_vector::outresults_shrinkage(const ST::string & pathresults)
  {

  if(is_fix==false)
    {
    optionsp->out("\n");
    optionsp->out("  MAIN_REGRESSION: linear effects (shrinkage) \n");
    optionsp->out("\n");
    ST::string shrinkage_pathresults = pathresults;
    vector<ST::string> vnames;
    vnames = Cp->datanames;

    ofstream out1;
    ofstream out2;
    ofstream out3;

    FC_shrinkage.outresults(out1,out2,out3,"");
    FC_shrinkage.outresults_help(out1,out2,shrinkage_pathresults,vnames);

    // unsigned nr;
    // if(is_adaptive==true)
      // {
      // FC_shrinkage.outresults_help(out1,out2,shrinkage_pathresults,vnames);
      // }

    // if(is_adaptive==false)
      // {
      // FC_shrinkage.outresults_help(out1,out2,shrinkage_pathresults,vnames);
      // }

    optionsp->out("\n");
    optionsp->out("    Results for shrinkage parameters are also stored in file\n");
    optionsp->out("    " + shrinkage_pathresults + "\n");
    optionsp->out("\n");
    }
  }


//______________________________________________________________________________
//
// FUNCTION: outoptions
// TASK: - write options to output window
//______________________________________________________________________________

void FC_variance_pen_vector::outoptions(void)
  {

//TEMP:BEGIN--------------------------------------------------------------------
// KOMMENTAR Susanne:  FC_shrinkage kann dann hier wieder rausgenommen weden
//   FC_shrinkage = FC(optionsp,"",shrinkagefix.size(),1,samplepath + ".shrinkage");
//TEMP:END----------------------------------------------------------------------

  int i;
//  int nrpen = beta.rows();
  vector<ST::string> vnames;   // = FC_shrinkage.datanames.getV();
  vnames = Cp->datanames;

  if(is_ridge == 0)
    {
    optionsp->out("  LINEAR EFFECTS WITH LASSO PENALTY: \n");
    }
  if(is_ridge == 1)
    {
    optionsp->out("  LINEAR EFFECTS WITH RIDGE PENALTY: \n");
    }

  if(is_adaptive==false)
    {
    optionsp->out("  Hyperparameter a for shrinkage: " +
                     ST::doubletostring(a_shrinkagegamma[0]) + "\n" );
    optionsp->out("  Hyperparameter b for shrinkage: " +
                     ST::doubletostring(b_shrinkagegamma[0]) + "\n" );
    if(is_fix==true)
      {
      optionsp->out("  Shrinkage is fixed at value = " +
                       ST::doubletostring(shrinkagestart[0]) + "\n" );
      }
    optionsp->out("\n");
    }

  if(is_adaptive==true)
    {
    for(i=0; i<nrpen; i++)
      {
      optionsp->out("  Hyperparameter a for shrinkage of " + vnames[i] + ": " +
                       ST::doubletostring(a_shrinkagegamma[i]) + "\n" );
      optionsp->out("  Hyperparameter b for shrinkage of " + vnames[i] + ": " +
                       ST::doubletostring(b_shrinkagegamma[i]) + "\n" );
      if(is_fix==true)
        {
        optionsp->out("  Shrinkage of " + vnames[i] + " is fixed at value = " +
                         ST::doubletostring(shrinkagestart[i]) + "\n" );
        }
      }
    optionsp->out("\n");

    }

  }




//--------------------------------------------------------------------------------------------------------------
//-------------------------------------- FC_variance_pen_vector_ssvs -------------------------------------------
//--------------------------------------------------------------------------------------------------------------



void FC_variance_pen_vector_ssvs::add_variable(datamatrix & x,vector<ST::string> & op,
                                               vector<ST::string> & vn)
  {

  int f;
  double a,b,rhelp;
  bool cpr;

  f = op[53].strtodouble(a);
  f = op[54].strtodouble(b);
  f = op[41].strtodouble(rhelp);

  if(op[57]=="true")
     NBPSS = true;
  else
     NBPSS = false;

  atau2.push_back(a);
  btau2_orig.push_back(b);
  btau2.push_back(b);
//  btau2.push_back(masterp->level1_likep[equationnr]->trmult*b);
  r.push_back(rhelp);

  if (op[61] == "true")
    {
    cpr = true;
    }
  else
    cpr = false;

  cprior.push_back(cpr);
  nrpen++;

  delta = FC(optionsp,"",nrpen,2,"");
  delta.setbeta(nrpen,1,0);
  delta.setbeta(nrpen,2,0.5);
  setbeta(nrpen,1,0.1);

  FC_psi2 = FC(optionsp,"",nrpen,1,"");
  FC_psi2.setbeta(nrpen,1,0.1);

  FC_delta = FC(optionsp,"",nrpen,2,"");
  FC_delta.setbeta(nrpen,1,0);
  FC_delta.setbeta(nrpen,2,0.5);

  Cp->tau2 = datamatrix(nrpen,1,0);
  Cp->tau2oldinv = datamatrix(nrpen,1,0);
  for(int i=0;i<nrpen;i++)
    {
    Cp->tau2(i,0) = beta(i,0);
    }
  }



FC_variance_pen_vector_ssvs::FC_variance_pen_vector_ssvs(MASTER_OBJ * mp,
                        unsigned & enr,
                        GENERAL_OPTIONS * o, FC_linear_pen * p,
                        DISTR * d,const ST::string & ti,
                        const ST::string & fp)
                        :FC(o,ti,1,1,fp)
  {

  masterp = mp;
  equationnr = enr;

  Cp = p;

  distrp = d;

  nrpen = 0;

  theta = FC(optionsp,"",1,1,"");
  theta.setbeta(1,1,0.5);

  FC_omega = FC(optionsp,"",1,1,"");
  FC_omega.setbeta(1,1,0.5);

  atheta=1;
  btheta=1;

  }

FC_variance_pen_vector_ssvs::FC_variance_pen_vector_ssvs(const FC_variance_pen_vector_ssvs & t)
    : FC(FC(t))
  {
  masterp = t.masterp;
  equationnr = t.equationnr;

  cprior = t.cprior;
  nrpen = t.nrpen;

  distrp = t.distrp;
  Cp = t.Cp;

  delta = t.delta;
  atau2 = t.atau2;
  btau2 = t.btau2;
  btau2_orig = t.btau2_orig;
  theta = t.theta;
  atheta = t.atheta;
  btheta = t.btheta;

  r = t.r;

  FC_omega = t.FC_omega;
  FC_psi2 =  t.FC_psi2;
  FC_delta = t.FC_delta;

  NBPSS = t.NBPSS;
  }


const FC_variance_pen_vector_ssvs & FC_variance_pen_vector_ssvs::operator=(const FC_variance_pen_vector_ssvs & t)
  {
  if (this == &t)
    return *this;
  FC::operator=(FC(t));
  masterp = t.masterp;
  equationnr = t.equationnr;

  cprior = t.cprior;
  nrpen = t.nrpen;

  distrp = t.distrp;
  Cp = t.Cp;

  delta = t.delta;
  atau2 = t.atau2;
  btau2 = t.btau2;
  btau2_orig = t.btau2_orig;
  theta = t.theta;
  atheta = t.atheta;
  btheta = t.btheta;

  r = t.r;

  FC_omega = t.FC_omega;
  FC_psi2 =  t.FC_psi2;
  FC_delta = t.FC_delta;

  NBPSS = t.NBPSS;

  return *this;
  }


bool FC_variance_pen_vector_ssvs::posteriormode(void)
  {

  /*
  shelp(1,0) = double(tau2.size());

  datamatrix tau2inv(tau2.size(),1,0);
  unsigned i;
  for (i=0;i<tau2inv.rows();i++)
    tau2inv(i,0) = 1/tau2[i];

  setbeta(tau2inv);
  */

  return true;
  }





void FC_variance_pen_vector_ssvs::update(void)
  {
  unsigned j;
  for(j=0; j<btau2.size(); j++)
    btau2[j] = masterp->level1_likep[equationnr]->trmult*btau2_orig[j];

//  cout << masterp->level1_likep[equationnr]->trmult << endl;
//  cout << distrp->family << endl;
//  cout << endl;

  acceptance++;
  if(NBPSS)
    {
    double p, q, c, r_delta;
    double anew, bnew;
    double sumdelta=0.0;
    for (j=0; j<nrpen; j++)
      {
      // update psi2_j
      if(FC_delta.beta(j,0) > 0)
        r_delta = 1.0;
      else
        r_delta = r[j];

      anew = atau2[j] + 0.5;
      bnew = btau2[j]+0.5*beta(j,0)/r_delta;
      FC_psi2.beta(j,0) = rand_invgamma(anew, bnew);

//cout << "r_delta: " << r_delta << endl;
//cout << "beta(j,0): " << beta(j,0) << endl;
//cout << "atau2[j]: " << atau2[j] << endl;
//cout << "btau2[j]: " << btau2[j] << endl;
//cout << "anew: " << anew << endl;
//cout << "bnew: " << bnew << endl;
//cout << "FC_psi2.beta(j,0): " << FC_psi2.beta(j,0) << endl;

      // update delta_j

      double u = uniform();
      double L = 1/sqrt(r[j])*exp(- beta(j,0)/(2*FC_psi2.beta(j,0))*(1/r[j]-1));
      double pr1 = 1/(1+ ((1-FC_omega.beta(0,0))/FC_omega.beta(0,0))*L);

      if (u <= pr1)
        {
        FC_delta.beta(j,0) = 1;
        r_delta = 1;
        }
      else
        {
        FC_delta.beta(j,0) = 0;
        r_delta = r[j];
        }

      FC_delta.beta(j,1) = pr1;
      sumdelta += FC_delta.beta(j,0);

      // update tau^2_j

      p = 0.0;
      q = 1.0/(r_delta * FC_psi2.beta(j,0));
      c = pow(Cp->beta(j,0),2);

//cout << "r_delta: " << r_delta << endl;
//cout << "p: " << p << endl;
//cout << "q: " << q << endl;
//cout << "c: " << c << endl;

      double tau2 = randnumbers::GIG2(p, q, c);
      beta(j,0) = tau2;

//cout << "tau2: " << tau2 << endl;
//cout << endl;

      Cp->tau2(j,0) = beta(j,0);
      if(cprior[j])
        Cp->tau2(j,0) = beta(j,0)*distrp->get_scale();
      }

    // update omega

    FC_omega.beta(0,0) = randnumbers::rand_beta(atheta+sumdelta,btheta+nrpen-sumdelta);

    FC_omega.update();
    FC_psi2.update();
    FC_delta.update();

    }
  else
  {
  unsigned j;
  double anew_tau2;
  double bnew_tau2;
  double beta2;
//  double thetanew;

  double sumdelta=0;


  for (j=0;j<nrpen;j++)
    {

    // update tau^2_j

    anew_tau2 = atau2[j] + 0.5;
    beta2 = pow(Cp->beta(j,0),2);

    if (delta.beta(j,0) > 0)
      {
      bnew_tau2 = btau2[j]+0.5*beta2;
      }
    else
      {
      bnew_tau2 = btau2[j]+0.5*beta2/r[j];
      }
    double tau2help = rand_invgamma(anew_tau2,bnew_tau2);

    if (delta.beta(j,0) > 0)
      {
      beta(j,0) = tau2help;
      }
    else
      {
      beta(j,0) = r[j]*tau2help;
      }

    Cp->tau2(j,0) = beta(j,0);
    if(cprior[j])
      {
  //    cout << "j: " << j << endl;
      Cp->tau2(j,0) = beta(j,0)*distrp->get_scale();
      }

 //   cout << "atau2: " << atau2[j] << endl;
 //   cout << "btau2: " << btau2[j] << endl;
    //cout << "sigma2: " << distrp->get_scale() << endl;

    // update delta_j

//    btau2_pow_atau2 = pow(btau2[j],atau2[j]);
//    nu0btau2_pow_atau2 = pow(nu0*btau2[j],atau2[j]);


//    helpIG1 = btau2_pow_atau2*exp(-btau2[j]/beta(j,0));
//    helpIG2 = nu0btau2_pow_atau2*exp(-nu0*btau2[j]/beta(j,0));

    double u = uniform();
    double L = 1/sqrt(r[j])*exp(- beta2/(2*tau2help)*(1/r[j]-1));
    double thetanew = 1/(1+ ((1-theta.beta(0,0))/theta.beta(0,0))*L);
  //  thetanew = (theta.beta(0,0) * helpIG1)/(theta.beta(0,0) * helpIG1 + (1-theta.beta(0,0))*helpIG2);

//    cout << "theta: " << thetanew << endl;
    delta.beta(j,1) = thetanew;
    if (u <= thetanew)
      delta.beta(j,0) = 1;
    else
      delta.beta(j,0) = 0;


    sumdelta += delta.beta(j,0);

    }

  delta.update();

  // update theta

  theta.beta(0,0) = randnumbers::rand_beta(atheta+sumdelta,btheta+nrpen-sumdelta);

  /*cout << "atheta: " << atheta << endl;
  cout << "btheta: " << btheta << endl;
  cout << "nrpen: " << nrpen << endl;
  cout << "theta: " << theta.beta(0,0) << endl;
*/

  theta.update();
  }


  FC::update();
}


void FC_variance_pen_vector_ssvs::get_samples(const ST::string & filename,ofstream & outg) const
  {
  FC::get_samples(filename, outg);

  if(NBPSS)
     {
     ST::string filename_delta = filename.substr(0,filename.length()-15) + "_delta_sample.raw";
     FC_delta.get_samples(filename_delta,outg);

     ST::string filename_omega = filename.substr(0,filename.length()-15) + "_omega_sample.raw";
     FC_omega.get_samples(filename_omega,outg);
     }
  else
     {
     ST::string filename_delta = filename.substr(0,filename.length()-15) + "_delta_sample.raw";
     delta.get_samples(filename_delta,outg);

     ST::string filename_theta = filename.substr(0,filename.length()-15) + "_omega_sample.raw";
     theta.get_samples(filename_theta,outg);
     }
  }


void FC_variance_pen_vector_ssvs::outresults(ofstream & out_stata, ofstream & out_R, ofstream & out_R2BayesX,
                                             const ST::string & pathresults)
  {

   optionsp->out("\n");
   optionsp->out("    Variance parameters:\n");

   FC::outresults(out_stata,out_R,out_R2BayesX,"");
   FC::outresults_help(out_stata,out_R,pathresults,Cp->datanames);


   optionsp->out("\n");
   optionsp->out("    Results for linear effects variance parameters are also stored in file\n");
   optionsp->out("    " + pathresults + "\n");
   optionsp->out("\n");

   optionsp->out("\n");
   optionsp->out("    Inclusion probabilities:\n");

   ST::string pathresults_delta = pathresults.substr(0,pathresults.length()-7)+"delta.res";
   ST::string pathresults_delta_prob = pathresults.substr(0,pathresults.length()-7)+"delta_prob.res";

   if(NBPSS)
     {
     FC_delta.outresults(out_stata,out_R,out_R2BayesX,"");
     FC_delta.outresults_help(out_stata,out_R,pathresults_delta,Cp->datanames);
     }
   else
     {
     delta.outresults(out_stata,out_R,out_R2BayesX,"");
     delta.outresults_help(out_stata,out_R,pathresults_delta,Cp->datanames);
     }

   optionsp->out("\n");
   optionsp->out("    Results for delta parameters are also stored in file\n");
   optionsp->out("    " + pathresults_delta + "\n");
   optionsp->out("\n");
   optionsp->out("\n");

   optionsp->out("    Rao-Blackwellised inclusion probabilities:\n");

   if(NBPSS)
     {
     FC_delta.outresults_help(out_stata,out_R,pathresults_delta_prob,Cp->datanames,1);
     }
   else
     {
     delta.outresults_help(out_stata,out_R,pathresults_delta_prob,Cp->datanames,1);
     }
   optionsp->out("\n");
   optionsp->out("    Results for Rao-Blackwellised inclusion probabilities are also stored in file\n");
   optionsp->out("    " + pathresults_delta + "\n");
   optionsp->out("\n");
   optionsp->out("\n");


/*   unsigned j;
   ofstream ou(pathresults_delta_prob.strtochar());

   ou << "pmean_prob" << endl;
   for(j=0;j<nrpen;j++)
     {
    ou  << " " << delta.betamean(j,1) << endl;
     }*/

   ST::string pathresults_theta = pathresults.substr(0,pathresults.length()-7)+"omega.res";

   optionsp->out("    Inclusion probability parameter omega:\n");
   optionsp->out("\n");

   if(NBPSS)
     {
     FC_omega.outresults(out_stata,out_R,out_R2BayesX,"");
     FC_omega.outresults_singleparam(out_stata,out_R,pathresults_theta);
     }
   else
     {
     theta.outresults(out_stata,out_R,out_R2BayesX,"");
     theta.outresults_singleparam(out_stata,out_R,pathresults_theta);
     }

   optionsp->out("    Results for linear effects omega parameters are also stored in file\n");
   optionsp->out("    " + pathresults_theta + "\n");
   optionsp->out("\n");

  }


void FC_variance_pen_vector_ssvs::outoptions(void)
  {
  FC::outoptions();

  int maxvarnamelength = 0;
  int len;

  for(unsigned i=0;i<Cp->datanames.size();i++)
    {
    len = Cp->datanames[i].length();

    if (len > maxvarnamelength)
      maxvarnamelength = len;
    }

  unsigned nsp, nsp2, nsp3;
  if (maxvarnamelength  > 10)
    nsp = 4 + maxvarnamelength - 8;
  else
    nsp = 4;
  ST::string l(' ',nsp);

  optionsp->out("  Hyperparameters for SSVS priors of linear effects:\n\n" );
  optionsp->out("    Variable" + l +
                "v1       " +
                "v2       " +
                "r        \n");

  for(unsigned i=0; i< Cp->datanames.size(); i++)
    {
    if (maxvarnamelength  > 10)
      nsp = 4 + maxvarnamelength - Cp->datanames[i].length();
    else
      nsp = 12-Cp->datanames[i].length();
    ST::string ls(' ',nsp);

    nsp2 = 9 - ST::doubletostring(atau2[i],3).length();
    ST::string ls2(' ',nsp2);

    nsp3 = 9 - ST::doubletostring(btau2[i],3).length();
    ST::string ls3(' ',nsp3);

    optionsp->out("    " + Cp->datanames[i] + ls +
                  ST::doubletostring(atau2[i],3) + ls2 +
                  ST::doubletostring(btau2[i],3) + ls3 +
                  ST::doubletostring(r[i],3) + "\n");
    }
  optionsp->out("\n");
  }



} // end: namespace MCMC






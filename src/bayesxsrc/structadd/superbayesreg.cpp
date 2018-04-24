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


#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#include<statwin_haupt.h>

#endif

#include"superbayesreg.h"
#include <mapobject.h>

// Vorschlag:
//#include<typeinfo.h>

#include<stddef.h>
using std::ifstream;


//------------------------------------------------------------------------------
//------------- CLASS superbayesreg: implementation of member functions --------
//------------------------------------------------------------------------------

bool superbayesreg::find_binomial(DISTR* & b)
  {
  bool found = false;

  if (distr_binomials.size()==1)
    {
    found = true;
    b = &distr_binomials[0];
    }
  else if (distr_binomialprobits.size()==1)
    {
    found = true;
    b = &distr_binomialprobits[0];
    }
  else if (distr_logit_fruehwirths.size()==1)
    {
    found = true;
    b = &distr_logit_fruehwirths[0];
    }
  else if (distr_cloglogs.size()==1)
    {
    found = true;
    b = &distr_cloglogs[0];
    }
  return found;
  }



bool superbayesreg::find_continuous_singleparam(DISTR* & m)
  {
  bool found = false;

  if (distr_gaussians.size()==1)
    {
    found = true;
    m = &distr_gaussians[0];
    }
  else if (distr_loggaussians.size()==1)
    {
    found = true;
    m = &distr_loggaussians[0];
    }

  return found;
  }


bool superbayesreg::find_continuous_multparam(vector<DISTR*> & m)
  {
  bool found = false;

  if (distr_lognormal_mus.size()==1)
    {
    found = true;
    m.push_back(&distr_lognormal_sigma2s[0]);
    m.push_back(&distr_lognormal_mus[0]);
    }
  if (distr_lognormal2_mus.size()==1)
    {
    found = true;
    m.push_back(&distr_lognormal2_sigmas[0]);
    m.push_back(&distr_lognormal2_mus[0]);
    }
  if (distr_invgaussian_mus.size()==1)
    {
    found = true;
    m.push_back(&distr_invgaussian_sigma2s[0]);
    m.push_back(&distr_invgaussian_mus[0]);
    }
  if (distr_gamma_mus.size()==1)
    {
    found = true;
    m.push_back(&distr_gamma_sigmas[0]);
    m.push_back(&distr_gamma_mus[0]);
    }
  else if (distr_hetgaussians.size()==1)
    {
    found = true;
    m.push_back(&distr_vargaussians[0]);
    m.push_back(&distr_hetgaussians[0]);
    }

  return found;
  }

  vector<DISTR_gamma_mu> distr_gamma_mus;
  vector<DISTR_gamma_sigma> distr_gamma_sigmas;


void superbayesreg::make_paths(ST::string & pathnonp,
                          ST::string & pathres,
                          ST::string & title,vector<ST::string> vn,
                          ST::string  endingraw,ST::string endingres,
                          ST::string  endingtitle)
  {

  unsigned modnr = equations.size()-1;

  ST::string h = equations[modnr].paths;

  ST::string varname1 = vn[0];
  ST::string varname2;

  if (vn.size() >=2)
    varname2 = vn[1];
  else
    varname2 = "";

  if (varname2=="")
    {
#if defined(__BUILDING_LINUX)
    pathnonp = defaultpath + "/temp/" + name + "_" + h + "_f_"
                                  + varname1 + endingraw;
#else
    pathnonp = defaultpath + "\\temp\\" + name + "_" + h + "_f_"
                                  + varname1 + endingraw;
#endif

    pathres = outfile.getvalue() + "_" + h + "_" +
                     endingres + "_" + varname1 + ".res";

    title = h + ": " + endingtitle + varname1;

    }
  else
    {
#if defined(__BUILDING_LINUX)
    pathnonp = defaultpath + "/temp/" + name + "_" + h + "_" + varname2
                 + "_f_"  + varname1 + endingraw;
#else
    pathnonp = defaultpath + "\\temp\\" + name + "_" + h + "_" + varname2
                 + "_f_"  + varname1 + endingraw;
#endif

    pathres = outfile.getvalue() + "_" + h + "_" +
                     endingres + "_" + varname1 + "_" + varname2 + ".res";


    title = h + ": " + endingtitle + varname1 + " and " + varname2;

    }

  }

// forward declaration
void hregressrun(superbayesreg & b);

void superbayesreg::create_hregress(void)
  {

  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // --------------------------- method hregress -------------------------------

  vector<ST::string> tnames;
  tnames.push_back("pspline");
  tnames.push_back("hrandom");
  tnames.push_back("random");
  tnames.push_back("spatial");
  tnames.push_back("kriging");
  tnames.push_back("geokriging");
  tnames.push_back("hrandom_pspline");
  tnames.push_back("hrandom_mrf");
  tnames.push_back("mrf_pspline");
  tnames.push_back("hrandomexp_pspline");
  tnames.push_back("ridge");
  tnames.push_back("lasso");
  tnames.push_back("ssvs");
  tnames.push_back("offset");
  tnames.push_back("userdefined");
  tnames.push_back("tensor");
  tnames.push_back("pspline_merror");
  tnames.push_back("multiplicative_pspline");


  tnonp = term_nonp(tnames);
  lineareffects = basic_termtype();
  termtypes.push_back(&tnonp);
  termtypes.push_back(&lineareffects);

  modreg = modelterm(&termtypes);

  udata = use();

  modeonly = simpleoption("modeonly",false);
  setseed = intoption("setseed",-1,0,MAXINT);

  iterations = intoption("iterations",52000,1,10000000);
  burnin = intoption("burnin",2000,0,500000);
  step = intoption("step",50,1,1000);
  level1 = doubleoption("level1",95,40,99);
  level2 = doubleoption("level2",80,40,99);

  standardize = simpleoption("standardize",false);

  families.reserve(80);
  families.push_back("gaussian");
  families.push_back("gaussian_mixture");
  families.push_back("loggaussian");
  families.push_back("quantreg");
  families.push_back("gaussian_re");
  families.push_back("gaussian_exp");
  families.push_back("gaussian_mult");
  families.push_back("binomial_logit");
  families.push_back("za_binomial_logit");
  families.push_back("za_binomial_cloglog");
  families.push_back("binomial_probit");
  families.push_back("binomial_svm");
  families.push_back("binomial_cloglog");
  families.push_back("poisson");
  families.push_back("poisson_ext");
  families.push_back("poisson_extlin");
  families.push_back("multinom_probit");
  families.push_back("multgaussian");
  families.push_back("binomial_logit_l1");
  families.push_back("multinom_logit");
  families.push_back("zip");
  families.push_back("hurdle");
  families.push_back("hurdle_poisson");
  families.push_back("hurdle_negbin");
  families.push_back("zinb");
  families.push_back("zip_cloglog");
  families.push_back("negbin");
  families.push_back("beta");
  families.push_back("lognormal");
  families.push_back("lognormal2");
  families.push_back("normal");
  families.push_back("normal2");
  families.push_back("cnormal");
  families.push_back("truncnormal");
  families.push_back("truncnormal2");
  families.push_back("gamma");
  families.push_back("za_gamma");
  families.push_back("za_lognormal");
  families.push_back("za_invgaussian");
  families.push_back("hetgaussian");
  families.push_back("pareto");
  families.push_back("invgaussian");
  families.push_back("gengamma");
  families.push_back("t");
  families.push_back("bivt");
  families.push_back("bivnormal");
  families.push_back("bivnormal_fz");
  families.push_back("bivprobit");
  families.push_back("bivprobit2");
  families.push_back("bivlogit");
  families.push_back("weibull");
  families.push_back("weibull2");
  families.push_back("dagum");
  families.push_back("betainf0");
  families.push_back("betainf1");
  families.push_back("betainf");
  families.push_back("za");
  families.push_back("BCCG");
  families.push_back("gumbel");
  families.push_back("gumbel2");
  families.push_back("gumbelcopula");
  families.push_back("gumbelcopula2_rho");
  families.push_back("gumbelcopula2_normal_rho");
  families.push_back("gumbelcopula2_normal_mu");
  families.push_back("gumbelcopula2_normal_sigma2");
  families.push_back("gumbelcopula2_normal");
  families.push_back("claytoncopula2_rho");
  families.push_back("claytoncopula");
  families.push_back("claytoncopula2_normal_rho");
  families.push_back("claytoncopula2_normal_mu");
  families.push_back("claytoncopula2_normal_sigma2");
  families.push_back("sfa0_sigma_v");
  families.push_back("sfa0_sigma_u");
  families.push_back("sfa0_mu_y");
  families.push_back("sfa_mu_y_id");
  families.push_back("sfa_mu_u_id");
  families.push_back("sfa_sigma_v");
  families.push_back("sfa_sigma_u");
  families.push_back("sfa_mu_y");
  families.push_back("sfa_mu_u");
  families.push_back("sfa_alpha");
  families.push_back("sfa2_mu_y");
  families.push_back("sfa2_sigma_v");
  families.push_back("sfa2_sigma_u");
  families.push_back("sfa2_mu_u");
  families.push_back("sfa2_mu_u_id");
  families.push_back("sfa2_mu_y_id");
  families.push_back("copulas");
  families.push_back("tcopula");
  families.push_back("gauss_copula");
  families.push_back("gauss_copula2");
  families.push_back("clayton_copula");
  families.push_back("gumbel_copula");
  families.push_back("gaussiancopula");
  families.push_back("frankcopula");
  families.push_back("frankcopula_rho");
  families.push_back("frankcopula2_rho");
  families.push_back("frankcopula2_exp_rho");
  families.push_back("frankcopula_exp");
  families.push_back("frankcopula2_normal_mu");
  families.push_back("frankcopula2_normal_sigma2");
  families.push_back("gaussiancopula_dagum");
  families.push_back("gaussiancopula_binary_normal");
  families.push_back("gaussiancopula_binary_dagum");
  families.push_back("dirichlet");
  families.push_back("gaussian_shared");
  families.push_back("quantreg_shared");
  families.push_back("gaussian_multeffect");

  family = stroption("family",families,"normal");
  aresp = doubleoption("aresp",0.001,-1.0,500);
  bresp = doubleoption("bresp",0.001,0.0,500);
  H = intoption("H",6,2,6);

  hlevel = intoption("hlevel",1,1,3);
  equationnr = intoption("equation",1,1,50);
  equationtypes.reserve(50);

  equationtypes.push_back("main");
  equationtypes.push_back("servant");
  equationtypes.push_back("mu");
  equationtypes.push_back("pi");
  equationtypes.push_back("tau");
  equationtypes.push_back("lambda");
  equationtypes.push_back("delta");
  equationtypes.push_back("sigma2");
  equationtypes.push_back("sigma");
  equationtypes.push_back("sigma_u");
  equationtypes.push_back("sigma_v");
  equationtypes.push_back("scale");
  equationtypes.push_back("scale1");
  equationtypes.push_back("scale2");
  equationtypes.push_back("nu");
  equationtypes.push_back("df");
  equationtypes.push_back("rho");
  equationtypes.push_back("alpha");
  equationtypes.push_back("alpha1");
  equationtypes.push_back("alpha2");
  equationtypes.push_back("alpha3");
  equationtypes.push_back("alpha4");
  equationtypes.push_back("alpha5");
  equationtypes.push_back("alpha6");
  equationtypes.push_back("alpha7");
  equationtypes.push_back("oddsratio");
  equationtypes.push_back("u1");
  equationtypes.push_back("u2");
  equationtypes.push_back("p");
  equationtypes.push_back("b");
  equationtypes.push_back("a");
  equationtypes.push_back("mean");
  equationtypes.push_back("scale");
  equationtypes.push_back("shape2");
  equationtypes.push_back("shape1");
  equationtypes.push_back("latent");
  equationtypes.push_back("mu_shared");

  equationtype = stroption("equationtype",equationtypes,"mu");

  predictop.reserve(20);
  predictop.push_back("no");
  predictop.push_back("nosamples");
  predictop.push_back("full");
  predictop.push_back("fulls");
  predictop.push_back("light");
  predictop.push_back("predictor");

  predict = stroption("predict",predictop,"no");

  WAICoff = simpleoption("WAICoff",false);

  pred_check = simpleoption("pred_check",false);

  cv = simpleoption("cv",false);

  MSEop.reserve(20);
  MSEop.push_back("no");
  MSEop.push_back("yes");
  MSEop.push_back("quadratic");
  MSEop.push_back("check");
  mse = stroption("MSE",MSEop,"no");

  mseparam = doubleoption("mseparam",0.5,-999999999,999999999);

  centerlinear = simpleoption("centerlinear",false);

  quantile = doubleoption("quantile",0.5,0.001,0.999);

  includescale = simpleoption("includescale",false);

  modemaxit = intoption("modemaxit",1000,0,5000);

  stopsum = doubleoption("stopsum",0.999,0.9,1);
  stoprmax = intoption("stoprmax",5,1,100000);

  fraclimit = doubleoption("fraclimit",0.001,0,1);

  scaleconst = simpleoption("scaleconst",false);

  utilities = simpleoption("utilities",false);


  aexp = doubleoption("aexp",0,0.00000000001,100000000);

  bexp = doubleoption("bexp",1,0.00000000001,100000000);

  adaptexp = simpleoption("adaptexp",false);


  slow = simpleoption("slow",false);

  nrbetween = intoption("nrbetween",1000000,1,2000000);

  changelinpredlimits = simpleoption("changelinpredlimits",false);
  linpredminlimit = doubleoption("linpredminlimit",-1000000000,-1000000000,1000000000);
  linpredmaxlimit = doubleoption("linpredmaxlimit",1000000000,-1000000000,1000000000);
  saveestimation = simpleoption("saveestimation",false);

  nrcat = intoption("nrcat",2,1,10);

  imeasures = simpleoption("imeasures",false);

  singleomega = simpleoption("singleomega",false);

  copula = simpleoption("copula",false);

  samplesel = simpleoption("samplesel",false);
  sampleselval = doubleoption("sampleselval",-1000000000,9999,1000000000);

  rotation = intoption("rotation",0,0,270);

  ssvsvarlimit = doubleoption("ssvsvarlimit",0.0000000001,0,0.1);
  IWLSlineff = simpleoption("IWLSlineff", false);
  forceIWLS = simpleoption("forceIWLS", false);

  importance = simpleoption("importance", false);


  regressoptions.reserve(200);

  regressoptions.push_back(&modeonly);
  regressoptions.push_back(&setseed);
  regressoptions.push_back(&iterations);
  regressoptions.push_back(&burnin);
  regressoptions.push_back(&step);
  regressoptions.push_back(&level1);
  regressoptions.push_back(&level2);
  regressoptions.push_back(&family);
  regressoptions.push_back(&aresp);
  regressoptions.push_back(&bresp);
  regressoptions.push_back(&H);
  regressoptions.push_back(&hlevel);
  regressoptions.push_back(&equationnr);
  regressoptions.push_back(&equationtype);
  regressoptions.push_back(&predict);
  regressoptions.push_back(&pred_check);
  regressoptions.push_back(&mse);
  regressoptions.push_back(&mseparam);
  regressoptions.push_back(&centerlinear);
  regressoptions.push_back(&quantile);
  regressoptions.push_back(&cv);
  regressoptions.push_back(&includescale);
  regressoptions.push_back(&standardize);
  regressoptions.push_back(&modemaxit);
  regressoptions.push_back(&stopsum);
  regressoptions.push_back(&stoprmax);
  regressoptions.push_back(&fraclimit);
  regressoptions.push_back(&scaleconst);
  regressoptions.push_back(&utilities);
  regressoptions.push_back(&aexp);
  regressoptions.push_back(&bexp);
  regressoptions.push_back(&adaptexp);
  regressoptions.push_back(&slow);
  regressoptions.push_back(&nrbetween);
  regressoptions.push_back(&changelinpredlimits);
  regressoptions.push_back(&linpredminlimit);
  regressoptions.push_back(&linpredmaxlimit);
  regressoptions.push_back(&saveestimation);
  regressoptions.push_back(&nrcat);
  regressoptions.push_back(&imeasures);
  regressoptions.push_back(&singleomega);
  regressoptions.push_back(&rotation);
  regressoptions.push_back(&copula);
  regressoptions.push_back(&samplesel);
  regressoptions.push_back(&sampleselval);
  regressoptions.push_back(&WAICoff);
  regressoptions.push_back(&ssvsvarlimit);
  regressoptions.push_back(&IWLSlineff);
  regressoptions.push_back(&forceIWLS);
  regressoptions.push_back(&importance);

  // methods 0
  methods.push_back(command("hregress",&modreg,&regressoptions,&udata,required,
			 optional,optional,optional,optional,required));

  functions[0] = hregressrun;
  }

// forward declaration
void autocorrrun(superbayesreg & b);

void superbayesreg::create_autocorr(void)
  {

  // --------------------------- method autocor --------------------------------

  maxlag = intoption("maxlag",250,1,500);

  autocorroptions.push_back(&maxlag);

  // methods 1
  methods.push_back(command("autocor",&ma,&autocorroptions,&ad,notallowed,
						  notallowed,notallowed,notallowed,optional,
                          notallowed));

  functions[1] = autocorrrun;

  }

// forward declaration
void getsamplerun(superbayesreg & b);

void superbayesreg::create_getsample(void)
  {

  // methods 2
  methods.push_back(command("getsample",&mgetsample,&getsampleoptions,
						  &usegetsample,notallowed,notallowed,notallowed,
                          notallowed,notallowed,notallowed));

  functions[2] = getsamplerun;

  }

void superbayesreg::create(void)
  {
#if defined(__BUILDING_LINUX)
  ST::string h = defaultpath+"/output/"+name;
#else
  ST::string h = defaultpath+"\\output\\"+name;
#endif

  outfile = fileoption("outfile",h,false);

  globaloptions.push_back(&outfile);

  create_hregress();
  create_autocorr();
  create_getsample();
  }


void superbayesreg::clear(void)
  {

  equations.erase(equations.begin(),equations.end());
  equations.reserve(200);

  distr_gaussians.erase(distr_gaussians.begin(),distr_gaussians.end());
  distr_gaussians.reserve(20);

  distr_hetgaussians.erase(distr_hetgaussians.begin(),distr_hetgaussians.end());
  distr_hetgaussians.reserve(20);

  distr_vargaussians.erase(distr_vargaussians.begin(),distr_vargaussians.end());
  distr_vargaussians.reserve(20);

  distr_gaussianmixtures.erase(distr_gaussianmixtures.begin(),
  distr_gaussianmixtures.end());
  distr_gaussianmixtures.reserve(20);

  distr_quantregs.erase(distr_quantregs.begin(),distr_quantregs.end());
  distr_quantregs.reserve(20);

  distr_loggaussians.erase(distr_loggaussians.begin(),distr_loggaussians.end());
  distr_loggaussians.reserve(20);

  distr_gaussian_res.erase(distr_gaussian_res.begin(),distr_gaussian_res.end());
  distr_gaussian_res.reserve(100);

  distr_gaussian_exps.erase(distr_gaussian_exps.begin(),distr_gaussian_exps.end());
  distr_gaussian_exps.reserve(20);

  distr_gaussian_mults.erase(distr_gaussian_mults.begin(),
  distr_gaussian_mults.end());
  distr_gaussian_mults.reserve(20);

  distr_binomials.erase(distr_binomials.begin(),distr_binomials.end());
  distr_binomials.reserve(20);

  distr_cloglogs.erase(distr_cloglogs.begin(),distr_cloglogs.end());
  distr_cloglogs.reserve(20);

  distr_poissons.erase(distr_poissons.begin(),distr_poissons.end());
  distr_poissons.reserve(20);

  distr_poisson_exts.erase(distr_poisson_exts.begin(),distr_poisson_exts.end());
  distr_poisson_exts.reserve(20);

  distr_poisson_extlins.erase(distr_poisson_extlins.begin(),
                              distr_poisson_extlins.end());
  distr_poisson_extlins.reserve(20);

  distr_binomialprobits.erase(distr_binomialprobits.begin(),
                              distr_binomialprobits.end());
  distr_binomialprobits.reserve(20);

  distr_binomialsvms.erase(distr_binomialsvms.begin(),
                              distr_binomialsvms.end());
  distr_binomialsvms.reserve(20);

  distr_multinomprobits.erase(distr_multinomprobits.begin(),
                              distr_multinomprobits.end());
  distr_multinomprobits.reserve(20);

  distr_multgaussians.erase(distr_multgaussians.begin(),
                              distr_multgaussians.end());
  distr_multgaussians.reserve(20);

  distr_multinomlogits.erase(distr_multinomlogits.begin(),
                              distr_multinomlogits.end());
  distr_multinomlogits.reserve(20);

  distr_logit_fruehwirths.erase(distr_logit_fruehwirths.begin(),
                               distr_logit_fruehwirths.end());
  distr_logit_fruehwirths.reserve(20);

  distr_ziplambdas.erase(distr_ziplambdas.begin(),distr_ziplambdas.end());
  distr_ziplambdas.reserve(20);

  distr_zippis.erase(distr_zippis.begin(),distr_zippis.end());
  distr_zippis.reserve(20);

  distr_hurdle_lambdas.erase(distr_hurdle_lambdas.begin(),distr_hurdle_lambdas.end());
  distr_hurdle_lambdas.reserve(20);

  distr_hurdle_pis.erase(distr_hurdle_pis.begin(),distr_hurdle_pis.end());
  distr_hurdle_pis.reserve(20);

  distr_hurdle_mus.erase(distr_hurdle_mus.begin(),distr_hurdle_mus.end());
  distr_hurdle_mus.reserve(20);

  distr_hurdle_deltas.erase(distr_hurdle_deltas.begin(),distr_hurdle_deltas.end());
  distr_hurdle_deltas.reserve(20);

  distr_negbinzip_mus.erase(distr_negbinzip_mus.begin(),distr_negbinzip_mus.end());
  distr_negbinzip_mus.reserve(20);

  distr_negbinzip_pis.erase(distr_negbinzip_pis.begin(),distr_negbinzip_pis.end());
  distr_negbinzip_pis.reserve(20);

  distr_negbinzip_deltas.erase(distr_negbinzip_deltas.begin(),distr_negbinzip_deltas.end());
  distr_negbinzip_deltas.reserve(20);

  distr_zip_cloglog_mus.erase(distr_zip_cloglog_mus.begin(),distr_zip_cloglog_mus.end());
  distr_zip_cloglog_mus.reserve(20);

  distr_zip_cloglog_pis.erase(distr_zip_cloglog_pis.begin(),distr_zip_cloglog_pis.end());
  distr_zip_cloglog_pis.reserve(20);

  distr_negbin_mus.erase(distr_negbin_mus.begin(),distr_negbin_mus.end());
  distr_negbin_mus.reserve(20);

  distr_negbin_deltas.erase(distr_negbin_deltas.begin(),distr_negbin_deltas.end());
  distr_negbin_deltas.reserve(20);

  distr_beta_mus.erase(distr_beta_mus.begin(),distr_beta_mus.end());
  distr_beta_mus.reserve(20);

  distr_beta_sigma2s.erase(distr_beta_sigma2s.begin(),distr_beta_sigma2s.end());
  distr_beta_sigma2s.reserve(20);

  distr_lognormal_mus.erase(distr_lognormal_mus.begin(),distr_lognormal_mus.end());
  distr_lognormal_mus.reserve(20);

  distr_lognormal_sigma2s.erase(distr_lognormal_sigma2s.begin(),distr_lognormal_sigma2s.end());
  distr_lognormal_sigma2s.reserve(20);

  distr_lognormal2_mus.erase(distr_lognormal2_mus.begin(),distr_lognormal2_mus.end());
  distr_lognormal2_mus.reserve(20);

  distr_lognormal2_sigmas.erase(distr_lognormal2_sigmas.begin(),distr_lognormal2_sigmas.end());
  distr_lognormal2_sigmas.reserve(20);

  distr_normal_mus.erase(distr_normal_mus.begin(),distr_normal_mus.end());
  distr_normal_mus.reserve(20);

  distr_cnormal_mus.erase(distr_cnormal_mus.begin(),distr_cnormal_mus.end());
  distr_cnormal_mus.reserve(20);

  distr_normal_sigma2s.erase(distr_normal_sigma2s.begin(),distr_normal_sigma2s.end());
  distr_normal_sigma2s.reserve(20);

  distr_normal2_mus.erase(distr_normal2_mus.begin(),distr_normal2_mus.end());
  distr_normal2_mus.reserve(20);

  distr_normal2_sigmas.erase(distr_normal2_sigmas.begin(),distr_normal2_sigmas.end());
  distr_normal2_sigmas.reserve(20);

  distr_truncnormal2_mus.erase(distr_truncnormal2_mus.begin(),distr_truncnormal2_mus.end());
  distr_truncnormal2_mus.reserve(20);

  distr_truncnormal2_sigmas.erase(distr_truncnormal2_sigmas.begin(),distr_truncnormal2_sigmas.end());
  distr_truncnormal2_sigmas.reserve(20);

  distr_invgaussian_mus.erase(distr_invgaussian_mus.begin(),distr_invgaussian_mus.end());
  distr_invgaussian_mus.reserve(20);

  distr_invgaussian_sigma2s.erase(distr_invgaussian_sigma2s.begin(),distr_invgaussian_sigma2s.end());
  distr_invgaussian_sigma2s.reserve(20);

  distr_gamma_mus.erase(distr_gamma_mus.begin(),distr_gamma_mus.end());
  distr_gamma_mus.reserve(20);

  distr_gamma_sigmas.erase(distr_gamma_sigmas.begin(),distr_gamma_sigmas.end());
  distr_gamma_sigmas.reserve(20);

  distr_pareto_bs.erase(distr_pareto_bs.begin(),distr_pareto_bs.end());
  distr_pareto_bs.reserve(20);

  distr_pareto_ps.erase(distr_pareto_ps.begin(),distr_pareto_ps.end());
  distr_pareto_ps.reserve(20);

  distr_gengamma_mus.erase(distr_gengamma_mus.begin(),distr_gengamma_mus.end());
  distr_gengamma_mus.reserve(20);

  distr_gengamma_sigmas.erase(distr_gengamma_sigmas.begin(),distr_gengamma_sigmas.end());
  distr_gengamma_sigmas.reserve(20);

  distr_gengamma_taus.erase(distr_gengamma_taus.begin(),distr_gengamma_taus.end());
  distr_gengamma_taus.reserve(20);

  distr_t_mus.erase(distr_t_mus.begin(),distr_t_mus.end());
  distr_t_mus.reserve(20);

  distr_t_sigma2s.erase(distr_t_sigma2s.begin(),distr_t_sigma2s.end());
  distr_t_sigma2s.reserve(20);

  distr_t_dfs.erase(distr_t_dfs.begin(),distr_t_dfs.end());
  distr_t_dfs.reserve(20);

  distr_BCCG_mus.erase(distr_BCCG_mus.begin(),distr_BCCG_mus.end());
  distr_BCCG_mus.reserve(20);

  distr_BCCG_sigmas.erase(distr_BCCG_sigmas.begin(),distr_BCCG_sigmas.end());
  distr_BCCG_sigmas.reserve(20);

  distr_BCCG_nus.erase(distr_BCCG_nus.begin(),distr_BCCG_nus.end());
  distr_BCCG_nus.reserve(20);

  distr_bivt_mus.erase(distr_bivt_mus.begin(),distr_bivt_mus.end());
  distr_bivt_mus.reserve(20);

  distr_gumbel_mus.erase(distr_gumbel_mus.begin(),distr_gumbel_mus.end());
  distr_gumbel_mus.reserve(20);

  distr_gumbel_sigmas.erase(distr_gumbel_sigmas.begin(),distr_gumbel_sigmas.end());
  distr_gumbel_sigmas.reserve(20);

  distr_gumbel2_mus.erase(distr_gumbel2_mus.begin(),distr_gumbel2_mus.end());
  distr_gumbel2_mus.reserve(20);

  distr_gumbel2_sigma2s.erase(distr_gumbel2_sigma2s.begin(),distr_gumbel2_sigma2s.end());
  distr_gumbel2_sigma2s.reserve(20);

  distr_bivt_sigmas.erase(distr_bivt_sigmas.begin(),distr_bivt_sigmas.end());
  distr_bivt_sigmas.reserve(20);

  distr_bivt_rhos.erase(distr_bivt_rhos.begin(),distr_bivt_rhos.end());
  distr_bivt_rhos.reserve(20);

  distr_bivt_dfs.erase(distr_bivt_dfs.begin(),distr_bivt_dfs.end());
  distr_bivt_dfs.reserve(20);

  distr_bivnormal_mufzs.erase(distr_bivnormal_mufzs.begin(),distr_bivnormal_mufzs.end());
  distr_bivnormal_mufzs.reserve(20);

  distr_bivnormal_mus.erase(distr_bivnormal_mus.begin(),distr_bivnormal_mus.end());
  distr_bivnormal_mus.reserve(20);

  distr_bivnormal_sigmas.erase(distr_bivnormal_sigmas.begin(),distr_bivnormal_sigmas.end());
  distr_bivnormal_sigmas.reserve(20);

  distr_bivnormal_rhos.erase(distr_bivnormal_rhos.begin(),distr_bivnormal_rhos.end());
  distr_bivnormal_rhos.reserve(20);

  distr_bivnormal_rhofzs.erase(distr_bivnormal_rhofzs.begin(),distr_bivnormal_rhofzs.end());
  distr_bivnormal_rhofzs.reserve(20);

  distr_bivprobit_mus.erase(distr_bivprobit_mus.begin(),distr_bivprobit_mus.end());
  distr_bivprobit_mus.reserve(20);

  distr_bivprobit_rhos.erase(distr_bivprobit_rhos.begin(),distr_bivprobit_rhos.end());
  distr_bivprobit_rhos.reserve(20);

  distr_bivprobit2_mus.erase(distr_bivprobit2_mus.begin(),distr_bivprobit2_mus.end());
  distr_bivprobit2_mus.reserve(20);

  distr_bivprobit2_rhos.erase(distr_bivprobit2_rhos.begin(),distr_bivprobit2_rhos.end());
  distr_bivprobit2_rhos.reserve(20);

  distr_bivlogit_mus.erase(distr_bivlogit_mus.begin(),distr_bivlogit_mus.end());
  distr_bivlogit_mus.reserve(20);

  distr_bivlogit_ors.erase(distr_bivlogit_ors.begin(),distr_bivlogit_ors.end());
  distr_bivlogit_ors.reserve(20);

  distr_zeroadjusteds.erase(distr_zeroadjusteds.begin(),distr_zeroadjusteds.end());
  distr_zeroadjusteds.reserve(5);

  distr_dirichlets.erase(distr_dirichlets.begin(),distr_dirichlets.end());
  distr_dirichlets.reserve(20);

  distr_zeroadjusted_mults.erase(distr_zeroadjusted_mults.begin(),distr_zeroadjusted_mults.end());
  distr_zeroadjusted_mults.reserve(5);

  distr_weibull_lambdas.erase(distr_weibull_lambdas.begin(),distr_weibull_lambdas.end());
  distr_weibull_lambdas.reserve(20);

  distr_weibull_alphas.erase(distr_weibull_alphas.begin(),distr_weibull_alphas.end());
  distr_weibull_alphas.reserve(20);

  distr_weibull2_lambdas.erase(distr_weibull2_lambdas.begin(),distr_weibull2_lambdas.end());
  distr_weibull2_lambdas.reserve(20);

  distr_weibull2_alphas.erase(distr_weibull2_alphas.begin(),distr_weibull2_alphas.end());
  distr_weibull2_alphas.reserve(20);

  distr_dagum_as.erase(distr_dagum_as.begin(),distr_dagum_as.end());
  distr_dagum_as.reserve(20);

  distr_dagum_bs.erase(distr_dagum_bs.begin(),distr_dagum_bs.end());
  distr_dagum_bs.reserve(20);

  distr_dagum_ps.erase(distr_dagum_ps.begin(),distr_dagum_ps.end());
  distr_dagum_ps.reserve(20);

  distr_betainf_mus.erase(distr_betainf_mus.begin(),distr_betainf_mus.end());
  distr_betainf_mus.reserve(20);

  distr_betainf_sigma2s.erase(distr_betainf_sigma2s.begin(),distr_betainf_sigma2s.end());
  distr_betainf_sigma2s.reserve(20);

  distr_betainf_nus.erase(distr_betainf_nus.begin(),distr_betainf_nus.end());
  distr_betainf_nus.reserve(20);

  distr_betainf_taus.erase(distr_betainf_taus.begin(),distr_betainf_taus.end());
  distr_betainf_taus.reserve(20);

  distr_betainf0_nus.erase(distr_betainf0_nus.begin(),distr_betainf0_nus.end());
  distr_betainf0_nus.reserve(20);

  distr_betainf1_taus.erase(distr_betainf1_taus.begin(),distr_betainf1_taus.end());
  distr_betainf1_taus.reserve(20);

  distr_gumbelcopula_rhos.erase(distr_gumbelcopula_rhos.begin(),distr_gumbelcopula_rhos.end());
  distr_gumbelcopula_rhos.reserve(20);

  distr_gaussiancopula_binary_normal_rhos.erase(distr_gaussiancopula_binary_normal_rhos.begin(),distr_gaussiancopula_binary_normal_rhos.end());
  distr_gaussiancopula_binary_normal_rhos.reserve(20);

  distr_gaussiancopula_binary_normal_sigmas.erase(distr_gaussiancopula_binary_normal_sigmas.begin(),distr_gaussiancopula_binary_normal_sigmas.end());
  distr_gaussiancopula_binary_normal_sigmas.reserve(20);

  distr_gaussiancopula_binary_normal_mus.erase(distr_gaussiancopula_binary_normal_mus.begin(),distr_gaussiancopula_binary_normal_mus.end());
  distr_gaussiancopula_binary_normal_mus.reserve(20);

  distr_gaussiancopula_binary_normal_latents.erase(distr_gaussiancopula_binary_normal_latents.begin(),distr_gaussiancopula_binary_normal_latents.end());
  distr_gaussiancopula_binary_normal_latents.reserve(20);

  distr_gaussiancopula_binary_dagum_rhos.erase(distr_gaussiancopula_binary_dagum_rhos.begin(),distr_gaussiancopula_binary_dagum_rhos.end());
  distr_gaussiancopula_binary_dagum_rhos.reserve(20);

  distr_gaussiancopula_binary_dagum_as.erase(distr_gaussiancopula_binary_dagum_as.begin(),distr_gaussiancopula_binary_dagum_as.end());
  distr_gaussiancopula_binary_dagum_as.reserve(20);

  distr_gaussiancopula_binary_dagum_bs.erase(distr_gaussiancopula_binary_dagum_bs.begin(),distr_gaussiancopula_binary_dagum_bs.end());
  distr_gaussiancopula_binary_dagum_bs.reserve(20);

  distr_gaussiancopula_binary_dagum_ps.erase(distr_gaussiancopula_binary_dagum_ps.begin(),distr_gaussiancopula_binary_dagum_ps.end());
  distr_gaussiancopula_binary_dagum_ps.reserve(20);

  distr_gaussiancopula_binary_dagum_latents.erase(distr_gaussiancopula_binary_dagum_latents.begin(),distr_gaussiancopula_binary_dagum_latents.end());
  distr_gaussiancopula_binary_dagum_latents.reserve(20);

  distr_gaussiancopula_dagum_rhos.erase(distr_gaussiancopula_dagum_rhos.begin(),distr_gaussiancopula_dagum_rhos.end());
  distr_gaussiancopula_dagum_rhos.reserve(20);

  distr_gaussiancopula_dagum_as.erase(distr_gaussiancopula_dagum_as.begin(),distr_gaussiancopula_dagum_as.end());
  distr_gaussiancopula_dagum_as.reserve(20);

  distr_gaussiancopula_dagum_bs.erase(distr_gaussiancopula_dagum_bs.begin(),distr_gaussiancopula_dagum_bs.end());
  distr_gaussiancopula_dagum_bs.reserve(20);

  distr_gaussiancopula_dagum_ps.erase(distr_gaussiancopula_dagum_ps.begin(),distr_gaussiancopula_dagum_ps.end());
  distr_gaussiancopula_dagum_ps.reserve(20);

  distr_gumbelcopula2_rhos.erase(distr_gumbelcopula2_rhos.begin(),distr_gumbelcopula2_rhos.end());
  distr_gumbelcopula2_rhos.reserve(20);

  distr_gumbelcopula2_normal_mus.erase(distr_gumbelcopula2_normal_mus.begin(),distr_gumbelcopula2_normal_mus.end());
  distr_gumbelcopula2_normal_mus.reserve(20);

  distr_gumbelcopula2_normal_sigma2s.erase(distr_gumbelcopula2_normal_sigma2s.begin(),distr_gumbelcopula2_normal_sigma2s.end());
  distr_gumbelcopula2_normal_sigma2s.reserve(20);

  distr_gumbelcopula2_normal_mu_2s.erase(distr_gumbelcopula2_normal_mu_2s.begin(),distr_gumbelcopula2_normal_mu_2s.end());
  distr_gumbelcopula2_normal_mu_2s.reserve(20);

  distr_gumbelcopula2_normal_sigma2_2s.erase(distr_gumbelcopula2_normal_sigma2_2s.begin(),distr_gumbelcopula2_normal_sigma2_2s.end());
  distr_gumbelcopula2_normal_sigma2_2s.reserve(20);

  distr_claytoncopula_rhos.erase(distr_claytoncopula_rhos.begin(),distr_claytoncopula_rhos.end());
  distr_claytoncopula_rhos.reserve(20);

  distr_claytoncopula2_rhos.erase(distr_claytoncopula2_rhos.begin(),distr_claytoncopula2_rhos.end());
  distr_claytoncopula2_rhos.reserve(20);

  distr_claytoncopula2_normal_mus.erase(distr_claytoncopula2_normal_mus.begin(),distr_claytoncopula2_normal_mus.end());
  distr_claytoncopula2_normal_mus.reserve(20);

  distr_claytoncopula2_normal_sigma2s.erase(distr_claytoncopula2_normal_sigma2s.begin(),distr_claytoncopula2_normal_sigma2s.end());
  distr_claytoncopula2_normal_sigma2s.reserve(20);

  distr_frankcopula_rhos.erase(distr_frankcopula_rhos.begin(),distr_frankcopula_rhos.end());
  distr_frankcopula_rhos.reserve(20);

  distr_frankcopula2_rhos.erase(distr_frankcopula2_rhos.begin(),distr_frankcopula2_rhos.end());
  distr_frankcopula2_rhos.reserve(20);

  distr_frankcopula2_exp_rhos.erase(distr_frankcopula2_exp_rhos.begin(),distr_frankcopula2_exp_rhos.end());
  distr_frankcopula2_exp_rhos.reserve(20);

  distr_frankcopula_exp_rhos.erase(distr_frankcopula_exp_rhos.begin(),distr_frankcopula_exp_rhos.end());
  distr_frankcopula_exp_rhos.reserve(20);

  distr_frankcopula2_normal_mus.erase(distr_frankcopula2_normal_mus.begin(),distr_frankcopula2_normal_mus.end());
  distr_frankcopula2_normal_mus.reserve(20);

  distr_frankcopula2_normal_sigma2s.erase(distr_frankcopula2_normal_sigma2s.begin(),distr_frankcopula2_normal_sigma2s.end());
  distr_frankcopula2_normal_sigma2s.reserve(20);

  distr_sfa0_mu_ys.erase(distr_sfa0_mu_ys.begin(),distr_sfa0_mu_ys.end());
  distr_sfa0_mu_ys.reserve(20);

  distr_sfa0_sigma_us.erase(distr_sfa0_sigma_us.begin(),distr_sfa0_sigma_us.end());
  distr_sfa0_sigma_us.reserve(20);

  distr_sfa0_sigma_vs.erase(distr_sfa0_sigma_vs.begin(),distr_sfa0_sigma_vs.end());
  distr_sfa0_sigma_vs.reserve(20);

  distr_sfa_mu_y_ids.erase(distr_sfa_mu_y_ids.begin(),distr_sfa_mu_y_ids.end());
  distr_sfa_mu_y_ids.reserve(20);

  distr_sfa_mu_u_ids.erase(distr_sfa_mu_u_ids.begin(),distr_sfa_mu_u_ids.end());
  distr_sfa_mu_u_ids.reserve(20);

  distr_sfa_mu_ys.erase(distr_sfa_mu_ys.begin(),distr_sfa_mu_ys.end());
  distr_sfa_mu_ys.reserve(20);

  distr_sfa_mu_us.erase(distr_sfa_mu_us.begin(),distr_sfa_mu_us.end());
  distr_sfa_mu_us.reserve(20);

  distr_sfa_sigma_us.erase(distr_sfa_sigma_us.begin(),distr_sfa_sigma_us.end());
  distr_sfa_sigma_us.reserve(20);

  distr_sfa_sigma_vs.erase(distr_sfa_sigma_vs.begin(),distr_sfa_sigma_vs.end());
  distr_sfa_sigma_vs.reserve(20);

  distr_sfa_alphas.erase(distr_sfa_alphas.begin(),distr_sfa_alphas.end());
  distr_sfa_alphas.reserve(20);

  distr_sfa2_mu_ys.erase(distr_sfa2_mu_ys.begin(),distr_sfa2_mu_ys.end());
  distr_sfa2_mu_ys.reserve(20);

  distr_sfa2_mu_u_ids.erase(distr_sfa2_mu_u_ids.begin(),distr_sfa2_mu_u_ids.end());
  distr_sfa2_mu_u_ids.reserve(20);

  distr_sfa2_mu_y_ids.erase(distr_sfa2_mu_y_ids.begin(),distr_sfa2_mu_y_ids.end());
  distr_sfa2_mu_y_ids.reserve(20);

  distr_sfa2_mu_us.erase(distr_sfa2_mu_us.begin(),distr_sfa2_mu_us.end());
  distr_sfa2_mu_us.reserve(20);

  distr_sfa2_sigma_us.erase(distr_sfa2_sigma_us.begin(),distr_sfa2_sigma_us.end());
  distr_sfa2_sigma_us.reserve(20);

  distr_sfa2_sigma_vs.erase(distr_sfa2_sigma_vs.begin(),distr_sfa2_sigma_vs.end());
  distr_sfa2_sigma_vs.reserve(20);

  distr_copulas.erase(distr_copulas.begin(),distr_copulas.end());
  distr_copulas.reserve(20);

  distr_tcopula_dfs.erase(distr_tcopula_dfs.begin(),distr_tcopula_dfs.end());
  distr_tcopula_dfs.reserve(20);

  distr_tcopula_rhos.erase(distr_tcopula_rhos.begin(),distr_tcopula_rhos.end());
  distr_tcopula_rhos.reserve(20);

  distr_gaussiancopula_rhos.erase(distr_gaussiancopula_rhos.begin(),distr_gaussiancopula_rhos.end());
  distr_gaussiancopula_rhos.reserve(20);

  distr_gaussiancopula_rhofzs.erase(distr_gaussiancopula_rhofzs.begin(),distr_gaussiancopula_rhofzs.end());
  distr_gaussiancopula_rhofzs.reserve(20);

  distr_gausscopulas.erase(distr_gausscopulas.begin(),distr_gausscopulas.end());
  distr_gausscopulas.reserve(20);

  distr_gausscopula2s.erase(distr_gausscopula2s.begin(),distr_gausscopula2s.end());
  distr_gausscopula2s.reserve(20);

  distr_clayton_copulas.erase(distr_clayton_copulas.begin(),distr_clayton_copulas.end());
  distr_clayton_copulas.reserve(20);

  distr_gumbel_copulas.erase(distr_gumbel_copulas.begin(),distr_gumbel_copulas.end());
  distr_gumbel_copulas.reserve(20);

  distr_binomialprobit_copulas.erase(distr_binomialprobit_copulas.begin(),
                              distr_binomialprobit_copulas.end());
  distr_binomialprobit_copulas.reserve(20);

  distr_JMs.erase(distr_JMs.begin(),distr_JMs.end());
  distr_JMs.reserve(20);

  distr_gaussian_multeffects.erase(distr_gaussian_multeffects.begin(),distr_gaussian_multeffects.end());
  distr_gaussian_multeffects.reserve(20);

  FC_linears.erase(FC_linears.begin(),FC_linears.end());
  FC_linears.reserve(200);

  design_psplines.erase(design_psplines.begin(),design_psplines.end());
  design_psplines.reserve(200);

  design_userdefineds.erase(design_userdefineds.begin(),design_userdefineds.end());
  design_userdefineds.reserve(200);

  design_userdefined_tensors.erase(design_userdefined_tensors.begin(),design_userdefined_tensors.end());
  design_userdefined_tensors.reserve(200);

  FC_tensor_omegas.erase(FC_tensor_omegas.begin(),FC_tensor_omegas.end());
  FC_tensor_omegas.reserve(200);

  design_mrfs.erase(design_mrfs.begin(),design_mrfs.end());
  design_mrfs.reserve(200);

  design_krigings.erase(design_krigings.begin(),design_krigings.end());
  design_krigings.reserve(200);

  design_hrandoms.erase(design_hrandoms.begin(),design_hrandoms.end());
  design_hrandoms.reserve(200);

  FC_merrors.erase(FC_merrors.begin(),FC_merrors.end());
  FC_merrors.reserve(20);

  FC_nonps.erase(FC_nonps.begin(),FC_nonps.end());
  FC_nonps.reserve(400);

  FC_hrandoms.erase(FC_hrandoms.begin(),FC_hrandoms.end());
  FC_hrandoms.reserve(200);

  FC_mults.erase(FC_mults.begin(),FC_mults.end());
  FC_mults.reserve(200);

  FC_nonp_variances.erase(FC_nonp_variances.begin(),FC_nonp_variances.end());
  FC_nonp_variances.reserve(400);

  FC_nonp_variance_varselections.erase(FC_nonp_variance_varselections.begin(),
                                       FC_nonp_variance_varselections.end());
  FC_nonp_variance_varselections.reserve(400);

  FC_varselection_omegas.erase(FC_varselection_omegas.begin(),
                                       FC_varselection_omegas.end());
  FC_varselection_omegas.reserve(200);

  FC_linear_pens.erase(FC_linear_pens.begin(),FC_linear_pens.end());
  FC_linear_pens.reserve(200);

  FC_variance_pen_vectors.erase(FC_variance_pen_vectors.begin(),
                                FC_variance_pen_vectors.end());
  FC_variance_pen_vectors.reserve(200);

  FC_variance_pen_vector_ssvss.erase(FC_variance_pen_vector_ssvss.begin(),
                                FC_variance_pen_vector_ssvss.end());
  FC_variance_pen_vector_ssvss.reserve(200);

  FC_hrandom_variances.erase(FC_hrandom_variances.begin(),
  FC_hrandom_variances.end());
  FC_hrandom_variances.reserve(200);

  FC_hrandom_variance_vecs.erase(FC_hrandom_variance_vecs.begin(),
  FC_hrandom_variance_vecs.end());
  FC_hrandom_variance_vecs.reserve(200);

  FC_hrandom_variance_vec_nmigs.erase(FC_hrandom_variance_vec_nmigs.begin(),
  FC_hrandom_variance_vec_nmigs.end());
  FC_hrandom_variance_vec_nmigs.reserve(200);

  FC_hrandom_variance_ssvss.erase(FC_hrandom_variance_ssvss.begin(),
  FC_hrandom_variance_ssvss.end());
  FC_hrandom_variance_ssvss.reserve(200);

  FC_predicts.erase(FC_predicts.begin(),FC_predicts.end());
  FC_predicts.reserve(200);

  FC_predicts_mult.erase(FC_predicts_mult.begin(),FC_predicts_mult.end());
  FC_predicts_mult.reserve(200);

  predict_mult_distrs.erase(predict_mult_distrs.begin(),predict_mult_distrs.end());
  predict_mult_distrs.reserve(200);

  FC_predict_predictors.erase(FC_predict_predictors.begin(),
                              FC_predict_predictors.end());
  FC_predict_predictors.reserve(200);


  FC_predictive_checks.erase(FC_predictive_checks.begin(),
  FC_predictive_checks.end());
  FC_predictive_checks.reserve(200);

  }


#if defined(JAVA_OUTPUT_WINDOW)
superbayesreg::superbayesreg(
administrator_basic * adb, administrator_pointer * adp,
const ST::string & n,ofstream * lo,istream * in,
						 ST::string p,vector<statobject*> * st)
						 : statobject(adb,n,"mcmcreg",lo,in,p)
  {
  clear();
  adminp_p = adp;
  statobj = st;
  master = MASTER_OBJ();
  create();
  resultsyesno = false;
  run_yes = false;
  posteriormode = false;
  computemodeforstartingvalues = true;
  describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
  }
#else
superbayesreg::superbayesreg(const ST::string & n,ofstream * lo,istream * in,
						 ST::string p,vector<statobject*> * st)
						 : statobject(n,"mcmcreg",lo,in,p)
  {
  clear();
  statobj = st;
  master =MASTER_OBJ();
  create();
  resultsyesno = false;
  run_yes = false;
  mainequation=false;
  shared=false;
  posteriormode = false;
  computemodeforstartingvalues = true;
  firstvarselection = true;
  describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
  }
#endif


superbayesreg::superbayesreg(const superbayesreg & b) : statobject(statobject(b))
  {
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = b.adminp_p;
  #endif

  firstvarselection = b.firstvarselection;

  pathres = b.pathres;
  title = b.title;
  pathnonp = b.pathnonp;

  statobj = b.statobj;

  D = b.D;

  modelvarnamesv = b.modelvarnamesv;

  equations=b.equations;
  mainequation=b.mainequation;
  countmarginal=b.countmarginal;
  nrlevel1 = b.nrlevel1;
  simobj = b.simobj;

  shared = b.shared;

  master = b.master;

  generaloptions = b.generaloptions;

  distr_gaussians = b.distr_gaussians;
  distr_gaussian_multeffects = b.distr_gaussian_multeffects;
  distr_hetgaussians = b.distr_hetgaussians;
  distr_vargaussians = b.distr_vargaussians;
  distr_quantregs = b.distr_quantregs;
  distr_gaussianmixtures = b.distr_gaussianmixtures;
  distr_loggaussians = b.distr_loggaussians;
  distr_gaussian_res = b.distr_gaussian_res;
  distr_gaussian_exps = b.distr_gaussian_exps;
  distr_gaussian_mults = b.distr_gaussian_mults;
  distr_binomials = b.distr_binomials;
  distr_cloglogs = b.distr_cloglogs;
  distr_poissons = b.distr_poissons;
  distr_poisson_exts = b.distr_poisson_exts;
  distr_poisson_extlins = b.distr_poisson_extlins;
  distr_binomialprobits = b.distr_binomialprobits;
  distr_binomialsvms = b.distr_binomialsvms;
  distr_multinomprobits = b.distr_multinomprobits;
  distr_multgaussians = b.distr_multgaussians;
  distr_multinomlogits = b.distr_multinomlogits;
  distr_logit_fruehwirths = b.distr_logit_fruehwirths;
  distr_ziplambdas = b.distr_ziplambdas;
  distr_zippis = b.distr_zippis;
  distr_hurdle_lambdas = b.distr_hurdle_lambdas;
  distr_hurdle_pis = b.distr_hurdle_pis;
  distr_hurdle_mus = b.distr_hurdle_mus;
  distr_hurdle_deltas = b.distr_hurdle_deltas;
  distr_negbinzip_mus = b.distr_negbinzip_mus;
  distr_negbinzip_pis = b.distr_negbinzip_pis;
  distr_negbinzip_deltas = b.distr_negbinzip_deltas;
  distr_zip_cloglog_mus = b.distr_zip_cloglog_mus;
  distr_zip_cloglog_pis = b.distr_zip_cloglog_pis;
  distr_negbin_mus = b.distr_negbin_mus;
  distr_negbin_deltas = b.distr_negbin_deltas;
  distr_beta_mus = b.distr_beta_mus;
  distr_beta_sigma2s = b.distr_beta_sigma2s;
  distr_lognormal_mus = b.distr_lognormal_mus;
  distr_lognormal_sigma2s = b.distr_lognormal_sigma2s;
  distr_lognormal2_mus = b.distr_lognormal2_mus;
  distr_lognormal2_sigmas = b.distr_lognormal2_sigmas;
  distr_normal_mus = b.distr_normal_mus;
  distr_normal_sigma2s = b.distr_normal_sigma2s;
  distr_cnormal_mus = b.distr_cnormal_mus;
  distr_cnormal_sigmas = b.distr_cnormal_sigmas;
  distr_normal2_mus = b.distr_normal2_mus;
  distr_normal2_sigmas = b.distr_normal2_sigmas;
  distr_truncnormal2_mus = b.distr_truncnormal2_mus;
  distr_truncnormal2_sigmas = b.distr_truncnormal2_sigmas;
  distr_invgaussian_mus = b.distr_invgaussian_mus;
  distr_invgaussian_sigma2s = b.distr_invgaussian_sigma2s;
  distr_gamma_mus = b.distr_gamma_mus;
  distr_gamma_sigmas = b.distr_gamma_sigmas;
  distr_pareto_bs = b.distr_pareto_bs;
  distr_pareto_ps = b.distr_pareto_ps;
  distr_gengamma_mus = b.distr_gengamma_mus;
  distr_gengamma_sigmas = b.distr_gengamma_sigmas;
  distr_gengamma_taus = b.distr_gengamma_taus;
  distr_t_mus = b.distr_t_mus;
  distr_t_sigma2s = b.distr_t_sigma2s;
  distr_t_dfs = b.distr_t_dfs;
  distr_BCCG_mus = b.distr_BCCG_mus;
  distr_BCCG_sigmas = b.distr_BCCG_sigmas;
  distr_BCCG_nus = b.distr_BCCG_nus;
  distr_bivt_mus = b.distr_bivt_mus;
  distr_bivt_sigmas = b.distr_bivt_sigmas;
  distr_bivt_dfs = b.distr_bivt_dfs;
  distr_bivt_rhos = b.distr_bivt_rhos;
  distr_bivnormal_mus = b.distr_bivnormal_mus;
  distr_bivnormal_mufzs = b.distr_bivnormal_mufzs;
  distr_bivnormal_sigmas = b.distr_bivnormal_sigmas;
  distr_bivnormal_rhos = b.distr_bivnormal_rhos;
  distr_bivnormal_rhofzs = b.distr_bivnormal_rhofzs;
  distr_bivprobit_mus = b.distr_bivprobit_mus;
  distr_bivprobit_rhos = b.distr_bivprobit_rhos;
  distr_bivprobit2_mus = b.distr_bivprobit2_mus;
  distr_bivprobit2_rhos = b.distr_bivprobit2_rhos;
  distr_bivlogit_mus = b.distr_bivlogit_mus;
  distr_bivlogit_ors = b.distr_bivlogit_ors;
  distr_zeroadjusteds = b.distr_zeroadjusteds;
  distr_zeroadjusted_mults = b.distr_zeroadjusted_mults;
  distr_dirichlets = b.distr_dirichlets;
  distr_weibull_lambdas = b.distr_weibull_lambdas;
  distr_weibull_alphas = b.distr_weibull_alphas;
  distr_weibull2_lambdas = b.distr_weibull2_lambdas;
  distr_weibull2_alphas = b.distr_weibull2_alphas;
  distr_dagum_as = b.distr_dagum_as;
  distr_dagum_bs = b.distr_dagum_bs;
  distr_dagum_ps = b.distr_dagum_ps;
  distr_betainf_mus = b.distr_betainf_mus;
  distr_betainf_sigma2s = b.distr_betainf_sigma2s;
  distr_betainf_nus = b.distr_betainf_nus;
  distr_betainf_taus = b.distr_betainf_taus;
  distr_betainf0_nus = b.distr_betainf0_nus;
  distr_betainf1_taus = b.distr_betainf1_taus;
  distr_gumbel_mus = b.distr_gumbel_mus;
  distr_gumbel_sigmas = b.distr_gumbel_sigmas;
  distr_gumbel2_mus = b.distr_gumbel2_mus;
  distr_gumbel2_sigma2s = b.distr_gumbel2_sigma2s;
  distr_gumbelcopula_rhos = b.distr_gumbelcopula_rhos;
  distr_gumbelcopula2_rhos = b.distr_gumbelcopula2_rhos;
  distr_gumbelcopula2_normal_mus = b.distr_gumbelcopula2_normal_mus;
  distr_gumbelcopula2_normal_sigma2s = b.distr_gumbelcopula2_normal_sigma2s;
  distr_gumbelcopula2_normal_mu_2s = b.distr_gumbelcopula2_normal_mu_2s;
  distr_gumbelcopula2_normal_sigma2_2s = b.distr_gumbelcopula2_normal_sigma2_2s;
  distr_gaussiancopula_binary_normal_rhos = b.distr_gaussiancopula_binary_normal_rhos;
  distr_gaussiancopula_binary_normal_sigmas = b.distr_gaussiancopula_binary_normal_sigmas;
  distr_gaussiancopula_binary_normal_mus = b.distr_gaussiancopula_binary_normal_mus;
  distr_gaussiancopula_binary_normal_latents = b.distr_gaussiancopula_binary_normal_latents;
  distr_gaussiancopula_binary_dagum_rhos = b.distr_gaussiancopula_binary_dagum_rhos;
  distr_gaussiancopula_binary_dagum_as = b.distr_gaussiancopula_binary_dagum_as;
  distr_gaussiancopula_binary_dagum_bs = b.distr_gaussiancopula_binary_dagum_bs;
  distr_gaussiancopula_binary_dagum_ps = b.distr_gaussiancopula_binary_dagum_ps;
  distr_gaussiancopula_binary_dagum_latents = b.distr_gaussiancopula_binary_dagum_latents;
  distr_gaussiancopula_dagum_rhos = b.distr_gaussiancopula_dagum_rhos;
  distr_gaussiancopula_dagum_as = b.distr_gaussiancopula_dagum_as;
  distr_gaussiancopula_dagum_bs = b.distr_gaussiancopula_dagum_bs;
  distr_gaussiancopula_dagum_ps = b.distr_gaussiancopula_dagum_ps;
  distr_claytoncopula_rhos = b.distr_claytoncopula_rhos;
  distr_claytoncopula2_rhos = b.distr_claytoncopula2_rhos;
  distr_claytoncopula2_normal_mus = b.distr_claytoncopula2_normal_mus;
  distr_claytoncopula2_normal_sigma2s = b.distr_claytoncopula2_normal_sigma2s;
  distr_copulas = b.distr_copulas;
  distr_tcopula_dfs = b.distr_tcopula_dfs;
  distr_tcopula_rhos = b.distr_tcopula_rhos;
  distr_gaussiancopula_rhos = b.distr_gaussiancopula_rhos;
  distr_gaussiancopula_rhofzs = b.distr_gaussiancopula_rhofzs;
  distr_frankcopula_rhos = b.distr_frankcopula_rhos;
  distr_frankcopula2_rhos = b.distr_frankcopula2_rhos;
  distr_frankcopula2_exp_rhos = b.distr_frankcopula2_exp_rhos;
  distr_frankcopula_exp_rhos = b.distr_frankcopula_exp_rhos;
  distr_frankcopula2_normal_mus = b.distr_frankcopula2_normal_mus;
  distr_frankcopula2_normal_sigma2s = b.distr_frankcopula2_normal_sigma2s;
  distr_sfa0_mu_ys = b.distr_sfa0_mu_ys;
  distr_sfa0_sigma_us = b.distr_sfa0_sigma_us;
  distr_sfa0_sigma_vs = b.distr_sfa0_sigma_vs;
  distr_sfa_mu_y_ids = b.distr_sfa_mu_y_ids;
  distr_sfa_mu_u_ids = b.distr_sfa_mu_u_ids;
  distr_sfa2_mu_ys = b.distr_sfa2_mu_ys;
  distr_sfa2_mu_us = b.distr_sfa2_mu_us;
  distr_sfa2_sigma_us = b.distr_sfa2_sigma_us;
  distr_sfa2_sigma_vs = b.distr_sfa2_sigma_vs;
  distr_sfa2_mu_y_ids = b.distr_sfa2_mu_y_ids;
  distr_sfa2_mu_u_ids = b.distr_sfa2_mu_u_ids;
  distr_sfa_mu_ys = b.distr_sfa_mu_ys;
  distr_sfa_mu_us = b.distr_sfa_mu_us;
  distr_sfa_sigma_us = b.distr_sfa_sigma_us;
  distr_sfa_sigma_vs = b.distr_sfa_sigma_vs;
  distr_sfa_alphas = b.distr_sfa_alphas;
  distr_gausscopulas = b.distr_gausscopulas;
  distr_gausscopula2s = b.distr_gausscopula2s;
  distr_clayton_copulas = b.distr_clayton_copulas;
  distr_gumbel_copulas = b.distr_gumbel_copulas;
  distr_binomialprobit_copulas = b.distr_binomialprobit_copulas;
  distr_JMs = b.distr_JMs;

  resultsyesno = b.resultsyesno;
  run_yes = b.run_yes;
  posteriormode = b.posteriormode;
  computemodeforstartingvalues = b.computemodeforstartingvalues;

  FC_linears = b.FC_linears;
  design_psplines = b.design_psplines;
  design_userdefineds = b.design_userdefineds;
  design_userdefined_tensors = b.design_userdefined_tensors;
  FC_tensor_omegas = b.FC_tensor_omegas;
  FC_nonps = b.FC_nonps;
  FC_merrors = b.FC_merrors;
  FC_nonp_variances = b.FC_nonp_variances;
  FC_nonp_variance_varselections = b.FC_nonp_variance_varselections;
  FC_varselection_omegas = b.FC_varselection_omegas;

  FC_predicts = b.FC_predicts;
  FC_predicts_mult = b.FC_predicts_mult;
  FC_predict_predictors = b.FC_predict_predictors;
  FC_predictive_checks = b.FC_predictive_checks;

  ridge = b.ridge;
  lasso = b.lasso;
  ssvs  = b.ssvs;
  ridge_linear = b.ridge_linear;
  lasso_linear = b.lasso_linear;
  ssvs_linear = b.ssvs_linear;
  FC_variance_pen_vectors = b.FC_variance_pen_vectors;
  FC_variance_pen_vector_ssvss = b.FC_variance_pen_vector_ssvss;
  FC_linear_pens = b.FC_linear_pens;

  design_hrandoms = b.design_hrandoms;
  FC_hrandom_variances = b.FC_hrandom_variances;
  FC_hrandom_variance_vecs = b.FC_hrandom_variance_vecs;
  FC_hrandom_variance_vec_nmigs = b.FC_hrandom_variance_vec_nmigs;
  FC_hrandom_variance_ssvss = b.FC_hrandom_variance_ssvss;

  design_krigings = b.design_krigings;
  design_mrfs = b.design_mrfs;

  errors = b.errors;

  }


const superbayesreg & superbayesreg::operator=(const superbayesreg & b)
  {
  if (this == & b)
	 return *this;
  statobject::operator=(statobject(b));
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = b.adminp_p;
  #endif

  firstvarselection = b.firstvarselection;

  pathres = b.pathres;
  title = b.title;
  pathnonp = b.pathnonp;

  statobj = b.statobj;

  D = b.D;

  modelvarnamesv = b.modelvarnamesv;

  equations=b.equations;
  mainequation=b.mainequation;
  countmarginal=b.countmarginal;
  nrlevel1 = b.nrlevel1;
  simobj = b.simobj;

  shared = b.shared;
  master = b.master;

  generaloptions = b.generaloptions;

  distr_gaussians = b.distr_gaussians;
  distr_gaussian_multeffects = b.distr_gaussian_multeffects;
  distr_hetgaussians = b.distr_hetgaussians;
  distr_vargaussians = b.distr_vargaussians;
  distr_quantregs = b.distr_quantregs;
  distr_gaussianmixtures = b.distr_gaussianmixtures;
  distr_loggaussians = b.distr_loggaussians;
  distr_gaussian_res = b.distr_gaussian_res;
  distr_gaussian_exps = b.distr_gaussian_exps;
  distr_gaussian_mults = b.distr_gaussian_mults;
  distr_binomials = b.distr_binomials;
  distr_cloglogs = b.distr_cloglogs;
  distr_poissons = b.distr_poissons;
  distr_poisson_exts = b.distr_poisson_exts;
  distr_poisson_extlins = b.distr_poisson_extlins;
  distr_binomialprobits = b.distr_binomialprobits;
  distr_binomialsvms = b.distr_binomialsvms;
  distr_multinomprobits = b.distr_multinomprobits;
  distr_multgaussians = b.distr_multgaussians;
  distr_multinomlogits = b.distr_multinomlogits;
  distr_logit_fruehwirths = b.distr_logit_fruehwirths;
  distr_ziplambdas = b.distr_ziplambdas;
  distr_zippis = b.distr_zippis;
  distr_hurdle_lambdas = b.distr_hurdle_lambdas;
  distr_hurdle_pis = b.distr_hurdle_pis;
  distr_hurdle_mus = b.distr_hurdle_mus;
  distr_hurdle_deltas = b.distr_hurdle_deltas;
  distr_negbinzip_mus = b.distr_negbinzip_mus;
  distr_negbinzip_pis = b.distr_negbinzip_pis;
  distr_negbinzip_deltas = b.distr_negbinzip_deltas;
  distr_zip_cloglog_mus = b.distr_zip_cloglog_mus;
  distr_zip_cloglog_pis = b.distr_zip_cloglog_pis;
  distr_negbin_mus = b.distr_negbin_mus;
  distr_negbin_deltas = b.distr_negbin_deltas;
  distr_beta_mus = b.distr_beta_mus;
  distr_beta_sigma2s = b.distr_beta_sigma2s;
  distr_lognormal_mus = b.distr_lognormal_mus;
  distr_lognormal_sigma2s = b.distr_lognormal_sigma2s;
  distr_lognormal2_mus = b.distr_lognormal2_mus;
  distr_lognormal2_sigmas = b.distr_lognormal2_sigmas;
  distr_normal_mus = b.distr_normal_mus;
  distr_normal_sigma2s = b.distr_normal_sigma2s;
  distr_cnormal_mus = b.distr_cnormal_mus;
  distr_cnormal_sigmas = b.distr_cnormal_sigmas;
  distr_normal2_mus = b.distr_normal2_mus;
  distr_truncnormal2_mus = b.distr_truncnormal2_mus;
  distr_truncnormal2_sigmas = b.distr_truncnormal2_sigmas;
  distr_normal2_sigmas = b.distr_normal2_sigmas;
  distr_invgaussian_mus = b.distr_invgaussian_mus;
  distr_invgaussian_sigma2s = b.distr_invgaussian_sigma2s;
  distr_gamma_mus = b.distr_gamma_mus;
  distr_gamma_sigmas = b.distr_gamma_sigmas;
  distr_pareto_bs = b.distr_pareto_bs;
  distr_pareto_ps = b.distr_pareto_ps;
  distr_gengamma_mus = b.distr_gengamma_mus;
  distr_gengamma_sigmas = b.distr_gengamma_sigmas;
  distr_gengamma_taus = b.distr_gengamma_taus;
  distr_t_mus = b.distr_t_mus;
  distr_t_sigma2s = b.distr_t_sigma2s;
  distr_t_dfs = b.distr_t_dfs;
  distr_BCCG_mus = b.distr_BCCG_mus;
  distr_BCCG_sigmas = b.distr_BCCG_sigmas;
  distr_BCCG_nus = b.distr_BCCG_nus;
  distr_bivnormal_mus = b.distr_bivnormal_mus;
  distr_bivnormal_sigmas = b.distr_bivnormal_sigmas;
  distr_bivnormal_rhos = b.distr_bivnormal_rhos;
  distr_bivnormal_rhofzs = b.distr_bivnormal_rhofzs;
  distr_bivnormal_mufzs = b.distr_bivnormal_mufzs;
  distr_bivt_mus = b.distr_bivt_mus;
  distr_bivt_sigmas = b.distr_bivt_sigmas;
  distr_bivt_dfs = b.distr_bivt_dfs;
  distr_bivt_rhos = b.distr_bivt_rhos;
  distr_bivprobit_mus = b.distr_bivprobit_mus;
  distr_bivprobit_rhos = b.distr_bivprobit_rhos;
  distr_bivprobit2_mus = b.distr_bivprobit2_mus;
  distr_bivprobit2_rhos = b.distr_bivprobit2_rhos;
  distr_bivlogit_mus = b.distr_bivlogit_mus;
  distr_bivlogit_ors = b.distr_bivlogit_ors;
  distr_zeroadjusteds = b.distr_zeroadjusteds;
  distr_dirichlets = b.distr_dirichlets;
  distr_zeroadjusted_mults = b.distr_zeroadjusted_mults;
  distr_weibull_lambdas = b.distr_weibull_lambdas;
  distr_weibull_alphas = b.distr_weibull_alphas;
  distr_weibull2_lambdas = b.distr_weibull2_lambdas;
  distr_weibull2_alphas = b.distr_weibull2_alphas;
  distr_dagum_as = b.distr_dagum_as;
  distr_dagum_bs = b.distr_dagum_bs;
  distr_dagum_ps = b.distr_dagum_ps;
  distr_betainf_mus = b.distr_betainf_mus;
  distr_betainf_sigma2s = b.distr_betainf_sigma2s;
  distr_betainf_nus = b.distr_betainf_nus;
  distr_betainf_taus = b.distr_betainf_taus;
  distr_betainf0_nus = b.distr_betainf0_nus;
  distr_betainf1_taus = b.distr_betainf1_taus;
  distr_gumbel_mus = b.distr_gumbel_mus;
  distr_gumbel_sigmas = b.distr_gumbel_sigmas;
  distr_gumbel2_mus = b.distr_gumbel2_mus;
  distr_gumbel2_sigma2s = b.distr_gumbel2_sigma2s;
  distr_gumbelcopula_rhos = b.distr_gumbelcopula_rhos;
  distr_gumbelcopula2_rhos = b.distr_gumbelcopula2_rhos;
  distr_gumbelcopula2_normal_mus = b.distr_gumbelcopula2_normal_mus;
  distr_gumbelcopula2_normal_sigma2s = b.distr_gumbelcopula2_normal_sigma2s;
  distr_gumbelcopula2_normal_mu_2s = b.distr_gumbelcopula2_normal_mu_2s;
  distr_gumbelcopula2_normal_sigma2_2s = b.distr_gumbelcopula2_normal_sigma2_2s;
  distr_gaussiancopula_binary_normal_rhos = b.distr_gaussiancopula_binary_normal_rhos;
  distr_gaussiancopula_binary_normal_sigmas = b.distr_gaussiancopula_binary_normal_sigmas;
  distr_gaussiancopula_binary_normal_mus = b.distr_gaussiancopula_binary_normal_mus;
  distr_gaussiancopula_binary_normal_latents = b.distr_gaussiancopula_binary_normal_latents;
  distr_gaussiancopula_binary_dagum_rhos = b.distr_gaussiancopula_binary_dagum_rhos;
  distr_gaussiancopula_binary_dagum_as = b.distr_gaussiancopula_binary_dagum_as;
  distr_gaussiancopula_binary_dagum_bs = b.distr_gaussiancopula_binary_dagum_bs;
  distr_gaussiancopula_binary_dagum_ps = b.distr_gaussiancopula_binary_dagum_ps;
  distr_gaussiancopula_binary_dagum_latents = b.distr_gaussiancopula_binary_dagum_latents;
  distr_gaussiancopula_dagum_rhos = b.distr_gaussiancopula_dagum_rhos;
  distr_gaussiancopula_dagum_as = b.distr_gaussiancopula_dagum_as;
  distr_gaussiancopula_dagum_bs = b.distr_gaussiancopula_dagum_bs;
  distr_gaussiancopula_dagum_ps = b.distr_gaussiancopula_dagum_ps;
  distr_claytoncopula_rhos = b.distr_claytoncopula_rhos;
  distr_claytoncopula2_rhos = b.distr_claytoncopula2_rhos;
  distr_claytoncopula2_normal_mus = b.distr_claytoncopula2_normal_mus;
  distr_claytoncopula2_normal_sigma2s = b.distr_claytoncopula2_normal_sigma2s;
  distr_copulas = b.distr_copulas;
  distr_tcopula_dfs = b.distr_tcopula_dfs;
  distr_tcopula_rhos = b.distr_tcopula_rhos;
  distr_gaussiancopula_rhos = b.distr_gaussiancopula_rhos;
  distr_gaussiancopula_rhofzs = b.distr_gaussiancopula_rhofzs;
  distr_frankcopula_rhos = b.distr_frankcopula_rhos;
  distr_frankcopula2_rhos = b.distr_frankcopula2_rhos;
  distr_frankcopula_exp_rhos = b.distr_frankcopula_exp_rhos;
  distr_frankcopula2_exp_rhos = b.distr_frankcopula2_exp_rhos;
  distr_frankcopula2_normal_mus = b.distr_frankcopula2_normal_mus;
  distr_frankcopula2_normal_sigma2s = b.distr_frankcopula2_normal_sigma2s;
  distr_sfa0_mu_ys = b.distr_sfa0_mu_ys;
  distr_sfa0_sigma_us = b.distr_sfa0_sigma_us;
  distr_sfa0_sigma_vs = b.distr_sfa0_sigma_vs;
  distr_sfa2_mu_y_ids = b.distr_sfa2_mu_y_ids;
  distr_sfa2_mu_u_ids = b.distr_sfa2_mu_u_ids;
  distr_sfa2_mu_ys = b.distr_sfa2_mu_ys;
  distr_sfa2_mu_us = b.distr_sfa2_mu_us;
  distr_sfa_mu_y_ids = b.distr_sfa_mu_y_ids;
  distr_sfa_mu_u_ids = b.distr_sfa_mu_u_ids;
  distr_sfa2_sigma_us = b.distr_sfa2_sigma_us;
  distr_sfa2_sigma_vs = b.distr_sfa2_sigma_vs;
  distr_sfa_mu_ys = b.distr_sfa_mu_ys;
  distr_sfa_mu_us = b.distr_sfa_mu_us;
  distr_sfa_sigma_us = b.distr_sfa_sigma_us;
  distr_sfa_sigma_vs = b.distr_sfa_sigma_vs;
  distr_sfa_alphas = b.distr_sfa_alphas;
  distr_gausscopulas = b.distr_gausscopulas;
  distr_gausscopula2s = b.distr_gausscopula2s;
  distr_clayton_copulas = b.distr_clayton_copulas;
  distr_gumbel_copulas = b.distr_gumbel_copulas;
  distr_binomialprobit_copulas = b.distr_binomialprobit_copulas;
  distr_JMs = b.distr_JMs;

  resultsyesno = b.resultsyesno;
  run_yes = b.run_yes;
  posteriormode = b.posteriormode;
  computemodeforstartingvalues = b.computemodeforstartingvalues;

  FC_linears = b.FC_linears;
  design_psplines = b.design_psplines;
  design_userdefineds = b.design_userdefineds;
  design_userdefined_tensors = b.design_userdefined_tensors;
  FC_tensor_omegas = b.FC_tensor_omegas;
  FC_nonps = b.FC_nonps;
  FC_merrors = b.FC_merrors;
  FC_nonp_variances = b.FC_nonp_variances;
  FC_nonp_variance_varselections = b.FC_nonp_variance_varselections;
  FC_varselection_omegas = b.FC_varselection_omegas;
  FC_predicts = b.FC_predicts;
  FC_predicts_mult = b.FC_predicts_mult;
  FC_predict_predictors = b.FC_predict_predictors;
  FC_predictive_checks = b.FC_predictive_checks;

  ridge = b.ridge;
  lasso = b.lasso;
  ssvs = b.ssvs;
  ridge_linear = b.ridge_linear;
  lasso_linear = b.lasso_linear;
  ssvs_linear = b.ssvs_linear;
  FC_variance_pen_vectors = b.FC_variance_pen_vectors;
  FC_variance_pen_vector_ssvss = b.FC_variance_pen_vector_ssvss;
  FC_linear_pens = b.FC_linear_pens;

  design_hrandoms = b.design_hrandoms;
  FC_hrandom_variances = b.FC_hrandom_variances;
  FC_hrandom_variance_vecs = b.FC_hrandom_variance_vecs;
  FC_hrandom_variance_vec_nmigs = b.FC_hrandom_variance_vec_nmigs;
  FC_hrandom_variance_ssvss = b.FC_hrandom_variance_ssvss;

  design_krigings = b.design_krigings;
  design_mrfs = b.design_mrfs;

  errors = b.errors;

  return *this;
  }


int superbayesreg::parse(const ST::string & c)
  {

  int u = 0;
  u = statobject::parse(c);

  int pos = statobject::parsecom(c,methods,globaloptions);

  if (pos >= 0)
	 (*functions[pos])(*this);

  return pos;
  }


bool superbayesreg::create_generaloptions(void)
  {
  if (equations.size()==1)
    {
    if (iterations.getvalue()- burnin.getvalue() < 1)
      {
      outerror("ERROR: number of iterations must exceed number of burnin iterations\n");
      return true;
      }

    if (step.getvalue() > iterations.getvalue() - burnin.getvalue())
      {
      outerror("ERROR: thinning parameter too large\n");
      return true;
      }


    generaloptions = MCMC::GENERAL_OPTIONS(
    #if defined(JAVA_OUTPUT_WINDOW)
    adminb_p,
    #endif
    iterations.getvalue(),burnin.getvalue(),
                                step.getvalue(),saveestimation.getvalue(),
                                copula.getvalue(),rotation.getvalue(),
                                samplesel.getvalue(), sampleselval.getvalue(),
                                forceIWLS.getvalue(),
                                logout,
                                level1.getvalue(),level2.getvalue());

    describetext.push_back("ESTIMATION OPTIONS:\n");
    describetext.push_back("\n");
    describetext.push_back("Number of Iterations: "
                              + ST::inttostring(iterations.getvalue()) + "\n");
    describetext.push_back("Burnin: " + ST::inttostring(burnin.getvalue()) + "\n");
    describetext.push_back("Thinning parameter: " +
                              ST::inttostring(step.getvalue()) + "\n");

  //  generaloptions.nrout(iterationsprint.getvalue());

    return false;
   }
  else
    return false;

  }


void superbayesreg::make_header(unsigned & modnr)
  {
  if (equations[modnr].hlevel == 1)
    {

    ST::string rn = equations[modnr].distrp->responsename;

    equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": MAIN " + equations[modnr].equationtype + " REGRESSION_"+ rn;
    equations[modnr].paths = "MAIN_" + equations[modnr].equationtype + "_REGRESSION_"+ rn;

    }
  else if (equations[modnr].hlevel == 2)
    {

    equations[modnr].header = "MCMCREG OBJECT " + name.to_bstr() +
                              ": RANDOM EFFECTS " + equations[modnr].equationtype + " REGRESSION";
    equations[modnr].paths = "RANDOM_EFFECTS_" + equations[modnr].equationtype;

    }

  }


void hregressrun(superbayesreg & b)
  {

  if ((b.equations.size()==0) || (b.run_yes==true))
    {
    b.clear();
    b.run_yes=false;
    b.generaloptions_yes = false;
    b.errors=false;
    b.nrlevel1 = 0;
    b.countmarginal = 0;
    b.mainequation=false;
    }

  if (b.errors==false)    // errors from level 2 equations
    {
    b.resultsyesno = false;
    if (b.modeonly.getvalue() == true)
      b.posteriormode = true;
    else
      b.posteriormode = false;

    b.terms = b.modreg.getterms();

    b.describetext.erase(b.describetext.begin(),b.describetext.end());
    b.describetext.push_back("LAST ESTIMATED MODEL: \n");
    b.describetext.push_back("\n");
    b.describetext.push_back(b.modreg.getModelText());
    b.describetext.push_back("\n");

    unsigned modnr = b.equations.size();
    b.equations.push_back(equation(modnr+1,b.hlevel.getvalue(),
                                   b.equationtype.getvalue()));

    b.firstvarselection = true;
    bool failure = false;

    if (!failure && b.generaloptions_yes==false)
      {
      failure = b.create_generaloptions();
      b.generaloptions_yes = true;
      }

    if (!failure)
      failure = b.create_distribution();

    if (!failure)
      b.make_header(modnr);


    if (!failure)
      failure = b.create_linear();
    if (!failure && b.terms.size() >= 1)
      failure = b.create_nonp();

    if (!failure)
      failure = b.create_predict();

    if (!failure)
      b.create_predictive_check();

    if (!failure)
      b.create_cv();

    if (b.hlevel.getvalue() == 1)
      b.nrlevel1++;

    if ((! failure) && (b.hlevel.getvalue() == 1) &&
        (b.mainequation==true))
      {

      if (!failure)
        failure = b.check_errors();

      if (!failure)
        {
        b.run_yes=true;

        unsigned mit = b.modemaxit.getvalue();
        b.simobj = MCMCsim(&b.generaloptions,b.equations,mit);

        ST::string pathgraphs = b.outfile.getvalue();

        if (b.modeonly.getvalue())
          {
          failure = b.simobj.posteriormode(pathgraphs,b.shared,false);
          }
        else
          failure = b.simobj.simulate(pathgraphs,b.setseed.getvalue(),
          b.computemodeforstartingvalues, b.shared);
        }

      if (!failure)
        b.resultsyesno = true;
      }


    if (failure)
      b.errors=true;

    } // end: if errors == false

  }


bool superbayesreg::check_errors(void)
  {

  bool err=false;

  unsigned j,k;

  for (k=0;k<equations.size();k++)
    {
//    bool test = equations[k].distrp->errors;
//    int test2 = k;
    if (equations[k].distrp->errors == true)
      {
      err = true;
      for (j=0;j<equations[k].distrp->errormessages.size();j++)
        {
//        ST::string test3 = equations[k].distrp->errormessages[j];
        outerror(equations[k].distrp->errormessages[j]);
        }
      }

    }


/*
  for (k=0;k<equations.size();k++)
    {
    for (i=0;i<equations[k].distrp.size();i++)
      {
      if (equations[k].distrp[i]->errors == true)
        {
        err = true;
        for (j=0;j<equations[k].distrp[i]->errormessages.size();j++)
          outerror(equations[k].distrp[i]->errormessages[j]);
        }

      }

    }
*/

  return err;
  }


bool superbayesreg::create_distribution(void)
  {

  unsigned i;
  bool failure=false;


  //--------------------- reading dataset information --------------------------

  dataobject * datap;               // pointer to dataset

  int objpos = findstatobject(*statobj,udata.getusingtext(),"dataset");
  statobject * s;
  if (objpos >= 0)
    {
    s = statobj->at(objpos);
    datap = dynamic_cast<dataobject*>(s);
    }
  else
    {
    if (objpos == -1)
      outerror("ERROR: " + udata.getusingtext() + " is not existing\n");
    else
      outerror("ERROR: " + udata.getusingtext() + " is not a dataset object\n");
    return true;
    }

  //------------------ end: reading dataset information ------------------------


  //---------------- reading data, creating designmatrices ---------------------

  ST::string rname;
  ST::string wn;
  ST::string ifexpression;
  unsigned weightpos=0;

  vector<ST::string> modelvarnamesh;
  modelvarnamesh = modreg.getModelVarnamesAsVector();
  modelvarnamesv.erase(modelvarnamesv.begin(),modelvarnamesv.end());
  for(i=0;i<modelvarnamesh.size();i++)
    if (modelvarnamesh[i] != "const")
      modelvarnamesv.push_back(modelvarnamesh[i]);

  rname = modelvarnamesv[0].to_bstr();
  wn = methods[0].get_weight_variable().to_bstr();

  if (wn.length() != 0)
    {
    modelvarnamesv.push_back(wn);
    weightpos = modelvarnamesv.size()-1;
    }

  ifexpression = methods[0].getexpression();

  // testing, wether all variables specified are already existing
  vector<ST::string> notex;
  if ((datap->allexisting(modelvarnamesv,notex)) == false)
    {
    for (i=0;i<notex.size();i++)
      {
      ST::string test = notex[i];
      if (notex[i] != "const" && notex[i] != "noconst")
        {
        outerror("ERROR: variable " + notex[i] + " is not existing\n");
        failure = true;
        }
      }
    if (failure)
      return true;
    } // end: if ((datap->allexisting(modelvarnamesv,notex)) == false)


  datap->makematrix(modelvarnamesv,D,ifexpression);

  errormessages = datap->geterrormessages();

  if (!errormessages.empty())
    return true;

//  datamatrix offs(1,1,0);
//  failure = create_offset(offs);


  datamatrix w;
  bool weightsdefined = false;

  if (wn.length() > 0)
    {
    w = D.getCol(weightpos);
    weightsdefined=true;
    }
  else
    w = datamatrix(1,1);



  describetext.push_back("Response distribution: "
                           + family.getvalue() + "\n");


  unsigned modnr = equations.size()-1;

  // ST::string f = family.getvalue();

//---------------------------- Gaussian response -------------------------------
  if (((family.getvalue() == "gaussian") ||  (family.getvalue() == "normal")) && (equationtype.getvalue()=="mu") &&
      (distr_vargaussians.size()==0) && (distr_normal_sigma2s.size()==0)
     )
    {
    computemodeforstartingvalues = true;

    mainequation=true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_gaussians.push_back(DISTR_gaussian(aresp.getvalue(),bresp.getvalue(),
                                      &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussians[distr_gaussians.size()-1];
    equations[modnr].pathd = outfile.getvalue() + "_scale.res";

    if(shared)
      {
      mainequation=false;
      distr_JMs[distr_JMs.size()-1].dist2 = &distr_gaussians[distr_gaussians.size()-1];
      }
    }
//-------------------------- END: Gaussian response ----------------------------

//------------------------- Gaussian multiplicative response -------------------
  else if ((family.getvalue() == "gaussian_multeffect") && (equationtype.getvalue()=="mu"))
    {
    computemodeforstartingvalues = true;
    mainequation=false;
    distr_gaussian_multeffects.push_back(DISTR_gaussian_multeffect(aresp.getvalue(),bresp.getvalue(),
                                      &generaloptions,D.getCol(0),"",w) );
    equations[modnr].distrp = &distr_gaussian_multeffects[distr_gaussian_multeffects.size()-1];
    }
//-------------------------- END: Gaussian response ----------------------------


//-------------------------------- normal sigma2 ---------------------------------
  else if ( ((family.getvalue() == "normal")) && equationtype.getvalue()=="sigma2")
    {

    mainequation=false;

    computemodeforstartingvalues = true;

    distr_normal_sigma2s.push_back(DISTR_normal_sigma2(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_normal_sigma2s[distr_normal_sigma2s.size()-1];
    equations[modnr].pathd = "";

    }
//---------------------------- END: normal sigma2 -------------------------------

//------------------------------- normal mu ------------------------------------
  else if ((family.getvalue() == "normal") && (equationtype.getvalue()=="mu") && (distr_normal_sigma2s.size()>=1))

    {
    if(generaloptions.copula)
      {
      mainequation=false;
      }
    else
      {
      mainequation=true;
      }


    computemodeforstartingvalues = true;

    distr_normal_mus.push_back(DISTR_normal_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_normal_mus[distr_normal_mus.size()-1];
    equations[modnr].pathd = "";

    if (distr_normal_sigma2s.size() != distr_normal_mus.size())
      {
      outerror("ERROR: Equation for sigma2 is missing");
      return true;
      }
    else
      {
      distr_normal_sigma2s[distr_normal_sigma2s.size()-1].distrp.push_back
      (&distr_normal_mus[distr_normal_mus.size()-1]);

      distr_normal_mus[distr_normal_mus.size()-1].distrp.push_back
      (&distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);

      countmarginal += 1;
      }
    if(generaloptions.copula){}
    else
      {
      predict_mult_distrs.push_back(&distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);
      predict_mult_distrs.push_back(&distr_normal_mus[distr_normal_mus.size()-1]);
      }

    if(generaloptions.copula)
      {
      if (countmarginal == 1)
        {
        distr_normal_mus[distr_normal_mus.size()-1].set_copulapos(0);
        distr_normal_sigma2s[distr_normal_sigma2s.size()-1].set_copulapos(0);
        }
      else if(countmarginal == 2)
        {
        distr_normal_mus[distr_normal_mus.size()-1].set_copulapos(1);
        distr_normal_sigma2s[distr_normal_sigma2s.size()-1].set_copulapos(1);
        }
      else
        {
        outerror("ERROR: currently only bivariate copula models implemented");
        return true;
        }
      }
    }
//------------------------------- END: normal mu -------------------------------


//------------- variance equation of heteroscedastic Gaussian ------------------
  else if (family.getvalue() == "hetgaussian" && equationtype.getvalue()=="sigma2")
    {

    distr_vargaussians.push_back(DISTR_vargaussian(&generaloptions,D.getCol(0)));

    equations[modnr].distrp = &distr_vargaussians[distr_vargaussians.size()-1];
    equations[modnr].pathd = "";

    }
//---------- END: variance equation of heteroscedastic Gaussian ----------------

//---------------------- heteroscedastic Gaussian response ---------------------
  else if ((family.getvalue() == "hetgaussian") && (equationtype.getvalue()=="mu") && (distr_vargaussians.size()==1))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_hetgaussians.push_back(DISTR_hetgaussian(aresp.getvalue(),
                                 bresp.getvalue(), &generaloptions,
                                 D.getCol(0),path,scaleconst.getvalue(),w) );

    equations[modnr].distrp = &distr_hetgaussians[distr_hetgaussians.size()-1];

#if defined(__BUILDING_LINUX)
    equations[modnr].pathd = defaultpath + "/temp/" + name  + "_scale.res";
#else
    equations[modnr].pathd = defaultpath + "\\temp\\" + name  + "_scale.res";
#endif

    distr_vargaussians[0].dgaussian = &distr_hetgaussians[distr_hetgaussians.size()-1];

    }
//----------------- END: heteroscedastic Gaussian response ---------------------


//-------------------------------- lognormal sigma2 ---------------------------------
  else if ((family.getvalue() == "lognormal" || family.getvalue() == "za_lognormal")
           && equationtype.getvalue()=="sigma2")
    {

    computemodeforstartingvalues = true;

    datamatrix dnew = D.getCol(0);
    if(family.getvalue() == "za_lognormal")
      {
      if(w.rows()==1)
        w = datamatrix(dnew.rows(),1,1);
      for(unsigned i=0; i<dnew.rows(); i++)
        {
        if(D(i,0)==0)
          {
          w(i,0) = 0;
          dnew(i,0) = 0.1;
          }
        }
      }

    distr_lognormal_sigma2s.push_back(DISTR_lognormal_sigma2(&generaloptions,dnew,w));

    equations[modnr].distrp = &distr_lognormal_sigma2s[distr_lognormal_sigma2s.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_lognormal_sigma2s[distr_lognormal_sigma2s.size()-1]);

    }
//---------------------------- END: lognormal sigma2 -------------------------------

//------------------------------- lognormal mu ------------------------------------
  else if ((family.getvalue() == "lognormal" || family.getvalue() == "za_lognormal")
           && ((equationtype.getvalue()=="mu"))
          )
    {
    datamatrix dnew= D.getCol(0);

    if(family.getvalue() == "lognormal")
      mainequation=true;
    else
      {
      mainequation=false;
      if(w.rows()==1)
        w = datamatrix(dnew.rows(),1,1);
      for(unsigned i=0; i<dnew.rows(); i++)
        {
        if(D(i,0)==0)
          {
          w(i,0) = 0;
          dnew(i,0) = 0.1;
          }
        }
      }
    computemodeforstartingvalues = true;

    distr_lognormal_mus.push_back(DISTR_lognormal_mu(&generaloptions,dnew,w));

    equations[modnr].distrp = &distr_lognormal_mus[distr_lognormal_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_lognormal_mus[distr_lognormal_mus.size()-1]);

    if (distr_lognormal_sigma2s.size() != 1)
      {
      outerror("ERROR: Equation for sigma2 is missing");
      return true;
      }
    else
      {
      distr_lognormal_sigma2s[distr_lognormal_sigma2s.size()-1].distrp.push_back
      (&distr_lognormal_mus[distr_lognormal_mus.size()-1]);

      distr_lognormal_mus[distr_lognormal_mus.size()-1].distrp.push_back
      (&distr_lognormal_sigma2s[distr_lognormal_sigma2s.size()-1]);
      }

    }
//------------------------------- END: lognormal mu -------------------------------



//-------------------------------- lognormal2 sigma ---------------------------------
  else if (family.getvalue() == "lognormal2" && equationtype.getvalue()=="sigma")
    {

    computemodeforstartingvalues = true;

    distr_lognormal2_sigmas.push_back(DISTR_lognormal2_sigma(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_lognormal2_sigmas[distr_lognormal2_sigmas.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_lognormal2_sigmas[distr_lognormal2_sigmas.size()-1]);

    }
//---------------------------- END: lognormal2 sigma -------------------------------

//------------------------------- lognormal2 mu ------------------------------------
  else if ((family.getvalue() == "lognormal2") && (equationtype.getvalue()=="mu"))

    {

     mainequation=true;

    computemodeforstartingvalues = true;

    distr_lognormal2_mus.push_back(DISTR_lognormal2_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_lognormal2_mus[distr_lognormal2_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_lognormal2_mus[distr_lognormal2_mus.size()-1]);

    if (distr_lognormal2_sigmas.size() != 1)
      {
      outerror("ERROR: Equation for sigma is missing");
      return true;
      }
    else
      {
      distr_lognormal2_sigmas[distr_lognormal2_sigmas.size()-1].distrp.push_back
      (&distr_lognormal2_mus[distr_lognormal2_mus.size()-1]);

      distr_lognormal2_mus[distr_lognormal2_mus.size()-1].distrp.push_back
      (&distr_lognormal2_sigmas[distr_lognormal2_sigmas.size()-1]);
      }

    }
//------------------------------- END: lognormal2 mu -------------------------------


//-------------------------------- cnormal sigma ---------------------------------
  else if (family.getvalue() == "cnormal" && equationtype.getvalue()=="sigma")
    {

    computemodeforstartingvalues = true;

    distr_cnormal_sigmas.push_back(DISTR_cnormal_sigma(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_cnormal_sigmas[distr_cnormal_sigmas.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_cnormal_sigmas[distr_cnormal_sigmas.size()-1]);

    }
//---------------------------- END: cnormal sigma -------------------------------

//------------------------------- cnormal mu ------------------------------------
  else if ((family.getvalue() == "cnormal") && (equationtype.getvalue()=="mu"))

    {

     mainequation=true;

    computemodeforstartingvalues = true;

    distr_cnormal_mus.push_back(DISTR_cnormal_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_cnormal_mus[distr_cnormal_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_cnormal_mus[distr_cnormal_mus.size()-1]);

    if (distr_cnormal_sigmas.size() != 1)
      {
      outerror("ERROR: Equation for sigma is missing");
      return true;
      }
    else
      {
      distr_cnormal_sigmas[distr_cnormal_sigmas.size()-1].distrp.push_back
      (&distr_cnormal_mus[distr_cnormal_mus.size()-1]);

      distr_cnormal_mus[distr_cnormal_mus.size()-1].distrp.push_back
      (&distr_cnormal_sigmas[distr_cnormal_sigmas.size()-1]);
      }
    }
//------------------------------- END: cnormal mu -------------------------------



//---------------------------------- ZINB pi -----------------------------------
  else if (family.getvalue() == "zinb" && equationtype.getvalue()=="pi")
    {

    computemodeforstartingvalues = true;

    distr_negbinzip_pis.push_back(DISTR_negbinzip_pi(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_negbinzip_pis[distr_negbinzip_pis.size()-1];
    equations[modnr].pathd = "";

    }
//------------------------------- END: ZINB pi ---------------------------------

//---------------------------------- ZINB delta --------------------------------
  else if (family.getvalue() == "zinb" && equationtype.getvalue()=="delta")
    {

    computemodeforstartingvalues = true;

    int strmax = stoprmax.getvalue();
    bool sl = slow.getvalue();
    int nb = nrbetween.getvalue();
    double sts = stopsum.getvalue();

    distr_negbinzip_deltas.push_back(DISTR_negbinzip_delta(&generaloptions,
                                     D.getCol(0),sts,strmax,nb,sl,w));

    equations[modnr].distrp = &distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1];
    equations[modnr].pathd = "";

    }
//------------------------------ END: ZINB delta -------------------------------


//----------------------------------- ZINB mu ----------------------------------
  else if (family.getvalue() == "zinb" && equationtype.getvalue()=="mu")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_negbinzip_mus.push_back(DISTR_negbinzip_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_negbinzip_mus[distr_negbinzip_mus.size()-1];
    equations[modnr].pathd = "";

    if (distr_negbinzip_pis.size() != 1)
      {
      outerror("ERROR: Equation for pi is missing");
      return true;
      }

    if (distr_negbinzip_deltas.size() != 1)
      {
      outerror("ERROR: Equation for delta is missing");
      return true;
      }

    predict_mult_distrs.push_back(&distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1]);
    predict_mult_distrs.push_back(&distr_negbinzip_pis[distr_negbinzip_pis.size()-1]);
    predict_mult_distrs.push_back(&distr_negbinzip_mus[distr_negbinzip_mus.size()-1]);

    distr_negbinzip_pis[distr_negbinzip_pis.size()-1].distrmu =
    &distr_negbinzip_mus[distr_negbinzip_mus.size()-1];

    distr_negbinzip_pis[distr_negbinzip_pis.size()-1].distrdelta =
    &distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1];


    distr_negbinzip_mus[distr_negbinzip_mus.size()-1].distrpi =
    &distr_negbinzip_pis[distr_negbinzip_pis.size()-1];

    distr_negbinzip_mus[distr_negbinzip_mus.size()-1].distrdelta =
    &distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1];


    distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1].distrpi =
    &distr_negbinzip_pis[distr_negbinzip_pis.size()-1];

    distr_negbinzip_deltas[distr_negbinzip_deltas.size()-1].distrmu =
    &distr_negbinzip_mus[distr_negbinzip_mus.size()-1];

    }
//------------------------------ END: ZINB mu ----------------------------------

//-------------------------------- negbin delta --------------------------------
  else if (family.getvalue() == "negbin" && equationtype.getvalue()=="delta")
    {

    computemodeforstartingvalues = true;

//    double flimit = fraclimit.getvalue();
    int strmax = stoprmax.getvalue();
    bool sl = slow.getvalue();
    int nb = nrbetween.getvalue();
    double sts = stopsum.getvalue();

    distr_negbin_deltas.push_back(DISTR_negbin_delta(&generaloptions,D.getCol(0),
    sts,strmax,nb,sl,w));

    equations[modnr].distrp = &distr_negbin_deltas[distr_negbin_deltas.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_negbin_deltas[distr_negbin_deltas.size()-1]);

    }
//---------------------------- END: negbin delta -------------------------------

//------------------------------- negbin mu ------------------------------------
  else if (family.getvalue() == "negbin" && equationtype.getvalue()=="mu")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_negbin_mus.push_back(DISTR_negbin_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_negbin_mus[distr_negbin_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_negbin_mus[distr_negbin_mus.size()-1]);

    if (distr_negbin_deltas.size() != 1)
      {
      outerror("ERROR: Equation for delta is missing");
      return true;
      }
    else
      {
      distr_negbin_deltas[distr_negbin_deltas.size()-1].distrp.push_back
      (&distr_negbin_mus[distr_negbin_mus.size()-1]);

      distr_negbin_mus[distr_negbin_mus.size()-1].distrp.push_back
      (&distr_negbin_deltas[distr_negbin_deltas.size()-1]);
      }

    }
//------------------------------- END: negbin mu -------------------------------

//-------------------------------- betainf_mu ---------------------------------
  else if (family.getvalue() == "betainf" && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;
    int i;
    if(w.rows()==1)
        w = datamatrix(D.rows(),1,1);
    double * wp = w.getV();
    double * rp = D.getV();
    for(i=0;i<w.rows();i++,wp++,rp+=D.cols())
      {
      if((*rp)==0)
        *wp = 0;
	  else if((*rp)==1)
	    *wp = 0;
      }
    distr_betainf_mus.push_back(DISTR_betainf_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_betainf_mus[distr_betainf_mus.size()-1];
    equations[modnr].pathd = "";


    }
//---------------------------- END: betainf_mu -------------------------------

//-------------------------------- betainf_sigma2 ---------------------------------
  else if (family.getvalue() == "betainf" && equationtype.getvalue()=="sigma2")
    {

    computemodeforstartingvalues = true;
    int i;
    if(w.rows()==1)
        w = datamatrix(D.rows(),1,1);
    double * wp = w.getV();
    double * rp = D.getV();
    for(i=0;i<w.rows();i++,wp++,rp+=D.cols())
      {
      if((*rp)==0)
        *wp = 0;
	  else if((*rp)==1)
	    *wp = 0;
      }
    distr_betainf_sigma2s.push_back(DISTR_betainf_sigma2(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_betainf_sigma2s[distr_betainf_sigma2s.size()-1];
    equations[modnr].pathd = "";

    }
//---------------------------- END: betainf_sigma2 -------------------------------

//-------------------------------- betainf_nu ---------------------------------
  else if (family.getvalue() == "betainf" && equationtype.getvalue()=="nu")
    {

    computemodeforstartingvalues = true;


    distr_betainf_nus.push_back(DISTR_betainf_nu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_betainf_nus[distr_betainf_nus.size()-1];
    equations[modnr].pathd = "";

    }
//---------------------------- END: betainf_nu -------------------------------

//------------------------------- betainf_tau ------------------------------------
  else if (family.getvalue() == "betainf" && equationtype.getvalue()=="tau")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_betainf_taus.push_back(DISTR_betainf_tau(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_betainf_taus[distr_betainf_taus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_betainf_mus[distr_betainf_mus.size()-1]);
    predict_mult_distrs.push_back(&distr_betainf_sigma2s[distr_betainf_sigma2s.size()-1]);
    predict_mult_distrs.push_back(&distr_betainf_nus[distr_betainf_nus.size()-1]);
    predict_mult_distrs.push_back(&distr_betainf_taus[distr_betainf_taus.size()-1]);


    if (distr_betainf_sigma2s.size() != 1)
      {
      outerror("ERROR: Equation for sigma2 is missing");
      return true;
      }

    if (distr_betainf_nus.size() != 1)
      {
      outerror("ERROR: Equation for nu is missing");
      return true;
      }

    if (distr_betainf_mus.size() != 1)
      {
      outerror("ERROR: Equation for mu is missing");
      return true;
      }

      distr_betainf_taus[distr_betainf_taus.size()-1].distrp.push_back
      (&distr_betainf_nus[distr_betainf_nus.size()-1]);

      distr_betainf_nus[distr_betainf_nus.size()-1].distrp.push_back
      (&distr_betainf_taus[distr_betainf_taus.size()-1]);

      distr_betainf_sigma2s[distr_betainf_sigma2s.size()-1].distrp.push_back
      (&distr_betainf_mus[distr_betainf_mus.size()-1]);

      distr_betainf_mus[distr_betainf_mus.size()-1].distrp.push_back
      (&distr_betainf_sigma2s[distr_betainf_sigma2s.size()-1]);


    }
//------------------------------- END: betainf_tau -------------------------------


//-------------------------------- betainf0_sigma2 ---------------------------------
  else if (family.getvalue() == "betainf0" && equationtype.getvalue()=="sigma2")
    {

    computemodeforstartingvalues = true;
    int i;
    if(w.rows()==1)
        w = datamatrix(D.rows(),1,1);
    double * wp = w.getV();
    double * rp = D.getV();
    for(i=0;i<w.rows();i++,wp++,rp+=D.cols())
      {
      if((*rp)==0)
        *wp = 0;
      }
    distr_beta_sigma2s.push_back(DISTR_beta_sigma2(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_beta_sigma2s[distr_beta_sigma2s.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_beta_sigma2s[distr_beta_sigma2s.size()-1]);

    }
//---------------------------- END: betainf0_sigma2 -------------------------------

//-------------------------------- betainf0_mu ---------------------------------
  else if (family.getvalue() == "betainf0" && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;
    int i;
    if(w.rows()==1)
        w = datamatrix(D.rows(),1,1);
    double * wp = w.getV();
    double * rp = D.getV();
    for(i=0;i<w.rows();i++,wp++,rp+=D.cols())
      {
      if((*rp)==0)
        *wp = 0;
      }
    distr_beta_mus.push_back(DISTR_beta_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_beta_mus[distr_beta_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_beta_mus[distr_beta_mus.size()-1]);

    }
//---------------------------- END: betainf_1mu -------------------------------

//------------------------------- betainf0_nu ------------------------------------
  else if (family.getvalue() == "betainf0" && equationtype.getvalue()=="nu")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_betainf0_nus.push_back(DISTR_betainf0_nu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_betainf0_nus[distr_betainf0_nus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_beta_mus[distr_beta_mus.size()-1]);
    predict_mult_distrs.push_back(&distr_beta_sigma2s[distr_beta_sigma2s.size()-1]);
    predict_mult_distrs.push_back(&distr_betainf0_nus[distr_betainf0_nus.size()-1]);


    if (distr_beta_sigma2s.size() != 1)
      {
      outerror("ERROR: Equation for sigma2 is missing");
      return true;
      }

    if (distr_beta_mus.size() != 1)
      {
      outerror("ERROR: Equation for mu is missing");
      return true;
      }

      distr_beta_sigma2s[distr_beta_sigma2s.size()-1].distrp.push_back
      (&distr_beta_mus[distr_beta_mus.size()-1]);

      distr_beta_mus[distr_beta_mus.size()-1].distrp.push_back
      (&distr_beta_sigma2s[distr_beta_sigma2s.size()-1]);


    }
//------------------------------- END: betainf0_nu -------------------------------

//-------------------------------- betainf1_sigma2 ---------------------------------
  else if (family.getvalue() == "betainf1" && equationtype.getvalue()=="sigma2")
    {

    computemodeforstartingvalues = true;
    int i;
    if(w.rows()==1)
        w = datamatrix(D.rows(),1,1);
    double * wp = w.getV();
    double * rp = D.getV();
    for(i=0;i<w.rows();i++,wp++,rp+=D.cols())
      {
	  if((*rp)==1)
	    *wp = 0;
      }
    distr_beta_sigma2s.push_back(DISTR_beta_sigma2(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_beta_sigma2s[distr_beta_sigma2s.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_beta_sigma2s[distr_beta_sigma2s.size()-1]);

    }
//---------------------------- END: betainf_1sigma2 -------------------------------

//-------------------------------- betainf1_mu ---------------------------------
  else if (family.getvalue() == "betainf1" && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;
    int i;
    if(w.rows()==1)
        w = datamatrix(D.rows(),1,1);
    double * wp = w.getV();
    double * rp = D.getV();
    for(i=0;i<w.rows();i++,wp++,rp+=D.cols())
      {
	  if((*rp)==1)
	    *wp = 0;
      }
    distr_beta_mus.push_back(DISTR_beta_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_beta_mus[distr_beta_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_beta_mus[distr_beta_mus.size()-1]);

    }
//---------------------------- END: betainf_1mu -------------------------------

//------------------------------- betainf1_tau ------------------------------------
  else if (family.getvalue() == "betainf1" && equationtype.getvalue()=="tau")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_betainf1_taus.push_back(DISTR_betainf1_tau(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_betainf1_taus[distr_betainf1_taus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_beta_mus[distr_beta_mus.size()-1]);
    predict_mult_distrs.push_back(&distr_beta_sigma2s[distr_beta_sigma2s.size()-1]);
    predict_mult_distrs.push_back(&distr_betainf1_taus[distr_betainf1_taus.size()-1]);


    if (distr_beta_sigma2s.size() != 1)
      {
      outerror("ERROR: Equation for sigma2 is missing");
      return true;
      }

    if (distr_beta_mus.size() != 1)
      {
      outerror("ERROR: Equation for mu is missing");
      return true;
      }

      distr_beta_sigma2s[distr_beta_sigma2s.size()-1].distrp.push_back
      (&distr_beta_mus[distr_beta_mus.size()-1]);

      distr_beta_mus[distr_beta_mus.size()-1].distrp.push_back
      (&distr_beta_sigma2s[distr_beta_sigma2s.size()-1]);


    }
//------------------------------- END: betainf1_tau -------------------------------


  //---------------------------------- dagum_p -----------------------------------
   else if (family.getvalue() == "dagum" && equationtype.getvalue()=="p")
     {

     computemodeforstartingvalues = true;

     distr_dagum_ps.push_back(DISTR_dagum_p(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_dagum_ps[distr_dagum_ps.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: dagum_p ---------------------------------

 //---------------------------------- dagum_b --------------------------------
   else if (family.getvalue() == "dagum" && equationtype.getvalue()=="b")
     {

     computemodeforstartingvalues = true;

     distr_dagum_bs.push_back(DISTR_dagum_b(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_dagum_bs[distr_dagum_bs.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------ END: dagum_b -------------------------------

 // ----------------------------------- dagum_a ----------------------------------
   else if (family.getvalue() == "dagum" && equationtype.getvalue()=="a")
     {

     if(generaloptions.copula)
      {
      mainequation=false;
      }
    else
      {
      mainequation=true;
      }

     computemodeforstartingvalues = true;

     distr_dagum_as.push_back(DISTR_dagum_a(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_dagum_as[distr_dagum_as.size()-1];
     equations[modnr].pathd = "";

     if (distr_dagum_ps.size() != distr_dagum_as.size())
       {
       outerror("ERROR: Equation for p is missing");
       return true;
       }

     if (distr_dagum_bs.size() != distr_dagum_as.size())
       {
       outerror("ERROR: Equation for b is missing");
       return true;
       }
     else
       {
       distr_dagum_ps[distr_dagum_ps.size()-1].distrp.push_back(
       &distr_dagum_bs[distr_dagum_bs.size()-1]);

	     distr_dagum_ps[distr_dagum_ps.size()-1].distrp.push_back(
       &distr_dagum_as[distr_dagum_as.size()-1]);

       distr_dagum_bs[distr_dagum_bs.size()-1].distrp.push_back(
       &distr_dagum_ps[distr_dagum_ps.size()-1]);

       distr_dagum_bs[distr_dagum_bs.size()-1].distrp.push_back(
       &distr_dagum_as[distr_dagum_as.size()-1]);

       distr_dagum_as[distr_dagum_as.size()-1].distrp.push_back(
       &distr_dagum_ps[distr_dagum_ps.size()-1]);

       distr_dagum_as[distr_dagum_as.size()-1].distrp.push_back(
       &distr_dagum_bs[distr_dagum_bs.size()-1]);

       countmarginal += 1;
       }
     if(generaloptions.copula){}
     else
       {
       predict_mult_distrs.push_back(&distr_dagum_ps[distr_dagum_ps.size()-1]);
       predict_mult_distrs.push_back(&distr_dagum_bs[distr_dagum_bs.size()-1]);
       predict_mult_distrs.push_back(&distr_dagum_as[distr_dagum_as.size()-1]);
       }

     if(generaloptions.copula)
      {
      if (countmarginal == 1)
        {
        distr_dagum_as[distr_dagum_as.size()-1].set_copulapos(0);
        distr_dagum_bs[distr_dagum_bs.size()-1].set_copulapos(0);
        distr_dagum_ps[distr_dagum_ps.size()-1].set_copulapos(0);
        }
      else if(countmarginal == 2)
        {
        distr_dagum_as[distr_dagum_as.size()-1].set_copulapos(1);
        distr_dagum_bs[distr_dagum_bs.size()-1].set_copulapos(1);
        distr_dagum_ps[distr_dagum_ps.size()-1].set_copulapos(1);
        }
      else
        {
        outerror("ERROR: currently only bivariate copula models implemented");
        return true;
        }
      }
     }
 //------------------------------ END: dagum_a ----------------------------------

  //---------------------------------- t_df -----------------------------------
   else if (family.getvalue() == "t" && equationtype.getvalue()=="df")
     {

     computemodeforstartingvalues = true;

     distr_t_dfs.push_back(DISTR_t_df(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_t_dfs[distr_t_dfs.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: t_df ---------------------------------

 //---------------------------------- t_sigma2 --------------------------------
   else if (family.getvalue() == "t" && equationtype.getvalue()=="sigma2")
     {

     computemodeforstartingvalues = true;

     distr_t_sigma2s.push_back(DISTR_t_sigma2(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_t_sigma2s[distr_t_sigma2s.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------ END: t_sigma2 -------------------------------


 // ----------------------------------- t_mu ----------------------------------
   else if (family.getvalue() == "t" && equationtype.getvalue()=="mu")
     {
     mainequation=true;

     computemodeforstartingvalues = true;

     distr_t_mus.push_back(DISTR_t_mu(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_t_mus[distr_t_mus.size()-1];
     equations[modnr].pathd = "";

     if (distr_t_dfs.size() != 1)
       {
       outerror("ERROR: Equation for degrees of freedom is missing");
       return true;
       }

     if (distr_t_sigma2s.size() != 1)
       {
       outerror("ERROR: Equation for sigma2 is missing");
       return true;
       }

     predict_mult_distrs.push_back(&distr_t_dfs[distr_t_dfs.size()-1]);
     predict_mult_distrs.push_back(&distr_t_sigma2s[distr_t_sigma2s.size()-1]);
     predict_mult_distrs.push_back(&distr_t_mus[distr_t_mus.size()-1]);

     distr_t_dfs[distr_t_dfs.size()-1].distrp.push_back(
     &distr_t_sigma2s[distr_t_sigma2s.size()-1]);

	 distr_t_dfs[distr_t_dfs.size()-1].distrp.push_back(
     &distr_t_mus[distr_t_mus.size()-1]);

     distr_t_sigma2s[distr_t_sigma2s.size()-1].distrp.push_back(
     &distr_t_dfs[distr_t_dfs.size()-1]);

     distr_t_sigma2s[distr_t_sigma2s.size()-1].distrp.push_back(
     &distr_t_mus[distr_t_mus.size()-1]);

     distr_t_mus[distr_t_mus.size()-1].distrp.push_back(
     &distr_t_dfs[distr_t_dfs.size()-1]);

     distr_t_mus[distr_t_mus.size()-1].distrp.push_back(
     &distr_t_sigma2s[distr_t_sigma2s.size()-1]);

     }
 //------------------------------ END: t_mu ----------------------------------

  //---------------------------------- BCCG_nu -----------------------------------
   else if (family.getvalue() == "BCCG" && equationtype.getvalue()=="nu")
     {

     computemodeforstartingvalues = true;

     distr_BCCG_nus.push_back(DISTR_BCCG_nu(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_BCCG_nus[distr_BCCG_nus.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: BCCG_nu ---------------------------------

 //---------------------------------- BCCG_sigma --------------------------------
   else if (family.getvalue() == "BCCG" && equationtype.getvalue()=="sigma")
     {

     computemodeforstartingvalues = true;

     distr_BCCG_sigmas.push_back(DISTR_BCCG_sigma(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------ END: BCCG_sigma -------------------------------


 // ----------------------------------- BCCG_mu ----------------------------------
   else if (family.getvalue() == "BCCG" && equationtype.getvalue()=="mu")
     {

     mainequation=true;

     computemodeforstartingvalues = true;

     distr_BCCG_mus.push_back(DISTR_BCCG_mu(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_BCCG_mus[distr_BCCG_mus.size()-1];
     equations[modnr].pathd = "";

     if (distr_BCCG_nus.size() != 1)
       {
       outerror("ERROR: Equation for degrees of freedom is missing");
       return true;
       }

     if (distr_BCCG_sigmas.size() != 1)
       {
       outerror("ERROR: Equation for sigma2 is missing");
       return true;
       }

     predict_mult_distrs.push_back(&distr_BCCG_nus[distr_BCCG_nus.size()-1]);
     predict_mult_distrs.push_back(&distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1]);
     predict_mult_distrs.push_back(&distr_BCCG_mus[distr_BCCG_mus.size()-1]);

     distr_BCCG_nus[distr_BCCG_nus.size()-1].distrp.push_back(
     &distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1]);

	 distr_BCCG_nus[distr_BCCG_nus.size()-1].distrp.push_back(
     &distr_BCCG_mus[distr_BCCG_mus.size()-1]);

     distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1].distrp.push_back(
     &distr_BCCG_nus[distr_BCCG_nus.size()-1]);

     distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1].distrp.push_back(
     &distr_BCCG_mus[distr_BCCG_mus.size()-1]);

     distr_BCCG_mus[distr_BCCG_mus.size()-1].distrp.push_back(
     &distr_BCCG_nus[distr_BCCG_nus.size()-1]);

     distr_BCCG_mus[distr_BCCG_mus.size()-1].distrp.push_back(
     &distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1]);

     }
 //------------------------------ END: BCCG_mu ----------------------------------


 //---------------------------------- bivt_rho -----------------------------------
   else if (family.getvalue() == "bivt" && equationtype.getvalue()=="rho")
     {

     //mainequation=true;

     computemodeforstartingvalues = true;

     distr_bivt_rhos.push_back(DISTR_bivt_rho(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_bivt_rhos[distr_bivt_rhos.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: bivt_rho ---------------------------------

 //---------------------------------- bivt_df -----------------------------------
   else if (family.getvalue() == "bivt" && equationtype.getvalue()=="df")
     {

     //mainequation=true;

     computemodeforstartingvalues = true;

     distr_bivt_dfs.push_back(DISTR_bivt_df(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_bivt_dfs[distr_bivt_dfs.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: bivt_df ---------------------------------

 //---------------------------------- bivt_sigma --------------------------------
   else if ((family.getvalue() == "bivt") && ((equationtype.getvalue()=="sigma") || (equationtype.getvalue()=="scale1") || (equationtype.getvalue()=="scale2")))
     {

     computemodeforstartingvalues = true;

      unsigned pos;
      if (distr_bivt_sigmas.size()==0)
        pos=0;
      else
        pos=1;

     distr_bivt_sigmas.push_back(DISTR_bivt_sigma(&generaloptions,D.getCol(0),pos,w));

     equations[modnr].distrp = &distr_bivt_sigmas[distr_bivt_sigmas.size()-1];
     equations[modnr].pathd = "";

     if (distr_bivt_sigmas.size() == 2)
       {
       distr_bivt_sigmas[distr_bivt_sigmas.size()-2].response2 = distr_bivt_sigmas[distr_bivt_sigmas.size()-1].response;
       distr_bivt_sigmas[distr_bivt_sigmas.size()-1].response2 = distr_bivt_sigmas[distr_bivt_sigmas.size()-2].response;
       }


     }
 //------------------------------ END: bivt_sigma -------------------------------


 // ----------------------------------- bivt_mu ----------------------------------
   else if ((family.getvalue() == "bivt") && ((equationtype.getvalue()=="mu"))
           )
     {

     //mainequation=true;

      computemodeforstartingvalues = true;

     unsigned pos;
      if (distr_bivt_mus.size()==0)
        pos=0;
      else
        {
        pos=1;
        mainequation=true;
        }


     distr_bivt_mus.push_back(DISTR_bivt_mu(&generaloptions,D.getCol(0),pos,w));

     equations[modnr].distrp = &distr_bivt_mus[distr_bivt_mus.size()-1];
     equations[modnr].pathd = "";

     if (distr_bivt_dfs.size() != 1)
       {
       outerror("ERROR: Equation for degrees of freedom is missing");
       return true;
       }

     if (distr_bivt_rhos.size() != 1)
       {
       outerror("ERROR: Equation for rho is missing");
       return true;
       }

     if (distr_bivt_sigmas.size() != 2)
       {
       outerror("ERROR: Equation for sigma is missing");
       return true;
       }


     if (pos==1)
       {

       predict_mult_distrs.push_back(&distr_bivt_dfs[distr_bivt_dfs.size()-1]);
       predict_mult_distrs.push_back(&distr_bivt_rhos[distr_bivt_rhos.size()-1]);
       predict_mult_distrs.push_back(&distr_bivt_sigmas[distr_bivt_sigmas.size()-2]);
       predict_mult_distrs.push_back(&distr_bivt_sigmas[distr_bivt_sigmas.size()-1]);
       predict_mult_distrs.push_back(&distr_bivt_mus[distr_bivt_mus.size()-2]);
       predict_mult_distrs.push_back(&distr_bivt_mus[distr_bivt_mus.size()-1]);

       distr_bivt_mus[distr_bivt_mus.size()-2].response2 = distr_bivt_mus[distr_bivt_mus.size()-1].response;
       distr_bivt_mus[distr_bivt_mus.size()-1].response2 = distr_bivt_mus[distr_bivt_mus.size()-2].response;
       distr_bivt_rhos[distr_bivt_rhos.size()-1].response2 = distr_bivt_mus[distr_bivt_mus.size()-2].response;
       distr_bivt_dfs[distr_bivt_dfs.size()-1].response2 = distr_bivt_mus[distr_bivt_mus.size()-2].response;
    //   distr_bivt_dfs[distr_bivt_dfs.size()-1].response = distr_bivt_mus[distr_bivt_mus.size()-1].response;

       distr_bivt_dfs[distr_bivt_dfs.size()-1].distrp.push_back(
       &distr_bivt_rhos[distr_bivt_rhos.size()-1]);

       distr_bivt_dfs[distr_bivt_dfs.size()-1].distrp.push_back(
       &distr_bivt_sigmas[distr_bivt_sigmas.size()-2]);

       distr_bivt_dfs[distr_bivt_dfs.size()-1].distrp.push_back(
       &distr_bivt_sigmas[distr_bivt_sigmas.size()-1]);

	   distr_bivt_dfs[distr_bivt_dfs.size()-1].distrp.push_back(
       &distr_bivt_mus[distr_bivt_mus.size()-2]);

       distr_bivt_dfs[distr_bivt_dfs.size()-1].distrp.push_back(
       &distr_bivt_mus[distr_bivt_mus.size()-1]);

       distr_bivt_rhos[distr_bivt_rhos.size()-1].distrp.push_back(
       &distr_bivt_dfs[distr_bivt_dfs.size()-1]);

       distr_bivt_rhos[distr_bivt_rhos.size()-1].distrp.push_back(
       &distr_bivt_sigmas[distr_bivt_sigmas.size()-2]);

       distr_bivt_rhos[distr_bivt_rhos.size()-1].distrp.push_back(
       &distr_bivt_sigmas[distr_bivt_sigmas.size()-1]);

	   distr_bivt_rhos[distr_bivt_rhos.size()-1].distrp.push_back(
       &distr_bivt_mus[distr_bivt_mus.size()-2]);

       distr_bivt_rhos[distr_bivt_rhos.size()-1].distrp.push_back(
       &distr_bivt_mus[distr_bivt_mus.size()-1]);

       distr_bivt_sigmas[distr_bivt_sigmas.size()-2].distrp.push_back(
       &distr_bivt_dfs[distr_bivt_dfs.size()-1]);

       distr_bivt_sigmas[distr_bivt_sigmas.size()-2].distrp.push_back(
       &distr_bivt_rhos[distr_bivt_rhos.size()-1]);

        distr_bivt_sigmas[distr_bivt_sigmas.size()-2].distrp.push_back(
       &distr_bivt_sigmas[distr_bivt_sigmas.size()-1]);

       distr_bivt_sigmas[distr_bivt_sigmas.size()-2].distrp.push_back(
       &distr_bivt_mus[distr_bivt_mus.size()-2]);

        distr_bivt_sigmas[distr_bivt_sigmas.size()-2].distrp.push_back(
       &distr_bivt_mus[distr_bivt_mus.size()-1]);

       distr_bivt_sigmas[distr_bivt_sigmas.size()-1].distrp.push_back(
       &distr_bivt_dfs[distr_bivt_dfs.size()-1]);

       distr_bivt_sigmas[distr_bivt_sigmas.size()-1].distrp.push_back(
       &distr_bivt_rhos[distr_bivt_rhos.size()-1]);

        distr_bivt_sigmas[distr_bivt_sigmas.size()-1].distrp.push_back(
       &distr_bivt_sigmas[distr_bivt_sigmas.size()-2]);

       distr_bivt_sigmas[distr_bivt_sigmas.size()-1].distrp.push_back(
       &distr_bivt_mus[distr_bivt_mus.size()-1]);

        distr_bivt_sigmas[distr_bivt_sigmas.size()-1].distrp.push_back(
       &distr_bivt_mus[distr_bivt_mus.size()-2]);

        distr_bivt_mus[distr_bivt_mus.size()-2].distrp.push_back(
       &distr_bivt_dfs[distr_bivt_dfs.size()-1]);

       distr_bivt_mus[distr_bivt_mus.size()-2].distrp.push_back(
       &distr_bivt_rhos[distr_bivt_rhos.size()-1]);

       distr_bivt_mus[distr_bivt_mus.size()-2].distrp.push_back(
       &distr_bivt_mus[distr_bivt_mus.size()-1]);

       distr_bivt_mus[distr_bivt_mus.size()-2].distrp.push_back(
       &distr_bivt_sigmas[distr_bivt_sigmas.size()-2]);

       distr_bivt_mus[distr_bivt_mus.size()-2].distrp.push_back(
       &distr_bivt_sigmas[distr_bivt_sigmas.size()-1]);

        distr_bivt_mus[distr_bivt_mus.size()-1].distrp.push_back(
       &distr_bivt_dfs[distr_bivt_dfs.size()-1]);

       distr_bivt_mus[distr_bivt_mus.size()-1].distrp.push_back(
       &distr_bivt_rhos[distr_bivt_rhos.size()-1]);

       distr_bivt_mus[distr_bivt_mus.size()-1].distrp.push_back(
       &distr_bivt_mus[distr_bivt_mus.size()-2]);

       distr_bivt_mus[distr_bivt_mus.size()-1].distrp.push_back(
       &distr_bivt_sigmas[distr_bivt_sigmas.size()-1]);

       distr_bivt_mus[distr_bivt_mus.size()-1].distrp.push_back(
       &distr_bivt_sigmas[distr_bivt_sigmas.size()-2]);
       }

     }
 //------------------------------ END: bivt_mu ----------------------------------

//---------------------------------- bivnormal_rhofz -----------------------------------
   else if (family.getvalue() == "bivnormal_fz" && equationtype.getvalue()=="rho")
     {

     computemodeforstartingvalues = true;

     distr_bivnormal_rhofzs.push_back(DISTR_bivnormal_rhofz(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: bivnormal_rhofz ---------------------------------

//---------------------------------- bivnormal_rho -----------------------------------
   else if (family.getvalue() == "bivnormal" && equationtype.getvalue()=="rho")
     {

     computemodeforstartingvalues = true;

     distr_bivnormal_rhos.push_back(DISTR_bivnormal_rho(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: bivnormal_rho ---------------------------------

 //---------------------------------- bivnormal_sigma --------------------------------
   else if ((family.getvalue() == "bivnormal") && ((equationtype.getvalue()=="sigma") || (equationtype.getvalue()=="scale1") || (equationtype.getvalue()=="scale2")))
     {

     computemodeforstartingvalues = true;

      unsigned pos;
      if (distr_bivnormal_sigmas.size()==0)
        pos=0;
      else
        pos=1;

     distr_bivnormal_sigmas.push_back(DISTR_bivnormal_sigma(&generaloptions,D.getCol(0),pos,w));

     equations[modnr].distrp = &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1];
     equations[modnr].pathd = "";

     if (distr_bivnormal_sigmas.size() == 2)
       {
       distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2].response2 = distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1].response;
       distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1].response2 = distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2].response;
       }


     }
 //------------------------------ END: bivnormal_sigma -------------------------------


 // ----------------------------------- bivnormal_mu ----------------------------------
   else if ((family.getvalue() == "bivnormal") && ((equationtype.getvalue()=="mu"))
           )
     {

    // computemodeforstartingvalues = true;

     unsigned pos;
     if (distr_bivnormal_mus.size()==0)
       {
       pos = 0;
       }
     else
       {
       pos = 1;
       mainequation=true;
       }

     distr_bivnormal_mus.push_back(DISTR_bivnormal_mu(&generaloptions,D.getCol(0),pos,w));

     equations[modnr].distrp = &distr_bivnormal_mus[distr_bivnormal_mus.size()-1];
     equations[modnr].pathd = "";

     if (distr_bivnormal_rhos.size() != 1)
       {
       outerror("ERROR: Equation for rho is missing");
       return true;
       }

     if (distr_bivnormal_sigmas.size() != 2)
       {
       outerror("ERROR: Equation for sigma is missing");
       return true;
       }


     if (pos==1)
       {

       predict_mult_distrs.push_back(&distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1]);
       predict_mult_distrs.push_back(&distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2]);
       predict_mult_distrs.push_back(&distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1]);
       predict_mult_distrs.push_back(&distr_bivnormal_mus[distr_bivnormal_mus.size()-2]);
       predict_mult_distrs.push_back(&distr_bivnormal_mus[distr_bivnormal_mus.size()-1]);

       distr_bivnormal_mus[distr_bivnormal_mus.size()-2].response2 = distr_bivnormal_mus[distr_bivnormal_mus.size()-1].response;
       distr_bivnormal_mus[distr_bivnormal_mus.size()-1].response2 = distr_bivnormal_mus[distr_bivnormal_mus.size()-2].response;
       distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1].response2 = distr_bivnormal_mus[distr_bivnormal_mus.size()-2].response;


       distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2]);

       distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1]);

	   distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1].distrp.push_back(
       &distr_bivnormal_mus[distr_bivnormal_mus.size()-2]);

       distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1].distrp.push_back(
       &distr_bivnormal_mus[distr_bivnormal_mus.size()-1]);

       distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2].distrp.push_back(
       &distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1]);

        distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1]);

       distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2].distrp.push_back(
       &distr_bivnormal_mus[distr_bivnormal_mus.size()-2]);

        distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2].distrp.push_back(
       &distr_bivnormal_mus[distr_bivnormal_mus.size()-1]);

       distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1].distrp.push_back(
       &distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1]);

        distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2]);

       distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1].distrp.push_back(
       &distr_bivnormal_mus[distr_bivnormal_mus.size()-1]);

        distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1].distrp.push_back(
       &distr_bivnormal_mus[distr_bivnormal_mus.size()-2]);

       distr_bivnormal_mus[distr_bivnormal_mus.size()-2].distrp.push_back(
       &distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1]);

       distr_bivnormal_mus[distr_bivnormal_mus.size()-2].distrp.push_back(
       &distr_bivnormal_mus[distr_bivnormal_mus.size()-1]);

       distr_bivnormal_mus[distr_bivnormal_mus.size()-2].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2]);

       distr_bivnormal_mus[distr_bivnormal_mus.size()-2].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1]);

       distr_bivnormal_mus[distr_bivnormal_mus.size()-1].distrp.push_back(
       &distr_bivnormal_rhos[distr_bivnormal_rhos.size()-1]);

       distr_bivnormal_mus[distr_bivnormal_mus.size()-1].distrp.push_back(
       &distr_bivnormal_mus[distr_bivnormal_mus.size()-2]);

       distr_bivnormal_mus[distr_bivnormal_mus.size()-1].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1]);

       distr_bivnormal_mus[distr_bivnormal_mus.size()-1].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2]);
       }

     }
 //------------------------------ END: bivnormal_mu ----------------------------------



 // ----------------------------------- bivnormal_mufz ----------------------------------
   else if ((family.getvalue() == "bivnormal_fz") && (equationtype.getvalue()=="mu")
           )
     {

    // computemodeforstartingvalues = true;

     unsigned pos;
     if (distr_bivnormal_mufzs.size()==0)
       pos = 0;
     else
       {
       pos = 1;
       mainequation=true;
       }

     distr_bivnormal_mufzs.push_back(DISTR_bivnormal_mufz(&generaloptions,D.getCol(0),pos,w));

     equations[modnr].distrp = &distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1];
     equations[modnr].pathd = "";

     if (distr_bivnormal_rhofzs.size() != 1)
       {
       outerror("ERROR: Equation for rho is missing");
       return true;
       }

     if (distr_bivnormal_sigmas.size() != 2)
       {
       outerror("ERROR: Equation for sigma is missing");
       return true;
       }


     if (pos==1)
       {


       predict_mult_distrs.push_back(&distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1]);
       predict_mult_distrs.push_back(&distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2]);
       predict_mult_distrs.push_back(&distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1]);
       predict_mult_distrs.push_back(&distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2]);
       predict_mult_distrs.push_back(&distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1]);

       distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2].response2 = distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1].response;
       distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1].response2 = distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2].response;
       distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1].response2 = distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2].response;


       distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2]);

       distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1]);

	   distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1].distrp.push_back(
       &distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2]);

       distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1].distrp.push_back(
       &distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1]);

       distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2].distrp.push_back(
       &distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1]);

        distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1]);

       distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2].distrp.push_back(
       &distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2]);

        distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2].distrp.push_back(
       &distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1]);

       distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1].distrp.push_back(
       &distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1]);

        distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2]);

       distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1].distrp.push_back(
       &distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1]);

        distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1].distrp.push_back(
       &distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2]);

       distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2].distrp.push_back(
       &distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1]);

       distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2].distrp.push_back(
       &distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1]);

       distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2]);

       distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1]);

       distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1].distrp.push_back(
       &distr_bivnormal_rhofzs[distr_bivnormal_rhofzs.size()-1]);

       distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1].distrp.push_back(
       &distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-2]);

       distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-1]);

       distr_bivnormal_mufzs[distr_bivnormal_mufzs.size()-1].distrp.push_back(
       &distr_bivnormal_sigmas[distr_bivnormal_sigmas.size()-2]);
       }

     }
 //------------------------------ END: bivnormal_mufz ----------------------------------


   //---------------------------------- bivprobit_rho -----------------------------------
   else if (family.getvalue() == "bivprobit" && equationtype.getvalue()=="rho")
     {

     computemodeforstartingvalues = true;

     distr_bivprobit_rhos.push_back(DISTR_bivprobit_rho(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_bivprobit_rhos[distr_bivprobit_rhos.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: bivprobit_rho ---------------------------------

  // ----------------------------------- bivprobit_mu ----------------------------------
   else if ((family.getvalue() == "bivprobit") && (equationtype.getvalue()=="mu"))
     {

     computemodeforstartingvalues = true;

     unsigned pos;
     if (distr_bivprobit_mus.size()==0)
       pos = 0;
     else
       {
       pos = 1;
       mainequation=true;
       }

     distr_bivprobit_mus.push_back(DISTR_bivprobit_mu(&generaloptions,D.getCol(0),pos,w));

     equations[modnr].distrp = &distr_bivprobit_mus[distr_bivprobit_mus.size()-1];
     equations[modnr].pathd = "";

     if (distr_bivprobit_rhos.size() != 1)
       {
       outerror("ERROR: Equation for rho is missing");
       return true;
       }


     if (pos==1)
       {


       predict_mult_distrs.push_back(&distr_bivprobit_rhos[distr_bivprobit_rhos.size()-1]);
       predict_mult_distrs.push_back(&distr_bivprobit_mus[distr_bivprobit_mus.size()-2]);
       predict_mult_distrs.push_back(&distr_bivprobit_mus[distr_bivprobit_mus.size()-1]);

       distr_bivprobit_mus[distr_bivprobit_mus.size()-2].workingresponse2p = &distr_bivprobit_mus[distr_bivprobit_mus.size()-1].response;
       distr_bivprobit_mus[distr_bivprobit_mus.size()-1].workingresponse2p = &distr_bivprobit_mus[distr_bivprobit_mus.size()-2].response;
       distr_bivprobit_rhos[distr_bivprobit_rhos.size()-1].workingresponse1p = &distr_bivprobit_mus[distr_bivprobit_mus.size()-1].response;
       distr_bivprobit_rhos[distr_bivprobit_rhos.size()-1].workingresponse2p = &distr_bivprobit_mus[distr_bivprobit_mus.size()-2].response;


	   distr_bivprobit_rhos[distr_bivprobit_rhos.size()-1].distrp.push_back(
       &distr_bivprobit_mus[distr_bivprobit_mus.size()-2]);

       distr_bivprobit_rhos[distr_bivprobit_rhos.size()-1].distrp.push_back(
       &distr_bivprobit_mus[distr_bivprobit_mus.size()-1]);

       distr_bivprobit_mus[distr_bivprobit_mus.size()-2].distrp.push_back(
       &distr_bivprobit_rhos[distr_bivprobit_rhos.size()-1]);

       distr_bivprobit_mus[distr_bivprobit_mus.size()-2].distrp.push_back(
       &distr_bivprobit_mus[distr_bivprobit_mus.size()-1]);

       distr_bivprobit_mus[distr_bivprobit_mus.size()-1].distrp.push_back(
       &distr_bivprobit_rhos[distr_bivprobit_rhos.size()-1]);

       distr_bivprobit_mus[distr_bivprobit_mus.size()-1].distrp.push_back(
       &distr_bivprobit_mus[distr_bivprobit_mus.size()-2]);

        }

     }
 //------------------------------ END: bivprobit_mu ----------------------------------


//---------------------------------- bivprobit2_rho -----------------------------------
   else if (family.getvalue() == "bivprobit2" && equationtype.getvalue()=="rho")
     {

     computemodeforstartingvalues = true;

     distr_bivprobit2_rhos.push_back(DISTR_bivprobit2_rho(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_bivprobit2_rhos[distr_bivprobit2_rhos.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: bivprobit2_rho ---------------------------------

  // ----------------------------------- bivprobit2_mu ----------------------------------
   else if ((family.getvalue() == "bivprobit2") && (equationtype.getvalue()=="mu"))
     {

     computemodeforstartingvalues = true;

     unsigned pos;
     if (distr_bivprobit2_mus.size()==0)
       pos = 0;
     else
       {
       pos = 1;
       mainequation=true;
       }

     distr_bivprobit2_mus.push_back(DISTR_bivprobit2_mu(&generaloptions,D.getCol(0),pos,w));

     equations[modnr].distrp = &distr_bivprobit2_mus[distr_bivprobit2_mus.size()-1];
     equations[modnr].pathd = "";

     if (distr_bivprobit2_rhos.size() != 1)
       {
       outerror("ERROR: Equation for rho is missing");
       return true;
       }


     if (pos==1)
       {


       predict_mult_distrs.push_back(&distr_bivprobit2_rhos[distr_bivprobit2_rhos.size()-1]);
       predict_mult_distrs.push_back(&distr_bivprobit2_mus[distr_bivprobit2_mus.size()-2]);
       predict_mult_distrs.push_back(&distr_bivprobit2_mus[distr_bivprobit2_mus.size()-1]);

       distr_bivprobit2_mus[distr_bivprobit2_mus.size()-2].workingresponse2p = &distr_bivprobit2_mus[distr_bivprobit2_mus.size()-1].response;
       distr_bivprobit2_mus[distr_bivprobit2_mus.size()-1].workingresponse2p = &distr_bivprobit2_mus[distr_bivprobit2_mus.size()-2].response;
       distr_bivprobit2_rhos[distr_bivprobit2_rhos.size()-1].workingresponse1p = &distr_bivprobit2_mus[distr_bivprobit2_mus.size()-1].response;
       distr_bivprobit2_rhos[distr_bivprobit2_rhos.size()-1].workingresponse2p = &distr_bivprobit2_mus[distr_bivprobit2_mus.size()-2].response;


	   distr_bivprobit2_rhos[distr_bivprobit2_rhos.size()-1].distrp.push_back(
       &distr_bivprobit2_mus[distr_bivprobit2_mus.size()-2]);

       distr_bivprobit2_rhos[distr_bivprobit2_rhos.size()-1].distrp.push_back(
       &distr_bivprobit2_mus[distr_bivprobit2_mus.size()-1]);

       distr_bivprobit2_mus[distr_bivprobit2_mus.size()-2].distrp.push_back(
       &distr_bivprobit2_rhos[distr_bivprobit2_rhos.size()-1]);

       distr_bivprobit2_mus[distr_bivprobit2_mus.size()-2].distrp.push_back(
       &distr_bivprobit2_mus[distr_bivprobit2_mus.size()-1]);

       distr_bivprobit2_mus[distr_bivprobit2_mus.size()-1].distrp.push_back(
       &distr_bivprobit2_rhos[distr_bivprobit2_rhos.size()-1]);

       distr_bivprobit2_mus[distr_bivprobit2_mus.size()-1].distrp.push_back(
       &distr_bivprobit2_mus[distr_bivprobit2_mus.size()-2]);

       }

     }
 //------------------------------ END: bivprobit2_mu ----------------------------------


    //---------------------------------- bivlogit_or -----------------------------------
   else if (family.getvalue() == "bivlogit" && equationtype.getvalue()=="oddsratio")
     {

     computemodeforstartingvalues = true;

     distr_bivlogit_ors.push_back(DISTR_bivlogit_or(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_bivlogit_ors[distr_bivlogit_ors.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: bivlogit_or ---------------------------------

  // ----------------------------------- bivlogit_mu ----------------------------------
   else if ((family.getvalue() == "bivlogit") && (equationtype.getvalue()=="mu"))

     {

     computemodeforstartingvalues = true;

     unsigned pos;
     if (distr_bivlogit_mus.size()==0)
       pos = 0;
     else
       {
       pos = 1;
       mainequation=true;
       }

     distr_bivlogit_mus.push_back(DISTR_bivlogit_mu(&generaloptions,D.getCol(0),pos,w));

     equations[modnr].distrp = &distr_bivlogit_mus[distr_bivlogit_mus.size()-1];
     equations[modnr].pathd = "";

     if (distr_bivlogit_ors.size() != 1)
       {
       outerror("ERROR: Equation for oddsratio is missing");
       return true;
       }


     if (pos==1)
       {


       predict_mult_distrs.push_back(&distr_bivlogit_ors[distr_bivlogit_ors.size()-1]);
       predict_mult_distrs.push_back(&distr_bivlogit_mus[distr_bivlogit_mus.size()-2]);
       predict_mult_distrs.push_back(&distr_bivlogit_mus[distr_bivlogit_mus.size()-1]);

       distr_bivlogit_mus[distr_bivlogit_mus.size()-2].response2 = distr_bivlogit_mus[distr_bivlogit_mus.size()-1].response;
       distr_bivlogit_mus[distr_bivlogit_mus.size()-1].response2 = distr_bivlogit_mus[distr_bivlogit_mus.size()-2].response;
       distr_bivlogit_ors[distr_bivlogit_ors.size()-1].response2 = distr_bivlogit_mus[distr_bivlogit_mus.size()-2].response;


	   distr_bivlogit_ors[distr_bivlogit_ors.size()-1].distrp.push_back(
       &distr_bivlogit_mus[distr_bivlogit_mus.size()-2]);

       distr_bivlogit_ors[distr_bivlogit_ors.size()-1].distrp.push_back(
       &distr_bivlogit_mus[distr_bivlogit_mus.size()-1]);

       distr_bivlogit_mus[distr_bivlogit_mus.size()-2].distrp.push_back(
       &distr_bivlogit_ors[distr_bivlogit_ors.size()-1]);

       distr_bivlogit_mus[distr_bivlogit_mus.size()-2].distrp.push_back(
       &distr_bivlogit_mus[distr_bivlogit_mus.size()-1]);

       distr_bivlogit_mus[distr_bivlogit_mus.size()-1].distrp.push_back(
       &distr_bivlogit_ors[distr_bivlogit_ors.size()-1]);

       distr_bivlogit_mus[distr_bivlogit_mus.size()-1].distrp.push_back(
       &distr_bivlogit_mus[distr_bivlogit_mus.size()-2]);

        }

     }
 //------------------------------ END: bivlogit_mu ----------------------------------

//----------------------------- dirichlet -------------------------------
  else if ((family.getvalue() == "dirichlet") && ((equationtype.getvalue()=="alpha")
                                || (equationtype.getvalue()=="alpha2")|| (equationtype.getvalue()=="alpha3")|| (equationtype.getvalue()=="alpha4")
                                || (equationtype.getvalue()=="alpha5")|| (equationtype.getvalue()=="alpha6")||
(equationtype.getvalue()=="alpha7")))
    {
    computemodeforstartingvalues = true;

    int nrc = nrcat.getvalue();

    unsigned pos;
    unsigned countpos;
    for (countpos=0;countpos<(nrc);countpos++)
    {
       if (distr_dirichlets.size()==countpos)
        pos=countpos;
    }



    distr_dirichlets.push_back(DISTR_dirichlet(&generaloptions,D.getCol(0),nrc,pos,w));

    equations[modnr].distrp = &distr_dirichlets[distr_dirichlets.size()-1];
    equations[modnr].pathd = "";


    if (distr_dirichlets.size()>10)
       {
       outerror("ERROR: Number of equations has to be between 2 and 10");
       return true;
       }


   if (distr_dirichlets.size() == nrc)
       {

       mainequation=true;

           unsigned i;
           for(i=(nrc);i>0;i--)
                predict_mult_distrs.push_back(&distr_dirichlets[distr_dirichlets.size()-i]);

           unsigned k;
           for(k=(nrc);k>0;k--) {
              unsigned j;
             for(j=(nrc);j>0;j--) {

                if(k==j) {

                }
                else {
                    distr_dirichlets[distr_dirichlets.size()-k].distrp.push_back(
                    &distr_dirichlets[distr_dirichlets.size()-j]);
                }
             }
           }

       }

    }
//-------------------------- END: dirichlet -----------------------------

 //---------------------------------- gengamma tau -----------------------------------
   else if (family.getvalue() == "gengamma" && equationtype.getvalue()=="tau")
     {

     computemodeforstartingvalues = true;

     distr_gengamma_taus.push_back(DISTR_gengamma_tau(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_gengamma_taus[distr_gengamma_taus.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: gengamma tau ---------------------------------

 //---------------------------------- gengamma sigma --------------------------------
   else if (family.getvalue() == "gengamma" && equationtype.getvalue()=="sigma")
     {

     computemodeforstartingvalues = true;

     distr_gengamma_sigmas.push_back(DISTR_gengamma_sigma(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_gengamma_sigmas[distr_gengamma_sigmas.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------ END: gengamma sigma -------------------------------


 // ----------------------------------- gengamma mu ----------------------------------
   else if (family.getvalue() == "gengamma" && equationtype.getvalue()=="mu")
     {

     mainequation=true;
     computemodeforstartingvalues = true;

     distr_gengamma_mus.push_back(DISTR_gengamma_mu(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_gengamma_mus[distr_gengamma_mus.size()-1];
     equations[modnr].pathd = "";

     if (distr_gengamma_taus.size() != 1)
       {
       outerror("ERROR: Equation for tau is missing");
       return true;
       }

     if (distr_gengamma_sigmas.size() != 1)
       {
       outerror("ERROR: Equation for sigma is missing");
       return true;
       }

     predict_mult_distrs.push_back(&distr_gengamma_taus[distr_gengamma_taus.size()-1]);
     predict_mult_distrs.push_back(&distr_gengamma_sigmas[distr_gengamma_sigmas.size()-1]);
     predict_mult_distrs.push_back(&distr_gengamma_mus[distr_gengamma_mus.size()-1]);

     distr_gengamma_taus[distr_gengamma_taus.size()-1].distrp.push_back(
     &distr_gengamma_sigmas[distr_gengamma_sigmas.size()-1]);

	 distr_gengamma_taus[distr_gengamma_taus.size()-1].distrp.push_back(
     &distr_gengamma_mus[distr_gengamma_mus.size()-1]);

     distr_gengamma_sigmas[distr_gengamma_sigmas.size()-1].distrp.push_back(
     &distr_gengamma_taus[distr_gengamma_taus.size()-1]);

     distr_gengamma_sigmas[distr_gengamma_sigmas.size()-1].distrp.push_back(
     &distr_gengamma_mus[distr_gengamma_mus.size()-1]);

     distr_gengamma_mus[distr_gengamma_mus.size()-1].distrp.push_back(
     &distr_gengamma_taus[distr_gengamma_taus.size()-1]);

     distr_gengamma_mus[distr_gengamma_mus.size()-1].distrp.push_back(
     &distr_gengamma_sigmas[distr_gengamma_sigmas.size()-1]);

     }
 //------------------------------ END: gengamma mu ----------------------------------


  //-------------------------------- pareto_p ---------------------------------
  else if (family.getvalue() == "pareto" && equationtype.getvalue()=="p")
    {

    computemodeforstartingvalues = true;

    distr_pareto_ps.push_back(DISTR_pareto_p(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_pareto_ps[distr_pareto_ps.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_pareto_ps[distr_pareto_ps.size()-1]);

    }
//---------------------------- END: pareto_p -------------------------------


//------------------------------- pareto_b ------------------------------------
  else if (family.getvalue() == "pareto" && equationtype.getvalue()=="b")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_pareto_bs.push_back(DISTR_pareto_b(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_pareto_bs[distr_pareto_bs.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_pareto_bs[distr_pareto_bs.size()-1]);

    if (distr_pareto_ps.size() != 1)
      {
      outerror("ERROR: Equation for p is missing");
      return true;
      }
    else
      {
      distr_pareto_ps[distr_pareto_ps.size()-1].distrp.push_back
      (&distr_pareto_bs[distr_pareto_bs.size()-1]);

      distr_pareto_bs[distr_pareto_bs.size()-1].distrp.push_back
      (&distr_pareto_ps[distr_pareto_ps.size()-1]);
      }

    }
//------------------------------- END: pareto_b -------------------------------


 //-------------------------------- weibull_alpha ---------------------------------
  else if (family.getvalue() == "weibull" && equationtype.getvalue()=="alpha")
    {

    computemodeforstartingvalues = true;

    distr_weibull_alphas.push_back(DISTR_weibull_alpha(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_weibull_alphas[distr_weibull_alphas.size()-1];
    equations[modnr].pathd = "";

    }
//---------------------------- END: weibull_alpha -------------------------------


//------------------------------- weibull_lambda ------------------------------------
  else if (family.getvalue() == "weibull" && equationtype.getvalue()=="lambda")
    {
    if(generaloptions.copula)
      {
      mainequation=false;
      }
    else
      {
      mainequation=true;
      }


    computemodeforstartingvalues = true;

    distr_weibull_lambdas.push_back(DISTR_weibull_lambda(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_weibull_lambdas[distr_weibull_lambdas.size()-1];
    equations[modnr].pathd = "";

    if (distr_weibull_alphas.size() != distr_weibull_lambdas.size())
      {
      outerror("ERROR: Equation for alpha is missing");
      return true;
      }
    else
      {
      distr_weibull_alphas[distr_weibull_alphas.size()-1].distrp.push_back
      (&distr_weibull_lambdas[distr_weibull_lambdas.size()-1]);

      distr_weibull_lambdas[distr_weibull_lambdas.size()-1].distrp.push_back
      (&distr_weibull_alphas[distr_weibull_alphas.size()-1]);

      countmarginal += 1;
      }
    if(generaloptions.copula){}
    else
      {
      predict_mult_distrs.push_back(&distr_weibull_alphas[distr_weibull_alphas.size()-1]);
      predict_mult_distrs.push_back(&distr_weibull_lambdas[distr_weibull_lambdas.size()-1]);
      }

    if(generaloptions.copula)
      {
      if (countmarginal == 1)
        {
        distr_weibull_lambdas[distr_weibull_lambdas.size()-1].set_copulapos(0);
        distr_weibull_alphas[distr_weibull_alphas.size()-1].set_copulapos(0);
        }
      else if(countmarginal == 2)
        {
        distr_weibull_lambdas[distr_weibull_lambdas.size()-1].set_copulapos(1);
        distr_weibull_alphas[distr_weibull_alphas.size()-1].set_copulapos(1);
        }
      else
        {
        outerror("ERROR: currently only bivariate copula models implemented");
        return true;
        }
      }

    }
//------------------------------- END: weibull_lambda -------------------------------


//-------------------------------- weibull2_alpha ---------------------------------
  else if (family.getvalue() == "weibull2" && equationtype.getvalue()=="alpha")
    {

    computemodeforstartingvalues = true;

    distr_weibull2_alphas.push_back(DISTR_weibull2_alpha(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_weibull2_alphas[distr_weibull2_alphas.size()-1];
    equations[modnr].pathd = "";

    }
//---------------------------- END: weibull2_alpha -------------------------------


//------------------------------- weibull2_lambda ------------------------------------
  else if (family.getvalue() == "weibull2" && equationtype.getvalue()=="lambda")
    {
    if(generaloptions.copula)
      {
      mainequation=false;
      }
    else
      {
      mainequation=true;
      }


    computemodeforstartingvalues = true;

    distr_weibull2_lambdas.push_back(DISTR_weibull2_lambda(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_weibull2_lambdas[distr_weibull2_lambdas.size()-1];
    equations[modnr].pathd = "";

    if (distr_weibull2_alphas.size() != distr_weibull2_lambdas.size())
      {
      outerror("ERROR: Equation for alpha is missing");
      return true;
      }
    else
      {
      distr_weibull2_alphas[distr_weibull2_alphas.size()-1].distrp.push_back
      (&distr_weibull2_lambdas[distr_weibull2_lambdas.size()-1]);

      distr_weibull2_lambdas[distr_weibull2_lambdas.size()-1].distrp.push_back
      (&distr_weibull2_alphas[distr_weibull2_alphas.size()-1]);

      countmarginal += 1;
      }
    if(generaloptions.copula){}
    else
      {
      predict_mult_distrs.push_back(&distr_weibull2_alphas[distr_weibull2_alphas.size()-1]);
      predict_mult_distrs.push_back(&distr_weibull2_lambdas[distr_weibull2_lambdas.size()-1]);
      }

    if(generaloptions.copula)
      {
      if (countmarginal == 1)
        {
        distr_weibull2_lambdas[distr_weibull2_lambdas.size()-1].set_copulapos(0);
        distr_weibull2_alphas[distr_weibull2_alphas.size()-1].set_copulapos(0);
        }
      else if(countmarginal == 2)
        {
        distr_weibull2_lambdas[distr_weibull2_lambdas.size()-1].set_copulapos(1);
        distr_weibull2_alphas[distr_weibull2_alphas.size()-1].set_copulapos(1);
        }
      else
        {
        outerror("ERROR: currently only bivariate copula models implemented");
        return true;
        }
      }

    }
//------------------------------- END: weibull2_lambda -------------------------------


 //-------------------------------- gumbel_sigma ---------------------------------
  else if (family.getvalue() == "gumbel" && equationtype.getvalue()=="sigma")
    {

    computemodeforstartingvalues = true;

    distr_gumbel_sigmas.push_back(DISTR_gumbel_sigma(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gumbel_sigmas[distr_gumbel_sigmas.size()-1];
    equations[modnr].pathd = "";

    }
//---------------------------- END: gumbel_sigma -------------------------------


//------------------------------- gumbel_mu ------------------------------------
  else if (family.getvalue() == "gumbel" && equationtype.getvalue()=="mu")
    {
    if(generaloptions.copula)
      {
      mainequation=false;
      }
    else
      {
      mainequation=true;
      }


    computemodeforstartingvalues = true;

    distr_gumbel_mus.push_back(DISTR_gumbel_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gumbel_mus[distr_gumbel_mus.size()-1];
    equations[modnr].pathd = "";

    if (distr_gumbel_sigmas.size() != distr_gumbel_mus.size())
      {
      outerror("ERROR: Equation for sigma is missing");
      return true;
      }
    else
      {
      distr_gumbel_sigmas[distr_gumbel_sigmas.size()-1].distrp.push_back
      (&distr_gumbel_mus[distr_gumbel_mus.size()-1]);

      distr_gumbel_mus[distr_gumbel_mus.size()-1].distrp.push_back
      (&distr_gumbel_sigmas[distr_gumbel_sigmas.size()-1]);

      countmarginal += 1;
      }
    if(generaloptions.copula){}
    else
      {
      predict_mult_distrs.push_back(&distr_gumbel_sigmas[distr_gumbel_sigmas.size()-1]);
      predict_mult_distrs.push_back(&distr_gumbel_mus[distr_gumbel_mus.size()-1]);
      }

    if(generaloptions.copula)
      {
      if (countmarginal == 1)
        {
        distr_gumbel_mus[distr_gumbel_mus.size()-1].set_copulapos(0);
        distr_gumbel_sigmas[distr_gumbel_sigmas.size()-1].set_copulapos(0);
        }
      else if(countmarginal == 2)
        {
        distr_gumbel_mus[distr_gumbel_mus.size()-1].set_copulapos(1);
        distr_gumbel_sigmas[distr_gumbel_sigmas.size()-1].set_copulapos(1);
        }
      else
        {
        outerror("ERROR: currently only bivariate copula models implemented");
        return true;
        }
      }

    }
//------------------------------- END: gumbel_mu -------------------------------


 //-------------------------------- gumbel2_sigma2 ---------------------------------
  else if (family.getvalue() == "gumbel2" && equationtype.getvalue()=="sigma2")
    {

    computemodeforstartingvalues = true;

    distr_gumbel2_sigma2s.push_back(DISTR_gumbel2_sigma2(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gumbel2_sigma2s[distr_gumbel2_sigma2s.size()-1];
    equations[modnr].pathd = "";

    }
//---------------------------- END: gumbel2_sigma2 -------------------------------


//------------------------------- gumbel2_mu ------------------------------------
  else if (family.getvalue() == "gumbel2" && equationtype.getvalue()=="mu")
    {
    if(generaloptions.copula)
      {
      mainequation=false;
      }
    else
      {
      mainequation=true;
      }


    computemodeforstartingvalues = true;

    distr_gumbel2_mus.push_back(DISTR_gumbel2_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gumbel2_mus[distr_gumbel2_mus.size()-1];
    equations[modnr].pathd = "";

    if (distr_gumbel2_sigma2s.size() != distr_gumbel2_mus.size())
      {
      outerror("ERROR: Equation for sigma2 is missing");
      return true;
      }
    else
      {
      distr_gumbel2_sigma2s[distr_gumbel2_sigma2s.size()-1].distrp.push_back
      (&distr_gumbel2_mus[distr_gumbel2_mus.size()-1]);

      distr_gumbel2_mus[distr_gumbel2_mus.size()-1].distrp.push_back
      (&distr_gumbel2_sigma2s[distr_gumbel2_sigma2s.size()-1]);

      countmarginal += 1;
      }
    if(generaloptions.copula){}
    else
      {
      predict_mult_distrs.push_back(&distr_gumbel2_sigma2s[distr_gumbel2_sigma2s.size()-1]);
      predict_mult_distrs.push_back(&distr_gumbel2_mus[distr_gumbel2_mus.size()-1]);
      }

    if(generaloptions.copula)
      {
      if (countmarginal == 1)
        {
        distr_gumbel2_mus[distr_gumbel2_mus.size()-1].set_copulapos(0);
        distr_gumbel2_sigma2s[distr_gumbel2_sigma2s.size()-1].set_copulapos(0);
        }
      else if(countmarginal == 2)
        {
        distr_gumbel2_mus[distr_gumbel2_mus.size()-1].set_copulapos(1);
        distr_gumbel2_sigma2s[distr_gumbel2_sigma2s.size()-1].set_copulapos(1);
        }
      else
        {
        outerror("ERROR: currently only bivariate copula models implemented");
        return true;
        }
      }

    }
//------------------------------- END: gumbel2_mu -------------------------------


//-------------------------------- gamma_sigma ---------------------------------
  else if ((family.getvalue() == "gamma" || family.getvalue() == "za_gamma")
           && equationtype.getvalue()=="sigma")
    {

    computemodeforstartingvalues = true;

    datamatrix dnew = D.getCol(0);
    if(family.getvalue() == "za_gamma")
      {
      if(w.rows()==1)
        w = datamatrix(dnew.rows(),1,1);
      for(unsigned i=0; i<dnew.rows(); i++)
        {
        if(D(i,0)==0)
          {
          w(i,0) = 0;
          dnew(i,0) = 0.1;
          }
        }
      }

    distr_gamma_sigmas.push_back(DISTR_gamma_sigma(&generaloptions,dnew,w));

    equations[modnr].distrp = &distr_gamma_sigmas[distr_gamma_sigmas.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_gamma_sigmas[distr_gamma_sigmas.size()-1]);

    }
//---------------------------- END: gamma_sigma -------------------------------


//------------------------------- gamma_mu ------------------------------------
  else if ((family.getvalue() == "gamma" || family.getvalue() == "za_gamma")
           && ((equationtype.getvalue()=="mu"))
          )
    {
    datamatrix dnew = D.getCol(0);
    if(generaloptions.copula)
      {
      mainequation=false;
      }
    else
      {
      if(family.getvalue() == "gamma")
        mainequation=true;
      else
        {
        mainequation=false;
        if(w.rows()==1)
          w = datamatrix(dnew.rows(),1,1);
        for(unsigned i=0; i<dnew.rows(); i++)
          {
          if(D(i,0)==0)
            {
            w(i,0) = 0;
            dnew(i,0) = 0.1;
            }
          }
        }
      }


    computemodeforstartingvalues = true;

    distr_gamma_mus.push_back(DISTR_gamma_mu(&generaloptions,dnew,w));

    equations[modnr].distrp = &distr_gamma_mus[distr_gamma_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_gamma_mus[distr_gamma_mus.size()-1]);

    if (distr_gamma_sigmas.size() != distr_gamma_mus.size())
      {
      outerror("ERROR: Equation for sigma is missing");
      return true;
      }
    else
      {
      distr_gamma_sigmas[distr_gamma_sigmas.size()-1].distrp.push_back
      (&distr_gamma_mus[distr_gamma_mus.size()-1]);

      distr_gamma_mus[distr_gamma_mus.size()-1].distrp.push_back
      (&distr_gamma_sigmas[distr_gamma_sigmas.size()-1]);

      countmarginal += 1;
      }
    if(generaloptions.copula){}
    else
      {
      predict_mult_distrs.push_back(&distr_gamma_sigmas[distr_gumbel_sigmas.size()-1]);
      predict_mult_distrs.push_back(&distr_gamma_mus[distr_gumbel_mus.size()-1]);
      }

    if(generaloptions.copula)
      {
      if (countmarginal == 1)
        {
        distr_gamma_mus[distr_gamma_mus.size()-1].set_copulapos(0);
        distr_gamma_sigmas[distr_gamma_sigmas.size()-1].set_copulapos(0);
        }
      else if(countmarginal == 2)
        {
        distr_gamma_mus[distr_gamma_mus.size()-1].set_copulapos(1);
        distr_gamma_sigmas[distr_gamma_sigmas.size()-1].set_copulapos(1);
        }
      else
        {
        outerror("ERROR: currently only bivariate copula models implemented");
        return true;
        }
      }
    }
//------------------------------- END: gamma_mu -------------------------------

//-------------------------------- invgaussian sigma2 ---------------------------------
  else if ((family.getvalue() == "invgaussian" || family.getvalue() == "za_invgaussian")
           && equationtype.getvalue()=="sigma2")
    {

    computemodeforstartingvalues = true;

    datamatrix dnew = D.getCol(0);
    if(family.getvalue() == "za_invgaussian")
      {
      if(w.rows()==1)
        w = datamatrix(dnew.rows(),1,1);
      for(unsigned i=0; i<dnew.rows(); i++)
        {
        if(D(i,0)==0)
          {
          w(i,0) = 0;
          dnew(i,0) = 0.1;
          }
        }
      }


    computemodeforstartingvalues = true;

    distr_invgaussian_sigma2s.push_back(DISTR_invgaussian_sigma2(&generaloptions,dnew,w));

    equations[modnr].distrp = &distr_invgaussian_sigma2s[distr_invgaussian_sigma2s.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_invgaussian_sigma2s[distr_invgaussian_sigma2s.size()-1]);

    }
//---------------------------- END: invgaussian sigma2 -------------------------------

//------------------------------- invgaussian mu ------------------------------------
  else if ((family.getvalue() == "invgaussian" || family.getvalue() == "za_invgaussian")
           && ((equationtype.getvalue()=="mu"))
          )
    {
    datamatrix dnew= D.getCol(0);

    if(family.getvalue() == "invgaussian")
      mainequation=true;
    else
      {
      mainequation=false;
      if(w.rows()==1)
        w = datamatrix(dnew.rows(),1,1);
      for(unsigned i=0; i<dnew.rows(); i++)
        {
        if(D(i,0)==0)
          {
          w(i,0) = 0;
          dnew(i,0) = 0.1;
          }
        }
      }

    computemodeforstartingvalues = true;

    distr_invgaussian_mus.push_back(DISTR_invgaussian_mu(&generaloptions,dnew,w));

    equations[modnr].distrp = &distr_invgaussian_mus[distr_invgaussian_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_invgaussian_mus[distr_invgaussian_mus.size()-1]);

    if (distr_invgaussian_sigma2s.size() != 1)
      {
      outerror("ERROR: Equation for sigma2 is missing");
      return true;
      }
    else
      {
      distr_invgaussian_sigma2s[distr_invgaussian_sigma2s.size()-1].distrp.push_back
      (&distr_invgaussian_mus[distr_invgaussian_mus.size()-1]);

      distr_invgaussian_mus[distr_invgaussian_mus.size()-1].distrp.push_back
      (&distr_invgaussian_sigma2s[distr_invgaussian_sigma2s.size()-1]);
      }

    }
//------------------------------- END: invgaussian mu -------------------------------


//-------------------------------- truncnormal2 sigma ---------------------------------
  else if (family.getvalue() == "truncnormal2" && equationtype.getvalue()=="sigma")
    {

    computemodeforstartingvalues = true;

    distr_truncnormal2_sigmas.push_back(DISTR_truncnormal2_sigma(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_truncnormal2_sigmas[distr_truncnormal2_sigmas.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_truncnormal2_sigmas[distr_truncnormal2_sigmas.size()-1]);

    }
//---------------------------- END: truncnormal2 sigma -------------------------------

//------------------------------- truncnormal2 mu ------------------------------------
  else if ((family.getvalue() == "truncnormal2") && (equationtype.getvalue()=="mu"))

    {

     mainequation=true;

    computemodeforstartingvalues = true;

    distr_truncnormal2_mus.push_back(DISTR_truncnormal2_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_truncnormal2_mus[distr_truncnormal2_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_truncnormal2_mus[distr_truncnormal2_mus.size()-1]);

    if (distr_truncnormal2_sigmas.size() != 1)
      {
      outerror("ERROR: Equation for sigma is missing");
      return true;
      }
    else
      {
      distr_truncnormal2_sigmas[distr_truncnormal2_sigmas.size()-1].distrp.push_back
      (&distr_truncnormal2_mus[distr_truncnormal2_mus.size()-1]);

      distr_truncnormal2_mus[distr_truncnormal2_mus.size()-1].distrp.push_back
      (&distr_truncnormal2_sigmas[distr_truncnormal2_sigmas.size()-1]);
      }

    }
//------------------------------- END: truncnormal2 mu -------------------------------


//-------------------------------- normal2 sigma ---------------------------------
  else if (family.getvalue() == "normal2" && equationtype.getvalue()=="sigma")
    {

    computemodeforstartingvalues = true;

    distr_normal2_sigmas.push_back(DISTR_normal2_sigma(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_normal2_sigmas[distr_normal2_sigmas.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_normal2_sigmas[distr_normal2_sigmas.size()-1]);

    }
//---------------------------- END: normal2 sigma -------------------------------

//------------------------------- normal2 mu ------------------------------------
  else if ((family.getvalue() == "normal2") && (equationtype.getvalue()=="mu"))

    {

     mainequation=true;

    computemodeforstartingvalues = true;

    distr_normal2_mus.push_back(DISTR_normal2_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_normal2_mus[distr_normal2_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_normal2_mus[distr_normal2_mus.size()-1]);

    if (distr_normal2_sigmas.size() != 1)
      {
      outerror("ERROR: Equation for sigma is missing");
      return true;
      }
    else
      {
      distr_normal2_sigmas[distr_normal2_sigmas.size()-1].distrp.push_back
      (&distr_normal2_mus[distr_normal2_mus.size()-1]);

      distr_normal2_mus[distr_normal2_mus.size()-1].distrp.push_back
      (&distr_normal2_sigmas[distr_normal2_sigmas.size()-1]);
      }

    }
//------------------------------- END: normal2 mu -------------------------------

//-------------------------------- beta sigma2 ---------------------------------
  else if (family.getvalue() == "beta" && equationtype.getvalue()=="sigma2")
    {

    computemodeforstartingvalues = true;

    distr_beta_sigma2s.push_back(DISTR_beta_sigma2(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_beta_sigma2s[distr_beta_sigma2s.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_beta_sigma2s[distr_beta_sigma2s.size()-1]);

    }
//---------------------------- END: beta sigma2 -------------------------------

//------------------------------- beta mu ------------------------------------
  else if ( (family.getvalue() == "beta") && (equationtype.getvalue()=="mu"))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_beta_mus.push_back(DISTR_beta_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_beta_mus[distr_beta_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_beta_mus[distr_beta_mus.size()-1]);

    if (distr_beta_sigma2s.size() != 1)
      {
      outerror("ERROR: Equation for sigma2 is missing");
      return true;
      }
    else
      {
      distr_beta_sigma2s[distr_beta_sigma2s.size()-1].distrp.push_back
       (&distr_beta_mus[distr_beta_mus.size()-1]);

      distr_beta_mus[distr_beta_mus.size()-1].distrp.push_back
       (&distr_beta_sigma2s[distr_beta_sigma2s.size()-1]);
      }

    }
//------------------------------- END: beta mu -------------------------------

//----------------------------- cloglog -------------------------------
  else if (family.getvalue() == "binomial_cloglog" || family.getvalue() == "za_binomial_cloglog")//&& equationtype.getvalue()=="mean"|| (equationtype.getvalue()=="meanservant"))
    {
    datamatrix dnew = D.getCol(0);
    if(family.getvalue() == "binomial_cloglog")
      mainequation=true;
    else
      {
      mainequation=false;
      for(unsigned i=0; i<dnew.rows(); i++)
        {
        if(dnew(i,0)>0)
          dnew(i,0) = 1;
        }
      }

    computemodeforstartingvalues = true;

    distr_cloglogs.push_back(DISTR_cloglog(&generaloptions,dnew,w));

    equations[modnr].distrp = &distr_cloglogs[distr_cloglogs.size()-1];
    equations[modnr].pathd = "";
    equations[modnr].equationtype = "pi";

    predict_mult_distrs.push_back(&distr_cloglogs[distr_cloglogs.size()-1]);

    }
//-------------------------- END: cloglog -----------------------------

//-------------------------------- claytoncopula2_normal sigma2 ---------------------------------
  else if (family.getvalue() == "claytoncopula2_normal_sigma2" && equationtype.getvalue()=="sigma2")
    {

	unsigned pos;
      if (distr_claytoncopula2_normal_sigma2s.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_claytoncopula2_normal_sigma2s.push_back(DISTR_claytoncopula2_normal_sigma2(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1];
    equations[modnr].pathd = "";

/*	if (distr_claytoncopula2_normal_sigma2s.size() == 2)
       {
       distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2].response2 = distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1].response;
       distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1].response2 = distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2].response;
       }*/


    }
//---------------------------- END: claytoncopula2_normal sigma2 -------------------------------

//------------------------------- claytoncopula2_normal mu ------------------------------------
  else if ((family.getvalue() == "claytoncopula2_normal_mu") && (equationtype.getvalue()=="mu"))
    {

	unsigned pos;
      if (distr_claytoncopula2_normal_mus.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_claytoncopula2_normal_mus.push_back(DISTR_claytoncopula2_normal_mu(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1];
    equations[modnr].pathd = "";

	/*if (distr_claytoncopula2_normal_mus.size() == 2)
       {
       distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2].response2 = distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1].response;
       distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1].response2 = distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2].response;
       }*/

    }
//------------------------------- END: claytoncopula2_normal mu -------------------------------

//----------------------------- claytoncopula2_rho ----------------------
  else if (((family.getvalue() == "claytoncopula2_normal_rho")||(family.getvalue() == "claytoncopula2_rho")) && ((equationtype.getvalue()=="rho")))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_claytoncopula2_rhos.push_back(DISTR_claytoncopula2_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1];
    equations[modnr].pathd = "";

    if (((distr_normal_mus.size() != 2) || (distr_normal_sigma2s.size() != 2)) && ((distr_claytoncopula2_normal_mus.size() != 2) || (distr_claytoncopula2_normal_sigma2s.size() != 2)))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }
		if ((distr_normal_mus.size() == 2))
		{
            predict_mult_distrs.push_back(&distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);
            predict_mult_distrs.push_back(&distr_normal_mus[distr_normal_mus.size()-2]);
            predict_mult_distrs.push_back(&distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);
            predict_mult_distrs.push_back(&distr_normal_mus[distr_normal_mus.size()-1]);
            predict_mult_distrs.push_back(&distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1]);

            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].response2 = distr_normal_mus[distr_normal_mus.size()-2].response;
          //  distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].response = distr_normal_mus[distr_normal_mus.size()-1].response;

            distr_normal_sigma2s[distr_normal_sigma2s.size()-2].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-2]);

            distr_normal_sigma2s[distr_normal_sigma2s.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-1]);

            distr_normal_mus[distr_normal_mus.size()-2].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);

            distr_normal_mus[distr_normal_mus.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);

            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);

            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-2]);

            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);


            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-1]);
		}
		else
		{
            predict_mult_distrs.push_back(&distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2]);
            predict_mult_distrs.push_back(&distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2]);
            predict_mult_distrs.push_back(&distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1]);
            predict_mult_distrs.push_back(&distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1]);
            predict_mult_distrs.push_back(&distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1]);

            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].response2 = distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2].response;
            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].response = distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1].response;
            distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2].response2 = distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1].response;
            distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1].response2 = distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2].response;
            distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2].response2 = distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1].response;
            distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1].response2 = distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2].response;


            distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1]);

            distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1]);

            distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2]);

            distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1]);

            distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1]);

            distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2]);

            distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1]);

            distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2]);

            distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2].distrp.push_back(
            &distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1]);

            distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2].distrp.push_back(
            &distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1]);

            distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2].distrp.push_back(
            &distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2]);

            distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2].distrp.push_back(
            &distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1]);

            distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1].distrp.push_back(
            &distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1]);

            distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1].distrp.push_back(
            &distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2]);

            distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1].distrp.push_back(
            &distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1]);

            distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1].distrp.push_back(
            &distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2]);

            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].distrp.push_back(
            &distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-2]);

            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].distrp.push_back(
            &distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-2]);

            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].distrp.push_back(
            &distr_claytoncopula2_normal_sigma2s[distr_claytoncopula2_normal_sigma2s.size()-1]);

            distr_claytoncopula2_rhos[distr_claytoncopula2_rhos.size()-1].distrp.push_back(
            &distr_claytoncopula2_normal_mus[distr_claytoncopula2_normal_mus.size()-1]);
		}

    }
//-------------------------- END: claytoncopula2_rho ---------------------


//----------------------------- claytoncopula_rho ----------------------
  else if ( (family.getvalue() == "claytoncopula") && ((equationtype.getvalue()=="rho")))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_claytoncopula_rhos.push_back(DISTR_claytoncopula_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_claytoncopula_rhos[distr_claytoncopula_rhos.size()-1];
    equations[modnr].pathd = "";

    if ((distr_copulas.size() != 2))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }

            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-2]);
            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-1]);
			predict_mult_distrs.push_back(&distr_claytoncopula_rhos[distr_claytoncopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-2].response2 = distr_copulas[distr_copulas.size()-1].response;
            distr_copulas[distr_copulas.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;
            distr_claytoncopula_rhos[distr_claytoncopula_rhos.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;
          //  distr_claytoncopula_rhos[distr_claytoncopula_rhos.size()-1].response = distr_copulas[distr_copulas.size()-1].response;

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_claytoncopula_rhos[distr_claytoncopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_claytoncopula_rhos[distr_claytoncopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

			distr_claytoncopula_rhos[distr_claytoncopula_rhos.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

            distr_claytoncopula_rhos[distr_claytoncopula_rhos.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);

    }
//-------------------------- END: claytoncopula_rho ---------------------


//-------------------------------- gumbelcopula2_normal sigma2 ---------------------------------
  else if (family.getvalue() == "gumbelcopula2_normal_sigma2" && equationtype.getvalue()=="sigma2")
    {

	unsigned pos;
      if (distr_gumbelcopula2_normal_sigma2s.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_gumbelcopula2_normal_sigma2s.push_back(DISTR_gumbelcopula2_normal_sigma2(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1];
    equations[modnr].pathd = "";

/*	if (distr_gumbelcopula2_normal_sigma2s.size() == 2)
       {
       distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2].response2 = distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1].response;
       distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1].response2 = distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2].response;
       }*/


    }
//---------------------------- END: gumbelcopula2_normal sigma2 -------------------------------

//------------------------------- gumbelcopula2_normal mu ------------------------------------
  else if ((family.getvalue() == "gumbelcopula2_normal_mu") && equationtype.getvalue()=="mu")
    {

	unsigned pos;
      if (distr_gumbelcopula2_normal_mus.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_gumbelcopula2_normal_mus.push_back(DISTR_gumbelcopula2_normal_mu(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1];
    equations[modnr].pathd = "";

	/*if (distr_gumbelcopula2_normal_mus.size() == 2)
       {
       distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].response2 = distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].response;
       distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].response2 = distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].response;
       }*/

    }
//------------------------------- END: gumbelcopula2_normal mu -------------------------------


//-------------------------------- gumbelcopula2_normal sigma2_2 ---------------------------------
  else if (family.getvalue() == "gumbelcopula2_normal" && equationtype.getvalue()=="sigma2")
    {

	unsigned pos;
      if (distr_gumbelcopula2_normal_sigma2_2s.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_gumbelcopula2_normal_sigma2_2s.push_back(DISTR_gumbelcopula2_normal_sigma2_2(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_gumbelcopula2_normal_sigma2_2s[distr_gumbelcopula2_normal_sigma2_2s.size()-1];
    equations[modnr].pathd = "";

/*	if (distr_gumbelcopula2_normal_sigma2s.size() == 2)
       {
       distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2].response2 = distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1].response;
       distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1].response2 = distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2].response;
       }*/


    }
//---------------------------- END: gumbelcopula2_normal sigma2_2 -------------------------------

//------------------------------- gumbelcopula2_normal mu_2 ------------------------------------
  else if ((family.getvalue() == "gumbelcopula2_normal") && equationtype.getvalue()=="mu")
    {

	unsigned pos;
      if (distr_gumbelcopula2_normal_mu_2s.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_gumbelcopula2_normal_mu_2s.push_back(DISTR_gumbelcopula2_normal_mu_2(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-1];
    equations[modnr].pathd = "";

	/*if (distr_gumbelcopula2_normal_mus.size() == 2)
       {
       distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].response2 = distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].response;
       distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].response2 = distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].response;
       }*/

    }
//------------------------------- END: gumbelcopula2_normal_mu_2 -------------------------------

//----------------------------- gumbelcopula2_rho ----------------------
  else if ((family.getvalue() == "gumbelcopula2_normal_rho") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_gumbelcopula2_rhos.push_back(DISTR_gumbelcopula2_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1];
    equations[modnr].pathd = "";

    if (((distr_gumbelcopula2_normal_mu_2s.size() != 2) && (distr_gumbelcopula2_normal_mus.size() != 2)))
       {
       outerror("ERROR: Two equations for marginal distributions of mu required");
       return true;
       }
	if (((distr_gumbelcopula2_normal_sigma2_2s.size() != 2) && (distr_gumbelcopula2_normal_sigma2s.size() != 2)))
       {
       outerror("ERROR: Two equations for marginal distributions of sigma2 required");
       return true;
       }
		if ((distr_gumbelcopula2_normal_mu_2s.size() == 2))
		{
            predict_mult_distrs.push_back(&distr_gumbelcopula2_normal_sigma2_2s[distr_gumbelcopula2_normal_sigma2_2s.size()-2]);
            predict_mult_distrs.push_back(&distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-2]);
            predict_mult_distrs.push_back(&distr_gumbelcopula2_normal_sigma2_2s[distr_gumbelcopula2_normal_sigma2_2s.size()-1]);
            predict_mult_distrs.push_back(&distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-1]);
            predict_mult_distrs.push_back(&distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].response2 = distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-2].response;
         //   distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].response = distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-1].response;

            distr_gumbelcopula2_normal_sigma2_2s[distr_gumbelcopula2_normal_sigma2_2s.size()-2].distrp.push_back(
            &distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-2]);

            distr_gumbelcopula2_normal_sigma2_2s[distr_gumbelcopula2_normal_sigma2_2s.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-1]);

            distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-2].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2_2s[distr_gumbelcopula2_normal_sigma2_2s.size()-2]);

            distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2_2s[distr_gumbelcopula2_normal_sigma2_2s.size()-1]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2_2s[distr_gumbelcopula2_normal_sigma2_2s.size()-2]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-2]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2_2s[distr_gumbelcopula2_normal_sigma2_2s.size()-1]);


            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_mu_2s[distr_gumbelcopula2_normal_mu_2s.size()-1]);
		}
		else
		{
            predict_mult_distrs.push_back(&distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2]);
            predict_mult_distrs.push_back(&distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2]);
            predict_mult_distrs.push_back(&distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1]);
            predict_mult_distrs.push_back(&distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1]);
            predict_mult_distrs.push_back(&distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].response2 = distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].response;
            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].response = distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].response;
            distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2].response2 = distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1].response;
            distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1].response2 = distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2].response;
            distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].response2 = distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].response;
            distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].response2 = distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].response;


            distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1]);

            distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1]);

            distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2]);

            distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1]);

            distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1]);

            distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2]);

            distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1]);

            distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2]);

            distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1]);

            distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1]);

            distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2]);

            distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1]);

            distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1]);

            distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2]);

            distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1]);

            distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-2]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-2]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_sigma2s[distr_gumbelcopula2_normal_sigma2s.size()-1]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_gumbelcopula2_normal_mus[distr_gumbelcopula2_normal_mus.size()-1]);
		}

    }
//-------------------------- END: gumbelcopula2_rho ---------------------



//----------------------------- gumbelcopula_rho ----------------------
  else if ((family.getvalue() == "gumbelcopula") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_gumbelcopula_rhos.push_back(DISTR_gumbelcopula_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gumbelcopula_rhos[distr_gumbelcopula_rhos.size()-1];
    equations[modnr].pathd = "";

    if ((distr_copulas.size() != 2))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }

            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-2]);
            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-1]);
			predict_mult_distrs.push_back(&distr_gumbelcopula_rhos[distr_gumbelcopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-2].response2 = distr_copulas[distr_copulas.size()-1].response;
            distr_copulas[distr_copulas.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;
            distr_gumbelcopula_rhos[distr_gumbelcopula_rhos.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;
      //      distr_gumbelcopula_rhos[distr_gumbelcopula_rhos.size()-1].response = distr_copulas[distr_copulas.size()-1].response;

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_gumbelcopula_rhos[distr_gumbelcopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_gumbelcopula_rhos[distr_gumbelcopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

			distr_gumbelcopula_rhos[distr_gumbelcopula_rhos.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

            distr_gumbelcopula_rhos[distr_gumbelcopula_rhos.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);

    }
//-------------------------- END: gumbelcopula_rho ---------------------

//----------------------------- gumbelcopula2_rho ----------------------
  else if ((family.getvalue() == "gumbelcopula2_rho") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_gumbelcopula2_rhos.push_back(DISTR_gumbelcopula2_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1];
    equations[modnr].pathd = "";

    if ((distr_normal_mus.size() != 2) || (distr_normal_sigma2s.size() != 2))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }

            predict_mult_distrs.push_back(&distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);
            predict_mult_distrs.push_back(&distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);
            predict_mult_distrs.push_back(&distr_normal_mus[distr_normal_mus.size()-2]);
            predict_mult_distrs.push_back(&distr_normal_mus[distr_normal_mus.size()-1]);
            predict_mult_distrs.push_back(&distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].response2 = distr_normal_mus[distr_normal_mus.size()-1].response;

            distr_normal_sigma2s[distr_normal_sigma2s.size()-2].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-2]);

            distr_normal_sigma2s[distr_normal_sigma2s.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-1]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-2]);

            distr_gumbelcopula2_rhos[distr_gumbelcopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-1]);

    }
//-------------------------- END: gumbelcopula2_rho ---------------------


//-------------------------------- gaussiancopula_dagum_p ---------------------------------
  else if (family.getvalue() == "gaussiancopula_dagum" && equationtype.getvalue()=="p")
    {

	unsigned pos;
      if (distr_gaussiancopula_dagum_ps.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_gaussiancopula_dagum_ps.push_back(DISTR_gaussiancopula_dagum_p(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1];
    equations[modnr].pathd = "";

/*	if (distr_gaussiancopula_dagum_bs.size() == 2)
       {
       distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].response2 = distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].response;
       distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].response2 = distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].response;
       }*/


    }
//---------------------------- END: gaussiancopula_dagum_p -------------------------------

//-------------------------------- gaussiancopula_dagum_b ---------------------------------
  else if (family.getvalue() == "gaussiancopula_dagum" && equationtype.getvalue()=="b")
    {

	unsigned pos;
      if (distr_gaussiancopula_dagum_bs.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_gaussiancopula_dagum_bs.push_back(DISTR_gaussiancopula_dagum_b(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1];
    equations[modnr].pathd = "";

/*	if (distr_gaussiancopula_dagum_bs.size() == 2)
       {
       distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].response2 = distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].response;
       distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].response2 = distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].response;
       }*/


    }
//---------------------------- END: gaussiancopula_dagum_b -------------------------------

//------------------------------- gaussiancopula_dagum_a ------------------------------------
  else if ((family.getvalue() == "gaussiancopula_dagum") && equationtype.getvalue()=="a")
    {

	unsigned pos;
      if (distr_gaussiancopula_dagum_as.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_gaussiancopula_dagum_as.push_back(DISTR_gaussiancopula_dagum_a(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1];
    equations[modnr].pathd = "";

	/*if (distr_gaussiancopula_dagum_as.size() == 2)
       {
       distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].response2 = distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].response;
       distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].response2 = distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].response;
       }*/

    }
//------------------------------- END: ggaussiancopula_dagum_a -------------------------------



//----------------------------- gaussiancopula_dagum_rho ----------------------
  else if (family.getvalue() == "gaussiancopula_dagum" && equationtype.getvalue()=="rho")
    {
	  mainequation = true;
    computemodeforstartingvalues = true;

    distr_gaussiancopula_dagum_rhos.push_back(DISTR_gaussiancopula_dagum_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1];
    equations[modnr].pathd = "";

    if (((distr_dagum_ps.size() != 2) || (distr_dagum_as.size() != 2) || (distr_dagum_bs.size() != 2)) && ((distr_gaussiancopula_dagum_ps.size() != 2) || (distr_gaussiancopula_dagum_as.size() != 2) || (distr_gaussiancopula_dagum_bs.size() != 2)))
       {
       outerror("ERROR: three equations for marginal distributions required");
       return true;
       }
		if ((distr_dagum_as.size() == 2))
		{
		        predict_mult_distrs.push_back(&distr_dagum_ps[distr_dagum_ps.size()-2]);
            predict_mult_distrs.push_back(&distr_dagum_bs[distr_dagum_bs.size()-2]);
            predict_mult_distrs.push_back(&distr_dagum_as[distr_dagum_as.size()-2]);
            predict_mult_distrs.push_back(&distr_dagum_ps[distr_dagum_ps.size()-1]);
            predict_mult_distrs.push_back(&distr_dagum_bs[distr_dagum_bs.size()-1]);
            predict_mult_distrs.push_back(&distr_dagum_as[distr_dagum_as.size()-1]);
            predict_mult_distrs.push_back(&distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].response2 = distr_dagum_as[distr_dagum_as.size()-2].response;
           // distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].response = distr_dagum_as[distr_dagum_as.size()-1].response;

            distr_dagum_ps[distr_dagum_ps.size()-2].distrp.push_back(
            &distr_dagum_bs[distr_dagum_bs.size()-2]);

            distr_dagum_ps[distr_dagum_ps.size()-2].distrp.push_back(
            &distr_dagum_as[distr_dagum_as.size()-2]);

            distr_dagum_ps[distr_dagum_ps.size()-1].distrp.push_back(
            &distr_dagum_bs[distr_dagum_bs.size()-1]);

            distr_dagum_ps[distr_dagum_ps.size()-1].distrp.push_back(
            &distr_dagum_as[distr_dagum_as.size()-1]);

            distr_dagum_bs[distr_dagum_bs.size()-2].distrp.push_back(
            &distr_dagum_ps[distr_dagum_ps.size()-2]);

            distr_dagum_bs[distr_dagum_bs.size()-2].distrp.push_back(
            &distr_dagum_as[distr_dagum_as.size()-2]);

            distr_dagum_bs[distr_dagum_bs.size()-1].distrp.push_back(
            &distr_dagum_ps[distr_dagum_ps.size()-1]);

            distr_dagum_bs[distr_dagum_bs.size()-1].distrp.push_back(
            &distr_dagum_as[distr_dagum_as.size()-1]);

            distr_dagum_as[distr_dagum_as.size()-2].distrp.push_back(
            &distr_dagum_ps[distr_dagum_ps.size()-2]);

            distr_dagum_as[distr_dagum_as.size()-2].distrp.push_back(
            &distr_dagum_bs[distr_dagum_bs.size()-2]);

            distr_dagum_as[distr_dagum_as.size()-1].distrp.push_back(
            &distr_dagum_ps[distr_dagum_ps.size()-1]);

            distr_dagum_as[distr_dagum_as.size()-1].distrp.push_back(
            &distr_dagum_bs[distr_dagum_bs.size()-1]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_dagum_ps[distr_dagum_ps.size()-2]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_dagum_bs[distr_dagum_bs.size()-2]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_dagum_as[distr_dagum_as.size()-2]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_dagum_ps[distr_dagum_ps.size()-1]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_dagum_bs[distr_dagum_bs.size()-1]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_dagum_as[distr_dagum_as.size()-1]);
		}
		else
		{
            predict_mult_distrs.push_back(&distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2]);
            predict_mult_distrs.push_back(&distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2]);
            predict_mult_distrs.push_back(&distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2]);
            predict_mult_distrs.push_back(&distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1]);
            predict_mult_distrs.push_back(&distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1]);
            predict_mult_distrs.push_back(&distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1]);
            predict_mult_distrs.push_back(&distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].response2 = distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].response;
            //
            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].response = distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].response;
            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2].response2 = distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1].response;
            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1].response2 = distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2].response;
            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].response2 = distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].response;
            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].response2 = distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].response;
            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].response2 = distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].response;
            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].response2 = distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].response;

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1]);

            distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1]);

            distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1]);

            distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-2]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-2]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-2]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_ps[distr_gaussiancopula_dagum_ps.size()-1]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_bs[distr_gaussiancopula_dagum_bs.size()-1]);

            distr_gaussiancopula_dagum_rhos[distr_gaussiancopula_dagum_rhos.size()-1].distrp.push_back(
            &distr_gaussiancopula_dagum_as[distr_gaussiancopula_dagum_as.size()-1]);
		}

    }
//-------------------------- END: gaussiancopula_dagum_rho ---------------------



//----------------------------- gaussiancopula_rho ----------------------
  else if ((family.getvalue() == "gaussiancopula") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_gaussiancopula_rhos.push_back(DISTR_gaussiancopula_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_rhos[distr_gaussiancopula_rhos.size()-1];
    equations[modnr].pathd = "";

    if ((distr_copulas.size() != 2))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }

            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-2]);
            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-1]);
            predict_mult_distrs.push_back(&distr_gaussiancopula_rhos[distr_gaussiancopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-2].response2 = distr_copulas[distr_copulas.size()-1].response;
            distr_copulas[distr_copulas.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;
            distr_gaussiancopula_rhos[distr_gaussiancopula_rhos.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_gaussiancopula_rhos[distr_gaussiancopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_gaussiancopula_rhos[distr_gaussiancopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

            distr_gaussiancopula_rhos[distr_gaussiancopula_rhos.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

            distr_gaussiancopula_rhos[distr_gaussiancopula_rhos.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);


    }
//-------------------------- END: gaussiancopula_rho ---------------------

//----------------------------- gaussiancopula_rhofz ----------------------
  else if ((family.getvalue() == "gaussiancopula") && (equationtype.getvalue()=="rhofz"))
    {

    mainequation=true;
    computemodeforstartingvalues = true;

    distr_gaussiancopula_rhofzs.push_back(DISTR_gaussiancopula_rhofz(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_rhofzs[distr_gaussiancopula_rhofzs.size()-1];
    equations[modnr].pathd = "";

    if ((distr_copulas.size() != 2))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }

            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-2]);
            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-1]);
			predict_mult_distrs.push_back(&distr_gaussiancopula_rhofzs[distr_gaussiancopula_rhofzs.size()-1]);

            distr_copulas[distr_copulas.size()-2].response2 = distr_copulas[distr_copulas.size()-1].response;
            distr_copulas[distr_copulas.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;
            distr_gaussiancopula_rhofzs[distr_gaussiancopula_rhofzs.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_gaussiancopula_rhofzs[distr_gaussiancopula_rhofzs.size()-1]);

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_gaussiancopula_rhofzs[distr_gaussiancopula_rhofzs.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

			distr_gaussiancopula_rhofzs[distr_gaussiancopula_rhofzs.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

            distr_gaussiancopula_rhofzs[distr_gaussiancopula_rhofzs.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);


    }
//-------------------------- END: gaussiancopula_rhofz ---------------------

//----------------------------- gausscopula --------------------------------
  else if ((family.getvalue() == "gauss_copula") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_gausscopulas.push_back(DISTR_gausscopula(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gausscopulas[distr_gausscopulas.size()-1];
    equations[modnr].pathd = "";

    if ((countmarginal != 2))
      {
      outerror("ERROR: Two equations for marginal distributions required or family not implemented yet");
      return true;
      }
    if ((generaloptions.rotation != 0))
      {
      out("WARNING: No meaningful rotation for Gauss copula. Copula will not be rotated.\n\n");
      generaloptions.rotation = 0;
      }

    // Fill distrp of the copula with some starting values to make it a vector of length=2
    distr_gausscopulas[distr_gausscopulas.size()-1].distrp.push_back(equations[modnr-1].distrp);
    distr_gausscopulas[distr_gausscopulas.size()-1].distrp.push_back(equations[modnr-2].distrp);

		vector<DISTR*> pred_mult_margin1;
		vector<DISTR*> pred_mult_margin0;

    if(distr_binomialprobit_copulas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_binomialprobit_copulas.size();coi++)
        {
   //     predict_mult_distrs.push_back(&distr_binomialprobit_copulas[coi]);
        distr_binomialprobit_copulas[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);

        if(distr_binomialprobit_copulas[coi].get_copulapos()==0)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[0]= (&distr_binomialprobit_copulas[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response2 = distr_binomialprobit_copulas[coi].response;
          distr_binomialprobit_copulas[coi].responsecopmat = &distr_gausscopulas[distr_gausscopulas.size()-1].response2;

					pred_mult_margin0.push_back(&distr_binomialprobit_copulas[coi]);
					}

        else if(distr_binomialprobit_copulas[coi].get_copulapos()==1)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[1]= (&distr_binomialprobit_copulas[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response1 = distr_binomialprobit_copulas[coi].response;
          distr_binomialprobit_copulas[coi].responsecopmat = &distr_gausscopulas[distr_gausscopulas.size()-1].response1;

					pred_mult_margin1.push_back(&distr_binomialprobit_copulas[coi]);
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_weibull_lambdas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_weibull_lambdas.size();coi++)
        {
/*         predict_mult_distrs.push_back(&distr_weibull_alphas[coi]);
        predict_mult_distrs.push_back(&distr_weibull_lambdas[coi]); */

        //distr_gausscopulas[distr_gausscopulas.size()-1].distrp.push_back(&distr_weibull_lambdas[coi]);

        distr_weibull_lambdas[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);
        distr_weibull_alphas[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);

        if(distr_weibull_lambdas[coi].get_copulapos()==0)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[0]= (&distr_weibull_lambdas[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response2 = distr_weibull_lambdas[coi].response;

					pred_mult_margin0.push_back(&distr_weibull_alphas[coi]);
					pred_mult_margin0.push_back(&distr_weibull_lambdas[coi]);

          }
        else if(distr_weibull_lambdas[coi].get_copulapos()==1)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[1]= (&distr_weibull_lambdas[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response1 = distr_weibull_lambdas[coi].response;

					pred_mult_margin1.push_back(&distr_weibull_alphas[coi]);
					pred_mult_margin1.push_back(&distr_weibull_lambdas[coi]);
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_weibull2_lambdas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_weibull2_lambdas.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_weibull2_alphas[coi]);
        predict_mult_distrs.push_back(&distr_weibull2_lambdas[coi]);

        //distr_gausscopulas[distr_gausscopulas.size()-1].distrp.push_back(&distr_weibull2_lambdas[coi]);
        distr_weibull2_lambdas[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);
        distr_weibull2_alphas[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);

        if(distr_weibull2_lambdas[coi].get_copulapos()==0)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[0]= (&distr_weibull2_lambdas[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response2 = distr_weibull2_lambdas[coi].response;
          }
        else if(distr_weibull2_lambdas[coi].get_copulapos()==1)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[1]= (&distr_weibull2_lambdas[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response1 = distr_weibull2_lambdas[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_gumbel_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_gumbel_mus.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_gumbel_mus[coi]);
        predict_mult_distrs.push_back(&distr_gumbel_sigmas[coi]);

        distr_gumbel_mus[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);
        distr_gumbel_sigmas[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);

        if(distr_gumbel_mus[coi].get_copulapos()==0)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[0]= (&distr_gumbel_mus[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response2 = distr_gumbel_mus[coi].response;
          }
        else if(distr_gumbel_mus[coi].get_copulapos()==1)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[1]= (&distr_gumbel_mus[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response1 = distr_gumbel_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_gamma_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_gamma_mus.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_gamma_mus[coi]);
        predict_mult_distrs.push_back(&distr_gamma_sigmas[coi]);

        distr_gamma_mus[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);
        distr_gamma_sigmas[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);

        if(distr_gamma_mus[coi].get_copulapos()==0)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[0]= (&distr_gamma_mus[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response2 = distr_gamma_mus[coi].response;
          }
        else if(distr_gamma_mus[coi].get_copulapos()==1)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[1]= (&distr_gamma_mus[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response1 = distr_gamma_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_normal_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_normal_mus.size();coi++)
        {
        // predict_mult_distrs.push_back(&distr_normal_sigma2s[coi]);
        // predict_mult_distrs.push_back(&distr_normal_mus[coi]);

        //distr_gausscopulas[distr_gausscopulas.size()-1].distrp.push_back(&distr_normal_mus[coi]);

        distr_normal_mus[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);
        distr_normal_sigma2s[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);

        if(distr_normal_mus[coi].get_copulapos()==0)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[0]= (&distr_normal_mus[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response2 = distr_normal_mus[coi].response;

					pred_mult_margin0.push_back(&distr_normal_sigma2s[coi]);
					pred_mult_margin0.push_back(&distr_normal_mus[coi]);
          }
        else if(distr_normal_mus[coi].get_copulapos()==1)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[1]= (&distr_normal_mus[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response1 = distr_normal_mus[coi].response;

					pred_mult_margin1.push_back(&distr_normal_sigma2s[coi]);
					pred_mult_margin1.push_back(&distr_normal_mus[coi]);
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_dagum_as.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_dagum_as.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_dagum_ps[coi]);
        predict_mult_distrs.push_back(&distr_dagum_bs[coi]);
        predict_mult_distrs.push_back(&distr_dagum_as[coi]);

        //distr_gausscopulas[distr_gausscopulas.size()-1].distrp.push_back(&distr_dagum_as[coi]);

        distr_dagum_as[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);
        distr_dagum_bs[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);
        distr_dagum_ps[coi].distrcopulap.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);

        if(distr_dagum_as[coi].get_copulapos()==0)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[0]= (&distr_dagum_as[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response2 = distr_dagum_as[coi].response;
          }
        else if(distr_dagum_as[coi].get_copulapos()==1)
          {
          distr_gausscopulas[distr_gausscopulas.size()-1].distrp[1]= (&distr_dagum_as[coi]);
          distr_gausscopulas[distr_gausscopulas.size()-1].response1 = distr_dagum_as[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }

    // initialize copulaoffset for the second marginal (required for marginal deviance etc.)
    int off = distr_gausscopulas[distr_gausscopulas.size()-1].distrp[0]->distrp.size()+1;
    distr_gausscopulas[distr_gausscopulas.size()-1].distrp[1]->copulaoffset = off;
    for(unsigned k=0; k<distr_gausscopulas[distr_gausscopulas.size()-1].distrp[1]->distrp.size(); k++)
      distr_gausscopulas[distr_gausscopulas.size()-1].distrp[1]->distrp[k]->copulaoffset = off;

    if(distr_binomialprobit_copulas.size()>0)
      {
      if(setseed.getvalue() >= 0)
        srand(setseed.getvalue());

      distr_gausscopulas[0].update_end();
      int coi;
      for(coi=0;coi<distr_binomialprobit_copulas.size();coi++)
        {
        distr_binomialprobit_copulas[coi].update();
        }
      }

		cerr  << "pred_mult_margin0.size " << pred_mult_margin0.size() << "\n\n";
		cerr  << "pred_mult_margin1.size " << pred_mult_margin1.size() << "\n\n";
         // abort();
		int lengi;
		for(lengi=0;lengi<pred_mult_margin0.size();lengi++)
        {
					predict_mult_distrs.push_back(pred_mult_margin0[lengi]);
				}
		for(lengi=0;lengi<pred_mult_margin1.size();lengi++)
        {
					predict_mult_distrs.push_back(pred_mult_margin1[lengi]);
				}
    predict_mult_distrs.push_back(&distr_gausscopulas[distr_gausscopulas.size()-1]);
    }
//-------------------------- END: gausscopula---------------------------


//----------------------------- gausscopula2 --------------------------------
  else if ((family.getvalue() == "gauss_copula2") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_gausscopula2s.push_back(DISTR_gausscopula2(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gausscopula2s[distr_gausscopula2s.size()-1];
    equations[modnr].pathd = "";

    if ((countmarginal != 2))
      {
      outerror("ERROR: Two equations for marginal distributions required or family not implemented yet");
      return true;
      }

    // Fill distrp of the copula with some starting values to make it a vector of length=2
    distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp.push_back(equations[modnr-1].distrp);
    distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp.push_back(equations[modnr-2].distrp);

    if(distr_weibull_lambdas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_weibull_lambdas.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_weibull_alphas[coi]);
        predict_mult_distrs.push_back(&distr_weibull_lambdas[coi]);

        //distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp.push_back(&distr_weibull_lambdas[coi]);

        distr_weibull_lambdas[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);
        distr_weibull_alphas[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);

        if(distr_weibull_lambdas[coi].get_copulapos()==0)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[0]= (&distr_weibull_lambdas[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response2 = distr_weibull_lambdas[coi].response;
          }
        else if(distr_weibull_lambdas[coi].get_copulapos()==1)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[1]= (&distr_weibull_lambdas[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response1 = distr_weibull_lambdas[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_weibull2_lambdas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_weibull2_lambdas.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_weibull2_alphas[coi]);
        predict_mult_distrs.push_back(&distr_weibull2_lambdas[coi]);

        //distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp.push_back(&distr_weibull2_lambdas[coi]);

        distr_weibull2_lambdas[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);
        distr_weibull2_alphas[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);

        if(distr_weibull2_lambdas[coi].get_copulapos()==0)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[0]= (&distr_weibull2_lambdas[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response2 = distr_weibull2_lambdas[coi].response;
          }
        else if(distr_weibull2_lambdas[coi].get_copulapos()==1)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[1]= (&distr_weibull2_lambdas[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response1 = distr_weibull2_lambdas[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_gumbel_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_gumbel_mus.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_gumbel_mus[coi]);
        predict_mult_distrs.push_back(&distr_gumbel_sigmas[coi]);

        distr_gumbel_mus[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);
        distr_gumbel_sigmas[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);

        if(distr_gumbel_mus[coi].get_copulapos()==0)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[0]= (&distr_gumbel_mus[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response2 = distr_gumbel_mus[coi].response;
          }
        else if(distr_gumbel_mus[coi].get_copulapos()==1)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[1]= (&distr_gumbel_mus[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response1 = distr_gumbel_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_gamma_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_gamma_mus.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_gamma_mus[coi]);
        predict_mult_distrs.push_back(&distr_gamma_sigmas[coi]);

        distr_gamma_mus[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);
        distr_gamma_sigmas[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);

        if(distr_gamma_mus[coi].get_copulapos()==0)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[0]= (&distr_gamma_mus[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response2 = distr_gamma_mus[coi].response;
          }
        else if(distr_gamma_mus[coi].get_copulapos()==1)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[1]= (&distr_gamma_mus[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response1 = distr_gamma_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_normal_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_normal_mus.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_normal_sigma2s[coi]);
        predict_mult_distrs.push_back(&distr_normal_mus[coi]);

        //distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp.push_back(&distr_normal_mus[coi]);

        distr_normal_mus[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);
        distr_normal_sigma2s[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);

        if(distr_normal_mus[coi].get_copulapos()==0)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[0]= (&distr_normal_mus[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response2 = distr_normal_mus[coi].response;
          }
        else if(distr_normal_mus[coi].get_copulapos()==1)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[1]= (&distr_normal_mus[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response1 = distr_normal_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_dagum_as.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_dagum_as.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_dagum_ps[coi]);
        predict_mult_distrs.push_back(&distr_dagum_bs[coi]);
        predict_mult_distrs.push_back(&distr_dagum_as[coi]);

        //distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp.push_back(&distr_dagum_as[coi]);

        distr_dagum_as[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);
        distr_dagum_bs[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);
        distr_dagum_ps[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);

        if(distr_dagum_as[coi].get_copulapos()==0)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[0]= (&distr_dagum_as[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response2 = distr_dagum_as[coi].response;
          }
        else if(distr_dagum_as[coi].get_copulapos()==1)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[1]= (&distr_dagum_as[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response1 = distr_dagum_as[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_binomialprobit_copulas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_binomialprobit_copulas.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_binomialprobit_copulas[coi]);

        //distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp.push_back(&distr_binomialprobit_copulas[coi]);

        distr_binomialprobit_copulas[coi].distrcopulap.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);

        if(distr_binomialprobit_copulas[coi].get_copulapos()==0)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[0]= (&distr_binomialprobit_copulas[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response2 = distr_binomialprobit_copulas[coi].response;
          distr_binomialprobit_copulas[coi].responsecopmat = &distr_gausscopula2s[distr_gausscopula2s.size()-1].response2;
          }
        else if(distr_binomialprobit_copulas[coi].get_copulapos()==1)
          {
          distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[1]= (&distr_binomialprobit_copulas[coi]);
          distr_gausscopula2s[distr_gausscopula2s.size()-1].response1 = distr_binomialprobit_copulas[coi].response;
          distr_binomialprobit_copulas[coi].responsecopmat = &distr_gausscopula2s[distr_gausscopula2s.size()-1].response1;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }

    // initialize copulaoffset for the second marginal (required for marginal deviance etc.)
    int off = distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[0]->distrp.size()+1;
    distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[1]->copulaoffset = off;
    for(unsigned k=0; k<distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[1]->distrp.size(); k++)
      distr_gausscopula2s[distr_gausscopula2s.size()-1].distrp[1]->distrp[k]->copulaoffset = off;

    if(distr_binomialprobit_copulas.size()>0)
      {
      if(setseed.getvalue() >= 0)
        srand(setseed.getvalue());

      distr_gausscopula2s[0].update_end();

      int coi;
      for(coi=0;coi<distr_binomialprobit_copulas.size();coi++)
        {
        distr_binomialprobit_copulas[coi].update();
        }
      }
    //cout << "size predict mults: " << predict_mult_distrs.size() << endl;
    predict_mult_distrs.push_back(&distr_gausscopula2s[distr_gausscopula2s.size()-1]);
    }
//-------------------------- END: gausscopula2---------------------------



//----------------------------- clayton_copula --------------------------------
  else if ((family.getvalue() == "clayton_copula") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_clayton_copulas.push_back(DISTR_clayton_copula(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_clayton_copulas[distr_clayton_copulas.size()-1];
    equations[modnr].pathd = "";

    if ((countmarginal != 2))
      {
      outerror("ERROR: Two equations for marginal distributions required or family not implemented yet");
      return true;
      }

    // Fill distrp of the copula with some starting values to make it a vector of length=2
    distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp.push_back(equations[modnr-1].distrp);
    distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp.push_back(equations[modnr-2].distrp);

    if(distr_weibull_lambdas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_weibull_lambdas.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_weibull_alphas[coi]);
        predict_mult_distrs.push_back(&distr_weibull_lambdas[coi]);

        //distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp.push_back(&distr_weibull_lambdas[coi]);

        distr_weibull_lambdas[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);
        distr_weibull_alphas[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);

        if(distr_weibull_lambdas[coi].get_copulapos()==0)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[0]= (&distr_weibull_lambdas[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response2 = distr_weibull_lambdas[coi].response;
          }
        else if(distr_weibull_lambdas[coi].get_copulapos()==1)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]= (&distr_weibull_lambdas[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response1 = distr_weibull_lambdas[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_weibull2_lambdas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_weibull2_lambdas.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_weibull2_alphas[coi]);
        predict_mult_distrs.push_back(&distr_weibull2_lambdas[coi]);

        //distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp.push_back(&distr_weibull2_lambdas[coi]);

        distr_weibull2_lambdas[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);
        distr_weibull2_alphas[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);

        if(distr_weibull2_lambdas[coi].get_copulapos()==0)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[0]= (&distr_weibull2_lambdas[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response2 = distr_weibull2_lambdas[coi].response;
          }
        else if(distr_weibull2_lambdas[coi].get_copulapos()==1)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]= (&distr_weibull2_lambdas[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response1 = distr_weibull2_lambdas[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_gumbel_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_gumbel_mus.size();coi++)
        {

        predict_mult_distrs.push_back(&distr_gumbel_sigmas[coi]);
        predict_mult_distrs.push_back(&distr_gumbel_mus[coi]);

        distr_gumbel_mus[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);
        distr_gumbel_sigmas[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);

        if(distr_gumbel_mus[coi].get_copulapos()==0)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[0]= (&distr_gumbel_mus[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response2 = distr_gumbel_mus[coi].response;
          }
        else if(distr_gumbel_mus[coi].get_copulapos()==1)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]= (&distr_gumbel_mus[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response1 = distr_gumbel_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
      if(distr_gumbel2_mus.size()>0)
        {
      int coi;
      for(coi=0;coi<distr_gumbel2_mus.size();coi++)
        {

        predict_mult_distrs.push_back(&distr_gumbel2_sigma2s[coi]);
        predict_mult_distrs.push_back(&distr_gumbel2_mus[coi]);

        distr_gumbel2_mus[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);
        distr_gumbel2_sigma2s[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);

        if(distr_gumbel2_mus[coi].get_copulapos()==0)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[0]= (&distr_gumbel2_mus[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response2 = distr_gumbel2_mus[coi].response;
          }
        else if(distr_gumbel2_mus[coi].get_copulapos()==1)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]= (&distr_gumbel2_mus[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response1 = distr_gumbel2_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_gamma_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_gamma_mus.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_gamma_mus[coi]);
        predict_mult_distrs.push_back(&distr_gamma_sigmas[coi]);

        distr_gamma_mus[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);
        distr_gamma_sigmas[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);

        if(distr_gamma_mus[coi].get_copulapos()==0)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[0]= (&distr_gamma_mus[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response2 = distr_gamma_mus[coi].response;
          }
        else if(distr_gamma_mus[coi].get_copulapos()==1)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]= (&distr_gamma_mus[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response1 = distr_gamma_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_normal_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_normal_mus.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_normal_sigma2s[coi]);
        predict_mult_distrs.push_back(&distr_normal_mus[coi]);

        //distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp.push_back(&distr_normal_mus[coi]);

        distr_normal_mus[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);
        distr_normal_sigma2s[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);

        if(distr_normal_mus[coi].get_copulapos()==0)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[0]= (&distr_normal_mus[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response2 = distr_normal_mus[coi].response;
          }
        else if(distr_normal_mus[coi].get_copulapos()==1)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]= (&distr_normal_mus[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response1 = distr_normal_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_dagum_as.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_dagum_as.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_dagum_ps[coi]);
        predict_mult_distrs.push_back(&distr_dagum_bs[coi]);
        predict_mult_distrs.push_back(&distr_dagum_as[coi]);

        //distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp.push_back(&distr_dagum_as[coi]);

        distr_dagum_as[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);
        distr_dagum_bs[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);
        distr_dagum_ps[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);

        if(distr_dagum_as[coi].get_copulapos()==0)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[0]= (&distr_dagum_as[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response2 = distr_dagum_as[coi].response;
          }
        else if(distr_dagum_as[coi].get_copulapos()==1)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]= (&distr_dagum_as[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response1 = distr_dagum_as[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_binomialprobit_copulas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_binomialprobit_copulas.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_binomialprobit_copulas[coi]);

       // distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp.push_back(&distr_binomialprobit_copulas[coi]);

        distr_binomialprobit_copulas[coi].distrcopulap.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);

        if(distr_binomialprobit_copulas[coi].get_copulapos()==0)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[0]= (&distr_binomialprobit_copulas[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response2 = distr_binomialprobit_copulas[coi].response;
          distr_binomialprobit_copulas[coi].responsecopmat = &distr_clayton_copulas[distr_clayton_copulas.size()-1].response2;
          }
        else if(distr_binomialprobit_copulas[coi].get_copulapos()==1)
          {
          distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]= (&distr_binomialprobit_copulas[coi]);
          distr_clayton_copulas[distr_clayton_copulas.size()-1].response1 = distr_binomialprobit_copulas[coi].response;
          distr_binomialprobit_copulas[coi].responsecopmat = &distr_clayton_copulas[distr_clayton_copulas.size()-1].response1;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    predict_mult_distrs.push_back(&distr_clayton_copulas[distr_clayton_copulas.size()-1]);

    // initialize copulaoffset for the second marginal (required for marginal deviance etc.)
    int off = distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[0]->distrp.size()+1;
    distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]->copulaoffset = off;
    for(unsigned k=0; k<distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]->distrp.size(); k++)
      distr_clayton_copulas[distr_clayton_copulas.size()-1].distrp[1]->distrp[k]->copulaoffset = off;

    if(distr_binomialprobit_copulas.size()>0)
      {
      if(setseed.getvalue() >= 0)
        srand(setseed.getvalue());

      distr_clayton_copulas[0].update_end();
      int coi;
      for(coi=0;coi<distr_binomialprobit_copulas.size();coi++)
        {
        distr_binomialprobit_copulas[coi].update();
        }
      }
    }
//-------------------------- END: clayton_copula---------------------------


//----------------------------- gumbel_copula --------------------------------
  else if ((family.getvalue() == "gumbel_copula") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_gumbel_copulas.push_back(DISTR_gumbel_copula(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gumbel_copulas[distr_gumbel_copulas.size()-1];
    equations[modnr].pathd = "";

    if ((countmarginal != 2))
      {
      outerror("ERROR: Two equations for marginal distributions required or family not implemented yet");
      return true;
      }

    // Fill distrp of the copula with some starting values to make it a vector of length=2
    distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp.push_back(equations[modnr-1].distrp);
    distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp.push_back(equations[modnr-2].distrp);

    if(distr_weibull_lambdas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_weibull_lambdas.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_weibull_alphas[coi]);
        predict_mult_distrs.push_back(&distr_weibull_lambdas[coi]);

        distr_weibull_lambdas[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);
        distr_weibull_alphas[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);

        if(distr_weibull_lambdas[coi].get_copulapos()==0)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[0]= (&distr_weibull_lambdas[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response2 = distr_weibull_lambdas[coi].response;
          }
        else if(distr_weibull_lambdas[coi].get_copulapos()==1)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[1]= (&distr_weibull_lambdas[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response1 = distr_weibull_lambdas[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_weibull2_lambdas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_weibull2_lambdas.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_weibull2_alphas[coi]);
        predict_mult_distrs.push_back(&distr_weibull2_lambdas[coi]);

        distr_weibull2_lambdas[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);
        distr_weibull2_alphas[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);

        if(distr_weibull2_lambdas[coi].get_copulapos()==0)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[0]= (&distr_weibull2_lambdas[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response2 = distr_weibull2_lambdas[coi].response;
          }
        else if(distr_weibull2_lambdas[coi].get_copulapos()==1)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[1]= (&distr_weibull2_lambdas[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response1 = distr_weibull2_lambdas[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_gumbel_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_gumbel_mus.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_gumbel_mus[coi]);
        predict_mult_distrs.push_back(&distr_gumbel_sigmas[coi]);

        distr_gumbel_mus[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);
        distr_gumbel_sigmas[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);

        if(distr_gumbel_mus[coi].get_copulapos()==0)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[0]= (&distr_gumbel_mus[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response2 = distr_gumbel_mus[coi].response;
          }
        else if(distr_gumbel_mus[coi].get_copulapos()==1)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[1]= (&distr_gumbel_mus[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response1 = distr_gumbel_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_gamma_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_gamma_mus.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_gamma_mus[coi]);
        predict_mult_distrs.push_back(&distr_gamma_sigmas[coi]);

        distr_gamma_mus[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);
        distr_gamma_sigmas[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);

        if(distr_gamma_mus[coi].get_copulapos()==0)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[0]= (&distr_gamma_mus[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response2 = distr_gamma_mus[coi].response;
          }
        else if(distr_gamma_mus[coi].get_copulapos()==1)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[1]= (&distr_gamma_mus[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response1 = distr_gamma_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_normal_mus.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_normal_mus.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_normal_sigma2s[coi]);
        predict_mult_distrs.push_back(&distr_normal_mus[coi]);

        distr_normal_mus[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);
        distr_normal_sigma2s[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);

        if(distr_normal_mus[coi].get_copulapos()==0)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[0]= (&distr_normal_mus[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response2 = distr_normal_mus[coi].response;
          }
        else if(distr_normal_mus[coi].get_copulapos()==1)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[1]= (&distr_normal_mus[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response1 = distr_normal_mus[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_dagum_as.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_dagum_as.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_dagum_ps[coi]);
        predict_mult_distrs.push_back(&distr_dagum_bs[coi]);
        predict_mult_distrs.push_back(&distr_dagum_as[coi]);

        distr_dagum_as[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);
        distr_dagum_bs[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);
        distr_dagum_ps[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);

        if(distr_dagum_as[coi].get_copulapos()==0)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[0]= (&distr_dagum_as[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response2 = distr_dagum_as[coi].response;
          }
        else if(distr_dagum_as[coi].get_copulapos()==1)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[1]= (&distr_dagum_as[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response1 = distr_dagum_as[coi].response;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    if(distr_binomialprobit_copulas.size()>0)
      {
      int coi;
      for(coi=0;coi<distr_binomialprobit_copulas.size();coi++)
        {
        predict_mult_distrs.push_back(&distr_binomialprobit_copulas[coi]);

        distr_binomialprobit_copulas[coi].distrcopulap.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);

        if(distr_binomialprobit_copulas[coi].get_copulapos()==0)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[0]= (&distr_binomialprobit_copulas[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response2 = distr_binomialprobit_copulas[coi].response;
          distr_binomialprobit_copulas[coi].responsecopmat = &distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response2;
          }
        else if(distr_binomialprobit_copulas[coi].get_copulapos()==1)
          {
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[1]= (&distr_binomialprobit_copulas[coi]);
          distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response1 = distr_binomialprobit_copulas[coi].response;
          distr_binomialprobit_copulas[coi].responsecopmat = &distr_gumbel_copulas[distr_gumbel_copulas.size()-1].response1;
          }
        else
          {
          outerror("ERROR: Two equations for marginal distributions required");
          return true;
          }
        }
      }
    predict_mult_distrs.push_back(&distr_gumbel_copulas[distr_gumbel_copulas.size()-1]);

    // initialize copulaoffset for the second marginal (required for marginal deviance etc.)
    int off = distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[0]->distrp.size()+1;
    distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[1]->copulaoffset = off;
    for(unsigned k=0; k<distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[1]->distrp.size(); k++)
      distr_gumbel_copulas[distr_gumbel_copulas.size()-1].distrp[1]->distrp[k]->copulaoffset = off;

    if(distr_binomialprobit_copulas.size()>0)
      {
      if(setseed.getvalue() >= 0)
        srand(setseed.getvalue());

      distr_gumbel_copulas[0].update_end();

      int coi;
      for(coi=0;coi<distr_binomialprobit_copulas.size();coi++)
        {
        distr_binomialprobit_copulas[coi].update();
        }
      }
    }
//-------------------------- END: clayton_copula---------------------------



// ----------------------------------- tcopula_df ----------------------
   else if ((family.getvalue() == "tcopula") && (equationtype.getvalue()=="df"))
     {

    computemodeforstartingvalues = true;

     unsigned pos;
     if (equationtype.getvalue()=="u1")
       pos = 1;
     else
       pos = 0;

     distr_tcopula_dfs.push_back(DISTR_tcopula_df(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_tcopula_dfs[distr_tcopula_dfs.size()-1];
     equations[modnr].pathd = "";


     }
 //------------------------------ END: tcopula_df ---------------------------

//----------------------------- tcopula_rho ----------------------
  else if ((family.getvalue() == "tcopula") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;
    computemodeforstartingvalues = true;

    distr_tcopula_rhos.push_back(DISTR_tcopula_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_tcopula_rhos[distr_tcopula_rhos.size()-1];
    equations[modnr].pathd = "";

    if ((distr_tcopula_dfs.size() != 1))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }

            predict_mult_distrs.push_back(&distr_tcopula_dfs[distr_tcopula_dfs.size()-1]);
			predict_mult_distrs.push_back(&distr_tcopula_rhos[distr_tcopula_rhos.size()-1]);

            distr_tcopula_dfs[distr_tcopula_dfs.size()-1].response2 = distr_tcopula_rhos[distr_tcopula_rhos.size()-1].response;
            distr_tcopula_rhos[distr_tcopula_rhos.size()-1].response2 = distr_tcopula_dfs[distr_tcopula_dfs.size()-1].response;

			distr_tcopula_dfs[distr_tcopula_dfs.size()-1].distrp.push_back(
            &distr_tcopula_rhos[distr_tcopula_rhos.size()-1]);

            distr_tcopula_rhos[distr_tcopula_rhos.size()-1].distrp.push_back(
            &distr_tcopula_dfs[distr_tcopula_dfs.size()-1]);


    }
//-------------------------- END: tcopula_rho ---------------------


//-------------------------------- frankcopula2_normal sigma2 ---------------------------------
  else if (family.getvalue() == "frankcopula2_normal_sigma2" && equationtype.getvalue()=="sigma2")
    {

	unsigned pos;
      if (distr_frankcopula2_normal_sigma2s.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_frankcopula2_normal_sigma2s.push_back(DISTR_frankcopula2_normal_sigma2(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1];
    equations[modnr].pathd = "";

/*	if (distr_frankcopula2_normal_sigma2s.size() == 2)
       {
       distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].response2 = distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].response;
       distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].response2 = distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].response;
       }*/


    }
//---------------------------- END: frankcopula2_normal sigma2 -------------------------------

//------------------------------- frankcopula2_normal mu ------------------------------------
  else if ((family.getvalue() == "frankcopula2_normal_mu") && equationtype.getvalue()=="mu")
    {

	unsigned pos;
      if (distr_frankcopula2_normal_mus.size()==0)
        pos=0;
      else
        pos=1;

    computemodeforstartingvalues = true;

    distr_frankcopula2_normal_mus.push_back(DISTR_frankcopula2_normal_mu(&generaloptions,D.getCol(0),pos,w));

    equations[modnr].distrp = &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1];
    equations[modnr].pathd = "";

	/*if (distr_frankcopula2_normal_mus.size() == 2)
       {
       distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].response2 = distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].response;
       distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].response2 = distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].response;
       }*/

    }
//------------------------------- END: frankcopula2_normal mu -------------------------------

//----------------------------- frankcopula2_rho ----------------------
  else if ((family.getvalue() == "frankcopula2_rho") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;
    computemodeforstartingvalues = true;

    distr_frankcopula2_rhos.push_back(DISTR_frankcopula2_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1];
    equations[modnr].pathd = "";

    if (((distr_normal_mus.size() != 2) || (distr_normal_sigma2s.size() != 2)) && ((distr_frankcopula2_normal_mus.size() != 2) || (distr_frankcopula2_normal_sigma2s.size() != 2)))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }
		if ((distr_normal_mus.size() == 2))
		{
            predict_mult_distrs.push_back(&distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);
            predict_mult_distrs.push_back(&distr_normal_mus[distr_normal_mus.size()-2]);
            predict_mult_distrs.push_back(&distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);
            predict_mult_distrs.push_back(&distr_normal_mus[distr_normal_mus.size()-1]);
            predict_mult_distrs.push_back(&distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1]);

            distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].response2 = distr_normal_mus[distr_normal_mus.size()-2].response;
       //     distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].response = distr_normal_mus[distr_normal_mus.size()-1].response;

            distr_normal_sigma2s[distr_normal_sigma2s.size()-2].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-2]);

            distr_normal_sigma2s[distr_normal_sigma2s.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-1]);

            distr_normal_mus[distr_normal_mus.size()-2].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);

            distr_normal_mus[distr_normal_mus.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);

            distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);

            distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-2]);

            distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);


            distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-1]);
		}
		else
		{
            predict_mult_distrs.push_back(&distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2]);
            predict_mult_distrs.push_back(&distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2]);
            predict_mult_distrs.push_back(&distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1]);
            predict_mult_distrs.push_back(&distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1]);
            predict_mult_distrs.push_back(&distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1]);

            distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].response2 = distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].response;
        //    distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].response = distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].response;
            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].response2 = distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].response;
            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].response2 = distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].response;
            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].response2 = distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].response;
            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].response2 = distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].response;


            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2]);

            distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2]);

            distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2]);

            distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1]);

            distr_frankcopula2_rhos[distr_frankcopula2_rhos.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1]);
		}

    }
//-------------------------- END: frankcopula2_rho ---------------------


//----------------------------- frankcopula2_exp_rho ----------------------
  else if (family.getvalue() == "frankcopula2_exp_rho" && equationtype.getvalue()=="mean")
    {
	mainequation = true;
    computemodeforstartingvalues = true;

    distr_frankcopula2_exp_rhos.push_back(DISTR_frankcopula2_exp_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1];
    equations[modnr].pathd = "";

    if (((distr_normal_mus.size() != 2) || (distr_normal_sigma2s.size() != 2)) && ((distr_frankcopula2_normal_mus.size() != 2) || (distr_frankcopula2_normal_sigma2s.size() != 2)))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }
		if ((distr_normal_mus.size() == 2))
		{
            predict_mult_distrs.push_back(&distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);
            predict_mult_distrs.push_back(&distr_normal_mus[distr_normal_mus.size()-2]);
            predict_mult_distrs.push_back(&distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);
            predict_mult_distrs.push_back(&distr_normal_mus[distr_normal_mus.size()-1]);
            predict_mult_distrs.push_back(&distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1]);

            distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].response2 = distr_normal_mus[distr_normal_mus.size()-2].response;
      //      distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].response = distr_normal_mus[distr_normal_mus.size()-1].response;

            distr_normal_sigma2s[distr_normal_sigma2s.size()-2].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-2]);

            distr_normal_sigma2s[distr_normal_sigma2s.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-1]);

            distr_normal_mus[distr_normal_mus.size()-2].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);

            distr_normal_mus[distr_normal_mus.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);

            distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-2]);

            distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-2]);

            distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].distrp.push_back(
            &distr_normal_sigma2s[distr_normal_sigma2s.size()-1]);


            distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].distrp.push_back(
            &distr_normal_mus[distr_normal_mus.size()-1]);
		}
		else
		{
            predict_mult_distrs.push_back(&distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2]);
            predict_mult_distrs.push_back(&distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2]);
            predict_mult_distrs.push_back(&distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1]);
            predict_mult_distrs.push_back(&distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1]);
            predict_mult_distrs.push_back(&distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1]);

            distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].response2 = distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].response;
         //   distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].response = distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].response;
            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].response2 = distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].response;
            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].response2 = distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].response;
            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].response2 = distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].response;
            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].response2 = distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].response;


            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1]);

            distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1]);

            distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2]);

            distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-2]);

            distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-2]);

            distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_sigma2s[distr_frankcopula2_normal_sigma2s.size()-1]);

            distr_frankcopula2_exp_rhos[distr_frankcopula2_exp_rhos.size()-1].distrp.push_back(
            &distr_frankcopula2_normal_mus[distr_frankcopula2_normal_mus.size()-1]);
		}

    }
//-------------------------- END: frankcopula2_exp_rho ---------------------



//----------------------------- frankcopula_rho ----------------------
  else if ((family.getvalue() == "frankcopula") && (equationtype.getvalue()=="rho"))
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_frankcopula_rhos.push_back(DISTR_frankcopula_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_frankcopula_rhos[distr_frankcopula_rhos.size()-1];
    equations[modnr].pathd = "";

    if ((distr_copulas.size() != 2))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }

            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-2]);
            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-1]);
			predict_mult_distrs.push_back(&distr_frankcopula_rhos[distr_frankcopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-2].response2 = distr_copulas[distr_copulas.size()-1].response;
            distr_copulas[distr_copulas.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;
            distr_frankcopula_rhos[distr_frankcopula_rhos.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;
        //    distr_frankcopula_rhos[distr_frankcopula_rhos.size()-1].response = distr_copulas[distr_copulas.size()-1].response;

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_frankcopula_rhos[distr_frankcopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_frankcopula_rhos[distr_frankcopula_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

			distr_frankcopula_rhos[distr_frankcopula_rhos.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

            distr_frankcopula_rhos[distr_frankcopula_rhos.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);


    }
//-------------------------- END: frankcopula_rho ---------------------


//----------------------------- frankcopula_exp_rho ----------------------
  else if (family.getvalue() == "frankcopula_exp"  && (equationtype.getvalue()=="rho"))
    {
	mainequation = true;
    computemodeforstartingvalues = true;

    distr_frankcopula_exp_rhos.push_back(DISTR_frankcopula_exp_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_frankcopula_exp_rhos[distr_frankcopula_exp_rhos.size()-1];
    equations[modnr].pathd = "";

    if ((distr_copulas.size() != 2))
       {
       outerror("ERROR: Two equations for marginal distributions required");
       return true;
       }

            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-2]);
            predict_mult_distrs.push_back(&distr_copulas[distr_copulas.size()-1]);
			predict_mult_distrs.push_back(&distr_frankcopula_exp_rhos[distr_frankcopula_exp_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-2].response2 = distr_copulas[distr_copulas.size()-1].response;
            distr_copulas[distr_copulas.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;
            distr_frankcopula_exp_rhos[distr_frankcopula_exp_rhos.size()-1].response2 = distr_copulas[distr_copulas.size()-2].response;
        //    distr_frankcopula_exp_rhos[distr_frankcopula_exp_rhos.size()-1].response = distr_copulas[distr_copulas.size()-1].response;

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_frankcopula_exp_rhos[distr_frankcopula_exp_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-2].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_frankcopula_exp_rhos[distr_frankcopula_exp_rhos.size()-1]);

            distr_copulas[distr_copulas.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

			distr_frankcopula_exp_rhos[distr_frankcopula_exp_rhos.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-2]);

            distr_frankcopula_exp_rhos[distr_frankcopula_exp_rhos.size()-1].distrp.push_back(
            &distr_copulas[distr_copulas.size()-1]);

    }
//-------------------------- END: frankcopula_exp_rho ---------------------


// ----------------------------------- copula ----------------------
   else if ((family.getvalue() == "copulas") &&
            ((equationtype.getvalue()=="u1") || (equationtype.getvalue()=="u2"))
           )
     {

    // computemodeforstartingvalues = true;

     unsigned pos;
     if (equationtype.getvalue()=="u1")
       pos = 1;
     else
       pos = 0;

     distr_copulas.push_back(DISTR_copula(&generaloptions,D.getCol(0),pos,w));

     equations[modnr].distrp = &distr_copulas[distr_copulas.size()-1];
     equations[modnr].pathd = "";


     }
 //------------------------------ END: copula ---------------------------

//---------------------------------- ZIP pi ------------------------------------
  else if (family.getvalue() == "zip" && equationtype.getvalue()=="pi")
    {

    computemodeforstartingvalues = true;

    distr_zippis.push_back(DISTR_zippi(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_zippis[distr_zippis.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_zippis[distr_zippis.size()-1]);

    }
//------------------------------- END: ZIP pi ----------------------------------

//-------------------------------- ZIP lambda ----------------------------------
  else if (family.getvalue() == "zip" && equationtype.getvalue()=="lambda")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_ziplambdas.push_back(DISTR_ziplambda(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_ziplambdas[distr_ziplambdas.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_ziplambdas[distr_ziplambdas.size()-1]);

    if (distr_zippis.size() != 1)
      {
      outerror("ERROR: Equation for pi is missing");
      return true;
      }
    else
      {
      distr_zippis[distr_zippis.size()-1].distrlambda =
      &distr_ziplambdas[distr_ziplambdas.size()-1];

      distr_ziplambdas[distr_ziplambdas.size()-1].distrpi =
      &distr_zippis[distr_zippis.size()-1];

      }

    }
//------------------------------- END: ZIP lambda ------------------------------


 //-------------------------------- hurdle_pi ---------------------------------
  else if (family.getvalue() == "hurdle" && equationtype.getvalue()=="pi")
    {

    computemodeforstartingvalues = true;

    distr_hurdle_pis.push_back(DISTR_hurdle_pi(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_hurdle_pis[distr_hurdle_pis.size()-1];
    equations[modnr].pathd = "";



    }
//---------------------------- END: hurdle_pi -------------------------------


//------------------------------- hurdle_lambda ------------------------------------
  else if (family.getvalue() == "hurdle" && equationtype.getvalue()=="lambda")
    {

    mainequation=true;

    computemodeforstartingvalues = true;
    int i;
    if(w.rows()==1)
        w = datamatrix(D.rows(),1,1);

    double * wp = w.getV();
    double * rp = D.getV();
    for(i=0;i<w.rows();i++,wp++,rp+=D.cols())
      {
      if((*rp)==0)
        *wp = 0;
      }

    distr_hurdle_lambdas.push_back(DISTR_hurdle_lambda(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_hurdle_lambdas[distr_hurdle_lambdas.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_hurdle_pis[distr_hurdle_pis.size()-1]);
    predict_mult_distrs.push_back(&distr_hurdle_lambdas[distr_hurdle_lambdas.size()-1]);

    if (distr_hurdle_pis.size() != 1)
      {
      outerror("ERROR: Equation for pi is missing");
      return true;
      }
    else
      {
   //   distr_hurdle_pis[distr_hurdle_pis.size()-1].distrp.push_back(
    //  &distr_hurdle_lambdas[distr_hurdle_lambdas.size()-1]);

      distr_hurdle_lambdas[distr_hurdle_lambdas.size()-1].distrp.push_back
      (&distr_hurdle_pis[distr_hurdle_pis.size()-1]);
      }

    }
//------------------------------- END: hurdle_lambda -------------------------------

//-------------------------------- hurdle delta --------------------------------
  else if (family.getvalue() == "hurdle" && equationtype.getvalue()=="delta")
    {

    computemodeforstartingvalues = true;

//    double flimit = fraclimit.getvalue();
    int strmax = stoprmax.getvalue();
    bool sl = slow.getvalue();
    int nb = nrbetween.getvalue();
    double sts = stopsum.getvalue();
    int i;
    if(w.rows()==1)
        w = datamatrix(D.cols(),1,1);
    double * wp = w.getV();
    double * rp = D.getV();
    for(i=0;i<w.rows();i++,wp++,rp+=D.cols())
      {
      if((*rp)==0)
        *wp = 0;
      }
    distr_hurdle_deltas.push_back(DISTR_hurdle_delta(&generaloptions,D.getCol(0),
    sts,strmax,nb,sl,w));

    equations[modnr].distrp = &distr_hurdle_deltas[distr_hurdle_deltas.size()-1];
    equations[modnr].pathd = "";

    }
//---------------------------- END: hurdle delta -------------------------------

//------------------------------- hurdle mu ------------------------------------
  else if (family.getvalue() == "hurdle" && equationtype.getvalue()=="mu")
    {

    mainequation=true;

    computemodeforstartingvalues = true;
    int i;
    if(w.rows()==1)
        w = datamatrix(D.rows(),1,1);
    double * wp = w.getV();
    double * rp = D.getV();
    for(i=0;i<w.rows();i++,wp++,rp+=D.cols())
      {
      if((*rp)==0)
        *wp = 0;
      }
    distr_hurdle_mus.push_back(DISTR_hurdle_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_hurdle_mus[distr_hurdle_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_hurdle_pis[distr_hurdle_pis.size()-1]);
    predict_mult_distrs.push_back(&distr_hurdle_deltas[distr_hurdle_deltas.size()-1]);
    predict_mult_distrs.push_back(&distr_hurdle_mus[distr_hurdle_mus.size()-1]);

    if (distr_hurdle_deltas.size() != 1)
      {
      outerror("ERROR: Equation for delta is missing");
      return true;
      }
    if (distr_hurdle_pis.size() != 1)
      {
      outerror("ERROR: Equation for pi is missing");
      return true;
      }

      distr_hurdle_deltas[distr_hurdle_deltas.size()-1].distrp.push_back
      (&distr_hurdle_mus[distr_hurdle_mus.size()-1]);

      distr_hurdle_mus[distr_hurdle_mus.size()-1].distrp.push_back
      (&distr_hurdle_deltas[distr_hurdle_deltas.size()-1]);


    }
//------------------------------- END: hurdle mu -------------------------------

//------------------------------  Zero adjusted --------------------------------

  else if (family.getvalue() == "za")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    DISTR * help_b;
    bool binomialfound = find_binomial(help_b);
    if (binomialfound==false)
      {
      outerror("ERROR: Binomial model for pi is missing\n");
      return true;
      }

    DISTR * help_m;
    bool continuous_single_found = find_continuous_singleparam(help_m);

    vector<DISTR*> help_m_vec;
    bool continuous_mult_found = find_continuous_multparam(help_m_vec);
    if ((continuous_single_found==false) && (continuous_mult_found==false))
      {
      outerror("ERROR: Continuous model is missing\n");
      return true;
      }

    if (continuous_single_found)
      {
      distr_zeroadjusteds.push_back(DISTR_zeroadjusted(&generaloptions,
      help_b,help_m));

      equations[modnr].distrp = &distr_zeroadjusteds[distr_zeroadjusteds.size()-1];
      equations[modnr].pathd = "";

      predict_mult_distrs.push_back(help_b);
      predict_mult_distrs.push_back(help_m);
      predict_mult_distrs.push_back(
      &distr_zeroadjusteds[distr_zeroadjusteds.size()-1]);
      }
    else
      {
      distr_zeroadjusted_mults.push_back(DISTR_zeroadjusted_mult(&generaloptions,
      help_b,help_m_vec));

      equations[modnr].distrp = &distr_zeroadjusted_mults[distr_zeroadjusted_mults.size()-1];
      equations[modnr].pathd = "";

      predict_mult_distrs.erase(predict_mult_distrs.begin(),predict_mult_distrs.end());
      predict_mult_distrs.push_back(help_b);
      unsigned i;
      for (i=0;i<help_m_vec.size();i++)
        predict_mult_distrs.push_back(help_m_vec[i]);
      predict_mult_distrs.push_back(
      &distr_zeroadjusted_mults[distr_zeroadjusted_mults.size()-1]);

      }

    }

//------------------------------  Zero adjusted --------------------------------


 //---------------------------------- BCCG_nu -----------------------------------
   else if (family.getvalue() == "BCCG" && equationtype.getvalue()=="nu")
     {

     computemodeforstartingvalues = true;

     distr_BCCG_nus.push_back(DISTR_BCCG_nu(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_BCCG_nus[distr_BCCG_nus.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------- END: BCCG_nu ---------------------------------

 //---------------------------------- BCCG_sigma --------------------------------
   else if (family.getvalue() == "BCCG" && equationtype.getvalue()=="sigma")
     {

     computemodeforstartingvalues = true;

     distr_BCCG_sigmas.push_back(DISTR_BCCG_sigma(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1];
     equations[modnr].pathd = "";

     }
 //------------------------------ END: BCCG_sigma -------------------------------


 // ----------------------------------- BCCG_mu ----------------------------------
   else if (family.getvalue() == "BCCG" && equationtype.getvalue()=="mu")
     {

      mainequation=true;


     computemodeforstartingvalues = true;

     distr_BCCG_mus.push_back(DISTR_BCCG_mu(&generaloptions,D.getCol(0),w));

     equations[modnr].distrp = &distr_BCCG_mus[distr_BCCG_mus.size()-1];
     equations[modnr].pathd = "";

     if (distr_BCCG_nus.size() != 1)
       {
       outerror("ERROR: Equation for degrees of freedom is missing");
       return true;
       }

     if (distr_BCCG_sigmas.size() != 1)
       {
       outerror("ERROR: Equation for sigma is missing");
       return true;
       }

     predict_mult_distrs.push_back(&distr_BCCG_nus[distr_BCCG_nus.size()-1]);
     predict_mult_distrs.push_back(&distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1]);
     predict_mult_distrs.push_back(&distr_BCCG_mus[distr_BCCG_mus.size()-1]);

     distr_BCCG_nus[distr_BCCG_nus.size()-1].distrp.push_back(
     &distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1]);

	 distr_BCCG_nus[distr_BCCG_nus.size()-1].distrp.push_back(
     &distr_BCCG_mus[distr_BCCG_mus.size()-1]);

     distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1].distrp.push_back(
     &distr_BCCG_nus[distr_BCCG_nus.size()-1]);

     distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1].distrp.push_back(
     &distr_BCCG_mus[distr_BCCG_mus.size()-1]);

     distr_BCCG_mus[distr_BCCG_mus.size()-1].distrp.push_back(
     &distr_BCCG_nus[distr_BCCG_nus.size()-1]);

     distr_BCCG_mus[distr_BCCG_mus.size()-1].distrp.push_back(
     &distr_BCCG_sigmas[distr_BCCG_sigmas.size()-1]);

     }
 //------------------------------ END: BCCG_mu ----------------------------------


 //-------------------------------- sfa0_sigma_v ---------------------------------
  else if (family.getvalue() == "sfa0_sigma_v" && equationtype.getvalue()=="sigma_v")
    {

    computemodeforstartingvalues = true;

    distr_sfa0_sigma_vs.push_back(DISTR_sfa0_sigma_v(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa0_sigma_vs[distr_sfa0_sigma_vs.size()-1];
    equations[modnr].pathd = "";

 //   predict_mult_distrs.push_back(&distr_sfa0_sigma_vs[distr_sfa0_sigma_vs.size()-1]);

    }
//---------------------------- END: sfa0_sigma_v -------------------------------

 //-------------------------------- sfa0_sigma_u ---------------------------------
  else if (family.getvalue() == "sfa0_sigma_u" && equationtype.getvalue()=="sigma_u")
    {

    computemodeforstartingvalues = true;

    distr_sfa0_sigma_us.push_back(DISTR_sfa0_sigma_u(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa0_sigma_us[distr_sfa0_sigma_us.size()-1];
    equations[modnr].pathd = "";

  //  predict_mult_distrs.push_back(&distr_sfa0_sigma_us[distr_sfa0_sigma_us.size()-1]);

    }
//---------------------------- END: sfa0_sigma_u -------------------------------


//------------------------------- sfa0_mu_y ------------------------------------
  else if (family.getvalue() == "sfa0_mu_y" && equationtype.getvalue()=="mu")
    {
    mainequation=true;

    computemodeforstartingvalues = true;

    distr_sfa0_mu_ys.push_back(DISTR_sfa0_mu_y(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa0_mu_ys[distr_sfa0_mu_ys.size()-1];
    equations[modnr].pathd = "";

 //   predict_mult_distrs.push_back(&distr_sfa0_mu_ys[distr_sfa0_mu_ys.size()-1]);


    if (distr_sfa0_sigma_vs.size() != 1)
      {
      outerror("ERROR: Equation for sigma_v is missing");
      return true;
      }
    if (distr_sfa0_sigma_us.size() != 1)
      {
      outerror("ERROR: Equation for sigma_u is missing");
      return true;
      }

    else
      {

     predict_mult_distrs.push_back(&distr_sfa0_sigma_vs[distr_sfa0_sigma_vs.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa0_sigma_us[distr_sfa0_sigma_us.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa0_mu_ys[distr_sfa0_mu_ys.size()-1]);


      distr_sfa0_sigma_vs[distr_sfa0_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa0_sigma_us[distr_sfa0_sigma_us.size()-1]);

      distr_sfa0_sigma_vs[distr_sfa0_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa0_mu_ys[distr_sfa0_mu_ys.size()-1]);

      distr_sfa0_sigma_us[distr_sfa0_sigma_us.size()-1].distrp.push_back
      (&distr_sfa0_sigma_vs[distr_sfa0_sigma_vs.size()-1]);

      distr_sfa0_sigma_us[distr_sfa0_sigma_us.size()-1].distrp.push_back
      (&distr_sfa0_mu_ys[distr_sfa0_mu_ys.size()-1]);

      distr_sfa0_mu_ys[distr_sfa0_mu_ys.size()-1].distrp.push_back
      (&distr_sfa0_sigma_vs[distr_sfa0_sigma_vs.size()-1]);

      distr_sfa0_mu_ys[distr_sfa0_mu_ys.size()-1].distrp.push_back
      (&distr_sfa0_sigma_us[distr_sfa0_sigma_us.size()-1]);

      }

    }
//------------------------------- END: sfa0_mu_y -------------------------------


 //-------------------------------- sfa_alpha ---------------------------------
  else if (family.getvalue() == "sfa_alpha" && equationtype.getvalue()=="alpha")
    {

    computemodeforstartingvalues = true;

    distr_sfa_alphas.push_back(DISTR_sfa_alpha(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa_alphas[distr_sfa_alphas.size()-1];
    equations[modnr].pathd = "";

  //  predict_mult_distrs.push_back(&distr_sfa_alphas[distr_sfa_alphas.size()-1]);

    }
//---------------------------- END: sfa_alpha -------------------------------

 //-------------------------------- sfa_sigma_v ---------------------------------
  else if (family.getvalue() == "sfa_sigma_v" && equationtype.getvalue()=="scale")
    {

    computemodeforstartingvalues = true;

    distr_sfa_sigma_vs.push_back(DISTR_sfa_sigma_v(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1];
    equations[modnr].pathd = "";

 //   predict_mult_distrs.push_back(&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);

    }
//---------------------------- END: sfa_sigma_v -------------------------------

 //-------------------------------- sfa_sigma_u ---------------------------------
  else if (family.getvalue() == "sfa_sigma_u" && equationtype.getvalue()=="scale")
    {

    computemodeforstartingvalues = true;

    distr_sfa_sigma_us.push_back(DISTR_sfa_sigma_u(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1];
    equations[modnr].pathd = "";

  //  predict_mult_distrs.push_back(&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);

    }
//---------------------------- END: sfa_sigma_u -------------------------------

 //-------------------------------- sfa_mu_u ---------------------------------
  else if (family.getvalue() == "sfa_mu_u" && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;

    distr_sfa_mu_us.push_back(DISTR_sfa_mu_u(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa_mu_us[distr_sfa_mu_us.size()-1];
    equations[modnr].pathd = "";

  //  predict_mult_distrs.push_back(&distr_sfa_mu_us[distr_sfa_mu_us.size()-1]);

    }
//---------------------------- END: sfa_mu_u -------------------------------


//------------------------------- sfa_mu_y ------------------------------------
  else if (family.getvalue() == "sfa_mu_y" && equationtype.getvalue()=="mu")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_sfa_mu_ys.push_back(DISTR_sfa_mu_y(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1];
    equations[modnr].pathd = "";

 //   predict_mult_distrs.push_back(&distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1]);

    if (distr_sfa_alphas.size() != 1)
      {
      outerror("ERROR: Equation for alpha is missing");
      return true;
      }
    if (distr_sfa_sigma_vs.size() != 1)
      {
      outerror("ERROR: Equation for sigma_v is missing");
      return true;
      }
    if (distr_sfa_sigma_us.size() != 1)
      {
      outerror("ERROR: Equation for sigma_u is missing");
      return true;
      }
    if (distr_sfa_mu_us.size() != 1)
      {
      outerror("ERROR: Equation for mu_u is missing");
      return true;
      }
    else
      {

     predict_mult_distrs.push_back(&distr_sfa_alphas[distr_sfa_alphas.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa_mu_us[distr_sfa_mu_us.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1]);


      distr_sfa_alphas[distr_sfa_alphas.size()-1].distrp.push_back
      (&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);

      distr_sfa_alphas[distr_sfa_alphas.size()-1].distrp.push_back
      (&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);

      distr_sfa_alphas[distr_sfa_alphas.size()-1].distrp.push_back
      (&distr_sfa_mu_us[distr_sfa_mu_us.size()-1]);

      distr_sfa_alphas[distr_sfa_alphas.size()-1].distrp.push_back
      (&distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1]);

      distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa_alphas[distr_sfa_alphas.size()-1]);

      distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);

      distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa_mu_us[distr_sfa_mu_us.size()-1]);

      distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1]);

      distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1].distrp.push_back
      (&distr_sfa_alphas[distr_sfa_alphas.size()-1]);

      distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1].distrp.push_back
      (&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);

      distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1].distrp.push_back
      (&distr_sfa_mu_us[distr_sfa_mu_us.size()-1]);

      distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1].distrp.push_back
      (&distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1]);

      distr_sfa_mu_us[distr_sfa_mu_us.size()-1].distrp.push_back
      (&distr_sfa_alphas[distr_sfa_alphas.size()-1]);

      distr_sfa_mu_us[distr_sfa_mu_us.size()-1].distrp.push_back
      (&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);

      distr_sfa_mu_us[distr_sfa_mu_us.size()-1].distrp.push_back
      (&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);

      distr_sfa_mu_us[distr_sfa_mu_us.size()-1].distrp.push_back
      (&distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1]);

      distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1].distrp.push_back
      (&distr_sfa_alphas[distr_sfa_alphas.size()-1]);

      distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1].distrp.push_back
      (&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);

      distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1].distrp.push_back
      (&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);

      distr_sfa_mu_ys[distr_sfa_mu_ys.size()-1].distrp.push_back
      (&distr_sfa_mu_us[distr_sfa_mu_us.size()-1]);
      }

    }
//------------------------------- END: sfa_mu_y -------------------------------

 //-------------------------------- sfa_mu_u_id ---------------------------------
  else if (family.getvalue() == "sfa_mu_u_id" && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;

    distr_sfa_mu_u_ids.push_back(DISTR_sfa_mu_u_id(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1];
    equations[modnr].pathd = "";

  //  predict_mult_distrs.push_back(&distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1]);

    }
//---------------------------- END: sfa_mu_u_id -------------------------------


//------------------------------- sfa_mu_y_id ------------------------------------
  else if (family.getvalue() == "sfa_mu_y_id" && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;

    distr_sfa_mu_y_ids.push_back(DISTR_sfa_mu_y_id(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1];
    equations[modnr].pathd = "";

 //   predict_mult_distrs.push_back(&distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1]);

    if (distr_sfa_alphas.size() != 1)
      {
      outerror("ERROR: Equation for alpha is missing");
      return true;
      }
    if (distr_sfa_sigma_vs.size() != 1)
      {
      outerror("ERROR: Equation for sigma_v is missing");
      return true;
      }
    if (distr_sfa_sigma_us.size() != 1)
      {
      outerror("ERROR: Equation for sigma_u is missing");
      return true;
      }
    if (distr_sfa_mu_u_ids.size() != 1)
      {
      outerror("ERROR: Equation for mu_u is missing");
      return true;
      }
    else
      {

     predict_mult_distrs.push_back(&distr_sfa_alphas[distr_sfa_alphas.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1]);


      distr_sfa_alphas[distr_sfa_alphas.size()-1].distrp.push_back
      (&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);

      distr_sfa_alphas[distr_sfa_alphas.size()-1].distrp.push_back
      (&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);

      distr_sfa_alphas[distr_sfa_alphas.size()-1].distrp.push_back
      (&distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1]);

      distr_sfa_alphas[distr_sfa_alphas.size()-1].distrp.push_back
      (&distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1]);

      distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa_alphas[distr_sfa_alphas.size()-1]);

      distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);

      distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1]);

      distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1]);

      distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1].distrp.push_back
      (&distr_sfa_alphas[distr_sfa_alphas.size()-1]);

      distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1].distrp.push_back
      (&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);

      distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1].distrp.push_back
      (&distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1]);

      distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1].distrp.push_back
      (&distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1]);

      distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1].distrp.push_back
      (&distr_sfa_alphas[distr_sfa_alphas.size()-1]);

      distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1].distrp.push_back
      (&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);

      distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1].distrp.push_back
      (&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);

      distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1].distrp.push_back
      (&distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1]);

      distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1].distrp.push_back
      (&distr_sfa_alphas[distr_sfa_alphas.size()-1]);

      distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1].distrp.push_back
      (&distr_sfa_sigma_vs[distr_sfa_sigma_vs.size()-1]);

      distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1].distrp.push_back
      (&distr_sfa_sigma_us[distr_sfa_sigma_us.size()-1]);

      distr_sfa_mu_y_ids[distr_sfa_mu_y_ids.size()-1].distrp.push_back
      (&distr_sfa_mu_u_ids[distr_sfa_mu_u_ids.size()-1]);
      }

    }
//------------------------------- END: sfa_mu_y_id -------------------------------


//-------------------------------- sfa2_sigma_v ---------------------------------
  else if (family.getvalue() == "sfa2_sigma_v" && equationtype.getvalue()=="scale")
    {

    computemodeforstartingvalues = true;

    distr_sfa2_sigma_vs.push_back(DISTR_sfa2_sigma_v(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1];
    equations[modnr].pathd = "";

 //   predict_mult_distrs.push_back(&distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1]);

    }
//---------------------------- END: sfa2_sigma_v -------------------------------

 //-------------------------------- sfa2_sigma_u ---------------------------------
  else if (family.getvalue() == "sfa2_sigma_u" && equationtype.getvalue()=="scale")
    {

    computemodeforstartingvalues = true;

    distr_sfa2_sigma_us.push_back(DISTR_sfa2_sigma_u(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1];
    equations[modnr].pathd = "";

  //  predict_mult_distrs.push_back(&distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1]);

    }
//---------------------------- END: sfa2_sigma_u -------------------------------

 //-------------------------------- sfa2_mu_u ---------------------------------
  else if (family.getvalue() == "sfa2_mu_u" && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;

    distr_sfa2_mu_us.push_back(DISTR_sfa2_mu_u(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa2_mu_us[distr_sfa2_mu_us.size()-1];
    equations[modnr].pathd = "";

  //  predict_mult_distrs.push_back(&distr_sfa2_mu_us[distr_sfa2_mu_us.size()-1]);

    }
//---------------------------- END: sfa2_mu_u -------------------------------


//------------------------------- sfa2_mu_y ------------------------------------
  else if (family.getvalue() == "sfa2_mu_y" && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;

    distr_sfa2_mu_ys.push_back(DISTR_sfa2_mu_y(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa2_mu_ys[distr_sfa2_mu_ys.size()-1];
    equations[modnr].pathd = "";

 //   predict_mult_distrs.push_back(&distr_sfa2_mu_ys[distr_sfa2_mu_ys.size()-1]);

    if (distr_sfa2_sigma_vs.size() != 1)
      {
      outerror("ERROR: Equation for sigma_v is missing");
      return true;
      }
    if (distr_sfa2_sigma_us.size() != 1)
      {
      outerror("ERROR: Equation for sigma_u is missing");
      return true;
      }
    if (distr_sfa2_mu_us.size() != 1)
      {
      outerror("ERROR: Equation for mu_u is missing");
      return true;
      }
    else
      {


     predict_mult_distrs.push_back(&distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa2_mu_us[distr_sfa2_mu_us.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa2_mu_ys[distr_sfa2_mu_ys.size()-1]);


      distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1]);

      distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa2_mu_us[distr_sfa2_mu_us.size()-1]);

      distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa2_mu_ys[distr_sfa2_mu_ys.size()-1]);

      distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1].distrp.push_back
      (&distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1]);

      distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1].distrp.push_back
      (&distr_sfa2_mu_us[distr_sfa2_mu_us.size()-1]);

      distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1].distrp.push_back
      (&distr_sfa2_mu_ys[distr_sfa2_mu_ys.size()-1]);

      distr_sfa2_mu_us[distr_sfa2_mu_us.size()-1].distrp.push_back
      (&distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1]);

      distr_sfa2_mu_us[distr_sfa2_mu_us.size()-1].distrp.push_back
      (&distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1]);

      distr_sfa2_mu_us[distr_sfa2_mu_us.size()-1].distrp.push_back
      (&distr_sfa2_mu_ys[distr_sfa2_mu_ys.size()-1]);

      distr_sfa2_mu_ys[distr_sfa2_mu_ys.size()-1].distrp.push_back
      (&distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1]);

      distr_sfa2_mu_ys[distr_sfa2_mu_ys.size()-1].distrp.push_back
      (&distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1]);

      distr_sfa2_mu_ys[distr_sfa2_mu_ys.size()-1].distrp.push_back
      (&distr_sfa2_mu_us[distr_sfa2_mu_us.size()-1]);
      }

    }
//------------------------------- END: sfa2_mu_y -------------------------------

 //-------------------------------- sfa2_mu_u_id ---------------------------------
  else if (family.getvalue() == "sfa2_mu_u_id" && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;

    distr_sfa2_mu_u_ids.push_back(DISTR_sfa2_mu_u_id(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa2_mu_u_ids[distr_sfa2_mu_u_ids.size()-1];
    equations[modnr].pathd = "";

  //  predict_mult_distrs.push_back(&distr_sfa2_mu_u_ids[distr_sfa2_mu_u_ids.size()-1]);

    }
//---------------------------- END: sfa2_mu_u_id -------------------------------


//------------------------------- sfa2_mu_y_id ------------------------------------
  else if (family.getvalue() == "sfa2_mu_y_id" && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;

    distr_sfa2_mu_y_ids.push_back(DISTR_sfa2_mu_y_id(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_sfa2_mu_y_ids[distr_sfa2_mu_y_ids.size()-1];
    equations[modnr].pathd = "";

 //   predict_mult_distrs.push_back(&distr_sfa2_mu_y_ids[distr_sfa2_mu_y_ids.size()-1]);

    if (distr_sfa2_sigma_vs.size() != 1)
      {
      outerror("ERROR: Equation for sigma_v is missing");
      return true;
      }
    if (distr_sfa2_sigma_us.size() != 1)
      {
      outerror("ERROR: Equation for sigma_u is missing");
      return true;
      }
    if (distr_sfa2_mu_u_ids.size() != 1)
      {
      outerror("ERROR: Equation for mu_u is missing");
      return true;
      }
    else
      {

     predict_mult_distrs.push_back(&distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa2_mu_u_ids[distr_sfa2_mu_u_ids.size()-1]);
     predict_mult_distrs.push_back(&distr_sfa2_mu_y_ids[distr_sfa2_mu_y_ids.size()-1]);

      distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1]);

      distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa2_mu_u_ids[distr_sfa2_mu_u_ids.size()-1]);

      distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1].distrp.push_back
      (&distr_sfa2_mu_y_ids[distr_sfa2_mu_y_ids.size()-1]);

      distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1].distrp.push_back
      (&distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1]);

      distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1].distrp.push_back
      (&distr_sfa2_mu_u_ids[distr_sfa2_mu_u_ids.size()-1]);

      distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1].distrp.push_back
      (&distr_sfa2_mu_y_ids[distr_sfa2_mu_y_ids.size()-1]);

      distr_sfa2_mu_u_ids[distr_sfa2_mu_u_ids.size()-1].distrp.push_back
      (&distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1]);

      distr_sfa2_mu_u_ids[distr_sfa2_mu_u_ids.size()-1].distrp.push_back
      (&distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1]);

      distr_sfa2_mu_u_ids[distr_sfa2_mu_u_ids.size()-1].distrp.push_back
      (&distr_sfa2_mu_y_ids[distr_sfa2_mu_y_ids.size()-1]);

      distr_sfa2_mu_y_ids[distr_sfa2_mu_y_ids.size()-1].distrp.push_back
      (&distr_sfa2_sigma_vs[distr_sfa2_sigma_vs.size()-1]);

      distr_sfa2_mu_y_ids[distr_sfa2_mu_y_ids.size()-1].distrp.push_back
      (&distr_sfa2_sigma_us[distr_sfa2_sigma_us.size()-1]);

      distr_sfa2_mu_y_ids[distr_sfa2_mu_y_ids.size()-1].distrp.push_back
      (&distr_sfa2_mu_u_ids[distr_sfa2_mu_u_ids.size()-1]);
      }

    }
//------------------------------- END: sfa2_mu_y_id -------------------------------

//-------------------------------- gaussiancopula_binary_normal_latent ---------------------------------
  else if ( ((family.getvalue() == "gaussiancopula_binary_normal"))
            && equationtype.getvalue()=="latent")
    {

    computemodeforstartingvalues = true;

    distr_gaussiancopula_binary_normal_latents.push_back(DISTR_gaussiancopula_binary_normal_latent(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1];
    equations[modnr].pathd = "";

    }
//---------------------------- END: gaussiancopula_binary_normal_latent -------------------------------

//------------------------------- gaussiancopula_binary_normal_sigma ------------------------------------
  else if ((family.getvalue() == "gaussiancopula_binary_normal") && equationtype.getvalue()=="sigma")
    {

    computemodeforstartingvalues = true;

    distr_gaussiancopula_binary_normal_sigmas.push_back(DISTR_gaussiancopula_binary_normal_sigma(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_binary_normal_sigmas[distr_gaussiancopula_binary_normal_sigmas.size()-1];
    equations[modnr].pathd = "";


    }
//------------------------------- END: ggaussiancopula_binary_normal_sigma -------------------------------


//------------------------------- gaussiancopula_binary_normal_mu ------------------------------------
  else if ( ((family.getvalue() == "gaussiancopula_binary_normal"))
          && equationtype.getvalue()=="mu")
    {

    computemodeforstartingvalues = true;

    distr_gaussiancopula_binary_normal_mus.push_back(DISTR_gaussiancopula_binary_normal_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_binary_normal_mus[distr_gaussiancopula_binary_normal_mus.size()-1];
    equations[modnr].pathd = "";


    }
//------------------------------- END: ggaussiancopula_binary_normal_mu -------------------------------


//----------------------------- gaussiancopula_binary_normal_rho ----------------------
  else if ( ((family.getvalue() == "gaussiancopula_binary_normal"))
              && equationtype.getvalue()=="rho")
    {
	  mainequation = true;
    computemodeforstartingvalues = true;

    distr_gaussiancopula_binary_normal_rhos.push_back(DISTR_gaussiancopula_binary_normal_rho(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_binary_normal_rhos[distr_gaussiancopula_binary_normal_rhos.size()-1];
    equations[modnr].pathd = "";

    if ( ((distr_gaussiancopula_binary_normal_sigmas.size() != 1) || (distr_gaussiancopula_binary_normal_mus.size() != 1)
            || (distr_gaussiancopula_binary_normal_latents.size() != 1))
        )
       {
       outerror("ERROR: three equations for marginal distribution required");
       return true;
       }
    predict_mult_distrs.push_back(&distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1]);
    predict_mult_distrs.push_back(&distr_gaussiancopula_binary_normal_sigmas[distr_gaussiancopula_binary_normal_sigmas.size()-1]);
    predict_mult_distrs.push_back(&distr_gaussiancopula_binary_normal_mus[distr_gaussiancopula_binary_normal_mus.size()-1]);
    predict_mult_distrs.push_back(&distr_gaussiancopula_binary_normal_rhos[distr_gaussiancopula_binary_normal_rhos.size()-1]);

    distr_gaussiancopula_binary_normal_rhos[distr_gaussiancopula_binary_normal_rhos.size()-1].workingresponse2p = &distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1].response;
    distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1].response2 = distr_gaussiancopula_binary_normal_mus[distr_gaussiancopula_binary_normal_mus.size()-1].response;
    distr_gaussiancopula_binary_normal_sigmas[distr_gaussiancopula_binary_normal_sigmas.size()-1].workingresponse2p = &distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1].response;
    distr_gaussiancopula_binary_normal_mus[distr_gaussiancopula_binary_normal_mus.size()-1].workingresponse2p = &distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1].response;

	  distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_sigmas[distr_gaussiancopula_binary_normal_sigmas.size()-1]);

	  distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_mus[distr_gaussiancopula_binary_normal_mus.size()-1]);

	  distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_rhos[distr_gaussiancopula_binary_normal_rhos.size()-1]);

	  distr_gaussiancopula_binary_normal_sigmas[distr_gaussiancopula_binary_normal_sigmas.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_mus[distr_gaussiancopula_binary_normal_mus.size()-1]);

	  distr_gaussiancopula_binary_normal_sigmas[distr_gaussiancopula_binary_normal_sigmas.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1]);

	  distr_gaussiancopula_binary_normal_sigmas[distr_gaussiancopula_binary_normal_sigmas.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_rhos[distr_gaussiancopula_binary_normal_rhos.size()-1]);

	  distr_gaussiancopula_binary_normal_mus[distr_gaussiancopula_binary_normal_mus.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_sigmas[distr_gaussiancopula_binary_normal_sigmas.size()-1]);

	  distr_gaussiancopula_binary_normal_mus[distr_gaussiancopula_binary_normal_mus.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1]);

	  distr_gaussiancopula_binary_normal_mus[distr_gaussiancopula_binary_normal_mus.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_rhos[distr_gaussiancopula_binary_normal_rhos.size()-1]);

    distr_gaussiancopula_binary_normal_rhos[distr_gaussiancopula_binary_normal_rhos.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_latents[distr_gaussiancopula_binary_normal_latents.size()-1]);

    distr_gaussiancopula_binary_normal_rhos[distr_gaussiancopula_binary_normal_rhos.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_sigmas[distr_gaussiancopula_binary_normal_sigmas.size()-1]);

	  distr_gaussiancopula_binary_normal_rhos[distr_gaussiancopula_binary_normal_rhos.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_normal_mus[distr_gaussiancopula_binary_normal_mus.size()-1]);
    }
//-------------------------- END: gaussiancopula_binary_normal_rho ---------------------


//-------------------------------- gaussiancopula_binary_normal_latent ---------------------------------
  else if ( ((family.getvalue() == "gaussiancopula_binary_dagum"))
            && equationtype.getvalue()=="latent")
    {

    computemodeforstartingvalues = true;

    distr_gaussiancopula_binary_dagum_latents.push_back(DISTR_gaussiancopula_binary_dagum_latent(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1];
    equations[modnr].pathd = "";

    }
//---------------------------- END: gaussiancopula_binary_normal_latent -------------------------------

//------------------------------- gaussiancopula_binary_dagum_p ------------------------------------
  else if ((family.getvalue() == "gaussiancopula_binary_dagum") && equationtype.getvalue()=="p")
    {

    computemodeforstartingvalues = true;

    distr_gaussiancopula_binary_dagum_ps.push_back(DISTR_gaussiancopula_binary_dagum_p(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1];
    equations[modnr].pathd = "";


    }
//------------------------------- END: ggaussiancopula_binary_dagum_p -------------------------------


//------------------------------- gaussiancopula_binary_dagum_b ------------------------------------
  else if ((family.getvalue() == "gaussiancopula_binary_dagum") && equationtype.getvalue()=="b")
    {

    computemodeforstartingvalues = true;

    distr_gaussiancopula_binary_dagum_bs.push_back(DISTR_gaussiancopula_binary_dagum_b(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1];
    equations[modnr].pathd = "";


    }
//------------------------------- END: ggaussiancopula_binary_dagum_b -------------------------------

//------------------------------- gaussiancopula_binary_dagum_a ------------------------------------
  else if ((family.getvalue() == "gaussiancopula_binary_dagum") && equationtype.getvalue()=="a")
    {

    computemodeforstartingvalues = true;

    distr_gaussiancopula_binary_dagum_as.push_back(DISTR_gaussiancopula_binary_dagum_a(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1];
    equations[modnr].pathd = "";


    }
//------------------------------- END: ggaussiancopula_binary_dagum_a -------------------------------

//----------------------------- gaussiancopula_binary_dagum_rho ----------------------
  else if ( ((family.getvalue() == "gaussiancopula_binary_dagum"))
              && equationtype.getvalue()=="rho")
    {
    mainequation = true;
    computemodeforstartingvalues = true;

    distr_gaussiancopula_binary_dagum_rhos.push_back(DISTR_gaussiancopula_binary_dagum_rho(&generaloptions,D.getCol(0),w));


    datamatrix sampleselweight;
    int i;
    if(generaloptions.samplesel)
      {
      sampleselweight = distr_gaussiancopula_binary_dagum_latents[0].response;
      double * wpa = (distr_gaussiancopula_binary_dagum_as[0].weight).getV();
      double * wpb = (distr_gaussiancopula_binary_dagum_bs[0].weight).getV();
      double * wpp = (distr_gaussiancopula_binary_dagum_ps[0].weight).getV();
      double * wprho = (distr_gaussiancopula_binary_dagum_rhos[0].weight).getV();
      double * sampleselwp = sampleselweight.getV();
      for(i=0; i<sampleselweight.rows(); i++,wpa++,wpb++,wpp++,wprho++,sampleselwp++)
        {
        if ((*sampleselwp)==0)
          {
          *wpa = 0;
          *wpb = 0;
          *wpp = 0;
          *wprho = 0;
          }
        }
      }

    equations[modnr].distrp = &distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1];
    equations[modnr].pathd = "";

    if ( ((distr_gaussiancopula_binary_dagum_bs.size() != 1) || (distr_gaussiancopula_binary_dagum_as.size() != 1)
           || (distr_gaussiancopula_binary_dagum_ps.size() != 1) || (distr_gaussiancopula_binary_dagum_latents.size() != 1))
        )
       {
       outerror("ERROR: three equations for marginal distribution required");
       return true;
       }
    predict_mult_distrs.push_back(&distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1]);
    predict_mult_distrs.push_back(&distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1]);
    predict_mult_distrs.push_back(&distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1]);
    predict_mult_distrs.push_back(&distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1]);
    predict_mult_distrs.push_back(&distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1]);

    distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1]);

    distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1]);

	distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1]);

	distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1]);

    distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1]);

	distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1]);

	distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1]);

	distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1]);

    distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1]);

	distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1]);

	distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1]);

	distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1]);

    distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1]);

	distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1]);

	distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1]);

	distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1]);

    distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1]);

    distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1]);

    distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1]);

	distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1].distrp.push_back(
    &distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1]);

    distr_gaussiancopula_binary_dagum_rhos[distr_gaussiancopula_binary_dagum_rhos.size()-1].workingresponse2p = &distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1].response;
    distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1].response2 = distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1].response;
    distr_gaussiancopula_binary_dagum_ps[distr_gaussiancopula_binary_dagum_ps.size()-1].workingresponse2p = &distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1].response;
    distr_gaussiancopula_binary_dagum_bs[distr_gaussiancopula_binary_dagum_bs.size()-1].workingresponse2p = &distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1].response;
    distr_gaussiancopula_binary_dagum_as[distr_gaussiancopula_binary_dagum_as.size()-1].workingresponse2p = &distr_gaussiancopula_binary_dagum_latents[distr_gaussiancopula_binary_dagum_latents.size()-1].response;

    if(setseed.getvalue() >= 0)
        srand(setseed.getvalue());
    distr_gaussiancopula_binary_dagum_rhos[0].update_end();
    distr_gaussiancopula_binary_dagum_latents[0].update();

    }
//-------------------------- END: gaussiancopula_binary_dagum_rho ---------------------


//-------------------------------- ZIP pi cloglog ------------------------------
  else if (family.getvalue() == "zip_cloglog" && equationtype.getvalue()=="pi")
    {

    computemodeforstartingvalues = true;

    distr_zip_cloglog_pis.push_back(DISTR_zip_cloglog_pi(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_zip_cloglog_pis[distr_zip_cloglog_pis.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_zip_cloglog_pis[distr_zip_cloglog_pis.size()-1]);

    }
//------------------------------ END: ZIP pi cloglog----------------------------


//------------------------------ ZIP lambda cloglog ----------------------------
  else if (family.getvalue() == "zip_cloglog" && equationtype.getvalue()=="lambda")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    distr_zip_cloglog_mus.push_back(DISTR_zip_cloglog_mu(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_zip_cloglog_mus[distr_zip_cloglog_mus.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_zip_cloglog_mus[distr_zip_cloglog_mus.size()-1]);

    if (distr_zip_cloglog_pis.size() != 1)
      {
      outerror("ERROR: Equation for pi is missing");
      return true;
      }
    else
      {
      distr_zip_cloglog_pis[distr_zip_cloglog_pis.size()-1].distrp.push_back
      (&distr_zip_cloglog_mus[distr_zip_cloglog_mus.size()-1]);

      distr_zip_cloglog_mus[distr_zip_cloglog_mus.size()-1].distrp.push_back
      (&distr_zip_cloglog_pis[distr_zip_cloglog_pis.size()-1]);
      }

    }
//---------------------------- END: ZIP lambda cloglog -------------------------


//----------------------- multinomial probit response --------------------------
  else if (family.getvalue() == "multinom_probit" && equationtype.getvalue()=="servant")
    {

    computemodeforstartingvalues = true;
    unsigned cat = distr_multinomprobits.size();

    distr_multinomprobits.push_back(DISTR_multinomprobit(&generaloptions,cat,false,D.getCol(0)));

    equations[modnr].distrp = &distr_multinomprobits[distr_multinomprobits.size()-1];
    equations[modnr].pathd = "";

    }
  else if (family.getvalue() == "multinom_probit" && equationtype.getvalue()=="main")
    {

    mainequation=true;

    computemodeforstartingvalues = true;
    unsigned cat = distr_multinomprobits.size();

    distr_multinomprobits.push_back(DISTR_multinomprobit(&generaloptions,cat,true,D.getCol(0),w) );

    equations[modnr].distrp = &distr_multinomprobits[distr_multinomprobits.size()-1];
    equations[modnr].pathd = "";

    if (distr_multinomprobits.size() > 1)
      {
      unsigned i;
      for (i=0;i<distr_multinomprobits.size()-1;i++)
        {
        distr_multinomprobits[distr_multinomprobits.size()-1].assign_othercat(&distr_multinomprobits[i]);
        }

      }

    distr_multinomprobits[distr_multinomprobits.size()-1].create_responsecat();

    }
//--------------------- END: multinomial probit response -----------------------

//----------------------- multivariate gaussian response -----------------------
  else if (family.getvalue() == "multgaussian" && equationtype.getvalue()=="servant")
    {

    computemodeforstartingvalues = true;
    unsigned cat = distr_multgaussians.size();

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_multgaussians.push_back(DISTR_multgaussian(aresp.getvalue(),
    bresp.getvalue(),cat,&generaloptions,false,D.getCol(0),path) );

    equations[modnr].distrp = &distr_multgaussians[distr_multgaussians.size()-1];
    equations[modnr].pathd = "";

    predict_mult_distrs.push_back(&distr_multgaussians[distr_multgaussians.size()-1]);

    }
  else if (family.getvalue() == "multgaussian" && equationtype.getvalue()=="main")
    {

    mainequation=true;

    computemodeforstartingvalues = true;
    unsigned cat = distr_multgaussians.size();

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif


    distr_multgaussians.push_back(DISTR_multgaussian(aresp.getvalue(),
    bresp.getvalue(),cat,&generaloptions,true,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_multgaussians[distr_multgaussians.size()-1];
    equations[modnr].pathd = outfile.getvalue() + "_scale.res";

    predict_mult_distrs.push_back(&distr_multgaussians[distr_multgaussians.size()-1]);

    unsigned i;
    vector<DISTR *> dp;
    for (i=0;i<distr_multgaussians.size();i++)
      {
      dp.push_back(&distr_multgaussians[i]);
      }

    for (i=0;i<distr_multgaussians.size();i++)
      distr_multgaussians[i].assign_distributions(dp);

    }
//--------------------- END: multivariate gaussian response --------------------


//----------------------- multinomial logit response ---------------------------
  else if (family.getvalue() == "multinom_logit" && equationtype.getvalue()=="servant")
    {

    computemodeforstartingvalues = true;
    unsigned cat = distr_multinomlogits.size();

    distr_multinomlogits.push_back(DISTR_multinomlogit(&generaloptions,cat,false,D.getCol(0)));

    equations[modnr].distrp = &distr_multinomlogits[distr_multinomlogits.size()-1];
    equations[modnr].pathd = "";

    }
  else if (family.getvalue() == "multinom_logit" && equationtype.getvalue()=="main")
    {

    mainequation=true;

    computemodeforstartingvalues = true;
    unsigned cat = distr_multinomlogits.size();

    distr_multinomlogits.push_back(DISTR_multinomlogit(&generaloptions,cat,true,D.getCol(0),w) );

    equations[modnr].distrp = &distr_multinomlogits[distr_multinomlogits.size()-1];
    equations[modnr].pathd = "";

    if (distr_multinomlogits.size() > 1)
      {
      unsigned i;
      for (i=0;i<distr_multinomlogits.size()-1;i++)
        {
        distr_multinomlogits[distr_multinomlogits.size()-1].assign_othercat(&distr_multinomlogits[i]);
        }

      distr_multinomlogits[distr_multinomlogits.size()-1].create_responsecat();

      }

    }
//--------------------- END: multinomial logit response ------------------------


//----------------------------- Poisson response -------------------------------
  else if (family.getvalue() == "poisson")
    {
    mainequation=true;

    computemodeforstartingvalues = true;

    equations[modnr].equationtype="lambda";

    distr_poissons.push_back(DISTR_poisson(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_poissons[distr_poissons.size()-1];
    equations[modnr].pathd = "";

    if(shared)
      {
      for(unsigned j=0; j<distr_JMs.size(); j++)
        distr_JMs[j].dpois = &distr_poissons[distr_poissons.size()-1];
      }
    }
//-------------------------- END: poisson response -----------------------------

//------------------------------- JM responses ---------------------------------
  else if (family.getvalue() == "gaussian_shared")
    {
    mainequation=false;
    shared=true;

    computemodeforstartingvalues = true;

    equations[modnr].equationtype="mu_shared";

    distr_JMs.push_back(DISTR_JM(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_JMs[distr_JMs.size()-1];
    equations[modnr].pathd = "";
    }

  else if (family.getvalue() == "quantreg_shared")
    {
    mainequation=false;
    shared=true;

    computemodeforstartingvalues = true;

    equations[modnr].equationtype="mu_shared";

    distr_JMs.push_back(DISTR_JM(&generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_JMs[distr_JMs.size()-1];
    equations[modnr].pathd = "";
    }
//---------------------------- END: JM responses -------------------------------

//-------------- Poisson response with extended response function---------------
  else if (family.getvalue() == "poisson_ext")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    equations[modnr].equationtype="lambda";

    distr_poisson_exts.push_back(DISTR_poisson_ext(
    &generaloptions,D.getCol(0),aexp.getvalue(),
    bexp.getvalue(),adaptexp.getvalue(),w));

    equations[modnr].distrp = &distr_poisson_exts[distr_poisson_exts.size()-1];
    equations[modnr].pathd = "";

    }
//--------------End: Poisson response with extended response function-----------

//---- Poisson response with extended response function, linear for eta > 0-----
  else if (family.getvalue() == "poisson_extlin")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    equations[modnr].equationtype="lambda";

    distr_poisson_extlins.push_back(DISTR_poisson_extlin(
    &generaloptions,D.getCol(0),w));

    equations[modnr].distrp =
    &distr_poisson_extlins[distr_poisson_extlins.size()-1];
    equations[modnr].pathd = "";

    }
//-- End: Poisson response with extended response function, linear for eta > 0--

//-------------------------- log-Gaussian response -----------------------------
  else if (((family.getvalue() == "loggaussian") || (family.getvalue() == "lognormal"))
           && (distr_lognormal_sigma2s.size()==0) )
    {

    mainequation=true;

    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_loggaussians.push_back(DISTR_loggaussian(aresp.getvalue(),
                                       bresp.getvalue(),
                                      &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_loggaussians[distr_loggaussians.size()-1];
    equations[modnr].pathd =  outfile.getvalue() + "_scale.res";

    }
//------------------------- END: log-Gaussian response -------------------------

//----------------------------- quantreg response ------------------------------
  else if (family.getvalue() == "quantreg")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    double quant = quantile.getvalue();
    distr_quantregs.push_back(DISTR_quantreg(aresp.getvalue(),
                                              bresp.getvalue(),
                                     &generaloptions,D.getCol(0),path,quant, w) );

    equations[modnr].distrp = &distr_quantregs[distr_quantregs.size()-1];
    equations[modnr].pathd = outfile.getvalue() + "_scale.res";

    if(shared)
      {
      mainequation=false;
      distr_JMs[distr_JMs.size()-1].dist2 = &distr_quantregs[distr_quantregs.size()-1];
      }
    }
//------------------------- END: quantreg response -----------------------------

//-------------------------- gaussian_mixture response -------------------------
  else if (family.getvalue() == "gaussian_mixture")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_gaussianmixtures.push_back(DISTR_gaussianmixture(aresp.getvalue(),
                                              bresp.getvalue(),
                                     &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussianmixtures[distr_gaussianmixtures.size()-1];
    equations[modnr].pathd = outfile.getvalue()  + "_scale.res";

    }
//------------------------- END: gaussian mixture response ---------------------

//--------------- Gaussian response, exponential response function -------------
  else if (family.getvalue() == "gaussian_exp")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_gaussian_exps.push_back(DISTR_gaussian_exp(
                                  aresp.getvalue(),bresp.getvalue(),
                                  &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussian_exps[distr_gaussian_exps.size()-1];
    equations[modnr].pathd =  outfile.getvalue() +  "_scale.res";

    }
//------------- END: Gaussian response, exponential response function ----------
//------------- Gaussian response, multiplicative random effects allowed -------
  else if (family.getvalue() == "gaussian_mult")
    {

#if defined(__BUILDING_LINUX)
    ST::string path = defaultpath + "/temp/" + name  + "_scale.raw";
#else
    ST::string path = defaultpath + "\\temp\\" + name  + "_scale.raw";
#endif

    distr_gaussian_mults.push_back(DISTR_gaussian_mult(
                                  aresp.getvalue(),bresp.getvalue(),
                                  &generaloptions,D.getCol(0),path,w) );

    equations[modnr].distrp = &distr_gaussian_mults[distr_gaussian_mults.size()-1];

#if defined(__BUILDING_LINUX)
    equations[modnr].pathd = defaultpath + "/temp/" + name  + "_scale.res";
#else
    equations[modnr].pathd = defaultpath + "\\temp\\" + name  + "_scale.res";
#endif

    }
//---------- END: Gaussian response, multiplicative random effects allowed -----
//-------------------------- Gaussian random effect ----------------------------
  else if (family.getvalue() == "gaussian_re")
    {

    distr_gaussian_res.push_back(DISTR_gaussian_re(
                                &generaloptions,D.getCol(0),w) );

    equations[modnr].distrp = &distr_gaussian_res[distr_gaussian_res.size()-1];
    equations[modnr].pathd = "";
    }
//------------------------- END: Gaussian random effect ------------------------
//---------------------------- Binomial response -------------------------------
  else if (family.getvalue() == "binomial_logit" || family.getvalue() == "za_binomial_logit")
    {
    datamatrix dnew = D.getCol(0);
    if(family.getvalue() == "binomial_logit")
      mainequation=true;
    else
      {
      mainequation=false;
      for(unsigned i=0; i<dnew.rows(); i++)
        {
        if(dnew(i,0)>0)
          dnew(i,0) = 1;
        }
      }

    computemodeforstartingvalues = true;
    equations[modnr].equationtype="pi";

    distr_binomials.push_back(DISTR_binomial(&generaloptions,dnew,w));
    equations[modnr].distrp = &distr_binomials[distr_binomials.size()-1];
    equations[modnr].pathd = "";

    }
//-------------------------- END: Binomial response ----------------------------
//---------------------------- Binomial response probit-------------------------
  else if (family.getvalue() == "binomial_probit")
    {

    if(generaloptions.copula)
      {
      mainequation=false;
      countmarginal += 1;
      }
    else
      {
      mainequation=true;
      }

    computemodeforstartingvalues = true;

    equations[modnr].equationtype="pi";

    if(generaloptions.copula)
      {
      distr_binomialprobit_copulas.push_back(DISTR_binomialprobit_copula(
      &generaloptions,D.getCol(0),w));

      equations[modnr].distrp = &distr_binomialprobit_copulas[distr_binomialprobit_copulas.size()-1];
      equations[modnr].pathd = "";

      if (countmarginal == 1)
        {
        distr_binomialprobit_copulas[distr_binomialprobit_copulas.size()-1].set_copulapos(0);
        }
      else if(countmarginal == 2)
        {
        distr_binomialprobit_copulas[distr_binomialprobit_copulas.size()-1].set_copulapos(1);
        }
      else
        {
        outerror("ERROR: currently only bivariate copula models implemented");
        return true;
        }
      }
    else
      {
      #if defined(__BUILDING_LINUX)
      ST::string path = defaultpath + "/temp/" + name  + "_latentutilities.raw";
      #else
      ST::string path = defaultpath + "\\temp\\" + name  + "_latentutilities.raw";
      #endif

      distr_binomialprobits.push_back(DISTR_binomialprobit(
      &generaloptions,D.getCol(0),utilities.getvalue(),path,w));

      equations[modnr].distrp = &distr_binomialprobits[distr_binomialprobits.size()-1];
      equations[modnr].pathd = outfile.getvalue() + "_latentvariables.res";
      }

    }
//-------------------------- END: Binomial response probit ---------------------
//---------------------------- Binomial response SVM ---------------------------
  else if (family.getvalue() == "binomial_svm")
    {

    mainequation=true;

    computemodeforstartingvalues = false;

    equations[modnr].equationtype="pi";

    distr_binomialsvms.push_back(DISTR_binomialsvm(
    &generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_binomialsvms[distr_binomialsvms.size()-1];
    equations[modnr].pathd = "";

    }
//-------------------------- END: Binomial response SVM ------------------------
//----------------------- Binomial response logit fruehwirth--------------------
  else if (family.getvalue() == "binomial_logit_l1")
    {

    mainequation=true;

    computemodeforstartingvalues = true;

    equations[modnr].equationtype="pi";

    distr_logit_fruehwirths.push_back(DISTR_logit_fruehwirth(H.getvalue(),
    &generaloptions,D.getCol(0),w));

    equations[modnr].distrp = &distr_logit_fruehwirths[distr_logit_fruehwirths.size()-1];
    equations[modnr].pathd = "";

    }
//-------------------------- END: Binomial response probit ---------------------
  else
    {
    outerror("ERROR: Invalid distribution specification (check family and/or equationtype option)\n");
    return true;
    }

  datamatrix sampleselweight;
  if(generaloptions.samplesel)
    {
    sampleselweight = distr_binomialprobit_copulas[0].responseorig;
    if(distr_dagum_as.size()>0)
      {
      distr_dagum_as[0].weight = sampleselweight;
      distr_dagum_bs[0].weight = sampleselweight;
      distr_dagum_ps[0].weight = sampleselweight;
      }
    if(distr_normal_mus.size()>0)
      {
      distr_normal_mus[0].weight = sampleselweight;
      distr_normal_sigma2s[0].weight = sampleselweight;
      }
    if(distr_gamma_mus.size()>0)
      {
      distr_gamma_mus[0].weight = sampleselweight;
      distr_gamma_sigmas[0].weight = sampleselweight;
      }
    if(distr_weibull_lambdas.size()>0)
      {
      distr_weibull_lambdas[0].weight = sampleselweight;
      distr_weibull_alphas[0].weight = sampleselweight;
      }
    if(distr_gausscopulas.size()>0)
      {
      distr_gausscopulas[0].weight = sampleselweight;
      }
    if(distr_clayton_copulas.size()>0)
      {
      distr_clayton_copulas[0].weight = sampleselweight;
      }
   if(distr_gumbel_copulas.size()>0)
      {
      distr_gumbel_copulas[0].weight = sampleselweight;
      }
    }



  equations[modnr].distrp->responsename=rname;
  equations[modnr].distrp->weightname=wn;
  equations[modnr].distrp->hlevel = hlevel.getvalue();
  equations[modnr].distrp->equationtype = equationtype.getvalue();
  equations[modnr].distrp->familyshort = family.getvalue();



  if (changelinpredlimits.getvalue() == true)
    {
    double min = linpredminlimit.getvalue();
    if (min == -1000000000)
      min = equations[modnr].distrp->linpredminlimit;
    double max = linpredmaxlimit.getvalue();
    if (max == 1000000000)
      max = equations[modnr].distrp->linpredmaxlimit;

    equations[modnr].distrp->changelimits(min,max);

    }

  if (equations[modnr].hlevel==1)
    master.level1_likep.push_back(equations[modnr].distrp);

/*
  if(mainequation)
    {
    datamatrix checkweight = equations[equations.size()-1].distrp->weight;
    unsigned j,k;
    for(j=0; j<equations.size(); j++)
      {
      if(equations[j].hlevel==1)
        {
        if(!generaloptions.samplesel && family.getvalue()!="betainf" && family.getvalue()!="betainf0" &&
           family.getvalue()!="betainf1" && family.getvalue()!="za_gamma" && family.getvalue()!="za_lognormal" &&
           family.getvalue()!="za_invgaussian" && family.getvalue()!="za_binomial_logit" && family.getvalue()!="za_binomial_cloglog" &&
           family.getvalue()!="za" && family.getvalue()!="hurdle")
          {
          datamatrix compareweight = equations[j].distrp->weight;
          for(k=0; k<checkweight.rows(); k++)
            {
            if(checkweight(k,0)!=compareweight(k,0))
              {
              outerror("ERROR: weights specified in the different equations are not equal!\n\n");
//              return true;
              }
            }
          }
        }
      }
    }
*/

  return false;

  }


void superbayesreg::extract_data(unsigned i, datamatrix & d,datamatrix & iv,
                                 unsigned dim_dm)
  {

  d = datamatrix(D.rows(),dim_dm);

  int j1,j2;
  unsigned p;
  for (p=0;p<dim_dm;p++)
    {
    j1 = terms[i].varnames[terms[i].varnames.size()-1-p].isinlist(modelvarnamesv);
    d.putCol(dim_dm-1-p,D.getCol(j1));
    }

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\d.res");
  // d.prettyPrint(out);
  // TEST

  if (terms[i].varnames.size() > dim_dm)
    {
    j2 = terms[i].varnames[0].isinlist(modelvarnamesv);
    iv = D.getCol(j2);
    }
  }


bool superbayesreg::create_predict(void)
  {
  if (predict.getvalue() != "no")
    {

    unsigned modnr = equations.size()-1;

    ST::string h = equations[modnr].paths;

#if defined(__BUILDING_LINUX)
    ST::string pathnonp = defaultpath + "/temp/" + name + "_" + h +
                            "_predict.raw";
#else
    ST::string pathnonp = defaultpath + "\\temp\\" + name + "_" + h +
                            "_predict.raw";
#endif

#if defined(__BUILDING_LINUX)
    ST::string pathnonp2 = defaultpath + "/temp/" + name + "_" + h +
                            "_deviance.raw";
#else
    ST::string pathnonp2 = defaultpath + "\\temp\\" + name + "_" + h +
                            "_deviance.raw";
#endif


    ST::string pathres = outfile.getvalue() +  "_" + h + "_predict.res";


    if ((predict.getvalue() == "full")   ||
        (predict.getvalue() == "fulls")  ||
        (predict.getvalue() == "light")
       )
      {

      if (equations[modnr].distrp->maindistribution == false)
        {
        outerror("ERROR: predict=full is allowed only for the main regression equation\n");
        outerror("HINT: Use predict=predictor to get predictions of the predictor for this equation\n");
        return true;
        }
      else
        {
        if (equations[modnr].distrp->predict_mult == false)
          {
          FC_predicts.push_back(FC_predict(&generaloptions,
                             equations[modnr].distrp,"",pathnonp,
                             pathnonp2,D,modelvarnamesv,WAICoff.getvalue()));

       //   cout << modnr << endl;
          if (predict.getvalue() == "fulls")
            FC_predicts[FC_predicts.size()-1].nosamples=false;

          if (predict.getvalue() == "light")
            FC_predicts[FC_predicts.size()-1].nosamplessave=true;

          if (mse.getvalue() ==  "yes")
            FC_predicts[FC_predicts.size()-1].MSE = MCMC::quadraticMSE;

          if (mse.getvalue() ==  "quadratic")
            FC_predicts[FC_predicts.size()-1].MSE = MCMC::quadraticMSE;

          if (mse.getvalue() ==  "check")
            {
            FC_predicts[FC_predicts.size()-1].MSE = MCMC::checkMSE;
            FC_predicts[FC_predicts.size()-1].MSEparam = mseparam.getvalue();
            }


          equations[modnr].add_FC(&FC_predicts[FC_predicts.size()-1],pathres);
          }
        else // predict_mult==true
          {

          FC_predicts_mult.push_back(FC_predict_mult(&generaloptions,
                             predict_mult_distrs,"",pathnonp,
                             pathnonp2,D,modelvarnamesv,WAICoff.getvalue()));

          if (predict.getvalue() == "fulls")
            FC_predicts_mult[FC_predicts_mult.size()-1].nosamples=false;

          if (predict.getvalue() == "light")
            FC_predicts_mult[FC_predicts_mult.size()-1].nosamplessave=true;


          equations[modnr].add_FC(&FC_predicts_mult[FC_predicts_mult.size()-1],pathres);


          }

        }
      }
    else if (predict.getvalue() == "predictor")
      {

      FC_predict_predictors.push_back(FC_predict_predictor(&generaloptions,
                           equations[modnr].distrp,"",pathnonp,pathnonp2,D,
                           modelvarnamesv));

      equations[modnr].add_FC(&FC_predict_predictors[FC_predict_predictors.size()-1],pathres);

      }


    }

  return false;
  }


void superbayesreg::create_predictive_check(void)
  {

  if (pred_check.getvalue() == true)
    {

    unsigned modnr = equations.size()-1;

    ST::string h = equations[modnr].paths;

    ST::string pathres = outfile.getvalue() +  "_" + h + "_predictive_check.res";

    FC_predictive_checks.push_back(FC_predictive_check(&generaloptions,
                         equations[modnr].distrp,"","",D,modelvarnamesv));

    equations[modnr].add_FC(&FC_predictive_checks[FC_predictive_checks.size()-1]
                           ,pathres);

    }

  }


void superbayesreg::create_cv(void)
  {

  if (cv.getvalue() == true)
    {

    unsigned modnr = equations.size()-1;

    ST::string h = equations[modnr].paths;


#if defined(__BUILDING_LINUX)
    ST::string pathnonp = defaultpath + "/temp/" + name + "_" + h +
                            "_cv_responses.raw";
#else
    ST::string pathnonp = defaultpath + "\\temp\\" + name + "_" + h +
                            "_cv_responses.raw";
#endif


    ST::string pathres = outfile.getvalue() +  "_" + h + "_cv.res";

    FCcv = FC_cv(&generaloptions,equations[modnr].distrp,"",pathnonp,&FC_hrandoms);

    equations[modnr].add_FC(&FCcv,pathres);

    }

  }

void superbayesreg::create_offset(unsigned i)
  {

  datamatrix d,iv;
  extract_data(i,d,iv,1);

  unsigned modnr = equations.size()-1;

  if (equations[modnr].distrp->linpred_current==1)
    equations[modnr].distrp->linearpred1.plus(d);
  else
    equations[modnr].distrp->linearpred2.plus(d);

  equations[modnr].distrp->offsetname = terms[i].varnames[0];

  }

bool superbayesreg::create_linear(void)
  {

  unsigned modnr = equations.size()-1;

  unsigned i;
  int j;

  vector<ST::string> varnames;
  vector<ST::string> varnamesh =  lineareffects.get_constvariables(terms);



  for(i=0;i<varnamesh.size();i++)
    varnames.push_back(varnamesh[i]);

  unsigned nr = varnames.size();


  ST::string title;
  ST::string pathconst;
  ST::string pathconstres;

  ST::string h = equations[modnr].paths;

  title = h + ": linear effects" ;

#if defined(__BUILDING_LINUX)
  pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                           + "_LinearEffects"  +
                           "_" + h + ".raw";
#else
  pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + "_LinearEffects"  +
                           "_" + h + ".raw";
#endif

  pathconstres = outfile.getvalue() + "_" + h + "_LinearEffects.res";

  if (pathconst.isvalidfile() == 1)
    {
    errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
    return true;
    }

  datamatrix X;
  if (nr > 0)
    X = datamatrix(D.rows(),nr,1);

  for (i=0;i<varnames.size();i++)
    {

    j = varnames[i].isinlist(modelvarnamesv);

    if (j != -1)
      {
      unsigned l;
      double * workX=X.getV()+i;
      double * workD=D.getV()+j;
      for (l=0;l<X.rows();l++,workX+=X.cols(),workD+=D.cols())
        *workX = *workD;
      }

    }

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\X.res");
  // X.prettyPrint(out);
  // TEST


  FC_linears.push_back(FC_linear(&master,nrlevel1,&generaloptions,equations[modnr].distrp,X,
                         varnames,title,pathconst,
                         centerlinear.getvalue(),IWLSlineff.getvalue()));

  equations[modnr].add_FC(&FC_linears[FC_linears.size()-1],pathconstres);

  if(shared && mainequation)
    {
    datamatrix help(X.rows(), 1, 0);
    for(j=0; j<distr_JMs.size(); j++)
      {

      ST::string helpname ="alpha" + ST::inttostring(j);
      FC_linears[FC_linears.size()-1].add_variable(help, helpname);
      distr_JMs[j].FClinp = &FC_linears[FC_linears.size()-1];
      distr_JMs[j].FClincol = X.cols()+j;
      }
    }

  return false;

  }

bool superbayesreg::create_merror(unsigned i)
  {
  datamatrix covdata, mevar;
  if (terms[i].options[72] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[72],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[72] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[72] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,mevar,expr);
    }

  if (terms[i].options[73] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[73],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[73] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[73] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,covdata,expr);
    }

  double minx = covdata.min();
  double maxx = covdata.max();
  int f;
  double nint;
  f = (terms[i].options[20]).strtodouble(nint);
  datamatrix datam = datamatrix(covdata.rows(), 1, 0.0);
  double sw = (maxx-minx)/nint;
  double h = minx + 0.5*sw;
  for(unsigned j=0; j<covdata.rows(); j++)
    {
    datam(j,0) = h;
    h += sw;
    if(h>maxx)
      h = minx + 0.5*sw;
    }
//  ofstream out1("c://temp//datam.raw");
//  datam.prettyPrint(out1);
//  out1.close();

  create_pspline(i, datam," with measurement error");

  unsigned modnr = equations.size()-1;
  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_merror.raw","merror"," Measurement Error Correction ");

/*  ofstream out1("c://temp//covdata.raw");
  covdata.prettyPrint(out1);
  out1.close();
  ofstream out2("c://temp//mevar.raw");
  mevar.prettyPrint(out2);
  out2.close();*/

  datamatrix d,iv;
  extract_data(i,d,iv,1);

  FC_merrors.push_back(FC_merror(&generaloptions,title,pathnonp,
                      terms[i].options,terms[i].varnames,covdata,mevar,d,
                      &FC_nonps[FC_nonps.size()-1]));

  equations[modnr].add_FC(&FC_merrors[FC_merrors.size()-1],pathres);

  return false;
  }

void superbayesreg::create_pspline(unsigned i, datamatrix & datam, ST::string t)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_pspline.raw","nonlinear_pspline_effect_of"," Nonlinear effect (P-spline"+t+") of ");

  datamatrix d,iv;
  extract_data(i,d,iv,1);

  if(datam.rows()>1)
    d = datam;

  if (terms[i].options[35] == "ssvs")
    iv = datamatrix(d.rows(),1,1);

  // TEST
  //  ofstream out("c:\\bayesx\\testh\\results\\d.res");
  //  d.prettyPrint(out);
  // TEST

  design_psplines.push_back(DESIGN_pspline(d,iv,&generaloptions,equations[modnr].distrp,
                            &FC_linears[FC_linears.size()-1],
                            terms[i].options,terms[i].varnames));

  FC_nonps.push_back(FC_nonp(&master,nrlevel1,&generaloptions,equations[modnr].distrp,title,
                     pathnonp,&design_psplines[design_psplines.size()-1],
                     terms[i].options,terms[i].varnames));


  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);


  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_pspline_var.raw","variance_of_nonlinear_pspline_effect_of","Variance of nonlinear effect of ");

  bool so = singleomega.getvalue();

  if (terms[i].options[35] == "iid")
    {
    FC_nonp_variances.push_back(FC_nonp_variance(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                                  title,pathnonp,&design_psplines[design_psplines.size()-1],
                                  &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                  terms[i].varnames));

    equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);

    }
  else if (terms[i].options[35] == "ssvs")
    {

    FC_nonp_variance_varselections.push_back(FC_nonp_variance_varselection(
                                  &master,nrlevel1,&generaloptions,equations[modnr].distrp,so,
                                  title,pathnonp,&design_psplines[design_psplines.size()-1],
                                  &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                  terms[i].varnames,ssvsvarlimit.getvalue()));

    equations[modnr].add_FC(&FC_nonp_variance_varselections[FC_nonp_variance_varselections.size()-1],pathres);

    if (so==true)
      {
      if (firstvarselection==true)
        {
        vector<ST::string> help;
        help.push_back("");
        make_paths(pathnonp,pathres,title,help,
        "_omega.raw","omega_of_nonlinear_effects","Prior inclusion probability for nonlinear terms");
        FC_varselection_omegas.push_back(FC_varselection_omega(&master,nrlevel1,&generaloptions,
                                         equations[modnr].distrp,title));
        equations[modnr].add_FC(&FC_varselection_omegas[FC_varselection_omegas.size()-1],pathres);
        firstvarselection=false;
        }

      FC_varselection_omegas[FC_varselection_omegas.size()-1].FC_tau2s.push_back
      (&FC_nonp_variance_varselections[FC_nonp_variance_varselections.size()-1] );

      }

    }

  }


bool superbayesreg::create_userdefined(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_userdefined.raw","nonlinear_userdefined_effect_of"," Nonlinear effect (user defined) of ");

  // adjust to cope with multiple covariates
  datamatrix d,iv;
  extract_data(i,d,iv,1);

  if (terms[i].options[35] == "ssvs")
    iv = datamatrix(d.rows(),1,1);

  datamatrix designmat, penmat, priormean,constrmat,betastart;
  if (terms[i].options[60] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[60],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[60] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[60] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,penmat,expr);
    }

  if (terms[i].options[62] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[62],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[62] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[62] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,designmat,expr);
    }

  if (terms[i].options[63] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[63],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[63] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[63] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,priormean,expr);
    }
  else
    {
    priormean = datamatrix(penmat.cols(),1,0);
    }


if (terms[i].options[71] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[71],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[71] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[71] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,betastart,expr);
    }
  else
    {
    betastart = datamatrix(penmat.cols(),1,0);
    }



  if (terms[i].options[67] != "")
    {
    terms[i].options[16] = "userdefined";
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[67],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[67] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[67] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,constrmat,expr);
    }
//  ST::string pathdesign = terms[i].options[60];
//  ST::string pathpenmat = pathdesign + "_penmat.raw";
//  ST::string pathdesignmat = pathdesign + "_designmat.raw";

//  cout << pathdesignmat << "\n";
//  cout << pathpenmat << "\n";

//  statmatrix<int> test;

/*  ifstream in("C:/temp/testmat.txt");

  if(in.fail())
    cout << "fail" << "\n";
  else
    cout << "success" << "\n";

  ST::string header;
  ST::getline(in,50000,header);
  header = header.eatallcarriagereturns();


  test.prettyScan(in);
  in.close();

  cout << header << "\n";*/

 /* ofstream out1("c:\\temp\\designmatrix_test.raw");
  designmat.prettyPrint(out1);
  out1.close();

  ofstream out2("c:\\temp\\penmatrix_test.raw");
  penmat.prettyPrint(out2);
  out2.close();

  cout << penmat.cols() << "\n";
  cout << penmat.rows() << "\n";

  cout << designmat.cols() << "\n";
  cout << designmat.rows() << "\n";

  cout << "test" << "\n";*/

  int p = designmat.cols();
  bool reduced = true;
  if (terms[i].options[79] == "false")
    reduced = false;

  if(!(priormean.rows()==p &&
       penmat.rows()==p &&
       penmat.cols()==p &&
       betastart.rows()==p &&
       (designmat.rows()==d.rows() || reduced) &&
       betastart.cols()==1 &&
       priormean.cols()==1))
    {
    outerror("ERROR: dimensions in term userdefined do not match:\n\n");
    outerror("       design matrix: (" + ST::inttostring(designmat.rows()) + " x " + ST::inttostring(designmat.cols()) + ")\n");
    outerror("       penalty matrix: (" + ST::inttostring(penmat.rows()) + " x " + ST::inttostring(penmat.cols()) + ")\n");
    outerror("       prior mean: (" + ST::inttostring(priormean.rows()) + " x " + ST::inttostring(priormean.cols()) + ")\n");
    outerror("       starting values: (" + ST::inttostring(betastart.rows()) + " x " + ST::inttostring(betastart.cols()) + ")\n\n");
    return true;
    }

  design_userdefineds.push_back(DESIGN_userdefined(d,iv,
                            designmat, penmat, priormean, constrmat,
                            &generaloptions,equations[modnr].distrp,
                            &FC_linears[FC_linears.size()-1],
                            terms[i].options,terms[i].varnames));

  FC_nonps.push_back(FC_nonp(&master,nrlevel1,&generaloptions,equations[modnr].distrp,title,
                     pathnonp,&design_userdefineds[design_userdefineds.size()-1],
                     terms[i].options,terms[i].varnames,betastart));

  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_userdefined_var.raw","variance_of_nonlinear_userdefined_effect_of","Variance of (user-defined) nonlinear effect of ");


  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_userdefined_var.raw","variance_of_nonlinear_userdefined_effect_of","Variance of nonlinear effect of ");

  bool so = singleomega.getvalue();

  if (terms[i].options[35] == "iid")
    {
    FC_nonp_variances.push_back(FC_nonp_variance(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                                  title,pathnonp,&design_userdefineds[design_userdefineds.size()-1],
                                  &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                  terms[i].varnames));

    equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);
    }
  else if (terms[i].options[35] == "ssvs")
    {

    FC_nonp_variance_varselections.push_back(FC_nonp_variance_varselection(
                                  &master,nrlevel1,&generaloptions,equations[modnr].distrp,so,
                                  title,pathnonp,&design_userdefineds[design_userdefineds.size()-1],
                                  &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                  terms[i].varnames,ssvsvarlimit.getvalue()));

    equations[modnr].add_FC(&FC_nonp_variance_varselections[FC_nonp_variance_varselections.size()-1],pathres);

    if (so==true)
      {
      if (firstvarselection==true)
        {
        vector<ST::string> help;
        help.push_back("");
        make_paths(pathnonp,pathres,title,help,
        "_omega.raw","omega_of_nonlinear_effects","Prior inclusion probability for nonlinear terms");
        FC_varselection_omegas.push_back(FC_varselection_omega(&master,nrlevel1,&generaloptions,
                                         equations[modnr].distrp,title));
        equations[modnr].add_FC(&FC_varselection_omegas[FC_varselection_omegas.size()-1],pathres);
        firstvarselection=false;
        }

      FC_varselection_omegas[FC_varselection_omegas.size()-1].FC_tau2s.push_back
      (&FC_nonp_variance_varselections[FC_nonp_variance_varselections.size()-1] );

      }

    }


  return false;
  }


bool superbayesreg::findREdistr(ST::string & na,equation & maine,unsigned & fnr)
  {
  bool found = false;
  unsigned i;
  for (i=0;i<equations.size()-1;i++)
    {

/*
    TEST
    int hl = equations[i].hlevel;
    ST::string typem = maine.equationtype;
    ST::string type = equations[i].equationtype;
    ST::string name = equations[i].distrp->responsename;
*/

    if ((equations[i].distrp->responsename==na) &&
        (equations[i].hlevel==2) &&
        (equations[i].equationtype==maine.equationtype)
        )
      {
      found = true;
      fnr = i;
      }
    }
  return found;
  }


bool superbayesreg::create_hrandom(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_hrandom.raw","random_effect_of","Random effect of ");

  ST::string pathnonp2 = pathnonp.substr(0,pathnonp.length()-4)  + "_2.res";


  datamatrix d,iv;
  extract_data(i,d,iv,1);

  unsigned fnr;
  ST::string na = terms[i].varnames[terms[i].varnames.size()-1];
  bool found = findREdistr(na,equations[modnr],fnr);

  if (found==false)
    {
    outerror("ERROR: level 2 equation for variable " + na + " not found \n");
    return true;
    }


  design_hrandoms.push_back(DESIGN_hrandom(d,iv,&generaloptions,
                            equations[modnr].distrp,
                            &FC_linears[FC_linears.size()-1],
                             equations[fnr].distrp,
                            terms[i].options,terms[i].varnames));
  if(terms[i].options[35] == "nn")
    {
    FC_hrandom_distributionals.push_back(FC_hrandom_distributional(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                        equations[fnr].distrp, title,pathnonp,pathnonp2,
                        &design_hrandoms[design_hrandoms.size()-1],
                        terms[i].options,terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_distributionals[FC_hrandom_distributionals.size()-1],pathres);
    }
  else
    {
    FC_hrandoms.push_back(FC_hrandom(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                        equations[fnr].distrp, title,pathnonp,pathnonp2,
                        &design_hrandoms[design_hrandoms.size()-1],
                        terms[i].options,terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandoms[FC_hrandoms.size()-1],pathres);
    }
  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_hrandom_var.raw","variance_of_random_effect_of","Variance of random effect of ");

  if (terms[i].options[35] == "iid")
    {
    FC_hrandom_variances.push_back(FC_hrandom_variance(&master,nrlevel1,
                                   &generaloptions,equations[modnr].distrp,
                                   equations[fnr].distrp,
                                  title,pathnonp,
                                  &design_hrandoms[design_hrandoms.size()-1],
                                  &FC_hrandoms[FC_hrandoms.size()-1],
                                  terms[i].options, terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_variances[FC_hrandom_variances.size()-1],pathres);
    }
  else if (terms[i].options[35] == "nn")
    {
    FC_hrandom_variances.push_back(FC_hrandom_variance(&master,nrlevel1,
                                   &generaloptions,equations[modnr].distrp,
                                   equations[fnr].distrp,
                                  title,pathnonp,
                                  &design_hrandoms[design_hrandoms.size()-1],
                                  &FC_hrandom_distributionals[FC_hrandom_distributionals.size()-1],
                                  terms[i].options, terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_variances[FC_hrandom_variances.size()-1],pathres);
    }
  else if (terms[i].options[35] == "lasso")
    {

    FC_hrandom_variance_vecs.push_back(FC_hrandom_variance_vec(&master,nrlevel1,
                                   &generaloptions,equations[modnr].distrp,
                                   equations[fnr].distrp,
                                  title,pathnonp,
                                  &design_hrandoms[design_hrandoms.size()-1],
                                  &FC_hrandoms[FC_hrandoms.size()-1],
                                  terms[i].options, terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_variance_vecs[FC_hrandom_variance_vecs.size()-1],pathres);

    }
  else if (terms[i].options[35] == "nmig")
    {

    FC_hrandom_variance_vec_nmigs.push_back(FC_hrandom_variance_vec_nmig(&master,
                                        nrlevel1,&generaloptions,equations[modnr].distrp,
                                        equations[fnr].distrp,
                                        title,pathnonp,
                                        &design_hrandoms[design_hrandoms.size()-1],
                                        &FC_hrandoms[FC_hrandoms.size()-1],
                                        terms[i].options, terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_variance_vec_nmigs[FC_hrandom_variance_vec_nmigs.size()-1],pathres);

    }
  else if ((terms[i].options[35] == "ssvs"))
    {

    FC_hrandom_variance_ssvss.push_back(FC_hrandom_variance_ssvs(&master,nrlevel1,
                                        &generaloptions,equations[modnr].distrp,
                                        equations[fnr].distrp,
                                        title,pathnonp,
                                        &design_hrandoms[design_hrandoms.size()-1],
                                        &FC_hrandoms[FC_hrandoms.size()-1],
                                        terms[i].options, terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_variance_ssvss[FC_hrandom_variance_ssvss.size()-1],pathres);

    }

  if (imeasures.getvalue() == true && FC_hrandoms.size() >= 1)
      FC_hrandoms[FC_hrandoms.size()-1].imeasures=true;

  return false;

  }


bool superbayesreg::create_random(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_hrandom.raw","random_effect_of","Random effect of ");

  ST::string pathnonp2 = pathnonp.substr(0,pathnonp.length()-4)  + "_2.res";

  datamatrix d,iv;
  extract_data(i,d,iv,1);

  design_hrandoms.push_back(DESIGN_hrandom(d,iv,&generaloptions,
                            equations[modnr].distrp,
                            &FC_linears[FC_linears.size()-1],terms[i].options, terms[i].varnames));

  if(terms[i].options[35] == "nn")
    {
    FC_hrandom_distributionals.push_back(FC_hrandom_distributional(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                         title,pathnonp,pathnonp2,
                        &design_hrandoms[design_hrandoms.size()-1],
                        terms[i].options,terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandom_distributionals[FC_hrandom_distributionals.size()-1],pathres);
    }
  else
    {
    FC_hrandoms.push_back(FC_hrandom(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                         title,pathnonp,pathnonp2,
                        &design_hrandoms[design_hrandoms.size()-1],
                        terms[i].options,terms[i].varnames));

    equations[modnr].add_FC(&FC_hrandoms[FC_hrandoms.size()-1],pathres);
    }

  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_hrandom_var.raw","variance_of_random_effect_of","Variance of random effect of ");

  if(terms[i].options[35] == "nn")
    {
    FC_hrandom_variances.push_back(FC_hrandom_variance(&master,nrlevel1,
                                   &generaloptions,equations[modnr].distrp,
                                  title,pathnonp,
                                  &design_hrandoms[design_hrandoms.size()-1],
                                  &FC_hrandom_distributionals[FC_hrandom_distributionals.size()-1],
                                  terms[i].options, terms[i].varnames));
    }
  else
    {
    FC_hrandom_variances.push_back(FC_hrandom_variance(&master,nrlevel1,
                                   &generaloptions,equations[modnr].distrp,
                                  title,pathnonp,
                                  &design_hrandoms[design_hrandoms.size()-1],
                                  &FC_hrandoms[FC_hrandoms.size()-1],
                                  terms[i].options, terms[i].varnames));
    }
  equations[modnr].add_FC(&FC_hrandom_variances[FC_hrandom_variances.size()-1],pathres);

  if (imeasures.getvalue() == true && FC_hrandoms.size() >= 1)
      FC_hrandoms[FC_hrandoms.size()-1].imeasures=true;

  return false;

  }

bool superbayesreg::create_multiplicative_pspline(unsigned i)
  {
  unsigned modnr = equations.size()-1;

  terms[i].options[12] = "true";

  datamatrix datadummy = datamatrix(1,1,0.0);
  create_pspline(i, datadummy);

  // FC_nonp * fcnp_pspline = &FC_nonps[FC_nonps.size()-1];
  MCMC::DESIGN * dp_pspline = &design_psplines[design_psplines.size()-1];

  dp_pspline->changingdesign=true;
  datamatrix effect(D.rows(),1,1);
  dp_pspline->set_intvar(effect,0);
  dp_pspline->meaneffectintvar = 1;

  term helpt = terms[i];
  terms[i].varnames.erase(terms[i].varnames.begin(),terms[i].varnames.end());
  terms[i].varnames.push_back(helpt.varnames[1]);
  terms[i].varnames.push_back(helpt.varnames[0]);

  terms[i].options[12] = "true";
  // bool error=false;

  unsigned fnr;
  ST::string na = terms[i].varnames[terms[i].varnames.size()-1];
  bool found = findREdistr(na,equations[modnr],fnr);

  if (found==false)
    {
    outerror("ERROR: level 2 equation for variable " + na + " not found \n");
    return true;
    }

  FC_nonps[FC_nonps.size()-1].set_multiplicative(equations[fnr].distrp);

  (distr_gaussian_multeffects[distr_gaussian_multeffects.size()-1]).dg = &(distr_gaussians[distr_gaussians.size()-1]);
  (distr_gaussian_multeffects[distr_gaussian_multeffects.size()-1]).dgexists = true;

  return false;
  }


bool superbayesreg::create_random_pspline(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  terms[i].options[12] = "true";
  bool multexp = false;
  if (terms[i].options[0] == "hrandomexp_pspline")
    {
    terms[i].options[17] = "true";
    multexp=true;
    }

  datamatrix datadummy = datamatrix(1,1,0.0);
  create_pspline(i, datadummy);
  FC_nonp * fcnp_pspline = &FC_nonps[FC_nonps.size()-1];
  MCMC::DESIGN * dp_pspline = &design_psplines[design_psplines.size()-1];
  dp_pspline->changingdesign=true;
  datamatrix effect(D.rows(),1,1);
  dp_pspline->set_intvar(effect,0);
  dp_pspline->meaneffectintvar = 1;

  FC_mults.push_back(FC_mult(true,multexp));
  equations[modnr].add_FC(&FC_mults[FC_mults.size()-1],"");

  term helpt = terms[i];
  terms[i].varnames.erase(terms[i].varnames.begin(),terms[i].varnames.end());
  terms[i].varnames.push_back(helpt.varnames[1]);
  terms[i].varnames.push_back(helpt.varnames[0]);

  terms[i].options[12] = "true";
  bool error=false;
  error = create_hrandom(i);
  if (error==true)
    return true;
  FC_nonp * fcnp_hrandom = &FC_hrandoms[FC_hrandoms.size()-1];
  MCMC::DESIGN * dp_hrandom = &design_hrandoms[design_hrandoms.size()-1];
  dp_hrandom->changingdesign=true;
  dp_hrandom->meaneffectintvar=0;
  dp_hrandom->center=true;

  FC_mults.push_back(FC_mult(false,multexp));

  FC_mults[FC_mults.size()-2].set_effectp(dp_pspline,fcnp_pspline);
  FC_mults[FC_mults.size()-2].set_intp(dp_hrandom,fcnp_hrandom);

  FC_mults[FC_mults.size()-1].set_effectp(dp_hrandom,fcnp_hrandom);
  FC_mults[FC_mults.size()-1].set_intp(dp_pspline,fcnp_pspline);

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_mult.raw","_multiplicative_effect_of",
             "Multiplicative effect of ");

  bool samplem;
  if (terms[i].options[13] == "false")
    samplem = false;
  else
    samplem = true;

  bool me=false ;
  if (terms[i].options[19] == "true")
    me = true;

  double mec;
  int f;
   f = (terms[i].options[34]).strtodouble(mec);


  FC_mults[FC_mults.size()-1].set_multeffects(&master,nrlevel1,&generaloptions,
                                              title,pathnonp,samplem,me,mec);

  equations[modnr].add_FC(&FC_mults[FC_mults.size()-1],pathres);

  return false;

  }



bool  superbayesreg::create_mrf_pspline(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  terms[i].options[12] = "true";

  datamatrix datadummy = datamatrix(1,1,0.0);
  create_pspline(i, datadummy);
  FC_nonp * fcnp_pspline = &FC_nonps[FC_nonps.size()-1];
  MCMC::DESIGN * dp_pspline = &design_psplines[design_psplines.size()-1];
  dp_pspline->changingdesign=true;
  datamatrix effect(D.rows(),1,1);
  dp_pspline->set_intvar(effect,0);
  dp_pspline->meaneffectintvar = 1;

  FC_mults.push_back(FC_mult(true,false));
  equations[modnr].add_FC(&FC_mults[FC_mults.size()-1],"");

  term helpt = terms[i];
  terms[i].varnames.erase(terms[i].varnames.begin(),terms[i].varnames.end());
  terms[i].varnames.push_back(helpt.varnames[1]);
  terms[i].varnames.push_back(helpt.varnames[0]);

  terms[i].options[12] = "true";
  bool error=false;
  error = create_mrf(i);
  if (error==true)
    return true;
  FC_nonp * fcnp_mrf = &FC_nonps[FC_nonps.size()-1];
  MCMC::DESIGN * dp_mrf = &design_mrfs[design_mrfs.size()-1];
  dp_mrf->changingdesign=true;
  dp_mrf->meaneffectintvar=0;
  dp_mrf->center=true;

  FC_mults.push_back(FC_mult(false,false));

  FC_mults[FC_mults.size()-2].set_effectp(dp_pspline,fcnp_pspline);
  FC_mults[FC_mults.size()-2].set_intp(dp_mrf,fcnp_mrf);

  FC_mults[FC_mults.size()-1].set_effectp(dp_mrf,fcnp_mrf);
  FC_mults[FC_mults.size()-1].set_intp(dp_pspline,fcnp_pspline);

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_mult.raw","_multiplicative_effect_of",
             "Multiplicative effect of ");

  bool samplem;
  if (terms[i].options[13] == "false")
    samplem = false;
  else
    samplem = true;

  bool me=false ;
  if (terms[i].options[19] == "true")
    me = true;

  double mec;
  int f;
   f = (terms[i].options[34]).strtodouble(mec);


  FC_mults[FC_mults.size()-1].set_multeffects(&master,nrlevel1,&generaloptions,
                                              title,pathnonp,samplem,me,mec);

  equations[modnr].add_FC(&FC_mults[FC_mults.size()-1],pathres);

  return false;

  }





bool  superbayesreg::create_random_mrf(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  terms[i].options[12] = "true";
  bool multexp = false;

  create_mrf(i);
  FC_nonp * fcnp_mrf = &FC_nonps[FC_nonps.size()-1];
  MCMC::DESIGN * dp_mrf = &design_mrfs[design_mrfs.size()-1];
  dp_mrf->changingdesign=true;
  datamatrix effect(D.rows(),1,1);
  dp_mrf->set_intvar(effect,0);
  dp_mrf->meaneffectintvar = 1;

  FC_mults.push_back(FC_mult(true,multexp));
  equations[modnr].add_FC(&FC_mults[FC_mults.size()-1],"");

  term helpt = terms[i];
  terms[i].varnames.erase(terms[i].varnames.begin(),terms[i].varnames.end());
  terms[i].varnames.push_back(helpt.varnames[1]);
  terms[i].varnames.push_back(helpt.varnames[0]);

  terms[i].options[12] = "true";
  bool error=false;
  error = create_hrandom(i);
  if (error==true)
    return true;
  FC_nonp * fcnp_hrandom = &FC_hrandoms[FC_hrandoms.size()-1];
  MCMC::DESIGN * dp_hrandom = &design_hrandoms[design_hrandoms.size()-1];
  dp_hrandom->changingdesign=true;
  dp_hrandom->meaneffectintvar=0;
  dp_hrandom->center=true;

  FC_mults.push_back(FC_mult(false,multexp));

  FC_mults[FC_mults.size()-2].set_effectp(dp_mrf,fcnp_mrf);
  FC_mults[FC_mults.size()-2].set_intp(dp_hrandom,fcnp_hrandom);

  FC_mults[FC_mults.size()-1].set_effectp(dp_hrandom,fcnp_hrandom);
  FC_mults[FC_mults.size()-1].set_intp(dp_mrf,fcnp_mrf);

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_mult.raw","_multiplicative_effect_of",
             "Multiplicative effect of ");

  bool samplem;
  if (terms[i].options[13] == "false")
    samplem = false;
  else
    samplem = true;

  bool me=false ;
  if (terms[i].options[19] == "true")
    me = true;

  double mec;
  int f;
   f = (terms[i].options[34]).strtodouble(mec);


  FC_mults[FC_mults.size()-1].set_multeffects(&master,nrlevel1,&generaloptions,
                                              title,pathnonp,samplem,me,mec);

  equations[modnr].add_FC(&FC_mults[FC_mults.size()-1],pathres);

  return false;

  }




bool  superbayesreg::find_map(unsigned i,MAP::map & m)
    {

    mapobject * mapp;

    int objpos = findstatobject(*statobj,terms[i].options[8],"map");

    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      mapp = dynamic_cast<mapobject*>(s);
      m = mapp->getmap();
      return false;
      }
    else
      {
      if (objpos == -1)
        {
        if ((terms[i].options[1] == "") || (terms[i].options[1] == " "))
          outerror("ERROR: map object must be specified to estimate a spatial effect\n");
        else
          outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
        }
      else
        outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
      return true;
      }

    }


bool superbayesreg::create_mrf(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_spatial.raw","spatial_MRF_effect_of","Spatial effect (MRF) of ");

  datamatrix d,iv;
  extract_data(i,d,iv,1);

  if (terms[i].options[35] == "ssvs")
    iv = datamatrix(d.rows(),1,1);

  /*
  mapobject * mapp;                           // pointer to mapobject

  int objpos = findstatobject(*statobj,terms[i].options[8],"map");

  if (objpos >= 0)
    {
    statobject * s = statobj->at(objpos);
    mapp = dynamic_cast<mapobject*>(s);
    }
  else
    {
    if (objpos == -1)
      {
      if ((terms[i].options[1] == "") || (terms[i].options[1] == " "))
        outerror("ERROR: map object must be specified to estimate a spatial effect\n");
      else
        outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
      }
    else
      outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
    return true;
    }

  MAP::map m = mapp->getmap();
  */

  MAP::map m;
  if (find_map(i,m) == true)
    return true;

  bool isconnected = m.isconnected();
  if (isconnected==false)
    {
    outerror("ERROR: map is disconnected, spatial effect cannot be estimated\n");
    return true;
    }

  design_mrfs.push_back(DESIGN_mrf(d,iv,&generaloptions,equations[modnr].distrp,
                           &FC_linears[FC_linears.size()-1],m,
                            terms[i].options,terms[i].varnames));

  FC_nonps.push_back(FC_nonp(&master,nrlevel1,&generaloptions,equations[modnr].distrp,title,
                     pathnonp,&design_mrfs[design_mrfs.size()-1],
                     terms[i].options,terms[i].varnames));

  if (FC_nonps[FC_nonps.size()-1].errors==true)
    return true;


  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);

  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_spatial_var.raw","variance_of_spatial_MRF_effect_of","Variance of spatial effect of ");

  bool so = singleomega.getvalue();

  if (terms[i].options[35] == "iid")
    {
    FC_nonp_variances.push_back(FC_nonp_variance(&master,nrlevel1,
                                &generaloptions,equations[modnr].distrp,
                                title,pathnonp,&design_mrfs[design_mrfs.size()-1],
                                &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                terms[i].varnames));

    equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);
    }
  else if (terms[i].options[35] == "ssvs")
    {

    FC_nonp_variance_varselections.push_back(FC_nonp_variance_varselection(
                                  &master,nrlevel1,&generaloptions,equations[modnr].distrp,so,
                                  title,pathnonp,&design_mrfs[design_mrfs.size()-1],
                                  &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                  terms[i].varnames,ssvsvarlimit.getvalue()));

    equations[modnr].add_FC(&FC_nonp_variance_varselections[FC_nonp_variance_varselections.size()-1],pathres);

    if (so==true)
      {
      if (firstvarselection==true)
        {
        vector<ST::string> help;
        help.push_back("");
        make_paths(pathnonp,pathres,title,help,
        "_omega.raw","omega_of_nonlinear_effects","Prior inclusion probability for nonlinear terms");
        FC_varselection_omegas.push_back(FC_varselection_omega(&master,nrlevel1,&generaloptions,
                                         equations[modnr].distrp,title));
        equations[modnr].add_FC(&FC_varselection_omegas[FC_varselection_omegas.size()-1],pathres);
        firstvarselection=false;
        }

      FC_varselection_omegas[FC_varselection_omegas.size()-1].FC_tau2s.push_back
      (&FC_nonp_variance_varselections[FC_nonp_variance_varselections.size()-1] );

      }
    }

  return false;
  }


bool superbayesreg::create_kriging(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_kriging.raw","2dim_kriging_effect_of","2dim effect (kriging) of ");

  datamatrix d,iv;
  extract_data(i,d,iv,2);

  //----------------------------------------------------------------------------

  datamatrix knotdata;
  if (terms[i].options[37] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[37],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[37] + " does not contain any data\n");
        return true;
        }
      else if (datap->getVarnames().size()>2)
        {
        outerror("ERROR: dataset object " + terms[i].options[37] + " contains more than two variables\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[37] + " is not existing\n");
      return true;
      }
    list<ST::string> knotnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(knotnames,knotdata,expr);

    }

  //----------------------------------------------------------------------------

  design_krigings.push_back(DESIGN_kriging(d,iv,&generaloptions,
                            equations[modnr].distrp,
                           &FC_linears[FC_linears.size()-1],
                            terms[i].options,terms[i].varnames,knotdata));

  FC_nonps.push_back(FC_nonp(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                     title, pathnonp,&design_krigings[design_krigings.size()-1],
                     terms[i].options,terms[i].varnames));

  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);

  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_kriging_var.raw","variance_of_2dim_kriging_effect_of",
  "Variance of 2dim effect of ");

  FC_nonp_variances.push_back(FC_nonp_variance(&master,nrlevel1,
                                &generaloptions,equations[modnr].distrp,
                                title,pathnonp,&design_krigings[design_krigings.size()-1],
                                &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                terms[i].varnames));

  equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);

  return false;
  }


bool superbayesreg::create_userdefined_tensor(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_tensor.raw","nonlinear_tensor_effect_of"," Nonlinear tensor effect of ");

  datamatrix d,iv;
  extract_data(i,d,iv,terms[i].varnames.size());

  datamatrix designmat, designmat2, penmat, penmat2, priormean, constrmat, betastart;

  if (terms[i].options[60] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[60],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[60] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[60] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,penmat,expr);
    }


  if (terms[i].options[65] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[65],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[65] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[65] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,penmat2,expr);
    }


  if (terms[i].options[62] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[62],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[62] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[62] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,designmat,expr);
    }


  if (terms[i].options[66] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[66],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[66] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[66] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,designmat2,expr);
    }
  else
    {
    designmat2=datamatrix(1,1,1);
    }


  if (terms[i].options[63] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[63],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[63] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[63] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,priormean,expr);
    }
  else
    {
    priormean = datamatrix(penmat.cols()*penmat2.cols(),1,0);
    }

  if (terms[i].options[67] != "")
    {
    terms[i].options[16] = "userdefined";

    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[67],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[67] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[67] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,constrmat,expr);
    }

  if (terms[i].options[71] != "")
    {
    dataobject * datap;                           // pointer to datasetobject
    int objpos = findstatobject(*statobj,terms[i].options[71],"dataset");
    if (objpos >= 0)
      {
      statobject * s = statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      if (datap->obs()==0 || datap->getVarnames().size()==0)
        {
        outerror("ERROR: dataset object " + terms[i].options[71] + " does not contain any data\n");
        return true;
        }
      }
    else
      {
      outerror("ERROR: dataset object " + terms[i].options[71] + " is not existing\n");
      return true;
      }
    list<ST::string> varnames = datap->getVarnames();
    ST::string expr = "";
    datap->makematrix(varnames,betastart,expr);
    }
  else
    {
    betastart = datamatrix(penmat.cols()*penmat2.cols(),1,0);
    }

//  ST::string pathdesign = terms[i].options[60];
//  ST::string pathpenmat = pathdesign + "_penmat.raw";
//  ST::string pathdesignmat = pathdesign + "_designmat.raw";

//  cout << pathdesignmat << "\n";
//  cout << pathpenmat << "\n";

//  statmatrix<int> test;

/*  ifstream in("C:/temp/testmat.txt");

  if(in.fail())
    cout << "fail" << "\n";
  else
    cout << "success" << "\n";

  ST::string header;
  ST::getline(in,50000,header);
  header = header.eatallcarriagereturns();


  test.prettyScan(in);
  in.close();

  cout << header << "\n";*/

 /* ofstream out1("c:\\temp\\designmatrix_test.raw");
  designmat.prettyPrint(out1);
  out1.close();

  ofstream out2("c:\\temp\\penmatrix_test.raw");
  penmat.prettyPrint(out2);
  out2.close();

  cout << penmat.cols() << "\n";
  cout << penmat.rows() << "\n";

  cout << designmat.cols() << "\n";
  cout << designmat.rows() << "\n";

  cout << "test" << "\n";*/

  if(designmat2.rows()==1 && designmat2.cols()==1)
    {
    int p = designmat.cols();
    int p1 = penmat.rows();
    int p2 = penmat2.rows();
    bool reduced = true;
    if (terms[i].options[79] == "false")
      reduced = false;

    if(!(priormean.rows()==p &&
         priormean.cols()==1 &&
         betastart.rows()==p &&
         betastart.cols()==1 &&
         (designmat.rows()==d.rows() || reduced) &&
         penmat.cols()==p1 &&
         penmat2.cols()==p2 &&
         p1 == p2 &&
         p1*p2 == p))
      {
      outerror("ERROR: dimensions in term userdefined (tensor) do not match:\n\n");
      outerror("       design matrix: (" + ST::inttostring(designmat.rows()) + " x " + ST::inttostring(designmat.cols()) + ")\n");
      outerror("       penalty matrix 1: (" + ST::inttostring(penmat.rows()) + " x " + ST::inttostring(penmat.cols()) + ")\n");
      outerror("       penalty matrix 2: (" + ST::inttostring(penmat2.rows()) + " x " + ST::inttostring(penmat2.cols()) + ")\n");
      outerror("       prior mean: (" + ST::inttostring(priormean.rows()) + " x " + ST::inttostring(priormean.cols()) + ")\n");
      outerror("       starting values: (" + ST::inttostring(betastart.rows()) + " x " + ST::inttostring(betastart.cols()) + ")\n\n");
      return true;
      }
    }
  else
    {
    int p = designmat.cols();
    int p1 = penmat.rows();
    int p2 = penmat2.rows();
    int p3 = designmat2.cols();
    if(!(priormean.rows()==p*p3 &&
         priormean.cols()==1 &&
         betastart.rows()==p*p3 &&
         betastart.cols()==1 &&
         designmat.rows()==d.rows() &&
         designmat2.rows()==d.rows() &&
         penmat.cols()==p1 &&
         penmat2.cols()==p2 &&
         p1 == p &&
         p2 == p3 &&
         p1*p2 == p*p3))
      {
      outerror("ERROR: dimensions in term userdefined (tensor) do not match:\n\n");
      outerror("       design matrix 1: (" + ST::inttostring(designmat.rows()) + " x " + ST::inttostring(designmat.cols()) + ")\n");
      outerror("       design matrix 2: (" + ST::inttostring(designmat2.rows()) + " x " + ST::inttostring(designmat2.cols()) + ")\n");
      outerror("       penalty matrix 1: (" + ST::inttostring(penmat.rows()) + " x " + ST::inttostring(penmat.cols()) + ")\n");
      outerror("       penalty matrix 2: (" + ST::inttostring(penmat2.rows()) + " x " + ST::inttostring(penmat2.cols()) + ")\n");
      outerror("       prior mean: (" + ST::inttostring(priormean.rows()) + " x " + ST::inttostring(priormean.cols()) + ")\n");
      outerror("       starting values: (" + ST::inttostring(betastart.rows()) + " x " + ST::inttostring(betastart.cols()) + ")\n\n");
      return true;
      }

    }
  design_userdefined_tensors.push_back(DESIGN_userdefined_tensor(d,iv,
                            designmat, designmat2, penmat, penmat2, priormean, constrmat,
                            &generaloptions,equations[modnr].distrp,
                            &FC_linears[FC_linears.size()-1],
                            terms[i].options,terms[i].varnames));

  FC_nonps.push_back(FC_nonp(&master,nrlevel1,&generaloptions,equations[modnr].distrp,title,
                     pathnonp,&design_userdefined_tensors[design_userdefined_tensors.size()-1],
                     terms[i].options,terms[i].varnames,betastart));

  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_tensor_var.raw","variance_of_nonlinear_tensor_effect_of","Variance of nonlinear tensor effect of ");


  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_tensor_var.raw","variance_of_nonlinear_tensor_effect_of","Variance of nonlinear tensor effect of ");


  FC_nonp_variances.push_back(FC_nonp_variance(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                              title,pathnonp,&design_userdefined_tensors[design_userdefined_tensors.size()-1],
                              &FC_nonps[FC_nonps.size()-1],terms[i].options,
                              terms[i].varnames));

  equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_tensor_anisotropy.raw","anisotropy_of_nonlinear_tensor_effect_of","Anisotropy of nonlinear tensor effect of ");

  long nro;
  int f;
  f = terms[i].options[68].strtolong(nro);
  unsigned nromegas = (unsigned)nro;
  FC_tensor_omegas.push_back(FC_tensor_omega(&design_userdefined_tensors[design_userdefined_tensors.size()-1],
                                             &FC_nonps[FC_nonps.size()-1],
                                             &generaloptions,
                                             title,nromegas,""));

  equations[modnr].add_FC(&FC_tensor_omegas[FC_tensor_omegas.size()-1],pathres);

  return false;
  }



bool superbayesreg::create_geokriging(unsigned i)
  {

  unsigned modnr = equations.size()-1;

  make_paths(pathnonp,pathres,title,terms[i].varnames,
             "_geokriging.raw","geokriging_effect_of",
             "spatial effect (kriging) of ");

  datamatrix d,iv;
  extract_data(i,d,iv,1);

  /*
  mapobject * mapp;                           // pointer to mapobject

  int objpos = findstatobject(*statobj,terms[i].options[8],"map");

  if (objpos >= 0)
    {
    statobject * s = statobj->at(objpos);
    mapp = dynamic_cast<mapobject*>(s);
    }
  else
    {
    if (objpos == -1)
      {
      if ((terms[i].options[1] == "") || (terms[i].options[1] == " "))
        outerror("ERROR: map object must be specified to estimate a geokriging effect\n");
      else
        outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
      }
    else
      outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
    return true;
    }

  MAP::map m = mapp->getmap();
  */

  MAP::map m;
  if (find_map(i,m) == true)
    return true;

  design_krigings.push_back(DESIGN_kriging(d,iv,m,&generaloptions,
                            equations[modnr].distrp,
                           &FC_linears[FC_linears.size()-1],
                            terms[i].options,terms[i].varnames));

  FC_nonps.push_back(FC_nonp(&master,nrlevel1,&generaloptions,equations[modnr].distrp,
                     title, pathnonp,&design_krigings[design_krigings.size()-1],
                     terms[i].options,terms[i].varnames));

  equations[modnr].add_FC(&FC_nonps[FC_nonps.size()-1],pathres);

  // variances

  make_paths(pathnonp,pathres,title,terms[i].varnames,
  "_geokriging_var.raw","variance_of_geokriging_effect_of",
  "Variance of geokriging effect of ");

  FC_nonp_variances.push_back(FC_nonp_variance(&master,nrlevel1,
                                &generaloptions,equations[modnr].distrp,
                                title,pathnonp,&design_krigings[design_krigings.size()-1],
                                &FC_nonps[FC_nonps.size()-1],terms[i].options,
                                terms[i].varnames));

  equations[modnr].add_FC(&FC_nonp_variances[FC_nonp_variances.size()-1],pathres);

  return false;
  }



bool superbayesreg::create_ridge_lasso(unsigned i)
  {
  MCMC::shrinktype sh = MCMC::ridge;
  if (terms[i].options[0]=="lasso")
    sh = MCMC::lasso;
  else if (terms[i].options[0]=="ssvs")
    sh = MCMC::ssvs;

  unsigned modnr = equations.size()-1;

  if ( ( (ridge == -1) && (terms[i].options[0] == "ridge")  ) ||
       ( (lasso == -1) && (terms[i].options[0] == "lasso")  ) ||
       ( (ssvs == -1) && (terms[i].options[0] == "ssvs")  )
     )
    {


    ST::string title;
    ST::string pathpen;
    ST::string pathpenres;

    ST::string titlevar;
    ST::string pathpenvar;
    ST::string pathpenresvar;


    ST::string h = equations[modnr].paths;

    if (terms[i].options[0] == "ridge")
      {

      title = h + ": linear effects with ridge penalty";
      pathpenres = outfile.getvalue() + "_" + h +
                   "_LinearEffects_ridgepenalty.res";
      pathpenresvar = outfile.getvalue() + "_" + h +
                   "_LinearEffects_ridgepenalty_var.res";
      titlevar = h + ": linear effects with ridge penalty (var)";

#if defined(__BUILDING_LINUX)
      pathpen = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                             + "_LinearEffects_ridgepenalty"  +
                             "_" + h + ".raw";
      pathpenvar = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                             + "_LinearEffects_ridgepenalty_var"  +
                             "_" + h + ".raw";
#else
      pathpen = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                             + "_LinearEffects_ridgepenalty"  +
                             "_" + h + ".raw";
      pathpenvar = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                             + "_LinearEffects_ridgepenalty_var"  +
                             "_" + h + ".raw";
#endif
      }
    else if (terms[i].options[0] == "lasso")
      {
      title = h + ": linear effects with lasso penalty";
      pathpenres = outfile.getvalue() + "_" + h +
                   "_LinearEffects_lassopenalty.res";
      titlevar = h + ": linear effects with lasso penalty (var)";
      pathpenresvar = outfile.getvalue() + "_" + h +
                   "_LinearEffects_lassopenalty_var.res";

#if defined(__BUILDING_LINUX)
      pathpen = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                             + "_LinearEffects_lassopenalty"  +
                             "_" + h + ".raw";
      pathpenvar = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                             + "_LinearEffects_lassopenalty_var"  +
                             "_" + h + ".raw";
#else
      pathpen = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                             + "_LinearEffects_lassopenalty"  +
                             "_" + h + ".raw";
      pathpenvar = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                             + "_LinearEffects_lassopenalty_var"  +
                             "_" + h + ".raw";
#endif
      }
    else if (terms[i].options[0] == "ssvs")
      {
      title = h + ": linear effects with ssvs prior";
      pathpenres = outfile.getvalue() + "_" + h +
                   "_LinearEffects_ssvs.res";
      pathpenresvar = outfile.getvalue() + "_" + h +
                   "_LinearEffects_ssvs_var.res";
      titlevar = h + ": variance of linear effects with ssvs prior";
//      titlevar = "";

#if defined(__BUILDING_LINUX)
      pathpen = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                             + "_LinearEffects_lassopenalty"  +
                             "_" + h + ".raw";
      pathpenvar = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                             + "_LinearEffects_lassopenalty_var"  +
                             "_" + h + ".raw";
#else
      pathpen = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                             + "_LinearEffects_ssvs"  +
                             "_" + h + ".raw";
      pathpenvar = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                             + "_LinearEffects_ssvs_var"  +
                             "_" + h + ".raw";
#endif
      }

    if (pathpen.isvalidfile() == 1)
      {
      errormessages.push_back("ERROR: unable to open file " + pathpen +
                                 " for writing\n");
      return true;
      }


    datamatrix d,iv;
    extract_data(i,d,iv,1);

    FC_linear_pens.push_back(FC_linear_pen(&master,nrlevel1,&generaloptions,
                         equations[modnr].distrp,d,terms[i].varnames,title,
                         pathpen, centerlinear.getvalue(), IWLSlineff.getvalue()));

    equations[modnr].add_FC(&FC_linear_pens[FC_linear_pens.size()-1],pathpenres);

    if (terms[i].options[0] != "ssvs")
      {

      FC_variance_pen_vectors.push_back(
      FC_variance_pen_vector(&master,&generaloptions,
                              &(FC_linear_pens[FC_linear_pens.size()-1])  ,
                              equations[modnr].distrp,titlevar, pathpenvar,
                              sh));
      }
    else
      {
      FC_variance_pen_vector_ssvss.push_back(
      FC_variance_pen_vector_ssvs(&master,nrlevel1,&generaloptions,
                              &(FC_linear_pens[FC_linear_pens.size()-1])  ,
                              equations[modnr].distrp,titlevar, pathpenvar));
      }


    if (terms[i].options[0] == "ridge")
      {
      ridge_linear = FC_linear_pens.size()-1;
      ridge = FC_variance_pen_vectors.size()-1;
      FC_variance_pen_vectors[ridge].add_variable(d,terms[i].options,
                                                  terms[i].varnames);
      }
    else if (terms[i].options[0] == "lasso")
      {
      lasso_linear = FC_linear_pens.size()-1;
      lasso = FC_variance_pen_vectors.size()-1;
      FC_variance_pen_vectors[lasso].add_variable(d,terms[i].options,
                                                  terms[i].varnames);
      }
    else
      {
      ssvs_linear = FC_linear_pens.size()-1;
      ssvs = FC_variance_pen_vector_ssvss.size()-1;
      FC_variance_pen_vector_ssvss[ssvs].add_variable(d,terms[i].options,
                                                  terms[i].varnames);
      }

    if (terms[i].options[0] != "ssvs")
      {
      equations[modnr].add_FC(&FC_variance_pen_vectors[
      FC_variance_pen_vectors.size()-1],pathpenresvar);
      }
    else
      {
//      cout << "begtest" << endl;
//      cout << "modnr: "  << modnr << endl;
//      int testsize = FC_variance_pen_vector_ssvss.size()-1;
//      cout << "size (ssvs): " << testsize << endl;
//      FC_variance_pen_vector_ssvss[testsize].outoptions();

      equations[modnr].add_FC(&FC_variance_pen_vector_ssvss[
      FC_variance_pen_vector_ssvss.size()-1],pathpenresvar);

//      cout << "size (FCs): " << equations[modnr].FCpointer.size()-1 << endl;
//      cout << "midtest" << endl;

//      testsize = equations[modnr].FCpointer.size()-1;
//      equations[modnr].FCpointer[testsize]->outoptions();

//      cout << "endtest" << endl;
      }

    }
  else
    {
    datamatrix d,iv;
    extract_data(i,d,iv,1);

    if (terms[i].options[0] == "ridge")
      {
      FC_linear_pens[ridge_linear].add_variable(d,terms[i].varnames[0]);
      FC_variance_pen_vectors[ridge].add_variable(d,terms[i].options,
                                                  terms[i].varnames);
      }
    else if (terms[i].options[0] == "lasso")
      {
      FC_linear_pens[lasso_linear].add_variable(d,terms[i].varnames[0]);
      FC_variance_pen_vectors[lasso].add_variable(d,terms[i].options,
                                                  terms[i].varnames);
      }
    else
      {
      FC_linear_pens[ssvs_linear].add_variable(d,terms[i].varnames[0]);
      FC_variance_pen_vector_ssvss[ssvs].add_variable(d,terms[i].options,
                                                  terms[i].varnames);

//      cout << "begtest" << endl;
//      cout << "modnr: "  << modnr << endl;
//      int testsize = FC_variance_pen_vector_ssvss.size()-1;
//      cout << "size: " << testsize << endl;
//      FC_variance_pen_vector_ssvss[testsize].outoptions();
//      cout << "size (FCs): " << equations[modnr].FCpointer.size()-1 << endl;

//      cout << "midtest" << endl;

//      testsize = equations[modnr].FCpointer.size()-1;
//      equations[modnr].FCpointer[testsize]->outoptions();

//      cout << "endtest" << endl;
      }

    }

  return false;
  }



bool superbayesreg::create_nonp(void)
  {

  unsigned i;
  bool error=false;

  lasso = -1;
  ridge = -1;
  ssvs = -1;
  lasso_linear = -1;
  ridge_linear = -1;
  ssvs_linear = -1;

  datamatrix datadummy = datamatrix(1,1,0.0);

  for(i=0;i<terms.size();i++)
    {
    if (terms[i].options.size() > 0)
      {
      if (terms[i].options[0] == "offset")
        create_offset(i);
      if (terms[i].options[0] == "pspline")
        create_pspline(i, datadummy);
      if (terms[i].options[0] == "pspline_merror")
        error = create_merror(i);
      if (terms[i].options[0] == "userdefined")
        error = create_userdefined(i);
      if (terms[i].options[0] == "tensor")
        error = create_userdefined_tensor(i);
      if (terms[i].options[0] == "hrandom")
        error = create_hrandom(i);
      if (terms[i].options[0] == "random")
        error = create_random(i);
      if (terms[i].options[0] == "spatial")
        error = create_mrf(i);
      if (terms[i].options[0] == "kriging")
        error = create_kriging(i);
      if (terms[i].options[0] == "geokriging")
        error = create_geokriging(i);
      if ((terms[i].options[0] == "hrandom_pspline") ||
          (terms[i].options[0] == "hrandomexp_pspline"))
        error = create_random_pspline(i);
      if ((terms[i].options[0] == "multiplicative_pspline"))
        error = create_multiplicative_pspline(i);
      if ((terms[i].options[0] == "mrf_pspline"))
        error = create_mrf_pspline(i);
      if ((terms[i].options[0] == "hrandom_mrf"))
        error = create_random_mrf(i);
      if (terms[i].options[0] == "ridge")
        error = create_ridge_lasso(i);
      if (terms[i].options[0] == "lasso")
        error = create_ridge_lasso(i);
      if (terms[i].options[0] == "ssvs")
        error = create_ridge_lasso(i);
      }

    if (error)
      return error;

    if (imeasures.getvalue() == true && FC_nonps.size() >= 1)
      FC_nonps[FC_nonps.size()-1].imeasures=true;
    }

  return false;
  }

void autocorrrun(superbayesreg & b)
  {

  if (b.resultsyesno==true)
	 {
     if (b.posteriormode == false)
       {
       ST::string path = b.outfile.getvalue()+  "_autocor" + ".raw";
       if (b.generaloptions.samplesize < unsigned(b.maxlag.getvalue()*4))
         b.outerror("ERROR: samplesize too small\n");
       else
         {
         ST::string pathgraphs = b.outfile.getvalue();
         b.simobj.autocorr(b.maxlag.getvalue(),pathgraphs);
         }
       }
     else
       b.outerror("ERROR: no MCMC simulation results\n");
	 }
  else
	 b.outerror("ERROR: no regression results\n");

  }


void getsamplerun(superbayesreg & b)
  {
  if (b.resultsyesno == true)
    {
    ST::string pathgraphs = b.outfile.getvalue();
    if (b.posteriormode == false)
      {
      #if defined(JAVA_OUTPUT_WINDOW)

// STEFAN: CHECKEN
// zweites Argument sollte ein Vektor sein.
      // b.simobj.get_samples(b.newcommands,b.outfile.getvalue() + "_");
      ST::string aString( b.outfile.getvalue() + "_" );
      b.simobj.get_samples(aString, b.newcommands);
      #else
      b.simobj.get_samples(pathgraphs);
      #endif
      }
    else
      b.outerror("ERROR: no MCMC simulation results\n");
    }
  else
    b.outerror("ERROR: no regression results\n");

  }


void superbayesreg::describe(const optionlist & globaloptions)
  {
  statobject::describe(globaloptions);
  }


#if defined(BORLAND_OUTPUT_WINDOW)
//------------------------------------------------------------------------------
#pragma package(smart_init)
#endif

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



#include "FC_merror.h"


//------------------------------------------------------------------------------
//-------------- CLASS: FC_hrandom implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{


FC_merror::FC_merror(void)
  {
  }

void FC_merror::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {
  int f;
  f = op[20].strtodouble(binning);

  f = op[74].strtodouble(a_tau2_x);
  f = op[75].strtodouble(b_tau2_x);
  f = op[76].strtodouble(m_mu_x);
  f = op[77].strtodouble(s_mu_x);
  s2_mu_x = s_mu_x*s_mu_x;
  f = op[78].strtodouble(mepropscale);
  }

FC_merror::FC_merror(GENERAL_OPTIONS * o, const ST::string & t,
            const ST::string & fp, vector<ST::string> & op,
            vector<ST::string> & vn, datamatrix & xo, datamatrix & mv,
            datamatrix & xd, FC_nonp * fcn)
     : FC(o,t,xo.rows(),1,fp)
  {
  read_options(op,vn);

  binning = fcn->beta.rows();

  xobs = xo;
  merror = (double)(xo.cols());
  FCp = fcn;

  unsigned i,j,k,l;
  mevar = mv;
  if(mevar.cols()==1)
    {
    mesd = mv;
    for(i=0; i<mesd.rows(); i++)
      mesd(i,0) = sqrt(mesd(i,0));
    }
  else
    {
    datamatrix help(xo.cols(),xo.cols(),0);
    mesd = datamatrix(xobs.rows(), 1, 0.0);
    for(i=0; i<mesd.rows(); i++)
      {
      l = xo.cols();
      for(j=0; j<xo.cols(); j++)
        {
        mesd(i,0) += mevar(i,j);
        help(j,j) = mevar(i,j);
        for(k=j+1; k<xo.cols(); k++)
          {
          help(j,k) = help(k,j) = mevar(i,l);
          l++;
          }
        }
      help = help.inverse();
      mecovinv.push_back(help);
      mesd(i,0) /= merror;
      mesd(i,0) = sqrt(mesd(i,0));
      }
    }
  xmean = datamatrix(xobs.rows(), 1, 0.0);
  for(i=0; i<xmean.rows(); i++)
    {
    for(j=0; j<merror; j++)
      xmean(i,0) += xobs(i,j);
    xmean(i,0) /= merror;
    }

//  ofstream out1("c://temp//xmean.raw");
//  xmean.prettyPrint(out1);
//  out1.close();

  FCp->designp->changingdesign=true;

  FC_tau2_x = FC(optionsp,"",1,1,samplepath + "_merror_tau2");
//  a_tau2_x = 0.001;
//  b_tau2_x = 0.001;

  FC_mu_x = FC(optionsp,"",1,1,samplepath + "_merror_mu");
//  m_mu_x = 0.0;
//  s_mu_x = 1000.0*1000.0;

  minbin = xobs.min();
  maxbin = xobs.max();
  deltabin = (maxbin-minbin)/binning;

  countmat = statmatrix<int>((unsigned)binning,1,0);

  datamatrix help = xmean;
  double h;
  double u1 = minbin+deltabin/2;
  for(unsigned i=0; i<help.rows(); i++)
    {
    if(help(i,0) < minbin)
      help(i,0) = minbin;
    if(help(i,0) > maxbin)
      help(i,0) = maxbin;
    h = floor((help(i,0) - minbin)/deltabin);
    if (h >= binning)
      {
      h -= 1.0;
      }
    help(i,0) = u1+h*deltabin;
    countmat((unsigned)h,0)++;
    }

  setbeta(help);

//  ofstream out1("c://temp//beta.raw");
//  beta.prettyPrint(out1);
//  out1.close();

//  ofstream out1("c://temp//countmat.raw");
//  countmat.prettyPrint(out1);
//  out1.close();

//  ofstream out1("c://temp//help.raw");
//  help.prettyPrint(out1);
//  out1.close();

//  ST::string help1 ="c://temp//merror0";
//  FCp->designp->test(help1);

  // 1. Indexsort of data
  FCp->designp->index_data.indexinit();
  help.indexsort(FCp->designp->index_data,0,help.rows()-1,0,0);

  //2. data = sorted observations
  double * workdata = FCp->designp->data.getV();
  int * workindex = FCp->designp->index_data.getV();
  for (j=0;j<help.rows();j++,workdata++,workindex++)
    {
    *workdata = help(*workindex,0);
    }

  // 3. Creates posbeg, posend
  int countsum=0;
  for(j=0; j<FCp->designp->posbeg.size(); j++)
    {
    if(countmat(j,0)!=0)
      {
      FCp->designp->posbeg[j] = countsum;
      FCp->designp->posend[j] = countsum + countmat(j,0)-1;
      countsum += countmat(j,0);
      }
    else
      {
      FCp->designp->posbeg[j] = -1;
      FCp->designp->posend[j] = -1;
      }
    }

  // 4. initializes ind
//  int k;
  workindex = FCp->designp->index_data.getV();
  for (j=0;j<FCp->designp->posend.size();j++)
    {
    if(FCp->designp->posbeg[j]!=-1)
      {
      for (k=FCp->designp->posbeg[j];k<=FCp->designp->posend[j];k++,workindex++)
        FCp->designp->ind(*workindex,0) = j;
      }
    }

  indexold = FCp->designp->ind;
  indexprop = indexold;

//  ST::string help2 ="c://temp//merror1";
//  FCp->designp->test(help2);

  }

FC_merror::FC_merror(const FC_merror & m)
  : FC(FC(m))
  {
  xobs = m.xobs;
  xmean = m.xmean;
  FCp = m.FCp;
  mevar = m.mevar;
  mesd = m.mesd;
  mecovinv = m.mecovinv;
  binning = m.binning;
  minbin = m.minbin;
  maxbin=m.maxbin;
  deltabin = m.deltabin;
  merror = m.merror;
  a_tau2_x = m.a_tau2_x;
  b_tau2_x = m.b_tau2_x;
  m_mu_x = m.m_mu_x;
  s_mu_x = m.s_mu_x;
  s2_mu_x = m.s2_mu_x;
  FC_mu_x = m.FC_mu_x;
  FC_tau2_x = m.FC_tau2_x;
  indexold = m.indexold;
  indexprop = m.indexprop;
  countmat = m.countmat;
  mepropscale = m.mepropscale;
  }


const FC_merror & FC_merror::operator=(const FC_merror & m)
  {
  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  xobs = m.xobs;
  xmean = m.xmean;
  FCp = m.FCp;
  mevar = m.mevar;
  mesd = m.mesd;
  mecovinv = m.mecovinv;
  binning = m.binning;
  minbin = m.minbin;
  maxbin=m.maxbin;
  deltabin = m.deltabin;
  merror = m.merror;
  a_tau2_x = m.a_tau2_x;
  b_tau2_x = m.b_tau2_x;
  m_mu_x = m.m_mu_x;
  s_mu_x = m.s_mu_x;
  s2_mu_x = m.s2_mu_x;
  FC_mu_x = m.FC_mu_x;
  FC_tau2_x = m.FC_tau2_x;
  indexold = m.indexold;
  indexprop = m.indexprop;
  countmat = m.countmat;
  mepropscale = m.mepropscale;
  return *this;
  }

void FC_merror::update(void)
  {
  unsigned i,j;
  double prop;
  double denominator = ((double)beta.rows())*s2_mu_x + FC_tau2_x.beta(0,0);
  double meanhelp = (((double)beta.rows()) * beta.mean(0) * s2_mu_x) / denominator;
  double sdhelp = sqrt((FC_tau2_x.beta(0,0) *s2_mu_x) / denominator);
  FC_mu_x.beta(0,0) = meanhelp + sdhelp * randnumbers::rand_normal();
  FC_mu_x.update();

  double ahelp = a_tau2_x + 0.5*(double(beta.rows()));
  double bhelp = b_tau2_x;
  double * betap = beta.getV();
  for(i=0; i<beta.rows(); i++, betap++)
    bhelp += 0.5*(*betap - FC_mu_x.beta(0,0))*(*betap - FC_mu_x.beta(0,0));
  FC_tau2_x.beta(0,0) = rand_invgamma(ahelp, bhelp);
  FC_tau2_x.update();
//  FC_mu_x.beta(0,0) = 0.0;
//  FC_tau2_x.beta(0,0) = 0.2;

  double sqrtM = sqrt(merror);
  double lognew, logold;
  double h, splinevalnew, splinevalold;
  double * linpredoldp;
  double * linprednewp;
  double linpred;

  if (FCp->likep->linpred_current==1)
    linpredoldp = FCp->likep->linearpred1.getV();
  else
    linpredoldp = FCp->likep->linearpred2.getV();

  double * resp = FCp->likep->response.getV();
  double * wp = FCp->likep->weight.getV();

  double priorold, priornew;
  double melikeold, melikenew;

  // reset countmat
  for(i=0; i<countmat.rows(); i++)
    countmat(i,0)=0;

  double u1 = minbin+0.5*deltabin;
  betap = beta.getV();
  unsigned * indexpropp = indexprop.getV();
  unsigned * indexoldp = indexold.getV();
  double * mesdp = mesd.getV();
  double * mevarp = mevar.getV();
  double * xobsp = xobs.getV();
  datamatrix help1(xobs.cols(), 1, 0.0);
  datamatrix help2(xobs.cols(), 1, 0.0);

  for(i=0; i<beta.rows(); i++, linpredoldp++, resp++, wp++, betap++, indexpropp++, indexoldp++, mesdp++, mevarp++)
    {
//    cout << "iter.: " << optionsp->nriter << endl;
//    cout << "obs.: " << i << endl;
    // generate proposal
    prop = *betap + mepropscale * *mesdp * randnumbers::rand_normal()/sqrtM;

//    cout << "current: " << beta(i,0) << endl;
//    cout << "proposal: " << prop << endl;

    // start binning of proposal
    if(prop < minbin)
      prop = minbin;
    if(prop > maxbin)
      prop = maxbin;
    h = floor((prop - minbin)/deltabin);
    if (h >= binning)
      h -= 1.0;
    *indexpropp = (unsigned)h;
    prop = u1+h*deltabin;

//    cout << "minbin: " << minbin << endl;
//    cout << "maxbin: " << maxbin << endl;
//    cout << "h: " << h << endl;
//    cout << "prop: " << prop << endl;
//    cout << "indexprop: " << indexprop(i,0) << endl;
//    cout << FCp->beta.rows() << " x " << FCp->beta.cols() << endl;
    // end binning of proposal

    // calculate log-likelihood
    splinevalnew = FCp->beta(*indexpropp,0);
    splinevalold = FCp->beta(*indexoldp,0);

/*    ofstream out1("c://temp//spline.raw");
    FCp->beta.prettyPrint(out1);
    out1.close();
    ofstream out2("c://temp//spline_x.raw");
    FCp->designp->data.prettyPrint(out2);
    out2.close();*/

    linpred = *linpredoldp + splinevalnew - splinevalold;
    linprednewp = &linpred;

    if(*wp != 0.0)
      {
      logold = FCp->likep->loglikelihood(resp, linpredoldp, wp);
      lognew = FCp->likep->loglikelihood(resp, linprednewp, wp);
      }
    else
      {
      logold = lognew = 0.0;
      }

/*    cout << "splinevalnew: " << splinevalnew << endl;
    cout << "splinevalold: " << splinevalold << endl;
    cout << "lognew: " << lognew << endl;
    cout << "logold: " << logold << endl;
    cout << "linprednew: " << linpred << endl;
    cout << "linpredold: " << *linpredoldp << endl;*/

     // calculate prior
    priornew = -0.5*(prop-FC_mu_x.beta(0,0))*(prop-FC_mu_x.beta(0,0))/(FC_tau2_x.beta(0,0));
    priorold = -0.5*(*betap-FC_mu_x.beta(0,0))*(*betap-FC_mu_x.beta(0,0))/(FC_tau2_x.beta(0,0));

    // calculate measurement error likelihood
    melikeold = 0.0;
    melikenew = 0.0;

    if(mevar.cols()==1)
      {
      for(j=0; j < merror; j++, xobsp++)
        {
        melikeold += (*xobsp-*betap)*(*xobsp-*betap);
        melikenew += (*xobsp-prop)*(*xobsp-prop);
        }
      melikeold *= -0.5 / *mevarp;
      melikenew *= -0.5 / *mevarp;
      }
    else
      {
      for(j=0; j<merror; j++, xobsp++)
        {
        help1(j,0) = *xobsp - *betap;
        help2(j,0) = *xobsp - prop;
        }
      melikeold = -0.5*(help1.transposed() * mecovinv[i] * help1)(0,0);
      melikenew = -0.5*(help2.transposed() * mecovinv[i] * help2)(0,0);
      }

    double logu = log(randnumbers::uniform());

    nrtrials++;
    if(logu <= lognew - logold + priornew - priorold + melikenew - melikeold)
      {
      acceptance++;
      *betap = prop;
      *indexoldp = *indexpropp;
      *linpredoldp = linpred;
      }
    else
      {

      }
//    cout << indexold(i,0) << endl;
    countmat(*indexoldp,0)++;
    }
  FC::update();

  // 1. Indexsort of data
  FCp->designp->index_data.indexinit();
  beta.indexsort(FCp->designp->index_data,0,beta.rows()-1,0,0);

  //2. data = sorted observations
  double * workdata = FCp->designp->data.getV();
  int * workindex = FCp->designp->index_data.getV();
  for (j=0;j<beta.rows();j++,workdata++,workindex++)
    {
    *workdata = beta(*workindex,0);
    }

  // 3. Creates posbeg, posend
  int countsum=0;
  for(j=0; j<FCp->designp->posbeg.size(); j++)
    {
    if(countmat(j,0)!=0)
      {
      FCp->designp->posbeg[j] = countsum;
      FCp->designp->posend[j] = countsum + countmat(j,0)-1;
      countsum += countmat(j,0);
      }
    else
      {
      FCp->designp->posbeg[j] = -1;
      FCp->designp->posend[j] = -1;
      }
    }

  // 4. initializes ind
  int k;
  workindex = FCp->designp->index_data.getV();
  for (j=0;j<FCp->designp->posend.size();j++)
    {
    if(FCp->designp->posbeg[j]!=-1)
      {
      for (k=FCp->designp->posbeg[j];k<=FCp->designp->posend[j];k++,workindex++)
        FCp->designp->ind(*workindex,0) = j;
      }
    }
  }

bool FC_merror::posteriormode(void)
  {
  return true;
  }

void FC_merror::outresults(ofstream & out_stata, ofstream & out_R, ofstream & out_R2BayesX,
                            const ST::string & pathresults)
  {
  FC::outresults(out_stata, out_R, out_R2BayesX, pathresults);
  FC::outresults_acceptance();
  optionsp->out("    Results for the true covariate values are stored in file\n");
  optionsp->out("    " +  pathresults + "\n");
  optionsp->out("\n");

  ST::string pathresults2 = pathresults.substr(0,pathresults.length()-4) +
                                 "_expectation_true_covariate_values.res";
  optionsp->out("    Expectation of true covariate values\n\n");
  FC_mu_x.outresults(out_stata, out_R, out_R2BayesX, pathresults2);
  FC_mu_x.outresults_singleparam(out_stata, out_R, pathresults2);
  optionsp->out("    Results for the variance component are stored in file\n");
  optionsp->out("    " +  pathresults2 + "\n");
  optionsp->out("\n");

  ST::string pathresults3 = pathresults.substr(0,pathresults.length()-4) +
                                 "_variance_true_covariate_values.res";
  optionsp->out("    Variance of true covariate values:\n\n");
  FC_tau2_x.outresults(out_stata, out_R, out_R2BayesX, pathresults3);
  FC_tau2_x.outresults_singleparam(out_stata, out_R, pathresults3);
  optionsp->out("    Results for the variance component are stored in file\n");
  optionsp->out("    " +  pathresults3 + "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }

void FC_merror::get_samples(const ST::string & filename, ofstream & outg) const
  {
  FC::get_samples(filename,outg);
  ST::string filed = filename.substr(0,filename.length()-11) +
                                 "_expectation_true_covariate_values_samples.raw";
  FC_mu_x.get_samples(filed, outg);
  filed = filename.substr(0,filename.length()-11) +
                                 "_variance_true_covariate_values_samples.raw";
  FC_tau2_x.get_samples(filed, outg);
  }

void FC_merror::outoptions(void)
  {
  optionsp->out("  " + title + "\n",true);
  optionsp->out("\n");
  optionsp->out("  Prior for true covariate values:\n");
  optionsp->out("    Normal distribution with mean mu and variance tau2\n");
  optionsp->out("  Hyperprior for mu:\n");
  optionsp->out("    Normal distribution with mean " + ST::doubletostring(m_mu_x,2) + " and variance " + ST::doubletostring(s2_mu_x,2) + "\n");
  optionsp->out("  Hyperprior for tau2:\n");
  optionsp->out("    Inverse gamma distribution with hyperparameters a=" + ST::doubletostring(a_tau2_x,4) + " and b=" + ST::doubletostring(b_tau2_x,4) + "\n");
  optionsp->out("  Scale parameter for the random walk proposal: "+ ST::doubletostring(mepropscale, 2) +"\n");
  optionsp->out("  Number of binning intervals: "+ ST::doubletostring(binning, 0) +"\n");
  optionsp->out("\n");
  }

} // end: namespace MCMC





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



#include"model_parameters.h"
#include<algorithm>


//------------------------------------------------------------------------------
//----------- class term_nonp: implementation of member functions --------------
//------------------------------------------------------------------------------

term_nonp::term_nonp(vector<ST::string> & na)
  {

  // DO NOT CHANGE ORDER!!!!
  termnames = na;
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  difforder =  intoption("difforder",2,1,3);
  lambda = doubleoption("lambda",10,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  nocenter = simpleoption("nocenter",false);
  map=stroption("map");
  lambda_re = doubleoption("lambda_re",10,0,10000000);
  a_re = doubleoption("a_re",0.001,-1.0,500);
  b_re = doubleoption("b_re",0.001,0,500);
  internal_mult = simpleoption("internal_mult",false);
  samplemult = simpleoption("samplemult",false);
  vector<ST::string> ctypes;
  ctypes.push_back("unconstrained");
  ctypes.push_back("increasing");
  ctypes.push_back("decreasing");
  constraints = stroption("constraints",ctypes,"unconstrained");
  round = doubleoption("round",-1,0,500);
  vector<ST::string> centermethods;
  centermethods.push_back("meancoeff");
  centermethods.push_back("meanintegral");
  centermethods.push_back("meaninvvar");
  centermethods.push_back("nullspace");
  centermethods.push_back("meansimple");      // subtract mean from parameters
  centermethods.push_back("integralsimple");      // subtract mean from parameters
  centermethods.push_back("meanf");           // sample centered f
  centermethods.push_back("meanfd");           // sample centered f
  centermethods.push_back("meansum2");
  centermethod = stroption("centermethod",centermethods,"meanfd");
  internal_multexp = simpleoption("internal_multexp",false);
  pvalue = simpleoption("pvalue",false);
  meaneffect = simpleoption("meaneffect",false);
  binning = doubleoption("binning",-1,10,1000);

  vector<ST::string> updatem;
  updatem.push_back("direct");
  updatem.push_back("orthogonal");
  update = stroption("update",updatem,"direct");

  vector<ST::string> nuv;
  nuv.push_back("0.5");
  nuv.push_back("1.5");
  nuv.push_back("2.5");
  nuv.push_back("3.5");
  nu = stroption("nu",nuv,"1.5");
  maxdist = doubleoption("maxdist",-1,0.00001,10000);

  ccovariate = simpleoption("ccovariate",false);

  sum2 = doubleoption("sum2",0,0,10000000);

  derivative = simpleoption("derivative",false);
  samplederivative = simpleoption("samplederivative",false);
  samplef = simpleoption("samplef",false);
  shrinkage = doubleoption("shrinkage",1,0,10000000);
  shrinkagefix = simpleoption("shrinkagefix",false);
  shrinkageweight = doubleoption("shrinkageweight",1,0,10000000);
  adaptiveshrinkage = simpleoption("adaptive",false);
  tau2 = doubleoption("tau2",1,0.000000000001,100000000);
  meaneffectconst = doubleoption("meaneffectconst",0,-10000000,10000000);

  vector<ST::string> priors;
  priors.push_back("iid");
  priors.push_back("lasso");
  priors.push_back("dirichlet");
  priors.push_back("nmig");
  priors.push_back("ssvs");
  priors.push_back("nn");

  prior = stroption("prior",priors,"iid");
  knotpath = fileoption("knotpath","");
  datasetref = stroption("datasetref","");

  lambdaconst = simpleoption("lambdaconst",false);

  abeta = doubleoption("abeta",1,0.00000001,500);
  bbeta = doubleoption("bbeta",1,0.00000001,500);
  r = doubleoption("r",0.000025,0.0000000001,1);
  v = doubleoption("v",5,0.0000000001,500);
  aQ = doubleoption("aQ",1,0.00000001,500);
  bQ = doubleoption("bQ",1,0.00000001,500);
  regiterates = intoption("regiterates",1000,0,1000000000);

  center = simpleoption("center",false);

  tildea = doubleoption("tildea",0.001,-1.0,500);
  tildeb = doubleoption("tildeb",0.001,0,500);
  cauchy = simpleoption("cauchy",false);
  wei = simpleoption("wei",false);
  scaletau2 = doubleoption("scaletau2",0,0,1000000000);
  r2 = doubleoption("r2",0.00456356,0,1);

  v1 = doubleoption("v1",5,0.0000000001,500);
  v2 = doubleoption("v2",25,0.0000000001,500);

  tildev1 = doubleoption("tildev1",10.5,0.0000000001,500);
  tildev2 = doubleoption("tildev2",1,0.0000000001,500);

  // gamma prior for tau2 instead of normal prior for tau in ssvs
  gig = simpleoption("gig",false);

  //test different options for proposal for scale dependent weibull  prior
  vector<ST::string> proposals;
  proposals.push_back("gamma");
  proposals.push_back("iwls_tau");
  proposals.push_back("iwls_logtau2");
  proposals.push_back("IG");

  proposal = stroption("proposal",proposals,"iwls_logtau2");

  rankK = intoption("rankK", -1, -1, 50000);
  penmatdata = stroption("penmatdata");

  //conjugate prior for effect fusion in Gaussian model (Paper Daniela Pauger and Helga Wagner, Linz)
  cprior = simpleoption("cprior",false);

  designmatdata = stroption("designmatdata");
  priormeandata = stroption("priormeandata","");

  vector<ST::string> hyperpriors;
  hyperpriors.push_back("invgamma"); // standard inverse gamma
  hyperpriors.push_back("scaledep"); // scale-dependent as proposed by Klein&Kneib, corresponds to weibull distribution
  hyperpriors.push_back("hcauchy"); // generalised beta prime based on half-cauchy for tau
  hyperpriors.push_back("hnormal");  // gamma prior based on half normal for tau
  hyperpriors.push_back("aunif");      // based on approximation for uniform prior for tau
  hyperprior = stroption("hyperprior",hyperpriors,"invgamma");

  penmatdata2 = stroption("penmatdata2");
  designmatdata2 = stroption("designmatdata2");

  constrmatdata = stroption("constrmatdata");

  nraniso = intoption("nraniso", 11, 1, 101);
  minaniso = doubleoption("minaniso", 0.05, 0.0000000001, 0.49);

  WAICoff = simpleoption("WAICoff",false);

  betastart = stroption("betastart");

  mevar = stroption("mevar");
  covdata = stroption("covdata");

  a_tau2_x = doubleoption("a_tau2_x", 0.001, -1, 500);
  b_tau2_x = doubleoption("b_tau2_x", 0.001, 0, 500);
  m_mu_x = doubleoption("m_mu_x", 0.0, -1000, 1000);
  s_mu_x = doubleoption("s_mu_x", 1000.0, 0.000000000001, 10000000000.0);
  mepropscale = doubleoption("mepropscale", 2.0, 0.000001, 1000000.0);

  reduceddesign = simpleoption("reduceddesign",false);

  vector<ST::string> ssvsupdates;
  ssvsupdates.push_back("regcoeff"); // IWLS update for tau based on interpretation as a regression coefficient
  ssvsupdates.push_back("sdev"); // IWLS update for tau based on interpretation as a standard deviation for beta
  ssvsupdates.push_back("gibbs"); // Gibbs update for tau^2 based GIG full conditional
  ssvsupdate = stroption("ssvsupdate",ssvsupdates,"regcoeff");

  }

void term_nonp::setdefault(void)
  {
  degree.setdefault();
  numberknots.setdefault();
  difforder.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  nocenter.setdefault();
  map.setdefault();
  lambda_re.setdefault();
  a_re.setdefault();
  b_re.setdefault();
  internal_mult.setdefault();
  samplemult.setdefault();
  constraints.setdefault();
  round.setdefault();
  centermethod.setdefault();
  internal_multexp.setdefault();
  pvalue.setdefault();
  meaneffect.setdefault();
  binning.setdefault();
  update.setdefault();
  nu.setdefault();
  maxdist.setdefault();
  ccovariate.setdefault();
  sum2.setdefault();
  derivative.setdefault();
  samplederivative.setdefault();
  samplef.setdefault();
  shrinkage.setdefault();
  shrinkagefix.setdefault();
  shrinkageweight.setdefault();
  adaptiveshrinkage.setdefault();
  tau2.setdefault();
  meaneffectconst.setdefault();
  prior.setdefault();
  knotpath.setdefault();
  datasetref.setdefault();
  lambdaconst.setdefault();

  abeta.setdefault();
  bbeta.setdefault();
  r.setdefault();
  v.setdefault();
  aQ.setdefault();
  bQ.setdefault();
  regiterates.setdefault();
  center.setdefault();

  tildea.setdefault();
  tildeb.setdefault();
  cauchy.setdefault();
  wei.setdefault();
  scaletau2.setdefault();
  r2.setdefault();

  v1.setdefault();
  v2.setdefault();
  tildev1.setdefault();
  tildev2.setdefault();
  gig.setdefault();
  proposal.setdefault();

  rankK.setdefault();
  penmatdata.setdefault();

  cprior.setdefault();

  designmatdata.setdefault();
  priormeandata.setdefault();

  hyperprior.setdefault();

  penmatdata2.setdefault();
  designmatdata2.setdefault();

  constrmatdata.setdefault();

  nraniso.setdefault();
  minaniso.setdefault();

  WAICoff.setdefault();

  betastart.setdefault();

  mevar.setdefault();
  covdata.setdefault();

  a_tau2_x.setdefault();
  b_tau2_x.setdefault();
  m_mu_x.setdefault();
  s_mu_x.setdefault();
  mepropscale.setdefault();

  reduceddesign.setdefault();

  ssvsupdate.setdefault();
  }


bool term_nonp::checkvector(const vector<term> & terms,const unsigned & i)
  {
  assert(i< terms.size());

  bool f = false;
  unsigned j=0;
  while ( (j<termnames.size()) && (f == false) )
    {
    if (terms[i].type == termnames[j])
      {
      f = true;
      }
    j ++;
    }

  return f;
  }


bool term_nonp::check(term & t)
  {

  if ( (t.varnames.size()<=2)  && (t.options.size() >= 1)
        && (t.options.size() <= 200) )
    {

    bool f = false;
    unsigned j=0;
    unsigned namespos=0;

    while ( (j<termnames.size()) && (f == false) )
      {
      if (t.options[0] == termnames[j])
        {
        namespos=j;
        f = true;
        }
      j ++;
      }

    if (f==false)
      {
      setdefault();
      return false;
      }

    optionlist optlist;
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&difforder);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&nocenter);
    optlist.push_back(&map);
    optlist.push_back(&lambda_re);
    optlist.push_back(&a_re);
    optlist.push_back(&b_re);
    optlist.push_back(&internal_mult);
    optlist.push_back(&samplemult);
    optlist.push_back(&constraints);
    optlist.push_back(&round);
    optlist.push_back(&centermethod);
    optlist.push_back(&internal_multexp);
    optlist.push_back(&pvalue);
    optlist.push_back(&meaneffect);
    optlist.push_back(&binning);
    optlist.push_back(&update);
    optlist.push_back(&nu);
    optlist.push_back(&maxdist);
    optlist.push_back(&ccovariate);
    optlist.push_back(&sum2);
    optlist.push_back(&derivative);
    optlist.push_back(&samplederivative);
    optlist.push_back(&samplef);
    optlist.push_back(&shrinkage);
    optlist.push_back(&shrinkagefix);
    optlist.push_back(&shrinkageweight);
    optlist.push_back(&adaptiveshrinkage);
    optlist.push_back(&tau2);
    optlist.push_back(&meaneffectconst);
    optlist.push_back(&prior);
    optlist.push_back(&knotpath);
    optlist.push_back(&datasetref);
    optlist.push_back(&lambdaconst);


    optlist.push_back(&abeta);
    optlist.push_back(&bbeta);
    optlist.push_back(&r);
    optlist.push_back(&v);
    optlist.push_back(&aQ);
    optlist.push_back(&bQ);
    optlist.push_back(&regiterates);
    optlist.push_back(&center);

    optlist.push_back(&tildea);
    optlist.push_back(&tildeb);
    optlist.push_back(&cauchy);
    optlist.push_back(&wei);
	  optlist.push_back(&scaletau2);
	  optlist.push_back(&r2);


    optlist.push_back(&v1);
    optlist.push_back(&v2);
    optlist.push_back(&tildev1);
    optlist.push_back(&tildev2);
    optlist.push_back(&gig);
	optlist.push_back(&proposal);

	optlist.push_back(&rankK);
	optlist.push_back(&penmatdata);

	optlist.push_back(&cprior);

	optlist.push_back(&designmatdata);
	optlist.push_back(&priormeandata);

	optlist.push_back(&hyperprior);

	optlist.push_back(&penmatdata2);
	optlist.push_back(&designmatdata2);

	optlist.push_back(&constrmatdata);
	optlist.push_back(&nraniso);
    optlist.push_back(&WAICoff);

	optlist.push_back(&minaniso);

    optlist.push_back(&betastart);

    optlist.push_back(&mevar);
    optlist.push_back(&covdata);

    optlist.push_back(&a_tau2_x);
    optlist.push_back(&b_tau2_x);
    optlist.push_back(&m_mu_x);
    optlist.push_back(&s_mu_x);
    optlist.push_back(&mepropscale);

    optlist.push_back(&reduceddesign);

	optlist.push_back(&ssvsupdate);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(200);
    t.options[0] = termnames[namespos];
    t.options[1] = ST::inttostring(degree.getvalue());
    t.options[2] = ST::inttostring(numberknots.getvalue());
    t.options[3] = ST::inttostring(difforder.getvalue());
    t.options[4] = ST::doubletostring(lambda.getvalue());
    t.options[5] = ST::doubletostring(a.getvalue());
    t.options[6] = ST::doubletostring(b.getvalue());

    if(nocenter.getvalue() == false)
      t.options[7] = "false";
    else
      t.options[7] = "true";

    t.options[8] = map.getvalue();

    t.options[9] = ST::doubletostring(lambda_re.getvalue());

    t.options[10] = ST::doubletostring(a_re.getvalue());

    t.options[11] = ST::doubletostring(b_re.getvalue());

    if(internal_mult.getvalue() == false)
      t.options[12] = "false";
    else
      t.options[12] = "true";

    if(samplemult.getvalue() == false)
      t.options[13] = "false";
    else
      t.options[13] = "true";

    t.options[14] = constraints.getvalue();

    t.options[15] = ST::doubletostring(round.getvalue());

    t.options[16] = centermethod.getvalue();

    if(internal_multexp.getvalue() == false)
      t.options[17] = "false";
    else
      t.options[17] = "true";

    if(pvalue.getvalue() == false)
      t.options[18] = "false";
    else
      t.options[18] = "true";

    if(meaneffect.getvalue() == false)
      t.options[19] = "false";
    else
      t.options[19] = "true";

    t.options[20] = ST::doubletostring(binning.getvalue());

    t.options[21] = update.getvalue();

    t.options[22] = nu.getvalue();

    t.options[23] = ST::doubletostring(maxdist.getvalue());

    if(ccovariate.getvalue() == false)
      t.options[24] = "false";
    else
      t.options[24] = "true";

    t.options[25] = ST::doubletostring(sum2.getvalue());

    if(derivative.getvalue() == false)
      t.options[26] = "false";
    else
      t.options[26] = "true";

    if(samplederivative.getvalue() == false)
      t.options[27] = "false";
    else
      t.options[27] = "true";

    if(samplef.getvalue() == false)
      t.options[28] = "false";
    else
      t.options[28] = "true";

    t.options[29] = ST::doubletostring(shrinkage.getvalue());

    if(shrinkagefix.getvalue() == false)
      t.options[30] = "false";
    else
      t.options[30] = "true";

    t.options[31] = ST::doubletostring(shrinkageweight.getvalue());

    if(adaptiveshrinkage.getvalue() == false)
      t.options[32] = "false";
    else
      t.options[32] = "true";

    t.options[33] = ST::doubletostring(tau2.getvalue());

    t.options[34] = ST::doubletostring(meaneffectconst.getvalue());

    t.options[35] = prior.getvalue();

    t.options[36] = knotpath.getvalue();
    t.options[37] = datasetref.getvalue();

    if(lambdaconst.getvalue() == false)
      t.options[38] = "false";
    else
      t.options[38] = "true";

    t.options[39] = ST::doubletostring(abeta.getvalue());
    t.options[40] = ST::doubletostring(bbeta.getvalue());
    t.options[41] = ST::doubletostring(r.getvalue());
    t.options[42] = ST::doubletostring(v.getvalue());
    t.options[43] = ST::doubletostring(aQ.getvalue());
    t.options[44] = ST::doubletostring(bQ.getvalue());
    t.options[45] = ST::inttostring(regiterates.getvalue());

    if (center.getvalue() == false)
      t.options[46] = "false";
    else
      t.options[46] = "true";

    t.options[47] = ST::doubletostring(tildea.getvalue());
    t.options[48] = ST::doubletostring(tildeb.getvalue());
    if (cauchy.getvalue() == false)
      t.options[49] = "false";
    else
      t.options[49] = "true";
    if (wei.getvalue() == false)
      t.options[50] = "false";
    else
      t.options[50] = "true";

    t.options[51] = ST::doubletostring(scaletau2.getvalue());
    t.options[52] = ST::doubletostring(r2.getvalue());

    // pararameters for ssvs psi2 (IG prior)
    t.options[53] = ST::doubletostring(v1.getvalue());
    t.options[54] = ST::doubletostring(v2.getvalue());

    // pararameters for ssvs if tau2 has gig FC
    t.options[55] = ST::doubletostring(tildev1.getvalue());
    t.options[56] = ST::doubletostring(tildev2.getvalue());
    //normal prior for tau or (if true gig prior for tau2 in class ssvs
    if (gig.getvalue() == false)
      t.options[57] = "false";
    else
      t.options[57] = "true";

    //type of proposal for weibull prior
    t.options[58] = proposal.getvalue();

    t.options[59] = ST::inttostring(rankK.getvalue());
    t.options[60] = penmatdata.getvalue();

	//conjugate prior for effect fusion ((paper Daniela Pauger and Helga Wagner, Linz)
    if (cprior.getvalue() == false)
      t.options[61] = "false";
    else
      t.options[61]= "true";

    t.options[62] = designmatdata.getvalue();
    t.options[63] = priormeandata.getvalue();

    t.options[64] = hyperprior.getvalue();

    t.options[65] = penmatdata2.getvalue();
    t.options[66] = designmatdata2.getvalue();

    t.options[67] = constrmatdata.getvalue();
    t.options[68] = ST::inttostring(nraniso.getvalue());
    if (WAICoff.getvalue() == false)
      t.options[69] = "false";
    else
      t.options[69] = "true";

    t.options[70] = ST::doubletostring(minaniso.getvalue());

    t.options[71] = betastart.getvalue();

    t.options[72] = mevar.getvalue();
    t.options[73] = covdata.getvalue();

    t.options[74] = ST::doubletostring(a_tau2_x.getvalue());
    t.options[75] = ST::doubletostring(b_tau2_x.getvalue());
    t.options[76] = ST::doubletostring(m_mu_x.getvalue());
    t.options[77] = ST::doubletostring(s_mu_x.getvalue());
    t.options[78] = ST::doubletostring(mepropscale.getvalue());

    if (reduceddesign.getvalue() == false)
      t.options[79] = "false";
    else
      t.options[79] = "true";

    t.options[80] = ssvsupdate.getvalue();

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


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



#include "design_userdefined.h"

using std::deque;

namespace MCMC
{


//------------------------------------------------------------------------------
//--------- CLASS: DESIGN_userdefined implementation of member functions -------
//------------------------------------------------------------------------------

void DESIGN_userdefined::compute_f(datamatrix & beta,datamatrix & betalin,
                       datamatrix & f, datamatrix & ftot)
  {
  unsigned i,j;

  double * workf = f.getV();

  for(i=0; i<Zout2.size(); i++, workf++)
    {
    *workf = 0.0;
    for(j=0; j<Zout2[i].size(); j++)
      {
      *workf += Zout2[i][j]*beta(index_Zout2[i][j],0);
      }
    }
  }

void DESIGN_userdefined::read_options(vector<ST::string> & op,
                                 vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  13      samplemult
  14      constraints
  15      round
  16      centermethod
  17      internal_mult
  18      pvalue
  19      meaneffect
  20      binning
  21      update
  22      nu
  23      maxdist
  24      ccovariate
  */


  int f;

  if (op[7] == "false")   //nocenter==false, i.e. center
    center = true;
  else
    center = false;



  f = op[15].strtodouble(round);

  if (op[16]=="meancoeff")
    centermethod = meancoeff;
  else if (op[16] == "meansimple")
    centermethod = meansimple;
  else if (op[16] == "integralsimple")
    centermethod = integralsimple;
  else if (op[16] == "nullspace")
    centermethod = nullspace;
  else if (op[16] == "meaninvvar")
    centermethod = cmeaninvvar;
  else if (op[16] == "meanintegral")
    centermethod = cmeanintegral;
  else if (op[16] == "meanf")
    centermethod = meanf;
  else if (op[16] == "meanfd")
    centermethod = meanfd;
  else if (op[16] == "meansum2")
    centermethod = meansum2;
  else if (op[16] == "userdefined")
    centermethod = userdefined;

  f = op[20].strtodouble(binning);

  f = op[59].strtodouble(rankK);

  datanames = vn;
  }


DESIGN_userdefined::DESIGN_userdefined(void) : DESIGN()
  {

  }


void DESIGN_userdefined::init_data(const datamatrix & dm, const datamatrix & iv)
  {

  unsigned j;

  // 1. Indexsort of data
  if (index_data.rows() <= 1)
    {
    index_data = statmatrix<int>(dm.rows(),1);
    index_data.indexinit();
    dm.indexsort(index_data,0,dm.rows()-1,0,0);
    }

  double dm_mean = dm.mean(0);

  //2. data = sorted observations, init intvar and intvar2
  data = datamatrix(dm.rows(),1);
  double * workdata = data.getV();
  int * workindex = index_data.getV();
  for (j=0;j<dm.rows();j++,workdata++,workindex++)
    {
    *workdata = dm(*workindex,0);
    }

  if (iv.rows() == dm.rows())
    {
    intvar = iv;
    intvar2 = datamatrix(iv.rows(),1);
    double * workintvar2 = intvar2.getV();
    double * workintvar = intvar.getV();

    for (j=0;j<iv.rows();j++,workintvar++,workintvar2++)
      {
      *workintvar2 = pow(*workintvar,2);
      }
    }

  // 3. Creates posbeg, posend

  posbeg.erase(posbeg.begin(),posbeg.end());
  posend.erase(posend.begin(),posend.end());
  posbeg.push_back(0);
  workdata = data.getV()+1;
  double help = data(0,0);
  for(j=1;j<data.rows();j++,workdata++)
    {
    if (  *workdata != help)
      {
      posend.push_back(j-1);
      if (j < data.rows())
        posbeg.push_back(j);
      }

    help = *workdata;

    }

  if (posend.size() < posbeg.size())
    posend.push_back(data.rows()-1);

  // 4. initializes ind
  int k;
  workindex = index_data.getV();
  ind = statmatrix<unsigned>(dm.rows(),1);
  for (j=0;j<posend.size();j++)
    {
    for (k=posbeg[j];k<=posend[j];k++,workindex++)
      ind(*workindex,0) = j;

    }

/*    j=0;
    ofstream out2("c:\\temp\\userdefined_spatial\\posend.res");
    for (j=0;j<posend.size();j++)
      out2 << posend[j] << endl;

    ofstream out2a("c:\\temp\\userdefined_spatial\\posbeg.res");
    for (j=0;j<posbeg.size();j++)
      out2a << posbeg[j] << endl;*/

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\ind.res");
  // ind.prettyPrint(out);
  // TEST


  // 5. Compute meaneffectnr, mclosest, effectvalues
  effectvalues.erase(effectvalues.begin(),effectvalues.end());
  double d;
  meaneffectnr = 0;
  double mclosest = data(posbeg[0],0);
  for(j=0;j<posbeg.size();j++)
    {
    d = data(posbeg[j],0);
    if ( fabs(d-dm_mean) < fabs(mclosest-dm_mean) )
      {
      meaneffectnr = j;
      mclosest = d;
      }

    effectvalues.push_back(ST::doubletostring(d,0));
    }

  compute_meaneffectintvar();
  }


void DESIGN_userdefined::compute_precision(double l)
  {
  if (precisiondeclared==false)
    {
    precision = envmatdouble(K.computeMaxXenv(XWX),1.0,nrpar);
    precisiondeclared = true;
    }


  precision.addto(XWX,K,1.0,l);

  //ofstream out("c:\\temp\\K.res");
  //K.print2(out);


  // ofstream out("c:\\temp\\precision.res");
  // precision.print2(out);


  }

void DESIGN_userdefined::compute_basisNull(void)
  {
  if (center==true)
  {
  unsigned i,j;

  if (centermethod==meancoeff || centermethod==meansimple)
    {
    basisNull = datamatrix(1,nrpar,1);
    position_lin = -1;
    }
  else if (centermethod==cmeaninvvar)
    {
    }
  else if ((centermethod==cmeanintegral) || (centermethod==integralsimple))  // integral f = 0
    {
    }
  else if (centermethod==meanf)            // sum of f's zero (over all observations)
    {

    basisNull = datamatrix(1,nrpar,1);

    unsigned k;
    for (k=0;k<nrpar;k++)
      basisNull(0,k) = compute_sumBk(k);

    position_lin = -1;

    }
  else if (centermethod==meanfd)            // sum of f's zero
    {

    basisNull = datamatrix(1,nrpar,1);

    unsigned k;
    for (k=0;k<nrpar;k++)
      basisNull(0,k) = compute_sumBk_different(k);

    // TEST
    // ofstream out("c:\\bayesx\\testh\\results\\basisnull.res");
    // basisNull.prettyPrint(out);
    // ende: TEST

    position_lin = -1;

    }
  else if (centermethod==nullspace)
    {
    if(nrpar-rankK > 0)
      {
      basisNull = datamatrix(nrpar-rankK,nrpar,1);
      datamatrix Kstat=Kdatamat;
      datamatrix vals(Kstat.rows(),1,0);
      // bool eigentest=eigen2(Kstat,vals);
      eigensort(vals,Kstat);
      unsigned j,k;
      for(j=0; j<vals.rows(); j++)
        {
        for(k=0; k<nrpar-rankK; k++)
          {
          basisNull(k,j) = Kstat(j,rankK+k);
          }
        }
      position_lin = -1;
      }
    }


  for(i=0;i<basisNull.rows();i++)
    {
    basisNullt.push_back(datamatrix(basisNull.cols(),1));
    for(j=0;j<basisNull.cols();j++)
      basisNullt[i](j,0) = basisNull(i,j);
    }


  if (basisNull.rows() > 1)
    {
    designlinear = datamatrix(posbeg.size(),basisNull.rows()-1);

    double * workdl = designlinear.getV();
    double h;
    for(i=0;i<posbeg.size();i++)
      for(j=0;j<designlinear.cols();j++,workdl++)
        {
        h = data(posbeg[i],0);
        *workdl =  pow(static_cast<double>(h),static_cast<double>(j+1));
        }
    }

  }

  // TEST
  /*
    ofstream out("c:\\bayesx\\test\\results\\data.res");
    data.prettyPrint(out);

    ofstream out2("c:\\bayesx\\test\\results\\designlin.res");
    designlinear.prettyPrint(out2);
   */
  // TEST

  }


void DESIGN_userdefined::compute_Zout(datamatrix & Z)
  {
  unsigned i,j;

  if(Z.rows()<data.rows())
    {
    for(i=0; i<posbeg.size() ; i++)
      {
      Zout2.push_back(vector<double>());
      index_Zout2.push_back(vector<int>());
      for(j=0; j<Z.cols(); j++)
        {
        if(Z(i,j)!=0)
          {
          Zout2[i].push_back(Z(i,j));
          index_Zout2[i].push_back(j);
          }
        }
      }
    }
  else
    {
    for(i=0; i<posbeg.size() ; i++)
      {
      Zout2.push_back(vector<double>());
      index_Zout2.push_back(vector<int>());
      for(j=0; j<Z.cols(); j++)
        {
        if(Z(index_data(posbeg[i],0),j)!=0)
          {
          Zout2[i].push_back(Z(index_data(posbeg[i],0),j));
          index_Zout2[i].push_back(j);
          }
        }
      }
    }

/*
  ofstream out("c:\\temp\\Zout.res");
  for (i=0;i<Zout2.size();i++)
    {
    for(j=0;j<Zout2[i].size();j++)
      out <<  Zout2[i][j] << "  ";
    out << endl;
    }

  ofstream out2("c:\\temp\\Zout_index.res");
  for (i=0;i<index_Zout2.size();i++)
    {
    for(j=0;j<index_Zout2[i].size();j++)
      out2 <<  index_Zout2[i][j] << "  ";
    out2 << endl;
    }


  ofstream out4("c:\\temp\\Zoutmat.res");
  datamatrix Zhelp(Zout2.size(),nrpar,0);
  for (i=0;i<Zout2.size();i++)
    {
    for(j=0;j<nrpar;j++)
      Zhelp(i,j) = Zout2[i][j];
    }

  Zhelp.prettyPrint(out4);
*/

  }


void DESIGN_userdefined::compute_Zout_transposed_vector(void)
  {

  vector<double> h;
  ZoutT = vector<vector<double> >(nrpar,h);

  vector<int> h2;
  index_ZoutT = vector<vector<int> >(nrpar,h2);

  unsigned i,j;

  for(i=0;i<Zout2.size();i++)
    for(j=0;j<Zout2[i].size();j++)
      {
      ZoutT[index_Zout2[i][j]].push_back(Zout2[i][j]);
      index_ZoutT[index_Zout2[i][j]].push_back(i);
      }


/*
  ofstream out("c:\\temp\\ZoutT.res");
  for (i=0;i<ZoutT.size();i++)
    {
    for(j=0;j<ZoutT[i].size();j++)
      out <<  ZoutT[i][j] << "  ";
    out << endl;
    }

  ofstream out2("c:\\temp\\ZoutT_index.res");
  for (i=0;i<index_ZoutT.size();i++)
    {
    for(j=0;j<index_ZoutT[i].size();j++)
      out2 <<  index_ZoutT[i][j] << "  ";
    out2 << endl;
    }

  ofstream out4("c:\\temp\\ZoutT.res");
  datamatrix Zhelp(Zout2.size(),ZoutT.size(),0);
  for (i=0;i<ZoutT.size();i++)
    {
    for(j=0;j<ZoutT[i].size();j++)
      Zhelp(index_ZoutT[i][j],i) = ZoutT[i][j];
    }

  Zhelp.prettyPrint(out4);
*/
  }

  // CONSTRUCTOR

DESIGN_userdefined::DESIGN_userdefined(GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl)
                      : DESIGN(o,dp,fcl)
  {
  }

DESIGN_userdefined::DESIGN_userdefined(datamatrix & dm,datamatrix & iv,
                       datamatrix & designmat, datamatrix & penmat, datamatrix & priormean,  datamatrix & constrmat,
                       GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
                       vector<ST::string> & op,
                       vector<ST::string> & vn)
                      : DESIGN(o,dp,fcl)
  {


  read_options(op,vn);

  // Adjust to multiple columns in dm
  init_data(dm,iv);

  nrpar = designmat.cols();

  mK = penmat*priormean;
  Kdatamat = penmat;

  if(centermethod==userdefined)
    {
    optionsp->out("\nNOTE: Effect will be centered according to the provided constraint matrix.\n\n");
    basisNull = constrmat;
    }

  compute_Zout(designmat);
  compute_Zout_transposed_vector();
  K = envmatdouble(penmat, 0.00000001);
  if(rankK==-1)
    {
    // compute rankK;
    datamatrix Kstat=penmat;
    datamatrix vals(Kstat.rows(),1,0);
    bool eigentest=eigen2(Kstat,vals);
    if(eigentest==false)
      {
      rankK = Kstat.rows();
      optionsp->out("WARNING: Unable to compute rank of penalty matrix.\n");
      optionsp->out("         Rank was set to the dimension of the penalty matrix: " + ST::inttostring(rankK) + "\n");
      optionsp->out("         Please specify argument rankK\n");
      if(centermethod==nullspace)
        {
        optionsp->out("         Option centermethod was changed to meanfd\n");
        centermethod = meanfd;
        }
      }
    else
      {
      rankK=0;
      unsigned j;
      for(j=0; j<vals.rows(); j++)
        {
        if(fabs(vals(j,0))>=0.000001)
          rankK++;
        }
      }
    }

  Wsum = datamatrix(posbeg.size(),1,1);

  datamatrix help = designmat.transposed()*designmat;
  XWX = envmatdouble(help, 0.0);
  compute_XtransposedWX();
  XWres = datamatrix(nrpar,1);

  compute_precision(1.0);
  compute_basisNull();
  }


  // COPY CONSTRUCTOR

DESIGN_userdefined::DESIGN_userdefined(const DESIGN_userdefined & m)
    : DESIGN(DESIGN(m))
  {
  round = m.round;
  binning = m.binning;
  Zout2 = m.Zout2;
  index_Zout2 = m.index_Zout2;
  mK = m.mK;
  Kdatamat = m.Kdatamat;
  }


  // OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_userdefined & DESIGN_userdefined::operator=(const DESIGN_userdefined & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));
  round = m.round;
  binning = m.binning;
  Zout2 = m.Zout2;
  index_Zout2 = m.index_Zout2;
  mK = m.mK;
  Kdatamat = m.Kdatamat;
  return *this;
  }


void DESIGN_userdefined::outbasis_R(ofstream & out)
  {
  }

void DESIGN_userdefined::outoptions(GENERAL_OPTIONS * op)
  {
  ST::string centerm;

  if (center==false)
    centerm = "uncentered sampling";
  else
    centerm = "centered sampling";

  op->out("  Prior: user defined\n");
  op->out("  Rank of penalty matrix: " + ST::inttostring(rankK) + "\n");
  op->out("  " + centerm + "\n" );

  op->out("\n");

  }

void DESIGN_userdefined::compute_XtransposedWres(datamatrix & partres, double l, double t2)
  {

  unsigned i,j;

  if (ZoutT.size() != nrpar)
    compute_Zout_transposed();

  if (consecutive_ZoutT == -1)
    {
    bool c = check_ZoutT_consecutive();
    consecutive_ZoutT = c;
    }

  double * workXWres = XWres.getV();
  double * workmK = mK.getV();
  unsigned size;
    vector<double>::iterator wZoutT;

  if (consecutive_ZoutT == 0)
    {

    vector<int>::iterator wZoutT_index;

    for(i=0;i<nrpar;i++,workXWres++,workmK++)
      {
      *workXWres=0;
      size = ZoutT[i].size();
      if(size>0)
        {
        wZoutT = ZoutT[i].begin();
        wZoutT_index = index_ZoutT[i].begin();
        for (j=0;j<size;j++,++wZoutT,++wZoutT_index)
          {
          *workXWres+= (*wZoutT)* partres(*wZoutT_index,0);
          }
        }
      *workXWres += *workmK/t2;
      }
    }
  else
    {
    double * wpartres;

    for(i=0;i<nrpar;i++,workXWres++,workmK++)
      {
      *workXWres=0;
      size = ZoutT[i].size();
      if(size>0)
        {
        wZoutT = ZoutT[i].begin();
        wpartres = partres.getV()+index_ZoutT[i][0];
        for (j=0;j<size;j++,++wZoutT,wpartres++)
          {
          *workXWres+= (*wZoutT)* (*wpartres);
          }
        }
      *workXWres += *workmK/t2;
      }
    }

  XWres_p = &XWres;
  XWX_p = &XWX;

  // TEST
  //ofstream out("c:\\temp\\XWres.res");
  //XWres.prettyPrint(out);
  // TEST
  }


void DESIGN_userdefined::compute_Zout_transposed(void)
  {

  vector<double> h;
  ZoutT = vector<vector<double> >(nrpar,h);

  vector<int> h2;
  index_ZoutT = vector<vector<int> >(nrpar,h2);


  unsigned i,j;

  for (i=0;i<Zout.rows();i++)
    for(j=0;j<Zout.cols();j++)
      {
      ZoutT[index_Zout(i,j)].push_back(Zout(i,j));
      index_ZoutT[index_Zout(i,j)].push_back(i);
      }
  }

//------------------------------------------------------------------------------
//---- CLASS: DESIGN_userdefined_tensor implementation of member functions -----
//------------------------------------------------------------------------------

void DESIGN_userdefined_tensor::read_options(vector<ST::string> & op,
                                 vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  13      samplemult
  14      constraints
  15      round
  16      centermethod
  17      internal_mult
  18      pvalue
  19      meaneffect
  20      binning
  21      update
  22      nu
  23      maxdist
  24      ccovariate
  */

  int f;

  if (op[7] == "false")   //nocenter==false, i.e. center
    center = true;
  else
    center = false;

  f = op[15].strtodouble(round);

  if (op[16]=="meancoeff")
    centermethod = meancoeff;
  else if (op[16] == "meansimple")
    centermethod = meansimple;
  else if (op[16] == "integralsimple")
    centermethod = integralsimple;
  else if (op[16] == "nullspace")
    centermethod = nullspace;
  else if (op[16] == "meaninvvar")
    centermethod = cmeaninvvar;
  else if (op[16] == "meanintegral")
    centermethod = cmeanintegral;
  else if (op[16] == "meanf")
    centermethod = meanf;
  else if (op[16] == "meanfd")
    centermethod = meanfd;
  else if (op[16] == "meansum2")
    centermethod = meansum2;
  else if (op[16] == "userdefined")
    centermethod = userdefined;

  f = op[20].strtodouble(binning);

  f = op[59].strtodouble(rankK);

  f = op[68].strtolong(nromega);

  f = op[70].strtodouble(minomega);

  datanames = vn;
  }


DESIGN_userdefined_tensor::DESIGN_userdefined_tensor(void) : DESIGN_userdefined()
  {

  }


void DESIGN_userdefined_tensor::init_data(const datamatrix & dm, const datamatrix & iv)
  {

  unsigned j;

  if(dm.cols()==1)
    {
      // 1. Indexsort of data
    if (index_data.rows() <= 1)
      {
      index_data = statmatrix<int>(dm.rows(),1);
      index_data.indexinit();
      dm.indexsort(index_data,0,dm.rows()-1,0,0);
      }

    double dm_mean = dm.mean(0);

    //2. data = sorted observations, init intvar and intvar2
    data = datamatrix(dm.rows(),1);
    double * workdata = data.getV();
    int * workindex = index_data.getV();
    for (j=0;j<dm.rows();j++,workdata++,workindex++)
      {
      *workdata = dm(*workindex,0);
      }

    if (iv.rows() == dm.rows())
      {
      intvar = iv;
      intvar2 = datamatrix(iv.rows(),1);
      double * workintvar2 = intvar2.getV();
      double * workintvar = intvar.getV();

      for (j=0;j<iv.rows();j++,workintvar++,workintvar2++)
        {
        *workintvar2 = pow(*workintvar,2);
        }
      }

    // 3. Creates posbeg, posend

    posbeg.erase(posbeg.begin(),posbeg.end());
    posend.erase(posend.begin(),posend.end());
    posbeg.push_back(0);
    workdata = data.getV()+1;
    double help = data(0,0);
    for(j=1;j<data.rows();j++,workdata++)
      {
      if (  *workdata != help)
        {
        posend.push_back(j-1);
        if (j < data.rows())
          posbeg.push_back(j);
        }

      help = *workdata;

      }

    if (posend.size() < posbeg.size())
      posend.push_back(data.rows()-1);

    // 4. initializes ind
    int k;
    workindex = index_data.getV();
    ind = statmatrix<unsigned>(dm.rows(),1);
    for (j=0;j<posend.size();j++)
      {
      for (k=posbeg[j];k<=posend[j];k++,workindex++)
        ind(*workindex,0) = j;

      }

  /*    j=0;
      ofstream out2("c:\\temp\\userdefined_spatial\\posend.res");
      for (j=0;j<posend.size();j++)
        out2 << posend[j] << endl;

      ofstream out2a("c:\\temp\\userdefined_spatial\\posbeg.res");
      for (j=0;j<posbeg.size();j++)
        out2a << posbeg[j] << endl;*/

    // TEST
    // ofstream out("c:\\bayesx\\testh\\results\\ind.res");
    // ind.prettyPrint(out);
    // TEST


    // 5. Compute meaneffectnr, mclosest, effectvalues
    effectvalues.erase(effectvalues.begin(),effectvalues.end());
    double d;
    meaneffectnr = 0;
    double mclosest = data(posbeg[0],0);
    for(j=0;j<posbeg.size();j++)
      {
      d = data(posbeg[j],0);
      if ( fabs(d-dm_mean) < fabs(mclosest-dm_mean) )
        {
        meaneffectnr = j;
        mclosest = d;
        }

      effectvalues.push_back(ST::doubletostring(d,0));
      }

    compute_meaneffectintvar();
    }
  else if(dm.cols()==2)
    {

    // 1. Indexsort of data
    if (index_data.rows() <= 1)
      {
      index_data = statmatrix<int>(dm.rows(),1);
      index_data.indexinit();
      dm.indexsort2d(index_data,0,dm.rows()-1,0,1,0);
      }

//  double dm_mean = dm.mean(0);

    //2. data = sorted observations, init intvar and intvar2

    data = datamatrix(dm.rows(),2);
    double * workdata = data.getV();
    int * workindex = index_data.getV();
    for (j=0;j<dm.rows();j++,workdata++,workindex++)
      {
      *workdata = dm(*workindex,0);
      workdata++;
      *workdata = dm(*workindex,1);
      }

    if (iv.rows() == dm.rows())
      {
      intvar = iv;
      intvar2 = datamatrix(iv.rows(),1);
      double * workintvar2 = intvar2.getV();
      double * workintvar = intvar.getV();

      for (j=0;j<iv.rows();j++,workintvar++,workintvar2++)
        {
        *workintvar2 = pow(*workintvar,2);
        }
      }


    // 3. Creates posbeg, posend

    posbeg.erase(posbeg.begin(),posbeg.end());
    posend.erase(posend.begin(),posend.end());
    xvalues.erase(xvalues.begin(),xvalues.end());
    yvalues.erase(yvalues.begin(),yvalues.end());
    posbeg.push_back(0);
    double help1 = data(0,0);
    double help2 = data(0,1);
    xvalues.push_back(help1);
    yvalues.push_back(help2);
    for(j=1;j<data.rows();j++)
      {
      if (  data(j,0) != help1 || data(j,1) != help2)
        {
        posend.push_back(j-1);
        if (j < data.rows())
          posbeg.push_back(j);
        xvalues.push_back(data(j,0));
        yvalues.push_back(data(j,1));
        }

      help1 = data(j,0);
      help2 = data(j,1);

      }

    if (posend.size() < posbeg.size())
      posend.push_back(data.rows()-1);


    // 4. initializes ind

    int k;
    workindex = index_data.getV();
    ind = statmatrix<unsigned>(dm.rows(),1);
    for (j=0;j<posend.size();j++)
      {
      for (k=posbeg[j];k<=posend[j];k++,workindex++)
        ind(*workindex,0) = j;
      }

    // TEST
    // ofstream out("c:\\bayesx\\testh\\results\\ind.res");
    // ind.prettyPrint(out);
    // TEST


    // 5. Compute meaneffectnr, mclosest, effectvalues

    double dm_mean1 = dm.mean(0);
    double dm_mean2 = dm.mean(1);
    effectvalues.erase(effectvalues.begin(),effectvalues.end());
    double d1,d2;
    meaneffectnr = 0;
    double distclosest,distcurrent;
    distclosest = pow(data(posbeg[0],0)-dm_mean1,2)+
                  pow(data(posbeg[0],1)-dm_mean2,2);

    for(j=0;j<posbeg.size();j++)
      {
      d1 = data(posbeg[j],0);
      d2 = data(posbeg[j],1);
      distcurrent = pow(d1-dm_mean1,2)+pow(d2-dm_mean2,2);
      if ( distcurrent < distclosest)
        {
        meaneffectnr = j;
        distclosest = distcurrent;
        }

      effectvalues.push_back(ST::doubletostring(d1,0) + "  "
                             + ST::doubletostring(d2,0));

      }

    compute_meaneffectintvar();
    }
  }


void DESIGN_userdefined_tensor::compute_precision(double l)
  {
  if (precisiondeclared==false)
    {
    precision = envmatdouble(K.computeMaxXenv(XWX),1.0,nrpar);
    precisiondeclared = true;
    }


  precision.addto(XWX,Ks[omegaindex],1.0,l);

  //ofstream out("c:\\temp\\K.res");
  //K.print2(out);


  // ofstream out("c:\\temp\\precision.res");
  // precision.print2(out);


  }

DESIGN_userdefined_tensor::DESIGN_userdefined_tensor(datamatrix & dm,datamatrix & iv,
                       datamatrix & designmat1, datamatrix & designmat2,
                       datamatrix & penmat1, datamatrix & penmat2, datamatrix & priormean, datamatrix & constrmat,
                       GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
                       vector<ST::string> & op,
                       vector<ST::string> & vn)
                      : DESIGN_userdefined(o,dp,fcl)
  {
  read_options(op,vn);

  // Adjust to multiple columns in dm
  init_data(dm,iv);

  nrpar = designmat1.cols()*designmat2.cols();

  mK = datamatrix(nrpar,1,0.0);

  if(centermethod==userdefined)
    {
    optionsp->out("\nNOTE: Effect will be centered according to the provided constraint matrix.\n\n");
    basisNull = constrmat;
    }

  unsigned i,j,k,l;
  datamatrix designmat;
  if(designmat2.cols()==1 && designmat2.rows()==1)
    {
    designmat.assign(designmat1);
    }
  else
    {
    designmat= datamatrix(designmat1.rows(), nrpar, 0.0);
    for(i=0; i<designmat1.rows(); i++)
      for(j=0, l=0; j<designmat1.cols(); j++)
         for(k=0; k<designmat2.cols(); k++, l++)
           designmat(i,l) = designmat1(i,j)*designmat2(i,k);
    }

/*  ofstream out("c:\\temp\\Z.raw");
  designmat.prettyPrint(out);
  out.close();*/

  if(nromega>1)
    {
    double rangeomega = 1.0-2.0*minomega;
    for(i=0; i<nromega; i++)
      {
      omegas.push_back(minomega + ((double)i)/((double)(nromega-1)) * rangeomega);
      }
    }
  else
    {
    omegas.push_back(0.5);
    }


  FC_omegas = FC(o,"",1,1,"");
  FC_omegas.setbeta(1,1,omegas[(int)((nromega-1)/2)]);
  omegaindex = (unsigned)((nromega-1)/2);

  datamatrix I1 = datamatrix::diag(penmat1.rows(), 1);
  datamatrix I2 = datamatrix::diag(penmat2.rows(), 1);
  datamatrix K1help = I2.kronecker(penmat1);
  datamatrix K2help = penmat2.kronecker(I1);

  datamatrix Khelp;
  logdets = datamatrix(nromega,1,0.0);

  for(i=0; i<nromega; i++)
    {
    Khelp = omegas[i]*K1help + (1-omegas[i])*K2help;
    Ks.push_back(envmatdouble(Khelp, 0.00000001));
    datamatrix vals(Khelp.rows(),1,0);
    bool eigentest=eigen2(Khelp,vals);
    if(eigentest==false)
      {
      optionsp->out("WARNING: Could not compute determinant\n");
      }
    else
      {
      unsigned j;
      for(j=0; j<vals.rows(); j++)
        {
        if(fabs(vals(j,0))>=0.000001)
          logdets(i,0) += log(vals(j,0));
        }
      }
    }

  compute_Zout(designmat);
  compute_Zout_transposed_vector();

  K = envmatdouble(omegas[omegaindex]*K1help+(1-omegas[omegaindex])*K2help, 0.00000001);
//  K = envmatdouble(K1help+K2help, 0.00000001);

  if(rankK==-1)
    {
    // compute rankK;
    datamatrix Kstat=K1help+K2help;
    datamatrix vals(Kstat.rows(),1,0);
    bool eigentest=eigen2(Kstat,vals);
    if(eigentest==false)
      {
      rankK = Kstat.rows();
      optionsp->out("WARNING: Unable to compute rank of penalty matrix.\n");
      optionsp->out("         Rank was set to the dimension of the penalty matrix: " + ST::inttostring(rankK) + "\n");
      optionsp->out("         Please specify argument rankK\n");
      }
    else
      {
      rankK=0;
      unsigned j;
      for(j=0; j<vals.rows(); j++)
        {
        if(fabs(vals(j,0))>=0.000001)
          rankK++;
        }
      }
    }

  Wsum = datamatrix(posbeg.size(),1,1);

  datamatrix help = designmat.transposed()*designmat;
  XWX = envmatdouble(help, 0.0);
  compute_XtransposedWX();
  XWres = datamatrix(nrpar,1);

  compute_precision(1.0);
  compute_basisNull();
  }


  // COPY CONSTRUCTOR

DESIGN_userdefined_tensor::DESIGN_userdefined_tensor(const DESIGN_userdefined_tensor & m)
    : DESIGN_userdefined(DESIGN_userdefined(m))
  {
  Ks = m.Ks;
  dets = m.dets;
  omegas = m.omegas;
  nromega = m.nromega;
  minomega = m.minomega;
  FC_omegas = m.FC_omegas;
  omegaindex = m.omegaindex;
  xvalues = m.xvalues;
  yvalues = m.yvalues;
  logdets = m.logdets;
  }


  // OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_userdefined_tensor & DESIGN_userdefined_tensor::operator=(const DESIGN_userdefined_tensor & m)
  {
  if (this == &m)
    return *this;
  DESIGN_userdefined::operator=(DESIGN_userdefined(m));
  Ks = m.Ks;
  dets = m.dets;
  omegas = m.omegas;
  nromega = m.nromega;
  minomega = m.minomega;
  FC_omegas = m.FC_omegas;
  omegaindex = m.omegaindex;
  xvalues = m.xvalues;
  yvalues = m.yvalues;
  logdets = m.logdets;
  return *this;
  }

void DESIGN_userdefined_tensor::outbasis_R(ofstream & out)
  {
  }

void DESIGN_userdefined_tensor::outoptions(GENERAL_OPTIONS * op)
  {
  ST::string centerm;

  if (center==false)
    centerm = "uncentered sampling";
  else
    centerm = "centered sampling";

  op->out("  Prior: user defined (tensor product)\n");
  op->out("  Rank of penalty matrix: " + ST::inttostring(rankK) + "\n");
  op->out("  " + centerm + "\n" );

  op->out("\n");

  }

double DESIGN_userdefined_tensor::penalty_compute_quadform(datamatrix & beta)
  {
//  beta.prettyPrint(cout);
  return Ks[omegaindex].compute_quadform(beta,0);
  }

} // end: namespace MCMC




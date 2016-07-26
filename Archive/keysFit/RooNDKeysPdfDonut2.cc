#include "TMath.h"
#include "TMatrixDSymEigen.h"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooArgList.hh" 
#include "RooFitCore/RooAbsReal.hh" 
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooRandom.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooHist.hh"
#include "RooFitCore/RooMsgService.hh"
#include "RooFitCore/RooAbsArg.hh"
#include "RooFitCore/RooAbsPdf.hh"
#include "RooFitCore/RooRealProxy.hh"
#include "RooFitCore/RooSetProxy.hh"
#include "RooFitCore/RooRealConstant.hh"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

class RooRealVar;
class RooArgList;
class RooArgSet;

using namespace std;
using std::cin;

#ifndef __CINT__
class VecVecDouble : public std::vector<std::vector<Double_t> >  { } ;
class VecTVecDouble : public std::vector<TVectorD> { } ;
typedef std::pair<Int_t, VecVecDouble::iterator > iiPair; 
typedef std::vector< iiPair > iiVec; 
typedef std::pair<Int_t, VecTVecDouble::iterator > itPair;
typedef std::vector< itPair > itVec;
#else
class itPair ;
#endif

class RooNDKeysPdfDonut : public RooAbsPdf {

public:


  enum Mirror {NoMirror, MirrorLeft, MirrorRight, MirrorBoth,
               MirrorAsymLeft, MirrorAsymLeftRight,
               MirrorAsymRight, MirrorLeftAsymRight,
               MirrorAsymBoth };

  RooNDKeysPdfDonut(const char *name, const char *title,
               const RooArgList& varList, RooDataSet& data, 
	       TString options="a", Double_t rho=1, Double_t nSigma=3, Bool_t rotate=kTRUE) ; 
  ClassDef(RooNDKeysPdfDonut,0) // General N-dimensional non-parametric kernel estimation p.d.f
};


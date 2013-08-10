
class RooDonutSemilep {
public:
  RooDonutSemilep(double _peak, double _muL, double _muPeak, double _muR, double _nuL, double _nuR,
		  double _biasOffset1, double _biasSlope1, double _biasCurv1,
		  double _biasOffset2, double _biasSlope2, double _biasCurv2,
		  double _sigLow1, double _sigHigh1, double _sigLow2,
		  double _sigHigh2, double _fracLow, double _fracHigh);

  inline virtual ~RooDonutSemilep() { }


  double peak        ;
  double muL         ;
  double muPeak      ;
  double muR         ;
  double nuL         ;
  double nuR         ;
  double biasOffset1 ;
  double biasSlope1  ;
  double biasCurv1   ;
  double biasOffset2 ;
  double biasSlope2  ;
  double biasCurv2   ;
  double sigLow1     ;
  double sigHigh1    ;
  double sigLow2     ;
  double sigHigh2    ;
  double fracLow     ;
  double fracHigh    ;
  
  Double_t evaluate(double pstarl, double mm2)  ;
  Double_t getNorm(Double_t y) const;
  Double_t getMean1(Double_t y) const;
  Double_t getMean2(Double_t y) const;
  Double_t getSigma1(Double_t y) const;
  Double_t getSigma2(Double_t y) const;
  Double_t getFrac(Double_t y) const;
};

class RooDonutCB {
public:
  RooDonutCB(double &_peak, double &_muL, double &_muPeak,
	     double &_muR, double &_nuL, double &_nuR,
	     double &_biasOffset1, double &_biasSlope1, double &_biasCurv1,
	     double &_biasOffset2, double &_biasSlope2, double &_biasCurv2,
	     double &_sigLow1, double &_sigHigh1, double &_sigLow2,
	     double &_sigHigh2, double &_bfrac, double &_bmean,
	     double &_bsig1, double &_bsig2, double &_bpow,
	     double &_bsigpow, double &_fracLow, double &_fracHigh);
  inline virtual ~RooDonutCB() { }



  double peak        ;
  double muL         ;
  double muPeak      ;
  double muR         ;
  double nuL         ;
  double nuR         ;
  double biasOffset1 ;
  double biasSlope1  ;
  double biasCurv1   ;
  double biasOffset2 ;
  double biasSlope2  ;
  double biasCurv2   ;
  double sigLow1     ;
  double sigHigh1    ;
  double sigLow2     ;
  double sigHigh2    ;
  double bfrac       ;
  double bmean       ;
  double bsig1       ;
  double bsig2       ;
  double bpow        ;
  double bsigpow     ;
  double fracLow     ;
  double fracHigh    ;
  
  Double_t evaluate(double pstarl, double mm2)  ;
  Double_t getNorm(Double_t y) const;
  Double_t getMean1(Double_t y) const;
  Double_t getMean2(Double_t y) const;
  Double_t getSigma1(Double_t y) const;
  Double_t getSigma2(Double_t y) const;
  Double_t getFrac(Double_t y) const;
  Double_t getBsigr(Double_t y) const;

};

void plotPl(int sam){                     

  bool isSL = true;
  if(sam==0||sam==1||sam==2||sam==10||sam==11||sam==12||sam==20||
     sam==23||sam==26||sam==29) isSL = false;
  TString fileName="babar_code/fit/ParFit";fileName+=sam;fileName+=".txt";
  fstream parFile;
  parFile.open(fileName,fstream::in);
  // Reading file with parameter range
  Float_t lowpeak, initpeak, highpeak, lowmuL, initmuL, highmuL, lowmuPeak, initmuPeak, highmuPeak;
  Float_t lowmuR, initmuR, highmuR, lownuL, initnuL, highnuL, lownuR, initnuR, highnuR;  
  Float_t lowbc1, initbc1, highbc1, lowbo1, initbo1, highbo1, lowbs1, initbs1, highbs1;
  Float_t lowbc2, initbc2, highbc2, lowbo2, initbo2, highbo2, lowbs2, initbs2, highbs2;
  Float_t lowfl, initfl, highfl, lowfh, initfh, highfh, lowsh1, initsh1, highsh1;
  Float_t lowsh2, initsh2, highsh2, lowsl1, initsl1, highsl1;
  Float_t lowsl2, initsl2, highsl2;
  Float_t lowbfrac, initbfrac, highbfrac, lowbmean, initbmean, highbmean;  
  Float_t lowbsigr, initbsigr, highbsigr, lowbsigl, initbsigl, highbsigl;  
  Float_t lowbpow, initbpow, highbpow, lowbsigpow,initbsigpow, highbsigpow;
  TString pName; int Nfit;
  parFile>>pName>>Nfit;
  if(Nfit!=sam){cout<<"The file "<<fileName<<" is not for sample "<<sam<<endl; return;}
  parFile>>pName>>lowpeak>>initpeak>>highpeak;
  parFile>>pName>>lowmuL>>initmuL>>highmuL;
  parFile>>pName>>lowmuPeak>>initmuPeak>>highmuPeak;
  parFile>>pName>>lowmuR>>initmuR>>highmuR;
  parFile>>pName>>lownuL>>initnuL>>highnuL;
  parFile>>pName>>lownuR>>initnuR>>highnuR;
  parFile>>pName>>lowbo1>>initbo1>>highbo1;
  parFile>>pName>>lowbs1>>initbs1>>highbs1;
  parFile>>pName>>lowbc1>>initbc1>>highbc1;
  parFile>>pName>>lowsl1>>initsl1>>highsl1;
  parFile>>pName>>lowsh1>>initsh1>>highsh1;
  parFile>>pName>>lowbo2>>initbo2>>highbo2;
  parFile>>pName>>lowbs2>>initbs2>>highbs2;
  parFile>>pName>>lowbc2>>initbc2>>highbc2;
  parFile>>pName>>lowsl2>>initsl2>>highsl2;
  parFile>>pName>>lowsh2>>initsh2>>highsh2;
  parFile>>pName>>lowfl>>initfl>>highfl;
  cout<<"The file "<<fileName<<" opened"<<endl;
  parFile>>pName>>lowfh>>initfh>>highfh;
  if(!isSL){
      parFile>>pName>>lowbfrac>>initbfrac>>highbfrac;
      parFile>>pName>>lowbmean>>initbmean>>highbmean;
      parFile>>pName>>lowbsigr>>initbsigr>>highbsigr;
      parFile>>pName>>lowbsigl>>initbsigl>>highbsigl;
      parFile>>pName>>lowbpow>>initbpow>>highbpow;
      parFile>>pName>>lowbsigpow>>initbsigpow>>highbsigpow;
  }


  RooDonutSemilep SL(initpeak, initmuL, initmuPeak,initmuR, initnuL, initnuR,initbo1,initbs1,initbc1,initbo2, 
		     initbs2, initbc2,initsl1,initsh1, initsl2,initsh2, initfl, initfh);
  RooDonutCB CB(initpeak, initmuL, initmuPeak,initmuR, initnuL, initnuR,initbo1,initbs1,initbc1,initbo2, 
		initbs2, initbc2,initsl1,initsh1, initsl2,initsh2, initbfrac, initbmean,
		initbsigr, initbsigl, initbpow,initbpow,initfl, initfh);

  TH1F h("h","p*_{l}",240,0,2.4);
  TCanvas c;
  for(double pl=0.005; pl<2.4; pl+=0.01){
    double val = 0;
    for(double mm=-3; mm<4; mm+=0.04){
      if(isSL) val += SL.evaluate(pl,mm);
      else val += CB.evaluate(pl,mm);
    }
    h.Fill(pl,val);
  }
  h.Draw();
  TString cName = "babar_code/fit/Plots/pl";cName+=sam;cName+=".eps";
  c.SaveAs(cName);
}

Double_t RooDonutSemilep::evaluate(double pstarl, double mm2) 
{
  Double_t norm = getNorm(pstarl);
  Double_t arg1 = mm2 - getMean1(pstarl);
  Double_t arg2 = mm2 - getMean2(pstarl);
  Double_t sigma1 = getSigma1(pstarl);
  Double_t sigma2 = getSigma2(pstarl);
  Double_t frac = getFrac(pstarl);
  if (frac < 0) frac = 0;
  if (frac > 1) frac = 1;
  Double_t eval = norm * (frac*exp(-0.5*arg1*arg1/(sigma1*sigma1)) +
			  (1.0-frac)*exp(-0.5*arg2*arg2/(sigma2*sigma2)));
  return eval;
}

Double_t RooDonutSemilep::getNorm(Double_t y) const
{
  Double_t arg,mu,nu;
  if (y > peak) {
    arg = y - peak;
    mu = muPeak + (muR-muPeak)*(y-peak)/(2.4-peak);
  } else {
    arg = peak - y;
    mu = muL + (muPeak-muL)*(y)/(peak);
  }
  nu = nuL + (nuR-nuL)*(y)/2.4;
  arg /= mu;

  return exp(-1.*pow(arg,nu));
}

Double_t RooDonutSemilep::getMean1(Double_t y) const
{
  return biasOffset1 + biasSlope1*y + biasCurv1*y*y;
}

Double_t RooDonutSemilep::getMean2(Double_t y) const
{
  return biasOffset2 + biasSlope2*y + biasCurv2*y*y;
}

Double_t RooDonutSemilep::getSigma1(Double_t y) const
{
  return sigLow1 + (sigHigh1-sigLow1)*(y/2.4);
}

Double_t RooDonutSemilep::getSigma2(Double_t y) const
{
  return sigLow2 + (sigHigh2-sigLow2)*(y/2.4);
}

Double_t RooDonutSemilep::getFrac(Double_t y) const
{
  return fracLow + (fracHigh-fracLow)*(y/2.4);
}

RooDonutSemilep::RooDonutSemilep(double _peak, double _muL, double _muPeak,
				 double _muR, double _nuL, double _nuR,
				 double _biasOffset1, double _biasSlope1, double _biasCurv1,
				 double _biasOffset2, double _biasSlope2, double _biasCurv2,
				 double _sigLow1, double _sigHigh1, double _sigLow2,
				 double _sigHigh2, double _fracLow, double _fracHigh){
  peak       =_peak        ;
  muL        =_muL         ;
  muPeak     =_muPeak      ;
  muR        =_muR         ;
  nuL        =_nuL         ;
  nuR        =_nuR         ;
  biasOffset1=_biasOffset1 ;
  biasSlope1 =_biasSlope1  ;
  biasCurv1  =_biasCurv1   ;
  biasOffset2=_biasOffset2 ;
  biasSlope2 =_biasSlope2  ;
  biasCurv2  =_biasCurv2   ;
  sigLow1    =_sigLow1     ;
  sigHigh1   =_sigHigh1    ;
  sigLow2    =_sigLow2     ;
  sigHigh2   =_sigHigh2    ;
  fracLow    =_fracLow     ;
  fracHigh   =_fracHigh    ;
}

Double_t RooDonutCB::evaluate(double pstarl, double mm2) 
{
  Double_t norm = getNorm(pstarl);
  Double_t arg1 = mm2 - getMean1(pstarl);
  Double_t arg2 = mm2 - getMean2(pstarl);
  Double_t sigma1 = getSigma1(pstarl);
  Double_t sigma2 = getSigma2(pstarl);
  Double_t frac = getFrac(pstarl);
  if (frac < 0) frac = 0;
  if (frac > 1) frac = 1;
  Double_t gaus  = exp(-0.5*arg1*arg1/(sigma1*sigma1));
  Double_t gaus2 = exp(-0.5*arg2*arg2/(sigma2*sigma2));
  Double_t eval = norm * (frac*gaus + (1.0-frac)*gaus2);

  Double_t bf2 = bfrac * pow((2.4-pstarl)/(2.2),bpow);
  if (bf2 > 1.0) bf2 = 1.0;
  eval *= (1.0-bf2);
  if (mm2>bmean)
    eval += bf2*norm*exp(-0.5*(bmean-mm2)*(bmean-mm2)/(getBsigr(pstarl)*getBsigr(pstarl)));
  else
    eval += bf2*norm*exp(-0.5*(bmean-mm2)*(bmean-mm2)/(bsig2*bsig2));

  return eval;
}

Double_t RooDonutCB::getNorm(Double_t pstarl) const
{
  Double_t arg,mu,nu;
  if (pstarl > peak) {
    arg = pstarl - peak;
    mu = muPeak + (muR-muPeak)*(pstarl-peak)/(2.4-peak);
  } else {
    arg = peak - pstarl;
    mu = muL + (muPeak-muL)*(pstarl)/(peak);
  }
  nu = nuL + (nuR-nuL)*(pstarl)/2.4;
  arg /= mu;

  return exp(-1.*pow(arg,nu));
}

Double_t RooDonutCB::getMean1(Double_t y) const
{
  return biasOffset1 + biasSlope1*y + biasCurv1*y*y;
}

Double_t RooDonutCB::getMean2(Double_t y) const
{
  return biasOffset2 + biasSlope2*y + biasCurv2*y*y;
}

Double_t RooDonutCB::getSigma1(Double_t y) const
{
  return sigLow1 + (sigHigh1-sigLow1)*(y/2.4);
}

Double_t RooDonutCB::getSigma2(Double_t y) const
{
  return sigLow2 + (sigHigh2-sigLow2)*(y/2.4);
}

Double_t RooDonutCB::getFrac(Double_t y) const
{
  return fracLow + (fracHigh-fracLow)*(y/2.4);
}

Double_t RooDonutCB::getBsigr(Double_t y) const
{
  return bsig1 * pow(1.0-y/2.4,bsigpow);
}

RooDonutCB::RooDonutCB(double &_peak, double &_muL, double &_muPeak,
		       double &_muR, double &_nuL, double &_nuR,
		       double &_biasOffset1, double &_biasSlope1, double &_biasCurv1,
		       double &_biasOffset2, double &_biasSlope2, double &_biasCurv2,
		       double &_sigLow1, double &_sigHigh1, double &_sigLow2,
		       double &_sigHigh2, double &_bfrac, double &_bmean,
		       double &_bsig1, double &_bsig2, double &_bpow,
		       double &_bsigpow, double &_fracLow, double &_fracHigh){
   peak        = _peak        ;
   muL         = _muL         ;
   muPeak      = _muPeak      ;
   muR         = _muR         ;
   nuL         = _nuL         ;
   nuR         = _nuR         ;
   biasOffset1 = _biasOffset1 ;
   biasSlope1  = _biasSlope1  ;
   biasCurv1   = _biasCurv1   ;
   biasOffset2 = _biasOffset2 ;
   biasSlope2  = _biasSlope2  ;
   biasCurv2   = _biasCurv2   ;
   sigLow1     = _sigLow1     ;
   sigHigh1    = _sigHigh1    ;
   sigLow2     = _sigLow2     ;
   sigHigh2    = _sigHigh2    ;
   bfrac       = _bfrac       ;
   bmean       = _bmean       ;
   bsig1       = _bsig1       ;
   bsig2       = _bsig2       ;
   bpow        = _bpow        ;
   bsigpow     = _bsigpow     ;
   fracLow     = _fracLow     ;
   fracHigh    = _fracHigh    ; 
}





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

void plotFunction(double _peak, double _muL, double _muPeak, double _muR, double _nuL, double _nuR,
		  double _biasOffset1, double _biasSlope1, double _biasCurv1,
		  double _biasOffset2, double _biasSlope2, double _biasCurv2,
		  double _sigLow1, double _sigHigh1, double _sigLow2,
		  double _sigHigh2, double _fracLow, double _fracHigh)
{                     // peak,  muL, muPeak, muR, nuL, nuR
  //RooDonutSemilep f(    0.6, 0.01,    0.5, 1.2, 0.1,   4,    6.4, -5, 1, 0, 0, 0, 0.01, 1.9, 1, 1, 1, 1);
  //plotFunction(6.1646e-01,2.9984e-02,4.3997e-01,1.2,1,6.9261,6.3045,-4.98,9.7955e-01,0,0,0,1.9508,0,1,1,1,1)
  RooDonutSemilep f(_peak, _muL, _muPeak,_muR, _nuL, _nuR,_biasOffset1, _biasSlope1, _biasCurv1,
		  _biasOffset2, _biasSlope2, _biasCurv2, _sigLow1,_sigHigh1, _sigLow2,_sigHigh2, _fracLow, _fracHigh);
  TH1F h("h","p*_{l}",240,0,2.4);
  TCanvas c;
  for(double pl=0.005; pl<2.4; pl+=0.01){
    double val = 0;
    for(double mm=-3; mm<4; mm+=0.04){
      val += f.evaluate(pl,mm);
    }
    h.Fill(pl,val);
  }
  h.Draw();
  c.SaveAs("pl.eps");
}
void plotCompare(double _peak, double _muL, double _muPeak, double _muR, double _nuL, double _nuR,
		  double _biasOffset1, double _biasSlope1, double _biasCurv1,
		  double _biasOffset2, double _biasSlope2, double _biasCurv2,
		  double _sigLow1, double _sigHigh1, double _sigLow2,
		  double _sigHigh2, double _fracLow, double _fracHigh, double v) {
  //plotCompare(6.1646e-01,2.9984e-02,4.3997e-01,1.2,1,6.9261,6.3045,-4.98,9.7955e-01,1.9508,0,12)  
  RooDonutSemilep f(_peak, _muL, _muPeak,_muR, _nuL, _nuR,_biasOffset1, _biasSlope1, _biasCurv1,
		  0, 0, 0, _sigLow1,_sigHigh1, _sigLow2,_sigHigh2, _fracLow, _fracHigh);
  RooDonutSemilep f2(_peak, _muL, _muPeak, _muR, v, _nuR,_biasOffset1, _biasSlope1, _biasCurv1,
		  0, 0, 0, _sigLow1,_sigHigh1, _sigLow2,_sigHigh2, _fracLow, _fracHigh);
  TH1F h("h","p*_{l}",240,0,2.4);
  TH1F h2("h2","p*_{l}",240,0,2.4);
  h2.SetLineColor(2);
  TCanvas c;
  for(double pl=0.005; pl<2.4; pl+=0.01){
    double val = 0,val2 = 0;
    for(double mm=-3; mm<8; mm+=0.04){
      val += f.evaluate(pl,mm);
      val2 += f2.evaluate(pl,mm);
    }
    h.Fill(pl,val);
    h2.Fill(pl,val2);
  }
  h.Draw();
  h2.Draw("same");
  c.SaveAs("plComp4.eps");
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


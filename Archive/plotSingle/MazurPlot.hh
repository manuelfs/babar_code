#ifndef MAZURPLOT_HH
#define MAZURPLOT_HH

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

#endif

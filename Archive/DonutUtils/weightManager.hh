#ifndef WEIGHTMANAGER_HH
#define WEIGHTMANAGER_HH
#include "Rtypes.h"
#include "TString.h"

class WeightManager {
public:
  WeightManager(const char *weightFile, int IsCoctail=0, Bool_t verb=false);
  ~WeightManager();

  Float_t getEventWeight(Int_t candType, Int_t candDstarType, Int_t MCType,Int_t MCSubmode,Int_t MCDssmode,
			 Int_t MCD,Int_t MCPions, Int_t isBzero, Int_t isSP6, Int_t MCTaumode,
			 Float_t candPstarLep, Float_t trueDmass,
			 Float_t ctl,Float_t ctv,Float_t chi,Float_t q2,Int_t lcharge, Float_t candM2);
  Float_t getCombWeight(Int_t MCCombB,Int_t MCCombDs,Int_t MCDoubleSL,Int_t candTruLep);
  Float_t getUnsimWeight(Int_t mode);
  Float_t getVariWeight(TString Variable, Float_t value, Int_t candType);
  Float_t getBDTBias();

  Bool_t verbose;
  Bool_t doBkgWeights;
  Bool_t doVariWeights;
  Bool_t doPlWeights;
  Bool_t doBFWeights;
  Bool_t doFFWeights;
  Bool_t doDlnFFWeights;
  Bool_t doPi0Weights;
  Bool_t doSigMC;
  Bool_t isZero;
  Int_t sigMode;
  Int_t isCocktail;

  void initializeBkgWeights(const char *folder);
  void initializeVariWeights(const char *folder);
  void initializePlWeights(const char *file);
  void initializeBFWeights(const char *file);
  void initializeFFWeights(const char *file);
  void initializeDlnFFWeights(const char *file);
  void initializePi0Weights(const char *file);
  void initializeSigMC(const char *msg);
  void Print();


  Float_t wD10,wD20,wD00,wD1prime0,wDstarpi,wDstar0pi0,wDpi,wD0pi0;
  Float_t wD1p,wD2p,wD0p,wD1primep,wDstarpi0,wDstar0pi,wDpi0,wD0pi;
  Float_t wDssGlobal;
  Float_t wD0,wDs0,wDp,wDsp;
  Float_t BFIterRatio[8];

  Float_t BkgWeight(Int_t MCType, Int_t candType);
  Float_t wBkg[4], wBkgmES[4];
  Float_t PlWeight(Float_t candPstarLep, Int_t candType);
  Float_t wPlD[4][1000];
  Int_t nbins[8], nVari;
  Float_t variWeights[8][4][40], ranges[8][41], BDTBias;
  TString Variables[8];
  Float_t BFWeight(Int_t candType, Int_t MCType,Int_t MCSubmode,Int_t MCD,Int_t MCPions,Int_t MCTaumode);

  Float_t HQET_R0,HQET_R1,HQET_R2,HQET_rho2, gSR;
  Float_t DssFFWeight(Int_t MCD, Float_t ctl, Float_t q2, Float_t mDss);
  Float_t FFWeight(Int_t MCType,Int_t MCSubmode,Float_t ctl,Float_t ctv,Float_t chi,Float_t q2,Int_t lcharge);

  Float_t ccbarTouds;
  Float_t rhoD2, Delta;
  Float_t DlnFFWeight(Int_t MCType,Float_t ctl,Float_t q2);

  Float_t pi0int,pi0slope;
  Float_t Pi0Weight(Int_t MCType, Int_t candType, Int_t candDstarType, Float_t ppi0, Float_t dssppi0);

  Float_t SigMCWeight(Int_t candType, Int_t MCType,Int_t MCSubmode,Int_t MCDssmode,
		      Int_t isBzero,Int_t isSP6,Int_t MCTaumode);

  Float_t mD0,mDp,mDs0,mDsp,mB;

  Float_t corrDssTau;
  Float_t corrSoftPi;

  Float_t wDoubleSL,wFakeLep,wUnknown;
  Float_t wDsTaunu,wDsEtalnu,wDsEtaprimelnu,wDsphilnu,wDsmunu;
  Float_t wDzbDs,wDzbDsstar,wDstarzbDs,wDstarzbDsstar,wDstarstarDs,wDpDzbKz,wDzbDstarpKz,
    wDstarzbDpKz,wDstarzbDstarpKz,wDzDzbKp,wDstarzbDstarzKp,wDstarzDzbKp,
    wDstarmDpKp,wDzba1p,wDstarzba1p,wDzbRhop,wDstarzbRhop;
  Float_t wDmDs,wDstarmDs,wDmDsstar,wDstarmDsstar,wDmDzKp,wDmDstarzKp,
    wDstarmDstarzKp,wDpDmKz,wDstarmDpKz,wDstarmDstarpKz,wDstarzDzbKz,wDstarpDstarm,
    wDpDstarm,wDpDm,wDma1p,wDstarma1p,wDzbRhoz,wDmRhop,wDstarzbRhoz,wDstarmRhop;
};


#endif

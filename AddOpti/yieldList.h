
#ifndef yieldList_h
#define yieldList_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class yieldList {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           MCType;
   Int_t           MCCombmode;
   Int_t           candType;
   Int_t           candBMode;
   Float_t         candEExtra;
   Float_t         candQ2;
   Float_t         candPMiss;
   Float_t         candMvaComb;
   Float_t         candMvaDl;
   Float_t         candM2;
   Int_t           candBntCha;
   Int_t           candBntNeu;
   Int_t           candDntCha;
   Int_t           candDntNeu;
   Int_t           candLepTru;

   // List of branches
   TBranch        *b_MCType;   //!
   TBranch        *b_MCCombmode;   //!
   TBranch        *b_candType;   //!
   TBranch        *b_candBMode;   //!
   TBranch        *b_candEExtra;   //!
   TBranch        *b_candQ2;   //!
   TBranch        *b_candPMiss;   //!
   TBranch        *b_candBntCha;   //!
   TBranch        *b_candBntNeu;   //!
   TBranch        *b_candDntCha;   //!
   TBranch        *b_candDntNeu;   //!
   TBranch        *b_candLepTru;   //!
   TBranch        *b_candMvaComb;   //!
   TBranch        *b_candMvaDl;   //!
   TBranch        *b_candM2;   //!

   yieldList(TString file="AWG82/ntuples/small/SkimR26_RunAll.root");
   virtual ~yieldList();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(bool doCut = false);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef yieldList_cxx

yieldList::yieldList(TString file) {
  TChain *ch = new TChain("ntp1");
  ch->Add(file);
  Init(ch);
}

yieldList::~yieldList()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t yieldList::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t yieldList::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void yieldList::Init(TChain *tree){
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("MCType", &MCType, &b_MCType);
   fChain->SetBranchAddress("candType", &candType, &b_candType);
   fChain->SetBranchAddress("candBMode", &candBMode, &b_candBMode);
   fChain->SetBranchAddress("candEExtra", &candEExtra, &b_candEExtra);
   fChain->SetBranchAddress("candQ2", &candQ2, &b_candQ2);
   fChain->SetBranchAddress("candPMiss", &candPMiss, &b_candPMiss);
/*    fChain->SetBranchAddress("MCCombmode", &MCCombmode, &b_MCCombmode); */
/*    fChain->SetBranchAddress("candBntCha", &candBntCha, &b_candBntCha); */
/*    fChain->SetBranchAddress("candBntNeu", &candBntNeu, &b_candBntNeu); */
/*    fChain->SetBranchAddress("candDntCha", &candDntCha, &b_candDntCha); */
/*    fChain->SetBranchAddress("candDntNeu", &candDntNeu, &b_candDntNeu); */
   fChain->SetBranchAddress("candLepTru", &candLepTru, &b_candLepTru);
   fChain->SetBranchAddress("candMvaComb", &candMvaComb, &b_candMvaComb);
   fChain->SetBranchAddress("candMvaDl", &candMvaDl, &b_candMvaDl);
   fChain->SetBranchAddress("candM2", &candM2, &b_candM2);
   Notify();
}

Bool_t yieldList::Notify(){
   return kTRUE;
}

void yieldList::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t yieldList::Cut(Long64_t entry){
   return 1;
}
#endif // #ifdef yieldList_cxx

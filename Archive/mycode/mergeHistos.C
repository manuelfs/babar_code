#include "TH2F.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TKey.h"
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
void MergeRootfile( TDirectory *target, TList *sourcelist );

void mergeHisto(){
  Target = TFile::Open("FitAll/MergedPulls.root", "RECREATE" );
  FileList = new TList();
  FileList->Add( TFile::Open("FitAll/PullsRest_Standard_GoodPoisson.root"));
  FileList->Add( TFile::Open("FitAll/PullsRest_Standard_GoodPoisson2.root"));
  FileList->Add( TFile::Open("FitAll/PullsRest_Standard_GoodPoisson3.root"));
  MergeRootfile( Target, FileList );
}

void MergeRootfile( TDirectory *target, TList *sourcelist ) {

  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );
  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd( path );
  TDirectory *current_sourcedir = gDirectory;

  // loop over all keys in this directory
  TChain *globChain = 0;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;
  while ( (key = (TKey*)nextkey())) {

    // read object from first source file
    first_source->cd( path );
    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // descendant of TH1 -> merge it

      cout << "Merging histogram " << obj->GetName() << endl;
      TH1 *h1 = (TH1*)obj;

      // loop over all source files and add the content of the
      // correspondant histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource ) {
        
        // make sure we are at the correct directory level by cd'ing to path
        nextsource->cd( path );
        TH1 *h2 = (TH1*)gDirectory->Get( h1->GetName() );
        if ( h2 ) {
          h1->Add( h2 );
          delete h2; // don't know if this is necessary, i.e. if 
                     // h2 is created by the call to gDirectory above.
        }

        nextsource = (TFile*)sourcelist->After( nextsource );
      }
    }
    else if ( obj->IsA()->InheritsFrom( "TTree" ) ) {
      
      // loop over all source files create a chain of Trees "globChain"
      const char* obj_name= obj->GetName();

      globChain = new TChain(obj_name);
      globChain->Add(first_source->GetName());
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      //      const char* file_name = nextsource->GetName();
      // cout << "file name  " << file_name << endl;
     while ( nextsource ) {
     	  
       globChain->Add(nextsource->GetName());
       nextsource = (TFile*)sourcelist->After( nextsource );
     }

    } else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      // it's a subdirectory

      cout << "Found subdirectory " << obj->GetName() << endl;

      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      MergeRootfile( newdir, sourcelist );

    } else {

      // object is of no type that we know or can handle
      cout << "Unknown object type, name: " 
           << obj->GetName() << " title: " << obj->GetTitle() << endl;
    }

    // now write the merged histogram (which is "in" obj) to the target file
    // note that this will just store obj in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
      target->cd();

      //!!if the object is a tree, it is stored in globChain...
	if(obj->IsA()->InheritsFrom( "TTree" ))
	  globChain->Write( key->GetName() );
	else
	obj->Write( key->GetName() );
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )

  // save modifications to target file
  target->Write();

}

void VaryFraction(int npdf, double fnon, double ftau=0){
  int npdf2, npdf3=-1;
  switch(npdf){
  case 0:
    npdf2 = 4; break;
  case 1:
    npdf2 = 5; break;
  case 2:
    npdf2 = 6; break;
  case 3:
    npdf2 = 7; break;
  case 12:
    npdf2 = 16; npdf3 = 86; break;
  case 13:
    npdf2 = 17; npdf3 = 86; break;
  case 14:
    npdf2 = 18; npdf3 = 86; break;
  case 15:
    npdf2 = 19; npdf3 = 86; break;
  case 20:
    npdf2 = 24; npdf3 = 87; break;
  case 22:
    npdf2 = 25; npdf3 = 87; break;
  case 26:
    npdf2 = 70; break;
  case 27:
    npdf2 = 71; break;
  case 28:
    npdf2 = 72; break;
  case 29:
    npdf2 = 73; break;
  case 30:
    npdf2 = 34; npdf3 = 74; break;
  case 31:
    npdf2 = 35; npdf3 = 75; break;
  case 32:
    npdf2 = 36; npdf3 = 76; break;
  case 33:
    npdf2 = 37; npdf3 = 77; break;
  case 38:
    npdf2 = 78; break;
  case 39:
    npdf2 = 79; break;
  case 40:
    npdf2 = 80; break;
  case 41:
    npdf2 = 81; break;
  case 42:
    npdf2 = 82; break;
  case 43:
    npdf2 = 83; break;
  case 44:
    npdf2 = 84; break;
  case 45:
    npdf2 = 85; break;
  default:
    cout<<"Wrong sample"<<endl;
    return;
  }
  if(ftau<=0 && npdf3 > 0){
    cout<<"Specify ftau"<<endl;
    return;
  }
  TString hname = "keys/root/fitSep/pdfKeys_"; hname += npdf; hname += "_Fit.root";
  TFile hfile(hname); hfile.cd();
  TString h2Name = "h2_"; h2Name += npdf;
  TH2F *h2Tru= (TH2F *)(hfile.Get("h2"))->Clone(h2Name);   

  hname = "keys/root/fitSep/pdfKeys_"; hname += npdf2; hname += "_Fit.root";
  TFile hfile2(hname); hfile2.cd();
  h2Name = "h2_"; h2Name += npdf2;
  TH2F *h2NonTru= (TH2F *)(hfile2.Get("h2"))->Clone(h2Name);   

  double scaleNon = h2Tru->Integral()/h2NonTru->Integral()/(1/fnon-1-ftau/fnon);
  cout<<"scaleNon = "<<scaleNon<<" and h2Tru->Integral() "<<h2Tru->Integral() <<endl;
  if(npdf3>0){
    hname = "keys/root/fitSep/pdfKeys_"; hname += npdf3; hname += "_Fit.root";
    TFile hfile3(hname); hfile3.cd();
    h2Name = "h2_"; h2Name += npdf3;
    TH2F *h2Tau= (TH2F *)(hfile3.Get("h2"))->Clone(h2Name);   
    double scaleTau = h2Tru->Integral()/h2Tau->Integral()/(1/ftau-1-fnon/ftau);
    h2Tru->Add(h2Tau, scaleTau);
    cout<<"scaleTau = "<<scaleTau<<" and h2Tru->Integral() "<<h2Tru->Integral() <<endl;
  }
  h2Tru->Add(h2NonTru, scaleNon);
  TString hnameFit = "keys/root/Fit/pdfKeys_"; hnameFit += npdf; hnameFit += "_Fit.root";
  TFile hfileFit(hnameFit,"RECREATE"); 
  hfileFit.cd();
  h2Tru->Write("h2");
  hfileFit.Close(); 

}


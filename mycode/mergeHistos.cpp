/*

  This macro will add histograms from a list of root files and write them
  to a target root file. The target file is newly created and must not be
  identical to one of the source files.

  Author: Sven A. Schmidt, sven.schmidt@cern.ch
  Date:   13.2.2001
  This code is based on the hadd.C example by Rene Brun and Dirk Geppert,
  which had a problem with directories more than one level deep.
  (see macro hadd_old.C for this previous implementation).
  
  The macro from Sven has been enhanced by 
     Anne-Sylvie Nicollerat <Anne-Sylvie.Nicollerat@cern.ch>
   to automatically add Trees (via a chain of trees).
  
  To use this macro, modify the file names in function hadd.
  
  NB: This macro is provided as a tutorial.
      Use $ROOTSYS/bin/hadd to merge many histogram files

This code is for P07 with global weight //mh

 */


#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "iostream.h"

TList *FileList;
TFile *Target;

void MergeRootfile( TDirectory *target, TList *sourcelist );
Int_t Add(const char* name, TList *sourcelist);

void hadd(const char* name) {
   
  cout << "opening output file" << endl;
  TString fname = name; fname += "/Merged.root";
  Target = TFile::Open( fname, "RECREATE" );
  
  cout << "creating file list " << endl;
  FileList = new TList();

  void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(name));
  int nfiles = 0;
  TString subdir;
  if (dir) {
    int ndir = 0;
     while ((subdir = gSystem->GetDirEntry(dir)) && ndir < 1) {
       if (subdir.Contains(".")) continue;
       if (subdir=="") break;
       subdir = name + subdir; subdir += "/*";  
       cout<<"manuelf: Opening directory "<<subdir<<endl;
       nfiles += Add(subdir, FileList);
       ndir++;
     }
  } else cout<<"manuelf: The directory does not exist"<<endl;

  cout << nfiles<<" files. Starting merge" << endl;
  MergeRootfile( Target, FileList );
  cout << "returning from merge" << endl;
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
      if(!obj->IsA()->InheritsFrom( "TTree" ))
	obj->Write( key->GetName() );
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )

  // save modifications to target file
  target->Write();

}

Int_t Add(const char* name, TList *sourcelist)
{
   // case with one single file
   if (!TString(name).MaybeWildcard()) {
     FileList->Add( TFile::Open(name));
     return 1;
   }

   // wildcarding used in name
   Int_t nf = 0;
   TString basename(name);

   Int_t slashpos = basename.Last('/');
   TString directory;
   if (slashpos>=0) {
      directory = basename(0,slashpos); // Copy the directory name
      basename.Remove(0,slashpos+1);      // and remove it from basename
   } else {
      directory = gSystem->UnixPathName(gSystem->WorkingDirectory());
   }

   const char *file;
   void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(directory.Data()));

   if (dir) {
     //create a TList to store the file names (not yet sorted)
     TList l;
     TRegexp re(basename,kTRUE);
     while ((file = gSystem->GetDirEntry(dir))) {
       if (!strcmp(file,".") || !strcmp(file,"..")) continue;
       TString s = file;
       if ( (basename!=file) && s.Index(re) == kNPOS) continue;
       l.Add(new TObjString(file));
     }
     gSystem->FreeDirectory(dir);
     //sort the files in alphanumeric order
     l.Sort();
     TIter next(&l);
     TObjString *obj;
     while ((obj = (TObjString*)next())) {
       file = obj->GetName();
       TString fopening =  Form("%s/%s",directory.Data(),file);
       FileList->Add( TFile::Open(fopening));
       cout<<"Opening "<<fopening<<endl;
       nf++;
     }
     l.Delete();
   }

   return nf;
}






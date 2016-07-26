// $Id: PidCorrectMesMean.cc,v 1.1 2009/09/03 01:31:56 manuelf Exp $
//
// The file mES_means_Runs12345.txt was obtained by splitting the R18
// Run 1-5 Breco sample (1000 events per ~1.3 fb^-1) into half-overlapping
// subsamples of 2500 events each and fitted the mES distribution in each of the
// subsamples with a sum of an Argus function and a single Gaussian. See
//
//   http://babar-hn.slac.stanford.edu:5090/HyperNews/get/recoTracking/1690.html
//
// The class PidCorrectMesMean, which is compatible both with ROOT and with  
// the BaBar software Framework, provides a correction that sets the Breco mES mean
// to a constant value, 5.2794 GeV/c^2 (PDF 2004). This correction does not change 
// the global average over Runs 1-4, which already have an mES correction applied, see
//
//   http://babar-hn.slac.stanford.edu:5090/HyperNews/get/condmgmt/81/2/1.html
//
// The selection of an appropriate correction is triggered by providing the run number. 
// If you wish use upperID or majorID instead, see 
//
//   int PidDchSvtDrcCalib/PidIdToRunNumber::runNumberfromUpperID(int upperID),
//   int PidDchSvtDrcCalib/PidIdToRunNumber::runNumberfromMajorID(unsigned int majorID).
//_________________________________________________________________________________________
//
// Author list: 
//   Alexandre (Sasha) Telnov, avtelnov@slac.stanford.edu
//
// (C) Princeton University, 2006 
//----------------------------------------------------------------------------------------
// From ROOT, do the following to make it find BaBar/BaBar.hh
// and other BaBar includes:
//
//     gSystem->AddIncludePath("-I../workdir/RELEASE -I../workdir/PARENT");
//
// Then
//
//     gROOT->ProcessLine(".L PidCorrectMesMean.cc++");
//     PidCorrectMesMean::initialize("mES_means_Runs12345.txt"); // R18 ("vivace"), or
//     PidCorrectMesMean::initialize("mES_means_Runs123456-R22c.txt"); // R22, summer 2007 (R22c skims)

#include "BaBar/BaBar.hh"

#include <assert.h>
#include <iostream>
#include <stdio.h>      // sprintf()
#include <fstream>
#include <sstream>
#include <string>

#include <vector>
using std::vector;
#include <map>

#include "DonutUtils/PidCorrectMesMean.hh"

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ifstream;
using std::istringstream; // istrstream is deprecated

using std::string;
using std::map;


//////////////////////////////////////////////////////////////
// Static variables have to be declared here
bool PidCorrectMesMean::_initialized; 
int PidCorrectMesMean::_nLinesRead;
int PidCorrectMesMean::_maxRunNumber;
int PidCorrectMesMean::_minRunNumber;
map<int,float>* PidCorrectMesMean::mapMesMeans; // std::map is not allowed as a static member in CINT, but a pointer to it is
string PidCorrectMesMean::_filename;

float PidCorrectMesMean::_MesReference; // this is what the Breco mES values are normalized to 


/////////////////////////////////////////////////////////////////
// Note that the constructor and the destructor are both private: this class cannot be instantiated
PidCorrectMesMean::PidCorrectMesMean() {}
PidCorrectMesMean::~PidCorrectMesMean() {}


//////////////////////////////////////////////////////////////////
// Load the dE/dx normalization constants from the specified file.
// Values do not have to be ordered, std:map takes care of that.

void PidCorrectMesMean::initialize(const char *filename) {

// First, check if initialize() has already been called; abort if it has been with a different file
  if (_initialized == true && _filename != filename) {
    cerr << "Error: PidCorrectMesMean::initialize() has already been called once, with a different configuration file, \n" << _filename << endl;
    cerr << "Fix your code! Aborting..." << endl;
    assert(false);
  }

  cout << "Loading the run number -> Breco mES mean correspondence table from file " << filename << endl;
  if (_initialized == true && _filename == filename) {
    cerr << "Warning: PidCorrectMesMean::initialize() has already been called once, with the same configuration file. Ignore it." << endl;
    return;
  }

  _initialized=true;
  _filename = filename;

  _MesReference = 5.2794;
  
  _minRunNumber = 200000;
  mapMesMeans = new map<int,float>;
  
// Then, check that the ascii file to load exists
  ifstream ifs(filename);
  if (!ifs.good()) {
    cerr << "Cannot open ASCII table file '" << filename << "', aborting" << endl;
    ifs.close();
    assert(false); // bye-bye
  }
  
// Now, load the constants from this file
  string buffer;
  int runNumber;
  float mESmean;
 
  while (ifs.good()) {

    // check first char (# means it is a comment)
    _nLinesRead++;
    char firstch;
    ifs.get(firstch);
    
    if (ifs.eof()) { // reached end-of-file
      break;
    }
    
    if (firstch == '#') {
      std::getline(ifs,buffer); // swallow the comment line 
      continue;
    }
    
  // OK, this is not a comment line, so put back first char
    ifs.putback(firstch);
  // get a line and turn it into stringstream
    std::getline(ifs, buffer);
    istringstream iss(buffer);
    
  // Parse the stringstream. Its format looks like this:
  // 9931    22-OCT-99       22:56:00
  // where the time and the date are local SLAC, not UTC

    iss >> runNumber;

// The mES correction is not implemented for MC
    assert (runNumber>=9931 && runNumber<=200000); 

    iss >> mESmean;

    (*mapMesMeans)[runNumber] = mESmean;    

    if (runNumber > _maxRunNumber) _maxRunNumber = runNumber;
    if (runNumber < _minRunNumber) _minRunNumber = runNumber;

  } // finished reading config file
  
  _maxRunNumber += 200; // a grace period
  (*mapMesMeans)[_maxRunNumber]=(*(--mapMesMeans->lower_bound(_maxRunNumber-200))).second;
  ++(mapMesMeans->lower_bound(_maxRunNumber-200));

  if (_minRunNumber-9931 > 3000) {
    cout << "Error in PidCorrectMesMean::initialize() : " << _filename 
	 << " does not include corrections for Run 1" << endl; 
    assert (_minRunNumber>12000);
  }

  (*mapMesMeans)[9930]=(*(--mapMesMeans->lower_bound(_minRunNumber+1))).second;
  ++(mapMesMeans->lower_bound(_minRunNumber+1));
  _minRunNumber = 9930; // a grace period

}


//////////////////////////////////////////////

float PidCorrectMesMean::correctMes(float mES, int runNumber) {
  if (_initialized != true) {
    cerr << "PidCorrectMesMean::initialize() must be called first, fix your code!" << endl; 
    assert(_initialized == true);
  }

  if (runNumber>_maxRunNumber || runNumber<_minRunNumber) {
    cerr << "PidCorrectMesMean fatal error: The table in " << _filename << " is valid only from run "
	 << _minRunNumber << " to run " << _maxRunNumber << endl;
    cerr << "You need to update it." << endl;
    assert(runNumber<_maxRunNumber);
  }

  float mESinBreco = (*(--mapMesMeans->lower_bound(runNumber))).second;
  ++(mapMesMeans->lower_bound(runNumber));

// cout << mES << "  " << _MesReference << "  " << mESinBreco << endl;

  return mES + _MesReference - mESinBreco;
}


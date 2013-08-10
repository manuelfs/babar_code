// $Id: PidCorrectMesMean.hh,v 1.1 2009/09/03 01:31:56 manuelf Exp $
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
//     PidCorrectMesMean::initialize("mES_means_Runs12345.txt");

#ifndef PidCorrectMesMean_HH
#define PidCorrectMesMean_HH

#include <vector>
#include <string>
#include <map>

class PidCorrectMesMean {

public:

//  PidCorrectMesMean();
//  ~PidCorrectMesMean();
// 
  static void initialize(const char *filename);
  static float correctMes(float mES, int runNumber);  

private:

// The constructor and the destructor, which should never be used because all members of this class are static
  PidCorrectMesMean();
  ~PidCorrectMesMean();

  static bool _initialized; // has ::initialize() been called?
  static int _nLinesRead; // the number of lines (including comments) read from the config file
  static std::string _filename; // the file used to ::initialize() 
  static int _maxRunNumber; // the table that we have is valid up to this runNumber only
  static int _minRunNumber; // the table that we have is valid from this runNumber only

  static std::map<int,float> *mapMesMeans; // std::map is not allowed as a static member in CINT, but a pointer to it is

  static float _MesReference; // this is what the mES values are Breco normalized to 

};

#endif

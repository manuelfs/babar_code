//----------------------------------------------------------------------------
// keys_utils: utilities for fitting with the keys functions
//
// Description:
//      UTILITIES FOR FITTING WITH THE KEYS
//      ParseName    - Parses the name of a histogram, finding sample number
//                     number of bins and smoothing
//      nameData     - Returns the filename of a pdf sample
//      FormatHisto  - Format for histograms with data points
//      setBins      - Returns plotting range and number of bins
//      WeightedTree - Returns a tree with a branch "weight"
//      RoundNumber  - Converts a number to a string with certain decimals
//      binHisto     - Projects a 2D histogram with certain binning
//      PlotFinalFit - Reads the results of a final fit a plots them
//      PlotEps      - Produces a ps/eps plot from a histogram file
//      IntG         - Integral of a gaussian in a range
//      GaussConv    - Convolution of a histogram and a gaussian
//      getNumberB   - Number of B (and uds,cc) in the files
//      Choleski     - Choleski decomposition
//      ReadFitFile  - Reads a text file with the results of a fit
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      12/04/25 manuelf -- Added ReadFitFile
//      12/03/24 manuelf -- Added Choleski
//      10/10/28 manuelf -- Added getNumberB
//      10/10/15 manuelf -- Added IntG and GaussConv
//      10/03/10 manuelf -- Added PlotFinalFit
//      10/03/02 manuelf -- Added binHisto
//      10/02/12 manuelf -- Added plotting code, PlotEps
//      10/01/18 manuelf -- Created
//----------------------------------------------------------------------------

#ifndef KEYS_UTILS_HH
#define KEYS_UTILS_HH


#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TMatrixT.h"
#include <fstream>
#include <iostream>

#include "styles.hpp"

using namespace TMath;
using namespace std;

double IntG(double mean, double sigma, double minX, double maxX);

double getWentries(TTree *tree, TString Wname="weight");

void ParseName(TString hname, TString &sam, TString &smoo, TString &scaleP, TString &Base, 
	       int &Sam, int &nM2bin, int &nPlbin);

TString nameData(int Sam, TString folder = "");

// Format for histograms with data points
void formatHisto(TH1F *h);

// Plotting range and number of bins
void setBins(int Sam, double &xlow, double &xhigh, int &nbinx, int &nbiny, int IsCocktail=1);

bool isSample(int sam, int MCType, int candType);

TString RoundNumber(double n, int e, double d=1);

void getNumberB(TString genName, TString runs, double &totMCB, double &totdata, 
		double &totuds, double &totccbar, double &totOffdata);

TH1F binHisto(TH2F *h2, int nbins, double xlow, double xhigh, double miny, double maxy, 
	      TString hname, TString var = "candM2");

void FillHisto(TH1F *h, TTree *tree, TString varName, int cand);

void PlotFinalFit(TString textName, TTree *tree, int nbins = 25, double minx=-0.4, double maxx=1.4,
		  TString typePlot="histo", TString var="candM2", TString sample="sig", 
		  double maxFactor = 1., TString psEPS = "EPS");;

TH1F *GaussConv(TH1F *h, float width, TString hname, int PointsBin, int Cnbins, double CminX, double CmaxX,
		double Goffset = 0) ;

void ConvKeys(TH2F *h2, double width);

void PlotEps(TString hname, TString dataFolder, TString epsFolder, TString weightName="1", TString psEPS = "EPS",
	     double maxFactor = 1., int nbins = -1, double minM2 = -99., double maxM2 = -99., 
	     double intFactor = 1., TString typePlot = "curve");

// Cholesky-Banachiewicz algorithm
TMatrixT<double> Choleski(TMatrixT<double> A, int nRows);

// Reads a text file with the results of a fit
void ReadFitFile(TString textName, double Yield[2][70], double Error[70]);
#endif  // KEYS_UTILS_HH

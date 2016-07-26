//------------------------------------------------------------------------
// Description:
//      ChiVector - Vector with chi2 functionality
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      12/10/11 manuelf -- Created 
//------------------------------------------------------------------------

#ifndef CHIVECTOR
#define CHIVECTOR

#include "TMath.h"

#define nMaxIndex 20
using namespace TMath;
using namespace std;

class  ChiVector {
public:
  double v[nMaxIndex][2];
  int nBins;
  ChiVector& operator+(const ChiVector &a);
};

#endif

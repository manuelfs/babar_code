#ifndef PLOT_RESULTS_HH
#define PLOT_RESULTS_HH

#include <vector>
#include <cmath>

#include "TString.h"

class Results {
public: 
  Results(TString iname, float ivalue, std::vector<float> istatError, std::vector<float> isystError={0.});
  TString name;
  float value;
  std::vector<float> statError, systError;

  float systUp()   {return systError[0];}
  float systDown() {return systError[1];}
  float statUp()   {return statError[0];}
  float statDown() {return statError[1];}

  float errUp()    {return sqrt(pow(statUp(),2) + pow(systUp(),2));}
  float errDown()  {return sqrt(pow(statDown(),2) + pow(systDown(),2));}

};


class PadResults {
public:
  PadResults(TString ititle, float iminX, float imaxX, std::vector<Results> ivresults, 
	     Results ism, Results iaverage);
  TString title;
  float minX, maxX;
  std::vector<Results> vresults;
  Results sm, average;
};

#endif  // PLOT_RESULTS_HH

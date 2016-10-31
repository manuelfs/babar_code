#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TEllipse.h"

#include "styles.hpp"
#include "results.hpp"
#include "keys_utils.hpp"

using namespace std;

void getEllipse(float s1, float s2, float corr, float &r1, float &r2, float &angle){
  float sum = pow(s1,2)+pow(s2,2);
  float diff = pow(s1,2)-pow(s2,2);
  float root = pow(s1,4) - 2*pow(s1,2)*pow(s2,2) + 4*pow(s1,2)*pow(corr,2)*pow(s2,2) + pow(s2,4);
  if(root<0) cout<<"Bad root "<<root<<endl;
  root = sqrt(root);

  if(corr==0) angle = 0;
  else angle = -180/3.141593*atan((diff-root)/(2*corr*s1*s2)); // Eigenvector for r1
  r1 = sqrt((sum+root)/2.); // Largest eigenvalue of the cov. matrix
  r2 = sqrt((sum-root)/2.); // Smallest eigenvalue of the cov. matrix
}

int main(){

  //////// Results
  // int cBabar = kGreen, cLHCb = kCyan+1, cBelleHT = kBlue, cBelleST = kOrange;
  // int cAverage = kRed+1, cSM = kGreen+2;

  Results RD_SM ("SM expectation", 0.300, {0.008}, {0}, kBlue-4);
  Results RDs_SM("SM expectation", 0.252, {0.003}, {0}, kViolet+2);

  Results RD_HFAG ("HFAG average", 0.397, {0.040}, {0.028}, kRed, -0.21);
  Results RDs_HFAG("HFAG average", 0.316, {0.016}, {0.010}, kRed, -0.21);

  Results RD_BABAR ("BABAR (HT)", 0.440, {0.058}, {0.042}, kGreen+1, -0.27);
  Results RDs_BABAR("BABAR (HT)", 0.332, {0.024}, {0.018}, kGreen+1, -0.27);

  Results RD_BelleHT ("Belle (HT)", 0.375, {0.064}, {0.026}, kOrange-7, -0.49);
  Results RDs_BelleHT("Belle (HT)", 0.293, {0.038}, {0.015}, kBlue, -0.49);

  Results RDs_BelleST("Belle (ST)", 0.302, {0.030}, {0.011}, kOrange);

  Results RDs_LHCb("LHCb", 0.336, {0.027}, {0.030}, kCyan-3);


  vector<vector<Results> > results;
  results.push_back({RDs_BelleST});
  results.push_back({RDs_LHCb});
  results.push_back({RD_BelleHT,  RDs_BelleHT});
  results.push_back({RD_BABAR,  RDs_BABAR});
  results.push_back({RD_HFAG,  RDs_HFAG});
  results.push_back({RD_SM,    RDs_SM});


  //////// Setting style
  float axisTextSize = 0.065;
  float lMargin = 0.15, rMargin = 0.05;
  float bMargin = 0.17, tMargin = 0.05;
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(bMargin);
  gStyle->SetPadTopMargin(tMargin);
  gStyle->SetPadLeftMargin(lMargin);
  gStyle->SetPadRightMargin(rMargin);
  gStyle->SetTitleSize(axisTextSize,"xy");     // Set the 2 axes title size
  gStyle->SetLabelSize(axisTextSize/1.1,"xy");     // Set the 2 axes label size
  gStyle->SetNdivisions(907, "xy");   
  gStyle->SetTitleOffset(1.17,"x");     
  gStyle->SetTitleOffset(1.17,"y");     
  gStyle->SetTitleFont(132,"xyz");          // Set the all 2 axes title font
  gStyle->SetLabelFont(132,"xyz");          // Set the all 2 axes label font
  gStyle->SetTextFont(132);                // Set global text font
  gStyle->SetPadTickX(1);             // Ticks at the top
  gStyle->SetPadTickY(1);             // Ticks at the right


  //////// Creating canvas and base histogram
  TCanvas can("can","", 500, 380);
  //TCanvas can("can","", 500, 410);
  float minX=0.22, maxX=0.56, minY=0.23, maxY=0.49;
  TH1D histo("histo", "", 10, minX, maxX);
  histo.SetMinimum(minY);
  histo.SetMaximum(maxY);
  histo.GetXaxis()->CenterTitle(true);
  histo.GetYaxis()->CenterTitle(true);
  histo.SetXTitle("R(D)");
  histo.SetYTitle("R(D#lower[-.1]{*})");
  histo.Draw("axis");
  


  //////// Drawing results
  float alpha = 0.04;
  int lWidth = 3;
  TEllipse ellipse; ellipse.SetLineWidth(lWidth);
  TBox box; box.SetLineWidth(0);
  TLine line;  line.SetLineWidth(lWidth+1);
  float maxR, minR, angle;
  for(auto &result : results){
    if(result.size()>1){
      float rd = result[0].value, rds = result[1].value;
      float erd = result[0].errUp(), erds = result[1].errUp();
      getEllipse(erd, erds, result[0].correl, maxR, minR, angle);
      if(result[0].name.Contains("aver") || result[0].name.Contains("SM")) 
	ellipse.SetFillColorAlpha(result[0].color, alpha+0.18);
      else ellipse.SetFillColorAlpha(result[0].color, alpha);
      ellipse.SetLineColor(result[0].color);
      ellipse.DrawEllipse(rd, rds, maxR, minR, 0, 360, angle, "c");
      if(result[0].name.Contains("SM")) ellipse.DrawEllipse(rd, rds, maxR, minR, 0, 360*10, angle, "c");
    } else {
      float rdsUp = result[0].value + result[0].errUp();
      float rdsDown = result[0].value - result[0].errDown();
      box.SetFillColorAlpha(result[0].color, alpha);
      box.DrawBox(minX, rdsDown, maxX, rdsUp);
      line.SetLineColor(result[0].color);
      line.DrawLine(minX, rdsUp, maxX, rdsUp);
      line.DrawLine(minX, rdsDown, maxX, rdsDown);
    }
  } // Loop over results

  //////// Legend
  int Ncols = 2;
  double legX(lMargin+0.05), legY = 1-tMargin-0.04, legSingle = 0.06;
  double legW = 0.28*Ncols, legH = legSingle*(results.size()+Ncols-1)/Ncols;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.045); leg.SetFillColor(0); 
  leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetNColumns(Ncols);
  vector<TH1D*> hLeg;
  for(size_t ind=results.size()-3; ind<=results.size()-3; ind--){
    TString hname = "hLeg"+to_string(ind);
    hLeg.push_back(new TH1D(hname,"",10,0,10));
    hLeg.back()->SetLineWidth(lWidth);
    hLeg.back()->SetLineColor(results[ind][0].color);
    hLeg.back()->SetFillColorAlpha(results[ind][0].color, alpha);
    leg.AddEntry((hLeg.back()), results[ind][0].name, "lf");
  }
  for(size_t ind=results.size()-2; ind<results.size(); ind++){
    TString hname = "hLeg"+to_string(ind);
    hLeg.push_back(new TH1D(hname,"",10,0,10));
    hLeg.back()->SetLineWidth(lWidth);
    hLeg.back()->SetLineColor(results[ind][0].color);
    hLeg.back()->SetFillColorAlpha(results[ind][0].color, alpha);
    leg.AddEntry((hLeg.back()), results[ind][0].name, "lf");
  }
  leg.Draw();


  histo.Draw("same axis");


  TString plotName = "plots/rdx_ellipses.pdf";
  can.SaveAs(plotName);

  return 0;
}


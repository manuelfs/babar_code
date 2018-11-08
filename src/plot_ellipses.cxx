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

  Results RD_SM ("SM pred.", 0.300, {0.008}, {0}, kBlue-4);
  Results RDs_SM("SM pred.", 0.252, {0.003}, {0}, kViolet+2);

  Results RD_HFLAV ("HFLAV aver.", 0.407, {0.039}, {0.024}, kRed, -0.20);
  Results RDs_HFLAV("HFLAV aver.", 0.304, {0.013}, {0.007}, kRed, -0.20);

  Results RD_BABAR ("BABAR (HT)", 0.440, {0.058}, {0.042}, kGreen+1, -0.27);
  Results RDs_BABAR("BABAR (HT)", 0.332, {0.024}, {0.018}, kGreen+1, -0.27);

  Results RD_BelleHT ("Belle (HT)", 0.375, {0.064}, {0.026}, kOrange-7, -0.49);
  Results RDs_BelleHT("Belle (HT)", 0.293, {0.038}, {0.015}, kBlue, -0.49);

  Results RDs_BelleST("Belle (ST)", 0.302, {0.030}, {0.011}, kOrange);

  Results RDs_Bellepi("Belle (#pi/#rho)", 0.270, {0.035}, {0.027}, kMagenta-4);

  Results RDs_LHCb("LHCb (#mu)", 0.336, {0.027}, {0.030}, kCyan-3);
  Results RDs_LHCb2("LHCb (#pi#pi#pi)", 0.285, {0.019}, {0.029}, kMagenta+2);


  vector<vector<Results> > results;
  results.push_back({RD_SM,    RDs_SM});
  results.push_back({RD_BABAR,  RDs_BABAR});
  results.push_back({RDs_LHCb});
  results.push_back({RD_BelleHT,  RDs_BelleHT});
  results.push_back({RDs_BelleST});
  results.push_back({RDs_Bellepi});
  results.push_back({RDs_LHCb2});
  results.push_back({RD_HFLAV,  RDs_HFLAV});


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
  float minX=0.20, maxX=0.54, minY=0.21, maxY=0.49;
  TH1D histo("histo", "", 10, minX, maxX);
  histo.SetMinimum(minY);
  histo.SetMaximum(maxY);
  histo.GetXaxis()->CenterTitle(true);
  histo.GetYaxis()->CenterTitle(true);
  histo.SetXTitle("R(D)");
  histo.SetYTitle("R(D#lower[-.1]{*})");
  histo.Draw("axis");
  

  double Nsingle = 0; // Number of 1D measurements

  //////// Drawing results
  float alpha = 0.04;
  int lWidth = 3;
  TEllipse ellipse; 
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
      ellipse.SetLineColor(result[0].color); ellipse.SetLineWidth(lWidth);
      ellipse.DrawEllipse(rd, rds, maxR, minR, 0, 360, angle, "c");
      if(result[0].name.Contains("SM")) ellipse.DrawEllipse(rd, rds, maxR, minR, 0, 360*10, angle, "c");
      //cout<<"RDs = "<<rds<<" +- "<<erds/rds*100<<endl;
      //cout<<"RD  = "<<rd<<" +- "<<erd/rd*100<<endl<<endl;

    } else {
      float dx = 0.003;
      float rd = 0.23+0.01*Nsingle;
      float rds = result[0].value;
      float rdsUp = result[0].value + result[0].errUp();
      float rdsDown = result[0].value - result[0].errDown();
      //cout<<"RDs = "<<rds<<" +- "<<(rdsUp-rds)/rds*100<<endl<<endl;

      ellipse.SetLineColor(result[0].color); ellipse.SetLineWidth(lWidth+2);
      ellipse.DrawEllipse(rd, rds, 0.0023, 0.0023, 0, 360, 0, "c");
      line.SetLineWidth(lWidth);
      line.SetLineColor(result[0].color);
      line.DrawLine(rd, rdsUp, rd, rdsDown);
      line.DrawLine(rd-dx, rdsUp, rd+dx, rdsUp);
      line.DrawLine(rd-dx, rdsDown, rd+dx, rdsDown);
      Nsingle++;
    }
  } // Loop over results

  //////// Legend
  int Ncols = 3;
  double legX(lMargin+0.05), legY = 1-tMargin-0.04, legSingle = 0.06;
  double legW = 0.25*Ncols, legH = legSingle*(8+Ncols-1)/Ncols;
  if(results.size()<=2) legW = 0.217*Ncols;
  if(results.size()<=3) {
    legH = legSingle*(8+Ncols-1)/Ncols/3;   
  } else if(results.size()<=6) legH = legSingle*(8+Ncols-1)/Ncols/3*2;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.045); leg.SetFillColor(0); 
  leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetNColumns(Ncols);
  vector<TH1D*> hLeg;
  // for(size_t ind=results.size()-3; ind<=results.size()-3; ind--){
  //   TString hname = "hLeg"+to_string(ind);
  //   hLeg.push_back(new TH1D(hname,"",10,0,10));
  //   hLeg.back()->SetLineWidth(lWidth);
  //   hLeg.back()->SetLineColor(results[ind][0].color);
  //   hLeg.back()->SetFillColorAlpha(results[ind][0].color, alpha);
  //   leg.AddEntry((hLeg.back()), results[ind][0].name, "lf");
  // }
  for(size_t ind=0; ind<results.size(); ind++){
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


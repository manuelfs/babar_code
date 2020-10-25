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

  if(corr==0) angle = 90;
  else angle = -180/3.141593*atan((diff-root)/(2*corr*s1*s2)); // Eigenvector for r1
  r1 = sqrt((sum+root)/2.); // Largest eigenvalue of the cov. matrix
  r2 = sqrt((sum-root)/2.); // Smallest eigenvalue of the cov. matrix
}

int main(){

  bool isWide = true;

  //////// Results
  // int cBabar = kGreen, cLHCb = kCyan+1, cBelleHT = kBlue, cBelleST = kOrange;
  // int cAverage = kRed+1, cSM = kGreen+2;

  // Numbers from HFLAV Spring 2019
  // https://hflav-eos.web.cern.ch/hflav-eos/semi/spring19/html/RDsDsstar/RDRDs.html
  Results RD_SM ("SM pred.", 0.299, {0.003}, {0}, kViolet-2);
  Results RDs_SM("SM pred.", 0.258, {0.005}, {0}, kViolet+2);

  Results RD_HFLAV ("HFLAV aver.", 0.340, {0.027}, {0.013}, kRed, -0.38);
  Results RDs_HFLAV("HFLAV aver.", 0.295, {0.011}, {0.008}, kRed, -0.38);

  Results RD_BABAR ("BABAR (HT)", 0.440, {0.058}, {0.042}, kGreen+2, -0.31);
  Results RDs_BABAR("BABAR (HT)", 0.332, {0.024}, {0.018}, kGreen+2, -0.31);

  Results RD_BABARe ("BABAR electron", 0.3466, {0.0823}, {0.04}, kGreen+1, -0.41);
  Results RDs_BABARe("BABAR HT electron", 0.3671, {0.0342}, {0.02}, kGreen+1, -0.41);

  Results RD_BABARmu ("BABAR muon", 0.5215, {0.0796}, {0.04}, kGreen-3, -0.41);
  Results RDs_BABARmu("BABAR HT muon", 0.2899, {0.0326}, {0.02}, kGreen-3, -0.41);

  Results RD_BelleST("Belle (ST)", 0.307, {0.037}, {0.016}, kBlue+3, -0.52);
  Results RDs_BelleST("Belle (ST)", 0.283, {0.018}, {0.014}, kBlue+3, -0.52);

  Results RD_BelleSTe("Belle electron", 0.281, {0.042}, {0.017}, kAzure+1, -0.62);
  Results RDs_BelleSTe("Belle ST electron", 0.304, {0.022}, {0.016}, kAzure+1, -0.62);

  Results RD_BelleSTmu("Belle muon", 0.373, {0.068}, {0.030}, kBlue-4, -0.42);
  Results RDs_BelleSTmu("Belle ST muon", 0.245, {0.035}, {0.020}, kBlue-4, -0.42);


  if(isWide){
    RD_BABAR.name = "BABAR, PRL #font[62]{109}, 101802 (2012)";
    //RDs_BelleST.name = "Belle, PRD #font[62]{94}, 072007 (2016)";
    RD_BelleST.name = "Belle, PRL #font[62]{124}, 161803 (2020)";
    RDs_BelleST.name = "Belle, PRL #font[62]{124}, 161803 (2020)";
    RD_HFLAV.name = "HFLAV average Spring 2019";
    RD_SM.name = "SM predictions";
  }

  vector<vector<Results> > results;
  results.push_back({RD_BABAR,  RDs_BABAR});
  results.push_back({RD_BelleST, RDs_BelleST});
  results.push_back({RD_BABARe,  RDs_BABARe});
  results.push_back({RD_BelleSTe, RDs_BelleSTe});
  results.push_back({RD_BABARmu,  RDs_BABARmu});
  results.push_back({RD_BelleSTmu, RDs_BelleSTmu});
  results.push_back({RD_HFLAV,  RDs_HFLAV});
  results.push_back({RD_SM,    RDs_SM});


  //////// Setting style
  float axisTextSize = 0.065;
  float lMargin = 0.1, rMargin = 0.02;
  float bMargin = 0.17, tMargin = 0.05;
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(bMargin);
  gStyle->SetPadTopMargin(tMargin);
  gStyle->SetPadLeftMargin(lMargin);
  gStyle->SetPadRightMargin(rMargin);
  gStyle->SetNdivisions(907, "xy");   
  gStyle->SetTitleOffset(1.17,"x");
  gStyle->SetTitleOffset(1.17,"y");     
  gStyle->SetTitleFont(132,"xyz");          // Set the all 2 axes title font
  gStyle->SetLabelFont(132,"xyz");          // Set the all 2 axes label font
  gStyle->SetTextFont(132);                // Set global text font
  gStyle->SetPadTickX(1);             // Ticks at the top
  gStyle->SetPadTickY(1);             // Ticks at the right

  if(isWide) {
    gStyle->SetTitleOffset(0.6,"y");
    gStyle->SetTitleOffset(0.98,"x");
    axisTextSize *= 1.3;
  }
  gStyle->SetTitleSize(axisTextSize,"xy");     // Set the 2 axes title size
  gStyle->SetLabelSize(axisTextSize/1.1,"xy");     // Set the 2 axes label size


  //////// Creating canvas and base histogram
  int cW = 500, cH = 380;
  if(isWide) cW = 780;
  TCanvas can("can","", cW, cH);
  //TCanvas can("can","", 500, 410);
  float minX=0.20, maxX=0.54, minY=0.21, maxY=0.49;
  if(isWide){
    minX = 0.21; maxX = 0.63; minY = 0.17; maxY = 0.62;
  }
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
  int lWidth = 2;
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
	ellipse.SetFillColorAlpha(result[0].color, alpha+0.32);
      else ellipse.SetFillColorAlpha(result[0].color, alpha);
      ellipse.SetLineColor(result[0].color); ellipse.SetLineWidth(lWidth);
      if(result[0].name.Contains("electron")) ellipse.SetLineStyle(2);
      else if(result[0].name.Contains("muon")) ellipse.SetLineStyle(3);
      else ellipse.SetLineStyle(1);
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
      // // Drawing bands for RD* results only
      // box.SetFillColorAlpha(result[0].color, alpha*2);
      // box.DrawBox(minX, rdsDown, maxX, rdsUp);
      // line.SetLineWidth(1);
      // line.DrawLine(minX, rdsUp, maxX, rdsUp);
      // line.DrawLine(minX, rdsDown, maxX, rdsDown);
      Nsingle++;
    }
  } // Loop over results

  if(isWide){
    Results RD_HFLAV4 ("noplot", RD_HFLAV.value, {static_cast<float>(3.4 * RD_HFLAV.errUp())}, {0}, kRed, RD_HFLAV.correl);
    Results RDs_HFLAV4("noplot", RDs_HFLAV.value, {static_cast<float>(3.4* RDs_HFLAV.errUp())}, {0}, kRed, RD_HFLAV.correl);
    vector<Results>  result({RD_HFLAV4,  RDs_HFLAV4});

    // float rd = result[0].value, rds = result[1].value;
    // float erd = result[0].errUp(), erds = result[1].errUp();
    // getEllipse(erd, erds, result[0].correl, maxR, minR, angle);
    // ellipse.SetFillStyle(0);
    // ellipse.SetLineWidth(1);
    // ellipse.SetLineStyle(2);
    // ellipse.SetLineColor(result[0].color); 
    // ellipse.DrawEllipse(rd, rds, maxR, minR, 0, 360, angle, "c");
    // TLatex label;  //label.SetNDC(kTRUE);
    // label.SetTextAlign(23); label.SetTextSize(0.06);
    // label.SetTextColor(kRed+1);
    // label.DrawLatex(0.3, 0.365, "#font[62]{3#sigma}");

  }

  //////// Legend
  int Ncols = 3; if(isWide) Ncols = 2;
  double legX(lMargin+0.05), legY = 1-tMargin-0.04, legSingle = 0.06, leg2X = legX;
  double legW = 0.25*Ncols, legH = legSingle*(8+Ncols-1)/Ncols, legFontSize = 0.045;
  if(results.size()<=2) legW = 0.217*Ncols;
  if(results.size()<=3) {
    legH = legSingle*(8+Ncols-1)/Ncols/3;   
  } else if(results.size()<=6) legH = legSingle*(8+Ncols-1)/Ncols/3*2;
  if(isWide) {
    Ncols = 1;
    legX -= 0.02;
    leg2X = legX+0.42;
    legW = 0.15;
    legH = legSingle * results.size()/1.8;
    legFontSize = 0.057;
  } 
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(legFontSize); leg.SetFillColor(0); 
  leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetNColumns(Ncols);
  TLegend leg2(leg);
  leg2.SetX1(leg2X); leg2.SetX2(leg2X+legW); 
  vector<TH1D*> hLeg;
  if(!isWide){
    for(size_t ind=0; ind<results.size(); ind++){
      TString hname = "hLeg"+to_string(ind);
      hLeg.push_back(new TH1D(hname,"",10,0,10));
      hLeg.back()->SetLineWidth(lWidth);
      hLeg.back()->SetLineColor(results[ind][0].color);
      hLeg.back()->SetFillColorAlpha(results[ind][0].color, alpha);
      leg.AddEntry((hLeg.back()), results[ind][0].name, "lf");
    }
    leg.Draw();
  }else {
    for(size_t ind=0; ind<results.size(); ind++){
      TString hname = "hLeg"+to_string(ind);
      hLeg.push_back(new TH1D(hname,"",10,0,10));
      hLeg.back()->SetLineWidth(lWidth);
      hLeg.back()->SetLineColor(results[ind][0].color);
      hLeg.back()->SetFillColorAlpha(results[ind][0].color, alpha);
      if(results[ind][0].name.Contains("electron")) hLeg.back()->SetLineStyle(2);
      else if(results[ind][0].name.Contains("muon")) hLeg.back()->SetLineStyle(3);
      else hLeg.back()->SetLineStyle(1);
      if(ind%2==0) leg.AddEntry((hLeg.back()), results[ind][0].name, "lf");
      else leg2.AddEntry((hLeg.back()), results[ind][0].name, "lf");
    }
    leg.Draw(); leg2.Draw();
  }

  histo.Draw("same axis");


  TString plotName = "plots/rdx_emu.pdf";
  can.SaveAs(plotName);

  return 0;
}


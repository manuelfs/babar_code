#include "mycode/cuts.hh"
TString round(double n, double d);

void CSample(int nvari=0, TString type = "pl", TString run = "All", TString cuts="Mva"){
  TString isMaz = "";
  if(nvari<0||nvari>8){
    cout<<"Variable must be between 0 and 8. Exiting"<<endl;
    return;
  }

  if(cuts=="Mva") basic += Mva;
  else if (nvari!=3) basic += ee;
  if(cuts.Contains("cosT"))
    basic += cosT;

  if(cuts.Contains("All"))
    dss += mpi0+dssee;
  else if(cuts.Contains("Mva"))
    dss += dssMva;
  if(cuts.Contains("cosT"))
    dss += cosT;

  TString dfile = "AWG82/ntuples/small/Adddata30_Run"; dfile += run; dfile +=".root";
  TString mfile = "AWG82/ntuples/small/Add_Truth30Rest_Run"; mfile += run; mfile +=".root";
  TChain d("ntp1");
  d.Add(dfile);
  double fracRest = 1.;
  if(mfile.Contains("Rest")) fracRest = 2/3.;
  TChain m("ntp1");
  m.SetLineColor(2);
  m.Add(mfile);
  cout<<"Data: "<<dfile<<" with entries: "<<d.GetEntries()<<endl;   
  cout<<"MC: "<<mfile<<" with entries: "<<m.GetEntries()<<endl;   
//   if(isMaz=="") {
//     TString directory = "AWG82/ntuples/small/Adddata30_Run"; 
//     //TString directory = "AWG82/ntuples/Adddata30/Merged/*Run"; 
//     //TString directory = "AWG82/ntuples/Adddata/Merged/*Run"; 
//     //TString directory = "AWG82/ntuples/dataNN/Merged/*Run"; 
//     for(int i=1; i<7; i++){
//       TString irun = ""; irun += i;
//       if(run.Contains(irun)){
// 	TString filename=directory; filename += irun;filename +="*";
// 	d.Add(filename);
// 	cout<<"Added "<<filename<<" with entries: "<<d.GetEntries()<<endl;   
//       }
//     }
//   }else{
//     if(run=="All")
//       d.Add("mdata/data*");
//     else
//       d.Add("mdata/data_run"+run+"*");
//   }
//   if(isMaz==""){
//     //TString directory = "AWG82/ntuples/small/Add_Truth30_Run"; 
//     //TString directory = "AWG82/ntuples/Add_Final/MergedRest/*Run"; 
//     //TString directory = "AWG82/ntuples/Add_Truth30/Merged/*Run"; 
//     //TString directory = "AWG82/ntuples/Add_BDT_KM/MergedRest/*Run"; 
//     //TString directory = "AWG82/ntuples/NN/Merged/*Run"; 
//     for(int i=1; i<7; i++){
//       TString irun = ""; irun += i;
//       if(run.Contains(irun)){
// 	TString filename=directory; filename += irun;filename +="*";
// 	m.Add(filename);
// 	cout<<"Added "<<filename<<" with entries: "<<m.GetEntries()<<endl;   
//       }
//     }
//     if(directory.Contains("Rest")) fracRest = 2/3.;
//     if(run=="Mva")
//       m.Add("AWG82/ntuples/small/Add_Truth30_Run5.root");      
//   }else{
//     if(run=="All")
//       m.Add("mdata/SP*");
//     else
//       m.Add("mdata/SP*"+run+"_maz.root");
//   }
  m.SetLineColor(2);
  TCut tot;
  int nlim = 0;
  if(type == "pl")
    tot = basic+"candPstarLep>1.5";
  else if(type == "q2"){
    tot = basic+"candQ2<5";
    nlim = 9;
  } else if(type == "dss"){
    tot = dss;
    nlim = 18;
  }else
    tot = basic+"candPstarLep>1.5&&candQ2>4";    

  TString vari[] = {"candM2","candPstarLep","candQ2","candEExtra","candMES","candMvaDl",
		    "candRejectedTracks","candRejectedPhotons","abs(candCosT)"};
  if(type=="dss"){
    vari[0] = "mm2pi0";
    vari[3] = "eextrapi0";
  }
  TString Namevari[] = {"M2","PstarLep","Q2","EExtra","MES","MvaDl",
		    "RejectedTracks","RejectedPhotons","CosT"};
  if(nvari==3) vari[3] = "eextrapi0";
  TString limits[] = {"(15,-.5,1.5)","(15,1.5,2.5)","(15,-2,14)","(25,0,2)","(15,5.27,5.29)","(15,-0.8,0.8)","(5,0,5)","(12,0,12)","(15,0,1)",
		      "(15,-.5,1.5)","(15,0,2.5)","(15,-1,5)","(25,0,.5)","(15,5.27,5.29)","(15,-0.8,0.8)","(5,0,5)","(12,0,12)","(15,0,1)",
		      "(20,-2.5,5)","(15,0,2.4)","(15,-1,5)","(15,0,.5)","(15,5.27,5.29)","(15,-0.,0.8)","(5,0,5)","(12,0,12)","(15,0,1)"};
  TString titles[] = {"D^{0}: MC/data = ","D*^{0}: MC/data = ","D^{+}: MC/data = ","D*^{+}: MC/data = "};
  TString xtitles[] = {"m^{2}_{miss} [GeV^{2}]","p*_{l} [GeV]","q^{2} [GeV^{2}]",
		       "E_{Extra} [GeV]","m_{ES} [GeV]","Dlnu MVA","Rejected tracks", "Rejected photons", "abs[cos(#theta_{T})]"};
  double MCBs[] = {36968000+37200000., 103498000+103124000., 50556000+49766000., 167994000+167466000., 
		   244322000+244812000., 68148000+68016000.};
  double dataBs[] = {22389980.4, 67394307.5, 35569248.8, 110449802.7, 147190396.5, 84767412.6};
  //double ratios[] = {3.31255, 3.06587, 2.82047, 3.03722, 3.32314, 1.60632};
  double nMCB = 0, ndataBs = 0;
  for(int i=1; i<7; i++){
    TString irun = ""; irun += i;
    if(run.Contains(irun)){
      nMCB += MCBs[i-1]; ndataBs += dataBs[i-1];
    }
  }
  double MCDlumi = 1;
  if(ndataBs) MCDlumi = nMCB/ndataBs;

  TLatex *label = new TLatex();
  label->SetNDC(kTRUE);
  TH1F *data[4];
  TH1F *mc[4];
  TH1F *data2[4];
  TH1F *mc2[4];
  TH1F *mcDs[4];
  TCanvas c("dataMC","data Vs MC",800,600);
  c.Divide(2,2);

  for(int i=0; i<4; i++){
    c.cd(i+1);
    TPad *p1 = new TPad("p1","",0,.25,1,1);
    p1->Draw();
    p1->cd();
    TString cand = "candType=="; cand+=i+1;
    TCut totcand = tot; totcand += cand;
    d.Draw(vari[nvari]+">>dat"+limits[nvari+nlim],totcand);
    TH1F *h = (TH1F*)gDirectory->Get("dat");
    data[i] = (TH1F*)h->Clone("data"+i);
    data[i]->Sumw2();
    data[i]->SetMarkerStyle(20);
    data[i]->SetMarkerSize(.6);
    data[i]->GetXaxis()->SetTitleOffset(1);
    data[i]->GetXaxis()->SetTitleSize(0.08);
    data[i]->SetXTitle("");
    gStyle->SetOptStat(0);

    TCut Dscut = "";
    m.Draw(vari[nvari]+">>MC"+limits[nvari+nlim],totcand+Dscut);
    h = (TH1F*)gDirectory->Get("MC");
    mc[i] = (TH1F*)h->Clone("MonteCarlo"+i);
    double ndata = (double)data[i]->GetEntries(); double nmc = (double)mc[i]->GetEntries();
    data[i]->SetTitle("");
    cout<<titles[i]<<": "<<ndata<<" Vs "<<nmc<<endl;
    if(nmc)
      mc[i]->Scale(ndata/nmc);

    float maxi = data[i]->GetMaximum();
    if(mc[i]->GetMaximum()>maxi) maxi = mc[i]->GetMaximum();
    data[i]->SetMaximum(1.15*maxi);
    data[i]->SetMinimum(0);
    data[i]->GetYaxis()->SetLabelSize(0.062);
    data[i]->GetXaxis()->SetLabelSize(0.08);
    data[i]->GetYaxis()->SetNdivisions(6+100*2);
    data[i]->GetXaxis()->SetNdivisions(7+100*2);
    data[i]->Draw();
    mc[i]->Draw("same");
    TString RatioTitle = titles[i]; RatioTitle += round(nmc/MCDlumi/fracRest,ndata);
    label->SetTextSize(0.085);
    label->DrawLatex(0.02,0.92,RatioTitle);
    TLegend *leg;
    if(nvari==-5)
      leg = new TLegend(0.1,.7,0.47,0.9);
    if(nvari==8){
      if(type=="dss")
	leg = new TLegend(0.1,.7,0.5,0.9);
      else
	leg = new TLegend(0.5,.1,0.9,0.3);
    }
    else
      leg = new TLegend(0.53,.7,0.9,0.9);
    leg->SetTextSize(0.068);
    //leg->SetBorderSize(0);
    leg->SetFillColor(0);
    TString dleg = "data ("; dleg+=(int)ndata;dleg+=" entries)";
    leg->AddEntry(data[i],dleg);
    TString mleg = "MC ("; mleg+=(int)nmc;mleg+=" entries)";
    leg->AddEntry(mc[i],mleg);
    leg->Draw();
    c.cd(i+1);
    TPad *p2 = new TPad("p2","",0,0,1,.2);
    p2->Draw();
    p2->cd();
    
    Float_t min = data[i]->GetXaxis()->GetBinLowEdge(data[i]->GetXaxis()->GetFirst());
    Float_t max = data[i]->GetXaxis()->GetBinLowEdge(data[i]->GetXaxis()->GetLast()+1);
    data2[i] = (TH1F *)data[i]->Clone();
    mc2[i] = (TH1F *)mc[i]->Clone();
    mc2[i]->Sumw2();
    data2[i]->Divide(mc2[i]);
    TLine l;
    l.SetLineColor(2);
    data2[i]->GetYaxis()->SetNdivisions(3+100*2);
    data2[i]->SetTitle("");
    data2[i]->SetXTitle("");
    data2[i]->SetMaximum(1.5);
    data2[i]->SetMinimum(0.5);
    data2[i]->GetYaxis()->SetLabelSize(0.25);
    data2[i]->GetXaxis()->SetLabelSize(0.01);
    data2[i]->SetMarkerSize(.5);
    data2[i]->Draw("E0 P");
    l.DrawLine(min, 1.0, max, 1.0);
    c.cd(i+1);
    label->SetTextSize(0.06);
    label->DrawLatex(0.68,0.21,xtitles[nvari]);
    
  }
  c.SaveAs("mycode/PlotsCS/"+isMaz+type+Namevari[nvari]+"Run"+run+cuts+".eps");
}


TString round(double n, double d){
  if(d==0) return " - ";
  double b = ((int)(n/d *100+.5));
  b/=100.;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(b<10){
    if(result.Length() == 1)
      result += ".00";
    if(result.Length() == 3)
      result += "0";
  }
  return result;
}



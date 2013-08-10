{
TChain c("ntp1");
c.Add("AWG82/ntuples/Add_BDT_KM/Merged/*");

TString title[] = {"Signal","Dss","Crossfeed","Combinatoric"};
TCut signal = "MCType>0&&MCType<7";
TCut cross = "MCType>6&&MCType<13";
TCut comb = "MCType==0";
TCut dss = "MCType>12";
TCut cuts[4];
cuts[0] = signal;
cuts[1] = dss;
cuts[2] = cross;
cuts[3] = comb;

for (int i=0; i<4; i++){
  cuts[i] += "candType==1";
  cuts[i] += "candBMode>11000&&candBMode<12000";
  //cuts[i] += "candEExtra<.3";
} 
TString sigd = "(candDmass-1.862)^2/0.01398^2+";
TString mes = "(candMES-5.279)^2/0.0034^2+";
TString de = "(candDeltaE)^2/0.03033^2+";
TString eex = "candEExtra^2/0.286^2+";
TString tagd = "(candBTagDmass-1.864)^2/.0099^2+";
TString vari = tagd+sigd+mes+de+eex; vari.Resize(vari.Sizeof()-2);
TString vari = "candEExtra";
TString vari = "candMES:candDeltaE";
TString vari = "candBTagDmass:candDmass";

TCanvas p("chi2","Chi^2 Vs EExtra",800,600);
c.Draw(vari+">>d0",cuts[0]);
TH2F *d0 =  (TH2F*)gDirectory->Get("d0");
p.Clear();
double mean = d0->GetMean();
p.Divide(2,2);
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

for(int i=0; i<4; i++){
  p.cd(i+1);
  TString histo = vari; histo += ">>"; histo+=title[i]; //histo+="(30,0,"; histo+=5*mean; histo+=")";
  c.Draw(histo,cuts[i],"cont4z");
  d0 =  (TH2F*)gDirectory->Get(title[i]);
  d0->SetTitle(title[i]);
//   d0->GetXaxis()->SetTitle("E_{Extra} [GeV]");
//   d0->GetYaxis()->SetTitle("m_{ES} [GeV]");
  d0->GetXaxis()->SetTitle("D_{sig} mass [GeV]");
  d0->GetYaxis()->SetTitle("D_{tag} mass [GeV]");
  d0->Draw("cont4z");
}

// p.cd(1);
// c.Draw("((candBTagDmass-1.864)^2/.0099^2+(candDeltaE)^2/0.03033^2+(candMES-5.279)^2/0.0034^2+(candDmass-1.862)^2/0.01398^2+candEExtra^2/0.286^2)/21.54>>d0(30,0,2.5)",cuts[1]);
// TH1F *d0 =  (TH1F*)gDirectory->Get("d0");
// p.cd(2);

// c.Draw("((candBTagDeltam-0.142)^2/.0049^2+(candBTagDmass-1.864)^2/.0086^2+(candDeltaE)^2/0.03033^2+(candMES-5.279)^2/0.0035^2+(candDmass-1.862)^2/0.01398^2+candEExtra^2/0.23^2)/25.07>>ds(30,0,2.5)",cuts[2]);
// TH1F *ds =  (TH1F*)gDirectory->Get("ds");
// p.cd(3);

// c.Draw("((candBTagDmass-3.093)^2/.0185^2+(candDeltaE)^2/0.0276^2+(candMES-5.279)^2/0.0034^2+(candDmass-1.863)^2/0.0137^2+candEExtra^2/0.237^2)/18.88>>J(30,0,2.5)",cuts[0]);
// TH1F *J =  (TH1F*)gDirectory->Get("J");

// d0->Add(ds); 
// d0->Add(J);
// p.cd(4);
// d0->Draw();
// c.SetLineColor(2);
// c.Draw("candEExtra>>ee(30,0,2.5)",cuts[3],"same");

// TH1F *ee =  (TH1F*)gDirectory->Get("ee");
//d0->Draw("same");
// c.Draw("(candBTagDmass-1.864)^2/.0099^2+(candDeltaE)^2/0.03033^2+(candMES-5.279)^2/0.0034^2+(candDmass-1.862)^2/0.01398^2+candEExtra^2/0.286^2>>d0(50,0,100)","candType==1&&candBMode>11000&&candBMode<12000&&MCType>0&&MCType<7&&candEExtra<2");
// TH1F *d0 =  (TH1F*)gDirectory->Get("d0");

// c.Draw("(candBTagDeltam-0.142)^2/.0049^2+(candBTagDmass-1.864)^2/.0086^2+(candDeltaE)^2/0.03033^2+(candMES-5.279)^2/0.0035^2+(candDmass-1.862)^2/0.01398^2+candEExtra^2/0.23^2>>ds(50,0,100)","candType==1&&candBMode>14000&&candBMode<16000&&MCType>0&&MCType<7&&candEExtra<2");
// TH1F *ds =  (TH1F*)gDirectory->Get("ds");
// c.Draw("(candBTagDmass-3.093)^2/.0185^2+(candDeltaE)^2/0.0276^2+(candMES-5.279)^2/0.0034^2+(candDmass-1.863)^2/0.0137^2+candEExtra^2/0.237^2>>J(50,0,100)","candType==1&&candBMode<11000&&MCType>0&&MCType<7&&candEExtra<2");
// TH1F *J =  (TH1F*)gDirectory->Get("J");

// cout<<"d0: "<<d0->GetMean()<<", ds: "<<ds->GetMean()<<", J: "<<J->GetMean()<<endl;

}


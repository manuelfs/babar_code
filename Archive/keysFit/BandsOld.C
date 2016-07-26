void bands(TString sam = "19", TString smoo = "100", TString rot = "Rot", int Times = 10){
  TCanvas mm("Bands","Bands");
  TH2F *h2i[500];
  double tot;
  double m2min = -4, m2max = 12, plmin = 0, plmax = 2.4;
  int ntotbin = 800;
  double limits[] = {0, 1, 1.4, 1.8, 2.4};
  int binlim[5];
  for(int i=0;i<5;i++) binlim[i] = limits[i]/2.4*ntotbin;
  int plbinlim[3] = {0,ntotbin*5./16,ntotbin};
  double x[801][801], x2[801][801];
  double y[5][801], y2[5][801];
  double p[5][801], p2[5][801];
  for(int k=1; k<ntotbin+1; k++){
    for(int i=0; i<5; i++){
      y[i][k] = 0;
      y2[i][k] = 0;
      p[i][k] = 0;
      p2[i][k] = 0;
    }
    for(int j=1; j<ntotbin+1; j++){
      x[k][j] = 0;
      x2[k][j] = 0;
    }
  }
  TH1F h("h","Band",50,0,6);
  for(int i=0; i<Times; i++){
    TString fname = "AWG82/results/keys/Archive/root/new/hKeys"; fname += rot; fname += sam; fname += "_";
    fname += smoo; fname += "_"; fname += i; fname += ".root";
    if(i%30==0) cout<<"Doing "<<fname<<" of "<<Times<<endl;
    TFile f(fname);
    f.cd();
    h2i[i] = (TH2F *)gDirectory->Get("h2");
    h.Fill(h2[i]->GetBinContent(190,500)*100000);
    for(int k=1; k<ntotbin+1; k++){
      double sumPlep[] = {0,0,0,0,0};
      double sumM2[] = {0,0,0};
      for(int j=1; j<ntotbin+1; j++){
	double val = h2i[i]->GetBinContent(k,j);
	double valp = h2i[i]->GetBinContent(j,k);
	x[k][j] += val;
	x2[k][j] += val*val;
	for(int s=0; s<4; s++){
	  if(j>binlim[s] && j<binlim[s+1]+1) {
	    y[s][k] += val;
	    sumPlep[s] += val;
	  }
	  if(j==binlim[s+1]) y2[s][k] += sumPlep[s]*sumPlep[s];
	}
	y[4][k] += val;
	sumPlep[4] += val;
	for(int s=0; s<2; s++){
	  if(j>plbinlim[s] && j<plbinlim[s+1]+1) {
	    p[s][k] += valp;
	    sumM2[s] += valp;
	  }
	  if(j==plbinlim[s+1]) p2[s][k] += sumM2[s]*sumM2[s];
	}
	p[2][k] += valp;
	sumM2[2] += valp;
      }
      p2[2][k] += sumM2[2]*sumM2[2];
      y2[4][k] += sumPlep[4]*sumPlep[4];
    }
    if(i==0) tot = h2i[i]->Integral();
  }
  cout<<"Mean "<<x[190][500]/Times<<" and rms "<<sqrt(x2[190][500]/Times)<<endl;
  TString fname = "hBootstrap"; fname += sam; fname += "_";fname += smoo; fname += ".root";
  TFile* hBoot = new TFile(fname,"RECREATE"); 
  hBoot->cd();
  h.Write();hBoot->Close(); hBoot->Delete();
  fname = "AWG82/results/keys/Archive/root/single/hKeys"; fname += rot; fname += sam; fname += "_";fname += smoo; fname += ".root";
  //TFile f(fname);
  //f.cd();
  TString inputfile = "fitSamples/pdfSample"; inputfile += sam; inputfile += ".root";
  cout << "File = " << inputfile << endl;	
  TChain c("ntp1");
  c.Add(inputfile);
  double entries = c.GetEntries();
  Float_t xlow,xhigh;
  Int_t nbinx,nbiny,Sam = sam.Atoi(), ntotbin=800;
  xlow = m2min;
  xhigh = m2max;
  nbinx = 80;
  nbiny = 80;
  if (Sam==0 || Sam==2 || Sam==10 || Sam==12 || Sam == 20 || Sam == 23 || Sam == 26 || Sam == 29) {
    xlow = -2; xhigh = 4;
    if (Sam > 12) {nbinx = 40; nbiny = 40;}
  }
  else if (Sam==1 || Sam==11) {
    xlow = -2; xhigh = 6;
  }
  if (Sam==6 || Sam==7 || Sam==16 || Sam==17) {
    nbinx = 40; nbiny = 40;
  }
  if (Sam==8 || Sam==18) {
    xhigh = 4; nbinx = 40; nbiny = 40;
  }
  if (Sam==9 || Sam==19) {
    nbinx = 40; nbiny = 40;
  }
  if (Sam==21 || Sam==22 || Sam==24 || Sam==25 || Sam==27 || Sam==28 || Sam==30 || Sam==31) {
    xhigh = 8; nbinx = 40; nbiny = 20;
  }
  if (Sam > 31) {
    nbinx = 40; nbiny = 20;
  }

  TString M2titles[] = {"0 < p*_{l} < 1 GeV","1 < p*_{l} < 1.4 GeV","1.4 < p*_{l} < 1.8 GeV",
		      "1.8 < p*_{l} < 2.4 GeV","0 < p*_{l} < 2.4 GeV"};
  TString Pltitles[] = {"-4 < m^{2}_{miss} < 1 GeV^{2}","1 < m^{2}_{miss} < 12 GeV^{2}",
			"-4 < m^{2}_{miss} < 12 GeV^{2}"};
  TString M2cuts[] = {"candPstarLep<1","candPstarLep>1&&candPstarLep<1.4",
		      "candPstarLep>1.4&&candPstarLep<1.8","candPstarLep>1.8&&candPstarLep<2.4", ""};
  TString Plcuts[] = {"candM2<1.","candM2>=1.",""};

  TString psname = "AWG82/results/keys/Bands/epsBand2Keys"; psname+=rot; psname+=sam; psname+="_";
  psname += smoo; psname += ".ps";
  mm.Print(psname+"[");
  TH1F *hm2[5], *hm2Max[5], *hm2Min[5], *m2[5], *hpl[3], *pl[3], *hplMax[3], *hplMin[3];
  TString M2names[5], Plnames[3];
  for(int i=0;i<5;i++){
    M2names[i] = "hm2_"; M2names[i] += i;
    hm2[i] = new TH1F(M2names[i],M2titles[i],ntotbin,m2min,m2max); 
    M2names[i] += "Max";
    hm2Max[i] = new TH1F(M2names[i],M2titles[i],ntotbin,m2min,m2max); 
    M2names[i] += "Min";
    hm2Min[i] = new TH1F(M2names[i],M2titles[i],ntotbin,m2min,m2max); 
    if(i<3) {
      Plnames[i] = "hpl_"; Plnames[i] += i;
      hpl[i] = new TH1F(Plnames[i],Pltitles[i],ntotbin,plmin,plmax); 
      Plnames[i] += "Max";
      hplMax[i] = new TH1F(Plnames[i],Pltitles[i],ntotbin,plmin,plmax); 
      Plnames[i] += "Min";
      hplMin[i] = new TH1F(Plnames[i],Pltitles[i],ntotbin,plmin,plmax); 
    }
  }    
  for(int i=0;i<5;i++){
    TString hname = "m2"; hname += i;
    TString vari = "candM2>>"; vari+=hname; vari+="("; vari+= nbinx; vari+=",";vari+= xlow; 
    vari+=",";vari+= xhigh; vari+=")";
    c.Draw(vari,M2cuts[i]);
    m2[i] = (TH1F*)gDirectory->Get(hname);
    m2[i]->SetXTitle("m^{2}_{miss} [GeV^{2}]");
    m2[i]->SetTitle(M2titles[i]);
    m2[i]->Sumw2();
    m2[i]->SetMarkerStyle(20);
    m2[i]->SetMarkerSize(1);
    gStyle->SetOptStat(0);
    if(i<4){
      for(int j=1; j<ntotbin+1; j++){
	double minVal = 0, maxVal = 0, mean, Square, variance, rms;
	double factor = entries*ntotbin*(xhigh-xlow)/nbinx/(m2max-m2min)/tot;
	for(int binp = binlim[i]+1; binp < binlim[i+1]+1; binp++){
	  mean = x[j][binp]*factor/Times;
	  Square = x2[j][binp]*factor*factor;
	  variance = Square/(Times-1)-mean*mean*Times/(Times-1);
	  rms = 0;
	  if(variance>=0) rms = sqrt(variance);
	  else cout<<"Bin "<<j<<" and "<<binp<<" is wrong. Square "<<Square<<" and mean "<<mean<<endl;
	  minVal += mean - rms;
	  maxVal += mean + rms;
	}
	mean = y[i][j]*factor/Times;
	Square = y2[i][j]*factor*factor;
	variance = Square/(Times-1)-mean*mean*Times/(Times-1);
	rms = 0;
	if(variance>=0) rms = sqrt(variance);
	else cout<<"Bin "<<j<<" is wrong. Square "<<Square<<" and mean "<<mean<<endl;
	hm2Min[i]->SetBinContent(j,mean-rms);
	hm2Max[i]->SetBinContent(j,mean+rms);
	hm2[i]->SetBinContent(j,mean);
	//if(j%50==0)cout<<"mean "<<mean<<", rms "<<rms<<" at "<<j<<". Square "<<Square<<", mean "<<mean<<endl;
      }
    } 
    hm2[i]->SetLineColor(4);
    hm2[i]->SetLineWidth(2);
    hm2Max[i]->SetLineColor(9);
    hm2Min[i]->SetLineColor(9);
    hm2Max[i]->SetFillColor(33);
    hm2Min[i]->SetFillColor(10);
    hm2Max[i]->SetLineWidth(2);
    hm2Min[i]->SetLineWidth(2);
    if(i<4) {
      hm2[4]->Add(hm2[i]);
      hm2Max[4]->Add(hm2Max[i]);
      hm2Min[4]->Add(hm2Min[i]);
    }
  }
  int plbinlim[3] = {0,ntotbin*5./16,ntotbin};
  for(int i=0;i<3;i++){
    TString hname = "pl"; hname += i;
    TString vari = "candPstarLep>>"; vari+=hname; vari+="("; vari+= nbiny; vari+=",0,2.4)";
    c.Draw(vari,Plcuts[i]);
    pl[i] = (TH1F*)gDirectory->Get(hname);
    pl[i]->SetXTitle("p*_{l} [GeV]");
    pl[i]->SetTitle(Pltitles[i]);
    pl[i]->Sumw2();
    pl[i]->SetMarkerStyle(20);
    pl[i]->SetMarkerSize(1);
    gStyle->SetOptStat(0);
    if(i<2){
      for(int j=1; j<ntotbin+1; j++){
	double minVal = 0, maxVal = 0, mean, Square, variance, rms;
	double factor = entries*ntotbin/nbiny/tot;
	for(int binp = plbinlim[i]+1; binp < plbinlim[i+1]+1; binp++){
	  mean = x[binp][j]*factor/Times;
	  Square = x2[binp][j]*factor*factor;
	  variance = Square/(Times-1)-mean*mean*Times/(Times-1);
	  rms = 0;
	  if(variance>=0) rms = sqrt(variance);
	  else cout<<"Bin "<<j<<" and "<<binp<<" is wrong. Square "<<Square<<" and mean "<<mean<<endl;
	  minVal += mean - rms;
	  maxVal += mean + rms;
	}
	mean = p[i][j]*factor/Times;
	Square = p2[i][j]*factor*factor;
	variance = Square/(Times-1)-mean*mean*Times/(Times-1);
	rms = 0;
	if(variance>=0) rms = sqrt(variance);
	else cout<<"Bin "<<j<<" is wrong. Square "<<Square<<" and mean "<<mean<<endl;
	hplMin[i]->SetBinContent(j,mean-rms);
	hplMax[i]->SetBinContent(j,mean+rms);
	hpl[i]->SetBinContent(j,mean);
      }
    }
    hpl[i]->SetLineColor(4);
    hpl[i]->SetLineWidth(2);
    hplMax[i]->SetLineColor(9);
    hplMin[i]->SetLineColor(9);
    hplMax[i]->SetFillColor(33);
    hplMin[i]->SetFillColor(10);
    hplMax[i]->SetLineWidth(2);
    hplMin[i]->SetLineWidth(2);
    if(i<2) {
      hpl[2]->Add(hpl[i]);
      hplMax[2]->Add(hplMax[i]);
      hplMin[2]->Add(hplMin[i]);
    }
  }
  m2[4]->Draw("e1"); hm2Max[4]->Draw("c same"); hm2Min[4]->Draw("c same"); 
  m2[4]->Draw("e1 same"); hm2[4]->Draw("c same"); mm.Print(psname);
  pl[2]->Draw("e1"); hplMax[2]->Draw("c same"); hplMin[2]->Draw("c same"); 
  pl[2]->Draw("e1 same"); hpl[2]->Draw("c same"); mm.Print(psname);
  for(int i=0;i<4;i++){
    m2[i]->Draw("e1");
    hm2Max[i]->Draw("c same"); hm2Min[i]->Draw("c same"); 
    m2[i]->Draw("e1 same"); hm2[i]->Draw("c same"); 
    mm.Print(psname);
  }
  for(int i=0;i<2;i++){
    pl[i]->Draw("e1");
    hplMax[i]->Draw("c same"); hplMin[i]->Draw("c same"); 
    pl[i]->Draw("e1 same"); hpl[i]->Draw("c same"); 
    mm.Print(psname);
  }
  mm.Print(psname+"]");
  return; 
}

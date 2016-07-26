#include "babar_code/keysFit/PlotSingle.C"


void Stitch(TString hname, TString Tail = "70", double stitch1=0.4, double stitch2=1.5){
  TString sam = hname; sam.Remove(0,sam.First('_')+1);
  sam.Remove(sam.First('_'),sam.Length());
  TString smoo = hname;smoo.Remove(0,smoo.Last('_')+1);
  smoo.Remove(smoo.First('.'),smoo.Length());
  TString binString = hname; binString.Remove(binString.Last('_'),binString.Length());
  binString.Remove(0,binString.Last('_')+1);
  int nM2bin = binString.Atoi(); int nPlbin = nM2bin%1000;
  nM2bin = nM2bin/1000;
  TString basename(hname);
  Int_t slashpos = basename.Last('D');
  if (slashpos>=0) {
    basename.Remove(0,slashpos+1); basename.Remove(basename.First('.'),basename.Length());
  } else { cout<<hname<<" is incorrect"<<endl; return;}

  TFile hfile(hname); 
  TH2F *hTemp = (TH2F *)gDirectory->Get("h2");
  TH2F *hSmooth = (TH2F *)hTemp->Clone("hSmooth");
  TH2F *hResult = (TH2F *)hTemp->Clone("hResult");
  TString hname2 = hname; hname2.Remove(hname2.First('_'),hname2.Length());
  hname2 += "_"; hname2 += sam; hname2 += "_"; hname2 += nM2bin*1000+nPlbin;
  hname2 += "_"; hname2 += Tail; hname2 += ".root";
  TFile hfile2(hname2); hfile2.cd();
  TH2F *hTemp = (TH2F *)gDirectory->Get("h2");
  if(hTemp==0) {cout<<"File "<<hname2<<" does not exist"<<endl;return;}
  TH2F *hTail = (TH2F *)hTemp->Clone("hTail");

  int binCut = (int)(stitch1/2.4*(double)nPlbin);
  int binRange = nPlbin/20;
  if(binRange+1>binCut) binRange = binCut-1;
  double tailBin,tailBin2,smooBin,smooBin2,mean,avgInc,Offset;
  for(int i=1; i<nM2bin+1; i++){
      tailBin = hTail->GetBinContent(i,binCut);
      tailBin2 = hTail->GetBinContent(i,binCut-1);
      smooBin = hResult->GetBinContent(i,binCut+1);
      smooBin2 = hResult->GetBinContent(i,binCut+2);
      mean = (tailBin+smooBin)/2;
      avgInc = (tailBin-tailBin2+smooBin2-smooBin)/2;
      Offset = tailBin-mean+avgInc/2;
    if(i==200){
      cout<<"smooBin: "<<smooBin<<", tailBin: "<<tailBin<<", Offset: "<<Offset<<endl;
    }      
    for(int j=1; j<binCut+binRange+1; j++){
      double val=0;
      if(j<binCut+1) {
	val = hTail->GetBinContent(i,j)*1.15;
	int dist = j-binCut+binRange;
	//if(dist>0) val -= Offset*(double)dist/(double)binRange;
      }else  {
	val = hResult->GetBinContent(i,j);
	int dist = binRange+binCut+1-j;
	//val += Offset*(double)dist/(double)binRange;
      }
      hResult->SetBinContent(i,j,val);
    }
  }
  binCut = (int)(stitch2/2.4*(double)nPlbin);
  if(binCut+binRange>nPlbin) binRange = nPlbin-binCut-1;
  for(int i=1; i<nM2bin+1; i++){
      tailBin = hTail->GetBinContent(i,binCut);
      tailBin2 = hTail->GetBinContent(i,binCut+1);
      smooBin = hResult->GetBinContent(i,binCut-1);
      smooBin2 = hResult->GetBinContent(i,binCut-2);
      mean = (tailBin+smooBin)/2;
      avgInc = (tailBin2-tailBin+smooBin-smooBin2)/2;
      Offset = smooBin-mean+avgInc/2;
    for(int j=binCut-binRange+1; j<nPlbin+1; j++){
      double val=0;
      if(j>=binCut) {
	val = hTail->GetBinContent(i,j)*1.15;
	int dist = binRange+binCut-j;
	//if(dist>0) val += Offset*(double)dist/(double)binRange;
      }else  {
	val = hResult->GetBinContent(i,j);
	int dist = j-binCut+1+binRange;
	//val -= Offset*(double)dist/(double)binRange;
      }
      hResult->SetBinContent(i,j,val);
    }
  }
  hname.Remove(hname.First('.'),hname.Length()); hname += "and"; 
  hname += Tail; hname += ".root";
  hname.Replace(hname.Index("single"),6,"Stitch");
  TFile* fResult = new TFile(hname,"RECREATE"); 
  hResult->Write();
  fResult->Close();
  fResult->Delete();
  cout<<hname<<" saved"<<endl; 
  PlotSingle(hname);

}


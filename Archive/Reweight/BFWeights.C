TString RoundNumber(double n, int e, double d=1);

void BFWeights(TString toWeights, TString fromWeights="SP10"){

  TString wNames[] = {"weightD0","weightDs0","weightDp","weightDsp",
		      "weightD10","weightD20","weightD00","weightD1prime0",
		      "weightD0pi0","weightDstar0pi0","weightDpi","weightDstarpi",
		      "weightD1p","weightD2p","weightD0p","weightD1primep",
		      "weightD0pi","weightDstar0pi","weightDpi0","weightDstarpi0"};

  int nWeights = 4;
  TString wSets[] = {"SP10","Newest","Oct2010","NoNR"};
  double weights[4][20] = {{2.24, 6.17, 2.07, 5.70, 0.56, 0.30, 0.49, 0.90, 0.10, 0.03, 0.19, 0.06, 
			    0.52, 0.23, 0.45, 0.83, 0.20, 0.07, 0.10, 0.03},
			   {2.32, 5.48, 2.17, 5.11, 0.77, 0.59, 0.88, 0.82, 0.00, 0.00, 0.00, 0.00,
			    0.69, 0.56, 0.81, 0.76, 0.00, 0.00, 0.00, 0.00},
			   {2.28, 5.53, 2.10, 5.19, 0.43, 0.42, 0.42, 0.43, 0.70, 0.32, 0.25, 0.12,
			    0.40, 0.39, 0.38, 0.40, 0.30, 0.64, 0.11, 0.22},
			   {2.28, 5.53, 2.10, 5.19, 0.43, 0.42, 0.42, 0.43, 0.00, 0.00, 0.00, 0.00,
			    0.40, 0.39, 0.38, 0.40, 0.00, 0.00, 0.00, 0.00}};

  int iTo = -1, iFrom = -1;
  for(int i=0; i<nWeights; i++){
    if(wSets[i] == toWeights) iTo = i;
    if(wSets[i] == fromWeights) iFrom = i;
  }
  if(iTo<0 || iFrom<0){cout<<toWeights<<" or "<<fromWeights<<" not found"<<endl;return;}
  double totTo = 0, totFrom = 0;
  for(int dec=0; dec<20; dec++){
    if(dec==0 || dec==4 || dec==12) cout<<endl;
    cout<<wNames[dec]<<"\t= "<<RoundNumber(weights[iTo][dec],4,weights[iFrom][dec])<<endl;
    totTo += weights[iTo][dec];
    totFrom += weights[iFrom][dec];
  }
  cout<<endl<<"Total "<<toWeights<<": "<<RoundNumber(totTo,2)<<"\t Total "<<fromWeights
      <<": "<<RoundNumber(totFrom,2)<<endl<<endl;
}

TString RoundNumber(double n, int e, double d){
  if(d==0) return " - ";
  double neg = 1; if(n*d<0) neg = -1;
  double b = (int)(neg*n/d*pow(10.,(double)e)+0.5);
  b /= pow(10.,(double)e)*neg;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}



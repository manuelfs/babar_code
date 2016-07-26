TString RoundNumber(double n, int e, double d=1);
void Aver2(double R1, double Esta1, double Esys1, double R2, double Esta2, double Esys2);
void CalcRatio(double BTau, double staTau, double sysTau, double Bl, double sl);

void AverErrors(){
  cout<<endl<<"Belle 2007"<<endl;
  CalcRatio(2.02, 0.39, 0.37, 5.01, 0.12);

  cout<<endl<<"Belle 2009"<<endl;
  Aver2(0.70, 0.19, 0.10, 0.48, 0.21, 0.06);
  Aver2(0.47, 0.11, 0.07, 0.48, 0.13, 0.05);

  cout<<endl<<"Belle 2010"<<endl;
  CalcRatio(0.77, 0.22, 0.12, 2.23, 0.11);
  CalcRatio(2.12, 0.28, 0.29, 5.68, 0.19);

  cout<<endl<<"BaBar 2011"<<endl;
  cout<<RoundNumber(sqrt(pow(0.054,2)+pow(0.054,2))/0.457*100,1)<<"% error"<<endl;
  cout<<RoundNumber(sqrt(pow(0.023,2)+pow(0.023,2))/0.324*100,1)<<"% error"<<endl<<endl;
}

void Aver2(double R1, double Esta1, double Esys1, double R2, double Esta2, double Esys2){
  double k = TMath::min(Esys1/R1, Esys2/R2);
  double E1 = sqrt(Esys1*Esys1+Esta1*Esta1);
  double E2 = sqrt(Esys2*Esys2+Esta2*Esta2);
  //double E = 1/sqrt(1/E1/E1+1/E2/E2);
  double E = 1/sqrt(1/Esta1/Esta1+1/Esta2/Esta2);
  double R = (R1/Esta1/Esta1+R2/Esta2/Esta2)*E*E;
  double totE = sqrt(E*E+pow(k*R,2));

  //cout<<RoundNumber(R,2)<<" \\pm "<<RoundNumber(E,2)<<" \\pm "<<RoundNumber(k*R,2)<<" - "<<RoundNumber(totE/R*100.,0)<<endl;
  cout<<RoundNumber(R,2)<<" \\pm "<<RoundNumber(totE,2)<<" - "<<RoundNumber(totE/R*100.,0)<<"% error"<<endl;
}
void CalcRatio(double BTau, double staTau, double sysTau, double Bl, double sl){
  double sTau = sqrt(staTau*staTau + sysTau*sysTau);
  double R = BTau/Bl;
  double E = sqrt(pow(sTau/Bl,2)+pow(BTau/Bl/Bl*sl,2));

  cout<<RoundNumber(R,2)<<" \\pm "<<RoundNumber(E,2)<<" - "<<RoundNumber(E/R*100.,0)<<"% error"<<endl;
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

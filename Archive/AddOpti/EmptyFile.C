//---------------------------------------------------------------------------------
// Description:
//      Fifth and last stage of the mode optimization
//      Macro that makes the text file with the BReco mode that are not useful
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/05/23 manuelf -- Creation
//---------------------------------------------------------------------------------

void EmptyFile(TString truthName = "babar_code/txt/truth12346.txt"){

  TString emptyName = "babar_code/txt/BSemiExclAdd_Empty.txt";
  fstream emptyFile;
  emptyFile.open(emptyName,fstream::out);
  fstream truthFile;
  truthFile.open(truthName,fstream::in);
  int SeedMode, YMode, type[90][54][6], S, N, Dss, Cr, Co;
  double Tmatched[90][54][6];
  for(int i=0; i<87; i++){
    if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
    if(i==49) i = 50;if(i==59) i = 60; if(i==62) i = 71;
    if(i==73) i = 80;if(i==85) i = 86;
    for(int j=1; j<54; j++){
      for(int k=0; k<6; k++){
	Tmatched[i][j][k] = 0;
      }
    }
  }
  TString rub;
  int Track, Neutral;
  double Tpuri, Tave, Npuri, Nave;
  while(truthFile){
    truthFile>>SeedMode>>YMode>>rub>>Track>>Tpuri>>Tave>>Neutral>>Npuri>>Nave;
    Tmatched[SeedMode][YMode][0] = (double)Track;
    Tmatched[SeedMode][YMode][1] = Tpuri;
    if(Track<3) Tmatched[SeedMode][YMode][1] = 0;
    Tmatched[SeedMode][YMode][2] = Tave;
    Tmatched[SeedMode][YMode][3] = (double)Neutral;
    Tmatched[SeedMode][YMode][4] = Npuri;
    Tmatched[SeedMode][YMode][5] = Nave;
  }
  for(int i=0; i<87; i++){
    int nModes = 0;
    if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
    if(i==49) i = 50;if(i==59) i = 60; if(i==62) i = 71;
    if(i==73) i = 80;if(i==85) i = 86;
    for(int j=1; j<54; j++){
      if(Tmatched[i][j][1] < 0.3){
	int mode = 10000 + i*100 + j;
	emptyFile<<mode<<" ";
	nModes++;
	if(nModes==12){
	  nModes=0;
	  emptyFile<<endl;
	}
      }
    }
    emptyFile<<endl;
    if(nModes>0) emptyFile<<endl;
  }

  cout<<"Empty modes stored in "<<emptyName<<endl;
}



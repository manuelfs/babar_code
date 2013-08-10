
#define nBDT 9

void Dchannels(int fixD = 1){
  TString BDTNames[] = {"030", "050", "075", "100", "120", "150", "200", "250", "300"};
  TString FixinFolder = "keys/root/fitOrix100", inFile, outFile, outFolder, inFolder;
  for(int iBDT = 0; iBDT < nBDT; iBDT++){
    outFolder = "keys/root/fit"; inFolder = "keys/root/fitOrix"; 
    if(fixD==1) outFolder += "Dfix";
    else outFolder += "Dsfx";
    outFolder += BDTNames[iBDT]; inFolder += BDTNames[iBDT];
    gSystem->mkdir(outFolder);
    for(int ind=1; ind<=68; ind++){
      if(ind==49) ind=51; if(ind==59) ind=61;   //Skipping non-existent numbers, and the signal combinatoric
      if(ind<=40&&ind%2==fixD || ind>40&&((ind%10)%2)==fixD) inFile = FixinFolder;
      else inFile = inFolder;
      inFile += "/pdfKeys_"; inFile += ind; inFile += "_Fit.root";
      outFile = outFolder; outFile += "/pdfKeys_"; outFile += ind; outFile += "_Fit.root";
      int result = gSystem->CopyFile(inFile, outFile, true);
      if(result!=0) cout<<inFile<<" copied into "<<outFolder<<" with result "<<result<<endl;
    }
    cout<<outFolder<<" done"<<endl;
  }
}


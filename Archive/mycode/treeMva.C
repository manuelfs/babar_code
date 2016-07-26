void treeMva(TString folder){
  Int_t slashpos = folder.Last('/');
  if (slashpos==folder.Sizeof()-2) {
    folder.Remove(slashpos,slashpos+1);     
  } 
  TString outdir = folder; outdir+="1/";
  gSystem->mkdir(outdir);
  folder += "/*";
  TChain c("ntp1");
  c.Add(folder);
  c.SetBranchStatus("*",0);
  c.SetBranchStatus("candEExtra",1);
  c.SetBranchStatus("candMES",1);
  c.SetBranchStatus("candDmass",1);
  c.SetBranchStatus("candRejectedTracks",1);
  c.SetBranchStatus("candRejectedPhotons",1);
  c.SetBranchStatus("candTagChargedMult",1);
  c.SetBranchStatus("candDeltam",1);
  c.SetBranchStatus("candBTagDeltam",1);
  c.SetBranchStatus("candBTagDmass",1);
  c.SetBranchStatus("candDeltaE",1);
  c.SetBranchStatus("MCType",1);
  c.SetBranchStatus("candType",1);
  c.SetBranchStatus("candLepTru",1);

  outdir += "mergedMva.root";
  cout<<"Merging "<<outdir<<endl;
  c.Merge(outdir);

}


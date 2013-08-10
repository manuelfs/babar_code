void CountModes(string Ename){
  ifstream Empty(Ename.c_str());
  int  mode, n=0;
  while(!Empty.eof() && Empty){
    Empty >> mode;
    n++;
  }
  int m = 0;
  for(int i=0; i<87; i++){
    if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
    if(i==49) i = 50;if(i==59) i = 60; if(i==62) i = 71;
    if(i==73) i = 80;if(i==85) i = 86;
    for(int j=1; j<54; j++){
      m++;
    }
  }

  cout<<"Killed "<<n<<" modes of "<<m<<" total"<<endl;
}

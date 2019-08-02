//TKey* loopdirGraph(TDirectory* dir, TString gname);

Double_t getSlope(TString kchan);



void adcPrinter(int strip){
  
  /*
   * usage: root 'adcPrinter(306)'
   * could be changed: 
   *  -- folder  KPiX_26
   */
  TH1::AddDirectory(kFALSE);

  // HGT, 
  TString fns = "/scratch/data/tracker_test/2018_09_06_15_48_39.bin.root";
  TString fnbg = "/scratch/data/tracker_test/2018_09_05_15_31_31.bin.root";
  
  /*// Aug: extTrig,
  TString fns = "/scratch/data/tracker_test/eudaq/testbeam_201808/2018_08_17_20_30_20.bin.root";
  TString fnbg = "/scratch/data/tracker_test/2018_08_19_21_34_17.bin.root";
  */
  
  TFile* fs = TFile::Open(fns.Data());
  TFile* fbg = TFile::Open(fnbg.Data());

  TObject *obj;
  TKey *key;

  fs->cd("KPiX_26/Strips_and_Channels");
  TIter next1( gDirectory->GetListOfKeys());
  TString hDir;
  while ((key = (TKey *) next1())) {
    // obj = fs->Get(key->GetName()); // copy object to memory
    // do something with obj
    hDir = (TString)key->GetName();
    if (hDir.Contains(Form("strip_%d_",strip))){
       printf(" found object:\t%s\n",key->GetName());
       break;
    }
    
  }
  
  printf(" found it:%s\n", hDir.Data());
  hDir = "KPiX_26/Strips_and_Channels/"+hDir+"/bucket_0/";

  fs->cd(hDir);
  TIter next2( gDirectory->GetListOfKeys());
  TString hist;
  while ((key = (TKey *) next2())) {
    // do something with obj
    hist = (TString)key->GetName();
    if (hist.Contains("hist_s")){
      printf(" found object:\t%s\n",key->GetName());
      break;
    }
  }

  TString kchan = hist;
  kchan.ReplaceAll("_b0_k26","");
  kchan.ReplaceAll(Form("hist_s%d_c",strip),"");
  printf(" kpix channel:\t%s\n",kchan.Data());
  double slope = getSlope(kchan);
  
  hist = hDir + hist;

  //  cout<< " debug: " << hist<< endl;
  TH1F* hs = nullptr;
  TH1F* hbg = nullptr;
  
  fs->GetObject(hist,hs);
  fbg->GetObject(hist,hbg);
   
  fs->Close();
  fbg->Close();
  
  if (hs==nullptr && hbg==nullptr) return;
  
  hs->Scale(1/hs->Integral());
  hbg->Scale(1/hbg->Integral());
  
  hs->SetMarkerColor(kRed);
  hbg->SetMarkerColor(kBlue);
  
  hs->SetLineColor(kRed);
  hbg->SetLineColor(kBlue);
  hs->Draw();
  hbg->Draw("same");
  
  hs->GetXaxis()->SetRange(500,800);
 
}



Double_t loopdirGraph(TDirectory* dir, TString gname) {
	
  TDirectory *dirsav = gDirectory;
  TIter keys_iter(dir->GetListOfKeys());
  TKey* key;
  Double_t slope=-5;
 
  while ((key = (TKey*)keys_iter())){   
    if (key->IsFolder()){
      // it is a subdir
      //printf(" found calib subdir:\t%s\n",key->GetName());
	
      dir->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      slope = loopdirGraph(subdir, gname);
      dirsav->cd();
      continue;
    }
    else{
      TString keyname = key->GetName();
      if (keyname.Contains(gname))  {
	printf(" found calib object:\t%s\n",key->GetName());
	
	TGraphErrors* graph = (TGraphErrors*)key->ReadObj();
	slope = graph->GetFunction("pol1")->GetParameter(1);
	cout<< "\t debug sub slope = "<< slope<< endl;
	break;
      }
    }
  }
  
  cout<< "\t debug func slope = "<< slope<< endl;

  return slope;
}

Double_t getSlope(TString kchan){

  // HGT:
  TFile* file = TFile::Open("/scratch/data/tracker_test/2018_09_05_11_06_13.bin.duoCalib.root");
  TString gname = "calib_1026_c"+kchan+"_b0_r0";

  double slope =  loopdirGraph(file, gname);
  
  cout<< " debug:\t" << slope<< endl;
  return slope;
}

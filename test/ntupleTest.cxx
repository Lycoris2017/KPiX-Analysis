#include "ntupleMaker.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TTree.h"
#include "TApplication.h"

using namespace std;

void treer(){
  TFile *f = TFile::Open("kpixTree.root");
  
  TTree*t1 = (TTree*)f->Get("General");
  TH1F *h1   = new TH1F("h1","CycleNr",7000,83000,90000);

  int cycleNr=0;
  t1->SetBranchAddress("cycleNr", &cycleNr);

  auto nentries = t1->GetEntries();
  for (Long64_t i =0; i<nentries; i++){

    t1->GetEntry(i);
    h1->Fill(cycleNr);
  }
  
  if (gROOT->IsBatch()) {
    // f->Close();
    return;
  }
  new TBrowser;
  t1->StartViewer();
}
  

int main( int argc, char **argv){
  lycoris::ntupleMaker kk(argv[1]);
  //kk.loopKpix();
  kk.setDebug(true);
  kk.CreateTree();

  //  kk.makeTreeEx();

  // TApplication *myapp=new TApplication("myapp",0,0);
  
  //treer();
  //myapp->Run();
  
  return(1);

}

/*
 *  Simple ROOT macro to check kpix channel w/ evts < 3.8 per Acq. Cycle
 *
 *  Usage: root -l 'debugger.C("xxxx.bin.root")'
 *
 */

void debugger(string infile){

  if (infile.empty()) infile ="/scratch/data/tracker_test/2018_03_22_11_16_41.bin.root";
  cout << " open file: "<< infile <<endl;
  
  TFile* f1=TFile::Open(infile.c_str());
  TH1F* h1 = (TH1F*)f1->Get("General/Channel_entries_total");
  if (!h1) {
    cout<< "Either your file broken or your histogram name is shitty!"<<endl;
    return;
  }

  // -- Check if a kpix channel has events less than 3.8 per Acq. Cycle:
  int nbins = h1->GetSize();
  TAxis *xaxis = h1->GetXaxis();
  for (int nbin = 1; nbin<nbins-1; nbin++){
    auto strip_evts=h1->GetBinContent(nbin);

    if (strip_evts<3.8 && strip_evts!=0) {
      auto strip_ID=xaxis->GetBinCenter(nbin);
      cout<< "[debug] stripId = "<< strip_ID
	  << ", entries = "<<strip_evts<<endl;
    }

  }


 return;
}

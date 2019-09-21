#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TCollection.h"
#include "TKey.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TLatex.h"

using namespace RooFit;

Double_t ScaleX(Double_t x);
void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t));

double Median(const TH1D * h1);
void StatBox(Double_t x,Double_t y,Color_t color, std::vector<TString> text);

std::vector<TH1F*> FindInHist(TFile *file);
std::vector<TH1F*> FindInHist(TFile *file, TFile *file2);
int PlotStack(TFile *file, TFile *file2);
int PlotSingle(TFile *file, TString hist, bool doFit, TString fitfunc );
int PlotFitres(TFile *file);

int main(int argc, char*argv[]){
	TH1::AddDirectory(kFALSE);
		
	printf("hello world, I have %d args\n", argc);
	TString filename;
	if (argc>0){
		printf("open file : %s\n", argv[1]);
		filename = TString(argv[1]);
	}
	else{
		TString filedir = "~/Documents/paper/vci2019proc/rooPlots/selfTrig/";
		filename = filedir + "2018_08_17_16_07_45.bin.root";
	}
	TFile *file = new TFile(filename, "READ");
	
	TFile *file2= 0;
	if (argc==3) {
		printf("open file2 : %s\n", argv[2]);
		file2 = new TFile(argv[2], "READ");
	}

	gROOT->LoadMacro("./rooMacros/AtlasStyle.C");
	gROOT->LoadMacro("./rooMacros/AtlasLabels.C");

	gROOT->ProcessLine("SetAtlasStyle()");
	// get rid of gDirectory dependence

	TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 500);

	//** !! START: variables to check !! **
	TString subdir = "";
	TString hname  = "";
	/*
	subdir = "KPiX_26/";
	hname = "timed_left_strip_entries_k26_b0";
	hname = "Strips_and_Channels/strip_490_channel_395/hist_timed_s490_c0395_b0_k26";
	*/
	//hname = "fc_response_median_made_CMmedian_subtracted_k26_total";
	//hname = "cluster_charge_cuts_sens0_b0";
	//PlotSingle(file, subdir+hname, true, "landau");
	//** !! END: variables to check !! ** 

	//PlotStack(file, "490");
	//PlotStack(file, file2);

	PlotFitres(file); 
	
	//gROOT->ProcessLine("LycorisLabel(0.48,0.75,\"Preliminary\",1);");
	//gROOT->ProcessLine("LycorisLabel(0.2,0.85,\"Preliminary\",1);");
	gROOT->ProcessLine("LycorisLabel(0.55,0.6,\"Preliminary\",1);");

	
	canvas->SaveAs("out.eps");
	file->Close();
	delete file;
	return(1);
}

int PlotSingle(TFile *file, TString hist, bool doFit=false, TString fitfunc="" ){
	TH1::SetDefaultSumw2();

	TH1F* h1=0;
	
	file->GetObject(hist, h1);
	if (h1 == 0) {
		std::cerr << "Error: did not find hist object called - "<< hist << std::endl;
		return (0);
	}

	h1->Scale(1./1514);
	
	if (doFit){
		if (fitfunc=="") fitfunc="gaus";
		if (h1->GetFunction("gaus")) h1->GetFunction("gaus")->Delete();
		if (h1->GetFunction("landau")) h1->GetFunction("landau")->Delete();
		
		double mean = h1->GetMean();
		TH1D h1D; h1->Copy(h1D);
		double median = Median(&h1D);
		double rms = h1->GetRMS();
		//h1->Fit(fitfunc.Data(), "RqW", "same", median-0.45*rms , median+0.45*rms);
		h1->Fit(fitfunc.Data(), "RWW", "same", 2 , 17);
	}
	
	h1->SetLineColor(1);
	h1->SetMarkerSize(1.0);
	h1->Draw("histe");
	TF1* fitres = h1->GetFunction(fitfunc.Data());
	if ( fitres ){
		fitres->SetLineColor(kRed);;
		fitres->Draw("same");
	}
	
	h1->GetYaxis()->SetTitle("Normalized / acq. cycle");
	h1->GetXaxis()->SetTitle("Charge [fC]");
	h1->GetXaxis()->SetRangeUser(-0.14, 17);
	h1->GetYaxis()->SetRangeUser(h1->GetMinimum(), 1.1*h1->GetMaximum());

	std::vector<TString> text;
	text.push_back("#mu = 0.00546 #pm 0.001742");
	text.push_back("#sigma = 0.2806 #pm 0.0019");
	//StatBox(0.58,0.75,1,text);
	
	return(1);
}

int PlotFitres(TFile *file){
	RooPlot* rplot=0;
	//TString rpname = "frame_fitTo_langauss";
	TString rpname = "frame_fitTo_lanOnly";
	rplot=(RooPlot*)file->Get(rpname);

	if (rplot == 0) {
		std::cerr << "Error: did not find hist object called - "<< rpname << std::endl;
		return (0);
	}

    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.16);
    gPad->SetTopMargin(0.05);
	rplot->Draw();
		
	return(1);
}




Double_t ScaleX(Double_t x){
	Double_t v;
	v= 10*x ;
	return v;
}
void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t))
{
  if (!a) return; // just a precaution
  if (a->GetXbins()->GetSize())
    {
      // an axis with variable bins
      // note: bins must remain in increasing order, hence the "Scale"
      // function must be strictly (monotonically) increasing
      TArrayD X(*(a->GetXbins()));
      for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
      a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
    }
  else
    {
      // an axis with fix bins
      // note: we modify Xmin and Xmax only, hence the "Scale" function
      // must be linear (and Xmax must remain greater than Xmin)
      a->Set( a->GetNbins(),
              Scale(a->GetXmin()), // new Xmin
              Scale(a->GetXmax()) ); // new Xmax
    }
  return;
}
double Median(const TH1D * h1) { 

   int n = h1->GetXaxis()->GetNbins();  
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray(); 
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]); 
}

void StatBox(Double_t x,Double_t y,Color_t color, std::vector<TString>text )
{
  TLatex* stat = new TLatex;
  stat->SetNDC();
  stat->SetTextFont(42);
  stat->SetTextSize(0.04);
  stat->SetTextColor(color);

  TString str="#splitline";
  for (auto && itex : text){
	  TString tex ="{" + itex + "}";
	  str+=tex;
  }
  printf("%s\n",str.Data());
  
  stat->DrawLatex(x,y,str.Data());
}

std::vector<TH1F*> FindInHist(TFile *file){
	TKey *key=0;
	TH1F* h1=0, *h2=0;
	//** !! START: variables to check !! **
	TString strip = "490";
	TString subdir = "KPiX_26/Strips_and_Channels/";
	TString hname1 = "hist_s"+strip;
	TString hname2 = "hist_timed_s"+strip;
	//** !! END: variables to check !! **

	//** find hist gDirectory
	TString hDir;
	file->cd(subdir);
	TIter next_d( gDirectory->GetListOfKeys() );
	while ( (key = (TKey*)next_d()) ){
		hDir = (TString)key->GetName();
		if( hDir.Contains(Form("strip_%s_", strip.Data()) ) ) break;
	}
	printf(" found!: %s\n", hDir.Data());
	hDir = subdir + hDir;

	//** find hists
	file->cd(hDir);
	TString hist;
	TIter next_h( gDirectory->GetListOfKeys() );
	while(( key = (TKey*)next_h()) ){
		hist = (TString)key->GetName();
		if ( hist.Contains(hname1.Data()) )
			h1 = (TH1F*)key->ReadObj();
		if ( hist.Contains(hname2.Data()) )
			h2 = (TH1F*)key->ReadObj();
	}
	if (h1!=0 && h2!=0)
		printf("found! %s\t%s\n", h1->GetName(), h2->GetName());

	std::vector<TH1F*> res;
	res.push_back(h1);
	res.push_back(h2);
	return res;
}

std::vector<TH1F*> FindInHist(TFile *file, TFile *file2){
	/*
	 * look for h1, h2 from two files in same structure
	 */

	TH1::AddDirectory(kFALSE);

	TString subdir = "KPiX_26/Buckets/Bucket_0/";
	TString hname  = "pedestalsRMS_fc_conn_k26_b0";

	TH1F* h1 = (TH1F*)file->Get(subdir+hname);
	//file2->GetObject(subdir+hname, h2);
	TH1F* h2 = (TH1F*)file2->Get(subdir+hname);
		
	if ( h1!=0 && h2!=0)
		printf("found from two files! %s\t%s\n", h1->GetName(), h2->GetName());
	std::vector<TH1F*> res;
	res.push_back(h1);
	res.push_back(h2);
	return res;
}



int PlotStack(TFile *file, TFile *file2=nullptr){
	THStack *hstack = new THStack("hstack", "stack histograms");
	std::vector<TH1F*> hvec;
	double ymax=105, ymin =0.;
	double xmax=1.02, xmin =0.12;
	
	if (file2==0)
		hvec = FindInHist(file);
	else
		hvec = FindInHist(file, file2);

	auto h1 = hvec.at(0);
	auto h2 = hvec.at(1);
	
	printf("found from two files! %s\t%s\n", h1->GetName(), h2->GetName());
	
	if (h1==0 || h2==0) {
		printf("!! NOT found histograms, check ur input.\n\n");
		return (0);
	}

	//** add to THStack
	h1->SetLineColor(1);
	h1->Rebin(4);
	h2->SetLineColor(2);
	h2->Rebin(4);
	hstack->Add(h1);
	hstack->Add(h2);

	//** if fancy axis used, you need space for top margin
	gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.20);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.10);
    
	hstack->Draw("histnostack");
	hstack->GetXaxis()->SetRangeUser(xmin, xmax);
	hstack->GetXaxis()->SetTitleOffset(1.2);
	hstack->GetYaxis()->SetTitleOffset(1.2);
	hstack->GetXaxis()->SetTitle("pedestal RMS [fC]");
	hstack->GetYaxis()->SetTitle("Channels");
	hstack->SetMinimum(ymin);
	hstack->SetMaximum(ymax);

	//** Draw fancy axes
	//hstack->GetXaxis();
	TGaxis *axis = new TGaxis(xmin, ymax, // xmin, ymin
	                          xmax, ymax, // xmax, ymax
	                          xmin, // wmin
	                          xmax,  // wmax
	                          002, "-");
	axis->SetLineColor(kBlue);
	axis->SetLabelFont(62);
	axis->SetLabelSize(0.03);
	axis->SetLabelColor(kBlue);
	axis->SetTitleFont(62);
	axis->SetTitleSize(0.03);
	axis->SetTitleOffset(1.3);
	axis->SetTitleColor(kBlue);
	axis->SetTitle("electrons");

	axis->ChangeLabel(-1, -1, -1, -1, -1, -1, "6242");
	axis->ChangeLabel(1, -1, -1, -1, -1, -1, "3121");

		
	axis->Draw();
	
	//** add legend
	TLegend* leg = new TLegend(0.52,0.58,0.76,0.70); //xmin, ymin, xmax, ymax
	leg->SetFillStyle(0); 
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->SetLineColor(0);
	
	leg->AddEntry(h1, "E-lab", "f");
	leg->AddEntry(h2, "Test Beam", "f");
	leg->Draw("same");

	return(1);
}

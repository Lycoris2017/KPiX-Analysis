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

using namespace RooFit;
int PlotSingle(TFile *file, TString hist, bool doFit=false ){
	TH1F* h1=0;
	
	file->GetObject(hist, h1);
	if (h1 == 0) {
		std::cerr << "Error: did not find hist object called - "<< hist << std::endl;
		return (0);
	}

	h1->SetLineColor(2);
	h1->SetMarkerSize(1.0);
	h1->Draw("e");
	h1->GetYaxis()->SetTitle("Normalized / acq. cycle");
	h1->GetXaxis()->SetTitle("Strip location");
	h1->GetXaxis()->SetRangeUser(150, 250);
	h1->GetYaxis()->SetRangeUser(h1->GetMinimum(), 1.4*h1->GetMaximum());

	return(1);
}

int PlotFitres(TFile *file){
	RooPlot* rplot=0;
	TString rpname = "frame_fitTo_langauss";
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
	h1->Rebin(2);
	h2->SetLineColor(2);
	h2->Rebin(2);
	hstack->Add(h1);
	hstack->Add(h2);
	
	hstack->Draw("histnostack");
	hstack->GetXaxis()->SetRangeUser(0.1,1.0);
	hstack->GetXaxis()->SetTitle("pedestal RMS [fC]");
	hstack->GetYaxis()->SetTitle("Channels");
	hstack->SetMinimum(0);
	hstack->SetMaximum(60);

	//** add legend
	
	TLegend* leg = new TLegend(0.52,0.68,0.76,0.80);	//TLegend(0.48,0.6,0.72,0.72);
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
	//TString subdir = "KPiX_26/";
	//TString hname = "timed_left_strip_entries_k26_b0";
	//TString hname = "Strips_and_Channels/strip_490_channel_395/hist_timed_s490_c0395_b0_k26";
	//PlotSingle(file, subdir+hname, true);
	//** !! END: variables to check !! ** 

	//PlotStack(file, "490");
	PlotStack(file, file2);

	//PlotFitres(file); 
	
	//gROOT->ProcessLine("LycorisLabel(0.48,0.85,\"Preliminary\",1);");
	//gROOT->ProcessLine("LycorisLabel(0.2,0.85,\"Preliminary\",1);");
	gROOT->ProcessLine("LycorisLabel(0.55,0.6,\"Preliminary\",1);");

	canvas->SaveAs("out.eps");
	file->Close();
	delete file;
	return(1);
}

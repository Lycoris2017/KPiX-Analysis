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

int PlotSingle(TFile *file){
	TH1F* h1=0;
	//** !! START: variables to check !! ** 
	TString subdir = "KPiX_26/";
	TString hname = "timed_left_strip_entries_k26_b0";
	//** !! END: variables to check !! ** 
	file->GetObject(subdir+hname, h1);
	if (h1 == 0) {
		std::cerr << "Error: did not find hist object called - "<< subdir+hname << std::endl;
		return (0);
	}

	h1->Draw("hist");
	h1->GetYaxis()->SetTitle("Normalized / acq. cycle");
	h1->GetXaxis()->SetTitle("Strip location");

	return(1);
}

int PlotStack(TFile *file, TString strip){
	THStack *hstack = new THStack("hstack", "stack histograms");
	
	TKey *key=0;
	TH1F *h1=0;
	TH1F *h2=0;
	//** !! START: variables to check !! ** 
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

	//** add to THStack
	h1->SetLineColor(1);
	h2->SetLineColor(2);
	hstack->Add(h1);
	hstack->Add(h2);

	hstack->Draw("histnostack");
	hstack->GetXaxis()->SetRangeUser(150, 230);
	hstack->GetXaxis()->SetTitle("Charge [ADC]");
	hstack->GetYaxis()->SetTitle("Normalized / acq. cycle");

	//** add legend
	
	TLegend* leg = new TLegend(0.48,0.68,0.72,0.80);	//TLegend(0.48,0.6,0.72,0.72);
	leg->SetFillStyle(0); 
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->SetLineColor(0);
	
	leg->AddEntry(h1, "All Self-Trigs", "f");
	leg->AddEntry(h2, "Matched Ext.Trig", "f");
	leg->Draw("same");
	return(1);
}

int main(){
	printf("hello world\n");
	gROOT->LoadMacro("./rooMacros/AtlasStyle.C");
	gROOT->LoadMacro("./rooMacros/AtlasLabels.C");
	gROOT->ProcessLine("SetAtlasStyle()");

	// get rid of gDirectory dependence
	TH1::AddDirectory(kFALSE);

	TString filedir = "~/Documents/paper/vci2019proc/rooPlots/selfTrig/";
	TString filename = "2018_08_17_16_07_45.bin.root";
	TFile *file = new TFile(filedir + filename, "READ");
	TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 500);

	//PlotSingle(file);
	PlotStack(file, "490");

	gROOT->ProcessLine("LycorisLabel(0.48,0.85,\"Preliminary\",1);");
		
	canvas->SaveAs("out.eps");
	file->Close();
	delete file;
	delete canvas;
	return(1);
}

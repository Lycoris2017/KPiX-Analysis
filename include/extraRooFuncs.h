/*
  Mengqing Wu <mengqing.wu@desy.de>
  @ 2020 Apr 15
  Project: Lycoris
  Target:  useful extra funcs for ploting with ROOT
*/
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

//! definition of funcs
double Median(const TH1D * h1);
void StatBox(Double_t x,Double_t y,Color_t color, TString text);
void StatBox(Double_t x,Double_t y,Color_t color, std::vector<TString> text);
Double_t ScaleX(Double_t x);
void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t));
void test(){
	std::cout << "hello world" << std::endl;
}
//! construction of funcs
double Median(const TH1D * h1) { 

   int n = h1->GetXaxis()->GetNbins();  
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray(); 
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]); 
}

void StatBox(Double_t x,Double_t y,Color_t color, TString text )
{
  TLatex* stat = new TLatex;
  stat->SetNDC();
  stat->SetTextFont(42);
  stat->SetTextSize(0.04);
  stat->SetTextColor(color);

  printf("%s\n",text.Data());
  stat->DrawLatex(x,y,text.Data());
}

void StatBox(Double_t x,Double_t y,Color_t color, std::vector<TString> text)
{
  TString str="#splitline";
  for (auto && itex : text){
	  TString tex ="{" + itex + "}";
	  str+=tex;
  }

  StatBox(x,y,color,str);
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

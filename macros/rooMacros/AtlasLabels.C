#include "AtlasLabels.h"

#include "TLatex.h"
#include "TLine.h"
#include "TPave.h"
#include "TPad.h"
#include "TMarker.h"


void LycorisLabel(Double_t x,Double_t y,const char* text,Color_t color) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);

  //double delx = 0.115*696*gPad->GetWh()/(472*gPad->GetWw());
  double delx = 0.150*696*gPad->GetWh()/(472*gPad->GetWw());
  
  l.DrawLatex(x,y,"LYCORIS");
  if (text) {
    TLatex p; 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+delx,y,text);
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }
}


void ATLASLabelOld(Double_t x,Double_t y,bool Internal,Color_t color) // replaced Preliminary to Internal
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x,y,"ATLAS");
  if (Internal) {
    TLatex p; 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+0.115,y,"Internal");
  }
}



void ATLASVersion(const char* version,Double_t x,Double_t y,Color_t color) 
{

  if (version) {
    char versionString[100];
    sprintf(versionString,"Version %s",version);
    TLatex l;
    l.SetTextAlign(22); 
    l.SetTextSize(0.04); 
    l.SetNDC();
    l.SetTextFont(72);
    l.SetTextColor(color);
    l.DrawLatex(x,y,versionString);
  }
}

void lumi2012(Double_t x,Double_t y,Color_t color, const char* text)
{
  TLatex* lumi = new TLatex;
  lumi->SetNDC();
  lumi->SetTextFont(42);
  lumi->SetTextSize(0.05);
  lumi->SetTextColor(color);
  if(text){
    TString str="#splitline{#surds = 8 TeV, #intL dt = 20.3 fb^{-1}}{";
    str+=text;
    str+="}";
    lumi->DrawLatex(x,y,str);
  }
  else  lumi->DrawLatex(x,y,"#splitline{#intL dt = 20.3 fb^{-1}}{#surds = 8 TeV}");

}

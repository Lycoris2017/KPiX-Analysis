//
// LycorisStyle, based on a style file from ATLAS Style
//   Mengqing Wu <mengqing.wu@desy.de>
//


#include <iostream>
#include "LycorisStyle.h"
#include "TROOT.h"

void SetMyGridStyle ()
{
  static TStyle* myGridStyle = 0;
  std::cout << "\n Applying MyGrid style settings...\n" << std::endl ;
  if ( myGridStyle==0 ) myGridStyle = MyGridStyle();
  gROOT->SetStyle("myGRID");
  gROOT->ForceStyle();
}

TStyle* MyGridStyle() 
{
  TStyle *myGridStyle = new TStyle("myGRID","my grid style");

  myGridStyle->SetPaintTextFormat("5.3f");//"5.2f");

  myGridStyle->SetPadTopMargin(0.05);
  myGridStyle->SetPadRightMargin(0.15);
  myGridStyle->SetPadBottomMargin(0.16);
  myGridStyle->SetPadLeftMargin(0.16);

  myGridStyle->SetMarkerSize(2);
  myGridStyle->SetMarkerStyle(8);

  myGridStyle->SetOptStat(kFALSE);
  myGridStyle->SetOptTitle(0); //set off the histogram title
  
  return myGridStyle;
}
 

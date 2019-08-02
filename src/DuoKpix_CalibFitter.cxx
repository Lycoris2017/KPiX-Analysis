//-----------------------------------------------------------------------------
// File          : calibrationFitter.cpp
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 05/30/2012
// Project       : Kpix Software Package
//-----------------------------------------------------------------------------
// Description :
// Application to process and fit kpix calibration data
//-----------------------------------------------------------------------------
// Copyright (c) 2012 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 05/30/2012: created
// 26/03/2018: modified for tracker sensor
//-----------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TStyle.h>
#include <stdarg.h>
#include <KpixEvent.h>
#include <KpixSample.h>
#include <Data.h>
#include <DataRead.h>
#include <math.h>
#include <fstream>
#include <XmlVariables.h>

#include "kpixmap.h"
#include "kpix_left_and_right.h"

using namespace std;

// Channel data
class ChannelData {
public:
  
  // Baseline Data
  uint         baseData[8192];
  uint         baseMin;
  uint         baseMax;
  double       baseCount;
  double       baseMean;
  double       baseSum;
  double       baseRms;
  
  // Baseline fit data
  double       baseFitMean;
  double       baseFitSigma;
  double       baseFitMeanErr;
  double       baseFitSigmaErr;
  double       baseFitChisquare;
  
  
  // Baseline hist data
  double       baseHistMean;
  double       baseHistRMS;
  
  // Calib Data
  double       calibCount[256];
  double       calibMean[256];
  double       calibSum[256];
  double       calibRms[256];
  double       calibError[256];
  double       calibOtherValue[1024];
  double       calibOtherDac[1024];
  
  ChannelData() {
    uint x;
    
    for (x=0; x < 8192; x++) baseData[x] = 0;
    baseMin          = 8192;
    baseMax          = 0;
    baseCount        = 0;
    baseMean         = 0;
    baseSum          = 0;
    baseRms          = 0;
    baseFitMean      = 0;
    baseFitSigma     = 0;
    baseFitMeanErr   = 0;
    baseFitSigmaErr  = 0;
    baseFitChisquare = 0;
    
    for (x=0; x < 256; x++) {
      calibCount[x]  = 0;
      calibMean[x]   = 0;
      calibSum[x]    = 0;
      calibRms[x]    = 0;
      calibError[x]  = 0;
    }
    for (x=0; x < 1024; x++) {
      calibOtherValue[x] = 0;
      calibOtherDac[x] = 0;
    }
  }
  
  void addBasePoint(uint data) {
    baseData[data]++;
    if ( data < baseMin ) baseMin = data;
    if ( data > baseMax ) baseMax = data;
    baseCount++;
    
    double tmpM = baseMean;
    double value = data;
    
    baseMean += (value - tmpM) / baseCount;
    baseSum  += (value - tmpM) * (value - baseMean);
  }
  
  void addCalibPoint(uint x, uint y) {
    calibCount[x]++;
    
    double tmpM = calibMean[x];
    double value = y;
    
    calibMean[x] += (value - tmpM) / calibCount[x];
    calibSum[x]  += (value - tmpM) * (value - calibMean[x]);
  }
  
  void addNeighborPoint(uint chan, uint x, uint y) {
    if ( y > calibOtherValue[chan] ) {
      calibOtherValue[chan] = y;
      calibOtherDac[chan] = x;
    }
  }
  
  void computeBase () {
    if ( baseCount > 0 ) baseRms = sqrt(baseSum / baseCount);
  }
  
  void computeCalib(double chargeError) {
    uint   x;
    double tmp;
    
    for (x=0; x < 256; x++) {
      if ( calibCount[x] > 0 ) {
	calibRms[x] = sqrt(calibSum[x] / calibCount[x]);
	tmp = calibRms[x] / sqrt(calibCount[x]);
	calibError[x] = sqrt((tmp * tmp) + (chargeError * chargeError));
      }
    }
  }
};


// Function to compute calibration charge
double calibCharge ( uint dac, bool positive, bool highCalib ) {
  double volt;
  double charge;
  
  if ( dac >= 0xf6 ) volt = 2.5 - ((double)(0xff-dac))*50.0*0.0001;
  else volt =(double)dac * 100.0 * 0.0001;
  
  if ( positive ) charge = (2.5 - volt) * 200e-15;
  else charge = volt * 200e-15;
  
  if ( highCalib ) charge *= 22.0;
  
  return(charge);
}

void addDoubleToXml ( ofstream *xml, uint indent, string variable, Double_t value ) {
  uint x;
  
  if ( ! isnan(value) ) {
    for (x=0; x < indent; x++){
      *xml << " ";
      *xml << "<" << variable << ">";
      *xml << value;
      *xml << "</" << variable << ">";
      *xml << endl;
    }
  }
}
  
void addStringToXml ( ofstream *xml, uint indent, string variable, string value ) {
  uint x;
  
  for (x=0; x < indent; x++) *xml << " ";
  *xml << "<" << variable << ">";
  *xml << value;
  *xml << "</" << variable << ">";
  *xml << endl;
}


// Process the data
int main ( int argc, char **argv ) {
  DataRead               dataRead;
  off_t                  fileSize;
  off_t                  filePos;
  KpixEvent              event;
  KpixSample             *sample;
  string                 calState;
  uint                   calChannel;
  uint                   calDac;
  uint                   lastPct;
  uint                   currPct;
  bool                   chanFound[32][1024];
  ChannelData            *chanData[32][1024][4][2];
  bool                   badMean[32][1024];
  bool                   badGain[32][1024];
  bool                   kpixFound[32];
  uint                   kpixMax;
  uint                   minDac;
  uint                   minChan;
  uint                   maxChan;
  uint                   x;
  uint                   range;
  uint                   value;
  uint                   kpix;
  uint                   channel;
  uint                   bucket;
  uint                   tstamp;
  string                 serial;
  KpixSample::SampleType type;
  TH1F                   *hist;
  //TH1F                   *resid[256];
  stringstream           tmp;
 // bool                   defresid=true;
  ofstream               xml;
  ofstream               csv;
  double                 grX[256];
  double                 grY[256];
  double                 grYErr[256];
  double                 grXErr[256];
  double                 grRes[256];
  uint                   grCount;
  TGraphErrors           *grCalib;
  TF1                    *fitCalib;
  TGraph                 *grResid;
  bool                   positive;
  bool                   b0CalibHigh;
  uint                   injectTime[5];
  uint                   eventCount;
  uint                   eventProcessed; // eventCount - skipped_cycles
  uint                   Baseline_eventCount;
  uint                   Inject_eventCount;
  
  uint                   skip_cycles_front;
  FILE*                  f_skipped_cycles;
  string                 outtxt;
  
  string                 outRoot;
  string                 outXml;
  string                 outCsv;
  TFile                  *rFile;
  TCanvas                *c1;
  char                   tstr[200];
  struct tm              *timeinfo;
  time_t                 tme;
  uint                   crChan;
  stringstream           crossString;
  stringstream           crossStringCsv;
  double                 crossDiff;
  uint                   badValue;
  XmlVariables           config;
  bool                   findBadMeanHist;
  bool                   findBadMeanFit;
  bool                   findBadMeanChisq;
  bool                   findBadGainFit;
  bool                   findBadGainChisq;
  double                 meanMin[2];
  double                 meanMax[2];
  double                 meanChisq;
  double                 gainMin[2];
  double                 gainMax[2];
  double                 chargeError[2];
  double                 gainChisq;
  double                 fitMin[2];
  double                 fitMax[2];
  double                 chisqNdf;
  ofstream               debug;
  ofstream				channel_file_bad;
  ofstream				channel_file_bad_fit;
  ofstream				channel_file_noise;
  ofstream				channel_file_calib;
  ofstream 				channel_file_adc_mean;
  uint                   badTimes;
  uint                   goodTimes;
  uint                   badMeanFitCnt;
  uint                   badMeanHistCnt;
  uint                   badMeanChisqCnt;
  uint                   badGainFitCnt;
  uint                   badGainChisqCnt;
  uint                   failedGainFit;
  uint                   failedMeanFit;
  uint                   badChannelCnt;
  uint					noiseSigmaCnt;
  uint                    errorSigmaCnt;
  stringstream			 FolderName;
  
  uint 					noise_cut = 1.0;

  uint nchann_calib=0;
  
 // bool                    printalot=false;
  
  unordered_map<uint, uint> kpix2strip_left;
  unordered_map<uint, uint> kpix2strip_right;
  kpix2strip_left = kpix_left();
  kpix2strip_right = kpix_right();
  
  // Init structure
  for (kpix=0; kpix < 32; kpix++) {
    for (channel=0; channel < 1024; channel++) {
      for (bucket=0; bucket < 4; bucket++) {
	chanData[kpix][channel][bucket][0] = NULL;
	chanData[kpix][channel][bucket][1] = NULL;
      }
      chanFound[kpix][channel] = false;
      badGain[kpix][channel] = false;
      badMean[kpix][channel] = false;
    }
    kpixFound[kpix] = false;
  }
  
  // Data file is the first and only arg
  if ( argc < 3 && argc >5 ) {
    cout << "Usage: new_calibrationFitter config_file data_file [skip_cycles_front] [debug]\n";
    return(1);
  }
  
   if ( argc >= 4 ){
    skip_cycles_front = atoi( argv[3] );
    cout<< " -- I am skipping first events: " << skip_cycles_front << endl;
    tmp.str("");
    tmp << argv[2] << ".printSkipped.txt";
    outtxt = tmp.str();
    f_skipped_cycles = fopen(outtxt.c_str(), "w");
   }
   else skip_cycles_front = 0;
   
   if ( argc == 5 ) debug.open(argv[4],ios::out | ios::trunc);

  
  // Read configuration
  if ( ! config.parseFile("config",argv[1]) ) {
    cout << "Failed to read configuration from " << argv[1] << endl;
    return(1);
  }

  // Extract configuration values
  findBadMeanHist  = config.getInt("FindBadMeanHist");
  findBadMeanFit   = config.getInt("FindBadMeanFit");
  meanMin[0]       = config.getDouble("GoodMeanMinR0");
  meanMax[0]       = config.getDouble("GoodMeanMaxR0");
  meanMin[1]       = config.getDouble("GoodMeanMinR1");
  meanMax[1]       = config.getDouble("GoodMeanMaxR1");
  findBadMeanChisq = config.getInt("FindBadMeanChisq");
  meanChisq        = config.getInt("GoodMeanChisqMax");
  findBadGainFit   = config.getInt("FindBadGainFit");
  gainMin[0]       = config.getDouble("GoodGainMinR0");
  gainMax[0]       = config.getDouble("GoodGainMaxR0");
  gainMin[1]       = config.getDouble("GoodGainMinR1");
  gainMax[1]       = config.getDouble("GoodGainMaxR1");
  findBadGainChisq = config.getInt("FindBadGainChisq");
  gainChisq        = config.getInt("GoodGainChisqMax");
  fitMin[0]        = config.getDouble("GainFitMinR0");
  fitMax[0]        = config.getDouble("GainFitMaxR0");
  fitMin[1]        = config.getDouble("GainFitMinR1");
  fitMax[1]        = config.getDouble("GainFitMaxR1");
  chargeError[0]   = config.getDouble("GainChargeErrorR0");
  chargeError[1]   = config.getDouble("GainChargeErrorR1");

  // Init a customized pol1, fit range will be re-range in fit()
  fitCalib = new TF1("fitCalib", "pol1",fitMin[0],fitMax[0] );
  fitCalib -> FixParameter( 0, 0 ); // offset to 0
  //fitCalib -> SetParameter( 0, 0 );
  fitCalib -> SetParameter( 1, 15); // slope
 

  // Open data file
  if ( ! dataRead.open(argv[2]) ) {
    cout << "Error opening data file " << argv[2] << endl;
    return(1);
  }
  
  // Create output names
  tmp.str("");
  if (skip_cycles_front==0){
    tmp << argv[2] << ".duoCalib.root";
    outRoot = tmp.str();
    
    tmp.str("");
    tmp << argv[2] << ".duoCalib.xml";
    outXml = tmp.str();
    
    tmp.str("");
    tmp << argv[2] << ".duoCalib.csv";
    outCsv = tmp.str();
  }
  else{
    tmp << argv[2] << ".skip" << skip_cycles_front <<".duoCalib.root";
    outRoot = tmp.str();

    tmp.str("");
    tmp << argv[2] << ".skip" << skip_cycles_front << ".duoCalib.xml";
    outXml = tmp.str();
    tmp.str("");
    tmp << argv[2] << ".skip" << skip_cycles_front << ".duoCalib.csv";
    outCsv = tmp.str();
  }
  
  
  //////////////////////////////////////////
  // Read Data
  //////////////////////////////////////////
  cout << "Opened data file: " << argv[2] << endl;
  fileSize = dataRead.size();
  filePos  = dataRead.pos();
  
  // Init
  currPct          	= 0;
  lastPct          	= 100;
  eventCount       	= 0;
  eventProcessed        = 0;
  Inject_eventCount     = 0;
  Baseline_eventCount   = 0;
  minChan          	= 0;
  maxChan          	= 0;
  badTimes         	= 0;
  badMeanFitCnt    	= 0;
  badMeanHistCnt   	= 0;
  badMeanChisqCnt  	= 0;
  badGainFitCnt    	= 0;
  badGainChisqCnt 	= 0;
  failedGainFit   	= 0;
  failedMeanFit   	= 0;
  badChannelCnt   	= 0;
  noiseSigmaCnt		= 0;
  errorSigmaCnt		= 0;
  cout << "\rReading File: 0 %" << flush;
  goodTimes       	= 0;

  // start - work in progress - wmq - Apr 11 2018
  /* KpixEvent etest;
  for (int it = 0; it<10; it++){
    dataRead.next(&etest);
    auto byte = etest.count();
    auto itrn = etest.eventNumber();
    printf(" index = %d , byte = %6d, train = %6d\n ", it, byte, itrn);

  }
  dataRead.close();
  dataRead.open(argv[2]);*/
  // end - work in progress - wmq - Apr 11 2018
  
  // Process each event
  while ( dataRead.next(&event) ) {
    
    
    // Get calibration state
    calState   = dataRead.getStatus("CalState");
    calChannel = dataRead.getStatusInt("CalChannel");
    calDac     = dataRead.getStatusInt("CalDac");

    // Get injection times
    /* 
     * 2018-04-11 by mengqing
     * - to check: if skipping event, should we get injection times from 1) #0 event 2) or 1st event to use?
     * - now using info from #0 evt;
     * - optional: if you want the alternative: change eventCount to eventProcessed
     */
    if ( eventCount == 0 ) { 
      minDac        = dataRead.getConfigInt("CalDacMin");
      minChan       = dataRead.getConfigInt("CalChanMin");
      maxChan       = dataRead.getConfigInt("CalChanMax");
      injectTime[0] = dataRead.getConfigInt("cntrlFpga:kpixAsic:Cal0Delay");
      injectTime[1] = dataRead.getConfigInt("cntrlFpga:kpixAsic:Cal1Delay") + injectTime[0] + 4;
      injectTime[2] = dataRead.getConfigInt("cntrlFpga:kpixAsic:Cal2Delay") + injectTime[1] + 4;
      injectTime[3] = dataRead.getConfigInt("cntrlFpga:kpixAsic:Cal3Delay") + injectTime[2] + 4;
      injectTime[4] = 8192;
    }

    if ( eventCount >= skip_cycles_front){
    
    // get each sample
    for (x=0; x < event.count(); x++) {
      
      // Get sample
      sample  = event.sample(x);
      kpix    = sample->getKpixAddress();
      channel = sample->getKpixChannel();
      bucket  = sample->getKpixBucket();
      value   = sample->getSampleValue();
      type    = sample->getSampleType();
      tstamp  = sample->getSampleTime();
      range   = sample->getSampleRange();
      
      // Only process real samples in the expected range
      if ( type == KpixSample::Data ) {
	
	// Create entry if it does not exist
	kpixFound[kpix]          = true;
	chanFound[kpix][channel] = true;
	if ( chanData[kpix][channel][bucket][range] == NULL ) chanData[kpix][channel][bucket][range] = new ChannelData;
	
	
	// Non calibration based run. Fill mean, ignore times
	if ( calState == "Idle" ) chanData[kpix][channel][bucket][range]->addBasePoint(value);
	
	// Filter for time
	else if ( tstamp > injectTime[bucket] && tstamp < injectTime[bucket+1] ) {
	  goodTimes++;
	  //cout << "Timestamp = " << tstamp << endl;
	  //cout << "Inject time of Bucket " << bucket << " = " << injectTime[bucket] << endl;
	  //cout << "Inject time of Bucket " << bucket+1 << " = " << injectTime[bucket+1] << endl << endl;
	  // Baseline
	  if ( calState == "Baseline" ) {
	    chanData[kpix][channel][bucket][range]->addBasePoint(value);
	  }
	  
	  // Injection
	  else if ( calState == "Inject" && calDac != minDac ) {
	    if ( channel == calChannel ) chanData[kpix][channel][bucket][range]->addCalibPoint(calDac, value);
	    else{
	      if ( chanData[kpix][calChannel][bucket][range] != NULL )
		chanData[kpix][calChannel][bucket][range]->addNeighborPoint(channel, calDac, value);
	      else cout<< "\n [Warn] Super Weird to check channel = "<< channel << endl;
	    }
	  }
	}
	else {
	  badTimes++;
	  //cout << "Timestamp = " << tstamp << endl;
	  //cout << "Inject time of Bucket " << bucket << " = " << injectTime[bucket] << endl;
	  //cout << "Inject time of Bucket " << bucket+1 << " = " << injectTime[bucket+1] << endl << endl;
	}
	
      }
    }
    eventProcessed++;
    } // skip cycle ends
    else {
      auto byte = event.count();
      auto train = event.eventNumber();
      if (f_skipped_cycles!=NULL)
	fprintf(f_skipped_cycles, " index = %d , byte = %6d, train = %6d, CalState = %s\n ", eventCount, byte, train, calState.c_str());
      
    }
    
    // Show progress
    filePos  = dataRead.pos();
    currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
    if ( currPct != lastPct ) {
      cout << "\rReading File: " << currPct << " %      " << flush;
      lastPct = currPct;
    }
    
    eventCount++;
    if (calState=="Baseline") 	    Baseline_eventCount++;
    if (calState == "Inject")       Inject_eventCount++;
    
  }
  cout << "\rReading File: Done.               " << endl;
  
  if (f_skipped_cycles!=NULL)  {
    fclose( f_skipped_cycles);
    cout << endl;
    cout << "Wrote skipped cycles to " << outtxt << endl;
    cout << endl;
  }
  //////////////////////////////////////////
  // Process Data
  //////////////////////////////////////////
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(111111111);
  
  // Default canvas
  c1 = new TCanvas("c1","c1");
	
  // Open root file
  rFile = new TFile(outRoot.c_str(),"recreate");
  
  // Open xml file
  xml.open(outXml.c_str(),ios::out | ios::trunc);
  xml << "<calibrationData>" << endl;
  
  // Open csv file
  csv.open(outCsv.c_str(),ios::out | ios::trunc);
  csv << "Kpix,Channel,Bucket,Range,BaseMean,BaseRms,BaseFitMean,BaseFitSigma";
  csv << ",BaseFitMeanErr,BaseFitSigmaErr,BaseFitChisquare";
  csv << ",CalibGain,CalibIntercept,CalibGainErr,CalibInterceptErr,CalibChisquare";
  csv << ",CalibGainRms,CrossTalkChan0,CrossTalkAmp0";
  csv << ",CrossTalkChan1,CrossTalkAmp1";
  csv << ",CrossTalkChan2,CrossTalkAmp2";
  csv << ",CrossTalkChan3,CrossTalkAmp3";
  csv << ",CrossTalkChan4,CrossTalkAmp4";
  csv << ",CrossTalkChan5,CrossTalkAmp5";
  csv << ",CrossTalkChan6,CrossTalkAmp6";
  csv << ",CrossTalkChan7,CrossTalkAmp7";
  csv << ",CrossTalkChan8,CrossTalkAmp8";
  csv << ",CrossTalkChan9,CrossTalkAmp9";
  csv << endl;
  
  // Open noise/bad channel list files
  channel_file_bad.open ("~/channel_list_bad.txt");
  channel_file_bad_fit.open ("~/channel_list_bad_fit.txt");
  channel_file_noise.open ("./channel_list_noise.txt");
  channel_file_calib.open("~/channel_list_calib.txt");
  channel_file_adc_mean.open("~/channel_adc_mean.txt");
  // Add notes
  xml << "   <sourceFile>" << argv[2] << "</sourceFile>" << endl;
  xml << "   <user>" <<  getlogin() << "</user>" << endl;

  time(&tme);
  timeinfo = localtime(&tme);
  strftime(tstr,200,"%Y_%m_%d_%H_%M_%S",timeinfo);
  xml << "   <timestamp>" << tstr << "</timestamp>" << endl;
  xml << "   <config>"<< endl;
  xml << config.getXml();
  xml << "   </config>"<< endl;
  
   
  
  // get calibration mode variables for charge computation
  positive    = (dataRead.getConfig("cntrlFpga:kpixAsic:CntrlPolarity") == "Positive");
  b0CalibHigh = (dataRead.getConfig("cntrlFpga:kpixAsic:CntrlCalibHigh") == "True");
  
  // Kpix count;
  for (kpix=0; kpix<32; kpix++) if ( kpixFound[kpix] ) kpixMax=kpix;
  
  //////////////////////////////////////////
  // Process Baselines 
  //////////////////////////////////////////
  
  // Process each kpix device
  
  rFile->cd(); // move into root folder base
  FolderName.str("");
  FolderName << "General";
  rFile->mkdir(FolderName.str().c_str()); // produce a sub folder with name of variable FolderName
  TDirectory *General_folder = rFile->GetDirectory(FolderName.str().c_str()); // get path to subdirectory
  General_folder->cd(); // move into subdirectory
  
  TH1F *pedestals = new TH1F("pedestals", "Pedestals distribution; Charge [ADC]; #entries", 9000, 0, 9000);
  TH1F *slope_hist = new TH1F("slope", "Slope distribution; Slope [ADC/fC]; #entries", 200, -20, 20);

  TH1F *slope_hist_k26 = new TH1F("slope_k26", "Slope distribution; Slope [ADC/fC]; #entries", 200, -20, 20);
  TH1F *slope_hist_k28 = new TH1F("slope_k28", "Slope distribution; Slope [ADC/fC]; #entries", 200, -20, 20);
  TH1F *slope_hist_k30 = new TH1F("slope_k30", "Slope distribution; Slope [ADC/fC]; #entries", 200, -20, 20);
   
  TH1F *slope_residual = new TH1F("slope_residual", "Slope_residual; Slope_residual[ADC/fC]; #entries", 2000, -100, 100);
  TH1F *pedestals_fc = new TH1F("pedestals_fc", "Pedestals distribution; Charge [fC]; #entries", 1000, -100, 100);
  
  TH1F *pedestalsRMS = new TH1F("pedestalsRMS", "Pedestal RMS; Charge [ADC]; #entries", 200, 0, 20);
  //  TH1F *slopeRMS = new TH1F("slopeRMS", "Slope RMS; Slope [ADC/fC]; #entries", 1000, 0, 20);
  //  TH1F *slope_residualRMS = new TH1F("slope_residualRMS", "Slope_residualRMS; Slope_residual[ADC/fC]; #entries", 1000, 0, 100);
  
  TH1F *slope_vs_channel_26 = new TH1F("slope_vs_channel_26", "; Channel ID; Slope [ADC/fC]", 1024, -0.5, 1023.5);
  TH1F *slope_vs_channel_28 = new TH1F("slope_vs_channel_28", "; Channel ID; Slope [ADC/fC]", 1024, -0.5, 1023.5);
  TH1F *slope_vs_channel_30 = new TH1F("slope_vs_channel_30", "; Channel ID; Slope [ADC/fC]", 1024, -0.5, 1023.5);
  
  TH1F *slope_vs_strip = new TH1F("slope_vs_strip", "; Strip_number; Slope [ADC/fC]", 1840, -0.5, 1839.5);
  TH1F *slope_vs_strip_26 = new TH1F("slope_vs_strip_26", "; Strip_number; Slope [ADC/fC]", 920, -0.5, 919.5);
  TH1F *slope_vs_strip_28 = new TH1F("slope_vs_strip_28", "; Strip_number; Slope [ADC/fC]", 920, 919.5, 1839.5);
    
  TH1F *RMSfc_v_channel = new TH1F("RMSfc_vs_channel", "; Channel_ID; RMS [fC]", 1024, -0.5, 1023.5);
  TH1F *RMSfc_v_channel_26 = new TH1F("RMSfc_vs_channel_26", "; Channel_ID; RMS [fC]", 1024, -0.5, 1023.5);
  TH1F *RMSfc_v_channel_28 = new TH1F("RMSfc_vs_channel_28", "; Channel_ID; RMS [fC]", 1024, -0.5, 1023.5);
  TH1F *RMSfc_v_channel_30 = new TH1F("RMSfc_vs_channel_30", "; Channel_ID; RMS [fC]", 1024, -0.5, 1023.5);
      
  TH1F *RMSfc_v_channel_connected = new TH1F("RMSfc_vs_channel_connected", "; Channel_ID; RMS [fC]", 1024, -0.5, 1023.5);
  
  TH1F *RMSfc_v_strip_left = new TH1F("RMSfc_vs_strip_left", "; Strip_number; RMS [fC]", 920, -0.5, 919.5);
  TH1F *RMSfc_v_strip_right = new TH1F("RMSfc_vs_strip_right", "; Strip_number; RMS [fC]", 920, 919.5, 1839.5);
  
  TH1F *PedestalsRMS_fc0 = new TH1F("pedestalsRMS_fc0", "Pedestals RMS, All Chn; [fC]; a.u.", 1000, 0, 4);
  TH1F *PedestalsRMS_fc0_26 = new TH1F("pedestalsRMS_fc0_26", "Pedestals RMS, All Chn; [fC]; a.u.", 1000, 0, 4);
  TH1F *PedestalsRMS_fc0_28 = new TH1F("pedestalsRMS_fc0_28", "Pedestals RMS, All Chn; [fC]; a.u.", 1000, 0, 4);
  TH1F *PedestalsRMS_fc0_30 = new TH1F("pedestalsRMS_fc0_30", "Pedestals RMS, All Chn; [fC]; a.u.", 1000, 0, 4);
  
  TH1F *PedestalsRMS_fc0_disc = new TH1F("pedestalsRMS_fc0_disc", "Pedestals RMS, disc. Chn; [fC]; a.u.", 1000, 0, 4);
  TH1F *PedestalsRMS_fc0_disc_26 = new TH1F("pedestalsRMS_fc0_disc_26", "Pedestals RMS, disc. Chn; [fC]; a.u.", 1000, 0, 4);
  TH1F *PedestalsRMS_fc0_disc_28 = new TH1F("pedestalsRMS_fc0_disc_28", "Pedestals RMS, disc. Chn; [fC]; a.u.", 1000, 0, 4);

  TH1F *PedestalsRMS_fc0_conn = new TH1F("pedestalsRMS_fc0_conn", "Pedestals RMS, conn. Chn; [fC]; a.u.", 1000, 0, 4);
  TH1F *PedestalsRMS_fc0_conn_26 = new TH1F("pedestalsRMS_fc0_conn_26", "Pedestals RMS, conn. Chn; [fC]; a.u.", 1000, 0, 4);
  TH1F *PedestalsRMS_fc0_conn_28 = new TH1F("pedestalsRMS_fc0_conn_28", "Pedestals RMS, conn. Chn; [fC]; a.u.", 1000, 0, 4);
    
  TH1F *PedestalsRMS_fc0_conn_120 = new TH1F("pedestalsRMS_fc0_conn_120", "Pedestals RMS, conn. Chn in 0-128; [fC]; a.u.", 1000, 0, 4);
  TH1F *PedestalsRMS_fc0_conn_800 = new TH1F("pedestalsRMS_fc0_conn_800", "Pedestals RMS, conn. Chn after 128; [fC]; a.u.", 1000, 0, 4);
  
  //TH1F *PedestalsRMS_fc2 = new TH1F("PedestalsRMS_fc2", "Pedestals RMS [fC]", 1000, 0, 10);
  //TH1F *PedestalsRMS_fc3 = new TH1F("PedestalsRMS_fc3", "Pedestals RMS [fC]", 1000, 0, 10);
  
  //TH1F *Pedestals_fc2 = new TH1F("Pedestals_fc2", "Pedestals [fC]", 1000, -100, 100);
  //TH1F *Pedestals_fc3 = new TH1F("Pedestals_fc3", "Pedestals [fC]", 1000, -100, 100);
  
  
  rFile->cd(); // move into root folder base
  FolderName.str("");
  FolderName << "Pedestals";
  rFile->mkdir(FolderName.str().c_str()); // produce a sub folder with name of variable FolderName
  TDirectory *pedestal_folder = rFile->GetDirectory(FolderName.str().c_str()); // get path to subdirectory
  pedestal_folder->cd(); // move into subdirectory
  
  for (kpix=0; kpix<32; kpix++) {
    if ( kpixFound[kpix] ) {
      
      // Get serial number
      tmp.str("");
      tmp << "cntrlFpga(0):kpixAsic(" << dec << kpix << "):SerialNumber";
      serial = dataRead.getConfig(tmp.str());
      
      // Process each channel
      for (channel=minChan; channel <= maxChan; channel++) {
	
	// Show progress
	cout << "\rProcessing baseline kpix " << dec << kpix << " / " << dec << kpixMax
	     << ", Channel " << channel << " / " << dec << maxChan
	     << "                 " << flush;
	
	// Channel is valid
	if ( chanFound[kpix][channel] ) {
	
	  // Each bucket
	  for (bucket = 0; bucket < 4; bucket++) {                   
	    
	    // Bucket is valid
	    if ( chanData[kpix][channel][bucket][0] != NULL || chanData[kpix][channel][bucket][1] != NULL ) {
	      
	      // Each range
	      for (range = 0; range < 2; range++) {
		
		// Range is valid
		if ( chanData[kpix][channel][bucket][range] != NULL ) {
		  chanData[kpix][channel][bucket][range]->computeBase();
		  
		  // Create histogram
		  tmp.str("");
		  tmp << "hist_" << serial << "_c" << dec << setw(4) << setfill('0') << channel;
		  tmp << "_b" << dec << bucket;
		  tmp << "_r" << dec << range;
		  
		  hist = new TH1F(tmp.str().c_str(),
				  (tmp.str()+";ADC; Entries / N_of_Cycles").c_str(),
				  8192,0,8192);
		  
		  //double num_of_entries; //ADDED
		  //double normed_bin_content;
		  // Fill histogram
		  //for (x=0; x < 8192; x++) num_of_entries += chanData[kpix][channel][bucket][range]->baseData[x]; //ADDED
		  //cout << num_of_entries << endl;
		  for (x=0; x < 8192; x++) {
		    hist->SetBinContent(x+1, double(chanData[kpix][channel][bucket][range]->baseData[x]));
		  }

		  //hist->Scale(1/double(eventCount)); // old version, scaled to all cycles;
		  hist->Scale(1/hist->Integral()); // normlized to 1;
		  
		  hist->GetXaxis()->SetRangeUser(chanData[kpix][channel][bucket][range]->baseMin,
						 chanData[kpix][channel][bucket][range]->baseMax);
		  hist->Fit("gaus","q");
		  // Double_t gausMin, gausMax;
		  // gausMin=hist->GetMean(1)-2*hist->GetRMS(1);
		  // gausMax=hist->GetMean(1)+2*hist->GetRMS(1);
		  // hist->Fit("gaus","q","",gausMin, gausMax);
		  
		  chanData[kpix][channel][bucket][range]->baseHistMean       = hist->GetMean();
		  chanData[kpix][channel][bucket][range]->baseHistRMS        = hist->GetRMS();
		  
		  if ( hist->GetFunction("gaus") ) {
		    chanData[kpix][channel][bucket][range]->baseFitMean      = hist->GetFunction("gaus")->GetParameter(1);
		    chanData[kpix][channel][bucket][range]->baseFitSigma     = hist->GetFunction("gaus")->GetParameter(2);
		    chanData[kpix][channel][bucket][range]->baseFitMeanErr   = hist->GetFunction("gaus")->GetParError(1);
		    chanData[kpix][channel][bucket][range]->baseFitSigmaErr  = hist->GetFunction("gaus")->GetParError(2);
		    
		    pedestals->Fill(chanData[kpix][channel][bucket][range]->baseFitMean);
		    channel_file_adc_mean << chanData[kpix][channel][bucket][range]->baseFitMean << " " << channel << " " << bucket << endl;
		    pedestalsRMS->Fill(chanData[kpix][channel][bucket][range]->baseFitSigma);
		    if ( hist->GetFunction("gaus")->GetNDF() == 0 ) {
		      chanData[kpix][channel][bucket][range]->baseFitChisquare = 0;
		    } else {
		      chanData[kpix][channel][bucket][range]->baseFitChisquare = 
			(hist->GetFunction("gaus")->GetChisquare() / hist->GetFunction("gaus")->GetNDF() );
		    }
		    
		    // Determine bad channel from fitted chisq
		    
		    if ( findBadMeanChisq && (chanData[kpix][channel][bucket][range]->baseFitChisquare >  meanChisq) ) {
		      debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
			    << " Range=" << dec << range 
			    << " Bad fit mean chisq=" << chanData[kpix][channel][bucket][range]->baseFitChisquare << endl;
		      badMean[kpix][channel] = true;
		      badMeanChisqCnt++;
		    }
		    
		    // Determine bad channel from fitted mean
		    if ( findBadMeanFit &&
			 ( (chanData[kpix][channel][bucket][range]->baseFitMean > meanMax[range]) ||  
			   (chanData[kpix][channel][bucket][range]->baseFitMean < meanMin[range]) ) ) {
		      debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
			    << " Range=" << dec << range 
			    << " Bad fit mean value=" << chanData[kpix][channel][bucket][range]->baseFitMean << endl;
		      badMean[kpix][channel] = true;
		      badMeanFitCnt++;
		    }
		    // Determine noisy channels from fitted sigma
		    if  ((chanData[kpix][channel][bucket][range]->baseFitSigma/chanData[kpix][channel][bucket][range]->baseFitMean) > 0.03) {
		      badMean[kpix][channel] = true;
		      noiseSigmaCnt++;
		      //cout << endl << "Noisy channel with " << "sigma = " << chanData[kpix][channel][bucket][range]->baseFitSigma << endl << endl;
		     // channel_file_noise << channel << endl;
		    }
		    
		    // Determine weird channels with only a single peak if fitted sigma is below the error of sigma
		    if (chanData[kpix][channel][bucket][range]->baseFitSigma < chanData[kpix][channel][bucket][range]->baseFitSigmaErr) {
		      badMean[kpix][channel] = true;
		      errorSigmaCnt++;
		      //cout << endl << "Sigma below error value with " << "sigma = " << chanData[kpix][channel][bucket][range]->baseFitSigma << endl << endl;
		     // channel_file_bad << channel << endl;
		    }
		  }
		  else if ( findBadMeanFit || findBadMeanChisq ) {
		    debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel 
			  << " Bucket=" << dec << bucket << " Range=" << dec << range
			  << " Failed to fit mean" << endl;
		    badMean[kpix][channel] = true;
		    failedMeanFit++;
		  }
		  
		  // Determine bad channel from histogram mean
		  if ( findBadMeanHist && 
		       ( (chanData[kpix][channel][bucket][range]->baseMean > meanMax[range]) ||
			 (chanData[kpix][channel][bucket][range]->baseMean < meanMin[range]) ) ) {
		    debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
			  << " Range=" << dec << range
			  << " Bad hist mean value=" << chanData[kpix][channel][bucket][range]->baseMean << endl;
		    badMeanHistCnt++;
		    badMean[kpix][channel] = true;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  cout << endl;
  
  //////////////////////////////////////////
  // Process Calibration
  //////////////////////////////////////////
  rFile->cd(); // move into root folder base
  FolderName.str("");
  FolderName << "Calibration_and_Residuals";
  rFile->mkdir(FolderName.str().c_str()); // produce a sub folder with name of variable FolderName
  TDirectory *calibration_folder = rFile->GetDirectory(FolderName.str().c_str()); // get path to subdirectory
  calibration_folder->cd(); // move into subdirectory
  // Process each kpix device
  for (kpix=0; kpix<32; kpix++){
    if ( kpixFound[kpix] ){ 
      
      // Get serial number
      tmp.str("");
      tmp << "cntrlFpga(0):kpixAsic(" << dec << kpix << "):SerialNumber";
      serial = dataRead.getConfig(tmp.str());
      xml << "   <kpixAsic id=\"" << serial << "\">" << endl;
      
      // Process each channel
      for (channel=minChan; channel <= maxChan; channel++) {
	
	// Show progress
	cout << "\rProcessing calibration kpix " << dec << kpix << " / " << dec << kpixMax
	     << ", Channel " << channel << " / " << dec << maxChan 
	     << "                 " << flush;
	
	// Channel is valid
	if ( chanFound[kpix][channel] ) {

	  // Start channel marker
	  xml << "      <Channel id=\"" << channel << "\">" << endl;
	  
	  // Each bucket
	  for (bucket = 0; bucket < 4; bucket++) {
	    
	    // Bucket is valid
	    if ( chanData[kpix][channel][bucket][0] != NULL || chanData[kpix][channel][bucket][1] != NULL ) {
	      xml << "         <Bucket id=\"" << bucket << "\">" << endl;
	      
	      // Each range
	      for (range = 0; range < 2; range++) {
		
		// Range is valid
		if ( chanData[kpix][channel][bucket][range] != NULL ) {
		  xml << "            <Range id=\"" << range << "\">" << endl;
		  chanData[kpix][channel][bucket][range]->computeCalib(chargeError[range]);
		  csv << serial << "," << dec << channel << "," << dec << bucket << "," << dec << range;
		  
		  // Add baseline data to xml
		  addDoubleToXml(&xml,15,"BaseMean",chanData[kpix][channel][bucket][range]->baseMean);
		  addDoubleToXml(&xml,15,"BaseRms",chanData[kpix][channel][bucket][range]->baseRms);
		  if ( chanData[kpix][channel][bucket][range]->baseFitMean != 0 ) {
		    addDoubleToXml(&xml,15,"BaseFitMean",chanData[kpix][channel][bucket][range]->baseFitMean);
		    addDoubleToXml(&xml,15,"BaseFitSigma",chanData[kpix][channel][bucket][range]->baseFitSigma);
		    addDoubleToXml(&xml,15,"BaseFitMeanErr",chanData[kpix][channel][bucket][range]->baseFitMeanErr);
		    addDoubleToXml(&xml,15,"BaseFitSigmaErr",chanData[kpix][channel][bucket][range]->baseFitSigmaErr);
		    addDoubleToXml(&xml,15,"BaseFitChisquare",chanData[kpix][channel][bucket][range]->baseFitChisquare);
		  }
		  
		  // Add baseline data to excel file
		  csv << "," << chanData[kpix][channel][bucket][range]->baseMean;
		  csv << "," << chanData[kpix][channel][bucket][range]->baseRms;
		  csv << "," << chanData[kpix][channel][bucket][range]->baseFitMean;
		  csv << "," << chanData[kpix][channel][bucket][range]->baseFitSigma;
		  csv << "," << chanData[kpix][channel][bucket][range]->baseFitMeanErr;
		  csv << "," << chanData[kpix][channel][bucket][range]->baseFitSigmaErr;
		  csv << "," << chanData[kpix][channel][bucket][range]->baseFitChisquare;
		  
		  // Create calibration graph
		  grCount = 0;
		  crossString.str("");
		  crossStringCsv.str("");
		  for (x=0; x < 256; x++) {
		    
		    // Calibration point is valid
		    if ( chanData[kpix][channel][bucket][range]->calibCount[x] > 0 ) {
		      grX[grCount]    = calibCharge ( x, positive, ((bucket==0)?b0CalibHigh:false));
		      grY[grCount]    = chanData[kpix][channel][bucket][range]->calibMean[x];
		      grYErr[grCount] = chanData[kpix][channel][bucket][range]->calibError[x];
		      grXErr[grCount] = 0;
		      debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
			    << " Range=" << dec << range
			    << " Adding point x=" << grX[grCount] 
			    << " Rms=" << chanData[kpix][channel][bucket][range]->calibRms[x]
			    << " Error=" << chanData[kpix][channel][bucket][range]->calibError[x] << endl;
		      grCount++;
		      
		      // Find crosstalk, value - base > 3 * sigma
		      for (crChan=0; crChan < 1024; crChan++ ) {
			
			if ( chanData[kpix][crChan][bucket][range] != NULL ) {
			  
			  crossDiff = chanData[kpix][channel][bucket][range]->calibOtherValue[crChan] - 
			    chanData[kpix][crChan][bucket][range]->baseMean;
			    
			  if ( (chanData[kpix][channel][bucket][range]->calibOtherDac[crChan] == x)  && 
			       (crChan != channel) && 
			       (chanData[kpix][channel][bucket][range] != NULL ) &&
			       (crossDiff > (10.0 * chanData[kpix][crChan][bucket][range]->baseRms))) 
			    {
			      
			      if ( crossString.str() != "" ) crossString << " ";
			      crossString << dec << crChan << ":" << dec << (uint)crossDiff;
			      crossStringCsv << "," << dec << crChan << "," << dec << (uint)crossDiff;
			    }
			}
		      }
		    }
		  }
		  
		  // Create graph
		  if ( grCount > 0 ) {
	    
		    grCalib = new TGraphErrors(grCount,grX,grY,grXErr,grYErr);
		    grCalib->Draw("Ap");
		    grCalib->GetXaxis()->SetTitle("Charge [C]");
		    grCalib->GetYaxis()->SetTitle("ADC");
		    grCalib->Fit("pol1","eq","",fitMin[range],fitMax[range]);
		    grCalib->GetFunction("pol1")->SetLineWidth(1);
		    
		    // Create name and write
		    tmp.str("");
		    tmp << "calib_" << serial << "_c" << dec << setw(4) << setfill('0') << channel;
		    tmp << "_b" << dec << bucket;
		    tmp << "_r" << dec << range;
		    //grCalib->SetTitle(tmp.str().c_str());
		    grCalib->Write(tmp.str().c_str());
		    nchann_calib++;
		    
		    // Create and store residual plot
		    for (x=0; x < grCount; x++) grRes[x] = (grY[x] - grCalib->GetFunction("pol1")->Eval(grX[x]));
		    grResid = new TGraph(grCount,grX,grRes);
		    grResid->Draw("Ap");
		    grResid->GetXaxis()->SetTitle("Charge [C]");
		    grResid->GetYaxis()->SetTitle("ADC");
		    
		    // Create name and write
		    tmp.str("");
		    tmp << "resid_" << serial << "_c" << dec << setw(4) << setfill('0') << channel;
		    tmp << "_b" << dec << bucket;
		    tmp << "_r" << dec << range;
		    //grResid->SetTitle(tmp.str().c_str());
		    grResid->Write(tmp.str().c_str());
		    
		    // Add to xml
		    if ( grCalib->GetFunction("pol1") ) {
		      chisqNdf = (grCalib->GetFunction("pol1")->GetChisquare() / grCalib->GetFunction("pol1")->GetNDF());
		      Double_t slope = grCalib->GetFunction("pol1")->GetParameter(1);
		      Double_t offset = offset;
		      
		      long double ped_charge = ( chanData[kpix][channel][bucket][range]->baseFitMean ) / slope;
		      long double ped_charge_err = (chanData[kpix][channel][bucket][range]->baseHistRMS) / slope ; // simple err
		      
		      //cout << "RMS in ADC " << chanData[kpix][channel][bucket][range]->baseHistRMS << endl;
		      //cout << "Slope " << slope << endl;
		      //cout << "RMS in fc " << ped_charge_err << endl;
		      
		      pedestals_fc->Fill( ped_charge * pow(10,15) );
		      
		      PedestalsRMS_fc0->Fill( ped_charge_err * pow(10,15) );
		      if (kpix==26) PedestalsRMS_fc0_26->Fill( ped_charge_err * pow(10,15) );
		      if (kpix==28) PedestalsRMS_fc0_28->Fill( ped_charge_err * pow(10,15) );
		      if (kpix==30) PedestalsRMS_fc0_30->Fill( ped_charge_err * pow(10,15) );
		      
		      //cout << ped_charge_err * pow(10,15) << endl;
		      if ( kpix2strip_left.at(channel) == 9999 ){
			PedestalsRMS_fc0_disc->Fill( ped_charge_err * pow(10,15) );

			if (kpix==26) PedestalsRMS_fc0_disc_26->Fill( ped_charge_err * pow(10,15) );
			if (kpix==28) PedestalsRMS_fc0_disc_28->Fill( ped_charge_err * pow(10,15) );
		      }
		      else{
			if (channel < 128) PedestalsRMS_fc0_conn_120 -> Fill(ped_charge_err * pow(10,15));
			else PedestalsRMS_fc0_conn_800 -> Fill(ped_charge_err * pow(10,15));
			
			PedestalsRMS_fc0_conn->Fill( ped_charge_err * pow(10,15) );

			if (kpix==26) PedestalsRMS_fc0_conn_26->Fill( ped_charge_err * pow(10,15) );
			if (kpix==28) PedestalsRMS_fc0_conn_28->Fill( ped_charge_err * pow(10,15) );
		      }
		      
		      if (ped_charge_err * pow(10,15) >= noise_cut){
			channel_file_noise << channel << endl ;
		      }
		      
		      slope_hist->Fill( slope / pow(10,15) );
		      if (kpix==26) {
				  slope_hist_k26->Fill(slope/pow(10,15));
				  slope_vs_channel_26->SetBinContent( channel+1, slope / pow(10,15));
				  if ( kpix2strip_left.at(channel)!=9999 ) {
				    slope_vs_strip_26->SetBinContent( kpix2strip_left.at(channel) + 1, slope / pow(10,15));
				    slope_vs_strip->SetBinContent( kpix2strip_left.at(channel) + 1, slope / pow(10,15) );
				  }
			  }
		      if (kpix==28) {
				  slope_hist_k28->Fill(slope/pow(10,15));
				  slope_vs_channel_28->SetBinContent( channel+1, slope / pow(10,15));
				  if ( kpix2strip_right.at(channel)!=9999 ){
				    slope_vs_strip_28->SetBinContent(  kpix2strip_right.at(channel)-919, slope / pow(10,15) );
				    slope_vs_strip->SetBinContent( kpix2strip_right.at(channel) + 1, slope / pow(10,15) ); 
				  }
		      } 

		      if (kpix==30){
			slope_hist_k30->Fill(slope/pow(10,15));
			slope_vs_channel_30->SetBinContent( channel+1, slope / pow(10,15) );
		      }
		      
		      if (abs(ped_charge_err * pow(10,15)) < 20) {
			RMSfc_v_channel->SetBinContent(channel, ped_charge_err * pow(10,15));
			if (kpix==26) RMSfc_v_channel_26->SetBinContent(channel, ped_charge_err * pow(10,15));
			if (kpix==28) RMSfc_v_channel_28->SetBinContent(channel, ped_charge_err * pow(10,15));
			if (kpix==30) RMSfc_v_channel_30->SetBinContent(channel, ped_charge_err * pow(10,15));
			
			if (kpix2strip_left.at(channel) != 9999) RMSfc_v_channel_connected->SetBinContent(channel, ped_charge_err * pow(10,15));
			RMSfc_v_strip_left->SetBinContent( kpix2strip_left.at(channel)+1, ped_charge_err * pow(10,15) );
			RMSfc_v_strip_right->SetBinContent( kpix2strip_right.at(channel)-919, ped_charge_err * pow(10,15) ); //as we set bin content need to move everything to the left by as channel 920 is bin 1, bin0 is underflow
			
			//cout << kpix2strip_right.at(channel) << endl;
		      }
		      //summary2_1->Fill( channel, slope / pow(10,15) );
		      
		      slope_residual->Fill( offset);
		      
		      //summary12->Fill( grCalib->GetFunction("pol1")->GetParError(1) / pow(10,15) );
		      //summary13->Fill( grCalib->GetFunction("pol1")->GetParError(0));
		      
		      addDoubleToXml(&xml,15,"CalibGain",grCalib->GetFunction("pol1")->GetParameter(1));
		      addDoubleToXml(&xml,15,"CalibIntercept",grCalib->GetFunction("pol1")->GetParameter(0));
		      addDoubleToXml(&xml,15,"CalibGainErr",grCalib->GetFunction("pol1")->GetParError(1));
		      addDoubleToXml(&xml,15,"CalibInterceptErr",grCalib->GetFunction("pol1")->GetParError(0));
		      addDoubleToXml(&xml,15,"CalibChisquare",chisqNdf);
		      csv << "," << grCalib->GetFunction("pol1")->GetParameter(1);
		      csv << "," << grCalib->GetFunction("pol1")->GetParameter(0);
		      csv << "," << grCalib->GetFunction("pol1")->GetParError(1);
		      csv << "," << grCalib->GetFunction("pol1")->GetParError(0);
		      csv << "," << chisqNdf;
		      
		      // Determine bad channel from fitted gain
		      if ( findBadGainFit && 
			   ( (grCalib->GetFunction("pol1")->GetParameter(1) > gainMax[range]) ||
			     (grCalib->GetFunction("pol1")->GetParameter(1) < gainMin[range]) ) ) {
			debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
			      << " Range=" << dec << range
			      << " Bad gain value=" << grCalib->GetFunction("pol1")->GetParameter(1) << endl;
			badGain[kpix][channel] = true;
			badGainFitCnt++;
		      }
		      
		      // Determine bad channel from fitted chisq
		      if ( findBadGainChisq && (chisqNdf >  gainChisq) ) {
			debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
			      << " Range=" << dec << range
			      << " Bad gain chisq=" << gainChisq << endl;
			badGain[kpix][channel] = true;
			badGainChisqCnt++;
		      }
		    }
		    else {
		      csv << ",0,0,0,0,0";
		      if ( findBadGainFit || findBadGainChisq )  {
			debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
			      << " Range=" << dec << range
			      << " Failed to fit gain" << endl;
			badGain[kpix][channel] = true;
			failedGainFit++;
		      }
		    }
		    
		    addDoubleToXml(&xml,15,"CalibGainRms",grCalib->GetRMS(2));
		    csv << "," << grCalib->GetRMS(2);
		    
		    if ( crossString.str() != "" ) addStringToXml(&xml,15,"CalibCrossTalk",crossString.str());
		    csv << crossStringCsv.str();
		  }
		  csv << endl; 
		  xml << "            </Range>" << endl;
		}
	      }
	      xml << "         </Bucket>" << endl;
	    }
	  }
	  
	  // Determine if the channel is bad
	  
	  badValue = 0;
	  if ( badMean[kpix][channel] ) badValue |= 0x1;
	  if ( badGain[kpix][channel] ) badValue |= 0x2;
	  
	  if ( badValue != 0 ) {
	    debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel
		  << " Marking channel bad." << endl;
	    badChannelCnt++;
	  }
	  
	  xml << "         <BadChannel>" << dec << badValue << "</BadChannel>" << endl;
	  xml << "      </Channel>" << endl;
	}
      }
         xml << "   </kpixAsic>" << endl;
    }
  }
  cout << endl;
  cout << "Wrote root plots to " << outRoot << endl;
  cout << "Wrote xml data to " << outXml << endl;
  cout << "Wrote csv data to " << outCsv << endl;
  cout << endl;
  
  cout << "Found " << dec << setw(10) << setfill(' ') << badTimes        << " events with bad times" << endl;
  cout << "Found " << dec << setw(10) << setfill(' ') << badMeanFitCnt   << " bad mean fit values" << endl;
  cout << "Found " << dec << setw(10) << setfill(' ') << badMeanChisqCnt << " bad mean fit chisq"  << endl;
  cout << "Found " << dec << setw(10) << setfill(' ') << badMeanHistCnt  << " bad mean hist values" << endl;
  cout << "Found " << dec << setw(10) << setfill(' ') << failedMeanFit   << " failed mean fits" << endl;
  cout << "Found " << dec << setw(10) << setfill(' ') << badGainFitCnt   << " bad gain fit values" << endl;
  cout << "Found " << dec << setw(10) << setfill(' ') << badGainChisqCnt << " bad gain fit chisq" << endl;
  cout << "Found " << dec << setw(10) << setfill(' ') << failedGainFit   << " failed gain fits" << endl;
  cout << "Found " << dec << setw(10) << setfill(' ') << badChannelCnt   << " bad channels" << endl;

  cout << "DEBUG : " << nchann_calib << endl ;
    
  xml << "</calibrationData>" << endl;
  xml.close();
  rFile->Write();
  delete rFile;
  
  // Cleanup
  for (kpix=0; kpix < 32; kpix++) {
    for (channel=0; channel < 1024; channel++) {
      for (bucket=0; bucket < 4; bucket++) {
	if ( chanData[kpix][channel][bucket][0] != NULL ) delete chanData[kpix][channel][bucket][0];
	if ( chanData[kpix][channel][bucket][1] != NULL ) delete chanData[kpix][channel][bucket][1];
      }
    }
  }
  
  // Close file
  dataRead.close();
  return(0);
}

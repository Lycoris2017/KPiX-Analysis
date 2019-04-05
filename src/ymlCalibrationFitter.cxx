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
#include "YmlVariables.h"


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
	
	if ( dac >= 0xf6 )
	{
			volt = 2.5 - ((double)(0xff-dac))*50.0*0.0001;
			//cout << "A " << volt << endl;
	}
	else
	{
		volt =(double)dac * 100.0 * 0.0001;
		//cout << "B " << volt << endl;
	}
	
	if ( positive )
	{
		 charge = (2.5 - volt) * 200e-15;
		 //cout << "C " << charge << endl;
	 }
	else
	{
		 charge = volt * 200e-15;
		//cout << "D " << charge << endl;
	}
	
	if ( highCalib ) charge *= 22.0;
	
	return(charge);
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
  bool                   chanFound[24][1024];
  ChannelData            *chanData[24][1024][4][2];
  bool                   badMean[24][1024];
  bool                   badGain[24][1024];
  bool                   kpixFound[24];
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
  TH1F                   *resid[256];
  stringstream           tmp;
 // bool                   defresid=true;
  double                 grX[256];
  double                 grDAC[256];
  double                 grY[256];
  double                 grYErr[256];
  double                 grXErr[256];
  double                 grRes[256];
  uint                   grCount;
  TGraphErrors           *grCalib;
  TGraphErrors           *grCalibDAC;
  TF1                    *fitCalib;
  TF1                    *fitCalibDAC;
  TGraph                 *grResid;
  bool                   positive[24];
  bool                   b0CalibHigh;
  uint                   injectTime[5];
  uint                   eventCount;
  uint                   eventProcessed; // eventCount - skipped_cycles
  uint                   Baseline_eventCount;
  uint                   Inject_eventCount;
  ofstream               logfile;
  uint                   skip_cycles_front;
  FILE*                  f_skipped_cycles;
  string                 outtxt;
  
  string                 outRoot;
  TFile                  *rFile;
  TCanvas                *c1;
  char                   tstr[200];
  struct tm              *timeinfo;
  time_t                 tme;
  uint                   crChan;
  double                 crossDiff;
  uint                   badValue;
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
  
	TH1F				*pedestals[24][4];
	TH1F				*pedestals_fc_0_127[24][4];
	TH1F				*slope_hist[24][4];
	TH1F				*slope_hist_conn[24][4];
	TH1F				*slope_hist_disc[24][4];
	TH1F				*slope_hist_0_127[24][4];
	TH1F				*slope_residual[24][4];
	TH1F				*pedestals_fc[24][4];
	TH1F				*pedestalsRMS_fc[24][4];
	TH1F				*pedestalsRMS_fc_0_127[24][4];
	TH1F				*pedestalsRMS_fc_disc[24][4];
	TH1F				*pedestalsRMS_fc_conn[24][4];
	TH1F				*slopeRMS[24][4];
	TH1F				*slope_vs_channel[24][4];
	TH1F				*slope_vs_right_strip[24][4];
	TH1F				*slope_vs_left_strip[24][4];
	TH1F				*RMSfC_vs_channel[24][4];
  
  
  
	uint 					noise_cut = 1.0;
	
	// bool                    printalot=false;
	
	unordered_map<uint, uint> kpix2strip_left;
	unordered_map<uint, uint> kpix2strip_right;
	kpix2strip_left = kpix_left();
	kpix2strip_right = kpix_right();
	
	// Init structure
	for (kpix=0; kpix < 24; kpix++) {
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
	
	
	//// Read configuration
	//if ( ! config.parseFile("config",argv[1]) ) {
		//cout << "Failed to read configuration from " << argv[1] << endl;
		//return(1);
	//}
	
	// Extract configuration values
	//findBadMeanHist  = config.getInt("FindBadMeanHist");
	//findBadMeanFit   = config.getInt("FindBadMeanFit");
	//meanMin[0]       = config.getDouble("GoodMeanMinR0");
	//meanMax[0]       = config.getDouble("GoodMeanMaxR0");
	//meanMin[1]       = config.getDouble("GoodMeanMinR1");
	//meanMax[1]       = config.getDouble("GoodMeanMaxR1");
	//findBadMeanChisq = config.getInt("FindBadMeanChisq");
	//meanChisq        = config.getInt("GoodMeanChisqMax");
	//findBadGainFit   = config.getInt("FindBadGainFit");
	//gainMin[0]       = config.getDouble("GoodGainMinR0");
	//gainMax[0]       = config.getDouble("GoodGainMaxR0");
	//gainMin[1]       = config.getDouble("GoodGainMinR1");
	//gainMax[1]       = config.getDouble("GoodGainMaxR1");
	//findBadGainChisq = config.getInt("FindBadGainChisq");
	//gainChisq        = config.getInt("GoodGainChisqMax");
	//fitMin[0]        = config.getDouble("GainFitMinR0");
	//fitMax[0]        = config.getDouble("GainFitMaxR0");
	//fitMin[1]        = config.getDouble("GainFitMinR1");
	//fitMax[1]        = config.getDouble("GainFitMaxR1");
	//chargeError[0]   = config.getDouble("GainChargeErrorR0");
	//chargeError[1]   = config.getDouble("GainChargeErrorR1");
	
	// Init a customized pol1, fit range will be re-range in fit()
	fitCalib = new TF1("fitCalib", "pol1",fitMin[0],fitMax[0] );
	fitCalib -> FixParameter( 0, 0 ); // offset to 0
	//fitCalib -> SetParameter( 0, 0 );
	fitCalib -> SetParameter( 1, 15); // slope
	
	
	////// Open data file +++++++++++++++++++++
	if ( argc != 2 ) {
		cout << "\nUsage: \ncalibrationFitter data_file \n\n";
		return(1);
	}
	
	// Open data file
	if ( ! dataRead.open(argv[1]) ) {
		cout << "Error opening data file " << argv[1] << endl;
		return(1);
	}
	
	// debug log:
	if (debug){
		cout << "open a debug log file"<< endl;
		logfile.open("debug.log", std::ofstream::out);
	}
	
	
	// Create output names
	tmp.str("");
	tmp << argv[1] << ".ymlcalib.root";
	outRoot = tmp.str();
		
	rFile = new TFile(outRoot.c_str(),"recreate");
	
	//////////////////////////////////////////
	// Read Data
	//////////////////////////////////////////
	//cout << "Opened data file: " << argv[2] << endl; ++++++++++++++++++++++++
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
	
	
	
	
	
	// Process each event +++++++++++++++++++++++
	while ( dataRead.next(&event) ) {
		
		
		if ( eventCount == 0 ) {
	
			minDac        = dataRead.getYmlStatusInt("CalDacMin");
			minChan       = dataRead.getYmlStatusInt("CalChanMin");
			maxChan       = dataRead.getYmlStatusInt("CalChanMax");
			
			for (int k = 0; k < 24; ++k && positive )
			{
				
				tmp.str("");
				tmp << "KpixDaqCore:KpixAsicArray:KpixAsic[" << k << "]:CntrlPolarity";
				positive[k] = (dataRead.getYmlConfig(tmp.str()) == "Positive");
			}
			
			
			logfile << "CalDacMin : " << minDac <<endl;
			logfile << "CalChanMin, Max : " << minChan << ", " << maxChan << endl;
			
			injectTime[0] = dataRead.getYmlConfigInt("KpixDaqCore:KpixAsicArray:KpixAsic[24]:Cal0Delay");
			injectTime[1] = dataRead.getYmlConfigInt("KpixDaqCore:KpixAsicArray:KpixAsic[24]:Cal1Delay") + injectTime[0] + 4;
			injectTime[2] = dataRead.getYmlConfigInt("KpixDaqCore:KpixAsicArray:KpixAsic[24]:Cal2Delay") + injectTime[1] + 4;
			injectTime[3] = dataRead.getYmlConfigInt("KpixDaqCore:KpixAsicArray:KpixAsic[24]:Cal3Delay") + injectTime[2] + 4;
			injectTime[4] = 8192;
			
			logfile << "CalDelay : [" << dec
				<< injectTime[0] << ", "	
				<< injectTime[1] << ", "
				<< injectTime[2] << ", "
				<< injectTime[3] << ", "
				<< injectTime[4] << "]" << endl;
		}
		
		// Get calibration state
		calState   = dataRead.getYmlStatus("CalState");
		calChannel = dataRead.getYmlStatusInt("CalChannel");
		calDac     = dataRead.getYmlStatusInt("CalDac");
	
		bool foundCalChannel = false;
		
		logfile << "CalState :" << calState << "\n"
			<< "CalChannel : " << calChannel << "\n"
			<< "CalDac :     " << calDac << "\n";
	
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
	
  
   
  
  // get calibration mode variables for charge computation
  //positive    = (dataRead.getConfig("cntrlFpga:kpixAsic:CntrlPolarity") == "Positive");
  b0CalibHigh = (dataRead.getConfig("cntrlFpga:kpixAsic:CntrlCalibHigh") == "True");

  
  // Kpix count;
  for (kpix=0; kpix<24; kpix++) if ( kpixFound[kpix] ) kpixMax=kpix;
  
  //////////////////////////////////////////
  // Process Baselines 
  //////////////////////////////////////////
  
  // Process each kpix device
  
	for (kpix = 0; kpix < 24; kpix++) //looping through all possible kpix
	{
		//
		//cout << "DEBUG test " << kpixFound[kpix] << endl;
		if (kpixFound[kpix]) //checking if kpix exists
		{
			rFile->cd(); //producing subfolder for kpix same as above for the event subfolder structure
			FolderName.str("");
			FolderName << "KPiX_" << kpix;
			rFile->mkdir(FolderName.str().c_str());
			TDirectory *kpix_folder = rFile->GetDirectory(FolderName.str().c_str());
			kpix_folder->cd();
	
		
			FolderName.str("");
			FolderName << "Buckets";
			kpix_folder->mkdir(FolderName.str().c_str());
			TDirectory *buckets_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
			rFile->cd(buckets_folder->GetPath());
			
			for (bucket = 0; bucket < 4; bucket++)
			{
				FolderName.str("");
				FolderName << "Bucket " << bucket;
				buckets_folder->mkdir(FolderName.str().c_str());
				TDirectory *bucket_folder = buckets_folder->GetDirectory(FolderName.str().c_str());
				rFile->cd(bucket_folder->GetPath());
				
				
				tmp.str("");
				tmp << "pedestals_k" << kpix << "_b" << bucket;
				pedestals[kpix][bucket] = new TH1F(tmp.str().c_str(), "Pedestals distribution; Charge [ADC]; #entries", 9000, 0, 9000);
				
				tmp.str("");
				tmp << "pedestals_fc_k" << kpix << "_b" << bucket;
				pedestals_fc[kpix][bucket] = new TH1F(tmp.str().c_str(), "Pedestals distribution; Charge [fC]; #entries", 1000, -100, 100);
				
				tmp.str("");
				tmp << "pedestals_fc_0_127_k" << kpix << "_b" << bucket;
				pedestals_fc_0_127[kpix][bucket] = new TH1F(tmp.str().c_str(), "Pedestals distribution channels 0 to 127; Charge [fC]; #entries", 1000, -100, 100);
				
				tmp.str("");
				tmp << "pedestalsRMS_fc_k" << kpix << "_b" << bucket;
				pedestalsRMS_fc[kpix][bucket] = new TH1F(tmp.str().c_str(), "Pedestals RMS, All Chn; [fC]; a.u.", 1000, 0, 4);
				
				tmp.str("");
				tmp << "pedestalsRMS_fc_0_127_k" << kpix << "_b" << bucket;
				pedestalsRMS_fc_0_127[kpix][bucket] = new TH1F(tmp.str().c_str(), "Pedestals RMS channels 0 to 127, All Chn; [fC]; a.u.", 1000, 0, 4);
				
				tmp.str("");
				tmp << "pedestalsRMS_fc_disc_k" << kpix << "_b" << bucket;
				pedestalsRMS_fc_disc[kpix][bucket] = new TH1F(tmp.str().c_str(), "Pedestals RMS, disc. Chn; [fC]; a.u.", 1000, 0, 4);
				
				tmp.str("");
				tmp << "pedestalsRMS_fc_conn_k" << kpix << "_b" << bucket;
				pedestalsRMS_fc_conn[kpix][bucket] = new TH1F(tmp.str().c_str(), "Pedestals RMS, conn. Chn; [fC]; a.u.", 1000, 0, 4);
				
				
				
				
				tmp.str("");
				tmp << "slope_k" << kpix << "_b" << bucket;
				slope_hist[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope distribution; Slope [ADC/fC]; #entries", 200, -20, 20);
				
				tmp.str("");
				tmp << "slope_conn_k" << kpix << "_b" << bucket;
				slope_hist_conn[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope distribution connected to sensor; Slope [ADC/fC]; #entries", 200, -20, 20);
				
				tmp.str("");
				tmp << "slope_disc_k" << kpix << "_b" << bucket;
				slope_hist_disc[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope distribution disconnected from sensor; Slope [ADC/fC]; #entries", 200, -20, 20);
				
				tmp.str("");
				tmp << "slope_0_127_k" << kpix << "_b" << bucket;
				slope_hist_0_127[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope distribution channels 0 to 127; Slope [ADC/fC]; #entries", 200, -20, 20);
				
				tmp.str("");
				tmp << "slope_residual_k" << kpix << "_b" << bucket;
				slope_residual[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope_residual; Slope_residual[ADC/fC]; #entries", 2000, -100, 100);
				
				tmp.str("");
				tmp << "slopeRMS_k" << kpix << "_b" << bucket;
				slopeRMS[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope RMS; Slope [ADC/fC]; #entries", 1000, 0, 20);
				
				tmp.str("");
				tmp << "slope_vs_channel_k" << kpix << "_b" << bucket;
				slope_vs_channel[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope [ADC/fC]; Channel ID; Slope [ADC/fC]", 1024, -0.5, 1023.5);
				
				tmp.str("");
				tmp << "slope_vs_right_strip_k" << kpix << "_b" << bucket;
				slope_vs_right_strip[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope [ADC/fC]; Strip ID; Slope [ADC/fC]", 920, 919.5, 1839.5);
				
				tmp.str("");
				tmp << "slope_vs_left_strip_k" << kpix << "_b" << bucket;
				slope_vs_left_strip[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope [ADC/fC]; Strip ID; Slope [ADC/fC]", 920, -0.5, 919.5);
				
				tmp.str("");
				tmp << "RMSfC_vs_channel_k" << kpix << "_b" << bucket;
				RMSfC_vs_channel[kpix][bucket] = new TH1F(tmp.str().c_str(), "Pedestal RMS [fC]; Channel ID; Pedestal RMS [fC]", 1024, -0.5, 1023.5);
				
				
			}
		}
	}
				
		
  
  
  rFile->cd(); // move into root folder base
  FolderName.str("");
  FolderName << "Pedestals";
  rFile->mkdir(FolderName.str().c_str()); // produce a sub folder with name of variable FolderName
  TDirectory *pedestal_folder = rFile->GetDirectory(FolderName.str().c_str()); // get path to subdirectory
  pedestal_folder->cd(); // move into subdirectory
  
for (kpix=0; kpix<24; kpix++) 
{
    if ( kpixFound[kpix] ) 
    {
      
		// Get serial number
		tmp.str("");
		tmp << "cntrlFpga(0):kpixAsic(" << dec << kpix << "):SerialNumber";
		serial = dataRead.getConfig(tmp.str());
		
		// Process each channel
		for (channel=minChan; channel <= maxChan; channel++) 
		{
	
			// Show progress
			cout << "\rProcessing baseline kpix " << dec << kpix << " / " << dec << kpixMax
				<< ", Channel " << channel << " / " << dec << maxChan
				<< "                 " << flush;
			
			// Channel is valid
			if ( chanFound[kpix][channel] ) 
			{
	
				// Each bucket
				for (bucket = 0; bucket < 4; bucket++) 
				{                   
					
					// Bucket is valid
					if ( chanData[kpix][channel][bucket][0] != NULL || chanData[kpix][channel][bucket][1] != NULL ) 
					{
					
						// Each range
						for (range = 0; range < 2; range++) 
						{
						
							// Range is valid
							if ( chanData[kpix][channel][bucket][range] != NULL ) 
							{
								chanData[kpix][channel][bucket][range]->computeBase();
						
								// Create histogram
								tmp.str("");
								tmp << "hist_" << serial << "_c" << dec << setw(4) << setfill('0') << channel;
								tmp << "_b" << dec << bucket;
								tmp << "_r" << dec << range;
								tmp << "_k" << dec << kpix;
								
								
								hist = new TH1F(tmp.str().c_str(),
										(tmp.str()+";ADC; Entries / N_of_Cycles").c_str(),
										8192,0,8192);
		  
								//double num_of_entries; //ADDED
								//double normed_bin_content;
								// Fill histogram
								//for (x=0; x < 8192; x++) num_of_entries += chanData[kpix][channel][bucket][range]->baseData[x]; //ADDED
								//cout << num_of_entries << endl;
								for (x=0; x < 8192; x++) 
								{
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
								
								if ( hist->GetFunction("gaus") ) 
								{
									chanData[kpix][channel][bucket][range]->baseFitMean      = hist->GetFunction("gaus")->GetParameter(1);
									chanData[kpix][channel][bucket][range]->baseFitSigma     = hist->GetFunction("gaus")->GetParameter(2);
									chanData[kpix][channel][bucket][range]->baseFitMeanErr   = hist->GetFunction("gaus")->GetParError(1);
									chanData[kpix][channel][bucket][range]->baseFitSigmaErr  = hist->GetFunction("gaus")->GetParError(2);
									
									pedestals[kpix][bucket]->Fill(chanData[kpix][channel][bucket][range]->baseFitMean);
									channel_file_adc_mean << chanData[kpix][channel][bucket][range]->baseFitMean << " " << channel << " " << bucket << endl;
									//pedestalsRMS->Fill(chanData[kpix][channel][bucket][range]->baseFitSigma);
									if ( hist->GetFunction("gaus")->GetNDF() == 0 ) 
									{
										chanData[kpix][channel][bucket][range]->baseFitChisquare = 0;
									} 
									else 
									{
										chanData[kpix][channel][bucket][range]->baseFitChisquare = 
										(hist->GetFunction("gaus")->GetChisquare() / hist->GetFunction("gaus")->GetNDF() );
									}
									
									// Determine bad channel from fitted chisq
									
									if ( findBadMeanChisq && (chanData[kpix][channel][bucket][range]->baseFitChisquare >  meanChisq) ) 
									{
										debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
											<< " Range=" << dec << range 
											<< " Bad fit mean chisq=" << chanData[kpix][channel][bucket][range]->baseFitChisquare << endl;
										badMean[kpix][channel] = true;
										badMeanChisqCnt++;
									}
		    
									// Determine bad channel from fitted mean
									if ( findBadMeanFit &&
										( (chanData[kpix][channel][bucket][range]->baseFitMean > meanMax[range]) ||  
										(chanData[kpix][channel][bucket][range]->baseFitMean < meanMin[range]) ) ) 
									{
										debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
											<< " Range=" << dec << range 
											<< " Bad fit mean value=" << chanData[kpix][channel][bucket][range]->baseFitMean << endl;
										badMean[kpix][channel] = true;
										badMeanFitCnt++;
									}
									// Determine noisy channels from fitted sigma
									if  ((chanData[kpix][channel][bucket][range]->baseFitSigma/chanData[kpix][channel][bucket][range]->baseFitMean) > 0.03) 
									{
										badMean[kpix][channel] = true;
										noiseSigmaCnt++;
										//cout << endl << "Noisy channel with " << "sigma = " << chanData[kpix][channel][bucket][range]->baseFitSigma << endl << endl;
										// channel_file_noise << channel << endl;
									}
		    
									// Determine weird channels with only a single peak if fitted sigma is below the error of sigma
									if (chanData[kpix][channel][bucket][range]->baseFitSigma < chanData[kpix][channel][bucket][range]->baseFitSigmaErr) 
									{
										badMean[kpix][channel] = true;
										errorSigmaCnt++;
										//cout << endl << "Sigma below error value with " << "sigma = " << chanData[kpix][channel][bucket][range]->baseFitSigma << endl << endl;
										// channel_file_bad << channel << endl;
									}
								}
								else if ( findBadMeanFit || findBadMeanChisq ) 
								{
									debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel 
									<< " Bucket=" << dec << bucket << " Range=" << dec << range
									<< " Failed to fit mean" << endl;
									badMean[kpix][channel] = true;
									failedMeanFit++;
								}
		  
								// Determine bad channel from histogram mean
								if ( findBadMeanHist && 
									( (chanData[kpix][channel][bucket][range]->baseMean > meanMax[range]) ||
									(chanData[kpix][channel][bucket][range]->baseMean < meanMin[range]) ) ) 
								{
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
for (kpix=0; kpix<24; kpix++)
{
	if ( kpixFound[kpix] )
	{ 
      
		// Get serial number
		tmp.str("");
		tmp << "cntrlFpga(0):kpixAsic(" << dec << kpix << "):SerialNumber";
		serial = dataRead.getConfig(tmp.str());
		
		// Process each channel
		for (channel=minChan; channel <= maxChan; channel++) 
		{
	
			// Show progress
			cout << "\rProcessing calibration kpix " << dec << kpix << " / " << dec << kpixMax
				<< ", Channel " << channel << " / " << dec << maxChan 
				<< "                 " << flush;
	
			// Channel is valid
			if ( chanFound[kpix][channel] ) 
			{

	  
				// Each bucket
				for (bucket = 0; bucket < 4; bucket++) 
				{
	    
					// Bucket is valid
					if ( chanData[kpix][channel][bucket][0] != NULL || chanData[kpix][channel][bucket][1] != NULL ) 
					{
	      
						// Each range
						for (range = 0; range < 2; range++) 
						{
		
							// Range is valid
							if ( chanData[kpix][channel][bucket][range] != NULL ) 
							{
								chanData[kpix][channel][bucket][range]->computeCalib(chargeError[range]);
								
					
								// Create calibration graph
								grCount = 0;
								for (x=0; x < 256; x++) 
								{
									
									// Calibration point is valid
									if ( chanData[kpix][channel][bucket][range]->calibCount[x] > 0 ) 
									{
										grX[grCount]    = calibCharge ( x, positive[kpix], ((bucket==0)?b0CalibHigh:false));
										grDAC[grCount]  = x;
										grY[grCount]    = chanData[kpix][channel][bucket][range]->calibMean[x];
										grYErr[grCount] = chanData[kpix][channel][bucket][range]->calibError[x];
										grXErr[grCount] = 0;
										debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
											<< " Range=" << dec << range
											<< " Adding point x=" << grX[grCount] 
											<< " Rms=" << chanData[kpix][channel][bucket][range]->calibRms[x]
											<< " Error=" << chanData[kpix][channel][bucket][range]->calibError[x] << endl;
										
										//cout << "Charge in fC : DAC = " << grX[grCount] << " : " << grDAC[grCount] << endl;
										grCount++;
										
										// Find crosstalk, value - base > 3 * sigma
										for (crChan=0; crChan < 1024; crChan++ ) 
										{
											
											if ( chanData[kpix][crChan][bucket][range] != NULL ) 
											{
												
												crossDiff = chanData[kpix][channel][bucket][range]->calibOtherValue[crChan] - 
												chanData[kpix][crChan][bucket][range]->baseMean;
													
												
											}
										}
									}
								}
		  
								// Create graph
								if ( grCount > 0 ) 
								{
	    
									grCalib = new TGraphErrors(grCount,grX,grY,grXErr,grYErr);
									grCalib->Draw("Ap");
									grCalib->GetXaxis()->SetTitle("Charge [C]");
									grCalib->GetYaxis()->SetTitle("ADC");
									grCalib->Fit("pol1","eq","",fitMin[range],fitMax[range]);
									grCalib->GetFunction("pol1")->SetLineWidth(1);
									
									grCalibDAC = new TGraphErrors(grCount,grDAC,grY,grXErr,grYErr);
									grCalibDAC->Draw("Ap");
									grCalibDAC->GetXaxis()->SetTitle("Charge [DAC]");
									grCalibDAC->GetYaxis()->SetTitle("ADC");
									grCalibDAC->Fit("pol1","eq","",fitMin[range],fitMax[range]);
									grCalibDAC->GetFunction("pol1")->SetLineWidth(1);
									
									
									// Create name and write
									tmp.str("");
									tmp << "calib_" << serial << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp << "_b" << dec << bucket;
									tmp << "_r" << dec << range;
									tmp << "_k" << dec << kpix;
									grCalib->SetTitle(tmp.str().c_str());
									grCalib->Write(tmp.str().c_str());
									
									tmp.str("");
									tmp << "calib_DAC" << serial << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp << "_b" << dec << bucket;
									tmp << "_r" << dec << range;
									tmp << "_k" << dec << kpix;
									grCalibDAC->SetTitle(tmp.str().c_str());
									grCalibDAC->Write(tmp.str().c_str());
				
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
									tmp << "_k" << dec << kpix;
									//grResid->SetTitle(tmp.str().c_str());
									grResid->Write(tmp.str().c_str());
				
									// Add to xml
									if ( grCalib->GetFunction("pol1") ) 
									{
										chisqNdf = (grCalib->GetFunction("pol1")->GetChisquare() / grCalib->GetFunction("pol1")->GetNDF());
										Double_t slope = grCalib->GetFunction("pol1")->GetParameter(1);
										Double_t offset = offset;
										
										long double ped_charge = ( chanData[kpix][channel][bucket][range]->baseFitMean ) / slope;
										long double ped_charge_err = (chanData[kpix][channel][bucket][range]->baseHistRMS) / slope ; // simple err
										
										
										pedestals_fc[kpix][bucket]->Fill( ped_charge * pow(10,15) );
										if (channel >= 0 && channel <= 127) pedestals_fc_0_127[kpix][bucket]->Fill( ped_charge * pow(10,15) );
										
										
										pedestalsRMS_fc[kpix][bucket]->Fill( ped_charge_err * pow(10,15) );
										if (channel >= 0 && channel <= 127) pedestalsRMS_fc_0_127[kpix][bucket]->Fill( ped_charge_err * pow(10,15) );

										if ( kpix2strip_left.at(channel) == 9999 ) 
										{
											pedestalsRMS_fc_disc[kpix][bucket]->Fill( ped_charge_err * pow(10,15) );
											slope_hist_disc[kpix][bucket]->Fill( slope / pow(10,15) );
										}
										else
										{
											slope_hist_conn[kpix][bucket]->Fill( slope / pow(10,15) );
											pedestalsRMS_fc_conn[kpix][bucket]->Fill( ped_charge_err * pow(10,15) );
										}
				
										if (ped_charge_err * pow(10,15) >= noise_cut)
										{
											channel_file_noise << channel << endl ;
										}
				
										slope_hist[kpix][bucket]->Fill( slope / pow(10,15) );
										if (channel >= 0 && channel <= 127) slope_hist_0_127[kpix][bucket]->Fill( slope / pow(10,15) );
				
										slope_vs_channel[kpix][bucket]->SetBinContent( channel+1, slope / pow(10,15));
										slope_vs_right_strip[kpix][bucket]->SetBinContent(kpix2strip_right.at(channel)-919, slope / pow(10,15));
										slope_vs_left_strip[kpix][bucket]->SetBinContent(kpix2strip_left.at(channel)+1, slope / pow(10,15));
										
										RMSfC_vs_channel[kpix][bucket]->SetBinContent( channel+1, ped_charge_err * pow(10,15));
										
				
										slope_residual[kpix][bucket]->Fill( offset);
									
										
								
				
										// Determine bad channel from fitted gain
										if ( findBadGainFit && 
											( (grCalib->GetFunction("pol1")->GetParameter(1) > gainMax[range]) ||
											(grCalib->GetFunction("pol1")->GetParameter(1) < gainMin[range]) ) ) 
										{
											debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
												<< " Range=" << dec << range
												<< " Bad gain value=" << grCalib->GetFunction("pol1")->GetParameter(1) << endl;
											badGain[kpix][channel] = true;
											badGainFitCnt++;
										}
				
		      // Determine bad c	hannel from fitted chisq
										if ( findBadGainChisq && (chisqNdf >  gainChisq) ) 
										{
											debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
												<< " Range=" << dec << range
												<< " Bad gain chisq=" << gainChisq << endl;
											badGain[kpix][channel] = true;
											badGainChisqCnt++;
										}
									}
									else 
									{
										if ( findBadGainFit || findBadGainChisq )  
										{
											debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel << " Bucket=" << dec << bucket
												<< " Range=" << dec << range
												<< " Failed to fit gain" << endl;
											badGain[kpix][channel] = true;
											failedGainFit++;
										}
									}
				
									

								}

							}	
						}
					}
				}
	  
				// Determine if the channel is bad
				
				badValue = 0;
				if ( badMean[kpix][channel] ) badValue |= 0x1;
				if ( badGain[kpix][channel] ) badValue |= 0x2;
				
				if ( badValue != 0 ) 
				{
					debug << "Kpix=" << dec << kpix << " Channel=" << dec << channel
					<< " Marking channel bad." << endl;
					badChannelCnt++;
				}
		
			}
		}
	}
}
  cout << endl;
  cout << "Wrote root plots to " << outRoot << endl;
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
  
  rFile->Write();
  delete rFile;
  
  // Cleanup
  for (kpix=0; kpix < 24; kpix++) 
  {
	for (channel=0; channel < 1024; channel++) 
	{
		for (bucket=0; bucket < 4; bucket++) 
		{
			if ( chanData[kpix][channel][bucket][0] != NULL ) delete chanData[kpix][channel][bucket][0];
			if ( chanData[kpix][channel][bucket][1] != NULL ) delete chanData[kpix][channel][bucket][1];
		}
    }
  }
  
  // Close file
  dataRead.close();
  return(0);
}

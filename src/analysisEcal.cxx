	//-----------------------------------------------------------------------------
	// File          : analysis.cpp
	// Author        : Uwe Kraemer (orig. Ryan Herbst) <uwe.kraemer@desy.de>
	// Created       : 06/28/2017
	// Project       : KPiX Analysis
	//-----------------------------------------------------------------------------
	// Description :
	// General analysis of KPiX data.
	//-----------------------------------------------------------------------------
	// Copyright (c) 2012 by SLAC. All rights reserved.
	// Proprietary and confidential to SLAC.
	//-----------------------------------------------------------------------------
	// Modification history :
	// 05/30/2012: created
	// 06/28/2017: large scale rewrite of original calibrationFitter.cpp
	//-----------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TTree.h>
#include <TROOT.h>
//#include "/afs/desy.de/user/k/kraemeru/public/map.cpp"
#include <TPolyLine3D.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TVector.h>
#include <TKey.h>
#include <TClass.h>
#include <stdarg.h>
#include <KpixEvent.h>
#include <KpixSample.h>
#include <Data.h>
#include <DataRead.h>
#include <math.h>
#include <fstream>
#include <XmlVariables.h>
#include <string.h>

#include "kpixmap.h"
using namespace std;
	
	
	
	

//////////////////////////////////////////
// Functions
//////////////////////////////////////////

// Coincidence function void coincidence(int* time_list1, int* time_list2, int* time_list3, int* channel_list1, int* channel_list2, int* channel_list3)

bool gtx_ltz ( int x, int y, int z) // simple function asking if x<y<z
{
	return ((x <= y) && (y <= z));
}

void addDoubleToXml ( ofstream *xml, uint indent, string variable, Double_t value ) { //xml output function
uint x;

if ( !(value != value) ) {
	for (x=0; x < indent; x++) *xml << " ";
	*xml << "<" << variable << ">";
	*xml << value;
	*xml << "</" << variable << ">";
	*xml << endl;
}
}

void addStringToXml ( ofstream *xml, uint indent, string variable, string value ) { //xml output function
uint x;

for (x=0; x < indent; x++) *xml << " ";
*xml << "<" << variable << ">";
*xml << value;
*xml << "</" << variable << ">";
*xml << endl;
}


//////////////////////////////////////////
// Begin of main analysis
//////////////////////////////////////////
int main ( int argc, char **argv )
{
	
TH1::SetDefaultSumw2();

	
//////////////////////////////////////////
// Class declaration and histogram initialization
//////////////////////////////////////////


DataRead               dataRead;  //kpix event classes used for analysis of binary date
off_t                  fileSize;  //
off_t                  filePos;   //
KpixEvent              event;    //
KpixSample             *sample;   //

string                 calState;
uint                   lastPct;
uint                   currPct;
bool                   bucketFound[32][1024][4];  // variable that gives true if bucket has an entry (32 possible number of KPiX, 1024 channels per KPiX, 4 buckets per channel)
bool                   chanFound[32][1024]; // variable that gives true if a channel has entries
bool                   kpixFound[32]; // variable that gives true if a kpix at the index n (0<=n<32) was found
uint                   x;
uint                   range;
uint                   value;
uint                   kpix;
uint                   channel;
uint                   bucket;
double                  tstamp;
string                 serial;
KpixSample::SampleType type;
TH1F                   	*hist[32][1024][4][2];  // #entries/ADC histograms per channel, bucket, kpix and histogram
TH1F					*hist_timed[32][1024][4][2]; //  #entries/time_of_event per channel, bucket, kpix and histogram
//TH1F					*hist_charge[32][1024][4][2]; //  #entries/charge in coulomb if calibration file was given per channel, bucket, kpix and histogram
TH1F					*channel_time[32][1024][4][2];
TH1F					*channel_entries[32][5]; // ADC distribution Total number of events differed per bucket and kpix
TH1F					*strip_entries[32][5];
TH1F					*channel_entries_timed[32][5]; // Time distribution Total number of events differed per bucket and kpix
TH1F					*channel_entries_no_strip[32][5];
TH1F					*kpix_x_ecal[32];
TH1F					*kpix_y_ecal[32];
TH1F					*trigger_difference[32];
TH1F					*channel_entries_no_monster[32][5];
TH1F					*times_kpix[32][5];
TH1F					*times_kpix_monster[32][5];
TH1F					*times_kpix_no_monster[32][5];
TH1F					*trig_count[32][5];
TH1F                   	*hist_buck_sum[32][1024];

TH2F					*kpix_map_ecal[32];

// Stringstream initialization for histogram naming
stringstream           tmp;
stringstream           tmp_units;

stringstream           tmp2;
stringstream           tmp2_units;

stringstream           tmp3;
stringstream           tmp3_units;

// Stringstream initialization for folder naming

stringstream			Folder1;
stringstream			Folder2;
stringstream			Folder3;

ofstream               xml;
ofstream               csv;
uint                   acquisitionCount;
string                 outRoot;
string                 outXml;
string                 outCsv;
TFile                  *rFile;
stringstream           crossString;
stringstream           crossStringCsv;
XmlVariables           config;
ofstream               debug;
ofstream				channel_file_bad;
ofstream				channel_file_bad_fit;
ofstream				channel_file_noise;
ofstream				channel_file_calib;
ofstream 				channel_file_adc_mean;
//double 					calib_slope[1024];
//double					calib_y0[1024];

int						strip[1024];
map_kpix_to_strip(strip);

pixel					pixel_kpix[1024];
pixel_mapping(pixel_kpix);


int 					num_of_channels[32] = {0};
TH1F*					acq_num_ext[32];
TH1F* 					AssignedChannelHist[32][1000];
TH1F* 					AssignedChannelHist_Total[32];
TH1F* 					trigger_difference_per_acq[32][1000];
//TH1F* 					AssignedNumberHist[32][1000];
TH1F*					event_time[32][1000];
TH1F*					event_time_ext[1000];
std::vector<int> monster_channels;



//for (int i = 0; i < 1000; i++)
//{
//	cout << "Strip mapping result, kpix #" << i << " is equal to strip #" << strip[i] << endl;
//}

// Data file is the first and only arg
//if ( argc != 3 && argc != 4 ) {
	//cout << "Usage: calibrationFitter config_file data_file [debug_file]\n";
	//return(1);
//}

//if ( argc == 4 ) debug.open(argv[3],ios::out | ios::trunc);
//bool calibration_check = 0;
//if ( argc == 4 ) calibration_check = 1;

// Read configuration
//if ( ! config.parseFile("config",argv[1]) ) {
	//cout << "Failed to read configuration from " << argv[1] << endl;
	//return(1);
//}
//Read root calibration file
//if ( argc == 4 && ! calibration_file.is_open()) {
	//cout << "Failed to read calibration from " << argv[3] << endl;
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


//////////////////////////////////////////
// Open Data file
//////////////////////////////////////////

if ( ! dataRead.open(argv[1])  ) {
	cout << "Error opening data file " << argv[1] << endl;
	return(1);
}

// Create output names
tmp.str("");
tmp << argv[1] << ".ecal.root";
outRoot = tmp.str();
tmp.str("");
tmp << argv[1] << ".ecal.xml";
outXml = tmp.str();
tmp.str("");
tmp << argv[1] << ".ecal.csv";
outCsv = tmp.str();

//================================================ Charge calibration, currently not working
//Start of change
//if (argc == 4)
//{
	//TFile *calibration_file = TFile::Open(argv[3]);
	//TIter cal_iter(calibration_file->GetListOfKeys());
	//TKey *cal_key;
	//int count = 0;
	//while ((cal_key = (TKey*)cal_iter()))
	//{
		////cout << "Current key1 = " << cal_key->GetClassName() << endl;
		//string title = cal_key->GetTitle();
		//if ((strcmp(cal_key->GetClassName(), "TGraphErrors") == 0) && (title.find("calib") != std::string::npos ))
		//{
			////cout << "test" << endl;
			//TGraphErrors *calib_graph = (TGraphErrors*)cal_key->ReadObj();
			//calib_graph->Fit("pol1","eq","",fitMin[range],fitMax[range]);
			//calib_slope[count] = calib_graph->GetFunction("pol1")->GetParameter(1);
			//calib_y0[count] = calib_graph->GetFunction("pol1")->GetParameter(0);
			////cout << "Title: " << cal_key->GetTitle() << " with Slope = " << calib_slope[count] << endl;
			//++count;
		//}
	//}
//}
//end of change

//////////////////////////////////////////
// Read Data
//////////////////////////////////////////
cout << "Opened data file: " << argv[1] << endl;
fileSize = dataRead.size();
filePos  = dataRead.pos();

// Init
currPct          	= 0;
lastPct          	= 100;
//int eventCount     	= 0;
//minChan          	= 0;
//maxChan          	= 0;
//badTimes         	= 0;
//badMeanFitCnt    	= 0;
//badMeanHistCnt   	= 0;
//badMeanChisqCnt  	= 0;
//badGainFitCnt    	= 0;
//badGainChisqCnt 	= 0;
//failedGainFit   	= 0;
//failedMeanFit   	= 0;
//badChannelCnt   	= 0;
//noiseSigmaCnt		= 0;
//errorSigmaCnt		= 0;
cout << "\rReading File: 0 %" << flush;  // Printout of progress bar
//goodTimes       	= 0;



// Open root file
rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
rFile->cd(); // move into root folder base
Folder1.str("");
Folder1 << "General";
rFile->mkdir(Folder1.str().c_str()); // produce a sub folder with name of variable Folder1
TDirectory *General_folder = rFile->GetDirectory(Folder1.str().c_str()); // get path to subdirectory
General_folder->cd(); // move into subdirectory


//////////////////////////////////////////
// Old histogram generation (needs to be cleaned up)
//////////////////////////////////////////

TH1F *channel_entries_total= new TH1F("Channel_entries_total", "Channel_Entries_total; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
TH1F *channel_entries_total_timed= new TH1F("Channel_entries_total_timed", "Channel_entries_total_timed; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);


TH1F *time_kpix= new TH1F("time_kpix", "time_kpix; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
TH1F *time_external= new TH1F("time external", "time external; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8192.5); // one higher because the accuracy is higher


TH1F *total= new TH1F("Total_response", "total_response; Charge (ADC); #entries/#acq.cycles", 300, -0.5, 8191.5);
TH1F *total_timed= new TH1F("Total_response_timed", "total_response_timed; Charge (ADC); #entries/#acq.cycles", 300, -0.5, 8191.5);


TH1F *beam_ext_time_diff = new TH1F("beam_ext_time_diff", "beam_ext_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 16384, -8192.5, 8191.5);

//TH2F *k28_k30_x_correlation = new TH2F ("x_correlation_k28_k30", "x_correlation_k28_k30; x28/column_width; x30/column_width ", 38,0.5,19,38,0.5,19);
//TH2F *k28_k30_y_correlation = new TH2F ("y_correlation_k28_k30", "y_correlation_k28_k30; x28/row_width; x30/row_width", 40,-18,13,40,-18,13);


//TH2F *k26_k30_x_correlation = new TH2F ("x_correlation_k26_k30", "x_correlation_k26_k30; x26/column_width; x30/column_width ", 38,0.5,19,38,0.5,19);
//TH2F *k26_k30_y_correlation = new TH2F ("y_correlation_k26_k30", "y_correlation_k26_k30; x26/row_width; x30/row_width", 40,-18,13,40,-18,13);


//TH2F *k26_k28_x_correlation = new TH2F ("x_correlation_k26_k28", "x_correlation_k28_k30; x26/column_width; x28/column_width ", 38,0.5,19,38,0.5,19);
//TH2F *k26_k28_y_correlation = new TH2F ("y_correlation_k26_k28", "y_correlation_k26_k28; x26/row_width; x28/row_width", 40,-18,13,40,-18,13);



//TH1F *hist_matched[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);

TH1F *ExtTrigPerCycle = new TH1F ("external_triggers_per_cycle", "ext_trig_per_acq.; #ext_triggers_per_acq.cycle; #entries/#acq.cycles",100,0.5,99.5);

//TH1F *three_coincidence_channel_entries= new TH1F("three_coincidence_channel_entries", "three_coincidence_channel_entries; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
TH1F *full_coincidence_channel_entries= new TH1F("full_coincidence_channel_entries", "full_coincidence_channel_entries; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
Folder1.str("");
Folder1 << "Events";
General_folder->mkdir(Folder1.str().c_str());
TDirectory *gen_event_folder = General_folder->GetDirectory(Folder1.str().c_str());
rFile->cd(gen_event_folder->GetPath());


for (int events = 0; events < 1000; events++) // produce subfolders per event
{
	Folder1.str("");
	Folder1 << "Event_" << events;
	gen_event_folder->mkdir(Folder1.str().c_str());
	TDirectory *events_folder = gen_event_folder->GetDirectory(Folder1.str().c_str());
	events_folder->cd();
	tmp.str("");
	tmp << "time_distribution_external" << "_evt_" << events;
	event_time_ext[events] = new TH1F(tmp.str().c_str(), "time_distribution_external; time/BunchClkCount; #entries/#acq.cycles", 8192, -0.5, 8191.5);
}

while ( dataRead.next(&event) ) // event read to check for filled channels and kpix to reduce number of empty histograms.
{
	acquisitionCount++;
	//if ( acquisitionCount < 10)
	//{
		//cout << "DEBUG: EVENTNUMBER " << event.eventNumber() << endl;
	//}
	for (x=0; x < event.count(); x++)
	{

		//// Get sample
		sample  = event.sample(x);
		kpix    = sample->getKpixAddress();
		channel = sample->getKpixChannel();
		bucket  = sample->getKpixBucket();
		type    = sample->getSampleType();
		if ( type == KpixSample::Data )
		{
			kpixFound[kpix]          = true;
			chanFound[kpix][channel] = true;
			bucketFound[kpix][channel][bucket] = true;
		}
		//cout << "KPIX: " << kpix << endl;
		//cout << "Channel: " << channel << endl;
		kpixFound[0] = false; // for some reason the system finds a kpix in slot 0
	}
}
dataRead.close();
range = 0;
//for (int u = 0; u<32; ++u)
//{
	//if (kpixFound[u])
	//{
		//for (int j =0; j<1024; ++j)
		//{
			//if (chanFound[u][j] == 0)
			//{
				//cout << "DEBUG  KPIX: " << u << " Channel: " << j << endl;
			//}
		//}
	//}
//}


double weight = 1.0/acquisitionCount; //normalization weight  #entries*weight = #entries/acq.cycle
int monster_counter[32] = {0}; // kpix monster counter


//////////////////////////////////////////
// New histogram generation within subfolder structure
//////////////////////////////////////////

for (kpix = 0; kpix < 32; kpix++) //looping through all possible kpix
{
	//
	//cout << "DEBUG test " << kpixFound[kpix] << endl;
	if (kpixFound[kpix]) //checking if kpix exists
	{
		rFile->cd(); //producing subfolder for kpix same as above for the event subfolder structure
		Folder1.str("");
		Folder1 << "KPiX_" << kpix;
		rFile->mkdir(Folder1.str().c_str());
		TDirectory *kpix_folder = rFile->GetDirectory(Folder1.str().c_str());
		kpix_folder->cd();
		tmp.str("");
		tmp << "Channel_entries_k_" << kpix << "_total";
		channel_entries[kpix][4] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
		tmp.str("");
		tmp << "Channel_entries_k_" << kpix << "_total_no_monster";
		channel_entries_no_monster[kpix][4] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
		tmp.str("");
		tmp << "Channel_entries_k_" << kpix << "_total_timed";
		channel_entries_timed[kpix][4] = new TH1F(tmp.str().c_str(), "Channel_Entries_timed; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
		tmp.str("");
		tmp << "Channel_entries_k_" << kpix << "_total_no_strip";
		channel_entries_no_strip[kpix][4] = new TH1F(tmp.str().c_str(), "Channel_Entries_no_strip; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
		tmp.str("");
		tmp << "Strip_entries_k_" << kpix << "_total";
		strip_entries[kpix][4] = new TH1F(tmp.str().c_str(), "Strip_Entries; Strip_address; #entries/#acq.cycles", 928,-0.5, 927.5);


		tmp.str("");
		tmp << "timestamp_kpix_k_" << kpix  << "_total";
		times_kpix[kpix][4] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time/#bunch_clk_count; #entries/#acq.cycles", 300,-0.5, 8191.5);
		tmp.str("");
		tmp << "timestamp_kpix_k_" << kpix  << "_monster_total";
		times_kpix_monster[kpix][4] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time/#bunch_clk_count; #entries/#acq.cycles", 300,-0.5, 8191.5);
		tmp.str("");
		tmp << "timestamp_kpix_k_" << kpix  << "_no_monster_total";
		times_kpix_no_monster[kpix][4] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time/#bunch_clk_count; #entries/#acq.cycles", 300,-0.5, 8191.5);
		tmp.str("");
		tmp << "trig_count_k" << kpix << "_total";
		trig_count[kpix][4]  = new TH1F (tmp.str().c_str(), "trig_count;  #triggered channels; #entries/#acq.cycles",1024, -0.5, 1023.5);
		tmp.str("");
		tmp << "beam_ext_time_diff_" << kpix;
		trigger_difference[kpix] = new TH1F (tmp.str().c_str(), "intern_extern_trig_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 16384, -8192.5, 8191.5);
		tmp.str("");
		tmp << "assigned_channels_k" << kpix << "_total";
		AssignedChannelHist_Total[kpix] = new TH1F (tmp.str().c_str(), "assigned_channels_per_ext_trig;   #assigned_channels; #entries/#acq.cycles",40, -0.5, 39.5);

		tmp.str("");
		tmp << "acq_num_ext_k_" << kpix;
		acq_num_ext[kpix] = new TH1F(tmp.str().c_str(), "acq_num_ext; #triggers/acq._cycle; #entries/#acq.cycles",5, -0.5, 4.5);
		
		tmp.str("");
		tmp << "k" << kpix << "_x";
		kpix_x_ecal[kpix] = new TH1F (tmp.str().c_str(), "kpix_x_ecal_hit_distribution; x/column_width; #entries/#acq.cycles",38,0.5,19);

		tmp.str("");
		tmp << "k" << kpix << "_y";
		kpix_y_ecal[kpix] = new TH1F (tmp.str().c_str(), "kpix_y_ecal_hit_distribution; y/row_width; #entries/#acq.cycles",40,-18,13);
		
		tmp.str("");
		tmp << "k" << kpix << "_map";
		kpix_map_ecal[kpix] = new TH2F (tmp.str().c_str(), "kpix_ecal_map; x/column_width; y/row_width",19,0.0,19,31,-18,13);

		for (bucket = 0; bucket< 4; bucket++)
		{
			tmp.str("");
			tmp << "Channel_entries_k_" << kpix << "_b" << bucket;
			channel_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
			tmp.str("");
			tmp << "Channel_entries_k_" << kpix << "_b" << bucket << "_no_strip";
			channel_entries_no_strip[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries_no_strip; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
			tmp.str("");
			tmp << "Strip_entries_k_" << kpix << "_b" << bucket;
			strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Strip_Entries; Strip_address; #entries/#acq.cycles", 928,-0.5, 927.5);

			tmp.str("");
			tmp << "Channel_entries_k_" << kpix <<  "_b" << bucket << "_timed";
			channel_entries_timed[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries_timed; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
			tmp.str("");
			tmp << "Channel_entries_k_" << kpix << "_b" << bucket << "_no_monster";
			channel_entries_no_monster[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
			tmp.str("");
			tmp << "timestamp_kpix_k_" << kpix << "_b" << bucket;
			times_kpix[kpix][bucket] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time/#bunch_clk_count; #entries/#acq.cycles", 300,-0.5, 8191.5);
			tmp.str("");
			tmp << "timestamp_kpix_monster_k_" << kpix << "_b" << bucket;
			times_kpix_monster[kpix][bucket] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time/#bunch_clk_count; #entries/#acq.cycles", 300,-0.5, 8191.5);
			tmp.str("");
			tmp << "trig_count_k" << kpix << "_b" << bucket ;
			trig_count[kpix][bucket]  = new TH1F(tmp.str().c_str(), "trig_count;  #triggered channels; #entries/#acq.cycles",1024, -0.5,1023.5);
			tmp.str("");
			tmp << "timestamp_kpix_k_" << kpix << "_b" << bucket << "_no_monster";
			times_kpix_no_monster[kpix][bucket] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time/#bunch_clk_count; #entries/#acq.cycles", 300,-0.5, 8191.5);
		}
		Folder1.str("");
		Folder1 << "Events";
		kpix_folder->mkdir(Folder1.str().c_str());
		TDirectory *event_folder = kpix_folder->GetDirectory(Folder1.str().c_str());
		rFile->cd(event_folder->GetPath());
		for (int events = 0; events < 1000; events++)
		{
			Folder1.str("");
			Folder1 << "Event_" << events;
			event_folder->mkdir(Folder1.str().c_str());
			TDirectory *events_folder = event_folder->GetDirectory(Folder1.str().c_str());
			events_folder->cd();
			tmp.str("");
			tmp << "time_distribution_k_" << kpix << "_evt_" << events;
			event_time[kpix][events] = new TH1F(tmp.str().c_str(), "time_distribution; time/BunchClkCount; #entries", 8192, -0.5, 8191.5);
			tmp.str("");
			tmp << "assigned_channels_k" << kpix << "_evt_" << events;
			AssignedChannelHist[kpix][events]  = new TH1F (tmp.str().c_str(), "assigned_channels_per_ext_trig;  external_trigger_number; #assigned_channels ",100, -0.5, 99.5);
			tmp.str("");
			tmp << "trigger_difference_k" << kpix << "_evt_" << events;
			trigger_difference_per_acq[kpix][events]  = new TH1F (tmp.str().c_str(), "trigger_difference;  #entries/#acq.cycles; #Delta T (BunchClkCount) ",16384, -8192.5, 8191.5);
			//tmp.str("");
			//tmp << "assigned_number_k" << kpix << "_evt_" << events;
			//AssignedNumberHist[kpix][events]  = new TH1F (tmp.str().c_str(), "assigned_NumberOfChannel_per_ext_trig;  #same_assignement; #entries/#acq.cycles",40,0,40);
		}
		Folder1.str("");
		Folder1 << "Channels";
		kpix_folder->mkdir(Folder1.str().c_str());
		TDirectory *channels_folder = kpix_folder->GetDirectory(Folder1.str().c_str());
		rFile->cd(channels_folder->GetPath());
		for (channel = 0; channel < 1024; channel++)
		{
			if (chanFound[kpix][channel])
			{
				Folder1.str("");
				Folder1 << "Channel_" << channel;
				channels_folder->mkdir(Folder1.str().c_str());
				TDirectory *channel_folder = channels_folder->GetDirectory(Folder1.str().c_str());
				rFile->cd(channel_folder->GetPath());

				tmp.str("");
				tmp << "hist" << "_c" << dec << setw(4) << setfill('0') << channel;
				tmp << "_k" << dec << kpix;

				tmp_units.str("");
				tmp_units << "hist" << "_c" << dec << setw(4) << setfill('0') << channel;
				tmp_units << "_k" << dec << kpix;
				tmp_units << "; Charge (ADC); #entries/#acq.cycles";

				hist_buck_sum[kpix][channel] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5, 8191.5);

				num_of_channels[kpix] = num_of_channels[kpix] + 1;

				for (bucket = 0; bucket < 4; bucket++)
				{
					if (bucketFound[kpix][channel][bucket])
					{
						// Naming of histograms and generating of histograms




						tmp.str("");  //set stringstream tmp to an empty string
						tmp << "hist" << "_c" << dec << setw(4) << setfill('0') << channel; // add hist_c0$channel
						tmp << "_b" << dec << bucket; // add _b$bucket
						tmp << "_k" << dec << kpix; // add _k$kpix to stringstream

						tmp_units.str(""); //set stringstream decribing histogram units to an empty string
						tmp_units << "hist" << "_c" << dec << setw(4) << setfill('0') << channel;  // add hist_c0$channel
						tmp_units << "_b" << dec << bucket; // add _b$bucket
						tmp_units << "_k" << dec << kpix; // add _k$kpix to stringstream
						tmp_units << "; Charge (ADC); #entries/#acq.cycles"; // add title: x label, y label to stringstream
					
						hist[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);

						tmp.str("");
						tmp << "hist_timed" << "_c" << dec << setw(4) << setfill('0') << channel;
						tmp << "_b" << dec << bucket;
						tmp << "_k" << dec << kpix;
						tmp << "_time_cut";

						tmp_units.str("");
						tmp_units << "hist_timed" << "_c" << dec << setw(4) << setfill('0') << channel;
						tmp_units << "_b" << dec << bucket;
						tmp_units << "_k" << dec << kpix;
						tmp_units << "; Charge (ADC); #entries/#acq.cycles";
					
						hist_timed[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);

						tmp.str("");
						tmp << "time" << "_c" << dec << setw(4) << setfill('0') << channel;
						tmp << "_b" << dec << bucket;
						tmp << "_r" << dec << range;
						tmp << "_k" << dec << kpix;

						tmp_units.str("");
						tmp_units << "time" << "_c" << dec << setw(4) << setfill('0') << channel;
						tmp_units << "_b" << dec << bucket;
						tmp_units << "_k" << dec << kpix;
						tmp_units << "; Time/bunchClkCount; #entries/#acq.cycles";

						channel_time[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),300, -0.5,8191.5);


						tmp.str("");
						tmp << "charge" << "_c" << dec << setw(4) << setfill('0') << channel;
						tmp << "_b" << dec << bucket;
						tmp << "_r" << dec << range;
						tmp << "_k" << dec << kpix;

						tmp_units.str("");
						tmp_units << "charge" << "_c" << dec << setw(4) << setfill('0') << channel;
						tmp_units << "_b" << dec << bucket;
						tmp_units << "_k" << dec << kpix;
						tmp_units << "; Charge/fC; #entries/#acq.cycles";

						//hist_charge[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),300, -0.5,8191.5);
						
						
					}
				}
			}
		}
	}
}

//////////////////////////////////////////
// Data read for first 1000 events for detailed look into single event structure
//////////////////////////////////////////
dataRead.open(argv[1]); //open binary file
int event_num = 0;
int event_num_ext = -1;
while ( dataRead.next(&event) ) //loop through binary file event structure until end of file
{

	double 	trigger_counter[32] = {0}; // fill the entire list of trigger_counter with 0
	
	for (x=0; x < event.count(); x++)  //within the binary file go through each event
	{
		//// Get sample
		sample  = event.sample(x);  // check event subtructure
		kpix    = sample->getKpixAddress();
		channel = sample->getKpixChannel();
		bucket  = sample->getKpixBucket();
		value   = sample->getSampleValue();
		type    = sample->getSampleType();
		tstamp  = sample->getSampleTime();
		range   = sample->getSampleRange();
		if (type == 2) //if type of event is ==2, the event is of type external timestamp
		{
			//cout << event_num << endl;
			//cout << event_num_ext << endl;
			if (x == 0) event_num_ext++;
			double time = tstamp + double(value * 0.125);
			if (event_num_ext < 1000)
			{
				//cout << event_num_ext <<  " " << time << endl;
				event_time_ext[event_num_ext]->Fill(time);
			}
		}
		if ( type == KpixSample::Data ) // If event is of type KPiX data
		{
			channel_entries[kpix][bucket]->Fill(channel, weight);
			channel_entries[kpix][4]->Fill(channel, weight);

			strip_entries[kpix][bucket]->Fill(strip[channel], weight);
			strip_entries[kpix][4]->Fill(strip[channel], weight);

			times_kpix[kpix][bucket]->Fill(tstamp, weight);
			times_kpix[kpix][4]->Fill(tstamp, weight);
			trigger_counter[kpix] = trigger_counter[kpix] + (1.0/num_of_channels[kpix]);
			if (strip[channel] == 9999)
			{
				channel_entries_no_strip[kpix][bucket]->Fill(channel,weight);
				channel_entries_no_strip[kpix][4]->Fill(channel,weight);
			}
			if (event.count() > 700)
			{
				times_kpix_monster[kpix][bucket]->Fill(tstamp, weight);
				times_kpix_monster[kpix][4]->Fill(tstamp, weight);
			}
			else
			{
				channel_entries_no_monster[kpix][bucket]->Fill(channel, weight);
				channel_entries_no_monster[kpix][4]->Fill(channel, weight);
				times_kpix_no_monster[kpix][bucket]->Fill(tstamp, weight);
				times_kpix_no_monster[kpix][4]->Fill(tstamp, weight);
			}
			if (event_num < 1000)
			{
				 event_time[kpix][event_num]->Fill(tstamp);
			 }
		}
	}
	event_num++;
	//if (trigger_counter[26] > 4) cout << trigger_counter[26] << endl;
	if (kpixFound[26]) acq_num_ext[26]->Fill(trigger_counter[26]); // trigger counting for monster check
	if (kpixFound[28]) acq_num_ext[28]->Fill(trigger_counter[28]);
	if (kpixFound[30]) acq_num_ext[30]->Fill(trigger_counter[30]);
}
dataRead.close(); // close file as we have looped through it and are now at the end
dataRead.open(argv[1]); //open file again to start from the beginning

int two_coincidence = 0;
int three_coincidence = 0;
int extern_trigger_id={0};


event_num = 0;


while ( dataRead.next(&event) )
{
	std::vector<double> time_ext;
	std::vector<int> channel_hits[32];
	std::vector<int> timestamp[32];
	std::vector<int> adc_value[32];
	std::vector<double> time_diff_kpix_ext[32];
	std::vector<int> AssignedTrigger[32];
	
	//std::vector<int> Assignment_number;
	int num_trig_count[32][5] = {0};
	//cout << " NEW EVENT " << endl;
	for (x=0; x < event.count(); x++)
	{
		//cout << "DEBUG: EVENT COUNT " << event.count() << endl;
		//// Get sample
		sample  = event.sample(x);
		kpix    = sample->getKpixAddress();
		channel = sample->getKpixChannel();
		bucket  = sample->getKpixBucket();
		value   = sample->getSampleValue();
		type    = sample->getSampleType();
		tstamp  = sample->getSampleTime();
		range   = sample->getSampleRange();





		if (type == 2)// If event is of type external timestamp
		{
			double time = tstamp + double(value * 0.125);
			time_external->Fill(time, weight);
			time_ext.push_back(time);
			//cout << "DEBUG: channel in timestmap = " << channel << endl;
			//cout << "DEBUG: bucket in timestmap = " << bucket << endl;
		}


		if ( type == KpixSample::Data ) // If event is of type KPiX data
		{

				channel_hits[kpix].push_back(channel);
				timestamp[kpix].push_back(tstamp);
				adc_value[kpix].push_back(value);
				time_kpix->Fill(tstamp, weight);
				
				std::vector<double> trig_diff_list;
				trig_diff_list.push_back(1);
				trig_diff_list.push_back(5);
				trig_diff_list.push_back(10);
				trig_diff_list.push_back(0.5);

				hist[kpix][channel][bucket][0]->Fill(value, weight);
				//hist_charge[kpix][channel][bucket][0]->Fill(double(value)/calib_slope[channel]*pow(10,15) , weight);
				hist_buck_sum[kpix][channel]->Fill(value,weight);
				channel_entries_total->Fill(channel, weight);
				channel_time[kpix][channel][bucket][0]->Fill(tstamp, weight);
				total->Fill(value, weight);


				kpix_map_ecal[kpix]->Fill(pixel_kpix[channel].x, pixel_kpix[channel].y);
				kpix_x_ecal[kpix]->Fill(pixel_kpix[channel].x, weight);
				kpix_y_ecal[kpix]->Fill(pixel_kpix[channel].y, weight);
				
				num_trig_count[kpix][bucket] += 1;
				num_trig_count[kpix][4] += 1;
				

				//cout << "DEBUG channel number: " << channel << endl;
				//cout << "DEBUG X coord: " << pixel_kpix[499].x << endl;
				//cout << "DEBUG y coord: " << pixel_kpix[499].y << endl;



				// Check for time difference between external time stamp and internal time stamp for noise filtering
				double trig_diff = 8200.0;
				int assigned_number;
				if (time_ext.size() > 0) //only calculate the time difference between triggers if there are some external triggers
				{
					for (unsigned int j = 0; j < time_ext.size(); ++j)
					{
						trig_diff_list.push_back(tstamp-time_ext.at(j));
						if (fabs(trig_diff) > fabs(tstamp-time_ext.at(j)))
						{
								trig_diff = tstamp-time_ext.at(j);
	
								assigned_number = j;
						}
						else
						{
							cout << "Difference not lower than before" << endl;
							cout << "Channel time stamp = " << tstamp << endl;
							cout << "External match = " << time_ext.at(j) << endl;
							cout << "Old difference = " << trig_diff << endl;
							cout << "New difference = " << tstamp-time_ext.at(j) << endl;
						}
						
					}
					//cout << "Trig diff old method = " <<  trig_diff << endl;
					//if (trig_diff_list.size() > 0) cout << "Trig diff new method = " <<  *std::min_element(trig_diff_list.begin(), trig_diff_list.end()) << endl; //seg fault when vector is empty
					//if (trig_diff_list.size() > 0) cout << "Trig diff new method position in vector = " <<  distance(trig_diff_list.begin(), min_element(trig_diff_list.begin(), trig_diff_list.end())) << endl; //seg fault when vector is empty
					assigned_number =  distance(trig_diff_list.begin(), min_element(trig_diff_list.begin(), trig_diff_list.end())); //position of smallest element in trigger difference vector
					time_diff_kpix_ext[kpix].push_back(trig_diff);
					//cout << assigned_number << endl;
					if (event_num < 1000)
					  {
					    AssignedChannelHist[kpix][event_num]->Fill(assigned_number);
					    trigger_difference_per_acq[kpix][event_num]->Fill(trig_diff);
					  }
					
					AssignedTrigger[kpix].push_back(assigned_number);
					//Assignment_number.push_back(assigned_number);
					if((trig_diff >= 0.0 )  && (trig_diff  <= 3.0) )
					  {
					    hist_timed[kpix][channel][bucket][0]->Fill(value, weight);
					    total_timed->Fill(value, weight);
					    channel_entries_total_timed->Fill(channel, weight);
					    channel_entries_timed[kpix][bucket]->Fill(channel, weight);
					    channel_entries_timed[kpix][4]->Fill(channel, weight);
					  }
					
					beam_ext_time_diff->Fill(trig_diff, weight);
					trigger_difference[kpix]->Fill(trig_diff, weight);
				}
				//if (kpix != 26 && kpix != 28 && kpix != 30) cout << "Weird..." << kpix << endl;
			//}
				
		}
	}
	for (kpix = 0; kpix < 32; kpix++)
	  {
	    if (kpixFound[kpix])
	      {
		for (uint trig = 0; trig < time_ext.size(); trig++)
		  {
		    int mycount = std::count(AssignedTrigger[kpix].begin(), AssignedTrigger[kpix].end(), trig);
		    AssignedChannelHist_Total[kpix]->Fill(mycount);
		  }
	      }
	  }
	//for (unsigned int k = 0; k < Assignment_number.size(); k++)
	//{
	//cout << "DEBUG: Assignment_Number " << k << " = " << Assignment_number.at(k) << endl;
		//if (event_num < 1000) AssignedNumberHist[kpix][event_num]->Fill(Assignment_number.at(k));
	//}
	ExtTrigPerCycle->Fill(time_ext.size());
	for (int h = 0; h<32; ++h)
	  {
	    if (kpixFound[h])
		{
			for (int q = 0; q<4; ++q)
			{
				if (num_trig_count[h][q] > 800)
				{
					monster_counter[h]++;
				}
				trig_count[h][q]->Fill(num_trig_count[h][q], weight);
				trig_count[h][4]->Fill(num_trig_count[h][q], weight);
			}
			
		}
	}	
	

	//////////////////////////////////////////
	// x_y correlation plots, take a lot of time and get screwed up when analysing external triggering data files
	//////////////////////////////////////////
	//for (int q = 0; q < timestamp[28].size(); q++)
	//{
		//for (int l = 0; l < timestamp[30].size(); l++)
		//{
			//if (timestamp[28].at(q) == timestamp[30].at(l))
			//{
				//k28_k30_x_correlation->Fill(pixel_kpix[channel_hits[28].at(q)].x, pixel_kpix[channel_hits[30].at(l)].x);
				//k28_k30_y_correlation->Fill(pixel_kpix[channel_hits[28].at(q)].y, pixel_kpix[channel_hits[30].at(l)].y);
			//}
		//}
	//}
	//for (int q = 0; q < timestamp[26].size(); q++)
	//{
		//for (int l = 0; l < timestamp[30].size(); l++)
		//{
			//if (timestamp[26].at(q) == timestamp[30].at(l))
			//{
				//k26_k30_x_correlation->Fill(pixel_kpix[channel_hits[26].at(q)].x, pixel_kpix[channel_hits[30].at(l)].x);
				//k26_k30_y_correlation->Fill(pixel_kpix[channel_hits[26].at(q)].y, pixel_kpix[channel_hits[30].at(l)].y);
			//}
		//}
	//}
	//for (int q = 0; q < timestamp[26].size(); q++)
	//{
		//for (int l = 0; l < timestamp[28].size(); l++)
		//{
			//if (timestamp[26].at(q) == timestamp[28].at(l))
			//{
				//k26_k28_x_correlation->Fill(pixel_kpix[channel_hits[26].at(q)].x, pixel_kpix[channel_hits[28].at(l)].x);
				//k26_k28_y_correlation->Fill(pixel_kpix[channel_hits[26].at(q)].y, pixel_kpix[channel_hits[28].at(l)].y);
				////for (int k = 0; k < timestamp[30].size(); k++)
				////{
					////if (timestamp[28].at(l) == timestamp[30].at(k))
					////{
						////track_plot->SetPoint(0, pixel_kpix[channel_hits[26].at(q)].x, pixel_kpix[channel_hits[26].at(q)].y, 0);
						////track_plot->SetPoint(1, pixel_kpix[channel_hits[28].at(l)].x, pixel_kpix[channel_hits[28].at(l)].y, 1);
						////track_plot->SetPoint(2, pixel_kpix[channel_hits[30].at(k)].x, pixel_kpix[channel_hits[30].at(k)].y, 2);
					////}
				////}
			//}
		//}
	//}


	//////////////////////////////////////////
	// Triggering efficiency and coincidence calculation, takes a lot of time.
	//////////////////////////////////////////
	std::vector<int> kpix_matched_time;
	std::vector<int> kpix_matched_channel;
	//int map_range = 1.0;
	//int time_range = 0;
	//for (int j = 0; j < adc_value[26].size(); ++j)
	//{
		//for (int i = 0; i < adc_value[28].size(); ++i)
		//{
			//for (int q = 0; q < adc_value[30].size(); ++q)
			//{
				//if ((timestamp[26].at(j) == timestamp[28].at(i)) && (timestamp[26].at(j) == timestamp[30].at(q)) && (channel_list[26].at(j) == channel_list[28].at(i)) && (channel_list[26].at(j) == channel_list[30].at(q))) three_coincidence++;
				//if ((timestamp[26].at(j) == timestamp[28].at(i)) && (timestamp[26].at(j) != timestamp[30].at(q)) && (channel_list[26].at(j) == channel_list[28].at(i)) && (channel_list[26].at(j) != channel_list[30].at(q))) two_coincidence++;
				//if ((timestamp[26].at(j) == timestamp[30].at(q)) && (timestamp[26].at(j) != timestamp[28].at(i)) && (channel_list[26].at(j) == channel_list[30].at(q)) && (channel_list[26].at(j) != channel_list[28].at(i))) two_coincidence++;
				//if ((timestamp[28].at(i) == timestamp[30].at(q)) && (timestamp[26].at(j) != timestamp[28].at(i)) && (channel_list[28].at(i) == channel_list[30].at(q)) && (channel_list[26].at(j) != channel_list[28].at(i))) two_coincidence++;
			//}
		//}
	//}
	//double x_0;
	//double y_0;
	//double dist_r;
	//for (unsigned int j = 0; j < adc_value[26].size(); ++j)
	//{
		//x_0 = pixel_kpix[channel_hits[26].at(j)].x;
		//y_0 = pixel_kpix[channel_hits[26].at(j)].y;
		//for (unsigned int i = 0; i < adc_value[30].size(); ++i)
		//{
			//dist_r = sqrt( pow(pixel_kpix[channel_hits[30].at(i)].x-x_0,2) + pow(pixel_kpix[channel_hits[30].at(i)].y-y_0,2) );
			//if (gtx_ltz(timestamp[30].at(i)-time_range, timestamp[26].at(j), timestamp[30].at(i)+time_range))
			//{
				//if (dist_r <= map_range)
				//{
					////if (gtx_ltz(0, time_diff_kpix_ext[26].at(j), 3))
					////{
						//kpix_matched_time.push_back(timestamp[26].at(j));
						//kpix_matched_channel.push_back(channel_hits[26].at(j));
					////}
				//}

			//}
		//}
	//}
	//for (unsigned int j = 0; j < kpix_matched_time.size(); j++)
	//{
		//x_0 = pixel_kpix[kpix_matched_channel.at(j)].x;
		//y_0 = pixel_kpix[kpix_matched_channel.at(j)].y;
		//for (unsigned int i = 0; i < adc_value[28].size(); ++i)
		//{
			//dist_r = sqrt( pow(pixel_kpix[channel_hits[28].at(i)].x-x_0,2) + pow(pixel_kpix[channel_hits[28].at(i)].y-y_0,2) );
			//if ((gtx_ltz(timestamp[28].at(i)-time_range, kpix_matched_time.at(j), timestamp[28].at(i)+time_range)))
			//{
				//if (dist_r <= map_range)
				//{
					//if (gtx_ltz(0, time_diff_kpix_ext[28].at(i), 3))
					//{
						//full_coincidence_channel_entries->Fill(channel_hits[28].at(i), weight);
					//}
					//three_coincidence_channel_entries->Fill(channel_hits[28].at(i), weight);


					//three_coincidence=three_coincidence+1;
					////cout <<
				//}
			//}
			//else two_coincidence=two_coincidence+1;
		//}
	//}
	
	
	extern_trigger_id=extern_trigger_id+time_ext.size();  // Counting which global external trigger was matched to a channel
	
  ////   Show progress
	event_num++;
	filePos  = dataRead.pos();
	currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
	if ( currPct != lastPct )
	{
		cout << "\rReading File: " << currPct << " %      " << flush;
		lastPct = currPct;
	}
}
//for (kpix = 0; kpix < 32; kpix++)
//{
	//if (kpixFound[kpix])
	//{

		//for (uint ext_trig = 0; ext_trig < AssignedChannel_Total[kpix].size(); ext_trig++)
		//{
			//int trigger_counter = 0;
			//for (int i =0; i < extern_trigger_id; ++i)
			//{
				//if (i == AssignedChannel_Total[kpix].at(ext_trig)) trigger_counter++;
			//}
			//AssignedChannelHist_Total[kpix]->Fill(trigger_counter);
		//}
	//}
//}

////// Writing of histograms
//double first = (hist_calib->FindFirstBinAbove(0))* (2000.0/8192) -1;
//double last = hist_calib->FindLastBinAbove(0)* (2000.0/8192) +1;

//hist_calib->GetXaxis()->SetRangeUser(first, last);



////TH2.SetOption("COLZ");
//k28_k30_x_correlation->SetOption("COLZ");
//k28_k30_y_correlation->SetOption("COLZ");
//k30_map->SetOption("COLZ");
//k26_k30_x_correlation->SetOption("COLZ");
//k26_k30_y_correlation->SetOption("COLZ");
//k26_map->SetOption("COLZ");
//k26_k28_x_correlation->SetOption("COLZ");
//k26_k28_y_correlation->SetOption("COLZ");
//k28_map->SetOption("COLZ");

//double r = sqrt(pow(pixel_kpix[771].x-pixel_kpix[717].x,2) + pow(pixel_kpix[771].y-pixel_kpix[717].y,2));
//cout << "pixel distance: " << r << endl;
//r = sqrt(pow(pixel_kpix[771].x-pixel_kpix[772].x,2) + pow(pixel_kpix[771].y-pixel_kpix[772].y,2));
//cout << "pixel distance: " << r << endl;




cout << "Full coincidence of sensors with external trigger: " << full_coincidence_channel_entries->GetEntries() << endl;
cout << "Three coincidence of sensors: " << three_coincidence << endl;
cout << "Two coincidence of sensors: " << two_coincidence << endl;
for (kpix = 0; kpix < 32; kpix++)
{
	if (kpixFound[kpix])
	{
		cout << endl << "Number of monster events in " << kpix << " = " << monster_counter[kpix] << endl;
		cout << "Number of normed monster events in " << kpix << " = " << monster_counter[kpix]*weight << endl;
		cout << "Number of entries in KPiX" << kpix << " = " << channel_entries[kpix][4]->GetEntries() << endl;
	}
}

//for (int k = 0; k < 1024; k++)
//{
	//cout << "DEBUG channel number: " << k << endl;
	//cout << "DEBUG X coord: " << pixel_kpix[k].x << endl;
	//cout << "DEBUG y coord: " << pixel_kpix[k].y << endl << endl;
//}

cout << endl;
cout << "Writing root plots to " << outRoot << endl;
cout << "Writing xml data to " << outXml << endl;
cout << "Writing csv data to " << outCsv << endl;
cout << endl;
rFile->Write();
gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever

rFile->Close();


//cout << "DEBUG" << pixel_kpix[872].x << endl;
//cout << "DEBUG" << pixel_kpix[872].y << endl;
//cout << "DEBUG" << pixel_kpix[810].x << endl;
//cout << "DEBUG" << pixel_kpix[810].y << endl;


dataRead.close();
return(0);
}

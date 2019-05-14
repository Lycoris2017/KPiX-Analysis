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
// 22/03/2018: clean up ecal plots and add strip plots by <mengqing.wu@desy.de>
//-----------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TTree.h>
#include <TROOT.h>
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
#include <algorithm>
#include <stdlib.h>

#include "kpixmap.h"
#include "kpix_left_and_right.h"
#include "TBFunctions.h"
using namespace std;


//////////////////////////////////////////
// Global Variables
//////////////////////////////////////////
vector<TGraphErrors*> calib_graphs; //needed for current loopdir	
//////////////////////////////////////////
// Functions
//////////////////////////////////////////


void loopdir(TDirectory* dir, string histname)
{
	
	TDirectory *dirsav = gDirectory;
	TIter keys_iter(dir->GetListOfKeys());
	TKey* key;
	
	while ((key = (TKey*)keys_iter()))
	{
		if (key->IsFolder())
		{
			dir->cd(key->GetName());
			TDirectory *subdir = gDirectory;
			//subfolder->cd();
			loopdir(subdir, histname);
			dirsav->cd();
			continue;
		}
		else
		{
			string keyname = key->GetName();
			size_t found = keyname.find(histname);
			if (found == 0) 
			{
				TGraphErrors *calib_graph = (TGraphErrors*)key->ReadObj();
				calib_graph->SetName(key->GetName());
				calib_graphs.push_back(calib_graph);
			}
		}
	}
}



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

double smallest_time_diff( vector<double> ext_list, int int_value)
{
	double trigger_diff = 8200.0;
	for (uint k = 0; k<ext_list.size(); ++k)
	{
		double delta_t = int_value-ext_list[k];
		if (fabs(trigger_diff) > fabs(delta_t) && delta_t > 0) 
		{
			trigger_diff = delta_t;
		}
	}
	return trigger_diff;
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
	
	// cycles to skip in front:
	long int                   skip_cycles_front;
	FILE*                  f_skipped_cycles;
	string                 outtxt;
	
	int kpix_checking = 24;
	int bucket_checking = 1;
	
	string                 calState;
	uint                   lastPct;
	uint                   currPct;
	bool                   bucketFound[kpix_checking][1024][bucket_checking];  // variable that gives true if bucket has an entry (24 possible number of KPiX, 1024 channels per KPiX, 4 buckets per channel)
	bool                   chanFound[kpix_checking][1024]; // variable that gives true if a channel has entries
	bool                   kpixFound[kpix_checking]; // variable that gives true if a kpix at the index n (0<=n<24) was found
	uint                   x;
	uint                   range;
	uint                   value;
	uint                   kpix;
	uint                   channel;
	uint                   bucket;
	double                  tstamp;
	string                 serial;
	
	double 					y;
	
	uint 					subCount;
	double 					bunchClk;
	
	KpixSample::SampleType type;
	TH1F                   	*hist[kpix_checking][1024][bucket_checking];  // #entries/ADC histograms per channel, bucket, kpix and histogram
	TH1F			*hist_timed[kpix_checking][1024][bucket_checking]; //  #entries/time_of_event per channel, bucket, kpix and histogram
	TH1F			*channel_time[kpix_checking][1024][bucket_checking];
	TH1F			*channel_entries[kpix_checking][bucket_checking+1]; // ADC distribution Total number of events differed per bucket and kpix
	TH1F			*channel_entries_timed[kpix_checking][bucket_checking+1]; // Time distribution Total number of events differed per bucket and kpix
	
	TH1F			*strip_entries[kpix_checking][bucket_checking+1];
	TH1F			*timed_strip_entries[kpix_checking][bucket_checking+1];
	
	TH1F			*y_entries[kpix_checking][bucket_checking+1];
	TH1F			*timed_y_entries[kpix_checking][bucket_checking+1];

	
	TH1D			*trigger_difference[kpix_checking]; //Time difference to an external timestamp
	TH1F			*times_kpix[kpix_checking][bucket_checking+1];
	TH1F			*trig_count[kpix_checking][bucket_checking+1];
	
	
	cout << "DEBUG : 1" << endl;
	
	// Stringstream initialization for histogram naming
	stringstream           tmp;
	stringstream           tmp_units;
	
	// Stringstream initialization for folder naming
	
	stringstream			FolderName;
	
	ofstream               xml;
	ofstream               csv;
	uint                   acqCount; // acquisitionCount
	uint                   acqProcessed;
	string                 outRoot;
	TFile                  *rFile;
	stringstream           crossString;
	stringstream           crossStringCsv;
	XmlVariables           config;
	ofstream               debug;
	ofstream		channel_file_bad;
	ofstream		channel_file_bad_fit;
	ofstream		channel_file_noise;
	ofstream		channel_file_calib;
	ofstream 		channel_file_adc_mean;
	
	// Calibration slope, is filled when 
	double					calib_slope[24][1024] = {1}; //ADD buckets later.
	int						calibration_check = 0;
	
	
	unordered_map<uint, uint> kpix2strip_left;
	unordered_map<uint, uint> kpix2strip_right;
	kpix2strip_left = kpix_left();
	kpix2strip_right = kpix_right();
	
	pixel					pixel_kpix[1024];
	pixel_mapping(pixel_kpix);
	
	
	int 					num_of_channels[kpix_checking] = {0};
	
	std::vector<int>        monster_channels;
	
	
	const int               monster_finder_limit = 100;
	std::vector<int>        monster_cycles[kpix_checking];
	
	int monster_counter[kpix_checking] = {0}; // kpix monster counter
	
	//////////////////////////////////////////
	// File Read
	//////////////////////////////////////////
	
	// Data file is the first and only arg
	if ( argc != 2 && argc!= 3) {
	cout << "Usage: ./analysis data_file [skip_cycles_front] \n";
	return(1);
	}
	
	char* end;
	// skip first few cycles or read in calibration data file:
	if ( argc == 3 ) {
		//cout << "Even more debug " << strtol(argv[2], &end, 10) << endl;
		if (strtol(argv[2], &end, 10) != 0 )
		{
			skip_cycles_front = strtol(argv[2], &end, 10);
			cout<< " -- I am skipping first events: " << skip_cycles_front << endl;
			tmp.str("");
			tmp << argv[1] << ".printSkipped.txt";
			outtxt = tmp.str();
			f_skipped_cycles = fopen(outtxt.c_str(), "w");
		}
		else 
		{
			cout << " -- Reading " << argv[2] << " as calibration input file." << endl;
			skip_cycles_front = 0;
			TFile *calibration_file = TFile::Open(argv[2]);
			calibration_check = 1;
			loopdir(calibration_file, "calib_");
			for (unsigned int i = 0; i<calib_graphs.size(); ++i)
			{
				//cout << "Current key1 = " << cal_key->GetClassName() << endl;
				
				string calib_name         = calib_graphs[i]->GetName();
				
				size_t kpix_num_start     = calib_name.find("_k")+2;
				size_t channel_num_start  = calib_name.find("_c")+2;
				size_t kpix_num_length       = calib_name.length() - kpix_num_start;
				size_t channel_num_length    = calib_name.find("_b") - channel_num_start;
				
			    string channel_string = calib_name.substr(channel_num_start, channel_num_length);
			    string kpix_string = calib_name.substr(kpix_num_start, kpix_num_length);
			    
			    int kpix = stoi(kpix_string);
			    int channel = stoi(channel_string);
				
				//cout << "KPiX Number = " << kpix << endl;
				//cout << "Channel Number = " << channel << endl;
				
				calib_slope[kpix][channel] = calib_graphs[i]->GetFunction("pol1")->GetParameter(1);
				//cout << "Slope of KPiX " << kpix << " and channel " << channel << " is " <<  calib_slope[kpix][channel] << endl;
				
			}
		}
		
	
	}
	
	
	
	//////////////////////////////////////////
	// Open Data file
	//////////////////////////////////////////
	
	if ( ! dataRead.open(argv[1])  ) {
		cout << "Error opening data file " << argv[1] << endl;
		return(1);
	}
	
	// Create output names
	tmp.str("");
	tmp << argv[1] << ".root";
	outRoot = tmp.str();

	//////////////////////////////////////////
	// Read Data
	//////////////////////////////////////////
	cout << "Opened data file: " << argv[1] << endl;
	fileSize = dataRead.size();
	filePos  = dataRead.pos();
	
	// Init
	currPct          	= 0;
	lastPct          	= 100;
	
	cout << "\rReading File: 0 %" << flush;  // Printout of progress bar
	//goodTimes       	= 0;
	
	cout << "DEBUG : 2" << endl;
	
	// Open root file
	rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
	rFile->cd(); // move into root folder base
	FolderName.str("");
	FolderName << "General";
	rFile->mkdir(FolderName.str().c_str()); // produce a sub folder with name of variable FolderName
	TDirectory *General_folder = rFile->GetDirectory(FolderName.str().c_str()); // get path to subdirectory
	General_folder->cd(); // move into subdirectory
	
	
	//////////////////////////////////////////
	// Old histogram generation (needs to be cleaned up)
	//////////////////////////////////////////
	
	TH1F *channel_entries_total= new TH1F("Channel_entries_total", "Channel_Entries_total; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
	TH1F *channel_entries_total_timed= new TH1F("Channel_entries_total_timed", "Channel_entries_total_timed; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
	
	
	TH1F *time_kpix= new TH1F("time_kpix", "time_kpix; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
	TH1F *time_kpix_b[3];
	time_kpix_b[0]= new TH1F("time_kpix_b0", "time_kpix_b0; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
	time_kpix_b[1]= new TH1F("time_kpix_b1", "time_kpix_b1; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
	time_kpix_b[2]= new TH1F("time_kpix_b2", "time_kpix_b2; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
	time_kpix_b[3]= new TH1F("time_kpix_b3", "time_kpix_b3; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
	
	TH1F *time_external= new TH1F("time_external", "time_external; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8192.5); // one higher because the accuracy is higher
	
	TH2F *strip_vs_kpix = new TH2F ("strip_vs_kpix", "strip_vs_kpix; KPiX_Address;  strip_address", 1024,-0.5,1023.5,1024,-0.5,1023.5);
	
	
	TH1F *total= new TH1F("Total_response", "total_response; Charge (ADC); #entries/#acq.cycles", 300, -0.5, 8191.5);
	TH1F *total_DisConnect= new TH1F("Total_response_disconnected", "total_respons of disconnected kpix channels; Charge (ADC); #entries/#acq.cycles", 300, -0.5, 8191.5);
	TH1F *total_Connect= new TH1F("Total_response_connected", "total response of connected kpix channels; Charge (ADC); #entries/#acq.cycles", 300, -0.5, 8191.5);
	
	
	
	TH1F *total_timed= new TH1F("Total_response_timed", "total_response_timed; Charge (ADC); #entries/#acq.cycles", 300, -0.5, 8191.5);
	
	
	TH1D *beam_ext_time_diff = new TH1D("beam_ext_time_diff", "beam_ext_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 8001, -0.5, 1000.5);
	
	
	TH1F *ExtTrigPerCycle = new TH1F ("external_triggers_per_cycle", "ext_trig_per_acq.; #ext_triggers_per_acq.cycle; #entries/#acq.cycles",100,0.5,99.5);
	
	//TH1F *three_coincidence_channel_entries= new TH1F("three_coincidence_channel_entries", "three_coincidence_channel_entries; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
	TH1F *full_coincidence_channel_entries= new TH1F("full_coincidence_channel_entries", "full_coincidence_channel_entries; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
	
	TH1F *hit_diff_timed = new TH1F("hit_diff_timed", "hit_diff_timed; #strips; #entries/#acq.cycles", 200, -100.5, 99.5);
	
	
	while ( dataRead.next(&event) ) // event read to check for filled channels and kpix to reduce number of empty histograms.
	{
		acqCount++;
		
		int cycle_time_local[24][8192] = {0}; //list for each kpix that is part of improved monster finder
		//if ( acqCount < 10)
		//{
		//cout << "DEBUG: EVENTNUMBER " << event.eventNumber() << endl;
		//}
		if (acqCount > skip_cycles_front)
		{
			acqProcessed++;
			for (x=0; x < event.count(); x++)
			{
		
				//// Get sample
				sample  = event.sample(x);
				kpix    = sample->getKpixAddress();
				tstamp  = sample->getSampleTime();
				channel = sample->getKpixChannel();
				bucket  = sample->getKpixBucket();
				type    = sample->getSampleType();
				if ( type == KpixSample::Data )
				{
					kpixFound[kpix]          = true;
					chanFound[kpix][channel] = true;
					bucketFound[kpix][channel][bucket] = true;
					cycle_time_local[kpix][int(tstamp)]+=1;
					
				}
			}
		}
		else 
		{
			auto byte = event.count();
			auto train = event.eventNumber();
			if (f_skipped_cycles!=NULL)
			fprintf(f_skipped_cycles, " index = %d , byte = %6d, train = %6d \n ", acqCount, byte, train);
		
		}
		
		
		
	}
	
	if (f_skipped_cycles!=NULL)  {
		fclose( f_skipped_cycles);
		cout << endl;
		cout << "Wrote skipped cycles to " << outtxt << endl;
		cout << endl;
	}
	range = 0;
	
	cout << "DEBUG : 3" << endl;
	
	int cycle_checking;  // Program crashes when more than ~1700 cycles are checked maybe a memory issues, therefore the checking will have a maximum of 1000
	
	
	if (acqCount < 1000) cycle_checking = acqCount/10.0;
	else cycle_checking = 100;
	TH1F*					cycle_time[24][cycle_checking];
	TH1F*					cycle_time_ext[cycle_checking];
	FolderName.str("");
	FolderName << "Acquisition_Cycles";
	General_folder->mkdir(FolderName.str().c_str());
	TDirectory *gen_cycle_folder = General_folder->GetDirectory(FolderName.str().c_str());
	rFile->cd(gen_cycle_folder->GetPath());
	for (int cycles = 0; cycles < cycle_checking; cycles++) // produce subfolders per cycle
	{
		FolderName.str("");
		FolderName << "Cycle_" << cycles;
		gen_cycle_folder->mkdir(FolderName.str().c_str());
		TDirectory *cycles_folder = gen_cycle_folder->GetDirectory(FolderName.str().c_str());
		cycles_folder->cd();
		tmp.str("");
		tmp << "time_distribution_external" << "_evt_" << cycles;
		cycle_time_ext[cycles] = new TH1F(tmp.str().c_str(), "time_distribution_external; time [#bunch_clk_count]; #entries/#acq.cycles", 8192, -0.5, 8191.5);
	}
	
	
	
	
	dataRead.close();
	
	//double weight = 1.0/acqCount; //normalization weight  #entries*weight = #entries/acq.cycle
	double weight = 1.0/acqProcessed;
	
	//////////////////////////////////////////
	// New histogram generation within subfolder structure
	//////////////////////////////////////////
	
	for (kpix = 0; kpix < kpix_checking; kpix++) //looping through all possible kpix
	{
		//
		//cout << "DEBUG test " << kpixFound[kpix] << endl;
		if (kpixFound[kpix]) //checking if kpix exists
		{
			
			
			if ( kpix%2 == 0 )
			{
				
				rFile->cd(); //producing subfolder for kpix same as above for the event subfolder structure
				FolderName.str("");
				FolderName << "KPiX_left_" << kpix;
				rFile->mkdir(FolderName.str().c_str());
				TDirectory *kpix_folder = rFile->GetDirectory(FolderName.str().c_str());
				kpix_folder->cd();
				tmp.str("");
				tmp << "Left_Strip_entries_k" << kpix << "_total";
				strip_entries[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Strip_Entries; Strip_address; #entries/#acq.cycles", 920,-0.5, 919.5);
				tmp.str("");
				tmp << "timed_left_strip_entries_k" << kpix << "_total";
				timed_strip_entries[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Timed_Strip_Entries; Strip_address; #entries/#ext_signals", 920,-0.5, 919.5);
				
				tmp.str("");
				tmp << "y_entries_k" << kpix << "_total";
				y_entries[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "y_positional_Entries; y/#mum; #entries/#acq.cycles", 1840,-0.5, 92000.5);
				
				tmp.str("");
				tmp << "timed_y_entries_k" << kpix << "_total";
				timed_y_entries[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "timed_y_positional_Entries; y/#mum; #entries/#acq.cycles", 1840,-0.5, 92000.5);
				
				
				tmp.str("");
				tmp << "Channel_entries_k" << kpix << "_total";
				channel_entries[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
				tmp.str("");
				tmp << "Channel_entries_k" << kpix << "_total_timed";
				channel_entries_timed[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Channel_Entries_timed; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
				tmp.str("");
				tmp << "timestamp_kpix_k" << kpix  << "_total";
				times_kpix[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);
				tmp.str("");
				tmp << "trig_count_k" << kpix << "_total";
				trig_count[kpix][bucket_checking]  = new TH1F (tmp.str().c_str(), "trig_count;  #triggered channels; #entries/#acq.cycles",1024, -0.5, 1023.5);
				tmp.str("");
				tmp << "ext_time_diff_k" << kpix;
				trigger_difference[kpix] = new TH1D (tmp.str().c_str(), "intern_extern_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 8001, -0.5, 1000.5);

				
				for (bucket = 0; bucket < bucket_checking; bucket++)
				{
					FolderName.str("");
					FolderName << "bucket_" << bucket;
					kpix_folder->mkdir(FolderName.str().c_str());
					TDirectory *gen_buckets_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
					rFile->cd(gen_buckets_folder->GetPath());
					tmp.str("");
					tmp << "Channel_entries_k" << kpix << "_b" << bucket;
					channel_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
					tmp.str("");
					tmp << "left_strip_entries_k" << kpix << "_b" << bucket;
					strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Strip_Entries; Strip_address; #entries/#acq.cycles", 920,-0.5, 919.5);
					tmp.str("");
					tmp << "timed_left_strip_entries_k" << kpix << "_b" << bucket;
					timed_strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Timed_Strip_Entries; Strip_address; #entries/#acq.cycles", 920,-0.5, 919.5);
					
					
					tmp.str("");
					tmp << "y_entries_k" << kpix << "_b" << bucket;
					y_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "y_positional_Entries; y/#mum; #entries/#acq.cycles", 1840,-0.5, 92000.5);
					
					tmp.str("");
					tmp << "timed_y_entries_k" << kpix << "_b" << bucket;
					timed_y_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "timed_y_positional_Entries; y/#mum; #entries/#acq.cycles", 1840,-0.5, 92000.5);
					
					
					tmp.str("");
					tmp << "Channel_entries_k" << kpix <<  "_b" << bucket << "_timed";
					channel_entries_timed[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries_timed; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
					tmp.str("");
					tmp << "timestamp_kpix_k" << kpix << "_b" << bucket;
					times_kpix[kpix][bucket] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);
					tmp.str("");
					tmp << "trig_count_k" << kpix << "_b" << bucket ;
					trig_count[kpix][bucket]  = new TH1F(tmp.str().c_str(), "trig_count;  #triggered channels; #entries/#acq.cycles",1024, -0.5,1023.5);
					
					
				}
				FolderName.str("");
				FolderName << "cycles";
				kpix_folder->mkdir(FolderName.str().c_str());
				TDirectory *cycle_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
				rFile->cd(cycle_folder->GetPath());
				FolderName.str("");
				FolderName << "Strips_and_Channels";
				kpix_folder->mkdir(FolderName.str().c_str());
				TDirectory *channels_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
				rFile->cd(channels_folder->GetPath());
				for (channel = 0; channel < 1024; channel++)
				{
					if (chanFound[kpix][channel])
					{
						FolderName.str("");
						FolderName << "strip_" << kpix2strip_left.at(channel) << "_channel_" << channel;
						channels_folder->mkdir(FolderName.str().c_str());
						TDirectory *channel_folder = channels_folder->GetDirectory(FolderName.str().c_str());
						rFile->cd(channel_folder->GetPath());
		
						num_of_channels[kpix] = num_of_channels[kpix] + 1;
		
						for (bucket = 0; bucket < bucket_checking; bucket++)
						{
							if (bucketFound[kpix][channel][bucket])
							{
								// Naming of histograms and generating of histograms
		
								FolderName.str("");
								FolderName << "bucket_" << bucket;
								channel_folder->mkdir(FolderName.str().c_str());
								TDirectory *buckets_folder = channel_folder->GetDirectory(FolderName.str().c_str());
								rFile->cd(buckets_folder->GetPath());
								if (calibration_check == 1)
								{
									tmp.str("");  //set stringstream tmp to an empty string	
									tmp << "hist_fc" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp << "_b" << dec << bucket; // add _b$bucket
									tmp << "_k" << dec << kpix; // add _k$kpix to stringstream
								
									tmp_units.str(""); //set stringstream decribing histogram units to an empty string
									tmp_units << "hist_fc" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp_units << "_b" << dec << bucket; // add _b$bucket
									tmp_units << "_k" << dec << kpix; // add _k$kpix to stringstream
									tmp_units << "; Charge (fC); #entries/#acq.cycles"; // add title: x label, y label to stringstream
									hist[kpix][channel][bucket] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),2000, -0.5,499.5);	
									
									tmp.str("");
									tmp << "hist_timed_fc" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp << "_b" << dec << bucket;
									tmp << "_k" << dec << kpix;
			
									tmp_units.str("");
									tmp_units << "hist_timed_fc" << "_s" << dec << kpix2strip_left.at(channel);
									tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp_units << "_b" << dec << bucket;
									tmp_units << "_k" << dec << kpix;
									tmp_units << "; Charge (fC); #entries/#acq.cycles";
								
									hist_timed[kpix][channel][bucket] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),2000, -0.5,499.5);
									
									
								}
								else
								{
									tmp.str("");  //set stringstream tmp to an empty string
									
									tmp << "hist" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp << "_b" << dec << bucket; // add _b$bucket
									tmp << "_k" << dec << kpix; // add _k$kpix to stringstream
			
									tmp_units.str(""); //set stringstream decribing histogram units to an empty string
									tmp_units << "hist" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp_units << "_b" << dec << bucket; // add _b$bucket
									tmp_units << "_k" << dec << kpix; // add _k$kpix to stringstream
									tmp_units << "; Charge (ADC); #entries/#acq.cycles"; // add title: x label, y label to stringstream
									hist[kpix][channel][bucket] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);	
									
									tmp.str("");
									tmp << "hist_timed" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp << "_b" << dec << bucket;
									tmp << "_k" << dec << kpix;
			
									tmp_units.str("");
									tmp_units << "hist_timed" << "_s" << dec << kpix2strip_left.at(channel);
									tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp_units << "_b" << dec << bucket;
									tmp_units << "_k" << dec << kpix;
									tmp_units << "; Charge (ADC); #entries/#acq.cycles";
								
									hist_timed[kpix][channel][bucket] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);
									
									
								}
														
								
		
								
		
								tmp.str("");
								tmp << "time" << "_s" << dec << kpix2strip_left.at(channel);
								tmp << "_c" << dec << setw(4) << setfill('0') << channel;
								tmp << "_b" << dec << bucket;
								tmp << "_r" << dec << range;
								tmp << "_k" << dec << kpix;
		
								tmp_units.str("");
								tmp_units << "time" << "_s" << dec << kpix2strip_left.at(channel);
								tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
								tmp_units << "_b" << dec << bucket;
								tmp_units << "_k" << dec << kpix;
								tmp_units << "; Time/bunchClkCount; #entries/#acq.cycles";
		
								channel_time[kpix][channel][bucket] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);
		
		
								
								
							}
						}
					}
				}
				
				
			}
			if ( kpix%2 == 1 )
			{
				rFile->cd(); //producing subfolder for kpix same as above for the event subfolder structure
				FolderName.str("");
				FolderName << "KPiX_right_" << kpix;
				rFile->mkdir(FolderName.str().c_str());
				TDirectory *kpix_folder = rFile->GetDirectory(FolderName.str().c_str());
				kpix_folder->cd();
				
				
				tmp.str("");
				tmp << "Right_Strip_entries_k" << kpix << "_total";
				strip_entries[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Strip_Entries; Strip_address; #entries/#acq.cycles", 920, 919.5, 1839.5);
				
				tmp.str("");
				tmp << "timed_right_strip_entries_k" << kpix << "_total";
				timed_strip_entries[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Timed_Strip_Entries; Strip_address; #entries/#ext_signals", 920, 919.5, 1839.5);
				
				tmp.str("");
				tmp << "y_entries_k" << kpix << "_total";
				y_entries[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "y_positional_Entries; y/#mum; #entries/#acq.cycles", 1840,-0.5, 92000.5);
				
				tmp.str("");
				tmp << "timed_y_entries_k" << kpix << "_total";
				timed_y_entries[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "timed_y_positional_Entries; y/#mum; #entries/#acq.cycles", 1840,-0.5, 92000.5);
				
				
				tmp.str("");
				tmp << "Channel_entries_k" << kpix << "_total";
				channel_entries[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
				tmp.str("");
				tmp << "Channel_entries_k" << kpix << "_total_timed";
				channel_entries_timed[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Channel_Entries_timed; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
				tmp.str("");
				tmp << "timestamp_kpix_k" << kpix  << "_total";
				times_kpix[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);
				tmp.str("");
				tmp << "trig_count_k" << kpix << "_total";
				trig_count[kpix][bucket_checking]  = new TH1F (tmp.str().c_str(), "trig_count;  #triggered channels; #entries/#acq.cycles",1024, -0.5, 1023.5);
				tmp.str("");
				tmp << "ext_time_diff_k" << kpix;
				trigger_difference[kpix] = new TH1D (tmp.str().c_str(), "intern_extern_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 8001, -0.5, 1000.5);

				for (bucket = 0; bucket < bucket_checking; bucket++)
				{
					FolderName.str("");
					FolderName << "bucket_" << bucket;
					kpix_folder->mkdir(FolderName.str().c_str());
					TDirectory *gen_buckets_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
					rFile->cd(gen_buckets_folder->GetPath());
					tmp.str("");
					tmp << "Channel_entries_k" << kpix << "_b" << bucket;
					channel_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
					tmp.str("");
					tmp << "timed_right_strip_entries_k" << kpix << "_b" << bucket;
					timed_strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Timed_Strip_Entries; Strip_address; #entries/#acq.cycles", 920, 919.5, 1839.5);
					tmp.str("");
					tmp << "right_strip_entries_k" << kpix << "_b" << bucket;
					strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Strip_Entries; Strip_address; #entries/#acq.cycles", 920, 919.5, 1839.5);
					
					
					tmp.str("");
					tmp << "y_entries_k" << kpix << "_b" << bucket;
					y_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "y_positional_Entries; y/#mum; #entries/#acq.cycles", 1840,-0.5, 92000.5);
					
					tmp.str("");
					tmp << "timed_y_entries_k" << kpix << "_b" << bucket;
					timed_y_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "timed_y_positional_Entries; y/#mum; #entries/#acq.cycles", 1840,-0.5, 92000.5);
					
					
					tmp.str("");
					tmp << "Channel_entries_k" << kpix <<  "_b" << bucket << "_timed";
					channel_entries_timed[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries_timed; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
					tmp.str("");
					tmp << "timestamp_kpix_k" << kpix << "_b" << bucket;
					times_kpix[kpix][bucket] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);
					tmp.str("");
					tmp << "trig_count_k" << kpix << "_b" << bucket ;
					trig_count[kpix][bucket]  = new TH1F(tmp.str().c_str(), "trig_count;  #triggered channels; #entries/#acq.cycles",1024, -0.5,1023.5);
					
					
				}
				FolderName.str("");
				FolderName << "cycles";
				kpix_folder->mkdir(FolderName.str().c_str());
				TDirectory *cycle_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
				rFile->cd(cycle_folder->GetPath());
				FolderName.str("");
				FolderName << "Strips_and_Channels";
				kpix_folder->mkdir(FolderName.str().c_str());
				TDirectory *channels_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
				rFile->cd(channels_folder->GetPath());
				for (channel = 0; channel < 1024; channel++)
				{
					if (chanFound[kpix][channel])
					{
						FolderName.str("");
						FolderName << "strip_" << kpix2strip_left.at(channel) << "_channel_" << channel;
						channels_folder->mkdir(FolderName.str().c_str());
						TDirectory *channel_folder = channels_folder->GetDirectory(FolderName.str().c_str());
						rFile->cd(channel_folder->GetPath());
		
		
						num_of_channels[kpix] = num_of_channels[kpix] + 1;
		
						for (bucket = 0; bucket < bucket_checking; bucket++)
						{
							if (bucketFound[kpix][channel][bucket])
							{
								// Naming of histograms and generating of histograms
		
								FolderName.str("");
								FolderName << "bucket_" << bucket;
								channel_folder->mkdir(FolderName.str().c_str());
								TDirectory *buckets_folder = channel_folder->GetDirectory(FolderName.str().c_str());
								rFile->cd(buckets_folder->GetPath());
								if (calibration_check == 1)
								{
									tmp.str("");  //set stringstream tmp to an empty string	
									tmp << "hist_fc" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp << "_b" << dec << bucket; // add _b$bucket
									tmp << "_k" << dec << kpix; // add _k$kpix to stringstream
								
									tmp_units.str(""); //set stringstream decribing histogram units to an empty string
									tmp_units << "hist_fc" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp_units << "_b" << dec << bucket; // add _b$bucket
									tmp_units << "_k" << dec << kpix; // add _k$kpix to stringstream
									tmp_units << "; Charge (fC); #entries/#acq.cycles"; // add title: x label, y label to stringstream
									hist[kpix][channel][bucket] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),2000, -0.5,499.5);	
									
									tmp.str("");
									tmp << "hist_timed_fc" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp << "_b" << dec << bucket;
									tmp << "_k" << dec << kpix;
			
									tmp_units.str("");
									tmp_units << "hist_timed_fc" << "_s" << dec << kpix2strip_left.at(channel);
									tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp_units << "_b" << dec << bucket;
									tmp_units << "_k" << dec << kpix;
									tmp_units << "; Charge (fC); #entries/#acq.cycles";
								
									hist_timed[kpix][channel][bucket] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),2000, -0.5,499.5);
									
									
								}
								else
								{
									tmp.str("");  //set stringstream tmp to an empty string
									
									tmp << "hist" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp << "_b" << dec << bucket; // add _b$bucket
									tmp << "_k" << dec << kpix; // add _k$kpix to stringstream
			
									tmp_units.str(""); //set stringstream decribing histogram units to an empty string
									tmp_units << "hist" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp_units << "_b" << dec << bucket; // add _b$bucket
									tmp_units << "_k" << dec << kpix; // add _k$kpix to stringstream
									tmp_units << "; Charge (ADC); #entries/#acq.cycles"; // add title: x label, y label to stringstream
									hist[kpix][channel][bucket] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);	
									
									tmp.str("");
									tmp << "hist_timed" << "_s" << dec <<  kpix2strip_left.at(channel);
									tmp << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp << "_b" << dec << bucket;
									tmp << "_k" << dec << kpix;
			
									tmp_units.str("");
									tmp_units << "hist_timed" << "_s" << dec << kpix2strip_left.at(channel);
									tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
									tmp_units << "_b" << dec << bucket;
									tmp_units << "_k" << dec << kpix;
									tmp_units << "; Charge (ADC); #entries/#acq.cycles";
								
									hist_timed[kpix][channel][bucket] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);
									
									
								}
														
								
		
								
		
								tmp.str("");
								tmp << "time" << "_s" << dec << kpix2strip_left.at(channel);
								tmp << "_c" << dec << setw(4) << setfill('0') << channel;
								tmp << "_b" << dec << bucket;
								tmp << "_r" << dec << range;
								tmp << "_k" << dec << kpix;
		
								tmp_units.str("");
								tmp_units << "time" << "_s" << dec << kpix2strip_left.at(channel);
								tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
								tmp_units << "_b" << dec << bucket;
								tmp_units << "_k" << dec << kpix;
								tmp_units << "; Time/bunchClkCount; #entries/#acq.cycles";
		
								channel_time[kpix][channel][bucket] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);
		
		
								
								
							}
						}
					}
				}
				
				//rFile->mkdir(FolderName.str().c_str());
				//TDirectory *kpix_folder = rFile->GetDirectory(FolderName.str().c_str());
				//kpix_folder->cd();
			}
		}
	}
	
	cout << "DEBUG : 4" << endl;
	
	////////////////////////////////////////////
	//// Data read for all events for detailed look into single event structure
	////////////////////////////////////////////
	dataRead.open(argv[1]); //open binary file
	int cycle_num = 0;
	int cycle_num_ext = -1;
	int ext_trig_count = 0;
	
	while ( dataRead.next(&event) ) //loop through binary file event structure until end of file
	{
		cycle_num++;
		if ( cycle_num > skip_cycles_front)
		{
			double 	trigger_counter[24] = {0}; // fill the entire list of trigger_counter with 0
			vector<int> trigger_times[24];
			
			for (x=0; x < event.count(); x++)  //within the binary file go through each event
			{
					////// Get sample
				sample   = event.sample(x);  // check event subtructure
				kpix     = sample->getKpixAddress();
				channel  = sample->getKpixChannel();
				bucket   = sample->getKpixBucket();
				value    = sample->getSampleValue();
				type     = sample->getSampleType();
				tstamp   = sample->getSampleTime();
				range    = sample->getSampleRange();
				         
				bunchClk = sample->getBunchCount();
				subCount = sample->getSubCount();
				
				int strip = 9999;
				if (kpix%2 == 0) // if left kpix
				{
					strip = kpix2strip_left.at(channel);
				}
				else  // if right kpix
				{
					strip  = kpix2strip_right.at(channel);
				}
				y = yParameter(strip, kpix); //testbeam parameter
				
				if (kpix < kpix_checking && bucket < bucket_checking)
				{
					if (type == 2) //if type of event is ==2, the event is of type external timestamp
					{
						
						if (x == 0) cycle_num_ext++;
						double time = bunchClk + double(subCount * 0.125);
						ext_trig_count++;
					
					}
					//if (kpixFound[kpix])
					//{
						//cout << "KPiX found at " << kpix << endl;
					//}
					if ( type == 0 ) // If event is of type KPiX data
					{
						channel_entries[kpix][bucket]->Fill(channel, weight);
						channel_entries[kpix][bucket_checking]->Fill(channel, weight);
						
						
						y_entries[kpix][bucket]->Fill(y, weight);
						y_entries[kpix][bucket_checking]->Fill(y, weight);
						
						if (kpix % 2 == 0)
						{
							strip_entries[kpix][bucket]->Fill(kpix2strip_left.at(channel), weight);
							strip_entries[kpix][bucket_checking]->Fill(kpix2strip_left.at(channel), weight);
						}
						
						if (kpix % 2 == 1)
						{
							strip_entries[kpix][bucket]->Fill(kpix2strip_right.at(channel), weight);
							strip_entries[kpix][bucket_checking]->Fill(kpix2strip_right.at(channel), weight);
							
						}
						times_kpix[kpix][bucket]->Fill(tstamp, weight);
						times_kpix[kpix][bucket_checking]->Fill(tstamp, weight);
						trigger_counter[kpix] = trigger_counter[kpix] + (1.0/num_of_channels[kpix]);
						
						//cout << find(trigger_times[kpix].begin(), trigger_times[kpix].end(), tstamp) << endl;
						if (find(trigger_times[kpix].begin(), trigger_times[kpix].end(), tstamp) == trigger_times[kpix].end() || trigger_times[kpix].empty())
						{
							//cout <<"test"<< endl;
							trigger_times[kpix].push_back(tstamp);
							
						}
					}
				}	
			}
			
		}
	}
	
	dataRead.close(); // close file as we have looped through it and are now at the end
	dataRead.open(argv[1]); //open file again to start from the beginning
	
	int two_coincidence = 0;
	int three_coincidence = 0;
	int extern_trigger_id={0};
	
	cout << "DEBUG : 5" << endl;
	cycle_num = 0;
	
	
	while ( dataRead.next(&event) )
	{
		cycle_num++;
		if ( cycle_num > skip_cycles_front)
		{
			std::vector<double> time_ext;
			std::vector<int> channel_hits[24];
			std::vector<int> timestamp[24];
			std::vector<int> adc_value[24];
			std::vector<double> time_diff_kpix_ext[24];
			std::vector<int> AssignedTrigger[24];
			std::vector<pair<double, double>> time_coincident_hits[24];
			
			//std::vector<int> Assignment_number;
			int num_trig_count[24][5] = {0};
		
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
		
				bunchClk = sample->getBunchCount();
				subCount = sample->getSubCount();
		
				int strip = 9999;
				if (kpix%2 == 0) // if left kpix
				{
					strip = kpix2strip_left.at(channel);
				}
				else  // if right kpix
				{
					strip  = kpix2strip_right.at(channel);
				}
				y = yParameter(strip, kpix); //Testbeam parameter
		
				if (type == 2)// If event is of type external timestamp
				{
					double time = bunchClk + double(subCount * 0.125);
					time_external->Fill(time, weight);
					time_ext.push_back(time);
					//cout << "DEBUG: channel in timestmap = " << channel << endl;
					//cout << "DEBUG: bucket in timestmap = " << bucket << endl;
				}
		
				if (kpix < kpix_checking && bucket < bucket_checking)
				{
					if ( type == KpixSample::Data ) // If event is of type KPiX data
					{
						if (!(kpix == 30 && (channel == 500 || channel == 501 || channel == 490 || channel == 491 || channel == 522 || channel == 523 || channel == 532 || channel == 533 )))
						{
							channel_hits[kpix].push_back(channel);
							timestamp[kpix].push_back(tstamp);
							adc_value[kpix].push_back(value);
							time_kpix->Fill(tstamp, weight);
							time_kpix_b[bucket]->Fill(tstamp, weight);
								
							if (calibration_check == 1)
							{
								hist[kpix][channel][bucket]->Fill(double(value)/calib_slope[kpix][channel]/* *pow(10,15) */ , weight);
							}
							else
							{
								hist[kpix][channel][bucket]->Fill(value, weight);
							}
							
							channel_entries_total->Fill(channel, weight);
							channel_time[kpix][channel][bucket]->Fill(tstamp, weight);
							total->Fill(value, weight);
					
							strip_vs_kpix->Fill(channel, kpix2strip_left.at(channel));
							if ( kpix2strip_left.at(channel)!=9999 ){
								total_DisConnect->Fill(value, weight);
							}
							else {
								total_Connect->Fill(value, weight);
							}
							
							
							num_trig_count[kpix][bucket] += 1;
							num_trig_count[kpix][bucket_checking] += 1;
							
							double trig_diff = smallest_time_diff(time_ext, tstamp); //Calculation of minimal difference is done in a function for cleanup
							
							
							int assigned_number;
							if (time_ext.size() > 0) //only calculate the time difference between triggers if there are some external triggers
							{

								time_diff_kpix_ext[kpix].push_back(trig_diff);
								AssignedTrigger[kpix].push_back(assigned_number);
								beam_ext_time_diff->Fill(trig_diff, weight);
								trigger_difference[kpix]->Fill(trig_diff, weight);
								if((trig_diff >= 0.0 )  && (trig_diff  <= 3.0) )
								{
									if (calibration_check == 1)
									{
										hist_timed[kpix][channel][bucket]->Fill(double(value)/calib_slope[kpix][channel]/* *pow(10,15) */, weight);
									}
									else
									{
										hist_timed[kpix][channel][bucket]->Fill(value, weight);
									}
									total_timed->Fill(value, weight);
									channel_entries_total_timed->Fill(channel, weight);
									channel_entries_timed[kpix][bucket]->Fill(channel, weight);
									channel_entries_timed[kpix][bucket_checking]->Fill(channel, weight);
									
									timed_y_entries[kpix][bucket]->Fill(y, weight);
									timed_y_entries[kpix][bucket_checking]->Fill(y, weight);
									 
									if (kpix%2 == 0)
									{
										timed_strip_entries[kpix][bucket]->Fill(kpix2strip_left.at(channel), 1.0/ext_trig_count);
										timed_strip_entries[kpix][bucket_checking]->Fill(kpix2strip_left.at(channel), 1.0/ext_trig_count);
									}
									else
									{
										timed_strip_entries[kpix][bucket_checking]->Fill(kpix2strip_right.at(channel), 1.0/ext_trig_count);
										timed_strip_entries[kpix][bucket]->Fill(kpix2strip_right.at(channel), 1.0/ext_trig_count);
									}
									time_coincident_hits[kpix].push_back(make_pair(kpix2strip_left.at(channel), tstamp));
								}
								//cout << "DEBUG " << trig_diff << endl;
			
							}
						}
					}
				}
				//cout << "DEBUG time size" << time_ext.size() << endl;
			}
			ExtTrigPerCycle->Fill(time_ext.size());
			for (int kpix = 0; kpix < kpix_checking; ++kpix)
			{
				if (kpixFound[kpix])
				{
					for (int bucket = 0; bucket < bucket_checking; ++bucket)
					{
						trig_count[kpix][bucket]->Fill(num_trig_count[kpix][bucket], weight);
						trig_count[kpix][bucket_checking]->Fill(num_trig_count[kpix][bucket], weight);
						
					}
					
				}
			}	
		cout << "DEBUG : 6" << endl;
		
			//////////////////////////////////////////
			// Triggering efficiency and coincidence calculation, takes a lot of time.
			// -- removed, but can be found in analysisEcal.cxx file
			//////////////////////////////////////////
			
			
			
			
			extern_trigger_id=extern_trigger_id+time_ext.size();  // Counting which global external trigger was matched to a channel
			
			
			for (int v = 0; v< time_coincident_hits[17].size(); ++v)
			{
				for (int j = 0; j < time_coincident_hits[19].size(); j++)
				{
					if (time_coincident_hits[17].at(v).second ==  time_coincident_hits[19].at(j).second)
					{
						hit_diff_timed->Fill(time_coincident_hits[17].at(v).first - time_coincident_hits[19].at(j).first);
					}
					
				}
			}
			
		}	
		////   Show progress
		filePos  = dataRead.pos();
		currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
		if ( currPct != lastPct ) 
		{
			cout << "\rReading File: " << currPct << " %      " << flush;
			lastPct = currPct;
		}
		
	}
	
	for (int kpix = 0; kpix < kpix_checking; ++kpix)
	{
		if (kpixFound[kpix])
		{
			for (int channel = 0; channel < 1024; ++channel)
			{
				if (chanFound[kpix][channel])
				{
					for (int bucket = 0; bucket < bucket_checking; ++bucket)
					{
						if (bucketFound[kpix][channel][bucket])
						{
							 hist[kpix][channel][bucket]->GetMean();
							 //hist[kpix][channel][bucket]->Fit("gaus","q");
						}
						
					}
				}
			}
			
		}
	}
	
	
	
	cout <<  endl << "Full coincidence of sensors with external trigger: " << full_coincidence_channel_entries->GetEntries() << endl;
	cout << "Three coincidence of sensors: " << three_coincidence << endl;
	cout << "Two coincidence of sensors: " << two_coincidence << endl;
	
	cout << endl <<  "An event is currently classified as a monster if the amount of triggers one acquisition clock within a cycle is above " << monster_finder_limit << endl;
	cout << "_______________________________________________________________" << endl;
	
	ofstream emptybinfile;
	emptybinfile.open ("emptybinfile.txt");
	
  
	
	for (kpix = 0; kpix < 24; kpix++)
	{
		if (kpixFound[kpix])
		{
			cout << "Number of monster events in " << kpix << " = " << monster_counter[kpix] << endl;
			cout << "Number of normed monster events in " << kpix << " = " << monster_counter[kpix]*weight << endl;
			cout << "Number of entries in KPiX" << kpix << " = " << channel_entries[kpix][bucket_checking]->GetEntries() << endl << endl;
			cout << "Number of NORMED entries in KPiX" << kpix << " = " << (channel_entries[kpix][bucket_checking]->GetEntries())*weight << endl << endl;
			
			for (int bin = 1; bin < 1025; bin++)  // bin 0 is the underflow, therefore need to start counting at 1. Bin == channel_address+1
			{
				if (channel_entries[kpix][bucket_checking]->GetBinContent(bin) == 0 && emptybinfile.is_open()) 
				{
					emptybinfile << bin-1 << endl;
				}
				else if (!emptybinfile.is_open()) cout << "error" << endl;
			}
		}	
		
	}	
	cout << "Saved the empty bins to file /home/lycoris-dev/Desktop/emptybinfile.txt" << endl ;
	emptybinfile.close();
	cout << "DEBUG : 7" << endl;
	
	//for (int k = 0; k < 1024; k++)
	//{
	//cout << "DEBUG channel number: " << k << endl;
	//cout << "DEBUG X coord: " << pixel_kpix[k].x << endl;
	//cout << "DEBUG y coord: " << pixel_kpix[k].y << endl << endl;
	//}
	
	cout << endl;
	cout << "Writing root plots to " << outRoot << endl;
	cout << endl;
	
	rFile->Write();
	gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
	rFile->Close();
	
	
	
	dataRead.close();
	return(0);
}

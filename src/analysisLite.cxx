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
// 18/09/2019: fix bugs and remove usage of multi-dimentional arrays by <mengqing.wu@desy.de>
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
	
	string                 calState;
	uint                   lastPct;
	uint                   currPct;
	bool                   bucketFound[32][1024][4];  // variable that gives true if bucket has an entry (32 possible number of KPiX, 1024 channels per KPiX, 4 buckets per channel)
	bool                   chanFound[32][1024]; // variable that gives true if a channel has entries
	bool                   kpixFound[32]; // variable that gives true if a kpix at the index n (0<=n<32) was found
	uint                   x;
	//uint                   range;
	uint                   value;
	uint                   kpix;
	uint                   channel;
	uint                   bucket;
	uint 					subCount;
	double 					bunchClk;
	
	uint64_t				frameruntime;
	uint64_t				runtime;
	
	double                  tstamp;
	string                 serial;
	KpixSample::SampleType type;
	TH1F                   	*hist[32][1024][4][2];  // #entries/ADC histograms per channel, bucket, kpix and histogram
	TH1F			*hist_timed[32][1024][4][2]; //  #entries/time_of_event per channel, bucket, kpix and histogram
	TH1F			*channel_time[32][1024][4][2];
	TH1F			*channel_entries[32][5]; // ADC distribution Total number of events differed per bucket and kpix
	TH1F			*left_strip_entries[32][5];
	TH1F			*right_strip_entries[32][5];
	TH1F			*timed_left_strip_entries[32][5];
	TH1F			*timed_right_strip_entries[32][5];
	TH1F			*channel_entries_timed[32][5]; // Time distribution Total number of events differed per bucket and kpix
	TH1F			*channel_entries_no_strip[32][5]; 
	TH1F            *fc_response_b0[32];
	
	TH1D			*trigger_difference[32]; //Time difference to an external timestamp
	TH1D			*trigger_diff_connected[32]; //Time difference to an external timestamp
	TH1D			*trigger_diff_disconnected[32]; //Time difference to an external timestamp
	TH1F			*channel_entries_no_monster[32][5];
	TH1F			*times_kpix[32][5];
	TH1F			*times_kpix_monster[32][5];
	TH1F			*times_kpix_no_monster[32][5];
	TH1F			*trig_count[32][5];
	TH1F                    *hist_buck_sum[32][1024];
	
	TH1F 			*unique_times[32];
	
	TH2F			*kpix_entries_left[32][5];
	TH2F			*kpix_entries_right[32][5];
	
	
	
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
	double					calib_slope[32][1024] = {1}; //ADD buckets later.
	int						calibration_check = 0;
	
	
	unordered_map<uint, uint> kpix2strip_left;
	unordered_map<uint, uint> kpix2strip_right;
	kpix2strip_left = kpix_left();
	kpix2strip_right = kpix_right();
	
	pixel					pixel_kpix[1024];
	pixel_mapping(pixel_kpix);
	
	
	int 					num_of_channels[32] = {0};
	TH1F*					acq_num_ext[32];
	
	TH1F* 					AssignedChannelHist_Total[32];
	std::vector<int>        monster_channels;
	
	
	const int               monster_finder_limit = 100;
	std::vector<int>        monster_cycles[32];
	
	int monster_counter[32] = {0}; // kpix monster counter
	//for (int i = 0; i < 1000; i++)
	//{
	//	cout << "Strip mapping result, kpix #" << i << " is equal to strip #" << kpix2strip.at(i) << endl;
	//}
	
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
			    
			    int kpix_int = stoi(kpix_string);
			    int channel_int = stoi(channel_string);
				
				//cout << "KPiX Number = " << kpix << endl;
				//cout << "Channel Number = " << channel << endl;
				
				calib_slope[kpix_int][channel_int] = calib_graphs[i]->GetFunction("pol1")->GetParameter(1);
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

	// event read to check for filled channels and kpix to reduce number of empty histograms.
	while ( dataRead.next(&event) ) {
		acqCount++;
		
		int cycle_time_local[32][8192] = {0}; //list for each kpix that is part of improved monster finder
		if ( acqCount < 10){
			cout << "DEBUG: EVENTNUMBER " << event.eventNumber() << endl;
		}
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
				if ( type == KpixSample::Data ){
					kpixFound[kpix]          = true;
					chanFound[kpix][channel] = true;
					bucketFound[kpix][channel][bucket] = true;
					cycle_time_local[kpix][int(tstamp)]+=1;
				}
				//cout << "KPIX: " << kpix << endl;
				//cout << "Channel: " << channel << endl;
				
			}
		}
		else {
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
	
	
	int cycle_checking;  // Program crashes when more than ~1700 cycles are checked maybe a memory issues, therefore the checking will have a maximum of 1000
	if (acqCount < 1000) cycle_checking = acqCount/10.0;
	else cycle_checking = 1000;
	TH1F* 					AssignedChannelHist[32][cycle_checking];
	TH1F* 					trigger_difference_per_acq[32][cycle_checking];
	TH1F*					cycle_time[32][cycle_checking];
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
	for (kpix = 0; kpix<32; kpix++){
		if (!kpixFound[kpix]) continue; // no such kpix
		
		rFile->cd();
		FolderName.str("");
		FolderName << "KPiX_" << kpix;
		rFile->mkdir(FolderName.str().c_str());
		TDirectory *kpix_folder = rFile->GetDirectory(FolderName.str().c_str());
		kpix_folder->cd();

		tmp.str("");
		tmp << "trig_count_k" << kpix << "_total";
		trig_count[kpix][4]  = new TH1F (tmp.str().c_str(), "trig_count;  #triggered channels; #entries/#acq.cycles",1024, -0.5, 1023.5);

		tmp.str("");
		tmp << "ext_time_diff_k" << kpix;
		trigger_difference[kpix] = new TH1D (tmp.str().c_str(), "intern_extern_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 8001, -0.5, 1000.5);

		tmp.str("");
		tmp << "ext_time_diff_con_k" << kpix;
		trigger_diff_connected[kpix] = new TH1D (tmp.str().c_str(), "intern_extern_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 8001, -0.5, 1000.5);

		tmp.str("");
		tmp << "ext_time_diff_discon_k" << kpix;
		trigger_diff_disconnected[kpix] = new TH1D (tmp.str().c_str(), "intern_extern_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 8001, -0.5, 1000.5);

		
		//-- START General Kpix plot at Bucket==0 ONLY
		bucket = 0 ;
		tmp.str("");
		tmp << "timestamp_kpix_k" << kpix << "_b" << bucket;
		times_kpix[kpix][bucket] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);

		tmp.str("");
		tmp << "timed_left_strip_entries_k" << kpix << "_b" << bucket;
		timed_left_strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Timed_Strip_Entries; Strip_address; #entries/#acq.cycles", 920,-0.5, 919.5);
		tmp.str("");
		tmp << "timed_right_strip_entries_k" << kpix << "_b" << bucket;
		timed_right_strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Timed_Strip_Entries; Strip_address; #entries/#acq.cycles", 920, 919.5, 1839.5);
		
		//-- END General Kpix plot at Bucket==0 ONLY
		

		FolderName.str("");
		FolderName << "Strips_and_Channels";
		kpix_folder->mkdir(FolderName.str().c_str());
		TDirectory *channels_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
		rFile->cd(channels_folder->GetPath());
		
		for (channel = 0; channel < 1024; channel++){
			if (!chanFound[kpix][channel]) continue;

			// to keep consistent with previous code: but hard coded left right kpix for S59
			FolderName.str("");
			if (kpix == 26)
				FolderName << "strip_" << kpix2strip_left.at(channel) << "_channel_" << channel;
			else if (kpix == 28)
				FolderName << "strip_" << kpix2strip_right.at(channel) << "_channel_" << channel;
			else
				FolderName << "channel_" << channel;
			channels_folder->mkdir(FolderName.str().c_str());
			TDirectory *channel_folder = channels_folder->GetDirectory(FolderName.str().c_str());
			rFile->cd(channel_folder->GetPath());

			for (bucket = 0; bucket < 4; bucket++){
				if (bucket!=0) continue;


				tmp.str("");
				tmp << "Channel_entries_k" << kpix << "_b" << bucket;
				channel_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
				
				tmp.str("");
				tmp << "Channel_entries_k" << kpix <<  "_b" << bucket << "_timed";
				channel_entries_timed[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries_timed; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);

				
				
				
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
				hist[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);	
				
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
				hist_timed[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);
				
			}// histograms per bucket [NB! but only bucket 0 defined]
			
			
		}// histograms per channel
		
	}// histograms per kpix
	
	
	////////////////////////////////////////////
	//// Data read for all events for detailed look into single event structure
	////////////////////////////////////////////
	dataRead.open(argv[1]); //open binary file
	int cycle_num = 0;
	int cycle_num_ext = -1;

	//loop through binary file event structure until end of file
	while ( dataRead.next(&event) ) {
		cycle_num++;
		if ( cycle_num > skip_cycles_front){
			double 	trigger_counter[32] = {0}; // fill the entire list of trigger_counter with 0
			vector<int> trigger_times[32];

			//if (kpixFound[26]) acq_num_ext[26]->Fill(trigger_counter[26]); // trigger counting for monster check
			//if (kpixFound[28]) acq_num_ext[28]->Fill(trigger_counter[28]);
			//if (kpixFound[30]) acq_num_ext[30]->Fill(trigger_counter[30]);
			
			for (x=0; x< event.count(); x++) {
				sample = event.sample(x);
				kpix    = sample->getKpixAddress();
				channel = sample->getKpixChannel();
				bucket  = sample->getKpixBucket();
				value   = sample->getSampleValue();
				type    = sample->getSampleType();
				tstamp  = sample->getSampleTime();

				if (bucket !=0 ) continue; // skip all non bucket0 events
				
				if (type == KpixSample::Timestamp) {
					
					if (x == 0) cycle_num_ext++;
					double time = tstamp + double(value * 0.125);
					
				}// case1: find Timestamp sample

				if (type == KpixSample::Data){
					channel_entries[kpix][bucket]->Fill(channel, weight);
					times_kpix[kpix][bucket]->Fill(tstamp, weight);
					
					
				}// case2: find Data sample
		  
				
				
			}// for-loop over all kpix samples
			
		} // if-end, skip cycles
	} // loop over all cycles
	

	dataRead.close(); // close file as we have looped through it and are now at the end
	dataRead.open(argv[1]); //open file again to start from the beginning
	
	int two_coincidence = 0;
	int three_coincidence = 0;
	int extern_trigger_id={0};
	
	cycle_num = 0;
	
	
	while ( dataRead.next(&event) )
		{
		cycle_num++;
		if ( cycle_num > skip_cycles_front)
		{
			std::vector<double> time_ext;
			std::vector<int> channel_hits[32];
			std::vector<int> timestamp[32];
			std::vector<int> adc_value[32];
			std::vector<double> time_diff_kpix_ext[32];
			std::vector<int> AssignedTrigger[32];
			std::vector<pair<double, double>> time_coincident_hits[32];
			
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
				if (bucket !=0 ) continue; // skip all non bucket0 events

				if (type == KpixSample::Timestamp){
					
					double time = tstamp + double(value * 0.125);
					
					time_external->Fill(time, weight);
					time_ext.push_back(time);
					//cout << "DEBUG: channel in timestmap = " << channel << endl;
					//cout << "DEBUG: bucket in timestmap = " << bucket << endl;
				}

				if (type == KpixSample::Data) {
					// NB: hard code to process only strip data
					//	if (kpix != 26 || kpix != 28) continue;

					hist[kpix][channel][bucket][0]->Fill(value, weight);
					
					if ( time_ext.size() > 0 ){
					double trig_diff = smallest_time_diff(time_ext, tstamp); //Calculation of minimal time difference 

					bool isconn = true;
					if (kpix == 26) isconn = (kpix2strip_left.at(channel) == 9999 );
					if (kpix == 28) isconn = (kpix2strip_right.at(channel) == 9999 );

					trigger_difference[kpix]->Fill(trig_diff, weight);
					if (isconn)
						trigger_diff_connected[kpix]->Fill(trig_diff,weight);
					else
						trigger_diff_disconnected[kpix]->Fill(trig_diff,weight);
					
					if ( trig_diff >= 0.0 && trig_diff <= 3.0 ){
						hist_timed[kpix][channel][bucket][0] -> Fill(value, weight);
						if (kpix == 26)
							timed_left_strip_entries[kpix][bucket] -> Fill(kpix2strip_left.at(channel), weight);
							
						if (kpix == 28)
							timed_right_strip_entries[kpix][bucket] -> Fill(kpix2strip_right.at(channel), weight);
					}
					
												
					}// if--end, algorithm with External Timestamp
				}
				
		
			}
		
		
			//////////////////////////////////////////
			// Triggering efficiency and coincidence calculation, takes a lot of time.
			// -- removed, but can be found in analysisEcal.cxx file
			//////////////////////////////////////////
			
			
			extern_trigger_id=extern_trigger_id+time_ext.size();  // Counting which global external trigger was matched to a channel
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
	
		
	cout <<  endl << "Full coincidence of sensors with external trigger: " << full_coincidence_channel_entries->GetEntries() << endl;
	cout << "Three coincidence of sensors: " << three_coincidence << endl;
	cout << "Two coincidence of sensors: " << two_coincidence << endl;
	
	cout << endl <<  "An event is currently classified as a monster if the amount of triggers one acquisition clock within a cycle is above " << monster_finder_limit << endl;
	cout << "_______________________________________________________________" << endl;
	

	
	// for (int k = 0; k < 1024; k++)
	// {
	// cout << "DEBUG channel number: " << k << endl;
	// cout << "DEBUG X coord: " << pixel_kpix[k].x << endl;
	// cout << "DEBUG y coord: " << pixel_kpix[k].y << endl << endl;
	// }
	
	cout << endl;
	cout << "Writing root plots to " << outRoot << endl;
	cout << endl;
	
	rFile->Write();
	gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
	rFile->Close();
	
	
	
	dataRead.close();
	return(0);
}

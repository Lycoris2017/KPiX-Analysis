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
using namespace std;


//////////////////////////////////////////
// Global Variables
//////////////////////////////////////////
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

double smallest_time_diff( vector<double> ext_list, int int_value)
{
    double trigger_diff = 8200.0;
    for (uint k = 0; k<ext_list.size(); ++k)
    {
        double delta_t = int_value-ext_list.at(k);
        if (fabs(trigger_diff) > fabs(delta_t))// && delta_t > 0)
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
	
//	TH1::SetDefaultSumw2();
	
//	//////////////////////////////////////////
//	// Class declaration and histogram initialization
//	//////////////////////////////////////////
	
//	DataRead               dataRead;  //kpix event classes used for analysis of binary date
//	off_t                  fileSize;  //
//	off_t                  filePos;   //
//	KpixEvent              event;    //
//	KpixSample             *sample;   //
	
//	// cycles to skip in front:
//	long int                   skip_cycles_front;
//	FILE*                  f_skipped_cycles;
//	string                 outtxt;
	
//	string                 calState;
//	uint                   lastPct;
//	uint                   currPct;
//	bool                   bucketFound[32][1024][4];  // variable that gives true if bucket has an entry (32 possible number of KPiX, 1024 channels per KPiX, 4 buckets per channel)
//	bool                   chanFound[32][1024]; // variable that gives true if a channel has entries
//	bool                   kpixFound[32]; // variable that gives true if a kpix at the index n (0<=n<32) was found
//	uint                   x;
//	uint                   range;
//	uint                   value;
//	uint                   kpix;
//	uint                   channel;
//	uint                   bucket;
//	uint 					subCount;
//	double 					bunchClk;
	
//	uint64_t				frameruntime;
//	uint64_t				runtime;
	
//	double                  tstamp;
//	string                 serial;
//	KpixSample::SampleType type;
//	TH1F                   	*hist[32][1024][4][2];  // #entries/ADC histograms per channel, bucket, kpix and histogram
//	TH1F			*hist_timed[32][1024][4][2]; //  #entries/time_of_event per channel, bucket, kpix and histogram
//	TH1F			*channel_time[32][1024][4][2];
//	TH1F			*channel_entries[32][5]; // ADC distribution Total number of events differed per bucket and kpix
//	TH1F			*left_strip_entries[32][5];
//	TH1F			*right_strip_entries[32][5];
//	TH1F			*timed_left_strip_entries[32][5];
//	TH1F			*timed_right_strip_entries[32][5];
//	TH1F			*channel_entries_timed[32][5]; // Time distribution Total number of events differed per bucket and kpix
//	TH1F			*channel_entries_no_strip[32][5];
//	TH1F            *fc_response_b0[32];
	
//	TH1D			*trigger_difference[32]; //Time difference to an external timestamp
//	TH1D			*trigger_diff_connected[32]; //Time difference to an external timestamp
//	TH1D			*trigger_diff_disconnected[32]; //Time difference to an external timestamp
//	TH1F			*channel_entries_no_monster[32][5];
//	TH1F			*times_kpix[32][5];
//	TH1F			*times_kpix_monster[32][5];
//	TH1F			*times_kpix_no_monster[32][5];
//	TH1F			*trig_count[32][5];
//	TH1F                    *hist_buck_sum[32][1024];
	
//	TH1F 			*unique_times[32];
	
//	TH2F			*kpix_entries_left[32][5];
//	TH2F			*kpix_entries_right[32][5];
	
	
	
//	// Stringstream initialization for histogram naming
//	stringstream           tmp;
//	stringstream           tmp_units;
	
//	// Stringstream initialization for folder naming
	
//	stringstream			FolderName;
	
//	ofstream               xml;
//	ofstream               csv;
//	uint                   acqCount; // acquisitionCount
//	uint                   acqProcessed;
//	string                 outRoot;
//	TFile                  *rFile;
//	stringstream           crossString;
//	stringstream           crossStringCsv;
//	XmlVariables           config;
//	ofstream               debug;
//	ofstream		channel_file_bad;
//	ofstream		channel_file_bad_fit;
//	ofstream		channel_file_noise;
//	ofstream		channel_file_calib;
//	ofstream 		channel_file_adc_mean;
	
//	// Calibration slope, is filled when
//	double					calib_slope[32][1024] = {1}; //ADD buckets later.
//	int						calibration_check = 0;
	
	
//	unordered_map<uint, uint> kpix2strip_left;
//	unordered_map<uint, uint> kpix2strip_right;
//	kpix2strip_left = kpix_left();
//	kpix2strip_right = kpix_right();
	
//	pixel					pixel_kpix[1024];
//	pixel_mapping(pixel_kpix);
	
	
//	int 					num_of_channels[32] = {0};
//	TH1F*					acq_num_ext[32];
	
//	TH1F* 					AssignedChannelHist_Total[32];
//	std::vector<int>        monster_channels;
	
	
//	const int               monster_finder_limit = 100;
//	std::vector<int>        monster_cycles[32];
	
//	int monster_counter[32] = {0}; // kpix monster counter
//	//for (int i = 0; i < 1000; i++)
//	//{
//	//	cout << "Strip mapping result, kpix #" << i << " is equal to strip #" << kpix2strip.at(i) << endl;
//	//}
	
//	//////////////////////////////////////////
//	// File Read
//	//////////////////////////////////////////
	
//	// Data file is the first and only arg
//	if ( argc != 2 && argc!= 3) {
//	cout << "Usage: ./analysis data_file [skip_cycles_front] \n";
//	return(1);
//	}
	
//	char* end;
//	// skip first few cycles or read in calibration data file:
//    if ( argc == 3 ) {
//        //cout << "Even more debug " << strtol(argv[2], &end, 10) << endl;
//        if (strtol(argv[2], &end, 10) != 0 )
//        {
//            skip_cycles_front = strtol(argv[2], &end, 10);
//            cout<< " -- I am skipping first events: " << skip_cycles_front << endl;
//            tmp.str("");
//            tmp << argv[1] << ".printSkipped.txt";
//            outtxt = tmp.str();
//            f_skipped_cycles = fopen(outtxt.c_str(), "w");
//        }
//        else
//        {
//            cout << " -- Reading " << argv[2] << " as calibration input file." << endl;
//            skip_cycles_front = 0;
//            TFile *calibration_file = TFile::Open(argv[2]);
//            calibration_check = 1;

////			cout << "Number of calibration slopes = " << calib_graphs.size() << endl;
//            TTree *calib_tree = (TTree*)calibration_file->Get("calibration_tree");
//            uint kpix_calib, channel_calib, bucket_calib, range_calib;
//            double slope_calib, slope_error_calib, pearsson_calib;
//            calib_tree->SetBranchAddress("kpix", &kpix_calib);
//            calib_tree->SetBranchAddress("channel", &channel_calib);
//            calib_tree->SetBranchAddress("bucket", &bucket_calib);
//            calib_tree->SetBranchAddress("range", &range_calib);
//            calib_tree->SetBranchAddress("calib_slope", &slope_calib);
//            calib_tree->SetBranchAddress("calib_error", &slope_error_calib);
//            calib_tree->SetBranchAddress("calib_pearsson", &pearsson_calib);
//            long int nEnTrees = calib_tree->GetEntries();
//            for (long int i = 0; i < nEnTrees; ++i){
//                calib_tree->GetEntry(i);
//                index = keybit(kpix_calib, channel_calib, bucket_calib);
//                std::pair<double,double> calib_entry = std::make_pair(slope_calib, pearsson_calib);
//                calibs.emplace(index,calib_entry);
//                calib_slope[kpix_calib][channel_calib][bucket_calib] = slope_calib;
////                calib_pearsson[kpix_calib][channel_calib][bucket_calib] = pearsson_calib;
//            }

////			cout << " -- Reading " << argv[3] << " as pedestal subtraction input file." << endl;
////			skip_cycles_front = 0;
////			TFile *pedestal_file = TFile::Open(argv[3]);

////			TTree* pedestal_tree = (TTree*)pedestal_file->Get("pedestal_tree");

////            int kpix_pedestal, channel_pedestal, bucket_pedestal;
////			double pedestal_median, pedestal_MAD;

////			pedestal_tree->SetBranchAddress("pedestal_median", &pedestal_median);
////            pedestal_tree->SetBranchAddress("kpix_num", &kpix_pedestal);
////            pedestal_tree->SetBranchAddress("channel_num", &channel_pedestal);
////            pedestal_tree->SetBranchAddress("bucket_num", &bucket_pedestal);
////			pedestal_tree->SetBranchAddress("pedestal_MAD", &pedestal_MAD);

////			long int nentries = pedestal_tree->GetEntries();
////			for (long int i = 0; i < nentries; ++i)
////			{
////				pedestal_tree->GetEntry(i);
////                pedestal_MedMAD[kpix_pedestal][channel_pedestal][bucket_pedestal][0] = pedestal_median;
////                pedestal_MedMAD[kpix_pedestal][channel_pedestal][bucket_pedestal][1] = pedestal_MAD;
////				//if (pedestal_MedMAD[kpix_num][channel_num][bucket_num][1] != 0)
////				//{
////					//cout << "Pedestal MAD " << pedestal_MedMAD[kpix_num][channel_num][bucket_num][1] << " kpix " << kpix_num << " channel " << channel_num << " bucket " << bucket_num << endl;
////				//}
////			}



//            gROOT->GetListOfFiles()->Remove(calibration_file);
//            calibration_file->Close();
////			gROOT->GetListOfFiles()->Remove(pedestal_file);
////			pedestal_file->Close();

//        }


//    }
	
	
	
//	//////////////////////////////////////////
//	// Open Data file
//	//////////////////////////////////////////
	
//	if ( ! dataRead.open(argv[1])  ) {
//		cout << "Error opening data file " << argv[1] << endl;
//		return(1);
//	}
	
//	// Create output names
//	tmp.str("");
//	tmp << argv[1] << ".root";
//	outRoot = tmp.str();
	
//	//================================================ Charge calibration, currently not working
//	//Start of change
//	//if (argc == 4)
//	//{
//		//TFile *calibration_file = TFile::Open(argv[3]);
//		//TIter cal_iter(calibration_file->GetListOfKeys());
//		//TKey *cal_key;
//		//int count = 0;
//		//while ((cal_key = (TKey*)cal_iter()))
//		//{
//			////cout << "Current key1 = " << cal_key->GetClassName() << endl;
//			//string title = cal_key->GetTitle();
//			//if ((strcmp(cal_key->GetClassName(), "TGraphErrors") == 0) && (title.find("calib") != std::string::npos ))
//			//{
//				////cout << "test" << endl;
//				//TGraphErrors *calib_graph = (TGraphErrors*)cal_key->ReadObj();
//				//calib_graph->Fit("pol1","eq","",fitMin[range],fitMax[range]);
//				//calib_slope[count] = calib_graph->GetFunction("pol1")->GetParameter(1);
//				//calib_y0[count] = calib_graph->GetFunction("pol1")->GetParameter(0);
//				////cout << "Title: " << cal_key->GetTitle() << " with Slope = " << calib_slope[count] << endl;
//				//++count;
//			//}
//		//}
//	//}
//	//end of change
	
//	//////////////////////////////////////////
//	// Read Data
//	//////////////////////////////////////////
//	cout << "Opened data file: " << argv[1] << endl;
//	fileSize = dataRead.size();
//	filePos  = dataRead.pos();
	
//	// Init
//	currPct          	= 0;
//	lastPct          	= 100;
	
//	cout << "\rReading File: 0 %" << flush;  // Printout of progress bar
//	//goodTimes       	= 0;
	
	
//	// Open root file
//	rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
//	rFile->cd(); // move into root folder base
//	FolderName.str("");
//	FolderName << "General";
//	rFile->mkdir(FolderName.str().c_str()); // produce a sub folder with name of variable FolderName
//	TDirectory *General_folder = rFile->GetDirectory(FolderName.str().c_str()); // get path to subdirectory
//	General_folder->cd(); // move into subdirectory
	
	
//	//////////////////////////////////////////
//	// Old histogram generation (needs to be cleaned up)
//	//////////////////////////////////////////
	
//	TH1F *channel_entries_total= new TH1F("Channel_entries_total", "Channel_Entries_total; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
//	TH1F *channel_entries_total_timed= new TH1F("Channel_entries_total_timed", "Channel_entries_total_timed; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
	
	
//	TH1F *time_kpix= new TH1F("time_kpix", "time_kpix; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
//	TH1F *time_kpix_b[3];
//	time_kpix_b[0]= new TH1F("time_kpix_b0", "time_kpix_b0; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
//	time_kpix_b[1]= new TH1F("time_kpix_b1", "time_kpix_b1; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
//	time_kpix_b[2]= new TH1F("time_kpix_b2", "time_kpix_b2; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
//	time_kpix_b[3]= new TH1F("time_kpix_b3", "time_kpix_b3; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8191.5);
	
//	TH1F *time_external= new TH1F("time_external", "time_external; Time/bunchClkCount; #entries/#acq.cycles", 300, -0.5, 8192.5); // one higher because the accuracy is higher
	
//	TH2F *strip_vs_kpix = new TH2F ("strip_vs_kpix", "strip_vs_kpix; KPiX_Address;  strip_address", 1024,-0.5,1023.5,1024,-0.5,1023.5);
	
	
//	TH1F *total= new TH1F("Total_response", "total_response; Charge (ADC); #entries/#acq.cycles", 300, -0.5, 8191.5);
//	TH1F *total_DisConnect= new TH1F("Total_response_disconnected", "total_respons of disconnected kpix channels; Charge (ADC); #entries/#acq.cycles", 300, -0.5, 8191.5);
//	TH1F *total_Connect= new TH1F("Total_response_connected", "total response of connected kpix channels; Charge (ADC); #entries/#acq.cycles", 300, -0.5, 8191.5);
	
	
	
//	TH1F *total_timed= new TH1F("Total_response_timed", "total_response_timed; Charge (ADC); #entries/#acq.cycles", 300, -0.5, 8191.5);
	
	
//	TH1D *beam_ext_time_diff = new TH1D("beam_ext_time_diff", "beam_ext_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 8001, -0.5, 1000.5);
	
	
//	TH1F *ExtTrigPerCycle = new TH1F ("external_triggers_per_cycle", "ext_trig_per_acq.; #ext_triggers_per_acq.cycle; #entries/#acq.cycles",100,0.5,99.5);
	
//	//TH1F *three_coincidence_channel_entries= new TH1F("three_coincidence_channel_entries", "three_coincidence_channel_entries; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
//	TH1F *full_coincidence_channel_entries= new TH1F("full_coincidence_channel_entries", "full_coincidence_channel_entries; KPiX_channel_address; #entries/#acq.cycles", 1024, -0.5, 1023.5);
	
//	TH1F *hit_diff_timed = new TH1F("hit_diff_timed", "hit_diff_timed; #strips; #entries/#acq.cycles", 200, -100.5, 99.5);
	
	
//	while ( dataRead.next(&event) ) // event read to check for filled channels and kpix to reduce number of empty histograms.
//	{
//		acqCount++;
//		//cout << "Event size " << event.size() << endl;
//		//cout << "Event count " << event.count() << endl;
//		//cout << "Timestamp of Event #" << event.eventNumber() << " = " << event.timestamp() << endl;
//		//cout << "Timestamp of Event #" << acqCount << " = " << event.timestamp() << endl;
		
//		int cycle_time_local[32][8192] = {0}; //list for each kpix that is part of improved monster finder
//		//if ( acqCount < 10)
//		//{
//		//cout << "DEBUG: EVENTNUMBER " << event.eventNumber() << endl;
//		//}
//		if (acqCount > skip_cycles_front)
//		{
//			acqProcessed++;
//			for (x=0; x < event.count(); x++)
//			{
		
//				//// Get sample
//				sample  = event.sample(x);
//				kpix    = sample->getKpixAddress();
//				tstamp  = sample->getSampleTime();
//				channel = sample->getKpixChannel();
//				bucket  = sample->getKpixBucket();
//				type    = sample->getSampleType();
//				if ( type == KpixSample::Data )
//				{
//					kpixFound[kpix]          = true;
//					chanFound[kpix][channel] = true;
//					bucketFound[kpix][channel][bucket] = true;
//					cycle_time_local[kpix][int(tstamp)]+=1;
					
//				}
//				//cout << "KPIX: " << kpix << endl;
//				//cout << "Channel: " << channel << endl;
				
//			}
//		}
//		else
//		{
//			auto byte = event.count();
//			auto train = event.eventNumber();
//			if (f_skipped_cycles!=NULL)
//			fprintf(f_skipped_cycles, " index = %d , byte = %6d, train = %6d \n ", acqCount, byte, train);
		
//		}
		
//		//for (int kpix = 0; kpix<32; kpix++){
//		  //if (kpixFound[kpix]) {
//		    //int monster_check = 0;
//		    //for (int i = 0; i < 8192; ++i){
//		      //if (cycle_time_local[kpix][i] > monster_finder_limit) {
			
//				//monster_check = 1; // if among the 8192 possible time points there exists one that has more than monster_finder_limit triggers the event is classified as a monster
//				//break; // no reason to keep looking, we found a monster.
//		      //}
//		    //}
//		    ////if (monster_check){
//		      ////cout << " There is a monster in house #" << kpix <<" (kpix slot number) with more than " << monster_finder_limit << " (finder limit) eyes  under bed #" << acqCount << " (cycle number)" << endl;
//		      ////monster_cycles[kpix].push_back(acqCount);
//		      ////monster_counter[kpix]++;
//		      ////}
//		  //}
		  
//		//}

		
//	}
	
//	if (f_skipped_cycles!=NULL)  {
//		fclose( f_skipped_cycles);
//		cout << endl;
//		cout << "Wrote skipped cycles to " << outtxt << endl;
//		cout << endl;
//	}
//	range = 0;
	
	
//	int cycle_checking;  // Program crashes when more than ~1700 cycles are checked maybe a memory issues, therefore the checking will have a maximum of 1000
//	if (acqCount < 1000) cycle_checking = acqCount/10.0;
//	else cycle_checking = 1000;
//	TH1F* 					AssignedChannelHist[32][cycle_checking];
//	TH1F* 					trigger_difference_per_acq[32][cycle_checking];
//	TH1F*					cycle_time[32][cycle_checking];
//	TH1F*					cycle_time_ext[cycle_checking];
//	FolderName.str("");
//	FolderName << "Acquisition_Cycles";
//	General_folder->mkdir(FolderName.str().c_str());
//	TDirectory *gen_cycle_folder = General_folder->GetDirectory(FolderName.str().c_str());
//	rFile->cd(gen_cycle_folder->GetPath());
//	for (int cycles = 0; cycles < cycle_checking; cycles++) // produce subfolders per cycle
//	{
//		FolderName.str("");
//		FolderName << "Cycle_" << cycles;
//		gen_cycle_folder->mkdir(FolderName.str().c_str());
//		TDirectory *cycles_folder = gen_cycle_folder->GetDirectory(FolderName.str().c_str());
//		cycles_folder->cd();
//		tmp.str("");
//		tmp << "time_distribution_external" << "_evt_" << cycles;
//		cycle_time_ext[cycles] = new TH1F(tmp.str().c_str(), "time_distribution_external; time [#bunch_clk_count]; #entries/#acq.cycles", 8192, -0.5, 8191.5);
//	}
	
	
	
	
//	dataRead.close();
	
//	//double weight = 1.0/acqCount; //normalization weight  #entries*weight = #entries/acq.cycle
//	double weight = 1.0/acqProcessed;
	
//	//////////////////////////////////////////
//	// New histogram generation within subfolder structure
//	//////////////////////////////////////////

//	for (kpix = 0; kpix < 32; kpix++) //looping through all possible kpix
//	{
//		//
//		//cout << "DEBUG test " << kpixFound[kpix] << endl;
//		if (kpixFound[kpix]) //checking if kpix exists
//		{
//			rFile->cd(); //producing subfolder for kpix same as above for the event subfolder structure
//			FolderName.str("");
//			FolderName << "KPiX_" << kpix;
//			rFile->mkdir(FolderName.str().c_str());
//			TDirectory *kpix_folder = rFile->GetDirectory(FolderName.str().c_str());
//			kpix_folder->cd();
//			//tmp.str("");
//			//tmp << "Channel_entries_k" << kpix << "_total";
//			//channel_entries[kpix][4] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
//			//tmp.str("");
//			//tmp << "Channel_entries_k" << kpix << "_total_no_monster";
//			//channel_entries_no_monster[kpix][4] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
//			//tmp.str("");
//			//tmp << "Channel_entries_k" << kpix << "_total_timed";
//			//channel_entries_timed[kpix][4] = new TH1F(tmp.str().c_str(), "Channel_Entries_timed; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
//			//tmp.str("");
//			//tmp << "Channel_entries_k" << kpix << "_total_no_strip";
//			//channel_entries_no_strip[kpix][4] = new TH1F(tmp.str().c_str(), "Channel_Entries_no_strip; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
//			//tmp.str("");
//			//tmp << "Left_Strip_entries_k" << kpix << "_total";
//			//left_strip_entries[kpix][4] = new TH1F(tmp.str().c_str(), "Strip_Entries; Strip_address; #entries/#acq.cycles", 920,-0.5, 919.5);
//			//tmp.str("");
//			//tmp << "Right_Strip_entries_k" << kpix << "_total";
//			//right_strip_entries[kpix][4] = new TH1F(tmp.str().c_str(), "Strip_Entries; Strip_address; #entries/#acq.cycles", 920, 919.5, 1839.5);
//			//tmp.str("");
//			//tmp << "timed_left_strip_entries_k" << kpix << "_total";
//			//timed_left_strip_entries[kpix][4] = new TH1F(tmp.str().c_str(), "Timed_Strip_Entries; Strip_address; #entries/#acq.cycles", 920,-0.5, 919.5);
//			//tmp.str("");
//			//tmp << "timed_right_strip_entries_k" << kpix << "_total";
//			//timed_right_strip_entries[kpix][4] = new TH1F(tmp.str().c_str(), "Timed_Strip_Entries; Strip_address; #entries/#acq.cycles", 920, 919.5, 1839.5);
//			tmp.str("");
//			tmp << "timestamp_kpix_k" << kpix  << "_total";
//			times_kpix[kpix][4] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);
//			tmp.str("");
//			tmp << "timestamp_kpix_k" << kpix  << "_monster_total";
//			times_kpix_monster[kpix][4] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);
//			tmp.str("");
//			tmp << "timestamp_kpix_k" << kpix  << "_no_monster_total";
//			times_kpix_no_monster[kpix][4] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);
//			tmp.str("");
//			tmp << "trig_count_k" << kpix << "_total";
//			trig_count[kpix][4]  = new TH1F (tmp.str().c_str(), "trig_count;  #triggered channels; #entries/#acq.cycles",1024, -0.5, 1023.5);
//			tmp.str("");
//			tmp << "ext_time_diff_k" << kpix;
//			trigger_difference[kpix] = new TH1D (tmp.str().c_str(), "intern_extern_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 8001, -0.5, 1000.5);
//			tmp.str("");
//			tmp << "ext_time_diff_con_k" << kpix;
//			trigger_diff_connected[kpix] = new TH1D (tmp.str().c_str(), "intern_extern_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 8001, -0.5, 1000.5);
//			tmp.str("");
//			tmp << "ext_time_diff_discon_k" << kpix;
//			trigger_diff_disconnected[kpix] = new TH1D (tmp.str().c_str(), "intern_extern_time_diff; #Delta T (BunchClkCount); #entries/#acq.cycles", 8001, -0.5, 1000.5);
//			tmp.str("");
//			tmp << "assigned_channels_k" << kpix << "_total";
//			AssignedChannelHist_Total[kpix] = new TH1F (tmp.str().c_str(), "assigned_channels_per_ext_trig;   #assigned_channels; #entries/#acq.cycles",40, -0.5, 39.5);
	
//			tmp.str("");
//			tmp << "acq_num_ext_k" << kpix;
//			acq_num_ext[kpix] = new TH1F(tmp.str().c_str(), "acq_num_ext; #triggers/acq._cycle; #entries/#acq.cycles",5, -0.5, 4.5);
			
//			tmp.str("");
//			tmp << "kpix_entries_left_k" << kpix << "_total";
//			kpix_entries_left[kpix][4] = new TH2F ("kpix_entries_left", "kpix_entries_left; column;  row", 32,-0.5,31.5, 32,-0.5,31.5);
			
//			tmp.str("");
//			tmp << "kpix_entries_right_k" << kpix << "_total";
//			kpix_entries_right[kpix][4] = new TH2F ("kpix_entries_right", "kpix_entries_right; column;  row", 32,-0.5,31.5, 32,-0.5,31.5);
			
//			tmp.str("");
//			tmp << "unique_times_k" << kpix;
//			unique_times[kpix] = new TH1F(tmp.str().c_str(), "unique_times; #triggers; #entries/#acq.cycles", 101,-0.5, 99.5);
			
			
//			for (bucket = 0; bucket< 4; bucket++)
//			{
//				FolderName.str("");
//				FolderName << "bucket_" << bucket;
//				kpix_folder->mkdir(FolderName.str().c_str());
//				TDirectory *gen_buckets_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
//				rFile->cd(gen_buckets_folder->GetPath());
//				tmp.str("");
//				tmp << "Channel_entries_k" << kpix << "_b" << bucket;
//				channel_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
//				tmp.str("");
//				tmp << "Channel_entries_k" << kpix << "_b" << bucket << "_no_strip";
//				channel_entries_no_strip[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries_no_strip; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
//				tmp.str("");
//				tmp << "left_strip_entries_k" << kpix << "_b" << bucket;
//				left_strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Strip_Entries; Strip_address; #entries/#acq.cycles", 920,-0.5, 919.5);
//				tmp.str("");
//				tmp << "timed_left_strip_entries_k" << kpix << "_b" << bucket;
//				timed_left_strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Timed_Strip_Entries; Strip_address; #entries/#acq.cycles", 920,-0.5, 919.5);
//				tmp.str("");
//				tmp << "timed_right_strip_entries_k" << kpix << "_b" << bucket;
//				timed_right_strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Timed_Strip_Entries; Strip_address; #entries/#acq.cycles", 920, 919.5, 1839.5);
//				tmp.str("");
//				tmp << "right_strip_entries_k" << kpix << "_b" << bucket;
//				right_strip_entries[kpix][bucket] = new TH1F(tmp.str().c_str(), "Strip_Entries; Strip_address; #entries/#acq.cycles", 920, 919.5, 1839.5);
				
				
				
//				tmp.str("");
//				tmp << "Channel_entries_k" << kpix <<  "_b" << bucket << "_timed";
//				channel_entries_timed[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel_Entries_timed; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
//				tmp.str("");
//				tmp << "Channel_entries_k" << kpix << "_b" << bucket << "_no_monster";
//				channel_entries_no_monster[kpix][bucket] = new TH1F(tmp.str().c_str(), "Channel Entries; KPiX_channel_address; #entries/#acq.cycles", 1024,-0.5, 1023.5);
//				tmp.str("");
//				tmp << "timestamp_kpix_k" << kpix << "_b" << bucket;
//				times_kpix[kpix][bucket] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);
//				tmp.str("");
//				tmp << "timestamp_kpix_monster_k" << kpix << "_b" << bucket;
//				times_kpix_monster[kpix][bucket] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);
//				tmp.str("");
//				tmp << "trig_count_k" << kpix << "_b" << bucket ;
//				trig_count[kpix][bucket]  = new TH1F(tmp.str().c_str(), "trig_count;  #triggered channels; #entries/#acq.cycles",1024, -0.5,1023.5);
//				tmp.str("");
//				tmp << "timestamp_kpix_k" << kpix << "_b" << bucket << "_no_monster";
//				times_kpix_no_monster[kpix][bucket] = new TH1F(tmp.str().c_str(), "timestamp_kpix; time [#bunch_clk_count]; #entries/#acq.cycles", 8192,-0.5, 8191.5);
				

//				tmp.str("");
//				tmp << "kpix_entries_left_k" << kpix << "_bucket_" << bucket;
//				kpix_entries_left[kpix][bucket] = new TH2F ("kpix_entries_left", "kpix_entries_left; column;  row", 32,-0.5,31.5, 32,-0.5,31.5);
				
//				tmp.str("");
//				tmp << "kpix_entries_right_k" << kpix << "_bucket_" << bucket;
//				kpix_entries_right[kpix][bucket] = new TH2F ("kpix_entries_right", "kpix_entries_right; column;  row", 32,-0.5,31.5, 32,-0.5,31.5);
				
				
//			}
//			FolderName.str("");
//			FolderName << "cycles";
//			kpix_folder->mkdir(FolderName.str().c_str());
//			TDirectory *cycle_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
//			rFile->cd(cycle_folder->GetPath());
//			for (int cycles = 0; cycles < cycle_checking; cycles++)
//			{
//				FolderName.str("");
//				FolderName << "cycle_" << cycles;
//				cycle_folder->mkdir(FolderName.str().c_str());
//				TDirectory *cycles_folder = cycle_folder->GetDirectory(FolderName.str().c_str());
//				cycles_folder->cd();
//				tmp.str("");
//				tmp << "time_distribution_k" << kpix << "_evt_" << cycles;
//				cycle_time[kpix][cycles] = new TH1F(tmp.str().c_str(), "time_distribution; time [#bunch_clk_count]; #entries", 8192, -0.5, 8191.5);
//				tmp.str("");
//				tmp << "assigned_channels_k" << kpix << "_evt_" << cycles;
//				AssignedChannelHist[kpix][cycles]  = new TH1F (tmp.str().c_str(), "assigned_channels_per_ext_trig;  external_trigger_number; #assigned_channels ",100, -0.5, 99.5);
//				tmp.str("");
//				tmp << "trigger_difference_k" << kpix << "_evt_" << cycles;
//				trigger_difference_per_acq[kpix][cycles]  = new TH1F (tmp.str().c_str(), "trigger_difference;  #entries/#acq.cycles; #Delta T [BunchClkCount] ",8001, -0.5, 1000.5);
//				//tmp.str("");
//				//tmp << "assigned_number_k" << kpix << "_evt_" << cycles;
//				//AssignedNumberHist[kpix][cycles]  = new TH1F (tmp.str().c_str(), "assigned_NumberOfChannel_per_ext_trig;  #same_assignement; #entries/#acq.cycles",40,0,40);
//			}
//			FolderName.str("");
//			FolderName << "Strips_and_Channels";
//			kpix_folder->mkdir(FolderName.str().c_str());
//			TDirectory *channels_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
//			rFile->cd(channels_folder->GetPath());
//			for (channel = 0; channel < 1024; channel++)
//			{
//				if (chanFound[kpix][channel])
//				{
//					FolderName.str("");
//					FolderName << "strip_" << kpix2strip_left.at(channel) << "_channel_" << channel;
//					channels_folder->mkdir(FolderName.str().c_str());
//					TDirectory *channel_folder = channels_folder->GetDirectory(FolderName.str().c_str());
//					rFile->cd(channel_folder->GetPath());
	
//					tmp.str("");
//					tmp << "hist" << "_s" << dec << kpix2strip_left.at(channel);
//					tmp << "_c" << dec << setw(4) << setfill('0') << channel;
//					tmp << "_k" << dec << kpix;
	
//					tmp_units.str("");
//					tmp_units << "hist" << "_s" << dec <<  kpix2strip_left.at(channel);
//					tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
//					tmp_units << "_k" << dec << kpix;
//					tmp_units << "; Charge (ADC); #entries/#acq.cycles";
	
//					hist_buck_sum[kpix][channel] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5, 8191.5);
	
//					num_of_channels[kpix] = num_of_channels[kpix] + 1;
	
//					for (bucket = 0; bucket < 4; bucket++)
//					{
//						if (bucketFound[kpix][channel][bucket])
//						{
//							// Naming of histograms and generating of histograms
	
//							FolderName.str("");
//							FolderName << "bucket_" << bucket;
//							channel_folder->mkdir(FolderName.str().c_str());
//							TDirectory *buckets_folder = channel_folder->GetDirectory(FolderName.str().c_str());
//							rFile->cd(buckets_folder->GetPath());
//							if (calibration_check == 1)
//							{
//								tmp.str("");  //set stringstream tmp to an empty string
//								tmp << "hist_fc" << "_s" << dec <<  kpix2strip_left.at(channel);
//								tmp << "_c" << dec << setw(4) << setfill('0') << channel;
//								tmp << "_b" << dec << bucket; // add _b$bucket
//								tmp << "_k" << dec << kpix; // add _k$kpix to stringstream
							
//								tmp_units.str(""); //set stringstream decribing histogram units to an empty string
//								tmp_units << "hist_fc" << "_s" << dec <<  kpix2strip_left.at(channel);
//								tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
//								tmp_units << "_b" << dec << bucket; // add _b$bucket
//								tmp_units << "_k" << dec << kpix; // add _k$kpix to stringstream
//								tmp_units << "; Charge (fC); #entries/#acq.cycles"; // add title: x label, y label to stringstream
//								hist[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),2000, -0.5,499.5);
								
//								tmp.str("");
//								tmp << "hist_timed_fc" << "_s" << dec <<  kpix2strip_left.at(channel);
//								tmp << "_c" << dec << setw(4) << setfill('0') << channel;
//								tmp << "_b" << dec << bucket;
//								tmp << "_k" << dec << kpix;
		
//								tmp_units.str("");
//								tmp_units << "hist_timed_fc" << "_s" << dec << kpix2strip_left.at(channel);
//								tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
//								tmp_units << "_b" << dec << bucket;
//								tmp_units << "_k" << dec << kpix;
//								tmp_units << "; Charge (fC); #entries/#acq.cycles";
							
//								hist_timed[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),2000, -0.5,499.5);
								
								
//							}
//							else
//							{
//								tmp.str("");  //set stringstream tmp to an empty string
								
//								tmp << "hist" << "_s" << dec <<  kpix2strip_left.at(channel);
//								tmp << "_c" << dec << setw(4) << setfill('0') << channel;
//								tmp << "_b" << dec << bucket; // add _b$bucket
//								tmp << "_k" << dec << kpix; // add _k$kpix to stringstream
		
//								tmp_units.str(""); //set stringstream decribing histogram units to an empty string
//								tmp_units << "hist" << "_s" << dec <<  kpix2strip_left.at(channel);
//								tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
//								tmp_units << "_b" << dec << bucket; // add _b$bucket
//								tmp_units << "_k" << dec << kpix; // add _k$kpix to stringstream
//								tmp_units << "; Charge (ADC); #entries/#acq.cycles"; // add title: x label, y label to stringstream
//								hist[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);
								
//								tmp.str("");
//								tmp << "hist_timed" << "_s" << dec <<  kpix2strip_left.at(channel);
//								tmp << "_c" << dec << setw(4) << setfill('0') << channel;
//								tmp << "_b" << dec << bucket;
//								tmp << "_k" << dec << kpix;
		
//								tmp_units.str("");
//								tmp_units << "hist_timed" << "_s" << dec << kpix2strip_left.at(channel);
//								tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
//								tmp_units << "_b" << dec << bucket;
//								tmp_units << "_k" << dec << kpix;
//								tmp_units << "; Charge (ADC); #entries/#acq.cycles";
							
//								hist_timed[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);
								
								
//							}
													
							
	
							
	
//							tmp.str("");
//							tmp << "time" << "_s" << dec << kpix2strip_left.at(channel);
//							tmp << "_c" << dec << setw(4) << setfill('0') << channel;
//							tmp << "_b" << dec << bucket;
//							tmp << "_r" << dec << range;
//							tmp << "_k" << dec << kpix;
	
//							tmp_units.str("");
//							tmp_units << "time" << "_s" << dec << kpix2strip_left.at(channel);
//							tmp_units << "_c" << dec << setw(4) << setfill('0') << channel;
//							tmp_units << "_b" << dec << bucket;
//							tmp_units << "_k" << dec << kpix;
//							tmp_units << "; Time/bunchClkCount; #entries/#acq.cycles";
	
//							//channel_time[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);
							
							
//						}
//					}
//				}
//			}
//		}
//	}
	
//	////////////////////////////////////////////
//	//// Data read for all events for detailed look into single event structure
//	////////////////////////////////////////////
//	dataRead.open(argv[1]); //open binary file
//	int cycle_num = 0;
//	int cycle_num_ext = -1;
	
//	while ( dataRead.next(&event) ) //loop through binary file event structure until end of file
//	{
//		cycle_num++;
//		if ( cycle_num > skip_cycles_front)
//		{
//			double 	trigger_counter[32] = {0}; // fill the entire list of trigger_counter with 0
//			vector<int> trigger_times[32];
			
//			for (x=0; x < event.count(); x++)  //within the binary file go through each event
//			{
//					////// Get sample
//				sample  = event.sample(x);  // check event subtructure
//				kpix    = sample->getKpixAddress();
//				channel = sample->getKpixChannel();
//				bucket  = sample->getKpixBucket();
//				value   = sample->getSampleValue();
//				type    = sample->getSampleType();
//				tstamp  = sample->getSampleTime();
//				range   = sample->getSampleRange();
//				bunchClk = sample->getBunchCount();
//				subCount = sample->getSubCount();
//				if (type == 2) //if type of event is ==2, the event is of type external timestamp
//				{
//					//cout << cycle_num << endl;
//					//cout << cycle_num_ext << endl;
//					if (x == 0) cycle_num_ext++;
//					double time = bunchClk + double(subCount * 0.125);
//					if (cycle_num_ext < cycle_checking)
//					{
//						//cout << cycle_num_ext <<  " " << time << endl;
//						cycle_time_ext[cycle_num_ext]->Fill(time);
//					}
//				}
//				if ( type == 0 ) // If event is of type KPiX data
//				{
//					//cout << "Bucket " << bucket << "KPiX " << kpix << endl;
//					channel_entries[kpix][bucket]->Fill(channel, weight);
//					channel_entries[kpix][4]->Fill(channel, weight);
					
//					kpix_entries_left[kpix][4]->Fill(31-(channel/32), channel%32, weight);
//					kpix_entries_left[kpix][bucket]->Fill(32-(channel/32), channel%32, weight);
					
//					kpix_entries_right[kpix][4]->Fill(channel/32, 31-channel%32, weight);
//					kpix_entries_right[kpix][bucket]->Fill(channel/32, 31-channel%32, weight);
//					//if (channel != 395 && channel != 375 && channel != 429 && channel != 333) //DIRTY FIX FOR PRC, GET RID OF IT AGAIN
//					//{
//						left_strip_entries[kpix][bucket]->Fill(kpix2strip_left.at(channel), weight);
//						left_strip_entries[kpix][4]->Fill(kpix2strip_left.at(channel), weight);
//					//}
					
//					right_strip_entries[kpix][bucket]->Fill(kpix2strip_right.at(channel), weight);
//					right_strip_entries[kpix][4]->Fill(kpix2strip_right.at(channel), weight);

					
//					times_kpix[kpix][bucket]->Fill(tstamp, weight);
//					times_kpix[kpix][4]->Fill(tstamp, weight);
//					trigger_counter[kpix] = trigger_counter[kpix] + (1.0/num_of_channels[kpix]);
//					if (kpix2strip_left.at(channel) == 9999)
//					{
//						channel_entries_no_strip[kpix][bucket]->Fill(channel,weight);
//						channel_entries_no_strip[kpix][4]->Fill(channel,weight);
//					}
//					if (find(monster_cycles[kpix].begin(), monster_cycles[kpix].end(), event.eventNumber()) != monster_cycles[kpix].end())
//					{
//						times_kpix_monster[kpix][bucket]->Fill(tstamp, weight);
//						times_kpix_monster[kpix][4]->Fill(tstamp, weight);
//					}
//					else
//					{
//						channel_entries_no_monster[kpix][bucket]->Fill(channel, weight);
//						channel_entries_no_monster[kpix][4]->Fill(channel, weight);
//						times_kpix_no_monster[kpix][bucket]->Fill(tstamp, weight);
//						times_kpix_no_monster[kpix][4]->Fill(tstamp, weight);
//					}
//					if (cycle_num < cycle_checking)
//					{
//						cycle_time[kpix][cycle_num]->Fill(tstamp);
//					}
					
//					//cout << find(trigger_times[kpix].begin(), trigger_times[kpix].end(), tstamp) << endl;
//					if (find(trigger_times[kpix].begin(), trigger_times[kpix].end(), tstamp) == trigger_times[kpix].end() || trigger_times[kpix].empty())
//					{
//						//cout <<"test"<< endl;
//						trigger_times[kpix].push_back(tstamp);
						
//					}
					
//				}
//			}
	
		
//			if (trigger_counter[26] > 4) cout << trigger_counter[26] << endl;
//			 if (kpixFound[26]) acq_num_ext[26]->Fill(trigger_counter[26]); // trigger counting for monster check
//			 if (kpixFound[28]) acq_num_ext[28]->Fill(trigger_counter[28]);
//			 if (kpixFound[30]) acq_num_ext[30]->Fill(trigger_counter[30]);
//			cout << trigger_times[26].size() << endl;
//			for (int i = 0; i<32; ++i)
//			{
//				if (kpixFound[i])
//				{
//					 unique_times[i]->Fill(trigger_times[i].size());
//					 //for (int j =0; j< trigger_times[i].size(); ++j)
//					 //{
//						 //cout << "Timestamp = " << trigger_times[i].at(j) << endl;
//					 //}
//					 //cout << endl;
//				}
//			}
			
//		}
//	}
	
//	dataRead.close(); // close file as we have looped through it and are now at the end
//	dataRead.open(argv[1]); //open file again to start from the beginning
	
//	int two_coincidence = 0;
//	int three_coincidence = 0;
//	int extern_trigger_id={0};
	
	
//	cycle_num = 0;
	
	
//	while ( dataRead.next(&event) )
//		{
//		cycle_num++;
//		if ( cycle_num > skip_cycles_front)
//		{
//			std::vector<double> time_ext;
//			std::vector<int> channel_hits[32];
//			std::vector<int> timestamp[32];
//			std::vector<int> adc_value[32];
//			std::vector<double> time_diff_kpix_ext[32];
//			std::vector<int> AssignedTrigger[32];
//			std::vector<pair<double, double>> time_coincident_hits[32];
			
//			//std::vector<int> Assignment_number;
//			int num_trig_count[32][5] = {0};
		
//			//cout << " NEW EVENT " << endl;
//			for (x=0; x < event.count(); x++)
//			{
//				//cout << "DEBUG: EVENT COUNT " << event.count() << endl;
//				//// Get sample
//				sample  = event.sample(x);
//				kpix    = sample->getKpixAddress();
//				channel = sample->getKpixChannel();
//				bucket  = sample->getKpixBucket();
//				value   = sample->getSampleValue();
//				type    = sample->getSampleType();
//				tstamp  = sample->getSampleTime();
//				range   = sample->getSampleRange();
//				bunchClk = sample->getBunchCount();
//				subCount = sample->getSubCount();
		
		
		
		
//				if (type == 2)// If event is of type external timestamp
//				{
//					double time = bunchClk + double(subCount * 0.125);
//					time_external->Fill(time, weight);
//					time_ext.push_back(time);
//					//cout << "DEBUG: channel in timestmap = " << channel << endl;
//					//cout << "DEBUG: bucket in timestmap = " << bucket << endl;
//				}
		
		
//                if ( type == KpixSample::Data ) // If event is of type KPiX data
//                {
//                    if (!(kpix == 30 && (channel == 500 || channel == 501 || channel == 490 || channel == 491 || channel == 522 || channel == 523 || channel == 532 || channel == 533 )))
//                    {
//                        channel_hits[kpix].push_back(channel);
//                        timestamp[kpix].push_back(tstamp);
//                        adc_value[kpix].push_back(value);
//                        time_kpix->Fill(tstamp, weight);
//                        time_kpix_b[bucket]->Fill(tstamp, weight);
							
//                        std::vector<double> trig_diff_list;
//                        trig_diff_list.push_back(1);
//                        trig_diff_list.push_back(5);
//                        trig_diff_list.push_back(10);
//                        trig_diff_list.push_back(0.5);
//                        if (calibration_check == 1)
//                        {
//                            hist[kpix][channel][bucket][0]->Fill(double(value)/calib_slope[kpix][channel]/* *pow(10,15) */ , weight);
							
							
//                        }
//                        else
//                        {
//                            hist[kpix][channel][bucket][0]->Fill(value, weight);
							
//                        }
						
//                        hist_buck_sum[kpix][channel]->Fill(value,weight);
//                        channel_entries_total->Fill(channel, weight);
//                        channel_time[kpix][channel][bucket][0]->Fill(tstamp, weight);
				
//                        total->Fill(value, weight);
				
//                        strip_vs_kpix->Fill(channel, kpix2strip_left.at(channel));
//                        if ( kpix2strip_left.at(channel)!=9999 ){
//                            total_DisConnect->Fill(value, weight);
//                        }
//                        else {
//                            total_Connect->Fill(value, weight);
//                        }
						
						
//                        num_trig_count[kpix][bucket] += 1;
//                        num_trig_count[kpix][4] += 1;
						
//                        double trig_diff = smallest_time_diff(time_ext, tstamp); //Calculation of minimal difference is done in a function for cleanup
						
//                        int assigned_number;
//                        if (time_ext.size() > 0) //only calculate the time difference between triggers if there are some external triggers
//                        {
//                            //for (unsigned int j = 0; j < time_ext.size(); ++j)
//                            //{
//                                //trig_diff_list.push_back(tstamp-time_ext.at(j));
//                                //double delta_t = tstamp-time_ext.at(j);
//                                //if (fabs(trig_diff) > fabs(delta_t))
//                                //{
//                                    //trig_diff = tstamp-time_ext.at(j);
//                                    //assigned_number = j;
//                                //}
//                                ////else
//                                ////{
//                                    ////if (x<10)
//                                    ////{
//                                        ////cout << "Difference not lower than before" << endl;
//                                        ////cout << "Channel time stamp = " << tstamp << endl;
//                                        ////cout << "External match = " << time_ext.at(j) << endl;
//                                        ////cout << "Old difference = " << trig_diff << endl;
//                                        ////cout << "New difference = " << tstamp-time_ext.at(j) << endl;
//                                    ////}
							
//                                ////}
							
//                            //}
//                            //cout << "Trig diff old method = " <<  trig_diff << endl;
//                            //if (trig_diff_list.size() > 0) cout << "Trig diff new method = " <<  *std::min_element(trig_diff_list.begin(), trig_diff_list.end()) << endl; //seg fault when vector is empty
//                            //if (trig_diff_list.size() > 0) cout << "Trig diff new method position in vector = " <<  distance(trig_diff_list.begin(), min_element(trig_diff_list.begin(), trig_diff_list.end())) << endl; //seg fault when vector is empty
						
//                            //assigned_number =  distance(trig_diff_list.begin(), min_element(trig_diff_list.begin(), trig_diff_list.end())); //position of smallest element in trigger difference vector
//                            //cout << "test" << endl;
//                            time_diff_kpix_ext[kpix].push_back(trig_diff);
//                            if (cycle_num < cycle_checking)
//                            {
//                                AssignedChannelHist[kpix][cycle_num]->Fill(assigned_number);
//                                trigger_difference_per_acq[kpix][cycle_num]->Fill(trig_diff);
//                            }
//                            AssignedTrigger[kpix].push_back(assigned_number);
//                            beam_ext_time_diff->Fill(trig_diff, weight);
//                            trigger_difference[kpix]->Fill(trig_diff, weight);
//                            if (kpix2strip_left.at(channel)!=9999)
//                            {
//                                trigger_diff_connected[kpix]->Fill(trig_diff,weight);
//                            }
//                            else
//                            {
//                                trigger_diff_disconnected[kpix]->Fill(trig_diff,weight);
//                            }
//                            if((trig_diff >= 0.0 )  && (trig_diff  <= 3.0) )
//                            {
//                                if (calibration_check == 1)
//                                {
//                                    hist_timed[kpix][channel][bucket][0]->Fill(double(value)/calib_slope[kpix][channel]/* *pow(10,15) */, weight);
//                                }
//                                else
//                                {
//                                    hist_timed[kpix][channel][bucket][0]->Fill(value, weight);
//                                }
//                                total_timed->Fill(value, weight);
//                                channel_entries_total_timed->Fill(channel, weight);
//                                channel_entries_timed[kpix][bucket]->Fill(channel, weight);
//                                channel_entries_timed[kpix][4]->Fill(channel, weight);
//                                timed_left_strip_entries[kpix][bucket]->Fill(kpix2strip_left.at(channel), weight);
//                                timed_right_strip_entries[kpix][bucket]->Fill(kpix2strip_right.at(channel), weight);
//                                timed_left_strip_entries[kpix][4]->Fill(kpix2strip_left.at(channel), weight);
//                                timed_right_strip_entries[kpix][4]->Fill(kpix2strip_right.at(channel), weight);
								
//                                time_coincident_hits[kpix].push_back(make_pair(kpix2strip_left.at(channel), tstamp));
//                            }
//                            //cout << "DEBUG " << trig_diff << endl;
		
//                        }
//                        //if (kpix != 26 && kpix != 28 && kpix != 30) cout << "Weird..." << kpix << endl;
//                        //}
//                    }
//                }
//            }
//            //for (kpix = 0; kpix < 32; kpix++)
//            //{
//                //if (kpixFound[kpix])
//                //{
//                //for (uint trig = 0; trig < time_ext.size(); trig++)
//                //{
//                    //int mycount = std::count(AssignedTrigger[kpix].begin(), AssignedTrigger[kpix].end(), trig);
//                    //AssignedChannelHist_Total[kpix]->Fill(mycount);
//                //}
//                //}
//            //}
//            //for (unsigned int k = 0; k < Assignment_number.size(); k++)
//            //{
//            //cout << "DEBUG: Assignment_Number " << k << " = " << Assignment_number.at(k) << endl;
//                //if (cycle_num < 1000) AssignedNumberHist[kpix][cycle_num]->Fill(Assignment_number.at(k));
//            //}
//            //ExtTrigPerCycle->Fill(time_ext.size());
//            //for (int kpix = 0; kpix<32; ++kpix)
//            //{
//                //if (kpixFound[kpix])
//                //{
//                    //for (int bucket = 0; bucket<4; ++bucket)
//                    //{
//                        //trig_count[kpix][bucket]->Fill(num_trig_count[kpix][bucket], weight);
//                        //trig_count[kpix][4]->Fill(num_trig_count[kpix][bucket], weight);
						
//                    //}
					
//                //}
//            //}
			
		
		
//			//////////////////////////////////////////
//			// Triggering efficiency and coincidence calculation, takes a lot of time.
//			// -- removed, but can be found in analysisEcal.cxx file
//			//////////////////////////////////////////
			
			
			
			
//			extern_trigger_id=extern_trigger_id+time_ext.size();  // Counting which global external trigger was matched to a channel
			
			
//			//for (int v = 0; v< time_coincident_hits[17].size(); ++v)
//			//{
//				//for (int j = 0; j < time_coincident_hits[19].size(); j++)
//				//{
//					//if (time_coincident_hits[17].at(v).second ==  time_coincident_hits[19].at(j).second)
//					//{
//						//hit_diff_timed->Fill(time_coincident_hits[17].at(v).first - time_coincident_hits[19].at(j).first);
//					//}
					
//				//}
//			//}
			
			
//		}
//		////   Show progress
//		filePos  = dataRead.pos();
//		currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
//		if ( currPct != lastPct )
//		{
//			cout << "\rReading File: " << currPct << " %      " << flush;
//			lastPct = currPct;
//		}
		
//	}
	
//	//for (int kpix = 0; kpix<32; ++kpix)
//		//{
//			//if (kpixFound[kpix])
//			//{
//				//for (int channel = 0; channel < 1024; ++channel)
//				//{
//					//if (chanFound[kpix][channel])
//					//{
//						//for (int bucket = 0; bucket<4; ++bucket)
//						//{
//							//if (bucketFound[kpix][channel][bucket])
//							//{
//								 //hist[kpix][channel][bucket][0]->GetMean();
//								 ////hist[kpix][channel][bucket][0]->Fit("gaus","q");
//							//}
							
//						//}
//					//}
//				//}
				
//			//}
//		//}
	
	
	
//	cout <<  endl << "Full coincidence of sensors with external trigger: " << full_coincidence_channel_entries->GetEntries() << endl;
//	cout << "Three coincidence of sensors: " << three_coincidence << endl;
//	cout << "Two coincidence of sensors: " << two_coincidence << endl;
	
//	cout << endl <<  "An event is currently classified as a monster if the amount of triggers one acquisition clock within a cycle is above " << monster_finder_limit << endl;
//	cout << "_______________________________________________________________" << endl;
	
//	ofstream emptybinfile;
//	emptybinfile.open ("emptybinfile.txt");
	
  
	
//	//for (kpix = 0; kpix < 32; kpix++)
//	//{
//		//if (kpixFound[kpix])
//		//{
//			//cout << "Number of monster events in " << kpix << " = " << monster_counter[kpix] << endl;
//			//cout << "Number of normed monster events in " << kpix << " = " << monster_counter[kpix]*weight << endl;
//			//cout << "Number of entries in KPiX" << kpix << " = " << channel_entries[kpix][4]->GetEntries() << endl << endl;
//			//cout << "Number of NORMED entries in KPiX" << kpix << " = " << (channel_entries[kpix][4]->GetEntries())*weight << endl << endl;
			
//			//for (int bin = 1; bin < 1025; bin++)  // bin 0 is the underflow, therefore need to start counting at 1. Bin == channel_address+1
//			//{
//				//if (channel_entries[kpix][4]->GetBinContent(bin) == 0 && emptybinfile.is_open())
//				//{
//					//emptybinfile << bin-1 << endl;
//				//}
//				//else if (!emptybinfile.is_open()) cout << "error" << endl;
//			//}
//		//}
		
//	//}
//	cout << "Saved the empty bins to file /home/lycoris-dev/Desktop/emptybinfile.txt" << endl ;
//	emptybinfile.close();
	
	
//	for (int k = 0; k < 1024; k++)
//	{
//	cout << "DEBUG channel number: " << k << endl;
//	cout << "DEBUG X coord: " << pixel_kpix[k].x << endl;
//	cout << "DEBUG y coord: " << pixel_kpix[k].y << endl << endl;
//	}
	
//	cout << endl;
//	cout << "Writing root plots to " << outRoot << endl;
//	cout << endl;
	
//	rFile->Write();
//	gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
//	rFile->Close();
	
	
	
//	dataRead.close();
    return(0);
}

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
#include <TTree.h>
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
#include "TBFunctions.h"
#include "kpix_left_and_right.h"
using namespace std;


//////////////////////////////////////////
// Global Variables
//////////////////////////////////////////
vector<TGraphErrors*> calib_graphs; //needed for current loopdir	
//////////////////////////////////////////
// Functions
//////////////////////////////////////////

double median2(vector<double>* v)
{
	if (v  == nullptr )
	{
		//cout << "Found a nullpointer" << endl;
		return 0;
	}
	else
	{
		//cout << "Calculating Median" << endl;
		if (v->empty()) cout << "We have a problem in the Median, it is empty!" << endl;
		size_t n = v->size() / 2;
		cout << n << endl;
		if (n != 0)
		{
			cout << "DEBUG X "  << v->at(n) << endl;
			if (v->size()%2 == 0)
			{
				nth_element(v->begin(), v->begin()+n, v->end());
				cout << "DEBUG Y" << endl;
				nth_element(v->begin(), v->begin()+n-1, v->end());
				return (v->at(n)+v->at(n-1))/2;
			}
			else
			{
				cout << "DEBUG Z "  << v->at(n) << endl;
				nth_element(v->begin(), v->begin()+n, v->end());
				return v->at(n);
			}
		}
		else return 0;
	}
}



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
			
			int found = keyname.find(histname);
			int found_bucket = keyname.find("b0");
			int found2 = keyname.find("DAC");
			if (found != -1 && found2 == -1 && found_bucket != -1) 
			{
				//cout << "Found histogram named " << keyname << endl; 
				TGraphErrors *calib_graph = (TGraphErrors*)key->ReadObj();
				calib_graph->SetName(key->GetName());
				calib_graphs.push_back(calib_graph);
			}
		}
	}
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
	
	const unsigned int n_kpix = 24;
	const unsigned int n_buckets = 4;
	const unsigned int n_channels = 1024;
	
	// cycles to skip in front:
	long int                   skip_cycles_front;
	FILE*                  f_skipped_cycles;
	string                 outtxt;
	
	string                 calState;
	uint                   lastPct;
	uint                   currPct;
	bool                   bucketFound[n_kpix][n_channels][n_buckets] = {false};  // variable that gives true if bucket has an entry (32 possible number of KPiX, n_channels channels per KPiX, 4 buckets per channel)
	bool                   channelFound[n_kpix][n_channels] = {false}; // variable that gives true if a channel has entries
	bool                   kpixFound[n_kpix] = {false}; // variable that gives true if a kpix at the index n (0<=n<32) was found
	uint                   range;
	uint                   value;
	uint                   kpix;
	uint                   channel;
	uint                   bucket;
	double                  tstamp;
	string                 serial;
	KpixSample::SampleType type;
	TTree*					pedestal;
	
	TH1F*					MAD0_v_channel;
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
	
	// Calibration slope, is filled when 
	double					calib_slope[n_kpix][n_channels] = {1}; //ADD buckets later.
	int						calibration_check = 0;

	
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
	tmp << argv[1] << ".tree_pedestal.root";
	outRoot = tmp.str();
	
	
	cout << "Opened data file: " << argv[1] << endl;
	fileSize = dataRead.size();
	filePos  = dataRead.pos();
	
	// Init
	currPct          	= 0;
	lastPct          	= 100;
	
	cout << "\rReading File: 0 %" << flush;  // Printout of progress bar
	//goodTimes       	= 0;
	
	
	// Open root file
	double pedestal_median, pedestal_MAD;
	int kpix_num, channel_num, bucket_num;
	rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
	rFile->cd(); // move into root folder base
	 
	pedestal = new TTree("pedestal_tree", "A ROOT Tree");
	pedestal->Branch("pedestal_median", &pedestal_median, "pedestal_median/D");
	pedestal->Branch("kpix_num", &kpix_num, "kpix_num/I");
	pedestal->Branch("channel_num", &channel_num, "channel_num/I:");
	pedestal->Branch("bucket_num", &bucket_num, "bucket_num/I");
	pedestal->Branch("pedestal_MAD", &pedestal_MAD, "pedestal_MAD/D");
		
	
	MAD0_v_channel = new TH1F("MAD0_v_channel", "MAD0_v_channel; Channel; #Entries", 1024, -0.5, 1023.5);

	range = 0;
	
	
	
	//double weight = 1.0/acqCount; //normalization weight  #entries*weight = #entries/acq.cycle
	double weight = 1.0;//acqProcessed;
	
	//////////////////////////////////////////
	// New histogram generation within subfolder structure
	//////////////////////////////////////////



	
	
	//////////////////////////////////////////
	// Data read for all events for detailed look into single event structure
	//////////////////////////////////////////
	dataRead.open(argv[1]); //open binary file
	int cycle_num = 0;
	int cycle_num_ext = -1;
	dataRead.open(argv[1]); //open file again to start from the beginning

	int extern_trigger_id={0};
	
	
	cycle_num = 0;
	ofstream myfile;
	vector<double>* pedestal_results[n_kpix][n_channels][n_buckets] = {new std::vector<double>};
	
	while ( dataRead.next(&event) )
	{
		cycle_num++;
		//cout << "KPiX event Number: " << event.eventNumber() << endl;
		//cout << "DEBUG" << endl;
		
		
		if ( cycle_num > skip_cycles_front)
		{
		
			uint channel_adc[n_kpix][n_channels];
			//cout << " NEW EVENT " << endl;
			for (uint x=0; x < event.count(); x++)
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
				
				
				if ( type == KpixSample::Data ) // If event is of type KPiX data
				{
					bucketFound[kpix][channel][bucket] = true;
					channelFound[kpix][channel] = true;
					kpixFound[kpix] = true;
					if (kpix == 0)
					cout <<"Weird!" << " k" << kpix << " c" << channel << " b" << bucket << " charge " << value << endl;
					if (calibration_check == 1)
					{
						double charge = value/calib_slope[kpix][channel];
						//cout << "Charge " << charge << endl;
						if (pedestal_results[kpix][channel][bucket] == nullptr)
						{
							pedestal_results[kpix][channel][bucket] = new std::vector<double>;
							if (pedestal_results[kpix][channel][bucket] == nullptr)
							{
								std::cerr << "Memory allocation error for vector kpix " <<
								"KPIX " << kpix << " CHANNEL " << channel << " BUCKET " << bucket << endl;
								exit(-1); // probably best to bail out
							}
							
						}
						
						pedestal_results[kpix][channel][bucket]->push_back(charge);
						//cout << pedestal_results[kpix][channel][bucket]->at(0) << endl;
						
						
					}
						
						
					
				}
				//if ( type == 1 )
				//{
					//cout << kpix << " " << value << endl;
				//}
				//cout << "DEBUG time size" << time_ext.size() << endl;
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
	}
	myfile.close();
	for (int kpix = 0; kpix < n_kpix ; kpix++)
	{
		if (kpixFound[kpix] == true)
		{
			for (int channel = 0; channel < n_channels ; channel++)
			{
				if (channelFound[kpix][channel] == true)
				{
					for (int bucket = 0; bucket < n_buckets ; bucket++)
					{
						if (bucketFound[kpix][channel][bucket] == true)
						{
							//cout << "DEBUG 1 " << endl;
							//if (pedestal_results[kpix][channel][bucket] == nullptr) cout << "Found a nullptr" << endl;
							//else cout << "Found a normal pointer" << endl;
							//cout << "Pedestal for k" << kpix << " c" << channel << " b" << bucket << endl ;
							pedestal_median = median(pedestal_results[kpix][channel][bucket]);
							//cout << "DEBUG 1.01 " << endl;
							
							//cout << pedestal_median << endl;
							pedestal_MAD = MAD(pedestal_results[kpix][channel][bucket]);
							channel_num = channel;
							kpix_num = kpix;
							bucket_num = bucket;
							//cout << "DEBUG 1.1 " << endl;
							//cout << "Median is " << pedestal_median << endl;
							pedestal->Fill();
							//cout << "DEBUG 2 " << endl;
							if (pedestal_MAD == 0)
							{
								MAD0_v_channel->Fill(channel);
							}
							//cout << "DEBUG 3 " << endl;
						}
					}
				}
			}
		}
	}
	for (int kpix = 0; kpix < n_kpix; ++kpix)
	{
		if (kpixFound[kpix] == true)
		{
			for (int channel = 0; channel < n_channels ; channel++)
			{
				if (channelFound[kpix][channel] == true)
				{
					for (int bucket = 0; bucket < n_buckets ; bucket++)
					{
						if (bucketFound[kpix][channel][bucket] == true)
						{
							delete pedestal_results[kpix][channel][bucket];
						}
					}
				}
			}
		}
	}
	
	cout << endl;
	cout << "Writing root plots to " << outRoot << endl;
	cout << endl;
	
	rFile->Write();
	gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
	rFile->Close();
	
	
	
	dataRead.close();
	return(0);
}

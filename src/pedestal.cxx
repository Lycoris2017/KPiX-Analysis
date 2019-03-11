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
	uint                   range;
	uint                   value;
	uint                   kpix;
	uint                   channel;
	uint                   bucket;
	double                  tstamp;
	string                 serial;
	KpixSample::SampleType type;
	TH1F                   	*hist[32][1024][4][2];  // #entries/ADC histograms per channel, bucket, kpix and histogram
	
	
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
	double					calib_slope[32][1024] = {1}; //ADD buckets later.
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
	tmp << argv[1] << ".pedestal.root";
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
	rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
	rFile->cd(); // move into root folder base
	
	
	
	
	while ( dataRead.next(&event) ) // event read to check for filled channels and kpix to reduce number of empty histograms.
	{
		acqCount++;
		if (acqCount > skip_cycles_front)
		{
			acqProcessed++;
			for (uint x=0; x < event.count(); x++)
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
					
				}
				kpixFound[0] = false; // for some reason the system finds a kpix in slot 0
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
	
	
	dataRead.close();
	
	//double weight = 1.0/acqCount; //normalization weight  #entries*weight = #entries/acq.cycle
	double weight = 1.0/acqProcessed;
	
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
			FolderName.str("");
			FolderName << "KPiX_" << kpix;
			rFile->mkdir(FolderName.str().c_str());
			TDirectory *kpix_folder = rFile->GetDirectory(FolderName.str().c_str());
			kpix_folder->cd();
			FolderName.str("");
			FolderName << "Channels";
			kpix_folder->mkdir(FolderName.str().c_str());
			TDirectory *channels_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
			rFile->cd(channels_folder->GetPath());
			for (channel = 0; channel < 1024; channel++)
			{
				if (chanFound[kpix][channel])
				{
					FolderName.str("");
					FolderName << "Channel_" << channel;
					channels_folder->mkdir(FolderName.str().c_str());
					TDirectory *channel_folder = channels_folder->GetDirectory(FolderName.str().c_str());
					rFile->cd(channel_folder->GetPath());

					for (bucket = 0; bucket < 4; bucket++)
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
								tmp << "hist_fc" << "_c" << dec << setw(4) << setfill('0') << channel;
								tmp << "_b" << dec << bucket; // add _b$bucket
								tmp << "_k" << dec << kpix; // add _k$kpix to stringstream
							
								tmp_units.str(""); //set stringstream decribing histogram units to an empty string
								tmp_units << "hist_fc" << "_c" << dec << setw(4) << setfill('0') << channel;
								tmp_units << "_b" << dec << bucket; // add _b$bucket
								tmp_units << "_k" << dec << kpix; // add _k$kpix to stringstream
								tmp_units << "; Charge (fC); #entries/#acq.cycles"; // add title: x label, y label to stringstream
								hist[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(), 16000, -0.5, 3999.5);
								
								
							}
							else
							{
								tmp.str("");  //set stringstream tmp to an empty string
								
								tmp << "hist" << "_c" << dec << setw(4) << setfill('0') << channel;
								tmp << "_b" << dec << bucket; // add _b$bucket
								tmp << "_k" << dec << kpix; // add _k$kpix to stringstream
		
								tmp_units.str(""); //set stringstream decribing histogram units to an empty string
								tmp_units << "hist" << "_c" << dec << setw(4) << setfill('0') << channel;
								tmp_units << "_b" << dec << bucket; // add _b$bucket
								tmp_units << "_k" << dec << kpix; // add _k$kpix to stringstream
								tmp_units << "; Charge (ADC); #entries/#acq.cycles"; // add title: x label, y label to stringstream
								hist[kpix][channel][bucket][0] = new TH1F(tmp.str().c_str(),tmp_units.str().c_str(),8192, -0.5,8191.5);	
							}
						}
					}
				}
			}
		}
	}
	
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
	myfile.open("/scratch/data/dieter.txt");
	
	
	while ( dataRead.next(&event) )
	{
		cycle_num++;
		cout << "KPiX event Number: " << event.eventNumber() << endl;
		if ( cycle_num > skip_cycles_front)
		{
		
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
					if (calibration_check == 1)
					{
						hist[kpix][channel][bucket][0]->Fill(double(value)/calib_slope[kpix][channel]*pow(10,15) , weight);
						
						
					}
					else
					{
						hist[kpix][channel][bucket][0]->Fill(value, weight);
						
					}
					if (bucket == 0 && kpix == 19)
					{
						if (cycle_num == 1)
						{
							myfile << "Itrn = " << event.eventNumber() << " j = " << setw(4) << channel << " k = " << bucket << "  ida=" << kpix << "  x= " << setw(4) << value << " z=" << tstamp << endl; 
						}
					}
					
				}
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
	for (int kpix = 0; kpix<32; ++kpix)
	{
		if (kpixFound[kpix])
		{
			for (int channel = 0; channel < 1024; ++channel)
			{
				if (chanFound[kpix][channel])
				{
					for (int bucket = 0; bucket<4; ++bucket)
					{
						if (bucketFound[kpix][channel][bucket])
						{
							 //cout << "KPIX=" << kpix << "  Channel=" << channel << "  Bucket=" << bucket << "       Mean = " <<  hist[kpix][channel][bucket][0]->GetMean() << endl;
							double mean = hist[kpix][channel][bucket][0]->GetMean();
							double RMS = hist[kpix][channel][bucket][0]->GetRMS();
							
							int firstbin = hist[kpix][channel][bucket][0]->FindFirstBinAbove(0);
							int lastbin = hist[kpix][channel][bucket][0]->FindLastBinAbove(0);
							
							if (calib_slope[kpix][channel]*pow(10,15) > 1 && calib_slope[kpix][channel]*pow(10,15) < 15) hist[kpix][channel][bucket][0]->GetXaxis()->SetRangeUser(firstbin, lastbin);
							double overflow = hist[kpix][channel][bucket][0]->GetBinContent(8001);
							if ( hist[kpix][channel][bucket][0]->GetBinContent(8001) > 0)
							{
								cout << "KPiX Number = " << kpix << "    Channel Number = " << channel << "     Bucket Number = " << bucket << "     overflow content = " << overflow << endl;
							}
							if (mean != 0 && RMS > 0.1)
							{
							 
								hist[kpix][channel][bucket][0]->Fit("gaus","Rq", "", mean-RMS, mean+RMS );
								
								//TF1 *gaussfit = hist[kpix][channel][bucket][0]->GetFunction("gaus");
								//cout << "KPiX Number = " << kpix << "    Channel Number = " << channel << "     Bucket Number = " << bucket << endl;
								//cout << "Chisquare = " << gaussfit->GetChisquare()/gaussfit->GetNDF() << endl;
							}
							
							
							
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

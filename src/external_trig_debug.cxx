//-----------------------------------------------------------------------------
// File          : analysisExternal.cxx
// Author        : Uwe Kraemer (orig. Ryan Herbst) <uwe.kraemer@desy.de>
// Created       : 06/28/2017
// Project       : KPiX Analysis
//-----------------------------------------------------------------------------
// Description :
// External trigger debugging.
//-----------------------------------------------------------------------------
// Copyright (c) 2012 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 05/30/2012: created
// 06/28/2017: large scale rewrite of original calibrationFitter.cpp
// 22/03/2018: clean up ecal plots and add strip plots by <mengqing.wu@desy.de>
// 24/09/2018: Branch off of original analysis.cxx to better focus on external trigger data.
// 06/08/2019: Branch off of cluster_analysis.cxx to focus on debugging of external trigger data.
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
#include "testbeam201905_noise_mask.h"
#include "clustr.h"
#include "PacMan.h"
#include "TBFunctions.h"
using namespace std;


//////////////////////////////////////////
// Global Variables
//////////////////////////////////////////
vector<TGraphErrors*> calib_graphs; //needed for current loopdir	
vector<TH1F*> pedestal_hists;
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
			string keytype = key->GetClassName();
			
			
			if (int(histname.find("pedestal") != -1))
			{
				int found2 = keyname.find("_b0");
				if (int(keytype.find("TH1") != -1))
				{
					if (int( keyname.find("hist_fc_") != -1) && found2 != -1)
					{
						//cout << "Pedestal name " << keyname << endl;
						TH1F *pedestal_hist = (TH1F*)key->ReadObj();
						pedestal_hist->SetName(key->GetName());
						pedestal_hists.push_back(pedestal_hist);
					}
				}
			}
			else if ( int(histname.find("calib") != -1))
			{
				if (int(keytype.find("TH1") == -1) && int(keyname.find("b0") != -1))
				{
					
					if (int(keyname.find("calib") != -1) && int(keyname.find("DAC") == -1)) 
					{
						//cout << "Calibration name " << keyname << endl;
						TGraphErrors *calib_graph = (TGraphErrors*)key->ReadObj();
						calib_graph->SetName(key->GetName());
						calib_graphs.push_back(calib_graph);
					}
				}
			
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
	
	cout << "CURRENTLY USING VERSION FOR NEWDAQ FOR OLD DAQ NEED TO EXCHANGE /Run WITH /20 !" << endl;
	
	
	DataRead               dataRead;  //kpix event classes used for analysis of binary date
	off_t                  fileSize;  //
	off_t                  filePos;   //
	KpixEvent              event;    //
	KpixSample             *sample;   //
	
	const unsigned int n_buckets = 1;
	const unsigned int n_kpix = 2;
	const unsigned int n_blocks = 32;
	const unsigned int n_channels = 1024;
	const unsigned int n_BCC = 8192;
	const unsigned int n_strips = 1840;
	
	const unsigned int n_input_kpix = 24;
	const unsigned int n_input_buckets = 4;
	const unsigned int n_input_channels = 1024;
	
	// cycles to skip in front:
	long int                   skip_cycles_front;
	FILE*                  f_skipped_cycles;
	string                 outtxt;
	
	string                 calState;
	uint                   lastPct;
	uint                   currPct;
	bool                   kpixFound[n_kpix] = {false}; // variable that gives true if a kpix at the index n (0<=n<32) was found
	bool                   channelFound[n_kpix][n_channels] = {false};
	bool                   bucketFound[n_kpix][n_channels][n_buckets] = {false};
	uint                   x;
	uint                   value;
	uint                   kpix;
	uint                   sensor;
	uint                   channel;
	uint                   bucket;
	double                  tstamp;
	uint 					subCount;
	double 					bunchClk;
	double 					y;
	string                 serial;
	KpixSample::SampleType type;	

	uint64_t				frameruntime;
	uint64_t				runtime;

	//TTree*					cluster_tree;
	
	TH1F					*fc_response_medCM_subtracted[n_kpix];
		
	TH1F					*noise_distribution[n_kpix];
	TH1F					*noise_distribution_sensor;
	TH1F					*noise_v_position;
	TH1F					*noise_v_channel[n_kpix];
	TGraphErrors 			*noise_v_time_graph;
	TH1F					*noise_v_time[n_kpix];
	TH1F					*noise_v_time_block[n_kpix][n_blocks];
	
	// Stringstream initialization for histogram naming
	stringstream           tmp;
	stringstream           tmp_units;
	
	// Stringstream initialization for folder naming
	
	stringstream			FolderName;
	
	ofstream				noise_file;
	ofstream               xml;
	ofstream               csv;
	unsigned int           acqCount; // acquisitionCount
	unsigned int           acqProcessed;
	string                 outRoot;
	TFile                  *rFile;
	stringstream           crossString;
	stringstream           crossStringCsv;
	XmlVariables           config;
	ofstream               debug;
	
	// Calibration slope, is filled when 
	double					calib_slope[n_input_kpix][n_input_channels] = {1}; //ADD buckets later.
	
	double					pedestal_MedMAD[n_input_kpix][n_input_channels][n_input_buckets][2] = {0};
	int						calibration_check = 0;
	double 					noise[n_kpix][n_channels];
	double 					MaximumSoN[n_strips] = {0};
	//int						pedestal_check = 0;
	
	unordered_map<uint, uint> kpix2strip_left;
	unordered_map<uint, uint> kpix2strip_right;
	
	
	kpix2strip_left = kpix_left();
	kpix2strip_right = kpix_right();
	
	pixel					pixel_kpix[n_channels];
	pixel_mapping(pixel_kpix);
	
	unordered_map<uint, uint> sensor2layer;
	
	sensor2layer.insert(make_pair(0, 15));
	sensor2layer.insert(make_pair(1, 14));
	sensor2layer.insert(make_pair(2, 13));
	sensor2layer.insert(make_pair(3, 10));
	sensor2layer.insert(make_pair(4, 11));
	sensor2layer.insert(make_pair(5, 12));
	sensor2layer.insert(make_pair(6, 9999));
	sensor2layer.insert(make_pair(7, 9999));
	sensor2layer.insert(make_pair(8, 9999));
	sensor2layer.insert(make_pair(9, 9999));
	sensor2layer.insert(make_pair(10, 9999));
	sensor2layer.insert(make_pair(11, 9999));
	
	
	//////////////////////////////////////////
	// Skip cycle/Read calibration file
	//////////////////////////////////////////
	
	// Data file is the first and only arg
	if (argc < 4) {
	cout << "Usage: ./analysis data_file [skip_cycles_front (int)]||[read calibration input file (char)] \n";
	return(1);
	}
	char* end;
	
	
	if ( argc >= 4 ) {
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
			loopdir(calibration_file, "calib");
			
			cout << "Number of calibration slopes = " << calib_graphs.size() << endl;
			
			for (unsigned int i = 0; i<calib_graphs.size(); ++i)
			{
				//cout << "Current key1 = " << cal_key->GetClassName() << endl;
				
				string calib_name         = calib_graphs[i]->GetName();
				//cout << calib_name << endl;
				size_t kpix_num_start     = calib_name.find("_k")+2;
				size_t channel_num_start  = calib_name.find("_c")+2;
				size_t kpix_num_length       = calib_name.length() - kpix_num_start;
				size_t channel_num_length    = calib_name.find("_b") - channel_num_start;
				
				//cout << kpix_num_start << endl;
				//cout << kpix_num_length << endl;
				//
				//cout << channel_num_start << endl;
				//cout << channel_num_length << endl;
				
			    string channel_string = calib_name.substr(channel_num_start, channel_num_length);
			    string kpix_string = calib_name.substr(kpix_num_start, kpix_num_length);
			    
			   //cout << "Channel Number = " <<  channel_string << endl;
			   //cout << "KPiX Number = " << kpix_string << endl;
			    
			    
			    int kpix_num = stoi(kpix_string)/6-1;
			    int channel_num = stoi(channel_string);
				
				//cout << "KPiX Number = " << kpix_num << endl;
				//cout << "Channel Number = " << channel_num << endl;
				
				calib_slope[kpix_num][channel_num] = calib_graphs[i]->GetFunction("pol1")->GetParameter(1);
				//cout << "Slope of KPiX " << kpix_num << " and channel " << channel_num << " is " <<  calib_slope[kpix_num][channel_num] << endl;
				
			}
			
			cout << " -- Reading " << argv[3] << " as pedestal subtraction input file." << endl;
			skip_cycles_front = 0;
			TFile *pedestal_file = TFile::Open(argv[3]);

			TTree* pedestal_tree = (TTree*)pedestal_file->Get("pedestal_tree");
			
			
			int kpix_num, channel_num, bucket_num;
			double pedestal_median, pedestal_MAD;
			
			pedestal_tree->SetBranchAddress("pedestal_median", &pedestal_median);
			pedestal_tree->SetBranchAddress("kpix_num", &kpix_num);
			pedestal_tree->SetBranchAddress("channel_num", &channel_num);
			pedestal_tree->SetBranchAddress("bucket_num", &bucket_num);
			pedestal_tree->SetBranchAddress("pedestal_MAD", &pedestal_MAD);
			
			long int nentries = pedestal_tree->GetEntries();
			for (long int i = 0; i < nentries; ++i)
			{
				pedestal_tree->GetEntry(i);
				if (kpix_num == 6 || kpix_num == 12)
				{
					pedestal_MedMAD[kpix_num/6-1][channel_num][bucket_num][0] = pedestal_median;
					pedestal_MedMAD[kpix_num/6-1][channel_num][bucket_num][1] = pedestal_MAD;
				}
			}
			
			
			
			gROOT->GetListOfFiles()->Remove(calibration_file);
			calibration_file->Close();
			gROOT->GetListOfFiles()->Remove(pedestal_file);
			pedestal_file->Close();
			
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
	string pedestalname = argv[3];
	string outname = argv[1];
	
	
	size_t name_start  = pedestalname.find("/Run") + 1;
	size_t name_length = pedestalname.find(".dat") - name_start;
	
	pedestalname = pedestalname.substr(name_start, name_length);
	
	cout << "Name of output file is " <<  pedestalname << endl;
	tmp.str("");
	tmp << argv[1] << "_" << pedestalname << ".debug.root";
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
	vector<double>* vec_corr_charge[n_kpix]; //delete does not work if I do not initialize all vectors
	std::map<int, double> common_modes_median[n_kpix];
	while ( dataRead.next(&event) ) // event read to check for filled channels and kpix to reduce number of empty histograms.
	{
		
		acqCount++;		
		if (acqCount > skip_cycles_front)
		{
			acqProcessed++;
			
			
			for (x=0; x < event.count(); x++)
			{
		
				//// Get sample
				sample  = event.sample(x);
				kpix    = (sample->getKpixAddress())/6-1;
				tstamp  = sample->getSampleTime();
				channel = sample->getKpixChannel();
				bucket  = sample->getKpixBucket();
				value   = sample->getSampleValue();
				type    = sample->getSampleType();
				
				//cout << type <<endl;
				//cout << "DEBUG 2" << endl;
				if ( type == KpixSample::Data )
				{
					//cout << kpix << endl;
					kpixFound[kpix]          = true;
					channelFound[kpix][channel] = true;
					bucketFound[kpix][channel][bucket] = true;
					//cout << "Found KPIX " << kpix << endl;
					if (bucket == 0)
					{
						if (calibration_check == 1)
						{
							//cout << "Pedestal MAD " << pedestal_MedMAD[kpix][channel][bucket][1] << endl;
							if (pedestal_MedMAD[kpix][channel][bucket][1] != 0 && calib_slope[kpix][channel] != 0)
							{
								
								if (vec_corr_charge[kpix] == NULL)
								{
									vec_corr_charge[kpix] = new std::vector<double>;
									if (vec_corr_charge[kpix]==NULL)
									{
										std::cerr << "Memory allocation error for vector kpix " <<
										kpix <<std::endl;
										exit(-1); // probably best to bail out
									}
								} 
								
								//cout << "DEBUG 2.1 " << kpix << endl;
								double charge_value = double(value)/calib_slope[kpix][channel];
								double corrected_charge_value_median = charge_value - pedestal_MedMAD[kpix][channel][bucket][0];
								vec_corr_charge[kpix]->push_back(corrected_charge_value_median);
							}
						}
					}
				}
			}
			for (unsigned int k = 0; k < n_kpix ; k++)
			{
				if (vec_corr_charge[k] != nullptr) 
				{
					//cout << "Debug size of vec: " << vec_corr_charge[k]->size() << endl; 
					common_modes_median[k].insert(std::pair<int, double>(event.eventNumber(), median(vec_corr_charge[k])));
					delete vec_corr_charge[k];
					vec_corr_charge[k] = NULL;
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
	//cout << tstamp << endl;
	dataRead.close();
	double weight = 1.0/acqProcessed;
	
	
	
	//////////////////////////////////////////
	// New histogram generation within subfolder structure
	//////////////////////////////////////////
	
	const unsigned int response_bins = 220;
	const double response_xmin = -20.5;
	const double response_xmax = 19.5;
	
	
	//clustr* clusttree = new clustr;
	//int sensr;
	//cluster_tree = new TTree("cluster_tree", "A ROOT Tree");
	//cluster_tree->Branch("cluster", &clusttree);
	//cluster_tree->Branch("sensor", &sensr, "sensr/I");

	//TH1F* mean_noise = new TH1F("mean_noise_left", "mean_noise; noise(fC); entries", 100, -0.05, 0.95);
	
	for (sensor = 0; sensor < 1; sensor++) //looping through all possible kpix
	{
		if (kpixFound[(sensor*2)] || kpixFound[(sensor*2+1)])
		{
			rFile->cd(); //producing subfolder for kpix same as above for the event subfolder structure
			FolderName.str("");
			FolderName << "Sensor_" << sensor;
			rFile->mkdir(FolderName.str().c_str());
			TDirectory *sensor_folder = rFile->GetDirectory(FolderName.str().c_str());
			sensor_folder->cd();
			
			
			
			tmp.str("");
			tmp << "noise_v_position_s" << sensor << "_b0";
			noise_v_position  = new TH1F(tmp.str().c_str(), "Noise; #mum; Noise (fC)", 1840,-46000, 46000);
			
			tmp.str("");
			tmp << "noise_distribution_s" << sensor << "_b0";
			noise_distribution_sensor = new TH1F(tmp.str().c_str(), "noise_distribution; Noise(fC);   #channels", 100,-0.005, 0.995);
			
			
			
			
			
			for (int k = 0; k < 2; k++) //looping through all possible kpix (left and right of each sensor)
			{
				kpix = (sensor*2)+k;
				if (kpixFound[kpix])
				{
					FolderName.str("");
					FolderName << "KPiX_" << kpix;
					sensor_folder->mkdir(FolderName.str().c_str());
					TDirectory *kpix_folder = sensor_folder->GetDirectory(FolderName.str().c_str());
					rFile->cd(kpix_folder->GetPath());
					
					tmp.str("");
					tmp << "fc_response_median_made_CMmedian_subtracted_k" << kpix << "_b0";
					fc_response_medCM_subtracted[kpix] = new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", response_bins, response_xmin, response_xmax);
					
					tmp.str("");
					tmp << "noise_distribution_k" << kpix << "_b0";
					noise_distribution[kpix] = new TH1F(tmp.str().c_str(), "noise_distribution; Noise(fC);   #channels", 100,-0.005, 4.995);
					cout << "Found the following KPiX in the file: " <<  kpix << endl;
					
					tmp.str("");
					tmp << "noise_v_channel_k" << kpix << "_b0";
					noise_v_channel[kpix]  = new TH1F(tmp.str().c_str(), "Noise; channel; Noise (fC)", 1024,0, 1023);
					
					tmp.str("");
					tmp << "noise_v_time_k" << kpix << "_b0";
					tmp_units.str("");
					tmp_units << "noise_v_time; Time (BCC); Noise (fC) ";
					noise_v_time[kpix]  = new TH1F(tmp.str().c_str(), tmp_units.str().c_str(), 8192, 0, 8191);
					
					
					for (int b = 0; b < 32; ++b)
					{
						tmp.str("");
						tmp << "noise_v_time_block" << b << "_k" << kpix << "_b0";
						tmp_units.str("");
						tmp_units << "noise_v_time_block; Time (BCC); Noise (fC) ";
						noise_v_time_block[kpix][b]  = new TH1F(tmp.str().c_str(), tmp_units.str().c_str(), 8192, 0, 8191);
					}
					
			
				}
			}
		}
		
	}
	//////////////////////////////////////////
	// Data read for all events for detailed look
	//////////////////////////////////////////
	dataRead.open(argv[1]); //open file again to start from the beginning


	int header = 1;
	
	std::vector<double>* corrected_charge_vec[n_kpix][n_channels] = {nullptr};
	std::vector<double>* corrected_charge_vec_time[n_kpix][n_BCC]  = {nullptr};
	std::vector<double>* corrected_charge_vec_time_block[n_kpix][n_BCC][n_blocks]  = {nullptr};
	//std::vector<double>* test[n_kpix][n_channels] = {nullptr}; // = { new std::vector<double> };
	
	while ( dataRead.next(&event) ) //preread to determine noise value of each channel in each KPiX.
	{
		int not_empty= 0;
		for (x=0; x < event.count(); x++)
		{
			//cout << "DEBUG: EVENT COUNT " << event.count() << endl;
			//// Get sample
			sample  = event.sample(x);
			kpix    = (sample->getKpixAddress())/6-1;
			channel = sample->getKpixChannel();
			bucket  = sample->getKpixBucket();
			value   = sample->getSampleValue();
			type    = sample->getSampleType();
			tstamp  = sample->getSampleTime();
			sensor = kpix/2;
			unsigned int block = channel/32;
			bunchClk = sample->getBunchCount();
			subCount = sample->getSubCount();
		
			if ( type == KpixSample::Data ) // If event is of type KPiX data
			{
				if (bucket == 0)
				{
					if (pedestal_MedMAD[kpix][channel][bucket][1] != 0 && calib_slope[kpix][channel] != 0 ) //ensuring we ignore 0 MAD channels
					{
						//cout << "DEBUG 1" << endl;
						double charge_value = double(value)/calib_slope[kpix][channel];
						double corrected_charge_value_median = charge_value - pedestal_MedMAD[kpix][channel][bucket][0];
						double charge_CM_corrected = corrected_charge_value_median - common_modes_median[kpix].at(event.eventNumber());
						if (isnan(charge_CM_corrected))
						{
							cout << "Slopes " << calib_slope[kpix][channel] << endl;
							cout << "Pedestal " << pedestal_MedMAD[kpix][channel][bucket][0] << endl;
							cout << "CommonMOde " << common_modes_median[kpix].at(event.eventNumber()) << endl;
						}
						if (corrected_charge_vec_time_block[kpix][int(tstamp)][block] == NULL)
						{
							corrected_charge_vec_time_block[kpix][int(tstamp)][block] = new std::vector<double>;
							if (corrected_charge_vec_time_block[kpix][int(tstamp)][block]==NULL)
							{
								std::cerr << "Memory allocation error for vector kpix " <<
								kpix << " tstamp "  << tstamp << " block "<< block <<" " <<std::endl;
								exit(-1); // probably best to bail out
							}
						}
						if (corrected_charge_vec_time[kpix][int(tstamp)] == NULL)
						{
							corrected_charge_vec_time[kpix][int(tstamp)] = new std::vector<double>;
							if (corrected_charge_vec_time[kpix][int(tstamp)] == NULL)
							{
								std::cerr << "Memory allocation error for vector kpix " <<
								kpix << " tstamp "  << tstamp <<" " <<std::endl;
								exit(-1); // probably best to bail out
							}
						} 
						if (corrected_charge_vec[kpix][channel] == NULL)
						{
							corrected_charge_vec[kpix][channel] = new std::vector<double>;
							if (corrected_charge_vec[kpix][channel] == NULL)
							{
								std::cerr << "Memory allocation error for vector kpix " <<
								kpix << " tstamp "  << tstamp <<" " <<std::endl;
								exit(-1); // probably best to bail out
							}
						} 
						corrected_charge_vec[kpix][channel]->push_back(charge_CM_corrected);
						corrected_charge_vec_time[kpix][int(tstamp)]->push_back(charge_CM_corrected);
						corrected_charge_vec_time_block[kpix][int(tstamp)][block]->push_back(charge_CM_corrected);
						
					}
					//else if (pedestal_MedMAD[kpix][channel][bucket][1] == 0 && kpix == 0) cout << "1KPIX " << kpix << " Channel " << channel << endl;
				}
			}
		}
		filePos  = dataRead.pos();
		currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
		if ( currPct != lastPct ) 
		{
			cout << "\rReading File for noise determination: " << currPct << " %      " << flush;
			lastPct = currPct;
		}
	}
	dataRead.close();
	
	 // END OF PREREAD
	 // BEGIN OF NOISE CALCULATION
	for (int k = 0; k < n_kpix; ++k)
	{
		if (kpixFound[k])
		{
			double x[n_BCC];
			double x_err[n_BCC];
			double y_graph[n_BCC];
			double y_graph_err[n_BCC];
			int sensor = k/2;
			for (int c = 0; c < n_channels; ++c)
			{
				if (channelFound[k][c])
				{
					int sensor = k/2;
					int strip = 9999;
					if (k%2 == 0) strip = kpix2strip_left.at(c);// if left kpix
					else strip  = kpix2strip_right.at(c); // if right kpix
					
					if (pedestal_MedMAD[k][c][0][1] != 0)
					{
						y = yParameterSensor(strip, sensor);
						noise[k][c] = 1.4826*MAD(corrected_charge_vec[k][c]);
//						if (isnan(noise[k][c]))
//						{
//							cout << k << " " << c << endl;
//							cout << calib_slope[k][c] << endl;
//							for (auto const i:(*corrected_charge_vec[k][c]))
//							{
//								cout << "Bla " << i << endl;
//							}
//						}
						if (strip != 9999)
						{
							noise_distribution[k]->Fill(noise[k][c]);
							noise_distribution_sensor->Fill(noise[k][c]);
							noise_v_position->Fill(y, noise[k][c]);
						}
						noise_v_channel[k]->SetBinContent(c+1, noise[k][c]);
						
					}
				}
				
			}
			for (int t = 0; t < n_BCC ; ++t)
			{
				for (int b = 0; b < n_blocks;  ++b)
				{
					noise_v_time_block[k][b]->SetBinContent(t+1, 1.4826*MAD(corrected_charge_vec_time_block[k][t][b]));
				}
				noise_v_time[k]->SetBinContent(t+1, 1.4826*MAD(corrected_charge_vec_time[k][t]));
				x[t] = t;
				x_err[t] = 0;
				y_graph[t] = 1.4826*MAD(corrected_charge_vec_time[k][t]);
				y_graph_err[t] = 0;
				//cout << x[t] << endl;
			}	
			noise_v_time_graph = new TGraphErrors(2*n_BCC ,x, y_graph, x_err, y_graph_err);
			noise_v_time_graph->Draw("Ap");
			noise_v_time_graph->GetXaxis()->SetTitle("Time (BCC)");
			noise_v_time_graph->GetYaxis()->SetTitle("Noise (fC)");
			tmp.str("");
			tmp << "noise_v_time_graph_k" << k;
			
			noise_v_time_graph->SetTitle(tmp.str().c_str());
			
			noise_v_time_graph->Write(tmp.str().c_str());
			
			
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

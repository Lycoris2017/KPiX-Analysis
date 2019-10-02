//-----------------------------------------------------------------------------
// File          : analysisExternal.cxx
// Author        : Uwe Kraemer (orig. Ryan Herbst) <uwe.kraemer@desy.de>
// Created       : 06/28/2017
// Project       : KPiX Analysis
//-----------------------------------------------------------------------------
// Description :
// Analysis of external triggering KPiX Data.
//-----------------------------------------------------------------------------
// Copyright (c) 2012 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 05/30/2012: created
// 06/28/2017: large scale rewrite of original calibrationFitter.cpp
// 22/03/2018: clean up ecal plots and add strip plots by <mengqing.wu@desy.de>
// 24/09/2018: Branch off of original analysis.cxx to better focus on external trigger data.
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
#include "testbeam201907_noise_mask.h"
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
	cout <<"DEBUG1" << endl;
	//////////////////////////////////////////
	// Class declaration and histogram initialization
	//////////////////////////////////////////
	
	cout << "CURRENTLY USING VERSION FOR NEWDAQ FOR OLD DAQ NEED TO EXCHANGE /Run WITH /20 !" << endl;
	
	
	DataRead               dataRead;  //kpix event classes used for analysis of binary date
	off_t                  fileSize;  //
	off_t                  filePos;   //
	KpixEvent              event;    //
	KpixSample             *sample;   //
	
	const unsigned int n_buckets = 4;
	const unsigned int n_kpix = 32;//24;
//	const unsigned int n_blocks = 32;
	const unsigned int n_channels = 1024;
	const unsigned int n_BCC = 8192;
//	const unsigned int n_strips = 1840;
	
	// cycles to skip in front:
	long int                   skip_cycles_front;
	FILE*                  f_skipped_cycles;
	string                 outtxt;
	
	string                 calState;
	uint                   lastPct;
	uint                   currPct;
	bool                   kpixFound[n_kpix] = {false}; // variable that gives true if a kpix at the index n (0<=n<32) was found
	bool                   channelFound[n_kpix][n_channels] = {false};
//	bool                   bucketFound[n_kpix][n_channels][n_buckets] = {false};
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
		
	TH1F 					*Max_SoN[n_kpix/2];
		
		
	TH1F 					*cluster_position_y[n_kpix/2][3];
	//TH1F 					*cluster_position_x[n_kpix/2][3];
	TH1F 					*clusters[n_kpix/2][2]; // only 2 because it makes no sense to ask for the number of clusters in the seed case cause it will always be 1.
	TH1F 					*cluster_charge[n_kpix][3];
	TH1F 					*cluster_significance[n_kpix][3];
	TH1F 					*cluster_significance2[n_kpix][3];
	TH1F 					*cluster_size[n_kpix][3];
	TH1F 					*cluster_sigma[n_kpix][3];
	
		
		
	TH2F					*cluster_correlation[3][n_kpix][n_kpix-1];
	TH1F					*cluster_offset_y[3][n_kpix/2][n_kpix/2-1];
	TH1F					*cluster_offset_x[3][n_kpix/2][n_kpix/2-1];

	TH1F					*cuts_entries[3][n_kpix/2][n_kpix/2-1];
	TH1F					*cuts_purity[3][n_kpix/2][n_kpix/2-1];
	TH1F					*cuts_performance[3][n_kpix/2][n_kpix/2-1];

	TH1F					*noise_distribution[n_kpix];
	TH1F					*noise_distribution_sensor[n_kpix/2];
	TH1F					*noise_v_position[n_kpix/2];
	TH1F					*noise_v_channel[n_kpix];
	
	TH1F					*noise_v_time[n_kpix/2];
	
	// Stringstream initialization for histogram naming
	stringstream           tmp;
	stringstream           tmp_units;
	
	// Stringstream initialization for folder naming
	
	stringstream			FolderName;
	
	ofstream				claus_file;
	ofstream				noise_file;
	ofstream               xml;
	ofstream               csv;
	uint                   acqCount = 0; // acquisitionCount
	uint                   acqProcessed;
	string                 outRoot;
	TFile                  *rFile;
	stringstream           crossString;
	stringstream           crossStringCsv;
	XmlVariables           config;
	ofstream               debug;
	
	// Calibration slope, is filled when 
	double					calib_slope[n_kpix][n_channels] = {1}; //ADD buckets later.
	
	double					pedestal_MedMAD[n_kpix][n_channels][n_buckets][2] = {0};
	int						calibration_check = 0;
	double 					noise[n_kpix][n_channels] = {0};
	double 					MaximumSoN[n_kpix/2][1840] = {0};
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
	//cout << "DEBUG: 1" << endl;
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
			    
			    
			    int kpix_num = stoi(kpix_string);
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
				pedestal_MedMAD[kpix_num][channel_num][bucket_num][0] = pedestal_median;
				pedestal_MedMAD[kpix_num][channel_num][bucket_num][1] = pedestal_MAD;
				//if (pedestal_MedMAD[kpix_num][channel_num][bucket_num][1] != 0)
				//{
					//cout << "Pedestal MAD " << pedestal_MedMAD[kpix_num][channel_num][bucket_num][1] << " kpix " << kpix_num << " channel " << channel_num << " bucket " << bucket_num << endl;
				//}
			}
			
			
			
			gROOT->GetListOfFiles()->Remove(calibration_file);
			calibration_file->Close();
			gROOT->GetListOfFiles()->Remove(pedestal_file);
			pedestal_file->Close();
			
		}
		
	
	}
	
	
	//cout << "DEBUG: 2" << endl;
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
	tmp << argv[1] << "_" << pedestalname << ".cluster.root";
	outRoot = tmp.str();

	tmp.str("");
	tmp << argv[1] << ".GBL_input.txt" ;
	cout << "Name of GBL input file is - " << tmp.str();
	//claus_file.open("claus_file_new.txt");
	claus_file.open(tmp.str());
	
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
	vector<double>* vec_corr_charge[n_kpix] = {nullptr}; //delete does not work if I do not initialize all vectors
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
				kpix    = sample->getKpixAddress();
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
//					bucketFound[kpix][channel][bucket] = true;
					//cout << "Found KPIX " << kpix << endl;
					if (bucket == 0)
					{
						if (calibration_check == 1)
						{
							//cout << "2nd Pedestal MAD " << pedestal_MedMAD[kpix][channel][bucket][1] << " kpix " << kpix << " channel " << channel << " bucket " << bucket << endl;
							if (pedestal_MedMAD[kpix][channel][bucket][1] != 0 && calib_slope[kpix][channel] != 0)
							{
								
								if (vec_corr_charge[kpix] == nullptr)
								{
									vec_corr_charge[kpix] = new std::vector<double>;
									if (vec_corr_charge[kpix]==nullptr)
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
								//cout << "Pointer of vec_corr_charge " << vec_corr_charge[kpix] << endl;
							}
						}
					}
				}
			}
			for (unsigned int k = 0; k < n_kpix ; k++)
			{
				//cout << "Pointer of vec_corr_charge " << vec_corr_charge[k] << endl;
				if (vec_corr_charge[k] != nullptr) 
				{
					//cout << "Debug size of vec: " << vec_corr_charge[k]->size() << endl; 
					common_modes_median[k].insert(std::pair<int, double>(event.eventNumber(), median(vec_corr_charge[k])));
					delete vec_corr_charge[k];
					//cout << "Common modes median of EventNumber " << event.eventNumber() << " entry " << common_modes_median[k].at(event.eventNumber()) << endl;
					vec_corr_charge[k] = nullptr;
					
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
	//cout << "DEBUG: 3" << endl;
	//cout << tstamp << endl;
	dataRead.close();
	double weight = 1.0/acqProcessed;
	;//acqProcessed;
	
	
	
	//////////////////////////////////////////
	// New histogram generation within subfolder structure
	//////////////////////////////////////////
	
	int response_bins = 220;
	double response_xmin = -20.5;
	double response_xmax = 19.5;
	
	
	//clustr* clusttree = new clustr;
	//int sensr;
	//cluster_tree = new TTree("cluster_tree", "A ROOT Tree");
	//cluster_tree->Branch("cluster", &clusttree);
	//cluster_tree->Branch("sensor", &sensr, "sensr/I");

	//TH1F* mean_noise = new TH1F("mean_noise_left", "mean_noise; noise(fC); entries", 100, -0.05, 0.95);
	
	for (sensor = 0; sensor < n_kpix/2; sensor++) //looping through all possible kpix
	{
		if (kpixFound[(sensor*2)] || kpixFound[(sensor*2+1)])
		{
			rFile->cd(); //producing subfolder for kpix same as above for the event subfolder structure
			FolderName.str("");
			FolderName << "Sensor_" << sensor;
			rFile->mkdir(FolderName.str().c_str());
			TDirectory *sensor_folder = rFile->GetDirectory(FolderName.str().c_str());
			sensor_folder->cd();
			
			for (int s = sensor+1; s < n_kpix/2; ++s)
			{
				if (kpixFound[(s*2)] || kpixFound[(s*2+1)])
				{
					tmp.str("");
					tmp << "cluster_offset_seed_y_sens" << sensor << "_to_sens" << s << "_b0";
					cluster_offset_y[0][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
					tmp.str("");
					tmp << "cluster_offset_seed_x_sens" << sensor << "_to_sens" << s << "_b0";
					cluster_offset_x[0][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
					tmp.str("");
					tmp << "cluster_correlation_seed_sens0" << sensor << "_to_sens" << s << "_b0";
					tmp_units.str("");
					tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
					cluster_correlation[0][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

					tmp.str("");
					tmp << "cluster_offset_all_y_sens" << sensor << "_to_sens" << s << "_b0";
					cluster_offset_y[1][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
					tmp.str("");
					tmp << "cluster_offset_all_x_sens" << sensor << "_to_sens" << s << "_b0";
					cluster_offset_x[1][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
					tmp.str("");
					tmp << "cluster_correlation_all_sens" << sensor << "_to_sens" << s << "_b0";
					tmp_units.str("");
					tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
					cluster_correlation[1][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

					tmp.str("");
					tmp << "cluster_offset_cuts_y_sens" << sensor << "_to_sens" << s << "_b0";
					cluster_offset_y[2][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
					tmp.str("");
					tmp << "cluster_offset_cuts_x_sens" << sensor << "_to_sens" << s << "_b0";
					cluster_offset_x[2][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
					tmp.str("");
					tmp << "cluster_correlation_cuts_sens" << sensor << "_to_sens" << s << "_b0";
					tmp_units.str("");
					tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
					cluster_correlation[2][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

//					tmp.str("");
//					tmp << "cluster_offset_HighSigmaCluster_y_sens" << sensor << "_to_sens" << s << "_b0";
//					cluster_offset_y[3][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
//					tmp.str("");
//					tmp << "cluster_offset_HighSigmaCluster_x_sens" << sensor << "_to_sens" << s << "_b0";
//					cluster_offset_x[3][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
//					tmp.str("");
//					tmp << "cluster_correlation_HighSigmaCluster_sens" << sensor << "_to_sens" << s << "_b0";
//					tmp_units.str("");
//					tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
//					cluster_correlation[3][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

//					tmp.str("");
//					tmp << "cluster_offset_4sigma_2fC_3size_sens_y_sens" << sensor << "_to_sens" << s << "_b0";
//					cluster_offset_y[4][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
//					tmp.str("");
//					tmp << "cluster_offset_4sigma_2fC_3size_sens_x_sens" << sensor << "_to_sens" << s << "_b0";
//					cluster_offset_x[4][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
//					tmp.str("");
//					tmp << "cluster_correlation_4sigma_2fC_3size_sens_sens" << sensor << "_to_sens" << s << "_b0";
//					tmp_units.str("");
//					tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
//					cluster_correlation[4][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

//					tmp.str("");
//					tmp << "cluster_offset_4sigma_3fC_3size_sens_y_sens" << sensor << "_to_sens" << s << "_b0";
//					cluster_offset_y[5][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
//					tmp.str("");
//					tmp << "cluster_offset_4sigma_3fC_3size_sens_x_sens" << sensor << "_to_sens" << s << "_b0";
//					cluster_offset_x[5][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
//					tmp.str("");
//					tmp << "cluster_correlation_4sigma_3fC_3size_sens_sens" << sensor << "_to_sens" << s << "_b0";
//					tmp_units.str("");
//					tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
//					cluster_correlation[5][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

//					tmp.str("");
//					tmp << "cluster_offset_4sigma_2_5fC_3size_sens_y_sens" << sensor << "_to_sens" << s << "_b0";
//					cluster_offset_y[6][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
//					tmp.str("");
//					tmp << "cluster_offset_4sigma_2_5fC_3size_sens_x_sens" << sensor << "_to_sens" << s << "_b0";
//					cluster_offset_x[6][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
//					tmp.str("");
//					tmp << "cluster_correlation_4sigma_2_5fC_3size_sens_sens" << sensor << "_to_sens" << s << "_b0";
//					tmp_units.str("");
//					tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
//					cluster_correlation[6][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

//					tmp.str("");
//					tmp << "cluster_offset_test_sens_y_sens" << sensor << "_to_sens" << s << "_b0";
//					cluster_offset_y[7][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
//					tmp.str("");
//					tmp << "cluster_offset_test_sens_x_sens" << sensor << "_to_sens" << s << "_b0";
//					cluster_offset_x[7][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
//					tmp.str("");
//					tmp << "cluster_correlation_test_sens_sens" << sensor << "_to_sens" << s << "_b0";
//					tmp_units.str("");
//					tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
//					cluster_correlation[7][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

//					tmp.str("");
//					tmp << "sigma_cut_purity_s" << sensor << "_s" << s << "_b0";
//					cuts_purity[0][sensor][s] = new TH1F(tmp.str().c_str(), "sigma_cut_puritye; cuts (#sigma); purity", 70, 0, 35);
//					tmp.str("");
//					tmp << "charge_cut_purity_s" << sensor << "_s" << s << "_b0";
//					cuts_purity[1][sensor][s] = new TH1F(tmp.str().c_str(), "charge_cut_purity; cuts (fC); purity", 70, 0, 35);
//					tmp.str("");
//					tmp << "size_cut_purity_s" << sensor << "_s" << s << "_b0";
//					cuts_purity[2][sensor][s] = new TH1F(tmp.str().c_str(), "size_cut_purity; cuts; purity", 35, 0, 35);

//					tmp.str("");
//					tmp << "sigma_cut_entries_s" << sensor << "_s" << s << "_b0";
//					cuts_entries[0][sensor][s] = new TH1F(tmp.str().c_str(), "sigma_cut_entries; cuts (#sigma); entries", 70, 0, 35);
//					tmp.str("");
//					tmp << "charge_cut_entries_s" << sensor << "_s" << s << "_b0";
//					cuts_entries[1][sensor][s] = new TH1F(tmp.str().c_str(), "charge_cut_entries; cuts (fc); entries", 70, 0, 35);
//					tmp.str("");
//					tmp << "size_cut_entries_s" << sensor << "_s" << s << "_b0";
//					cuts_entries[2][sensor][s] = new TH1F(tmp.str().c_str(), "size_cut_entries; cuts; entries", 35, 0, 35);

//					tmp.str("");
//					tmp << "sigma_cut_performance_s" << sensor << "_s" << s << "_b0";
//					cuts_performance[0][sensor][s] = new TH1F(tmp.str().c_str(), "sigma_cut_performance; cuts (#sigma); entries*purity", 70, 0, 35);
//					tmp.str("");
//					tmp << "charge_cut_performance_s" << sensor << "_s" << s << "_b0";
//					cuts_performance[1][sensor][s] = new TH1F(tmp.str().c_str(), "charge_cut_performance; cuts (fC); entries*purity", 70, 0, 35);
//					tmp.str("");
//					tmp << "size_cut_performance_s" << sensor << "_s" << s << "_b0";
//					cuts_performance[2][sensor][s] = new TH1F(tmp.str().c_str(), "size_cut_performance; cuts; entries*purity", 35, 0, 35);
				}
				
			}
			
			tmp.str("");
			tmp << "Max_SoN_sens" << sensor << "_b0";
			Max_SoN[sensor] = new TH1F(tmp.str().c_str(), "cluster position y0; strip; #Entries", 1840,-0.5, 1839.5);
			
			tmp.str("");
			tmp << "cluster_position_y_seed_sens" << sensor << "_b0";
			cluster_position_y[sensor][0] = new TH1F(tmp.str().c_str(), "cluster position y0; #mum; #Entries", 1840,-46000, 46000);
			tmp.str("");
			tmp << "cluster_charge_seed_sens" << sensor << "_b0";
			cluster_charge[sensor][0] = new TH1F(tmp.str().c_str(), "cluster charge0; Charge (fC); #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance_seed_sens" << sensor << "_b0";
			cluster_significance[sensor][0] = new TH1F(tmp.str().c_str(), "cluster significance_0; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance2_seed_sens" << sensor << "_b0";
			cluster_significance2[sensor][0] = new TH1F(tmp.str().c_str(), "cluster significance2_0; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_sigma_seed_sens" << sensor << "_b0";
			cluster_sigma[sensor][0] = new TH1F(tmp.str().c_str(), "cluster sigmae0; #sigma; #Entries", 200,-0.5, 4.95);
			tmp.str("");
			tmp << "cluster_size_seed_sens" << sensor << "_b0";
			cluster_size[sensor][0] = new TH1F(tmp.str().c_str(), "cluster size0; Size; #Entries", 10,-0.5, 9.5);
			
			
			tmp.str("");
			tmp << "cluster_position_y_all_sens" << sensor << "_b0";
			cluster_position_y[sensor][1] = new TH1F(tmp.str().c_str(), "cluster position y1; #mum; #Entries", 1840,-46000, 46000);
			tmp.str("");
			tmp << "cluster_charge_all_sens" << sensor << "_b0";
			cluster_charge[sensor][1] = new TH1F(tmp.str().c_str(), "cluster charge1; Charge (fC); #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance_all_sens" << sensor << "_b0";
			cluster_significance[sensor][1] = new TH1F(tmp.str().c_str(), "cluster significance_1; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance2_all_sens" << sensor << "_b0";
			cluster_significance2[sensor][1] = new TH1F(tmp.str().c_str(), "cluster significance2_1; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_sigma_all_sens" << sensor << "_b0";
			cluster_sigma[sensor][1] = new TH1F(tmp.str().c_str(), "cluster sigma1; #sigma; #Entries", 200,-0.5, 4.95);
			tmp.str("");
			tmp << "cluster_size_all_sens" << sensor << "_b0";
			cluster_size[sensor][1] = new TH1F(tmp.str().c_str(), "cluster size1; Size; #Entries", 10,-0.5, 9.5);
			tmp.str("");
			tmp << "clusters_all_sens" << sensor << "_b0";
			clusters[sensor][0] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; #Entries", 100,-0.5, 99.5);

			
			
			tmp.str("");
			tmp << "cluster_position_y_cuts_sens" << sensor << "_b0";
			cluster_position_y[sensor][2] = new TH1F(tmp.str().c_str(), "cluster position y2; #mum; #Entries", 1840,-46000, 46000);
			tmp.str("");
			tmp << "cluster_charge_cuts_sens" << sensor << "_b0";
			cluster_charge[sensor][2] = new TH1F(tmp.str().c_str(), "cluster charge2; Charge (fC); #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance_cuts_sens" << sensor << "_b0";
			cluster_significance[sensor][2] = new TH1F(tmp.str().c_str(), "cluster significance_2; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance2_cuts_sens" << sensor << "_b0";
			cluster_significance2[sensor][2] = new TH1F(tmp.str().c_str(), "cluster significance2_2; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_sigma_cuts_sens" << sensor << "_b0";
			cluster_sigma[sensor][2] = new TH1F(tmp.str().c_str(), "cluster sigma2; #sigma; #Entries", 200,-0.5, 4.95);
			tmp.str("");
			tmp << "cluster_size_cuts_sens" << sensor << "_b0";
			cluster_size[sensor][2] = new TH1F(tmp.str().c_str(), "cluster size2; Size; #Entries", 10,-0.5, 9.5);
			tmp.str("");
			tmp << "clusters_cuts_sens" << sensor << "_b0";
			clusters[sensor][1] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; #Entries", 100,-0.5, 99.5);
			
			
//			tmp.str("");
//			tmp << "cluster_position_y_HighSigmaCluster_sens" << sensor << "_b0";
//			cluster_position_y[sensor][3] = n4sigma_1fC_3sizeew TH1F(tmp.str().c_str(), "cluster position y3; #mum; #Entries x Charge", 1840,-46000, 46000);
//			tmp.str("");
//			tmp << "cluster_charge_HighSigmaCluster_sens" << sensor << "_b0";
//			cluster_charge[sensor][3] = new TH1F(tmp.str().c_str(), "cluster charge3; Charge (fC); #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_significance_HighSigmaCluster_sens" << sensor << "_b0";
//			cluster_significance[sensor][3] = new TH1F(tmp.str().c_str(), "cluster significance_3; S/N; #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_significance2_HighSigmaCluster_sens" << sensor << "_b0";
//			cluster_significance2[sensor][3] = new TH1F(tmp.str().c_str(), "cluster significance2_3; S/N; #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_sigma_HighSigmaCluster_sens" << sensor << "_b0";
//			cluster_sigma[sensor][3] = new TH1F(tmp.str().c_str(), "cluster sigma3; #sigma; #Entries", 200,-0.5, 4.95);
//			tmp.str("");
//			tmp << "cluster_size_HighSigmaCluster_sens" << sensor << "_b0";
//			cluster_size[sensor][3] = new TH1F(tmp.str().c_str(), "cluster size3; Size; #Entries", 10,-0.5, 9.5);
			
			
//			tmp.str("");
//			tmp << "cluster_position_y_4sigma_2fC_3size_sens" << sensor << "_b0";
//			cluster_position_y[sensor][4] = new TH1F(tmp.str().c_str(), "cluster position y4; #mum; #Entries x Charge", 1840,-46000, 46000);
//			tmp.str("");
//			tmp << "cluster_charge_4sigma_2fC_3size_sens" << sensor << "_b0";
//			cluster_charge[sensor][4] = new TH1F(tmp.str().c_str(), "cluster charge4; Charge (fC); #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_significance_4sigma_2fC_3size_sens" << sensor << "_b0";
//			cluster_significance[sensor][4] = new TH1F(tmp.str().c_str(), "cluster significance_4; S/N; #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_significance2_4sigma_2fC_3size_sens" << sensor << "_b0";
//			cluster_significance2[sensor][4] = new TH1F(tmp.str().c_str(), "cluster significance2_4; S/N; #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_sigma_4sigma_2fC_3size_sens" << sensor << "_b0";
//			cluster_sigma[sensor][4] = new TH1F(tmp.str().c_str(), "cluster sigma4; #sigma; #Entries", 200,-0.5, 4.95);
//			tmp.str("");
//			tmp << "cluster_size_4sigma_2fC_3size_sens" << sensor << "_b0";
//			cluster_size[sensor][4] = new TH1F(tmp.str().c_str(), "cluster size4; Size; #Entries", 10,-0.5, 9.5);
			
			
//			tmp.str("");
//			tmp << "cluster_position_y_4sigma_3fC_3size_sens" << sensor << "_b0";
//			cluster_position_y[sensor][5] = new TH1F(tmp.str().c_str(), "cluster position y5; #mum; #Entries x Charge", 1840,-46000, 46000);
//			tmp.str("");
//			tmp << "cluster_charge_4sigma_3fC_3size_sens" << sensor << "_b0";
//			cluster_charge[sensor][5] = new TH1F(tmp.str().c_str(), "cluster charge5; Charge (fC); #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_significance_4sigma_3fC_3size_sens" << sensor << "_b0";
//			cluster_significance[sensor][5] = new TH1F(tmp.str().c_str(), "cluster significance_5; S/N; #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_significance2_4sigma_3fC_3size_sens" << sensor << "_b0";
//			cluster_significance2[sensor][5] = new TH1F(tmp.str().c_str(), "cluster significance2_5; S/N; #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_sigma_4sigma_3fC_3size_sens" << sensor << "_b0";
//			cluster_sigma[sensor][5] = new TH1F(tmp.str().c_str(), "cluster sigma5; #sigma; #Entries", 200,-0.5, 4.95);
//			tmp.str("");
//			tmp << "cluster_size_4sigma_3fC_3size_sens" << sensor << "_b0";
//			cluster_size[sensor][5] = new TH1F(tmp.str().c_str(), "cluster size5; Size; #Entries", 10,-0.5, 9.5);
			
			
//			tmp.str("");
//			tmp << "cluster_position_y_4sigma_2_5fC_3size_sens" << sensor << "_b0";
//			cluster_position_y[sensor][6] = new TH1F(tmp.str().c_str(), "cluster position y6; #mum; #Entries x Charge", 1840,-46000, 46000);
//			tmp.str("");
//			tmp << "cluster_charge_4sigma_2_5fC_3size_sens" << sensor << "_b0";
//			cluster_charge[sensor][6] = new TH1F(tmp.str().c_str(), "cluster charge5; Charge (fC); #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_significance_4sigma_2_5fC_3size_sens" << sensor << "_b0";
//			cluster_significance[sensor][6] = new TH1F(tmp.str().c_str(), "cluster significance_6; S/N; #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_significance2_4sigma_2_5fC_3size_sens" << sensor << "_b0";
//			cluster_significance2[sensor][6] = new TH1F(tmp.str().c_str(), "cluster significance2_6; S/N; #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_sigma_4sigma_2_5fC_3size_sens" << sensor << "_b0";
//			cluster_sigma[sensor][6] = new TH1F(tmp.str().c_str(), "cluster sigma5; #sigma; #Entries", 200,-0.5, 4.95);
//			tmp.str("");
//			tmp << "cluster_size_4sigma_2_5fC_3size_sens" << sensor << "_b0";
//			cluster_size[sensor][6] = new TH1F(tmp.str().c_str(), "cluster size5; Size; #Entries", 10,-0.5, 9.5);
			
			
//			tmp.str("");
//			tmp << "cluster_position_y_test_sens" << sensor << "_b0";
//			cluster_position_y[sensor]3[7] = new TH1F(tmp.str().c_str(), "cluster position y6; #mum; #Entries x Charge", 1840,-46000, 46000);
//			tmp.str("");
//			tmp << "cluster_charge_test_sens" << sensor << "_b0";
//			cluster_charge[sensor][7] = new TH1F(tmp.str().c_str(), "cluster charge5; Charge (fC); #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_significance_test_sens" << sensor << "_b0";
//			cluster_significance[sensor][7] = new TH1F(tmp.str().c_str(), "cluster significance_7; S/N; #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_significance2_test_sens" << sensor << "_b0";
//			cluster_significance2[sensor][7] = new TH1F(tmp.str().c_str(), "cluster significance2_7; S/N; #Entries", 200,-0.5, 49.5);
//			tmp.str("");
//			tmp << "cluster_sigma_test_sens" << sensor << "_b0";
//			cluster_sigma[sensor][7] = new TH1F(tmp.str().c_str(), "cluster sigma5; #sigma; #Entries", 200,-0.5, 4.95);
//			tmp.str("");
//			tmp << "cluster_size_test_sens" << sensor << "_b0";
//			cluster_size[sensor][7] = new TH1F(tmp.str().c_str(), "cluster size5; Size; #Entries", 10,-0.5, 9.5);
			
			
			
			//tmp.str("");
			//tmp << "cluster_position_x0_sens" << sensor << "_b0";
			//cluster_position_x[sensor][0] = new TH1F(tmp.str().c_str(), "cluster position x0; #mum; #Entries", 1840,-0.5, 91999.5);
			//tmp.str("");
			//tmp << "cluster_position_x1_sens" << sensor << "_b0";
			//cluster_position_x[sensor][1] = new TH1F(tmp.str().c_str(), "cluster position x1; #mum; #Entries", 1840,-0.5, 91999.5);
			//tmp.str("");
			//tmp << "cluster_position_x2_sens" << sensor << "_b0";
			//cluster_position_x[sensor][2] = new TH1F(tmp.str().c_str(), "cluster position x2; #mum; #Entries", 1840,-0.5, 91999.5);
			//tmp.str("");
			//tmp << "cluster_position_x3_sens" << sensor << "_b0";
			//cluster_position_x[sensor][3] = new TH1F(tmp.str().c_str(), "cluster position x3; #mum; #Entries x Charge", 1840,-0.5, 91999.5);
			


			
			
			tmp.str("");
			tmp << "noise_v_position_s" << sensor << "_b0";
			noise_v_position[sensor]  = new TH1F(tmp.str().c_str(), "Noise; #mum; Noise (fC)", 1840,-46000, 46000);
			
			tmp.str("");
			tmp << "noise_distribution_s" << sensor << "_b0";
			noise_distribution_sensor[sensor] = new TH1F(tmp.str().c_str(), "noise_distribution; Noise(fC);   #channels", 100,-0.005, 0.995);
			
			tmp.str("");
			tmp << "noise_v_time_s" << sensor << "_b0";
			tmp_units.str("");
			tmp_units << "noise_v_time; Time (BCC); Noise (fC) ";
			noise_v_time[sensor]  = new TH1F(tmp.str().c_str(), tmp_units.str().c_str(), 8192, 0, 8191);
			
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
			
				}
			}
		}
		
	}
	//////////////////////////////////////////
	// Data read for all events for detailed look
	//////////////////////////////////////////
	dataRead.open(argv[1]); //open file again to start from the beginning
	
	int header = 1;
	

	int clstrcounter[n_kpix][1840] = {0};
	
	
	std::vector<double>* corrected_charge_vec[n_kpix][n_channels] = {nullptr};
	std::vector<double>* corrected_charge_vec_time[n_kpix][n_BCC]  = {nullptr};
	
	unordered_map<uint, uint> noise_mask[6];
	noise_mask[0] = noise_sensor_0();
	noise_mask[1] = noise_sensor_1();
	noise_mask[2] = noise_sensor_2();
	noise_mask[3] = noise_sensor_3();
	noise_mask[4] = noise_sensor_4();
	noise_mask[5] = noise_sensor_5();
	
	while ( dataRead.next(&event) ) //preread to determine noise value of each channel in each KPiX.
	{
		int not_empty= 0;
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
			sensor = kpix/2;
			unsigned int block = channel/32;
			bunchClk = sample->getBunchCount();
			subCount = sample->getSubCount();
		
			if ( type == KpixSample::Data ) // If event is of type KPiX data
			{
				if (bucket == 0)
				{
					if (pedestal_MedMAD[kpix][channel][bucket][1] != 0 && calib_slope[kpix][channel] != 0) //ensuring we ignore 0 MAD channels
					{
						//cout << "DEBUG 1" << endl;
						double charge_value = double(value)/calib_slope[kpix][channel];
						double corrected_charge_value_median = charge_value - pedestal_MedMAD[kpix][channel][bucket][0];
						double charge_CM_corrected = corrected_charge_value_median - common_modes_median[kpix].at(event.eventNumber());
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
	for (unsigned int k = 0; k < n_kpix; ++k)
	{
		if (kpixFound[k])
		{
			int sensor = k/2;
			for (unsigned int c = 0; c < n_channels; ++c)
			{
				if (channelFound[k][c])
				{
					int strip = 9999;
					if (k%2 == 0) strip = kpix2strip_left.at(c);// if left kpix
					else strip  = kpix2strip_right.at(c); // if right kpix
					if (pedestal_MedMAD[k][c][0][1] != 0)
					{
						y = yParameterSensor(strip, sensor);
						noise[k][c] = 1.4826*MAD(corrected_charge_vec[k][c]);
						//cout << "Noise of kpix " << k << " and channel " << c << " is " << noise[k][c] << endl;
						if (strip != 9999)
						{
							noise_distribution[k]->Fill(noise[k][c]);
							noise_distribution_sensor[sensor]->Fill(noise[k][c]);
							noise_v_position[sensor]->Fill(y, noise[k][c]);
						}
						noise_v_channel[k]->SetBinContent(c+1, noise[k][c]);
						
					}
				}
				
			}
			for (unsigned int t = 0; t < n_BCC ; ++t)
			{
				noise_v_time[sensor]->SetBinContent(t+1, 1.4826*MAD(corrected_charge_vec_time[k][t]));
			}	
		}
	}
	
	/////////////////////////////////////////////////////////
	//// BEGIN OF CLUSTER READ
	dataRead.open(argv[1]);
	
	int global_trig_counter = 0;
//	std::vector<clustr> all_clusters[n_kpix/2];
	std::vector<clustr>* all_clusters_pointer[n_kpix/2][acqCount] = {nullptr};
	while ( dataRead.next(&event) )
	{
		int not_empty= 0;
		std::vector<double> time_ext;
		std::map<int, double> cluster_Events_after_cut[n_kpix/2];
		std::map<int, double> cluster_Noise_after_cut[n_kpix/2];
		
		frameruntime = event.runtime();

		//cout << "Beginning a new EVENT" << endl;
		//cout << " NEW EVENT " << endl;
		int trig_counter = 0;
		
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
			sensor = kpix/2;
			bunchClk = sample->getBunchCount();
			subCount = sample->getSubCount();
			
			//channel to strip assignment.
			int strip = 9999;
			if (kpix%2 == 0) strip = kpix2strip_left.at(channel);// if left kpix
			else strip  = kpix2strip_right.at(channel); // if right kpix

			y = yParameter(strip, kpix); //Testbeam position parameter. Already adjusted for flipped sensors and kpix 1/2 position. NOT for the stereo angle
	
			if (type == KpixSample::Timestamp)// If event is of type external timestamp
			{
				trig_counter++;
				global_trig_counter++;
				double time = bunchClk + double(subCount * 0.125);
				time_ext.push_back(time);
				if (frameruntime!=0)
					runtime = sample->getSampleRuntime64(frameruntime);
				else
					cerr<< "Error: frameruntime is ZEROs!"<< endl;
				if (trig_counter == 1)
				{
					//cout << "DEBUG: " << event.eventNumber() << " ," << time << " ," << trig_counter << endl;
					tmp.str("");
					tmp << setw(7) << time  << " ," << setw(5) << global_trig_counter << runtime ;
					
				}
			}
			if ( type == KpixSample::Data ) // If event is of type KPiX data
			{
				
				if (sample->getEmpty()) cout << "Event is empty and I am analysing it" << endl; //debug output whether an empty event is still getting read
				else not_empty = 1; //checkmark that event is not empty

				//cout << tstamp << endl;
				if (bucket == 0)
				{
					if (pedestal_MedMAD[kpix][channel][bucket][1] != 0 && calib_slope[kpix][channel] != 0) //ensuring we ignore 0 MAD channels and channels with bad 0 slopes
					{
						
						//// ====== Calculation of Charge values, with pedestal and common mode subtraction  =============
						double charge_value = double(value)/calib_slope[kpix][channel];
						double corrected_charge_value_median = charge_value - pedestal_MedMAD[kpix][channel][bucket][0];
						double charge_CM_corrected = corrected_charge_value_median - common_modes_median[kpix].at(event.eventNumber());
						//// ========= Event cut ============
                
						fc_response_medCM_subtracted[kpix]->Fill(charge_CM_corrected);
//						if ( charge_CM_corrected > 3*noise[kpix][channel] && strip != 9999)
						if ( charge_CM_corrected > 3*noise[kpix][channel] && strip != 9999 && noise_mask[sensor].at(strip) == 1)  //only events with charge higher than 2 sigma of the noise are taken and with their charge being lower than 10 fC (to cut out weird channels), in addition no noise masked channels and no disconnected channels.
						{
							cluster_Events_after_cut[sensor].insert(std::pair<int, double>(strip, charge_CM_corrected));
							cluster_Noise_after_cut[sensor].insert(std::pair<int, double>(strip, noise[kpix][channel]));
							if (noise[kpix][channel] == 0) cout << "Something is going wrong here" << endl;
						}
						//if (kpix == 0 && (channel == 9 || channel == 105)) cout << "Something is going wrong here" << endl;
						
					}
				}
			}
			
			

		}
		if (not_empty == 1) //if event is not empty
		{
			//cout << "DEBUG 1" << endl;
			clustr Cluster[n_kpix/2];  // Another Cluster class variable
			std::vector<clustr> multi_cluster[8][n_kpix/2];
			for (sensor = 0; sensor < n_kpix/2; sensor++)
			{
				if (kpixFound[2*sensor] || kpixFound[2*sensor+1])
				{
					std::vector<clustr>::iterator it_SoN;
					if ( cluster_Events_after_cut[sensor].size() != 0 )
					{
						//cout << "===================" << endl;
						//cout << "Starting new PacMan" << endl;
						//cout << "===================" << endl;
						//cout << endl;
						//cout << "DEBUG 2" << endl;
						clustr Input;
						Input.Elements = cluster_Events_after_cut[sensor];
						Input.Noise = cluster_Noise_after_cut[sensor];
						int num_of_clusters = 0;
						if (all_clusters_pointer[sensor][event.eventNumber()] == nullptr)
						{
							all_clusters_pointer[sensor][event.eventNumber()] = new std::vector<clustr>;
							if (all_clusters_pointer[sensor][event.eventNumber()]==nullptr)
							{
								std::cerr << "Memory allocation error for vector kpix " <<
											 kpix <<std::endl;
								exit(-1); // probably best to bail out
							}
						}

						while (Input.Elements.size() != 0) // Keep repeating the clustering until either there are no valid candidates left
						{
							PacMan NomNom;
							//cout << "Maximum Signal over Noise Strip " << Input.MaxSoN() << endl;
							//cout << "Maximum Charge Strip " << Input.MaxCharge_w_Noise() << endl;
							double SoN_order = 0;
							Max_SoN[sensor]->Fill(Input.MaxSoN(), weight);
							MaximumSoN[sensor][Input.MaxSoN()]+=weight;
							NomNom.Eater(Input, Input.MaxSoN(), 9999, 99999);

							if (num_of_clusters == 0)
							{
								cluster_position_y[sensor][0]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), weight);
								cluster_charge[sensor][0]->Fill(NomNom.getClusterCharge(), weight);
								cluster_size[sensor][0]->Fill(NomNom.getClusterElementssize(), weight);
								cluster_significance[sensor][0]->Fill(NomNom.getClusterSignificance(), weight);
								cluster_significance2[sensor][0]->Fill(NomNom.getClusterSignificance2(), weight);
								cluster_sigma[sensor][0]->Fill(NomNom.getClusterSigma(), weight);
								multi_cluster[0][sensor].push_back(NomNom.getCluster());

								//cout << "Size: " << NomNom.getClusterElementssize() << endl;
							}

							cluster_position_y[sensor][1]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), weight);

							cluster_size[sensor][1]->Fill(NomNom.getClusterElementssize(), weight);
							cluster_significance[sensor][1]->Fill(NomNom.getClusterSignificance(), weight);
							cluster_significance2[sensor][1]->Fill(NomNom.getClusterSignificance2(), weight);
							cluster_charge[sensor][1]->Fill(NomNom.getClusterCharge(), weight);
							cluster_sigma[sensor][1]->Fill(NomNom.getClusterSigma(), weight);
							multi_cluster[1][sensor].push_back(NomNom.getCluster());

							//						all_clusters[sensor].push_back(NomNom.getCluster());
//							cout << "Sensor number: " << sensor << endl;
//							all_clusters_pointer[sensor][event.eventNumber()]->push_back(NomNom.getCluster());


							if (NomNom.getClusterSignificance2() >= 7.0 && NomNom.getClusterCharge() >= 2.0 && NomNom.getClusterElementssize() <= 2)
							{
								cluster_position_y[sensor][2]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), weight);
								cluster_charge[sensor][2]->Fill(NomNom.getClusterCharge(), weight);
								cluster_size[sensor][2]->Fill(NomNom.getClusterElementssize(), weight);
								cluster_significance[sensor][2]->Fill(NomNom.getClusterSignificance(), weight);
								cluster_significance2[sensor][2]->Fill(NomNom.getClusterSignificance2(), weight);
								cluster_sigma[sensor][2]->Fill(NomNom.getClusterSigma(), weight);
								multi_cluster[2][sensor].push_back(NomNom.getCluster());
							}
							//						if (NomNom.getClusterSignificance2() > 4.0 && NomNom.getClusterCharge() > 2.0 && NomNom.getClusterElementssize() < 3)
							//						{
							//							cluster_position_y[sensor][4]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), weight);
							//							cluster_charge[sensor][4]->Fill(NomNom.getClusterCharge(), weight);
							//							cluster_size[sensor][4]->Fill(NomNom.getClusterElementssize(), weight);
							//							cluster_significance[sensor][4]->Fill(NomNom.getClusterSignificance(), weight);
							//							cluster_significance2[sensor][4]->Fill(NomNom.getClusterSignificance2(), weight);
							//							cluster_sigma[sensor][4]->Fill(NomNom.getClusterSigma(), weight);
							//							multi_cluster[4][sensor].push_back(NomNom.getCluster());
							//						}
							//						if (NomNom.getClusterSignificance2() > 4.0 && NomNom.getClusterCharge() > 3.0 && NomNom.getClusterElementssize() < 3)
							//						{
							//							cluster_position_y[sensor][5]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), weight);
							//							cluster_charge[sensor][5]->Fill(NomNom.getClusterCharge(), weight);
							//							cluster_size[sensor][5]->Fill(NomNom.getClusterElementssize(), weight);
							//							cluster_significance[sensor][5]->Fill(NomNom.getClusterSignificance(), weight);
							//							cluster_significance2[sensor][5]->Fill(NomNom.getClusterSignificance2(), weight);
							//							cluster_sigma[sensor][5]->Fill(NomNom.getClusterSigma(), weight);
							//							multi_cluster[5][sensor].push_back(NomNom.getCluster());
							//						}
							//						if (NomNom.getClusterSignificance2() > 4.0 && NomNom.getClusterCharge() > 2.5 && NomNom.getClusterElementssize() < 3)
							//						{
							//							cluster_position_y[sensor][6]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), weight);
							//							cluster_charge[sensor][6]->Fill(NomNom.getClusterCharge(), weight);
							//							cluster_size[sensor][6]->Fill(NomNom.getClusterElementssize(), weight);
							//							cluster_significance[sensor][6]->Fill(NomNom.getClusterSignificance(), weight);
							//							cluster_significance2[sensor][6]->Fill(NomNom.getClusterSignificance2(), weight);
							//							cluster_sigma[sensor][6]->Fill(NomNom.getClusterSigma(), weight);
							//							multi_cluster[6][sensor].push_back(NomNom.getCluster());
							//						}
							//						if (NomNom.getClusterSignificance2() > 4.0 && NomNom.getClusterCharge() > 2.0)
							//						{
							//							cluster_position_y[sensor][7]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), weight);
							//							cluster_charge[sensor][7]->Fill(NomNom.getClusterCharge(), weight);
							//							cluster_size[sensor][7]->Fill(NomNom.getClusterElementssize(), weight);
							//							cluster_significance[sensor][7]->Fill(NomNom.getClusterSignificance(), weight);
							//							cluster_significance2[sensor][7]->Fill(NomNom.getClusterSignificance2(), weight);
							//							cluster_sigma[sensor][7]->Fill(NomNom.getClusterSigma(), weight);
							//							multi_cluster[7][sensor].push_back(NomNom.getCluster());
							//						}
							if (header == 1){
								header = 0;
								claus_file <<"Event Number,Layer,position,Significance,Significance2,Size,Charge,time, #trig" << endl;
							}
							claus_file << setw(5) << event.eventNumber()  << ","
							           << setw(1) << sensor2layer.at(sensor)  << ","
							           << setw(7) << yParameterSensor(NomNom.getClusterCoG(), sensor)  << ","
							           << setw(7) << NomNom.getClusterSignificance() << ","
							           << setw(2) << NomNom.getClusterElementssize() << ","
							           << setw(7) << NomNom.getClusterSignificance2() << ","
							           << setw(7) << NomNom.getClusterCharge() << ","
							           << tmp.str().c_str()
							           << endl;

							num_of_clusters++;
							//cout << "Cluster Position is " << yParameterSensor(NomNom.getClusterCoG(), sensor) << endl;
						}
						clusters[sensor][0]->Fill(num_of_clusters);
						clusters[sensor][1]->Fill(multi_cluster[2][sensor].size());


						sort( multi_cluster[1][sensor].begin(), multi_cluster[1][sensor].end(), [ ]( const clustr& lhs, const clustr& rhs )  // lambda expression sorting of the clusters by their significance, first element is most significant (very experimental as I do not fully understand lambda expressions, I only interpret them...)
						{
							return lhs.Significance > rhs.Significance;
						});
						//should now only fill the most significant cluster
						//					if ( multi_cluster[0][sensor].size() != 0)
						//					{
						//						cluster_position_y[sensor][3]->Fill(yParameterSensor(multi_cluster[1][sensor].at(0).CoG, sensor), weight);
						//						cluster_size[sensor][3]->Fill(multi_cluster[1][sensor].at(0).Elements.size(), weight);
						//						cluster_significance[sensor][3]->Fill(multi_cluster[1][sensor].at(0).Significance, weight);
						//						cluster_significance2[sensor][3]->Fill(multi_cluster[1][sensor].at(0).Significance2, weight);
						//						cluster_charge[sensor][3]->Fill(multi_cluster[1][sensor].at(0).Charge, weight);
						//						cluster_sigma[sensor][3]->Fill(multi_cluster[1][sensor].at(0).Sigma, weight);
						//					}
					}
				}
				for (int sensor1 = 0; sensor1 < n_kpix/2; sensor1++)
				{
					for (int sensor2 = sensor1+1; sensor2 < n_kpix/2; sensor2++)
					{
						for (auto const& s1 : multi_cluster[0][sensor1])
						{
							for (auto const& s2 : multi_cluster[0][sensor2])
							{
								double y1 = yParameterSensor(s1.CoG, sensor1);
								double y2 = yParameterSensor(s2.CoG, sensor2);
								double clstroffset_y  = y1 - y2;
								cluster_offset_y[0][sensor1][sensor2]->Fill(clstroffset_y, weight);
								cluster_correlation[0][sensor1][sensor2]->Fill(y1,y2);
							}
						}
						for (auto const& s1 : multi_cluster[1][sensor1])
						{
							for (auto const& s2 : multi_cluster[1][sensor2])
							{
								double y1 = yParameterSensor(s1.CoG, sensor1);
								double y2 = yParameterSensor(s2.CoG, sensor2);
								double clstroffset_y  = y1 - y2;
								cluster_offset_y[1][sensor1][sensor2]->Fill(clstroffset_y, weight);
								cluster_correlation[1][sensor1][sensor2]->Fill(y1,y2);
							}
						}
						for (auto const& s1 : multi_cluster[2][sensor1])
						{
							for (auto const& s2 : multi_cluster[2][sensor2])
							{
								double y1 = yParameterSensor(s1.CoG, sensor1);
								double y2 = yParameterSensor(s2.CoG, sensor2);
								double clstroffset_y  = y1 - y2;
								cluster_offset_y[2][sensor1][sensor2]->Fill(clstroffset_y, weight);
								cluster_correlation[2][sensor1][sensor2]->Fill(y1,y2);
							}
						}
						//					for (auto const& s1 : multi_cluster[3][sensor1])
						//					{
						//						for (auto const& s2 : multi_cluster[3][sensor2])
						//						{
						//							double y1 = yParameterSensor(s1.CoG, sensor1);
						//							double y2 = yParameterSensor(s2.CoG, sensor2);
						//							double clstroffset_y  = y1 - y2;
						//							cluster_offset_y[3][sensor1][sensor2]->Fill(clstroffset_y, weight);
						//							cluster_correlation[3][sensor1][sensor2]->Fill(y1,y2);
						//						}
						//					}
						//					for (auto const& s1 : multi_cluster[4][sensor1])
						//					{
						//						for (auto const& s2 : multi_cluster[4][sensor2])
						//						{
						//							double y1 = yParameterSensor(s1.CoG, sensor1);
						//							double y2 = yParameterSensor(s2.CoG, sensor2);
						//							double clstroffset_y  = y1 - y2;
						//							cluster_offset_y[4][sensor1][sensor2]->Fill(clstroffset_y, weight);
						//							cluster_correlation[4][sensor1][sensor2]->Fill(y1,y2);
						//						}
						//					}
						//					for (auto const& s1 : multi_cluster[5][sensor1])
						//					{
						//						for (auto const& s2 : multi_cluster[5][sensor2])
						//						{
						//							double y1 = yParameterSensor(s1.CoG, sensor1);
						//							double y2 = yParameterSensor(s2.CoG, sensor2);
						//							double clstroffset_y  = y1 - y2;
						//							cluster_offset_y[5][sensor1][sensor2]->Fill(clstroffset_y, weight);
						//							cluster_correlation[5][sensor1][sensor2]->Fill(y1,y2);
						//						}
						//					}
						//					for (auto const& s1 : multi_cluster[6][sensor1])
						//					{
						//						for (auto const& s2 : multi_cluster[6][sensor2])
						//						{
						//							double y1 = yParameterSensor(s1.CoG, sensor1);
						//							double y2 = yParameterSensor(s2.CoG, sensor2);
						//							double clstroffset_y  = y1 - y2;
						//							cluster_offset_y[6][sensor1][sensor2]->Fill(clstroffset_y, weight);
						//							cluster_correlation[6][sensor1][sensor2]->Fill(y1,y2);
						//						}
						//					}
						//					for (auto const& s1 : multi_cluster[7][sensor1])
						//					{
						//						for (auto const& s2 : multi_cluster[7][sensor2])
						//						{
						//							double y1 = yParameterSensor(s1.CoG, sensor1);
						//							double y2 = yParameterSensor(s2.CoG, sensor2);
						//							double clstroffset_y  = y1 - y2;
						//							cluster_offset_y[7][sensor1][sensor2]->Fill(clstroffset_y, weight);
						//							cluster_correlation[7][sensor1][sensor2]->Fill(y1,y2);
						//						}
						//					}
					}

				}
			}
		}
	////   Show progress
		filePos  = dataRead.pos();
		currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
		if ( currPct != lastPct ) 
		{
			cout << "\rReading File for cluster analysis: " << currPct << " %      " << flush;
			lastPct = currPct;
		}
		
	}
	//// Noise mask generation. ONLY USE WHEN NO MASK IS PUT IN!
//	noise_file.open("include/testbeam201907_noise_mask.h");
//	for (int s = 0; s < 6; ++s)
//	{
//		int badentries = 0;
//		noise_file << "unordered_map<uint, uint> noise_sensor_" << s << "()" << endl;
//		noise_file << "{" << endl;
//		noise_file << "    unordered_map<uint, uint> m1;" << endl;
//		for (int strips = 0; strips < 1840; strips++)
//		{
//			if (MaximumSoN[s][strips] > 0.08 ) //
//			{
//				cout << "Huge number of entries with no adjacent entries in strip " << strips << " with " << MaximumSoN[s][strips] << " entries in Sensor " << s << endl;
//				badentries ++;
//				noise_file << "    m1.insert(make_pair(" << strips << ",0));" << endl;
//			}
//			else
//			{
//				noise_file << "    m1.insert(make_pair(" << strips << ",1));" << endl;
//			}
//			//if (Max_SoN[s]->GetBinContent(strips) > 25)
//			//{
//				//cout << "Huge number of entries with no adjacent entries in strip in MaxSoN " << strips << " with " << Max_SoN[s]->GetBinContent(strips) << " entries in Sensor " << s << endl;
//			//}
//		}
//		noise_file << "    return m1;" << endl;
//		noise_file << "}" << endl;
//		cout << "Number of bad entries in this sensor is " << badentries << endl;
//	}
//	noise_file.close();



//	for (unsigned int sens = 0; sens<n_kpix/2; sens++ )
//	{
//		for (unsigned int sens2 = sens+1; sens2<n_kpix/2; sens2++ )
//		{
//			if (kpixFound[sens*2] && kpixFound[sens2*2])
//			{
////				TH1F* cut_performance = nullptr;
////				TH1F* cut_entries = nullptr;
////				TCanvas *canvas2 = nullptr;
////				TCanvas *canvas3 = nullptr;
////				tmp.str("");
////				tmp << "sigma_cut_performance_s" << sens << "_s" << sens2 << "_b0";
////				cut_performance = new TH1F(tmp.str().c_str(), "cut_performance; cuts; efficiency", 416, 0, 416);
////				canvas2 = new TCanvas(tmp.str().c_str(), "cut_performance", 3200, 1800);

////				tmp.str("");
////				tmp << "sigma_cut_entries_s" << sens << "_s" << sens2 << "_b0";
////				cut_entries = new TH1F(tmp.str().c_str(), "cut_entries; cuts; entries", 60, 0, 30);
////				canvas3 = new TCanvas(tmp.str().c_str(), "cut_entries", 3200, 1800);
//				for (int s = 0; s < 30; ++s)
//				{
//					double sigma_cut = double(s)/2;
//					double charge_cut = 2.0;
//					int size_cut =3;
//					tmp.str("");
//					tmp << "cluster_offset_y_s" << sens << "_s" << sens2 << "_b0" << "_SigmaCutG_" << sigma_cut << "_ChargeCutG_" << charge_cut << "_SizeCutL_" << size_cut;
//					TH1F* test = new TH1F(tmp.str().c_str(), "cluster offset test; #mum; #Entries", 100,-10000, 10000);
////					TCanvas *canvas1 = new TCanvas(tmp.str().c_str(), "cluster_offset", 3200, 1800);
////					canvas1->cd();
//					for (unsigned int evt = 0; evt < acqCount; evt++)
//					{
//						//cout << "Evt " << evt << endl;
//						if (all_clusters_pointer[sens][evt] != nullptr && all_clusters_pointer[sens2][evt] != nullptr)
//						{
//							for (auto const s1 : *all_clusters_pointer[sens][evt])
//							{
//								if (s1.Charge >= charge_cut && s1.Significance2 >= sigma_cut && s1.Elements.size() <= size_cut )
//								{
//									for (auto const s2 : *all_clusters_pointer[sens2][evt])
//									{
//										if (s2.Charge >= charge_cut && s2.Significance2 >= sigma_cut && s2.Elements.size() <= size_cut )
//										{
//											double y1 = yParameterSensor(s1.CoG, sens);
//											double y2 = yParameterSensor(s2.CoG, sens2);
//											double clstroffset_y  = y1 - y2;
//											test->Fill(clstroffset_y);
//										}
//									}
//								}
//							}
//						}
//					}
//					test->Draw("pe");

//					int binmax = test->GetMaximumBin();
//					double signal = test->Integral(binmax-4, binmax+4);
//					double entries = test->GetEntries();
//					double efficiency = 0;
//					if (entries != 0)
//						efficiency = signal/entries;
////					cut_performance->Fill(sigma_cut, efficiency);
////					cut_entries->Fill(sigma_cut, entries);
//					cuts_entries[0][sens][sens2]->SetBinContent(s, entries);
//					cuts_purity[0][sens][sens2]->SetBinContent(s, efficiency);
//					cuts_performance[0][sens][sens2]->SetBinContent(s, efficiency*entries);
////					stringstream tmp2;
////					tmp2.str("");
////					tmp2 << "/home/lycoris-dev/Documents/cut_test/" << tmp.str() << ".png";
////					canvas1->SaveAs(tmp2.str().c_str());
//					test = nullptr;
//					delete test;




//				}
////				canvas2->cd();
////				cut_performance->Draw("p");
////				canvas2->SetGridx();
////				tmp.str("");
////				tmp << "/home/lycoris-dev/Documents/cut_test/sigma_cut_performance_s" << sens << "_s" << sens2 <<".pdf";
////				canvas2->SaveAs(tmp.str().c_str());

////				canvas3->cd();
////				cut_entries->Draw("p");
////				canvas3->SetGridx();
////				tmp.str("");
////				tmp << "/home/lycoris-dev/Documents/cut_test/sigma_cut_entries_s" << sens << "_s" << sens2 <<".pdf";
////				canvas3->SaveAs(tmp.str().c_str());


////				tmp.str("");
////				tmp << "charge_cut_performance_s" << sens << "_s" << sens2 << "_b0";
////				cut_performance = new TH1F(tmp.str().c_str(), "cut_performance; cuts; efficiency", 40, 0, 20);
////				cut_performance->SetCanExtend(TH1::kAllAxes);
////				canvas2 = new TCanvas(tmp.str().c_str(), "cut_performance", 3200, 1800);

////				tmp.str("");
////				tmp << "charge_cut_entries_s" << sens << "_s" << sens2 << "_b0";
////				cut_entries = new TH1F(tmp.str().c_str(), "cut_entries; cuts; entries", 40, 0, 20);
////				canvas3 = new TCanvas(tmp.str().c_str(), "cut_entries", 3200, 1800);


//				for (int c = 0; c < 20; ++c)
//				{
//					double charge_cut = double(c)/2;
//					double sigma_cut = 7.0;
//					int size_cut = 3;

//					tmp.str("");
//					tmp << "cluster_offset_y_s" << sens << "_s" << sens2 << "_b0" << "_SigmaCutG_" << sigma_cut << "_ChargeCutG_" << charge_cut << "_SizeCutL_" << size_cut;
//					TH1F* test = new TH1F(tmp.str().c_str(), "cluster offset test; #mum; #Entries", 100,-10000, 10000);
////					TCanvas *canvas1 = new TCanvas(tmp.str().c_str(), "cluster_offset", 3200, 1800);
////					canvas1->cd();
//					for (unsigned int evt = 0; evt < acqCount; evt++)
//					{
//						//cout << "Evt " << evt << endl;
//						if (all_clusters_pointer[sens][evt] != nullptr && all_clusters_pointer[sens2][evt] != nullptr)
//						{
//							for (auto const s1 : *all_clusters_pointer[sens][evt])
//							{
//								if (s1.Charge >= charge_cut && s1.Significance2 >= sigma_cut && s1.Elements.size() <= size_cut )
//								{
//									for (auto const s2 : *all_clusters_pointer[sens2][evt])
//									{
//										if (s2.Charge >= charge_cut && s2.Significance2 >= sigma_cut && s2.Elements.size() <= size_cut )
//										{
//											double y1 = yParameterSensor(s1.CoG, sens);
//											double y2 = yParameterSensor(s2.CoG, sens2);
//											double clstroffset_y  = y1 - y2;
//											test->Fill(clstroffset_y);
//										}
//									}
//								}
//							}
//						}
//					}
//					test->Draw("pe");

//					int binmax = test->GetMaximumBin();
//					double signal = test->Integral(binmax-4, binmax+4);
//					double entries = test->GetEntries();
//					double efficiency = 0;
//					if (entries != 0)
//						efficiency = signal/entries;
////					cut_performance->Fill(charge_cut, efficiency);
////					cut_entries->Fill(charge_cut, entries);
//					cuts_entries[1][sens][sens2]->SetBinContent(c, entries);
//					cuts_purity[1][sens][sens2]->SetBinContent(c, efficiency);
//					cuts_performance[1][sens][sens2]->SetBinContent(c, efficiency*entries);
////					stringstream tmp2;
////					tmp2.str("");
////					tmp2 << "/home/lycoris-dev/Documents/cut_test/" << tmp.str() << ".png";
////					canvas1->SaveAs(tmp2.str().c_str());
//					test = nullptr;
//					delete test;
//				}
////				canvas2->cd();
////				cut_performance->Draw("p");
////				canvas2->SetGridx();
////				tmp.str("");
////				tmp << "/home/lycoris-dev/Documents/cut_test/charge_cut_performance_s" << sens << "_s" << sens2 <<".pdf";
////				canvas2->SaveAs(tmp.str().c_str());

////				canvas3->cd();
////				cut_entries->Draw("p");
////				canvas3->SetGridx();
////				tmp.str("");
////				tmp << "/home/lycoris-dev/Documents/cut_test/charge_cut_entries_s" << sens << "_s" << sens2 <<".pdf";
////				canvas3->SaveAs(tmp.str().c_str());




////				tmp.str("");
////				tmp << "size_cut_performance_s" << sens << "_s" << sens2 << "_b0";
////				cut_performance = new TH1F(tmp.str().c_str(), "cut_performance; cuts; efficiency", 416, 0, 416);
////				cut_performance->SetCanExtend(TH1::kAllAxes);
////				canvas2 = new TCanvas(tmp.str().c_str(), "cut_performance", 3200, 1800);

////				tmp.str("");
////				tmp << "size_cut_entries_s" << sens << "_s" << sens2 << "_b0";
////				cut_entries = new TH1F(tmp.str().c_str(), "cut_entries; cuts; entries", 12, 0, 12);
////				canvas3 = new TCanvas(tmp.str().c_str(), "cut_entries", 3200, 1800);

//				for (int size_cut = 0; size_cut < 12; ++size_cut)
//				{
//					double charge_cut = 2.0;
//					double sigma_cut = 7.0;
//					tmp.str("");
//					tmp << "cluster_offset_y_s" << sens << "_s" << sens2 << "_b0" << "_SigmaCutG_" << sigma_cut << "_ChargeCutG_" << charge_cut << "_SizeCutL_" << size_cut;
//					TH1F* test = new TH1F(tmp.str().c_str(), "cluster offset test; #mum; #Entries", 100,-10000, 10000);
////					TCanvas *canvas1 = new TCanvas(tmp.str().c_str(), "cluster_offset", 3200, 1800);
////					canvas1->cd();
//					for (unsigned int evt = 0; evt < acqCount; evt++)
//					{
//						//cout << "Evt " << evt << endl;
//						if (all_clusters_pointer[sens][evt] != nullptr && all_clusters_pointer[sens2][evt] != nullptr)
//						{
//							for (auto const s1 : *all_clusters_pointer[sens][evt])
//							{
//								if (s1.Charge >= charge_cut && s1.Significance2 >= sigma_cut && s1.Elements.size() <= size_cut )
//								{
//									for (auto const s2 : *all_clusters_pointer[sens2][evt])
//									{
//										if (s2.Charge >= charge_cut && s2.Significance2 >= sigma_cut && s2.Elements.size() <= size_cut )
//										{
//											double y1 = yParameterSensor(s1.CoG, sens);
//											double y2 = yParameterSensor(s2.CoG, sens2);
//											double clstroffset_y  = y1 - y2;
//											test->Fill(clstroffset_y);
//										}
//									}
//								}
//							}
//						}
//					}
//					test->Draw("pe");

//					int binmax = test->GetMaximumBin();
//					double signal = test->Integral(binmax-1, binmax+1);
//					double entries = test->GetEntries();
//					double efficiency = 0;
//					if (entries != 0)
//						efficiency = signal/entries;
////					cut_performance->Fill(size_cut, efficiency);
////					cut_entries->Fill(size_cut, entries);
//					cuts_entries[2][sens][sens2]->SetBinContent(size_cut, entries);
//					cuts_purity[2][sens][sens2]->SetBinContent(size_cut, efficiency);
//					cuts_performance[2][sens][sens2]->SetBinContent(size_cut, efficiency*entries);
////					stringstream tmp2;
////					tmp2.str("");
////					tmp2 << "/home/lycoris-dev/Documents/cut_test/" << tmp.str() << ".png";
////					canvas1->SaveAs(tmp2.str().c_str());
//					test = nullptr;
//					delete test;
//				}
////				canvas2->cd();
////				cut_performance->Draw("p");
////				canvas2->SetGridx();
////				tmp.str("");
////				tmp << "/home/lycoris-dev/Documents/cut_test/size_cut_performance_s" << sens << "_s" << sens2 <<".pdf";
////				canvas2->SaveAs(tmp.str().c_str());

////				canvas3->cd();
////				cut_entries->Draw("p");
////				canvas3->SetGridx();
////				tmp.str("");
////				tmp << "/home/lycoris-dev/Documents/cut_test/size_cut_entries_s" << sens << "_s" << sens2 <<".pdf";
////				canvas3->SaveAs(tmp.str().c_str());
//			}
//		}
//	}

	
	claus_file.close();
	
	cout << endl;
	cout << "Writing root plots to " << outRoot << endl;
	cout << endl;
	
	rFile->Write();
	gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
	rFile->Close();
	
	

	dataRead.close();
	return(0);
}

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
	
	int bucket_checking = 1;
	int  kpix_checking = 12;
	
	// cycles to skip in front:
	long int                   skip_cycles_front;
	FILE*                  f_skipped_cycles;
	string                 outtxt;
	
	string                 calState;
	uint                   lastPct;
	uint                   currPct;
	bool                   kpixFound[kpix_checking]; // variable that gives true if a kpix at the index n (0<=n<32) was found
	bool                   chanFound[kpix_checking][1024];
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



	//TTree*					cluster_tree;
	
	TH1F					*fc_response_medCM_subtracted[kpix_checking];
		
	TH1F 					*Max_SoN[kpix_checking/2];
		
		
	TH1F 					*cluster_position_y[kpix_checking/2][4];
	//TH1F 					*cluster_position_x[kpix_checking/2][4];
	TH1F 					*clusters[kpix_checking/2];
	TH1F 					*cluster_charge[kpix_checking][3];
	TH1F 					*cluster_significance[kpix_checking][3];
	TH1F 					*cluster_size[kpix_checking][3];
	
		
		
	TH2F					*cluster_correlation[kpix_checking][kpix_checking-1];
	TH1F					*cluster_offset_y[kpix_checking/2][kpix_checking/2-1];
	TH1F					*cluster_offset_x[kpix_checking/2][kpix_checking/2-1];
		
	TH1F					*noise_distribution[kpix_checking];
	TH1F					*noise_v_position[kpix_checking/2];
	// Stringstream initialization for histogram naming
	stringstream           tmp;
	stringstream           tmp_units;
	
	// Stringstream initialization for folder naming
	
	stringstream			FolderName;
	
	ofstream				claus_file;
	ofstream				noise_file;
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
	double					calib_slope[24][1024] = {1}; //ADD buckets later.
	
	double					pedestal_MedMAD[24][1024][4][2] = {0};
	int						calibration_check = 0;
	double 					noise[kpix_checking][1024];
	double 					MaximumSoN[kpix_checking/2][1840] = {0};
	//int						pedestal_check = 0;
	
	unordered_map<uint, uint> kpix2strip_left;
	unordered_map<uint, uint> kpix2strip_right;
	
	
	kpix2strip_left = kpix_left();
	kpix2strip_right = kpix_right();
	
	pixel					pixel_kpix[1024];
	pixel_mapping(pixel_kpix);
	
	unordered_map<uint, uint> sensor2layer;
	
	sensor2layer.insert(make_pair(0, 15));
	sensor2layer.insert(make_pair(1, 14));
	sensor2layer.insert(make_pair(2, 13));
	sensor2layer.insert(make_pair(3, 10));
	sensor2layer.insert(make_pair(4, 11));
	sensor2layer.insert(make_pair(5, 12));
	
	
	
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
	
	std::map<int, double> common_modes_median[kpix_checking];
	while ( dataRead.next(&event) ) // event read to check for filled channels and kpix to reduce number of empty histograms.
	{
		acqCount++;		
		int cycle_time_local[kpix_checking][8192] = {0}; //list for each kpix that is part of improved monster finder
		
		

		if (acqCount > skip_cycles_front)
		{
			acqProcessed++;
			vector<double> vec_corr_charge[kpix_checking];
			
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
				
				if ( type == KpixSample::Data )
				{
					if (kpix == 13 || kpix ==16) cout << "Found KPIX " << kpix << endl;
					kpixFound[kpix]          = true;
					chanFound[kpix][channel] = true;
					cycle_time_local[kpix][int(tstamp)]+=1;
					
					if (bucket == 0)
					{
						//cout << "Test" << endl;
						if (calibration_check == 1)
						{
							if (pedestal_MedMAD[kpix][channel][bucket][1] != 0)
							{
								double charge_value = double(value)/calib_slope[kpix][channel];
								double corrected_charge_value_median = charge_value - pedestal_MedMAD[kpix][channel][bucket][0];
							
							
							//if (calib_slope[kpix][channel] > 1 && calib_slope[kpix][channel] < 30)
							//{
								
								vec_corr_charge[kpix].push_back(corrected_charge_value_median);
								
								
								
								

								//if (CM2[kpix]/hitcounter[kpix] > 2)
								//{
								//cout << "KPiX Number = " << kpix << "        cycle number = " << acqProcessed << "             channel = " << channel << "                  charge value = " << charge_value << "                 pedestal_gauss = " << pedestal_gauss[kpix][channel][0] <<  "            common mode = " <<  CM2[kpix]/hitcounter[kpix] << endl;
								//if (channel == 373) cout << "KPiX Number = " << kpix << "       channel = " << channel << "                  pedestal median = " << pedestal_median[kpix][channel][0] << endl;
								
								//}
								
							}
						}
					}
					
					
				}
				//cout << "KPIX: " << kpix << endl;
				//cout << "Channel: " << channel << endl;
				
			}
			for (int k = 0; k < kpix_checking ; k++)
			{
				if (vec_corr_charge[k].size() != 0)
				{
					common_modes_median[k].insert(std::pair<int, double>(event.eventNumber(), median(vec_corr_charge[k])));
					//cout << "Common mode median of KPiX " << k << " is " << median(vec_corr_charge[k]) << endl;
					//}
					
					//cout << "KPiX Number = " << kpix <<  "            common mode mean calculated = " << CM2[k]/hitcounter[k] <<  "            common mode median calculated = " << median(vec_corr_charge[k]) << endl;
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
	
	
	for (sensor = 0; sensor < kpix_checking/2; sensor++) //looping through all possible kpix
	{
		//
		//cout << "DEBUG test " << kpixFound[kpix] << endl;
		if (kpixFound[kpix]) //checking if kpix exists
		{
			rFile->cd(); //producing subfolder for kpix same as above for the event subfolder structure
			FolderName.str("");
			FolderName << "Sensor_" << sensor;
			rFile->mkdir(FolderName.str().c_str());
			TDirectory *sensor_folder = rFile->GetDirectory(FolderName.str().c_str());
			sensor_folder->cd();
			
			for (int s = sensor+1; s < kpix_checking/2; ++s)
			{
				tmp.str("");
				tmp << "cluster_offset_y_sens" << sensor << "_to_sens" << s << "_b0";
				cluster_offset_y[sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
				
				tmp.str("");
				tmp << "cluster_offset_x_sens" << sensor << "_to_sens" << s << "_b0";
				cluster_offset_x[sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
				
				tmp.str("");
				tmp << "cluster_correlation_sens" << sensor << "_to_sens" << s << "_b0";
				tmp_units.str("");
				tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
				cluster_correlation[sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-0.5, 91999.5, 230,-0.5, 91999.5);
			}
			
			tmp.str("");
			tmp << "Max_SoN_sens" << sensor << "_b0";
			Max_SoN[sensor] = new TH1F(tmp.str().c_str(), "cluster position y0; strip; #Entries", 1840,-0.5, 1839.5);
			
			tmp.str("");
			tmp << "cluster_position_y0_sens" << sensor << "_b0";
			cluster_position_y[sensor][0] = new TH1F(tmp.str().c_str(), "cluster position y0; #mum; #Entries", 1840,-0.5, 91999.5);
			tmp.str("");
			tmp << "cluster_position_y1_sens" << sensor << "_b0";
			cluster_position_y[sensor][1] = new TH1F(tmp.str().c_str(), "cluster position y1; #mum; #Entries", 1840,-0.5, 91999.5);
			tmp.str("");
			tmp << "cluster_position_y2_sens" << sensor << "_b0";
			cluster_position_y[sensor][2] = new TH1F(tmp.str().c_str(), "cluster position y2; #mum; #Entries", 1840,-0.5, 91999.5);
			tmp.str("");
			tmp << "cluster_position_y3_sens" << sensor << "_b0";
			cluster_position_y[sensor][3] = new TH1F(tmp.str().c_str(), "cluster position y3; #mum; #Entries x Charge", 1840,-0.5, 91999.5);
			
			
			
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
			tmp << "clusters_sens" << sensor << "_b0";
			clusters[sensor] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; #Entries", 100,-0.5, 99.5);
			
			tmp.str("");
			tmp << "cluster_charge0_sens" << sensor << "_b0";
			cluster_charge[sensor][0] = new TH1F(tmp.str().c_str(), "cluster charge0; Charge (fC); #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_charge1_sens" << sensor << "_b0";
			cluster_charge[sensor][1] = new TH1F(tmp.str().c_str(), "cluster charge1; Charge (fC); #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_charge2_sens" << sensor << "_b0";
			cluster_charge[sensor][2] = new TH1F(tmp.str().c_str(), "cluster charge1; Charge (fC); #Entries", 200,-0.5, 49.5);
			
			
			tmp.str("");
			tmp << "cluster_significance0_sens" << sensor << "_b0";
			cluster_significance[sensor][0] = new TH1F(tmp.str().c_str(), "cluster significance0; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance1_sens" << sensor << "_b0";
			cluster_significance[sensor][1] = new TH1F(tmp.str().c_str(), "cluster significance1; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance2_sens" << sensor << "_b0";
			cluster_significance[sensor][2] = new TH1F(tmp.str().c_str(), "cluster significance1; S/N; #Entries", 200,-0.5, 49.5);
			
			tmp.str("");
			tmp << "cluster_size0_sens" << sensor << "_b0";
			cluster_size[sensor][0] = new TH1F(tmp.str().c_str(), "cluster size0; Size; #Entries", 10,-0.5, 9.5);
			tmp.str("");
			tmp << "cluster_size1_sens" << sensor << "_b0";
			cluster_size[sensor][1] = new TH1F(tmp.str().c_str(), "cluster size1; Size; #Entries", 10,-0.5, 9.5);
			tmp.str("");
			tmp << "cluster_size2_sens" << sensor << "_b0";
			cluster_size[sensor][2] = new TH1F(tmp.str().c_str(), "cluster size1; Size; #Entries", 10,-0.5, 9.5);
			
			tmp.str("");
			tmp << "noise_v_position_s" << sensor << "_b0";
			noise_v_position[sensor]  = new TH1F(tmp.str().c_str(), "Noise; #mum; Noise (fC)", 1840,-0.5, 91999.5);

			for (int k = 0; k < 2; k++) //looping through all possible kpix (left and right of each sensor)
			{
				kpix = (sensor*2)+k;
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
			}
		}
	}
	//////////////////////////////////////////
	// Data read for all events for detailed look
	//////////////////////////////////////////
	dataRead.open(argv[1]); //open file again to start from the beginning
	
	claus_file.open("claus_file_new.txt");

	int header = 1;
	

	int clstrcounter[kpix_checking][1840] = {0};
	
	
	std::vector<double> corrected_charge_vec[kpix_checking][1024];
	
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
			
				if ( type == KpixSample::Data ) // If event is of type KPiX data
				{
					if (bucket == 0)
					{
						if (pedestal_MedMAD[kpix][channel][bucket][1] != 0) //ensuring we ignore 0 MAD channels
						{
							double charge_value = double(value)/calib_slope[kpix][channel];
							double corrected_charge_value_median = charge_value - pedestal_MedMAD[kpix][channel][bucket][0];
							double charge_CM_corrected = corrected_charge_value_median - common_modes_median[kpix].at(event.eventNumber());
							corrected_charge_vec[kpix][channel].push_back(charge_CM_corrected);
						}
						//else if (pedestal_MedMAD[kpix][channel][bucket][1] == 0 && kpix == 0) cout << "1KPIX " << kpix << " Channel " << channel << endl;
					}
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
	for (int k = 0; k < kpix_checking; ++k)
	{
		for (int c = 0; c < 1024; ++c)
		{
			int sensor = k/2;
			int strip = 9999;
			if (k%2 == 0) strip = kpix2strip_left.at(c);// if left kpix
			else strip  = kpix2strip_right.at(c); // if right kpix
	    
			y = yParameter(strip, kpix);
			
			noise[k][c] = 1.4826*MAD(corrected_charge_vec[k][c]);
			//if (MAD(corrected_charge_vec[k][c]) == 0)
			//{
				////cout << "2KPIX " << k << " Channel " << c << endl;
				////for (auto const& i : corrected_charge_vec[k][c])
				////{
					////cout << "Value is " << i << endl;
				////}
			//}
			noise_distribution[k]->Fill(noise[k][c]);
			noise_v_position[sensor]->Fill(y, noise[k][c]);
		}
	}
	
	// BEGIN OF CLUSTER READ
	dataRead.open(argv[1]);
	
	
	int global_trig_counter = 0;
	
	while ( dataRead.next(&event) )
	{
		int not_empty= 0;
		
		std::vector<double> time_ext;
		std::map<int, double> cluster_Events_after_cut[kpix_checking/2];
		std::map<int, double> cluster_Noise_after_cut[kpix_checking/2];


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
	
			if (type == 2)// If event is of type external timestamp
			{
				trig_counter++;
				global_trig_counter++;
				double time = bunchClk + double(subCount * 0.125);
				time_ext.push_back(time);
				if (trig_counter == 1)
				{
					//cout << "DEBUG: " << event.eventNumber() << " ," << time << " ," << trig_counter << endl;
					tmp.str("");
					tmp << setw(7) << time  << " ," << setw(5) << global_trig_counter ;
					
				}
			}
			if ( type == KpixSample::Data ) // If event is of type KPiX data
			{
				
				if (sample->getEmpty()) cout << "Event is empty and I am analysing it" << endl; //debug output whether an empty event is still getting read
				else not_empty = 1; //checkmark that event is not empty

				//cout << tstamp << endl;
				if (bucket == 0)
				{
					if (pedestal_MedMAD[kpix][channel][bucket][1] != 0 && calib_slope[kpix][channel] >= 0.5) //ensuring we ignore 0 MAD channels and channels with bad 0 slopes
					{
						
						//// ====== Calculation of Charge values, with pedestal and common mode subtraction  =============
						double charge_value = double(value)/calib_slope[kpix][channel];
						double corrected_charge_value_median = charge_value - pedestal_MedMAD[kpix][channel][bucket][0];
						double charge_CM_corrected = corrected_charge_value_median - common_modes_median[kpix].at(event.eventNumber());
                
						//// ========= Event cut ============
                
						fc_response_medCM_subtracted[kpix]->Fill(charge_CM_corrected);
						//if ( charge_CM_corrected > 2*noise[kpix][channel] && charge_CM_corrected < 10 && strip != 9999)
						if ( charge_CM_corrected > 2*noise[kpix][channel] && charge_CM_corrected < 10 && strip != 9999 && noise_mask[sensor].at(strip) == 1)  //only events with charge higher than 2 sigma of the noise are taken and with their charge being lower than 10 fC (to cut out weird channels), in addition no noise masked channels and no disconnected channels.
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
			clustr Cluster[kpix_checking/2];  // Another Cluster class variable
			std::vector<clustr> multi_cluster[kpix_checking/2];
			for (sensor = 0; sensor < kpix_checking/2; sensor++)
			{
				
				if ( cluster_Events_after_cut[sensor].size() != 0 )
				{
					//cout << "===================" << endl;
					//cout << "Starting new PacMan" << endl;
					//cout << "===================" << endl;
					//cout << endl;
					clustr Input;
					Input.Elements = cluster_Events_after_cut[sensor];
					Input.Noise = cluster_Noise_after_cut[sensor];
					int num_of_clusters = 0;
					while (Input.Elements.size() != 0 && num_of_clusters < 9999) // Keep repeating the clustering until either there are no valid candidates left or the number of clusters is higher than X (currently 5)
					{
						PacMan NomNom;
						//cout << "Maximum Signal over Noise Strip " << Input.MaxSoN() << endl;
						//cout << "Maximum Charge Strip " << Input.MaxCharge_w_Noise() << endl;
						
						NomNom.Eater(Input, Input.MaxSoN(), 9999);
						if (num_of_clusters == 0)
						{
							cluster_position_y[sensor][0]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), weight);
							cluster_charge[sensor][0]->Fill(NomNom.getClusterCharge(), weight);
							cluster_size[sensor][0]->Fill(NomNom.getClusterElementssize(), weight);
							cluster_significance[sensor][0]->Fill(NomNom.getClusterSignificance(), weight);
							//cout << "Size: " << NomNom.getClusterElementssize() << endl;
							
						}
						Max_SoN[sensor]->Fill(Input.MaxSoN(), weight);
						MaximumSoN[sensor][Input.MaxSoN()]+=weight;
						cluster_position_y[sensor][1]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), weight);
						cluster_position_y[sensor][3]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), NomNom.getClusterSignificance()*weight);
						multi_cluster[sensor].push_back(NomNom.getCluster());
						cluster_size[sensor][1]->Fill(NomNom.getClusterElementssize(), weight);
						cluster_significance[sensor][1]->Fill(NomNom.getClusterSignificance(), weight);
						cluster_charge[sensor][1]->Fill(NomNom.getClusterCharge(), weight);
						//cout << "Significance: " << NomNom.getClusterSignificance() << endl;
						if (NomNom.getClusterSignificance() > 4 && NomNom.getClusterCharge() > 2)
						{
							cluster_position_y[sensor][2]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor), weight);
							cluster_charge[sensor][2]->Fill(NomNom.getClusterCharge(), weight);
							cluster_size[sensor][2]->Fill(NomNom.getClusterElementssize(), weight);
							cluster_significance[sensor][2]->Fill(NomNom.getClusterSignificance(), weight);
							if (header == 1)
							{
								header = 0;
								claus_file <<"Event Number, Layer, position, Significance, Size, Charge, time, #trig" << endl;
							}
							claus_file << setw(5) << event.eventNumber()  << ", " << setw(1) << sensor2layer.at(sensor)  << ", " <<  setw(7) << yParameterSensor(NomNom.getClusterCoG(), sensor)  << " ," << setw(7) << NomNom.getClusterSignificance() << " ," << setw(2) << NomNom.getClusterElementssize() << " ," << setw(7) << NomNom.getClusterCharge() << " ," << tmp.str().c_str() << endl;
						
						}
						num_of_clusters++;
						//cout << "Cluster Position is " << yParameterSensor(NomNom.getClusterCoG(), sensor) << endl;
					}
					clusters[sensor]->Fill(num_of_clusters);
				}
			
			}
			
			for (int sensor1 = 0; sensor1 < kpix_checking/2; sensor1++)
			{
				for (int sensor2 = sensor1+1; sensor2 < kpix_checking/2; sensor2++)
				{
					for (auto const& s1 : multi_cluster[sensor1])
					{
						int s1_count = 0;
						for (auto const& s2 : multi_cluster[sensor2])
						{
							double y1 = yParameterSensor(s1.CoG, sensor1);
							double y2 = yParameterSensor(s2.CoG, sensor2);
							double clstroffset_y  = y1 - y2;
							//cout << "Cluster offset is " << clstroffset << endl;
							cluster_offset_y[sensor1][sensor2]->Fill(clstroffset_y, weight);
							cluster_correlation[sensor1][sensor2]->Fill(y1,y2);
							if (10000 > fabs(clstroffset_y))
							{
								
								cluster_position_y[sensor2][2]->Fill(y2, weight);
								if (s1_count = 0) cluster_position_y[sensor1][2]->Fill(y1, weight);
							}
							s1_count++;
								
						}
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
	
	// Noise mask generation. ONLY USE WHEN NO MASK IS PUT IN!
	//noise_file.open("include/testbeam201905_noise_mask.h");
	//for (int s = 0; s < 6; ++s)
	//{
		//int badentries = 0;
		//noise_file << "unordered_map<uint, uint> noise_sensor_" << s << "()" << endl;
		//noise_file << "{" << endl;
		//noise_file << "    unordered_map<uint, uint> m1;" << endl;
		//for (int strips = 0; strips < 1840; strips++)
		//{
			//if (MaximumSoN[s][strips] > 0.03 ) //
			//{
				//cout << "Huge number of entries with no adjacent entries in strip " << strips << " with " << MaximumSoN[s][strips] << " entries in Sensor " << s << endl;
				//badentries += MaximumSoN[s][strips];
				//noise_file << "    m1.insert(make_pair(" << strips << ",0));" << endl;
			//}
			//else
			//{
				//noise_file << "    m1.insert(make_pair(" << strips << ",1));" << endl;
			//}
			////if (Max_SoN[s]->GetBinContent(strips) > 25)
			////{
				////cout << "Huge number of entries with no adjacent entries in strip in MaxSoN " << strips << " with " << Max_SoN[s]->GetBinContent(strips) << " entries in Sensor " << s << endl;
			////}
		//}
		//noise_file << "    return m1;" << endl;
		//noise_file << "}" << endl;
		//cout << "Number of bad entries in this sensor is " << badentries << endl;
	//}
	//noise_file.close();
	
	
	
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

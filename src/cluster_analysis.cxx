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

//double median(vector<double> &v)
//{
    //size_t n = v.size() / 2;
    //if (v.size()%2 == 0)
    //{
		//nth_element(v.begin(), v.begin()+n, v.end());
		//nth_element(v.begin(), v.begin()+n-1, v.end());
		//return (v[n]+v[n-1])/2;
	//}
	//else
	//{
		 //nth_element(v.begin(), v.begin()+n, v.end());
		 //return v[n];
	 //}
//}

//double MAD(vector<double> &v)
//{
	//double med = median(v);
	//vector<double> deviation;
	//for (auto const i:v)
	//{
		//deviation.push_back(fabs(i - med));
	//}
	//return median(deviation);
//}


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
	
	int 					datacounter	= 0;
	
	string                 serial;
	KpixSample::SampleType type;	



	
	TH1F				*fc_response_medCM_subtracted[kpix_checking];
	
	TH1F 				*cluster_position[kpix_checking/2][3];
	//TH1F 				*cluster_charge[kpix_checking];
	//TH1F 				*cluster_size[kpix_checking];
	
	
	//TH2F				*cluster_correlation[kpix_checking][kpix_checking-1];
	TH1F				*cluster_offset[kpix_checking/2][kpix_checking/2-1];
	
	TH1F				*noise_distribution[kpix_checking];
	
	// Stringstream initialization for histogram naming
	stringstream           tmp;
	stringstream           tmp_units;
	
	// Stringstream initialization for folder naming
	
	stringstream			FolderName;
	
	ofstream				claus_file;
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
	//int						pedestal_check = 0;
	
	unordered_map<uint, uint> kpix2strip_left;
	unordered_map<uint, uint> kpix2strip_right;
	kpix2strip_left = kpix_left();
	kpix2strip_right = kpix_right();
	
	pixel					pixel_kpix[1024];
	pixel_mapping(pixel_kpix);
	
	
	
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
	tmp << argv[1] << "_" << pedestalname << ".pedestal.tree_external.root";
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
	double weight = 1.0;//acqProcessed;
	
	
	
	//////////////////////////////////////////
	// New histogram generation within subfolder structure
	//////////////////////////////////////////
	
	int response_bins = 440;
	double response_xmin = -20.5;
	double response_xmax = 19.5;
	
	
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
				tmp << "cluster_offset_sens" << sensor << "_to_sens" << s << "_b0";
				cluster_offset[sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -5000.5, 4999.5);
			}
			
			
			tmp.str("");
			tmp << "cluster_position0_k" << kpix << "_b0";
			cluster_position[sensor][0] = new TH1F(tmp.str().c_str(), "cluster position0; #mum; #Entries", 1840,-0.5, 92000.5);
			tmp.str("");
			tmp << "cluster_position1_k" << kpix << "_b0";
			cluster_position[sensor][1] = new TH1F(tmp.str().c_str(), "cluster position1; #mum; #Entries", 1840,-0.5, 92000.5);
			tmp.str("");
			tmp << "cluster_position2_k" << kpix << "_b0";
			cluster_position[sensor][2] = new TH1F(tmp.str().c_str(), "cluster position2; #mum; #Entries", 1840,-0.5, 92000.5);
			
			
			

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
			}
		}
	}
	//////////////////////////////////////////
	// Data read for all events for detailed look
	//////////////////////////////////////////
	dataRead.open(argv[1]); //open file again to start from the beginning
	
	claus_file.open("claus_file.txt");
	
	
	int cycle_num = 0;
	
	double cluster_cut = 1.0;
	
	std::vector<double> corrected_charge_vec[kpix_checking][1024];
	while ( dataRead.next(&event) )
	{
		cycle_num++;
		int check = 0;
		if ( cycle_num > skip_cycles_front)
		{
			
			
			std::vector<double> time_ext;
			
			
			std::map<int, double> cluster_events_after_cut[kpix_checking/2];

			double channel_charge[kpix_checking][1024];

			//cout << "Beginning a new EVENT" << endl;
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
				int sensor = kpix/2;
				bunchClk = sample->getBunchCount();
				subCount = sample->getSubCount();
				
				//channel to strip assignment.
				int strip = 9999;
				if (kpix%2 == 0) // if left kpix
				{
					strip = kpix2strip_left.at(channel);
				}
				else  // if right kpix
				{
					strip  = kpix2strip_right.at(channel);
				}
				y = yParameter(strip, kpix); //Testbeam position parameter. Already adjusted for flipped sensors and kpix 1/2 position. Not for the stereo angle
		
				if (type == 2)// If event is of type external timestamp
				{
					double time = bunchClk + double(subCount * 0.125);
					time_ext.push_back(time);

				}
				
				

				
				if ( type == KpixSample::Data ) // If event is of type KPiX data
				{
					
					if (sample->getEmpty()) cout << "Event is empty and I am analysing it" << endl;
					else check = 1;

					//cout << tstamp << endl;
					if (bucket == 0)
					{
						//cout << "Test" << endl;
						if (calibration_check == 1)
						{
							if (pedestal_MedMAD[kpix][channel][bucket][1] != 0) //ensuring we ignore 0 MAD channels
							{
								
								//// ====== Calculation of Charge values, with pedestal and common mode subtraction  =============
								double charge_value = double(value)/calib_slope[kpix][channel];
								
								double corrected_charge_value_median = charge_value - pedestal_MedMAD[kpix][channel][bucket][0];
								
								double charge_CM_corrected = corrected_charge_value_median - common_modes_median[kpix].at(event.eventNumber());
								
								channel_charge[kpix][channel] = charge_CM_corrected;
								corrected_charge_vec[kpix][channel].push_back(charge_CM_corrected);
								
								//// ========= Event cut ============

								if (kpix == 0) claus_file << cycle_num << " " << channel << " " << value << " " << charge_value << " " << corrected_charge_value_median << " " << charge_CM_corrected << endl;

								fc_response_medCM_subtracted[kpix]->Fill(charge_CM_corrected, weight);
								
								if ( charge_CM_corrected > cluster_cut) 
								{
									cluster_events_after_cut[sensor].insert(std::pair<int, double>(strip, charge_CM_corrected));
								}
							}
						}
					}
				}
				
				

			}
			if (check == 1) //if event is not empty
			{
				clustr Cluster[kpix_checking/2];  // Another Cluster class variable
				std::vector<clustr> multi_cluster[kpix_checking/2];
				for (int sensor = 0; sensor < kpix_checking/2; sensor++)
				{
					if ( cluster_events_after_cut[sensor].size() != 0)
					{
						
						//cout << "===================" << endl;
						//cout << "Starting new PacMan" << endl;
						//cout << "===================" << endl;
						//cout << endl;
						clustr Input;
						Input.Elements = cluster_events_after_cut[sensor];
						int num_of_clusters = 0;
						while (Input.Elements.size() != 0 && num_of_clusters < 5) // Keep repeating the clustering until either there are no valid candidates left or the number of clusters is higher than X (currently 4)
						{
							PacMan NomNom;
							NomNom.Eater(Input, Input.MaxCharge(), 9999);
							if (num_of_clusters == 0)
							{
								cluster_position[sensor][0]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor));
								//cluster_charge[sensor]->Fill(NomNom.getClusterCharge());
								//cluster_size[sensor]->Fill(NomNom.getElementssize());
							}
							cluster_position[sensor][1]->Fill(NomNom.getClusterCoG());
							//cluster_charge[sensor][1]->Fill(NomNom.getClusterCharge()); //currently misusing the buckets
							//cluster_size[sensor][1]->Fill(NomNom.getElementssize());
							multi_cluster[sensor].push_back(NomNom.getCluster());
							num_of_clusters++;
						}
						//cout << "Number of Clusters in the event for KPIX " << KPIX << " = " << num_of_clusters << endl;
						//cout << "Cluster CoG = " << Cluster[KPIX].CoG << " : NomNom CoG = " << NomNom.getClusterCoG() << endl;
						//cout << "Cluster Charge = " << Cluster[KPIX].Charge << " : NomNom Charge = " << NomNom.getClusterCharge() << endl;
					}
					
				
				}
				
				for (int sensor1 = 0; sensor1 < kpix_checking/2; sensor1++)
				{
					for (int sensor2 = sensor1+1; sensor2 < kpix_checking/2; sensor2++)
					{
						for (auto const& s1 : multi_cluster[sensor1])
						{
							for (auto const& s2 : multi_cluster[sensor2])
							{
								double clstroffset  = yParameterSensor(s1.CoG, sensor1) - yParameterSensor(s2.CoG, sensor2);
								//cout << "Cluster offset is " << clstroffset << endl;
								cluster_offset[sensor1][sensor2]->Fill(clstroffset);
							}
						}
					}
				}
				
			}
		datacounter++;
			
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
	cout << endl;
	cout << "Writing root plots to " << outRoot << endl;
	cout << endl;
	
	rFile->Write();
	gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
	rFile->Close();
	claus_file.close();
	
	
	dataRead.close();
	return(0);
}

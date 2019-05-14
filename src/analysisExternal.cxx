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
	uint                   channel;
	uint                   bucket;
	double                  tstamp;
	uint 					subCount;
	double 					bunchClk;
	
	double 					y;
	
	string                 serial;
	KpixSample::SampleType type;	

	TH1F				*fc_response[kpix_checking][bucket_checking+1];
	
	vector<TH1F*> 		fc_response_cycle;
	vector<TH1F*> 		fc_response_cycle_median_subtracted;

	TH1F				*fc_response_median_subtracted[kpix_checking][bucket_checking+1];
	TH1F				*fc_response_channel[kpix_checking][1024];
	TH1F				*fc_response_median_subtracted_channel[kpix_checking][1024];
	
	TH1F				*fc_response_medCM_subtracted[kpix_checking][bucket_checking+1];
	TH1F				*fc_response_medCM_subtracted_channel[kpix_checking][1024];
	
	TH1F				*median_charge[kpix_checking][bucket_checking+1];
	
	
	TH1F				*common_mode_hist[kpix_checking];
	
	
	TH1F 				*hit_position[kpix_checking][bucket_checking+1][3];
	TH1F 				*cluster_position[kpix_checking][bucket_checking+1];
	TH1F 				*cluster_charge[kpix_checking][bucket_checking+1];
	TH1F 				*cluster_size[kpix_checking][bucket_checking+1];
	
	TH2F 				*position_vs_charge_corrected[kpix_checking][bucket_checking+1];
	TH2F 				*position_vs_charge_CM_corrected[kpix_checking][bucket_checking+1];
	TH2F 				*kpix_position_vs_charge_corrected[kpix_checking][bucket_checking+1];
	TH2F 				*position_vs_charge[kpix_checking][bucket_checking+1];
	
	//TH2F				*cluster_correlation[kpix_checking][kpix_checking-1][bucket_checking+1];
	TH1F				*cluster_offset[kpix_checking][kpix_checking-1][bucket_checking+1];
	
	
	
	TH1F				*slopes[kpix_checking][bucket_checking+1];
	
	TH1F				*slopes_vs_channel[kpix_checking][bucket_checking+1];
	TH1F				*slopes_vs_strip[kpix_checking][bucket_checking+1];
	TH1F				*pedestals_mean[kpix_checking][bucket_checking+1];
	TH1F				*pedestals_median[kpix_checking][bucket_checking+1];
	TH1F				*pedestals_gauss[kpix_checking][bucket_checking+1];
	
	TH2F				*noise_correlation12[kpix_checking];
	TH2F				*noise_correlation13[kpix_checking];
	TH2F				*noise_correlation23[kpix_checking];
	TH2F				*noise_correlation1k[kpix_checking];
	TH2F				*noise_correlation12_3[kpix_checking];
	
	TH1F				*noise_distribution[kpix_checking][bucket_checking];
	
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
	
	double					pedestal_median[24][1024][4] = {0};
	double					pedestal_mean[24][1024][4] = {0};
	double					pedestal_gauss[24][1024][4] = {0};
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
			//pedestal_check = 1;
			loopdir(pedestal_file, "pedestal");
			
			cout << "Number of pedestal histograms = " << pedestal_hists.size() << endl;
			for (unsigned int i = 0; i<pedestal_hists.size(); ++i)
			{
				//cout << "Current key1 = " << cal_key->GetClassName() << endl;
				
				string pedestal_name         = pedestal_hists[i]->GetName();
				
				//cout << pedestal_name << endl;
				
				size_t kpix_num_start     = pedestal_name.find("_k")+2;
				size_t channel_num_start  = pedestal_name.find("_c")+2;
				size_t kpix_num_length       = pedestal_name.length() - kpix_num_start;
				size_t channel_num_length    = pedestal_name.find("_b") - channel_num_start;
				
				size_t bucket_num_start 		= pedestal_name.find("_b")+2;
				size_t bucket_num_length      	= pedestal_name.length();
				
			    string channel_string = pedestal_name.substr(channel_num_start, channel_num_length);
			    string kpix_string = pedestal_name.substr(kpix_num_start, kpix_num_length);
			    string bucket_string = pedestal_name.substr(bucket_num_start, bucket_num_length);
			    
			    //cout << "Channel Number = " <<  channel_string << endl;
			    //cout << "KPiX Number = " << kpix_string << endl;
			    
			    int kpix_num = stoi(kpix_string);
			    int channel_num = stoi(channel_string);
			    int bucket_num = stoi(bucket_string);
			    
				
				pedestal_mean[kpix_num][channel_num][bucket_num] = pedestal_hists[i]->GetMean();
				//cout << "Pedestal of KPiX " << kpix << " and channel " << channel << " is " <<  pedestal[kpix_num][channel_num] << endl;
				
				
				//cout << "Mean " << pedestal_hists[i]->GetMean() << endl;
				if ( pedestal_mean[kpix_num][channel_num][bucket_num] != 0  )
				{
					//pedestal_hists[i]->Fit("gaus","q");
					//double test = pedestal_hists[i]->GetFunction("gaus")->GetParameter(1);
					//cout << "Parameter 0 " << pedestal_hists[i]->GetFunction("gaus")->GetParameter(0) << endl;
					//cout << "Parameter 1 " << pedestal_hists[i]->GetFunction("gaus")->GetParameter(1) << endl;
					//cout << "Parameter 2 " << pedestal_hists[i]->GetFunction("gaus")->GetParameter(2) << endl;
					//pedestal_gauss[kpix_num][channel_num][bucket_num] = pedestal_hists[i]->GetFunction("gaus")->GetParameter(1);
					
					Double_t k,q;
					q = 0.5;
					pedestal_hists[i]->GetQuantiles(1, &k, &q);
					//cout << "Median " << k << endl;
					pedestal_median[kpix_num][channel_num][bucket_num] = k;
					
					
				}
				
				
				
				//cout << "KPiX Number = " << kpix_num << "    Channel Number = " << channel_num << "     Bucket Number = " << bucket_num << endl;
				//cout << "Mean " << pedestal_mean[kpix_num][channel_num][bucket_num] << endl;
				//cout << "Median " << pedestal_median[kpix_num][channel_num][bucket_num] << endl;
				//cout << "Gauss " << pedestal_gauss[kpix_num][channel_num][bucket_num] << endl;
				//cout << endl;
				
				
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
	tmp << argv[1] << "_" << pedestalname << ".pedestal.external.root";
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
	
	
	vector<double> common_modes_median[kpix_checking];
	
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
							double charge_value = double(value)/calib_slope[kpix][channel];
							double corrected_charge_value_median = charge_value - pedestal_median[kpix][channel][0];
							
							//if (calib_slope[kpix][channel] > 1 && calib_slope[kpix][channel] < 30)
							//{
								
								vec_corr_charge[kpix].push_back(corrected_charge_value_median);
								
								
								
								
								//cout << "cycle " << setw(4) << acqProcessed << "    kpix " << setw(4) << kpix << "    bucket " << setw(4) << bucket << "    ADC value " << setw(4) << value << "    charge value " << setw(4) << charge_value << "    corrected charge value " << setw(4) <<  corrected_charge_value_median << "    pedestal value " << setw(4) <<  pedestal_median[kpix][channel][0] << endl;
								//if (CM2[kpix]/hitcounter[kpix] > 2)
								//{
								//cout << "KPiX Number = " << kpix << "        cycle number = " << acqProcessed << "             channel = " << channel << "                  charge value = " << charge_value << "                 pedestal_gauss = " << pedestal_gauss[kpix][channel][0] <<  "            common mode = " <<  CM2[kpix]/hitcounter[kpix] << endl;
								//if (channel == 373) cout << "KPiX Number = " << kpix << "       channel = " << channel << "                  pedestal median = " << pedestal_median[kpix][channel][0] << endl;
								
								//}
								
							//}
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
					common_modes_median[k].push_back(median(vec_corr_charge[k]));
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
	
	
	
	for (kpix = 0; kpix < kpix_checking; kpix++) //looping through all possible kpix
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
			
			
			tmp.str("");
			tmp << "noise_correlation12_k" << kpix << "_b0";
			noise_correlation12[kpix] = new TH2F(tmp.str().c_str(), "noise_correlation; channel_{1};  channel_{2}", 100, -10.5, 9.5, 100, -10.5, 9.5);
			tmp.str("");
			tmp << "noise_correlation13_k" << kpix << "_b0";
			noise_correlation13[kpix] = new TH2F(tmp.str().c_str(), "noise_correlation; channel_{1};  channel_{3}", 100, -10.5, 9.5, 100, -10.5, 9.5);
			tmp.str("");
			tmp << "noise_correlation23_k" << kpix << "_b0";
			noise_correlation23[kpix] = new TH2F(tmp.str().c_str(), "noise_correlation; channel_{2};  channel_{3}", 100, -10.5, 9.5, 100, -10.5, 9.5);
			
			tmp.str("");
			tmp << "noise_correlation1k_k" << kpix << "_b0";
			noise_correlation1k[kpix] = new TH2F(tmp.str().c_str(), "noise_correlation; channel_{1};  channel_{k}", 100, -10.5, 9.5, 100, -10.5, 9.5);
			
			tmp.str("");
			tmp << "noise_correlation12_3_k" << kpix << "_b0";
			noise_correlation12_3[kpix] = new TH2F(tmp.str().c_str(), "noise_correlation; channel_{1}-channel_{3};  channel_{2}-channel_{3}", 100, -10.5, 9.5, 100, -10.5, 9.5);
			
			
			tmp.str("");
			tmp << "median_charge_k" << kpix << "_total";
			median_charge[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "median_charge; Charge (fC);   #channels", 800,-50.5, 199.5);
			
			tmp.str("");
			tmp << "noise_distribution_k" << kpix << "_total";
			noise_distribution[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "noise_distribution; Charge (fC);   #channels", 100,-0.005, 4.995);
			
			for (int k = kpix+1; k < kpix_checking; ++k)
			{
				tmp.str("");
				tmp << "cluster_offset_k" << kpix << "_to_k" << k << "_total";
				cluster_offset[kpix][k][bucket_checking] = new TH1F(tmp.str().c_str(), "strip_offset; #mum; #entries", 200, -5000.5, 4999.5);
			}
			
			
			
			//tmp.str("");
			//tmp << "cluster_correlation_k" << kpix << "_total";
			//new TH2F("cluster_correlation", "strip offset Layer 1 to Layer 2; strips of module 1; strips of module 2",230, 919.5, 1819.5, 230, 919.5, 1819.5);
			
			
			if (kpix%2 == 0)
			{
				
				tmp.str("");
				tmp << "left_strip_position_vs_charge_corrected_k" << kpix << "_total";
				position_vs_charge_corrected[kpix][bucket_checking] = new TH2F(tmp.str().c_str(), "position vs corrected_charge; Strip position;  Charge (fC)", 230, -0.5, 919.5, 20,-0.5,19.5);
				
				tmp.str("");
				tmp << "left_strip_position_vs_charge_CM_corrected_k" << kpix << "_total";
				position_vs_charge_CM_corrected[kpix][bucket_checking] = new TH2F(tmp.str().c_str(), "position vs CM corrected charge; Strip position;  Charge (fC)", 230, -0.5, 919.5, 440,-20.5,19.5);
				
				tmp.str("");
				tmp << "left_strip_position_vs_charge_k" << kpix << "_total";
				position_vs_charge[kpix][bucket_checking] = new TH2F(tmp.str().c_str(), "position vs charge; Strip position;  Charge (fC)", 230,-0.5,919.5, 200,-0.5,199.5);
			
				
				tmp.str("");
				tmp << "hit_position_singlestrip_l_k" << kpix << "_total";
				hit_position[kpix][bucket_checking][0] = new TH1F(tmp.str().c_str(), "singlestrip hit position; Strip Position; #Entries", 920, -0.5, 919.5);
				tmp.str("");
				tmp << "hit_position_doublestrip_l_k" << kpix << "_total";
				hit_position[kpix][bucket_checking][1] = new TH1F(tmp.str().c_str(), "doublestrip hit position; Strip Position; #Entries", 920, -0.5, 919.5);
				tmp.str("");
				tmp << "hit_position_single+double_strip_l_k" << kpix << "_total";
				hit_position[kpix][bucket_checking][2] = new TH1F(tmp.str().c_str(), "total strip hit position; Strip Position; #Entries", 920, -0.5, 919.5);
				
				tmp.str("");
				tmp << "cluster_position_k" << kpix << "_total";
				cluster_position[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "cluster position; #mum; #Entries", 1840,-0.5, 92000.5);
				
				tmp.str("");
				tmp << "slopes_vs_strip_l_k" << kpix << "_total";
				slopes_vs_strip[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Slope [ADC/fC]; Strip Number; Slope [ADC/fC]", 920, -0.5, 919.5);
			
				
			}
			else
			{
				tmp.str("");
				tmp << "right_strip_position_vs_charge_corrected_k" << kpix << "_total";
				position_vs_charge_corrected[kpix][bucket_checking] = new TH2F(tmp.str().c_str(), "position vs corrected_charge; Strip position;  Charge (fC)", 230, 919.5, 1839.5, 20,-0.5,19.5);
				
				tmp.str("");
				tmp << "right_strip_position_vs_charge_CM_corrected_k" << kpix << "_total";
				position_vs_charge_CM_corrected[kpix][bucket_checking] = new TH2F(tmp.str().c_str(), "position vs CM corrected charge; Strip position;  Charge (fC)", 230, 919.5, 1839.5, 440,-20.5,19.5);
				
				tmp.str("");
				tmp << "right_strip_position_vs_charge_k" << kpix << "_total";
				position_vs_charge[kpix][bucket_checking] = new TH2F(tmp.str().c_str(), "position vs charge; Strip position;  Charge (fC)", 230, 919.5, 1839.5, 200,-0.5,199.5);
			
				
			
				tmp.str("");
				tmp << "hit_position_singlestrip_r_k" << kpix << "_total";
				hit_position[kpix][bucket_checking][0] = new TH1F(tmp.str().c_str(), "singlestrip hit position; Strip Position; #Entries", 920, 919.5, 1839.5);
				tmp.str("");
				tmp << "hit_position_doublestrip_r_k" << kpix << "_total";
				hit_position[kpix][bucket_checking][1] = new TH1F(tmp.str().c_str(), "doublestrip hit position; Strip Position; #Entries", 920, 919.5, 1839.5);
				tmp.str("");
				tmp << "hit_position_single+double_strip_r_k" << kpix << "_total";
				hit_position[kpix][bucket_checking][2] = new TH1F(tmp.str().c_str(), "total strip hit position; Strip Position; #Entries", 920, 919.5, 1839.5);
			
				tmp.str("");
				tmp << "cluster_position_k" << kpix << "_total";
				cluster_position[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "cluster position; #mum; #Entries", 1840,-0.5, 92000.5);
				
				tmp.str("");
				tmp << "slopes_vs_strip_r_k" << kpix << "_total";
				slopes_vs_strip[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Slope [ADC/fC]; Strip Number; Slope [ADC/fC]", 920, 919.5, 1839.5);
				
			}
			
			
			tmp.str("");
			tmp << "kpix_position_vs_charge_corrected_k" << kpix << "_total";
			kpix_position_vs_charge_corrected[kpix][bucket_checking] = new TH2F(tmp.str().c_str(), "kpix position vs corrected_charge; KPiX position;  Charge (fC)", 256,-0.5,1023.5, 20,-0.5,19.5);
			
			
			
			tmp.str("");
			tmp << "fc_response_k" << kpix << "_total";
			fc_response[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", 800,-0.5, 399.5);
			
			tmp.str("");
			tmp << "fc_response_median_subtracted_k" << kpix << "_total";
			fc_response_median_subtracted[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", response_bins, response_xmin, response_xmax);
			
			tmp.str("");
			tmp << "fc_response_median_made_CMmedian_subtracted_k" << kpix << "_total";
			fc_response_medCM_subtracted[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", response_bins, response_xmin, response_xmax);
			
			
			
			tmp.str("");
			tmp << "cluster_charge_k" << kpix << "_total";
			cluster_charge[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "cluster charge; Charge (fC); #Entries", 200,-0.5, 49.5);
			
			tmp.str("");
			tmp << "cluster_size_k" << kpix << "_total";
			cluster_size[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "cluster size; #Strips in cluster; N/N_{total}", 10,-0.5, 9.5);
			
			
			tmp.str("");
			tmp << "common_mode_hist_k" << kpix << "_b0";
			common_mode_hist[kpix] = new TH1F(tmp.str().c_str(), "common_mode; Charge (fC);   #Entries", 400,-9.5, 10.5);
			
			
			tmp.str("");
			tmp << "pedestals_mean_" << kpix << "_total";
			pedestals_mean[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "pedestals; Charge (fC); #Entries", 800,-0.5, 399.5);
			
			tmp.str("");
			tmp << "pedestals_median_" << kpix << "_total";
			pedestals_median[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "pedestals; Charge (fC); #Entries", 800,-0.5, 399.5);
			
			tmp.str("");
			tmp << "pedestals_gauss_" << kpix << "_total";
			pedestals_gauss[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "pedestals; Charge (fC); #Entries", 800,-0.5, 399.5);
			
			tmp.str("");
			tmp << "slopes_" << kpix << "_total";
			slopes[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "slopes; Charge (fC); #Entries", 250,-0.5, 50);
			
			
			tmp.str("");
			tmp << "slopes_vs_channel_k" << kpix << "_total";
			slopes_vs_channel[kpix][bucket_checking] = new TH1F(tmp.str().c_str(), "Slope [ADC/fC]; KPiX ID; Slope [ADC/fC]", 1024, -0.5, 1023.5);
			
			FolderName.str("");
			FolderName << "cycles";
			kpix_folder->mkdir(FolderName.str().c_str());
			TDirectory *cycle_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
			rFile->cd(cycle_folder->GetPath());
			for (int cycles = 0; cycles < acqProcessed; cycles++)
			{
				FolderName.str("");
				FolderName << "cycle_" << cycles;
				cycle_folder->mkdir(FolderName.str().c_str());
				TDirectory *cycles_folder = cycle_folder->GetDirectory(FolderName.str().c_str());
				cycles_folder->cd();
				
				tmp.str("");
				tmp << "fc_response_k" << kpix << "_cycle_" << cycles << "_b0";
				fc_response_cycle.push_back(new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", 800,-0.5, 399.5));
				
				
				tmp.str("");
				tmp << "fc_response_median_subtracted_k" << kpix << "_cycle_" << cycles << "_b0";
				fc_response_cycle_median_subtracted.push_back(new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", 500,-50.5, 199.5));
				
				
			}
			
			FolderName.str("");
			FolderName << "Buckets";
			kpix_folder->mkdir(FolderName.str().c_str());
			TDirectory *buckets_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
			rFile->cd(buckets_folder->GetPath());
			
			for (bucket = 0; bucket < bucket_checking; bucket++)
			{
				FolderName.str("");
				FolderName << "Bucket " << bucket;
				buckets_folder->mkdir(FolderName.str().c_str());
				TDirectory *bucket_folder = buckets_folder->GetDirectory(FolderName.str().c_str());
				rFile->cd(bucket_folder->GetPath());

				
				tmp.str("");
				tmp << "median_charge_k" << kpix << "_b" << bucket;
				median_charge[kpix][bucket] = new TH1F(tmp.str().c_str(), "median_charge; Charge (fC);   #channels", 800,-50.5, 199.5);
				
				for (int k = kpix; k < kpix_checking; ++k)
				{
				tmp.str("");
				tmp << "cluster_offset_k" << kpix << "_to_k" << k << "_b" << bucket;
				cluster_offset[kpix][k][bucket] = new TH1F(tmp.str().c_str(), "strip_offset; #mum; #entries", 200, -5000.5, 4999.5);
				}
				
				
				tmp.str("");
				tmp << "noise_distribution_k" << kpix << "_b" << bucket;
				noise_distribution[kpix][bucket] = new TH1F(tmp.str().c_str(), "noise_distribution; Charge (fC);   #channels", 100,-0.005, 4.995);
				
				if (kpix%2 == 0)
				{
					tmp.str("");
					tmp << "left_strip_position_vs_charge_corrected_k" << kpix << "_b" << bucket;
					position_vs_charge_corrected[kpix][bucket] = new TH2F(tmp.str().c_str(), "position vs corrected_charge; Strip position;  Charge (fC)", 230,-0.5,919.5, 20,-0.5,19.5);
					
					
					tmp.str("");
					tmp << "left_strip_position_vs_charge_CM_corrected_k" << kpix << "_b" << bucket;
					position_vs_charge_CM_corrected[kpix][bucket] = new TH2F(tmp.str().c_str(), "position vs CM corrected charge; Strip position;  Charge (fC)", 230,-0.5,919.5, 440,-20.5,19.5);
				
					tmp.str("");
					tmp << "left_strip_position_vs_charge_k" << kpix << "_b" << bucket;
					position_vs_charge[kpix][bucket] = new TH2F(tmp.str().c_str(), "position vs charge; Strip position;  Charge (fC)", 230,-0.5,919.5, 200,-0.5,199.5);
				
					tmp.str("");
					tmp << "cluster_position_k" << kpix << "_b" << bucket;
					cluster_position[kpix][bucket] = new TH1F(tmp.str().c_str(), "cluster position; #mum; #Entries", 1840,-0.5, 92000.5);
				
					tmp.str("");
					tmp << "hit_position_singlestrip_l_k" << kpix << "_b" << bucket;
					hit_position[kpix][bucket][0] = new TH1F(tmp.str().c_str(), "singlestrip hit position; Strip Position; #Entries", 920, -0.5, 919.5);
					tmp.str("");
					tmp << "hit_position_doublestrip_l_k" << kpix << "_b" << bucket;
					hit_position[kpix][bucket][1] = new TH1F(tmp.str().c_str(), "doublestrip hit position; Strip Position; #Entries", 920, -0.5, 919.5);
					tmp.str("");
					tmp << "hit_position_single+double_strip_l_k" << kpix << "_b" << bucket;
					hit_position[kpix][bucket][2] = new TH1F(tmp.str().c_str(), "total strip hit position; Strip Position; #Entries", 920, -0.5, 919.5);
					
					tmp.str("");
					tmp << "slopes_vs_stripl_k" << kpix << "_b" << bucket;
					slopes_vs_strip[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope [ADC/fC]; Strip Number; Slope [ADC/fC]", 920, -0.5, 919.5);
					
					
				}
				
				
				
				
				else
				{
					tmp.str("");
					tmp << "right_strip_position_vs_charge_corrected_k" << kpix << "_b" << bucket;
					position_vs_charge_corrected[kpix][bucket] = new TH2F(tmp.str().c_str(), "position vs corrected_charge; Strip position;  Charge (fC)", 230,919.5,1839.5, 20,-0.5,19.5);
					
					
					tmp.str("");
					tmp << "right_strip_position_vs_charge_CM_corrected_k" << kpix << "_b" << bucket;
					position_vs_charge_CM_corrected[kpix][bucket] = new TH2F(tmp.str().c_str(), "position vs CM corrected charge; Strip position;  Charge (fC)", 230,919.5,1839.5, 440,-20.5,19.5);
					
					tmp.str("");
					tmp << "right_strip_position_vs_charge_k" << kpix << "_b" << bucket;
					position_vs_charge[kpix][bucket] = new TH2F(tmp.str().c_str(), "position vs charge; Strip position;  Charge (fC)", 230,919.5, 1839.5, 200,-0.5,199.5);
					
					
				
					tmp.str("");
					tmp << "cluster_position_k" << kpix << "_b" << bucket;
					cluster_position[kpix][bucket] = new TH1F(tmp.str().c_str(), "cluster position; #mum; #Entries", 1840,-0.5, 92000.5);
				
					tmp.str("");
					tmp << "hit_position_singlestrip_r_k" << kpix << "_b" << bucket;
					hit_position[kpix][bucket][0] = new TH1F(tmp.str().c_str(), "singlestrip hit position; Strip Position; #Entries", 920, 919.5, 1839.5);
					tmp.str("");
					tmp << "hit_position_doublestrip_r_k" << kpix << "_b" << bucket;
					hit_position[kpix][bucket][1] = new TH1F(tmp.str().c_str(), "doublestrip hit position; Strip Position; #Entries", 920, 919.5, 1839.5);
					tmp.str("");
					tmp << "hit_position_single+double_strip_r_k" << kpix << "_b" << bucket;
					hit_position[kpix][bucket][2] = new TH1F(tmp.str().c_str(), "total strip hit position; Strip Position; #Entries", 920, 919.5, 1839.5);
					
					
					tmp.str("");
					tmp << "slopes_vs_stripr_k" << kpix << "_b" << bucket;
					slopes_vs_strip[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope [ADC/fC]; Strip Number; Slope [ADC/fC]", 920, 919.5, 1839.5);
				}
				
				tmp.str("");
				tmp << "kpix_position_vs_charge_corrected_k" << kpix << "_b" << bucket;
				kpix_position_vs_charge_corrected[kpix][bucket] = new TH2F(tmp.str().c_str(), "kpix position vs corrected_charge; KPiX position;  Charge (fC)", 256,-0.5,1023.5, 20,-0.5,19.5);
				
				
				
				tmp.str("");
				tmp << "fc_response_k" << kpix << "_b" << bucket;
				fc_response[kpix][bucket] = new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", 800,-0.5, 399.5);
				
				
				
				tmp.str("");
				tmp << "cluster_charge_k" << kpix << "_b" << bucket;
				cluster_charge[kpix][bucket] = new TH1F(tmp.str().c_str(), "cluster charge; Charge (fC); #Entries", 200,-0.5, 49.5);
				
				tmp.str("");
				tmp << "cluster_size_k" << kpix << "_b" << bucket;
				cluster_size[kpix][bucket] = new TH1F(tmp.str().c_str(), "cluster size; #Strips in cluster; N/N_{total}", 10,-0.5, 9.5);
				
				
				tmp.str("");
				tmp << "fc_response_median_made_CMmedian_subtracted_k" << kpix << "_b" << bucket;
				fc_response_medCM_subtracted[kpix][bucket] = new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", response_bins, response_xmin, response_xmax);
				
				
				tmp.str("");
				tmp << "fc_response_median_subtracted_k" << kpix << "_b" << bucket;
				fc_response_median_subtracted[kpix][bucket] = new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", response_bins, response_xmin, response_xmax);
				
				


				
				
				
				
				
				
				
				tmp.str("");
				tmp << "pedestals_mean_" << kpix << "_b" << bucket;
				pedestals_mean[kpix][bucket] = new TH1F(tmp.str().c_str(), "pedestals; Charge (fC); #Entries", 800,-0.5, 399.5);
				
				tmp.str("");
				tmp << "pedestals_median_" << kpix << "_b" << bucket;
				pedestals_median[kpix][bucket] = new TH1F(tmp.str().c_str(), "pedestals; Charge (fC); #Entries", 800,-0.5, 399.5);
				
				tmp.str("");
				tmp << "pedestals_gauss_" << kpix << "_b" << bucket;
				pedestals_gauss[kpix][bucket] = new TH1F(tmp.str().c_str(), "pedestals; Charge (fC); #Entries", 800,-0.5, 399.5);
				
				tmp.str("");
				tmp << "slopes_" << kpix << "_b" << bucket;
				slopes[kpix][bucket] = new TH1F(tmp.str().c_str(), "slopes; Charge (fC); #Entries", 250,-0.5, 50);
				
				
				
				
				tmp.str("");
				tmp << "slopes_vs_channel_k" << kpix << "_b" << bucket;
				slopes_vs_channel[kpix][bucket] = new TH1F(tmp.str().c_str(), "Slope [ADC/fC]; KPiX ID; Slope [ADC/fC]", 1024, -0.5, 1023.5);
			}
			
			FolderName.str("");
			FolderName << "Strips_and_Channels";
			kpix_folder->mkdir(FolderName.str().c_str());
			TDirectory *channels_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
			rFile->cd(channels_folder->GetPath());
			
			for (channel = 0; channel < 1024; channel++) //looping through all possible kpix
			{
				if (chanFound[kpix]) //checking if kpix exists
				{
					int strip = 9999;
					if (kpix%2 == 0) strip = kpix2strip_left.at(channel);
					else  strip = kpix2strip_right.at(channel);
					
					FolderName.str("");
					FolderName << "strip_" << strip << "_channel_" << channel;
					channels_folder->mkdir(FolderName.str().c_str());
					TDirectory *channel_folder = channels_folder->GetDirectory(FolderName.str().c_str());
					rFile->cd(channel_folder->GetPath());
					
					
					tmp.str("");
					tmp << "fc_response_c" << channel << setw(4) << "_k" << kpix << "_b0";
					fc_response_channel[kpix][channel] = new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", 800,-0.5, 399.5);
					
					
					tmp.str("");
					tmp << "fc_response_median_subtracted_c" << channel << setw(4) << "_k" << kpix << "_b0";
					fc_response_median_subtracted_channel[kpix][channel] = new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", 500,-50.5, 199.5);
					
					tmp.str("");
					tmp << "fc_response_median_made_CMmedian_subtracted_k" << kpix << "_c" << channel << "_b0";
					fc_response_medCM_subtracted_channel[kpix][channel] = new TH1F(tmp.str().c_str(), "fc_response; Charge (fC); #Entries", response_bins, response_xmin, response_xmax);
					
				}
			}
		}
	}
	//////////////////////////////////////////
	// Data read for all events for detailed look
	//////////////////////////////////////////
	dataRead.open(argv[1]); //open file again to start from the beginning
	
	claus_file.open("claus_file.txt");
	
	
	int cycle_num = 0;
	int datacounter = 0;
	double sstrip_cut = 1.4;
	double dstrip_cut = 1.0;
	//double cluster_cut[kpix_checking];
	//cluster_cut[17]  = 0.5;
	//cluster_cut[19]= 1.4;
	
	double cluster_cut = 1.4;
	
	int ssignal = 0;
	int dsignal = 0;
	
	int sbkgrnd = 0;
	int dbkgrnd = 0;
	std::vector<double> corrected_charge_vec[kpix_checking][1024][bucket_checking+1];
	while ( dataRead.next(&event) )
	{
		cycle_num++;
		
		if ( cycle_num > skip_cycles_front)
		{
			
			
			std::vector<double> time_ext;
			
			
			std::vector<pair<int, double>> singlestrip_events_after_cut[kpix_checking][bucket_checking];
			std::vector<pair<int, double>> doublestrip_events_after_cut[kpix_checking][bucket_checking];
			std::map<int, double> cluster_events_after_cut[kpix_checking][bucket_checking];

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
				
				bunchClk = sample->getBunchCount();
				subCount = sample->getSubCount();
				
				tmp.str("");
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
						
					
					//cout << tstamp << endl;
					if (bucket == 0)
					{
						//cout << "Test" << endl;
						if (calibration_check == 1)
						{
							//if (calib_slope[kpix][channel] > 3) //ensuring we ignore 0 slope channels
							//{
								
								//// ====== Calculation of Charge values, with pedestal and common mode subtraction  =============
								double charge_value = double(value)/calib_slope[kpix][channel];
								
								double corrected_charge_value_median = charge_value - pedestal_median[kpix][channel][bucket];
								
								
								double charge_CM_corrected = corrected_charge_value_median - common_modes_median[kpix][datacounter];
								
								channel_charge[kpix][channel] = charge_CM_corrected;
								corrected_charge_vec[kpix][channel][bucket].push_back(charge_CM_corrected);
								
								//// ========= Event cut ============

								if (kpix == 0) claus_file << cycle_num << " " << channel << " " << value << " " << charge_value << " " << corrected_charge_value_median << " " << charge_CM_corrected << endl;
								
								

								


								fc_response[kpix][bucket_checking]->Fill(charge_value, weight);
								
								
								fc_response_median_subtracted[kpix][bucket_checking]->Fill(corrected_charge_value_median, weight);
								fc_response_medCM_subtracted[kpix][bucket_checking]->Fill(charge_CM_corrected, weight);
								

								position_vs_charge[kpix][bucket_checking]->Fill(strip, charge_value, 1.0);
								position_vs_charge_corrected[kpix][bucket_checking]->Fill(strip, corrected_charge_value_median, 1.0);
								position_vs_charge_CM_corrected[kpix][bucket_checking]->Fill(strip, charge_CM_corrected, 1.0);
								if ( charge_CM_corrected > cluster_cut) 
								{
									cluster_events_after_cut[kpix][bucket].insert(std::pair<int, double>(strip, charge_CM_corrected));
								}
								if ( charge_CM_corrected > sstrip_cut )
								{
									singlestrip_events_after_cut[kpix][bucket].push_back(make_pair(strip, charge_CM_corrected));
								}
								if (charge_CM_corrected > dstrip_cut)
								{
									doublestrip_events_after_cut[kpix][bucket].push_back(make_pair(strip, charge_CM_corrected));
								}
								position_vs_charge[kpix][bucket]->Fill(strip, charge_value, 1.0);
								position_vs_charge_corrected[kpix][bucket]->Fill(strip, corrected_charge_value_median, 1.0);
								position_vs_charge_CM_corrected[kpix][bucket]->Fill(strip, charge_CM_corrected, 1.0);
								

								kpix_position_vs_charge_corrected[kpix][bucket_checking]->Fill(channel, corrected_charge_value_median, 1.0);
								kpix_position_vs_charge_corrected[kpix][bucket]->Fill(channel, corrected_charge_value_median, 1.0);
								
								
								//bucket separated
								
								fc_response[kpix][bucket]->Fill(charge_value, weight);
								fc_response_median_subtracted[kpix][bucket]->Fill(corrected_charge_value_median, weight);
								fc_response_medCM_subtracted[kpix][bucket]->Fill(charge_CM_corrected, weight);
								
								
								//channel separated
								fc_response_channel[kpix][channel]->Fill(charge_value, weight);
								fc_response_median_subtracted_channel[kpix][channel]->Fill(corrected_charge_value_median, weight);
								fc_response_medCM_subtracted_channel[kpix][channel]->Fill(charge_CM_corrected, weight);
								
							
							//}
							
							
						}
					}
				}
				
				

			}
			
			for (int k = 0; k < kpix_checking; ++k) 
			{
				common_mode_hist[k]->Fill(common_modes_median[k][datacounter], weight);
				noise_correlation12[k]->Fill(channel_charge[k][312], channel_charge[k][624]);
				noise_correlation13[k]->Fill(channel_charge[k][312], channel_charge[k][899]);
				noise_correlation23[k]->Fill(channel_charge[k][624], channel_charge[k][899]);
				noise_correlation1k[k]->Fill(channel_charge[k][312], channel_charge[k][394]);
				noise_correlation12_3[k]->Fill(channel_charge[k][312]-channel_charge[k][899], channel_charge[k][624]-channel_charge[k][899]);
				//cout << "Is it a strip? " << kpix2strip_left.at(312) << " " << kpix2strip_left.at(624) << " " << kpix2strip_left.at(899) << " " << endl;
				
				
			}
			
			
			clustr Cluster[kpix_checking];  // Another Cluster class variable
			std::vector<clustr> multi_cluster[kpix_checking];
			for (int KPIX = 0; KPIX < kpix_checking; KPIX++)
			{
				if (kpixFound[KPIX] == true)
				{
					if ( cluster_events_after_cut[KPIX][0].size() != 0)
					{
						
						//cout << "===================" << endl;
						//cout << "Starting new PacMan" << endl;
						//cout << "===================" << endl;
						//cout << endl;
						clustr Input;
						Input.Elements = cluster_events_after_cut[KPIX][0];
						int num_of_clusters = 0;
						while (Input.Elements.size() != 0 && num_of_clusters < 5) // Keep repeating the clustering until either there are no valid candidates left or the number of clusters is higher than X (currently 4)
						{
							PacMan NomNom;
							NomNom.Eater(Input, Input.MaxCharge(), 9999);
							if (num_of_clusters == 0)
							{
								cluster_position[KPIX][0]->Fill(yParameter(NomNom.getClusterCoG(), KPIX));
								cluster_charge[KPIX][0]->Fill(NomNom.getClusterCharge());
								cluster_size[KPIX][0]->Fill(NomNom.getElementssize());
								cluster_position[KPIX][bucket_checking]->Fill(yParameter(NomNom.getClusterCoG(), KPIX));
								cluster_charge[KPIX][bucket_checking]->Fill(NomNom.getClusterCharge());
								cluster_size[KPIX][bucket_checking]->Fill(NomNom.getElementssize());
								Cluster[KPIX] = NomNom.getCluster();
								Cluster[KPIX].SetParameters();
							}
							cluster_position[KPIX][1]->Fill(NomNom.getClusterCoG());
							cluster_charge[KPIX][1]->Fill(NomNom.getClusterCharge()); //currently misusing the buckets
							cluster_size[KPIX][1]->Fill(NomNom.getElementssize());
							multi_cluster[KPIX].push_back(NomNom.getCluster());
							num_of_clusters++;
						}
						//cout << "Number of Clusters in the event for KPIX " << KPIX << " = " << num_of_clusters << endl;
						//cout << "Cluster CoG = " << Cluster[KPIX].CoG << " : NomNom CoG = " << NomNom.getClusterCoG() << endl;
						//cout << "Cluster Charge = " << Cluster[KPIX].Charge << " : NomNom Charge = " << NomNom.getClusterCharge() << endl;
					}
					
				}
			}
			
			for (int KPIX1 = 0; KPIX1 < kpix_checking; KPIX1++)
			{
				for (int KPIX2 = KPIX1+1; KPIX2 < kpix_checking; KPIX2++)
				{
					for (auto const& k1 : multi_cluster[KPIX1])
					{
						for (auto const& k2 : multi_cluster[KPIX2])
						{
							double clstroffset  = yParameter(k1.CoG, KPIX1) - yParameter(k2.CoG, KPIX2);
							//cout << "Cluster offset is " << clstroffset << endl;
							cluster_offset[KPIX1][KPIX2][0]->Fill(clstroffset);
							cluster_offset[KPIX1][KPIX2][bucket_checking]->Fill(clstroffset);
						}
					}
				}
			}
			
			
			
			for (int kpix = 0; kpix < kpix_checking; kpix++)
			{
				if (kpixFound[kpix] == true)
				{
					
					
					if (singlestrip_events_after_cut[kpix][0].size() == 1)  // if there is exactly 1 strip with charge above sstrip_cut, this is classified as a signal event and therefore filled into our histogram
					{
						hit_position[kpix][0][0]->Fill(singlestrip_events_after_cut[kpix][0].at(0).first, weight);
						hit_position[kpix][0][2]->Fill(singlestrip_events_after_cut[kpix][0].at(0).first, weight);
						hit_position[kpix][bucket_checking][0]->Fill(singlestrip_events_after_cut[kpix][0].at(0).first, weight);
						hit_position[kpix][bucket_checking][2]->Fill(singlestrip_events_after_cut[kpix][0].at(0).first, weight);
						
						hit_position[kpix][0][0]->Fill(singlestrip_events_after_cut[kpix][0].at(0).first, weight);
						hit_position[kpix][0][2]->Fill(singlestrip_events_after_cut[kpix][0].at(0).first, weight);
						hit_position[kpix][bucket_checking][0]->Fill(singlestrip_events_after_cut[kpix][0].at(0).first, weight);
						hit_position[kpix][bucket_checking][2]->Fill(singlestrip_events_after_cut[kpix][0].at(0).first, weight);
						
						//cout << "TEEEEST" << singlestrip_events_after_cut[kpix][0].at(0).first << endl;
						if ( singlestrip_events_after_cut[kpix][0].at(0).first > 400 && singlestrip_events_after_cut[kpix][0].at(0).first < 700) ssignal++;
						else sbkgrnd++;
						
						
					}
					
					// doublestrip signal is more complicated
					
					sort(doublestrip_events_after_cut[kpix][0].begin(), doublestrip_events_after_cut[kpix][0].end()); // sort the vector for their strip number to ensure that [strip n < strip n+1]
					int oldchannel = -9999;
					int doublestrip_count = 0;
					double doublestrip_charge = 0;
					double doublestrip_channel = 0;
					for (int v = 0; v< doublestrip_events_after_cut[kpix][0].size(); ++v) // loop over all doublestrip candidate events (charge above dstrip_cut)
					{
						//cout << "channel number = " << doublestrip_events_after_cut[kpix][0].at(v).first << endl;
						if (oldchannel+1 == doublestrip_events_after_cut[kpix][0].at(v).first) // check if the channel in the last loop is adjacent to the channel in the current loop.
						{
							doublestrip_charge = ( doublestrip_events_after_cut[kpix][0].at(v).second + doublestrip_events_after_cut[kpix][0].at(v-1).second); // add charge of this channel and the previous one
							//doublestrip_channel = double(doublestrip_events_after_cut[kpix][0].at(v).first + doublestrip_events_after_cut[kpix][0].at(v-1).first)/2; // calculate channel position (no charge weight)
							doublestrip_channel = double(doublestrip_events_after_cut[kpix][0].at(v).second * doublestrip_events_after_cut[kpix][0].at(v).first + doublestrip_events_after_cut[kpix][0].at(v-1).second * doublestrip_events_after_cut[kpix][0].at(v-1).first)/doublestrip_charge; // calculate charge weighted channel position (charge1*pos1 + charge2*pos2)/(charge1+charge2)
							doublestrip_count++;
						}
						oldchannel = doublestrip_events_after_cut[kpix][0].at(v).first;
					}
					//cout << "DStrip_event number" <<  doublestrip_count << endl;
					if (doublestrip_count == 1)
					{
						hit_position[kpix][0][1]->Fill(doublestrip_channel, weight);
						hit_position[kpix][0][2]->Fill(doublestrip_channel, weight);
						hit_position[kpix][bucket_checking][1]->Fill(doublestrip_channel, weight);
						hit_position[kpix][bucket_checking][2]->Fill(doublestrip_channel, weight);
						
						if ( doublestrip_channel > 400 && doublestrip_channel < 700) dsignal++;
						else dbkgrnd++;
						
						
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
	
	for (int kpix = 0 ; kpix < kpix_checking; kpix++)
	{
		if (kpixFound[kpix] == true)
		{
			for (int channel = 0; channel < 1024 ; channel++)
			{
				if (chanFound[kpix][channel] == true)
				{
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
					
					double mean = fc_response_median_subtracted_channel[kpix][channel]->GetMean();
					//fc_response_medCM_subtracted_channel[kpix][channel]->Fit("gaus","Rq", "", mean-0.8, mean+0.8);
					for (int bucket = 0; bucket < bucket_checking; bucket++)
					{
						pedestals_mean[kpix][bucket]->Fill(pedestal_mean[kpix][channel][bucket]);
						pedestals_median[kpix][bucket]->Fill(pedestal_median[kpix][channel][bucket]);
						pedestals_gauss[kpix][bucket]->Fill(pedestal_gauss[kpix][channel][bucket]);
						
						slopes[kpix][bucket]->Fill(calib_slope[kpix][channel]);						
						slopes_vs_channel[kpix][bucket]->SetBinContent( channel+1, calib_slope[kpix][channel]);
						
						slopes[kpix][bucket_checking]->Fill(calib_slope[kpix][channel]);
						if (kpix%2 == 0) slopes_vs_strip[kpix][bucket_checking]->SetBinContent( strip+1, calib_slope[kpix][channel]);
						else slopes_vs_strip[kpix][bucket_checking]->SetBinContent( strip-919, calib_slope[kpix][channel]);
						slopes_vs_channel[kpix][bucket_checking]->SetBinContent( channel+1, calib_slope[kpix][channel]);
						
						double mean = fc_response_median_subtracted_channel[kpix][channel]->GetMean();
						Double_t median,q;
						q = 0.5;
						fc_response_median_subtracted_channel[kpix][channel]->GetQuantiles(1, &median, &q);
						median_charge[kpix][bucket]->Fill(median);
						if (corrected_charge_vec[kpix][channel][bucket].size() != 0)
						{
							noise_distribution[kpix][bucket]->Fill(MAD(corrected_charge_vec[kpix][channel][bucket]));
							noise_distribution[kpix][bucket_checking]->Fill(MAD(corrected_charge_vec[kpix][channel][bucket]));
						}
					}
				
				}
			}
			for (int bucket = 0; bucket < bucket_checking; bucket++)
			{
				double mean = fc_response_medCM_subtracted[kpix][bucket]->GetMean();
				double RMS = fc_response_medCM_subtracted[kpix][bucket]->GetRMS();
				
				fc_response_medCM_subtracted[kpix][bucket]->Fit("gaus","Rq", "", mean-0.8, mean+0.8);
				
				
				mean = fc_response_median_subtracted[kpix][bucket]->GetMean();
				RMS = fc_response_median_subtracted[kpix][bucket]->GetRMS();
				fc_response_median_subtracted[kpix][bucket]->Fit("gaus","Rq", "", mean-0.8, mean+0.8);
				cluster_charge[kpix][bucket]->Fit("landau","RqW", "", -0.14, 17);
				if (bucket < 3) cluster_size[kpix][bucket]->Scale(1/cluster_size[kpix][bucket]->Integral(), "width");
			}
		}
	}
	//cout <<  endl << "Full coincidence of sensors with external trigger: " << full_coincidence_channel_entries->GetEntries() << endl;
	//cout << "Three coincidence of sensors: " << three_coincidence << endl;
	//cout << "Two coincidence of sensors: " << two_coincidence << endl;
	
	//cout << endl <<  "An event is currently classified as a monster if the amount of triggers one acquisition clock within a cycle is above " << monster_finder_limit << endl;
	//cout << "_______________________________________________________________" << endl;

	
	cout << "Singlestrip Signal/Background = " << ssignal << "/" << sbkgrnd << endl;
	cout << "Singlestrip Signal/Background = " << dsignal << "/" << dbkgrnd << endl;
	
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

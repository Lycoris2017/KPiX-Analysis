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

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_RESET   "\x1b[0m"

//////////////////////////////////////////
// Global Variables
//////////////////////////////////////////
//////////////////////////////////////////
// Functions
//////////////////////////////////////////

struct tree_cluster_input
{
    double CoG ;
    double Charge;
    double Significance2; //Adding noises in quadrature and later dividing chargesum by this noise sum
    double Sigma;
    int Size;
};

//void loopdir(TDirectory* dir, string histname)
//{
	
//	TDirectory *dirsav = gDirectory;
//	TIter keys_iter(dir->GetListOfKeys());
//	TKey* key;
	
//	while ((key = (TKey*)keys_iter()))
//	{
//		if (key->IsFolder())
//		{
//			dir->cd(key->GetName());
//			TDirectory *subdir = gDirectory;
//			//subfolder->cd();
//			loopdir(subdir, histname);
//			dirsav->cd();
//			continue;
//		}
//		else
//		{
//			string keyname = key->GetName();
//			string keytype = key->GetClassName();
			
			
//			if (int(histname.find("pedestal") != -1))
//			{
//				int found2 = keyname.find("_b0");
//				if (int(keytype.find("TH1") != -1))
//				{
//					if (int( keyname.find("hist_fc_") != -1) && found2 != -1)
//					{
//						//cout << "Pedestal name " << keyname << endl;
//						TH1F *pedestal_hist = (TH1F*)key->ReadObj();
//						pedestal_hist->SetName(key->GetName());
//						pedestal_hists.push_back(pedestal_hist);
//					}
//				}
//			}
//			else if ( int(histname.find("calib") != -1))
//			{
//				if (int(keytype.find("TH1") == -1) && int(keyname.find("b0") != -1))
//				{
					
//					if (int(keyname.find("calib") != -1) && int(keyname.find("DAC") == -1))
//					{
//						//cout << "Calibration name " << keyname << endl;
//						TGraphErrors *calib_graph = (TGraphErrors*)key->ReadObj();
//						calib_graph->SetName(key->GetName());
//						calib_graphs.push_back(calib_graph);
//					}
//				}
			
//			}
//		}
//	}
//}




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
    TH1F                    *ext_trigs;
    TH1F                    *ext_trigs_cycle;
    TH1F                    *trig_diff_hist;
    TH1F *entries_v_channel[n_kpix][n_buckets];
    TH2F *entries_mapped[n_kpix][n_buckets];
    TH1F *entries_v_strip[n_kpix/2][n_buckets];
    TH1F *entries_v_position[n_kpix/2][n_buckets];
    TH1F *entries_v_time[n_kpix][n_buckets];

    TH1F *entries_v_channel_timed[n_kpix][n_buckets];
    TH1F *entries_v_strip_timed[n_kpix/2][n_buckets];
    TH1F *entries_v_position_timed[n_kpix/2][n_buckets];
    TH1F *entries_v_time_timed[n_kpix][n_buckets];
	
	// Stringstream initialization for histogram naming
	stringstream           tmp;
	stringstream           tmp_units;
	
	// Stringstream initialization for folder naming
	
	stringstream			FolderName;
	
	ofstream				claus_file;
    ofstream				CM_file;
	ofstream				noise_file;
	ofstream               xml;
	ofstream               csv;
	uint                   acqCount = 0; // acquisitionCount
	uint                   acqProcessed;
	string                 outRoot;
	TFile					*rFile;


	stringstream           crossString;
	stringstream           crossStringCsv;
	XmlVariables           config;
	ofstream               debug;
	
	// Calibration slope, is filled when 
    double					calib_slope[n_kpix][n_channels][n_buckets] = {1}; //ADD buckets later.
    std::bitset<18> index;
    std::unordered_map<std::bitset<18>, std::pair<double, double>> calibs;
    double                  pearsson_cut = 0.8;
	
	double					pedestal_MedMAD[n_kpix][n_channels][n_buckets][2] = {0};
	int						calibration_check = 0;
	double 					noise[n_kpix][n_channels] = {0};
	double 					MaximumSoN[n_kpix/2][1840] = {0};
	//int						pedestal_check = 0;
	
    int maxAcquisitions = 25000;

	unordered_map<uint, uint> kpix2strip_left;
	unordered_map<uint, uint> kpix2strip_right;
	
	
	kpix2strip_left = kpix_left();
	kpix2strip_right = kpix_right();
	
	pixel	pixel_kpix[n_channels];
	pixel_mapping(pixel_kpix);
	
	unordered_map<uint, uint> sensor2layer;
	
	sensor2layer.insert(make_pair(0, 10));
	sensor2layer.insert(make_pair(1, 11));
	sensor2layer.insert(make_pair(2, 12));
	sensor2layer.insert(make_pair(3, 15));
	sensor2layer.insert(make_pair(4, 14));
	sensor2layer.insert(make_pair(5, 13));
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
    if (argc < 3) {
	cout << "Usage: ./analysis data_file [skip_cycles_front (int)]||[read calibration input file (char)] \n";
	return(1);
	}
	//cout << "DEBUG: 1" << endl;
	char* end;
    if (argc == 3){
        cout << " -- Reading " << argv[2] << " as calibration input file." << endl;
        skip_cycles_front = 0;
        TFile *calibration_file = TFile::Open(argv[2]);
        calibration_check = 1;

//			cout << "Number of calibration slopes = " << calib_graphs.size() << endl;
        TTree *calib_tree = (TTree*)calibration_file->Get("calibration_tree");
        uint kpix_calib, channel_calib, bucket_calib, range_calib;
        double slope_calib, slope_error_calib, pearsson_calib;
        calib_tree->SetBranchAddress("kpix", &kpix_calib);
        calib_tree->SetBranchAddress("channel", &channel_calib);
        calib_tree->SetBranchAddress("bucket", &bucket_calib);
        calib_tree->SetBranchAddress("range", &range_calib);
        calib_tree->SetBranchAddress("calib_slope", &slope_calib);
        calib_tree->SetBranchAddress("calib_error", &slope_error_calib);
        calib_tree->SetBranchAddress("calib_pearsson", &pearsson_calib);
        long int nEnTrees = calib_tree->GetEntries();
        for (long int i = 0; i < nEnTrees; ++i){
            calib_tree->GetEntry(i);
            index = keybit(kpix_calib, channel_calib, bucket_calib);
            std::pair<double,double> calib_entry = std::make_pair(slope_calib, pearsson_calib);
            calibs.emplace(index,calib_entry);
            calib_slope[kpix_calib][channel_calib][bucket_calib] = slope_calib;
//                calib_pearsson[kpix_calib][channel_calib][bucket_calib] = pearsson_calib;
        }
    }
	
	
    cout << "DEBUG: 2" << endl;
	//////////////////////////////////////////
	// Open Data file
	//////////////////////////////////////////
	
	
	
	if ( ! dataRead.open(argv[1])  ) {
		cout << "Error opening data file " << argv[1] << endl;
		return(1);
	}
    CM_file.open("test.txt");
	// Create output names
	

    outRoot = argv[1];
    outRoot =  outRoot + ".analysis.root";

	tmp.str("");
	tmp << argv[1] << ".GBL_input.txt" ;
	cout << "Write to GBL input file : " << tmp.str() << endl;
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
    printf(ANSI_COLOR_YELLOW "(Warning) Maximum number of cycles that are being considered is set to %i" ANSI_COLOR_RESET "\n", maxAcquisitions);
	cout << "\rReading File: 0 %" << flush;  // Printout of progress bar
	//goodTimes       	= 0;
	
	// Open root file
    rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
	rFile->cd(); // move into root folder base


    cout <<" DEEEEEBUG 2.5" << endl;
    while ( dataRead.next(&event)  &&  event.eventNumber() <= maxAcquisitions) // event read to check for filled channels and kpix to reduce number of empty histograms.
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
                bunchClk = sample->getBunchCount();
                subCount = sample->getSubCount();
                index = keybit(kpix, channel, bucket);



                //cout << type <<endl;
                //cout << "DEBUG 2" << endl;
                if ( type == KpixSample::Data )
                {
                    //cout << kpix << endl;
                    kpixFound[kpix]          = true;
                    channelFound[kpix][channel] = true;
                }
            }
        }
    }


	double weight = 1.0/acqProcessed;
	;//acqProcessed;
	
	
	

    //////////////////////////////////////////
    // New histogram generation within subfolder structure
    //////////////////////////////////////////

    int response_bins = 220;
    double response_xmin = -20.5;
    double response_xmax = 19.5;


    tmp.str("");
    tmp << "external_trigger_time";
    ext_trigs = new TH1F(tmp.str().c_str(), "ext trigger time; T (BCC); Nr. of Entries", 8192, 0, 8192);

    tmp.str("");
    tmp << "external_triggers_per_cycle";
    ext_trigs_cycle = new TH1F(tmp.str().c_str(), "ext triggers; Nr. of Trigger per cycle; Nr. of Entries", 50, -0.5, 49.5);


    tmp.str("");
    tmp << "trig_diff";
    trig_diff_hist = new TH1F(tmp.str().c_str(), "trig_diff; #Delta T (BCC); Nr. of Entries", 2001, -1000, 1000);

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

            for (bucket = 0; bucket < n_buckets; bucket++){
                tmp.str("");
                tmp << "entries_v_strip_s" << sensor << "_b" << bucket;
                cout << kpix << " " << bucket << endl;
                entries_v_strip[sensor][bucket] = new TH1F(tmp.str().c_str(), "entries_v_strip; Strip; Nr. of Entries", 1840, 0, 1839);

                tmp.str("");
                tmp << "entries_v_position_s" << sensor << "_b" << bucket;
                cout << kpix << " " << bucket << endl;
                entries_v_position[sensor][bucket] = new TH1F(tmp.str().c_str(), "entries_v_position; Position (um); Nr. of Entries", 1840, -46000, 46000);

                tmp.str("");
                tmp << "entries_v_strip_timed_s" << sensor << "_b" << bucket;
                cout << kpix << " " << bucket << endl;
                entries_v_strip_timed[sensor][bucket] = new TH1F(tmp.str().c_str(), "entries_v_strip_timed; channel; Nr. of Entries", 1840, 0, 1839);

                tmp.str("");
                tmp << "entries_v_position_timed_s" << sensor << "_b" << bucket;
                cout << kpix << " " << bucket << endl;
                entries_v_position_timed[sensor][bucket] = new TH1F(tmp.str().c_str(), "entries_v_position_timed; Position (um); Nr. of Entries", 1840, -46000, 46000);
            }
            for (int k = 0; k < 2; k++) //looping through all possible kpix (left and right of each sensor)
            {
                kpix = (sensor*2)+k;
                if (kpixFound[kpix])
                {

                    for (bucket = 0; bucket < n_buckets; bucket++){
                        tmp.str("");
                        tmp << "entries_v_channel_k" << kpix << "_b" << bucket;
                        cout << kpix << " " << bucket << endl;
                        entries_v_channel[kpix][bucket] = new TH1F(tmp.str().c_str(), "entries_v_channel; channel; Nr. of Entries", 1024, 0, 1024);



                        tmp.str("");
                        tmp << "entries_v_time_k" << kpix << "_b" << bucket;
                        cout << kpix << " " << bucket << endl;
                        entries_v_time[kpix][bucket] = new TH1F(tmp.str().c_str(), "entries_v_time; time (BCC); Nr. of Entries", 8192, 0, 8191);

                        tmp.str("");
                        tmp << "entries_v_channel_timed_k" << kpix << "_b" << bucket;
                        cout << kpix << " " << bucket << endl;
                        entries_v_channel_timed[kpix][bucket] = new TH1F(tmp.str().c_str(), "entries_v_channel_timed; channel; Nr. of Entries", 1024, 0, 1024);



                        tmp.str("");
                        tmp << "entries_v_time_timed_k" << kpix << "_b" << bucket;
                        cout << kpix << " " << bucket << endl;
                        entries_v_time_timed[kpix][bucket] = new TH1F(tmp.str().c_str(), "entries_v_time_timed; time (BCC); Nr. of Entries", 8192, 0, 8191);

                        tmp.str("");
                        tmp << "entries_mapped_k" << kpix << "_b" << bucket;
                        entries_mapped[kpix][bucket] = new TH2F(tmp.str().c_str(), "entries_mapped; kpix_x; kpix_y; Nr. of Entries", 32, -0.5, 31.5, 32, -0.5, 31.5);

                    }


                }
            }
        }

    }
	

	//TH1F* mean_noise = new TH1F("mean_noise_left", "mean_noise; noise(fC); entries", 100, -0.05, 0.95);
	

	//////////////////////////////////////////
	// Data read for all events for detailed look
	//////////////////////////////////////////
	dataRead.open(argv[1]); //open file again to start from the beginning
	
	int header = 1;
	

	int clstrcounter[n_kpix][1840] = {0};
    cout <<" DEEEEEBUG 3" << endl;
	
	std::vector<double>* corrected_charge_vec[n_kpix][n_channels] = {nullptr};
	std::vector<double>* corrected_charge_vec_time[n_kpix][n_BCC]  = {nullptr};
	
	unordered_map<uint, uint> noise_mask[6];
	noise_mask[0] = noise_sensor_0();
	noise_mask[1] = noise_sensor_1();

	noise_mask[2] = noise_sensor_2();
	noise_mask[3] = noise_sensor_3();
	noise_mask[4] = noise_sensor_4();
	noise_mask[5] = noise_sensor_5();
	
    rFile->cd();

	
	/////////////////////////////////////////////////////////
    //// BEGIN OF ANALYSIS
	dataRead.open(argv[1]);
	
	int global_trig_counter = 0;
//	std::vector<clustr> all_clusters[n_kpix/2];



    //cout << "DEBUG 3" << endl;
    while ( dataRead.next(&event) &&  event.eventNumber() <= maxAcquisitions)
	{
		int not_empty= 0;
		std::vector<double> time_ext;
        int eventSample;
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
            index = keybit(kpix, channel, bucket);
            //cout <<"DEBUG3.1 " << kpix << " " <<  channel << endl;
            int kpix_x = channel/32;
            int kpix_y = channel%32;
			//channel to strip assignment.
			int strip = 9999;
			if (kpix%2 == 0) strip = kpix2strip_left.at(channel);// if left kpix
			else strip  = kpix2strip_right.at(channel); // if right kpix

			y = yParameter(strip, kpix); //Testbeam position parameter. Already adjusted for flipped sensors and kpix 1/2 position. NOT for the stereo angle
	
			if (type == KpixSample::Timestamp){
				trig_counter++;
				global_trig_counter++;
				double time = bunchClk + double(subCount * 0.125);
				time_ext.push_back(time);
                ext_trigs->Fill(time);
				runtime = sample->getSampleRuntime64(frameruntime);
				if (frameruntime==0)
					cerr<< "Warning: frameruntime is ZEROs!"<< endl;
				
				if (trig_counter == 1){ // only take the first trigger, which suppose to relate to bucket==0
					//cout << "DEBUG: " << event.eventNumber() << " ," << time << " ," << trig_counter << endl;
					tmp.str("");
					tmp << runtime    << ","
					    << runtime*5  << " ,"
					    << global_trig_counter ;
				}
			}
			
			if ( type == KpixSample::Data ) // If event is of type KPiX data
			{
				
				if (sample->getEmpty()) cout << "Event is empty and I am analysing it" << endl; //debug output whether an empty event is still getting read
				else not_empty = 1; //checkmark that event is not empty

				//cout << tstamp << endl;
//				if (bucket == 0)
//				{
//                    if (calibs.at(index).first != 0 && calibs.at(index).second > 0.85) //ensuring we ignore 0 MAD channels and channels with bad 0 slopes
//					{
                        //cout <<"DEBUG3.2 " << kpix << " " <<  channel << " " << bucket << endl;

                double trig_diff = smallest_time_diff(time_ext, tstamp);
                trig_diff_hist->Fill(trig_diff);
//                if (fabs(trig_diff) > 10){
//                    cout << "trig diff"   << trig_diff << endl;
//                    cout << "tstamp " << tstamp << endl;
//                    for (auto const& i : time_ext){
//                        cout << i << endl;
//                    }
//                }
                        entries_v_channel[kpix][bucket]->Fill(channel,weight);
                        entries_mapped[kpix][bucket]->Fill(kpix_x, kpix_y);
                        entries_v_strip[sensor][bucket]->Fill(strip,weight);
                        entries_v_time[kpix][bucket]->Fill(tstamp);
                        entries_v_position[sensor][bucket]->Fill(yParameterSensor(strip, sensor));
                        if (trig_diff >= 2 && trig_diff < 4){
                            entries_v_channel_timed[kpix][bucket]->Fill(channel,weight);
                            entries_v_strip_timed[sensor][bucket]->Fill(strip,weight);
                            entries_v_time_timed[kpix][bucket]->Fill(tstamp);
                            entries_v_position_timed[sensor][bucket]->Fill(yParameterSensor(strip, sensor));
                        }

						
//					}
//				}
			}
        }
        ext_trigs_cycle->Fill(trig_counter);

	////   Show progress
		filePos  = dataRead.pos();
		currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
		if ( currPct != lastPct ) 
		{
            cout << "\rReading File for analysis: " << currPct << " %      " << flush;
			lastPct = currPct;
		}
		
	}



    //cout <<"DEBUG6" << endl;
	claus_file.close();
    CM_file.close();
	
	cout << endl;
	cout << "Writing root plots to " << outRoot << endl;
	cout << endl;
	
	rFile->Write();
	gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
	rFile->Close();
	
	

	dataRead.close();
	return(0);
}

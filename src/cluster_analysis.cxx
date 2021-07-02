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
#include <TProfile.h>
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
	//cout <<"DEBUG1" << endl;
	//////////////////////////////////////////
	// Class declaration and histogram initialization
	//////////////////////////////////////////
	
	cout << "CURRENTLY USING VERSION FOR NEWDAQ FOR OLD DAQ NEED TO EXCHANGE /Run WITH /20 !" << endl;
	
	
    DataRead               dataRead;  //kpix event classes used for analysis of binary data
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

    TTree*					clusterTree;
    TTree*					stripTree;
    TTree*					extData;
	
    TH1F					*QTrue[n_kpix][4];
    TProfile				*ADC_v_Time[n_kpix];
    TH1F					*SoNTrue[n_kpix];
    TH1F					*QTrue_channel[n_kpix][n_channels][4];
		
	TH1F 					*Max_SoN[n_kpix/2];
		
	TH1F 					*cluster_position_y[n_kpix/2][4];
	//TH1F 					*cluster_position_x[n_kpix/2][3];
	TH1F 					*clusters[n_kpix/2][3]; // only 2 because it makes no sense to ask for the number of clusters in the seed case cause it will always be 1.
	TH1F 					*cluster_charge[n_kpix][4];
	TH1F 					*cluster_significance[n_kpix][4];
	TH1F 					*cluster_significance2[n_kpix][4];
	TH1F 					*cluster_size[n_kpix][4];
	TH1F 					*cluster_sigma[n_kpix][4];

    TH2F                    *charge_v_position[n_kpix/2];
    TH2F                    *charge_v_strip[n_kpix/2];

    TH1F 					*trig_diff[n_kpix];
    TH1F 					*trigger_time[n_kpix];
    TH1F                    *ext_trigs;
    TH1F                    *ext_trigs2;
		
		
	TH2F					*cluster_correlation[4][n_kpix][n_kpix-1];
	TH1F					*cluster_offset_y[4][n_kpix/2][n_kpix/2-1];
	TH1F					*cluster_offset_x[4][n_kpix/2][n_kpix/2-1];

	TH1F					*cuts_entries[3][n_kpix/2][n_kpix/2-1];
	TH1F					*cuts_purity[3][n_kpix/2][n_kpix/2-1];
	TH1F					*cuts_performance[3][n_kpix/2][n_kpix/2-1];

	TH1F					*noise_distribution[n_kpix];
	TH1F					*noise_distribution_strips_only[n_kpix];
	TH1F					*noise_distribution_sensor[n_kpix/2];
    TProfile				*noise_v_position[n_kpix/2];
    TProfile				*noise_v_strip[n_kpix/2];
    TProfile				*slope_v_strip[n_kpix/2];
    TProfile				*noise_v_channel[n_kpix];
    TH2F					*noise_v_slope[n_kpix/2];

    TH1F                    *kpix_temp[n_kpix];
    TH2F                    *kpix_temp_v_runNumber[n_kpix];

    TH2F                    *noise_on_kpix[n_kpix];
//    TH1F					*noise_v_column[n_kpix];
//    TH1F					*noise_v_row[n_kpix];

	
    TProfile					*noise_v_time[n_kpix/2];
    TProfile                *noise_v_bits_high[n_kpix/2];
    TProfile                *noise_v_Graybits_high[n_kpix/2];
    TProfile                *noise_v_high_graybit[n_kpix/2];

    TH1F                    *common_mode_kpix[n_kpix];
    TProfile                *common_mode_v_event_kpix[n_kpix];
	
	TGraphErrors			*runTimeDiff;

	// Stringstream initialization for histogram naming
	stringstream           tmp;
	stringstream           tmp_units;
	
	// Stringstream initialization for folder naming
	
	stringstream			FolderName;
	
	ofstream				claus_file;
	ofstream    pixelFile;
//    ofstream				marcel_file;
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
    //std::bitset<18> index;
    //std::unordered_map<std::bitset<18>, std::pair<double, double>> calibs;
    //std::unordered_map<std::bitset<18>, std::pair<double, double>> pedestals;

    uint index;
    std::unordered_map<uint, std::pair<double, double>> calibs;
    std::unordered_map<uint, std::pair<double, double>> pedestals;

	
	int						calibration_check = 0;
	double 					noise[n_kpix][n_channels] = {0};
	double 					MaximumSoN[n_kpix/2][1840] = {0};
	//int						pedestal_check = 0;
	
	int maxAcquisitions = 100000;

	unordered_map<uint, uint> kpix2strip_left;
	unordered_map<uint, uint> kpix2strip_right;
	
	
	kpix2strip_left = kpix_left();
	kpix2strip_right = kpix_right();
	
	pixel					pixel_kpix[n_channels];
	pixel_mapping(pixel_kpix);
	
//	unordered_map<uint, uint> sensor2layer;
	
//	sensor2layer.insert(make_pair(0, 10));
//	sensor2layer.insert(make_pair(1, 11));
//	sensor2layer.insert(make_pair(2, 12));
//	sensor2layer.insert(make_pair(3, 15));
//	sensor2layer.insert(make_pair(4, 14));
//	sensor2layer.insert(make_pair(5, 13));
//	sensor2layer.insert(make_pair(6, 9999));
//	sensor2layer.insert(make_pair(7, 9999));
//	sensor2layer.insert(make_pair(8, 9999));
//	sensor2layer.insert(make_pair(9, 9999));
//	sensor2layer.insert(make_pair(10, 9999));
//	sensor2layer.insert(make_pair(11, 9999));
	
	
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
                index = keyhash(kpix_calib, channel_calib, bucket_calib);
                std::pair<double,double> calib_entry = std::make_pair(slope_calib, pearsson_calib);
                calibs.emplace(index,calib_entry);
            }
			cout << " -- Reading " << argv[3] << " as pedestal subtraction input file." << endl;
			skip_cycles_front = 0;
			TFile *pedestal_file = TFile::Open(argv[3]);

			TTree* pedestal_tree = (TTree*)pedestal_file->Get("pedestal_tree");
			
            uint kpix_pedestal, channel_pedestal, bucket_pedestal;
			double pedestal_median, pedestal_MAD;
			
			pedestal_tree->SetBranchAddress("pedestal_median", &pedestal_median);
            pedestal_tree->SetBranchAddress("kpix_num", &kpix_pedestal);
            pedestal_tree->SetBranchAddress("channel_num", &channel_pedestal);
            pedestal_tree->SetBranchAddress("bucket_num", &bucket_pedestal);
			pedestal_tree->SetBranchAddress("pedestal_MAD", &pedestal_MAD);
			
			long int nentries = pedestal_tree->GetEntries();
			for (long int i = 0; i < nentries; ++i)
			{
				pedestal_tree->GetEntry(i);
                index = keyhash(kpix_pedestal, channel_pedestal, bucket_pedestal);
                std::pair<double,double> pedestal_entry = std::make_pair(pedestal_median, pedestal_MAD);
                pedestals.emplace(index, pedestal_entry);
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
	
	double SoNCut = 2;
	string NameAdd = to_string(SoNCut);
	NameAdd = NameAdd.substr(0,3);
	
	if ( ! dataRead.open(argv[1])  ) {
		cout << "Error opening data file " << argv[1] << endl;
		return(1);
	}
    CM_file.open("test.txt");
	// Create output names
	string pedestalname = argv[3];
    string outname = argv[1];
	
	
	size_t name_start  = pedestalname.find("/Run") + 1;

    if (name_start == 0){
        printf(ANSI_COLOR_RED "DEBUG: Name start not found, seeing whether input is /Calibration" ANSI_COLOR_RESET "\n");
        name_start = pedestalname.find("/Calibration")+1;
    }
	size_t name_length = pedestalname.find(".dat") - name_start;
	
	pedestalname = pedestalname.substr(name_start, name_length);
	
	cout << "Name of output file is " <<  pedestalname << endl;
	tmp.str("");
	tmp << argv[1] << "_" << pedestalname << "_SoNG" << NameAdd << ".cluster.root";
	outRoot = tmp.str();

	tmp.str("");
	tmp << argv[1] << "_SoNG" << NameAdd << ".GBL_input.txt" ;
	cout << "Write to GBL input file : " << tmp.str() << endl;
	//claus_file.open("claus_file_new.txt");
	claus_file.open(tmp.str());
	tmp.str("");
	tmp << argv[1] << "_SoNG" << NameAdd << ".Pixel_debug.txt" ;
	cout << "Write to pixel debug file : " << tmp.str() << endl;
	pixelFile.open(tmp.str());
//    tmp.str("");
//    tmp << argv[1] << "_SoNG" << NameAdd << ".marcel.csv" ;
//    marcel_file.open(tmp.str());


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



	cout << "DEEEEEEEEEBUG 1" << endl;

	vector<double>* vec_corr_charge[n_kpix] = {nullptr}; //delete does not work if I do not initialize all vectors
	std::map<int, double> common_modes_median[n_kpix];
    while ( dataRead.next(&event)  &&  event.eventNumber() <= maxAcquisitions) // event read to check for filled channels and kpix to reduce number of empty histograms.
	{
		
		acqCount++;		
//        cout << "New acquisition" << endl;
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
                index = keyhash(kpix, channel, bucket);
                if (bucket >= n_buckets) continue;


				//cout << type <<endl;
				//cout << "DEBUG 2" << endl;
//                if (type == KpixSample::Temperature)
//                {
//                    cout << " KPiX " << kpix << " Temperature " << value << " Time " << tstamp << " eventNum " << eventNum << endl;
//                }

				if ( type == KpixSample::Data )
				{
					//cout << kpix << endl;
					kpixFound[kpix]          = true;
					channelFound[kpix][channel] = true;
//					bucketFound[kpix][channel][bucket] = true;
					//cout << "Found KPIX " << kpix << endl;
					//cout << "DEEEEEEEEEEEBUG 1.1" << endl;
					if (bucket == 0)
					{
//                        cout << "Just testing something: Index " << index << " KPiX " << kpix << " Channel " << channel << " Bucket " << bucket << endl;
//                        cout << "Previous calib slope " << calib_slope[kpix][channel][bucket] << endl;
//                        cout << "New calib slope " << calibs.at(index).first << endl;
						if (calibration_check == 1)
						{
//                            cout << "DEBUG 1.1.1 " << index << " " << kpix << " " << channel << " "<< bucket <<  endl;
//                            cout << "DEBUG 1.1.2 " << pedestals.at(index).second << endl;
//                            cout << "DEBUG 1.1.3 " << calibs.at(index).first << endl;
//                            cout << "DEBUG 1.1.4 " << calibs.at(index).second << endl;
                            //cout << "2nd Pedestal MAD " << pedestals.at(index).second << " kpix " << kpix << " channel " << channel << " bucket " << bucket << endl;
                            if (pedestals.at(index).second != 0 && calibs.at(index).first != 0  && calibs.at(index).second > 0.85)
							{
								//cout << "DEEEEEEEEEEEBUG 1.2" << endl;

								
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
                                double corrected_charge_value_median = double(value - pedestals.at(index).first)/calibs.at(index).first;
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

                    CM_file << k << " " << median(vec_corr_charge[k]) << " " <<  event.eventNumber() << endl;
					common_modes_median[k].insert(std::pair<int, double>(event.eventNumber(), median(vec_corr_charge[k])));
					delete vec_corr_charge[k];
					//cout << "Common modes median of EventNumber " << event.eventNumber()  << " kpix " << k  << " entry " << common_modes_median[k].at(event.eventNumber()) << endl;
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
	maxAcquisitions = acqCount;
	//cout << "DEBUG: 3" << endl;
	//cout << tstamp << endl;
	dataRead.close();
    double weight = 1.0;//acqProcessed;
    //acqProcessed;
	
	
	//cout << "DEEEEEEEEEEEBUG 2" << endl;

    //////////////////////////////////////////
    // New histogram generation within subfolder structure
    //////////////////////////////////////////

    int response_bins = 500;
    double response_xmin = -10.5;
    double response_xmax = 9.5;

    tmp.str("");
    tmp << "bucket_trigger_diff";
    TH1F* bucket_trigger_diff = new TH1F(tmp.str().c_str(), "bucket_trigger diff; t (BCC); Entries", 2000, 0, 2000);
    tmp.str("");
    tmp << "trigger_diff";
    TH1F* trigger_diff = new TH1F(tmp.str().c_str(), "trigger diff; t (BCC); Entries", 2000, 0, 2000);
    tmp.str("");
    tmp << "external_trigger_time";
    ext_trigs = new TH1F(tmp.str().c_str(), "ext trigger time; t (BCC); Entries", 8192, 0, 8192);
    tmp.str("");
    tmp << "external_trigger_time2";
    ext_trigs2 = new TH1F(tmp.str().c_str(), "ext trigger time2; t (BCC); Entries", 8192, 0, 8192);
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

            for (uint s = sensor+1; s < n_kpix/2; ++s)
            {
//                if (kpixFound[(s*2)] || kpixFound[(s*2+1)])
//                {
//                    tmp.str("");
//                    tmp << "cluster_offset_seed_y_sens" << sensor << "_to_sens" << s << "_b0";
//                    cluster_offset_y[0][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; Entries", 100, -10000, 10000);
//                    tmp.str("");
//                    tmp << "cluster_offset_seed_x_sens" << sensor << "_to_sens" << s << "_b0";
//                    cluster_offset_x[0][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; Entries", 200, -4000, 4000);
//                    tmp.str("");
//                    tmp << "cluster_correlation_seed_sens0" << sensor << "_to_sens" << s << "_b0";
//                    tmp_units.str("");
//                    tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (mm); Sensor " << s << " | Position (mm)";
//                    cluster_correlation[0][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46.0005, 45.9995, 230,-46.000, 46.000);

//                    tmp.str("");
//                    tmp << "cluster_offset_all_y_sens" << sensor << "_to_sens" << s << "_b0";
//                    cluster_offset_y[1][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; Entries", 100, -10000, 10000);
//                    tmp.str("");
//                    tmp << "cluster_offset_all_x_sens" << sensor << "_to_sens" << s << "_b0";
//                    cluster_offset_x[1][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; Entries", 200, -4000, 4000);
//                    tmp.str("");
//                    tmp << "cluster_correlation_all_sens" << sensor << "_to_sens" << s << "_b0";
//                    tmp_units.str("");
//                    tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
//                    cluster_correlation[1][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

//                    tmp.str("");
//                    tmp << "cluster_offset_CUTS_y_sens" << sensor << "_to_sens" << s << "_b0";
//                    cluster_offset_y[2][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; Entries", 100, -10000, 10000);
//                    tmp.str("");
//                    tmp << "cluster_offset_CUTS_x_sens" << sensor << "_to_sens" << s << "_b0";
//                    cluster_offset_x[2][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; Entries", 200, -4000, 4000);
//                    tmp.str("");
//                    tmp << "cluster_correlation_CUTS_sens" << sensor << "_to_sens" << s << "_b0";
//                    tmp_units.str("");
//                    tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
//                    cluster_correlation[2][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);
//                }

            }

            tmp.str("");
            tmp << "Max_SoN_sens" << sensor << "_b0";
            Max_SoN[sensor] = new TH1F(tmp.str().c_str(), "cluster position y0; strip; Entries", 1840,-0.5, 1839.5);

            tmp.str("");
            tmp << "cluster_position_y_seed_sens" << sensor << "_b0";
            cluster_position_y[sensor][0] = new TH1F(tmp.str().c_str(), "cluster position y0; mm; No. of Clusters", 1840,-46.000, 46.000);
            tmp.str("");
            tmp << "cluster_charge_seed_sens" << sensor << "_b0";
            cluster_charge[sensor][0] = new TH1F(tmp.str().c_str(), "cluster charge0; Charge (fC); Entries", 500,-0.5, 19.5);
            tmp.str("");
            tmp << "cluster_significance_seed_sens" << sensor << "_b0";
            cluster_significance[sensor][0] = new TH1F(tmp.str().c_str(), "cluster significance_0; S/N; Entries", 500,-0.5, 49.5);
            tmp.str("");
            tmp << "cluster_significance2_seed_sens" << sensor << "_b0";
            cluster_significance2[sensor][0] = new TH1F(tmp.str().c_str(), "cluster significance2_0; S/N; Entries", 500,-0.5, 49.5);
            tmp.str("");
            tmp << "cluster_sigma_seed_sens" << sensor << "_b0";
            cluster_sigma[sensor][0] = new TH1F(tmp.str().c_str(), "cluster sigmae0; #sigma; Entries", 200,-0.5, 4.95);
            tmp.str("");
            tmp << "cluster_size_seed_sens" << sensor << "_b0";
            cluster_size[sensor][0] = new TH1F(tmp.str().c_str(), "cluster size0; Size(Strips); Entries", 10,-0.5, 9.5);

            tmp.str("");
            tmp << "charge_v_position_s" << sensor << "_b0";
            charge_v_position[sensor] = new TH2F(tmp.str().c_str(), "charge_v_position; position y; #Charge", 1840,-46000, 46000, 200, -20.5, 19.5);

            tmp.str("");
            tmp << "charge_v_strip_s" << sensor << "_b0";
            charge_v_strip[sensor] = new TH2F(tmp.str().c_str(), "charge_v_strip; strip; #Charge", 1840,-0.5, 1839.5, 200, -20.5, 19.5);

            tmp.str("");
            tmp << "cluster_position_y_all_sens" << sensor << "_b0";
            cluster_position_y[sensor][1] = new TH1F(tmp.str().c_str(), "cluster position y1; #mum; Entries", 1840,-46000, 46000);
            tmp.str("");
            tmp << "cluster_charge_all_sens" << sensor << "_b0";
            cluster_charge[sensor][1] = new TH1F(tmp.str().c_str(), "cluster charge1; Charge (fC); Entries", 500,-0.5, 19.5);
            tmp.str("");
            tmp << "cluster_significance_all_sens" << sensor << "_b0";
            cluster_significance[sensor][1] = new TH1F(tmp.str().c_str(), "cluster significance_1; S/N; Entries", 500,-0.5, 49.5);
            tmp.str("");
            tmp << "cluster_significance2_all_sens" << sensor << "_b0";
            cluster_significance2[sensor][1] = new TH1F(tmp.str().c_str(), "cluster significance2_1; S/N; Entries", 500,-0.5, 49.5);
            tmp.str("");
            tmp << "cluster_sigma_all_sens" << sensor << "_b0";
            cluster_sigma[sensor][1] = new TH1F(tmp.str().c_str(), "cluster sigma1; #sigma; Entries", 200,-0.5, 4.95);
            tmp.str("");
            tmp << "cluster_size_all_sens" << sensor << "_b0";
            cluster_size[sensor][1] = new TH1F(tmp.str().c_str(), "cluster size1; Size(Strips); Entries", 10,-0.5, 9.5);
            tmp.str("");
            tmp << "clusters_all_sens" << sensor << "_b0";
            clusters[sensor][0] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; Entries", 400,-0.5, 399.5);



            tmp.str("");
            tmp << "cluster_position_y_CUTS_sens" << sensor << "_b0";
            cluster_position_y[sensor][2] = new TH1F(tmp.str().c_str(), "cluster position y2; mm; No. of Entries", 1840,-46.000, 46.000);
            tmp.str("");
            tmp << "cluster_charge_CUTS_sens" << sensor << "_b0";
            cluster_charge[sensor][2] = new TH1F(tmp.str().c_str(), "cluster charge2; Charge (fC); Entries", 500,-0.5, 19.5);
            tmp.str("");
            tmp << "cluster_significance_CUTS_sens" << sensor << "_b0";
            cluster_significance[sensor][2] = new TH1F(tmp.str().c_str(), "cluster significance_2; S/N; Entries", 500,-0.5, 49.5);
            tmp.str("");
            tmp << "cluster_significance2_CUTS_sens" << sensor << "_b0";
            cluster_significance2[sensor][2] = new TH1F(tmp.str().c_str(), "cluster significance2_2; S/N; Entries", 500,-0.5, 49.5);
            tmp.str("");
            tmp << "cluster_sigma_CUTS_sens" <<sensor << "_b0";
            cluster_sigma[sensor][2] = new TH1F(tmp.str().c_str(), "cluster sigma2; #sigma; Entries", 200,-0.5, 4.95);
            tmp.str("");
            tmp << "cluster_size_CUTS_sens" << sensor << "_b0";
            cluster_size[sensor][2] = new TH1F(tmp.str().c_str(), "cluster size2; Size (Strips); Entries", 10,-0.5, 9.5);
            tmp.str("");
            tmp << "clusters_CUTS_sens" << sensor << "_b0";
            clusters[sensor][1] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; Entries", 400,-0.5, 399.5);





            tmp.str("");
            tmp << "noise_v_position_s" << sensor << "_b0";
            noise_v_position[sensor]  = new TProfile(tmp.str().c_str(), "Noise; #mum; Noise (fC)", 1840,-46000, 46000);

            tmp.str("");
            tmp << "noise_v_strip_s" << sensor << "_b0";
            noise_v_strip[sensor]  = new TProfile(tmp.str().c_str(), "Noise; strip number; Noise (fC)", 1840,-0.5, 1839.5);
            
            tmp.str("");
            tmp << "slope_v_strip_s" << sensor << "_b0";
            slope_v_strip[sensor]  = new TProfile(tmp.str().c_str(), "Calib Slope; strip number; slope (ADC/fC)", 1840,-0.5, 1839.5);
            
            
            tmp.str("");
            tmp << "noise_v_slope_s" << sensor << "_b0";
            noise_v_slope[sensor]  = new TH2F(tmp.str().c_str(), "Noise_v_slope; noise (fC); slope (ADC/fC)", 50,0,10, 200, -5, 35);

            tmp.str("");
            tmp << "noise_distribution_s" << sensor << "_b0";
            noise_distribution_sensor[sensor] = new TH1F(tmp.str().c_str(), "noise_distribution; Noise(fC);   #channels", 100,-0.005, 0.995);

            tmp.str("");
            tmp << "noise_v_time_s" << sensor << "_b0";
            tmp_units.str("");
            tmp_units << "noise_v_time; Time (BCC); Noise (fC) ";
            noise_v_time[sensor]  = new TProfile(tmp.str().c_str(), tmp_units.str().c_str(), 8192, 0, 8191);

            tmp.str("");
            tmp << "noise_v_bits_high_s" << sensor << "_b0";
            tmp_units.str("");
            tmp_units << "noise_v_bits_high; High Bits; Noise (fC) ";
            noise_v_bits_high[sensor]  = new TProfile(tmp.str().c_str(), tmp_units.str().c_str(), 14, -0.5, 13.5);

            tmp.str("");
            tmp << "noise_v_Graybits_high_s" << sensor << "_b0";
            tmp_units.str("");
            tmp_units << "noise_v_Graybits_high; High Gray Bits; Noise (fC) ";
            noise_v_Graybits_high[sensor]  = new TProfile(tmp.str().c_str(), tmp_units.str().c_str(), 14, -0.5, 13.5);

            tmp.str("");
            tmp << "noise_v_high_graybit_s" << sensor << "_b0";
            tmp_units.str("");
            tmp_units << "noise_v_high_graybit; Gray Bit High; Noise (fC) ";
            noise_v_high_graybit[sensor]  = new TProfile(tmp.str().c_str(), tmp_units.str().c_str(), 14, -0.5, 13.5);

            for (uint k = 0; k < 2; k++) //looping through all possible kpix (left and right of each sensor)
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
                    tmp << "QTrue_k" << kpix << "_b0";
                    QTrue[kpix][0] = new TH1F(tmp.str().c_str(), "True charge; Charge (fC); Entries", response_bins, response_xmin, response_xmax);

                    tmp.str("");
                    tmp << "Q_ADC_k" << kpix << "_b0";
                    QTrue[kpix][1] = new TH1F(tmp.str().c_str(), "ADC charge; Charge (fC); Entries", 8192, 0, 8192);

                    tmp.str("");
                    tmp << "Q_Ped_sub_k" << kpix << "_b0";
                    QTrue[kpix][2] = new TH1F(tmp.str().c_str(), "ped sub ADC charge; Charge (fC); Entries", 8192, -4096, 4096);

                    tmp.str("");
                    tmp << "Q_fC_k" << kpix << "_b0";
                    QTrue[kpix][3] = new TH1F(tmp.str().c_str(), "ped sub fC charge; Charge (fC); Entries", response_bins, response_xmin, response_xmax);

                    tmp.str("");
                    tmp << "SoNTrue_k" << kpix << "_b0";
                    SoNTrue[kpix] = new TH1F(tmp.str().c_str(), "True S/N; S/N; Entries", 4*response_bins, 4*response_xmin, 4*response_xmax);

                    tmp.str("");
                    tmp << "ADC_v_Time" << kpix << "_b0";
                    ADC_v_Time[kpix] = new TProfile(tmp.str().c_str(), "ADC_v_Time; Time (BCC); Charge (fC)", 8192, 0, 8191);


                    tmp.str("");
                    tmp << "noise_distribution_strips_only_k" << kpix << "_b0";
                    noise_distribution_strips_only[kpix] = new TH1F(tmp.str().c_str(), "noise_distribution_strips_only; Noise(fC);   Nr. of channels", 500,-0.005, 4.995);

                    tmp.str("");
                    tmp << "noise_distribution_k" << kpix << "_b0";
                    noise_distribution[kpix] = new TH1F(tmp.str().c_str(), "noise_distribution; Noise(fC);   Nr. of channels", 500,-0.005, 4.995);
                    cout << "Found the following KPiX in the file: " <<  kpix << endl;

                    tmp.str("");
                    tmp << "noise_v_channel_k" << kpix << "_b0";
                    noise_v_channel[kpix]  = new TProfile(tmp.str().c_str(), "Noise; channel; Noise (fC)", 1024,0, 1023);

                    tmp.str("");
                    tmp << "noise_on_kpix_k" << kpix << "_b0";
                    noise_on_kpix[kpix]  = new TH2F(tmp.str().c_str(), "Noise (fc); KPiX x; KPiX y", 32,0,32 ,32,0,32);

//                    tmp.str("");
//                    tmp << "noise_v_column_k" << kpix << "_b0";
//                    noise_v_column[kpix]  = new TH1F(tmp.str().c_str(), "Noise; column; Noise (fC)", 32,0, 31);

//                    tmp.str("");
//                    tmp << "noise_v_row_k" << kpix << "_b0";
//                    noise_v_row[kpix]  = new TH1F(tmp.str().c_str(), "Noise; row; Noise (fC)", 32,0, 31);

                    tmp.str("");
                    tmp << "trig_diff_k" << kpix << "_b0";
                    trig_diff[kpix] = new TH1F(tmp.str().c_str(), "Trigger difference; #Delta t (BCC); Entries", 80, 0, 10);


                    tmp.str("");
                    tmp << "trigger_time_k" << kpix << "_b0";
                    trigger_time[kpix] = new TH1F(tmp.str().c_str(), "Trigger_time; t (BCC); Entries", 8192, 0, 8191);


                    tmp.str("");
                    tmp << "common_mode_k" << kpix << "_b0";
                    common_mode_kpix[kpix] = new TH1F(tmp.str().c_str(), "common mode; Common Mode(fC);   Entries", 240,-12, 12);

                    tmp.str("");
                    tmp << "common_mode_v_event_k" << kpix << "_b0";
                    common_mode_v_event_kpix[kpix] = new TProfile(tmp.str().c_str(), "common mode; EvetNumber;   Common Mode(fC)", 10000, 0, 10000);
					
                    tmp.str("");
                    tmp << "kpix_temp_k" << kpix;
                    kpix_temp[kpix] = new TH1F("temp_hist", "temp_hist; temperature; Entries", 256,-62, 90 );
                    tmp.str("");
                    tmp << "kpix_temp_v_runNumber_k" << kpix;
                    kpix_temp_v_runNumber[kpix] = new TH2F(tmp.str().c_str(), "kpix_temp_v_runNumber; runNumber; Temperature(C)", 10000,0,100000, 123,-24, 52);


                    for (channel = 0; channel < 1024; channel++){
						if (channelFound) {
							FolderName.str("");
                            int strip = 9999;
                            if (kpix%2 == 0) strip = kpix2strip_left.at(channel);// if left kpix
                            else strip  = kpix2strip_right.at(channel); // if right kpix
                            FolderName << "Channel_" << channel << "_Strip_" << strip;
                            kpix_folder->mkdir(FolderName.str().c_str());
                            TDirectory *channel_folder = kpix_folder->GetDirectory(FolderName.str().c_str());
                            rFile->cd(channel_folder->GetPath());
							
							tmp.str("");
                            tmp << "Q_true_k" << kpix << "_c" << channel << "_s" << strip << "_b0";
                            QTrue_channel[kpix][channel][0] = new TH1F(tmp.str().c_str(), "True charge; Charge (fC); Entries", response_bins, response_xmin, response_xmax);

                            tmp.str("");
                            tmp << "Q_ADC_k" << kpix << "_c" << channel << "_s" << strip << "_b0";
                            QTrue_channel[kpix][channel][1] = new TH1F(tmp.str().c_str(), "ADC charge; Charge (fC); Entries", 8192, 0, 8192);

                            tmp.str("");
                            tmp << "Q_Ped_sub_k" << kpix << "_c" << channel << "_s" << strip << "_b0";
                            QTrue_channel[kpix][channel][2] = new TH1F(tmp.str().c_str(), "ped sub ADC charge; Charge (fC); Entries", 8192, -4096, 4096);

                            tmp.str("");
                            tmp << "Q_fC_k" << kpix << "_c" << channel << "_s" << strip << "_b0";
                            QTrue_channel[kpix][channel][3] = new TH1F(tmp.str().c_str(), "ped sub fC charge; Charge (fC); Entries", response_bins, response_xmin, response_xmax);



						}
					}
                }
            }
        }

    }
	

	//TH1F* mean_noise = new TH1F("mean_noise_left", "mean_noise; noise(fC); entries", 100, -0.05, 0.95);
	

	//////////////////////////////////////////
	// Data read for all events for detailed look
	/////////////////////////////////////////SoNCut/
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
	
    rFile->cd();

    uint grCount = 0;
    double prevTime = 0;
    double time = 0;
    while ( dataRead.next(&event)  &&  event.eventNumber() <= maxAcquisitions) //preread to determine noise value of each channel in each KPiX.
	{
		int not_empty= 0;
        std::vector<double> time_ext;
        uint TBucket0;
        uint TBucket1;

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
            index = keyhash(kpix, channel, bucket);
            if (bucket >= n_buckets) continue;
//            cout << "DEBUG: " << kpix << " " << channel << " " << value << " " << type << endl;

            if (type == KpixSample::Temperature)
            {
                uint tempADC = value%256; // supposedly only the first 8 bit are temperature values
//                uint graytemp = time%256;
//                test = bin(temp, test);
//                cout << test << endl;
//                std::string greytest = binarytoGray(test);
//                std::string greytstamp = "";
//                greytstamp = bin(graytemp, greytstamp);
//                cout << greytest << " : " << greytstamp << endl;
                double temp=0.598*(255-tempADC)-62;
//                cout << kpix << endl;
//                cout << "test" << endl;
//                cout << tmp<< endl;

                //kpix_temp[kpix]->Fill(tmp);
                //kpix_temp_v_runNumber[kpix]->Fill(event.eventNumber(),tmp);
            }

            if (type == 2)// If event is of type external timestamp
            {
                prevTime = time;
                time = bunchClk + double(subCount * 0.125);
                time_ext.push_back(time);
                trigger_diff->Fill(time-prevTime);
                ext_trigs->Fill(time);
                ext_trigs2->Fill(tstamp);
                grCount++;
                //cout << "DEBUG: channel in timestmap = " << channel << endl;
                //cout << "DEBUG: bucket in timestmap = " << bucket << endl;
            }
			if ( type == KpixSample::Data ) // If event is of type KPiX data
			{
//				cout << "Bucket " << bucket << " Time " << tstamp << endl;
				if (bucket == 0)
				{
                    TBucket0 = tstamp;
                    if (pedestals.at(index).second != 0 && calibs.at(index).first != 0  && calibs.at(index).second > 0.85 ) //ensuring we ignore 0 MAD channels
					{
//                         cout << "DEBUG: " << kpix << " " << channel << " " << value << " " << type << endl;
                        //if (kpix == 0 && channel == 188) cout << "Just testing things: " << calibs.at(index).first << " " << calibs.at(index).second;
                        trigger_time[kpix]->Fill(tstamp);
                        double time_diff_triggers = smallest_time_diff(time_ext, tstamp);
                        double charge_ped_corrected = value - pedestals.at(index).first;
                        trig_diff[kpix]->Fill(time_diff_triggers);
                        QTrue_channel[kpix][channel][1]->Fill(value, weight);
                        QTrue_channel[kpix][channel][2]->Fill(charge_ped_corrected, weight);
                        ADC_v_Time[kpix]->Fill(tstamp, value);
						//cout << "DEBUG 1" << endl;
                        double corrected_charge_value_median = charge_ped_corrected/calibs.at(index).first;
                        QTrue_channel[kpix][channel][3]->Fill(corrected_charge_value_median, weight);
                        double true_charge = corrected_charge_value_median - common_modes_median[kpix].at(event.eventNumber());
                        QTrue_channel[kpix][channel][0]->Fill(true_charge, weight);
                        if (channel == 0){
                            common_mode_v_event_kpix[kpix]->Fill(event.eventNumber(),common_modes_median[kpix].at(event.eventNumber()));
                            common_mode_kpix[kpix]->Fill(common_modes_median[kpix].at(event.eventNumber()), weight);
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
//                        cout << "CM corrected charge is " << true_charge << " for KPiX " << kpix << " and channel " << channel << endl;
                        corrected_charge_vec[kpix][channel]->push_back(true_charge);
                        corrected_charge_vec_time[kpix][int(tstamp)]->push_back(true_charge);

					}
                    //else if (pedestals.at(index).second == 0 && kpix == 0) cout << "1KPIX " << kpix << " Channel " << channel << endl;
				}
                if (bucket == 1)
                {
                    TBucket1 = tstamp;
//                    cout <<  tstamp - TBucket0 << endl;
                }
			}
		}
		bucket_trigger_diff->Fill(TBucket1-TBucket0);
		filePos  = dataRead.pos();
		currPct = (uint)(((double)filePos / (double)fileSize) * 100.0);
		if ( currPct != lastPct ) 
		{
			cout << "\rReading File for noise determination: " << currPct << " %      " << flush;
			lastPct = currPct;
		}


	}
	dataRead.close();
//    marcel_file << "Sensor, Time (BCC), Noise (fC)" << endl;

    //cout << "DEBUG 1" << endl;
	 // END OF PREREAD
	 // BEGIN OF NOISE CALCULATION
	for (unsigned int k = 0; k < n_kpix; ++k)
	{
		if (kpixFound[k])
		{
            sensor = k/2;
			for (unsigned int c = 0; c < n_channels; ++c)
			{
				if (channelFound[k][c])
				{
                    index = keyhash(k,c,0);
					int strip = 9999;
					if (k%2 == 0) strip = kpix2strip_left.at(c);// if left kpix
					else strip  = kpix2strip_right.at(c); // if right kpix
					if (pedestals.at(index).second != 0 && calibs.at(index).first != 0  && calibs.at(index).second > 0.85 )
					{
						y = yParameterSensor(strip, sensor);
						noise[k][c] = 1.4826*MAD(corrected_charge_vec[k][c]);
						noise_distribution[k]->Fill(noise[k][c]);
//                        cout << "Noise of kpix " << k << " and channel " << c << " is " << noise[k][c] << endl;
						if (strip != 9999)
						{
							noise_distribution_strips_only[k]->Fill(noise[k][c]);
							noise_distribution_sensor[sensor]->Fill(noise[k][c]);
							noise_v_slope[sensor]->Fill(noise[k][c], calibs.at(index).first);
							noise_v_position[sensor]->Fill(y, noise[k][c]);
                            noise_v_strip[sensor]->Fill(strip, noise[k][c]);
                            slope_v_strip[sensor]->Fill(strip, calibs.at(index).first);
						}
						noise_v_channel[k]->SetBinContent(c+1, noise[k][c]);
						noise_on_kpix[k]->SetBinContent((c+1)/32, (c+1)%32, noise[k][c]);
						//noise_v_column[k]->SetBinContent((c+1)/32, noise[k][c]);
						//noise_v_row[k]->SetBinContent((c+1)%32, noise[k][c]);
						
					}
				}
				
			}
			for (unsigned int t = 0; t < n_BCC ; ++t)
			{
                double time_noise =  1.4826*MAD(corrected_charge_vec_time[k][t]);
                noise_v_time[sensor]->Fill(t, time_noise);
                string TimeBin = "";
                TimeBin = bin(t, TimeBin);
                uint bitsHigh = 0;
                for (char const &c : TimeBin) {
                    if (c == '1'){
//                        cout << "Test 1" << endl;
                        bitsHigh++;
                    }
//                    if (c == '0'){
////                        cout << "Test 0" << endl;
//                    }
                }
//                std::string TimeGray = binarytoGray(TimeBin);
//                bitsHigh =
                noise_v_bits_high[sensor]->Fill(bitsHigh, time_noise, 1);
                uint graybitsHigh = 0;
                string TimeGray = binarytoGray(TimeBin);
                uint bitnum = 0;
                for (char const &c : TimeGray){
                    if (c == '1'){
                        graybitsHigh++;
                        noise_v_high_graybit[sensor]->Fill(bitnum, time_noise, 1);

                    }
                    bitnum++;
                }
                noise_v_Graybits_high[sensor]->Fill(graybitsHigh, time_noise, 1);

//                marcel_file << setw(3) << sensor  << ","
//                            << setw(5) << t << ","
//                            << setw(5) << time_noise
//                            << endl;
			}	
		}
	}
    //cout << "DEBUG 2" << endl;
	/////////////////////////////////////////////////////////
	//// BEGIN OF CLUSTER READ
	dataRead.open(argv[1]);
	
	int global_trig_counter = 0;
//	std::vector<clustr> all_clusters[n_kpix/2];

    uint clusterID = 0;

    std::vector<double> vector_charge, vector_sigma, vector_size, vector_pos, strip_pos, strip_signi, strip_charge;
    std::vector<int> vector_sensor, vector_ID,  strip_sensor, strip_ID;

    clusterTree = new TTree("vectorCluster", "vector cluster tree");
    clusterTree->Branch("sensor", &vector_sensor);
    clusterTree->Branch("charge", &vector_charge);
    clusterTree->Branch("signi", &vector_sigma);
    clusterTree->Branch("size", &vector_size);
    clusterTree->Branch("pos", &vector_pos);
    clusterTree->Branch("ID", &vector_ID);

    stripTree = new TTree("vectorClusterStrip", "vector strip tree");
    stripTree->Branch("sensor", &strip_sensor);
    stripTree->Branch("charge", &strip_charge);
    stripTree->Branch("signi", &strip_signi);
    stripTree->Branch("pos", &strip_pos);
    stripTree->Branch("ID", &strip_ID);
    uint64_t prev_runtime = 0;

    uint diffTime, TriggerNumber;
    extData = new TTree("extData", "extData");
    extData->Branch("runtimeDiff", &diffTime);
    extData->Branch("TrigNumber", &TriggerNumber);
    extData->Branch("runtime", &runtime, "runtime/l");
    grCount = 0;
    //cout << "DEBUG 3" << endl;
    pixelFile <<"Event Number,Layer,position,Significance2,Size,Charge,runtime,runtime_ns,trigN,ClusterID" << endl;
    while ( dataRead.next(&event) &&  event.eventNumber() <= maxAcquisitions)
	{
		int not_empty= 0;
		std::vector<double> time_ext;
		std::map<int, double> cluster_Events_after_cut[n_kpix/2];
		std::map<int, double> cluster_Noise_after_cut[n_kpix/2];
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
            index = keyhash(kpix, channel, bucket);
            //cout <<"DEBUG3.1 " << kpix << " " <<  channel << endl;
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
                runtime = sample->getSampleRuntime64(frameruntime);
				//cout << "Runtime output "<< runtime << endl;
				if (prev_runtime != 0){
					diffTime = runtime - prev_runtime;
					TriggerNumber = grCount;
				//cout << "Runtime diff to prev " << runtime - prev_runtime << endl;
					extData->Fill();

				}
				//cout << "FrameRuntime output "<< frameruntime << endl;
				prev_runtime = runtime;
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
//            if (trig_counter > 4){
//                cout << "Found " << trig_counter << " triggers in the event!" << endl;
//            }
			
			if ( type == KpixSample::Data ) // If event is of type KPiX data
			{
				
				if (sample->getEmpty()) cout << "Event is empty and I am analysing it" << endl; //debug output whether an empty event is still getting read
				else not_empty = 1; //checkmark that event is not empty

				//cout << tstamp << endl;
				if (bucket == 0)
				{
                    if (pedestals.at(index).second != 0 && calibs.at(index).first != 0 && calibs.at(index).second > 0.85) //ensuring we ignore 0 MAD channels and channels with bad 0 slopes
					{
                        //cout <<"DEBUG3.2 " << kpix << " " <<  channel << endl;
						//// ====== Calculation of Charge values, with pedestal and common mode subtraction  =============
                        QTrue[kpix][1]->Fill(value, weight);
                        double ped_sub_ADC = double(value) - pedestals.at(index).first;

                        QTrue[kpix][2]->Fill(ped_sub_ADC, weight);
                        double ped_sub_fC = ped_sub_ADC/calibs.at(index).first;

                        QTrue[kpix][3]->Fill(ped_sub_fC, weight);


                        double true_charge = ped_sub_fC - common_modes_median[kpix].at(event.eventNumber());
                        double true_SoN = true_charge/noise[kpix][channel];
						//// ========= Event cut ============
                        QTrue[kpix][0]->Fill(true_charge, weight);
                        SoNTrue[kpix]->Fill(true_SoN, weight);

                        charge_v_position[sensor]->Fill(yParameterSensor(strip, sensor), true_charge);
                        charge_v_strip[sensor]->Fill(strip, true_charge);
                        //cout <<"DEBUG3.2.1 " << kpix << " " <<  channel << endl;
                        if ( true_charge > SoNCut*noise[kpix][channel] && strip != 9999)
//                        if ( true_charge > 3*noise[kpix][channel] && strip != 9999 && noise_mask[sensor].at(strip) == 1)  //only events with charge higher than 3 sigma of the noise are taken and with their charge being lower than 10 fC (to cut out weird channels), in addition no noise masked channels and no disconnected channels.
						{
						  cluster_Events_after_cut[sensor].insert(std::pair<int, double>(strip, true_charge));
						  cluster_Noise_after_cut[sensor].insert(std::pair<int, double>(strip, noise[kpix][channel]));
						  eventSample = x;
//                            cout << "Time check " << tstamp << endl;
						  if (noise[kpix][channel] == 0) cout << "Something is going wrong here" << endl;
						  
						  pixelFile << setw(5) << event.eventNumber()  << "," << setw(1) << sensor2layer(sensor)  << ","  << setw(7) << strip  << ","
                                                                   << setw(7) << true_charge/noise[kpix][channel] << ","
                                                                   << setw(2) << 1 << ","
                                                                   << setw(7) << true_charge*6241.5 << ","
                                       << tmp.str().c_str() << ","
							    << setw(7) << "irr" << setw(6) << value << "," << setw(6) <<  ped_sub_ADC << "," << setw(6) << ped_sub_fC << "," << setw(6) << true_charge << "," << setw(6) << true_SoN << setw(6) << common_modes_median[kpix].at(event.eventNumber());
                                                                   << endl;

						}
                        //cout <<"DEBUG3.3 " << kpix << " " <<  channel << endl;
						//if (kpix == 0 && (channel == 9 || channel == 105)) cout << "Something is going wrong here" << endl;
						
					}
				}
			}

            //cout <<"DEBUG4 " << kpix << " " <<  channel << endl;

		}
//        cout << "Next event" << endl;
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
						while (Input.Elements.size() != 0) // Keep repeating the clustering until either there are no valid candidates left
						{
							PacMan NomNom;
							//cout << "Maximum Signal over Noise Strip " << Input.MaxSoN() << endl;
							//cout << "Maximum Charge Strip " << Input.MaxCharge_w_Noise() << endl;
//                            cout << "Size of input element: " << Input.Elements.size() << endl;
//                            cout << "Max input element: " << Input.MaxSoN() << endl;
							double SoN_order = 0;
							Max_SoN[sensor]->Fill(Input.MaxSoN(), weight);
							MaximumSoN[sensor][Input.MaxSoN()]+=weight;
                            NomNom.Eater(Input, Input.MaxSoN(), 9999, 99999, clusterID);
                            clusterID++;
                            clustr Cluster = NomNom.getCluster();
							if (num_of_clusters == 0)
							{
								cluster_position_y[sensor][0]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor)/1000, weight);
								cluster_charge[sensor][0]->Fill(NomNom.getClusterCharge(), weight);
								cluster_size[sensor][0]->Fill(NomNom.getClusterElementssize(), weight);
								cluster_significance[sensor][0]->Fill(NomNom.getClusterSignificance(), weight);
								cluster_significance2[sensor][0]->Fill(NomNom.getClusterSignificance2(), weight);
								cluster_sigma[sensor][0]->Fill(NomNom.getClusterSigma(), weight);
								multi_cluster[0][sensor].push_back(NomNom.getCluster());

								//cout << "Size: " << NomNom.getClusterElementssize() << endl;
							}

							cluster_position_y[sensor][1]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor)/1000, weight);

							cluster_size[sensor][1]->Fill(NomNom.getClusterElementssize(), weight);
							cluster_significance[sensor][1]->Fill(NomNom.getClusterSignificance(), weight);
							cluster_significance2[sensor][1]->Fill(NomNom.getClusterSignificance2(), weight);
							cluster_charge[sensor][1]->Fill(NomNom.getClusterCharge(), weight);
							cluster_sigma[sensor][1]->Fill(NomNom.getClusterSigma(), weight);
                            multi_cluster[1][sensor].push_back(NomNom.getCluster());
//                            cout << "Cross checking whether time assumption is correct " << eventTime << endl;

							//						all_clusters[sensor].push_back(NomNom.getCluster());
//							cout << "Sensor number: " << sensor << endl;
//							all_clusters_pointer[sensor][event.eventNumber()]->push_back(NomNom.getCluster());
                            vector_charge.push_back(NomNom.getClusterCharge());
                            vector_sensor.push_back(sensor);
                            vector_sigma.push_back(NomNom.getClusterSignificance2());
                            vector_size.push_back(NomNom.getClusterElementssize());
                            vector_pos.push_back(yParameterSensor(NomNom.getClusterCoG(), sensor));
                            vector_ID.push_back(NomNom.getClusterID());


                            for (auto const& i : Cluster.Elements)
                            {
                                strip_pos.push_back(yParameterSensor(i.first,sensor));
                                strip_charge.push_back(i.second);
                                strip_signi.push_back(i.second/Cluster.Noise.at(i.first));
                                strip_sensor.push_back(sensor);
                                strip_ID.push_back(Cluster.ID);
                            }

                            stripTree->Fill();
                            strip_pos.clear();
                            strip_charge.clear();
                            strip_signi.clear();
                            strip_sensor.clear();
                            strip_ID.clear();

                            if (NomNom.getClusterSignificance2() >= 6.0)
							{
								cluster_position_y[sensor][2]->Fill(yParameterSensor(NomNom.getClusterCoG(), sensor)/1000, weight);
								cluster_charge[sensor][2]->Fill(NomNom.getClusterCharge(), weight);
								cluster_size[sensor][2]->Fill(NomNom.getClusterElementssize(), weight);
								cluster_significance[sensor][2]->Fill(NomNom.getClusterSignificance(), weight);
								cluster_significance2[sensor][2]->Fill(NomNom.getClusterSignificance2(), weight);
								cluster_sigma[sensor][2]->Fill(NomNom.getClusterSigma(), weight);
								multi_cluster[2][sensor].push_back(NomNom.getCluster());
							}

							if (header == 1){
								header = 0;
								claus_file <<"Event Number,Layer,position,Significance2,Size,Charge,runtime,runtime_ns,trigN,ClusterID" << endl;
								
							}
							claus_file << setw(5) << event.eventNumber()  << ","
                                       << setw(1) << sensor2layer(sensor)  << ","
							           << setw(7) << yParameterSensor(NomNom.getClusterCoG(), sensor)  << ","
							           << setw(7) << NomNom.getClusterSignificance2() << ","
							           << setw(2) << NomNom.getClusterElementssize() << ","
							           << setw(7) << NomNom.getClusterCharge() << ","
                                       << tmp.str().c_str() << ","
                                       << setw(7) << NomNom.getClusterID()
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
//				for (int sensor1 = 0; sensor1 < n_kpix/2; sensor1++)
//				{
//					for (int sensor2 = sensor1+1; sensor2 < n_kpix/2; sensor2++)
//					{
//						for (auto const& s1 : multi_cluster[0][sensor1])
//						{
//							for (auto const& s2 : multi_cluster[0][sensor2])
//							{
//								double y1 = yParameterSensor(s1.CoG, sensor1);
//								double y2 = yParameterSensor(s2.CoG, sensor2);
//								double clstroffset_y  = y1 - y2;
//								cluster_offset_y[0][sensor1][sensor2]->Fill(clstroffset_y, weight);
//								cluster_correlation[0][sensor1][sensor2]->Fill(y1,y2);
//							}
//						}
//						for (auto const& s1 : multi_cluster[1][sensor1])
//						{
//							for (auto const& s2 : multi_cluster[1][sensor2])
//							{
//								double y1 = yParameterSensor(s1.CoG, sensor1);
//								double y2 = yParameterSensor(s2.CoG, sensor2);
//								double clstroffset_y  = y1 - y2;
//								cluster_offset_y[1][sensor1][sensor2]->Fill(clstroffset_y, weight);
//								cluster_correlation[1][sensor1][sensor2]->Fill(y1,y2);
//							}
//						}
//						for (auto const& s1 : multi_cluster[2][sensor1])
//						{
//							for (auto const& s2 : multi_cluster[2][sensor2])
//							{
//								double y1 = yParameterSensor(s1.CoG, sensor1);
//								double y2 = yParameterSensor(s2.CoG, sensor2);
//								double clstroffset_y  = y1 - y2;
//								cluster_offset_y[2][sensor1][sensor2]->Fill(clstroffset_y, weight);
//								cluster_correlation[2][sensor1][sensor2]->Fill(y1,y2);
//							}
//						}
//                        //cout <<"DEBUG5" << endl;
//					}

//				}
			}
            clusterTree->Fill();
            vector_charge.clear();
            vector_sigma.clear();
            vector_size.clear();
            vector_pos.clear();
            vector_sensor.clear();
            vector_ID.clear();
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


    //cout <<"DEBUG6" << endl;
	claus_file.close();
	pixelFile.close();
//    marcel_file.close();
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

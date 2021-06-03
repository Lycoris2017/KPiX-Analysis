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
	
    KpixEvent              event;    //
	
    DataRead               dataRead;  //kpix event classes used for analysis of binary data
	off_t                  fileSize;  //
	off_t                  filePos;   //
	
	const unsigned int n_buckets = 4;
	const unsigned int n_kpix = 32;//24;
	const unsigned int n_channels = 1024;
	
	
	uint                   lastPct;
	uint                   currPct;
//	bool                   bucketFound[n_kpix][n_channels][n_buckets] = {false};
	string                 serial;
	


	uint                   acqCount = 0; // acquisitionCount
	uint                   acqProcessed;

    stringstream           tmp_;

    bool                   kpixFound_[n_kpix] = {false}; // variable that gives true if a kpix at the index n (0<=n<32) was found
    bool                   channelFound_[n_kpix][n_channels] = {false};
	
	// Calibration slope, is filled when 
    //std::bitset<18> index;
    //std::unordered_map<std::bitset<18>, std::pair<double, double>> calibs;
    //std::unordered_map<std::bitset<18>, std::pair<double, double>> pedestals;

    std::unordered_map<uint, std::pair<double, double>> calibs;
    std::unordered_map<uint, std::pair<double, double>> pedestals;

	
	int						calibration_check = 0;
	
	int maxAcquisitions = 100000;

	
	
	//////////////////////////////////////////
	// Skip cycle/Read calibration file
	//////////////////////////////////////////
	
	// Data file is the first and only arg
	if (argc < 4) {
	cout << "Usage: ./analysis data_file [skip_cycles_front (int)]||[read calibration input file (char)] \n";
	return(1);
	}
    else {
		//cout << "Even more debug " << strtol(argv[2], &end, 10) << endl;

        cout << " -- Reading " << argv[2] << " as calibration input file." << endl;
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
            auto index = keyhash(kpix_calib, channel_calib, bucket_calib);
            std::pair<double,double> calib_entry = std::make_pair(slope_calib, pearsson_calib);
            calibs.emplace(index,calib_entry);
        }
        cout << " -- Reading " << argv[3] << " as pedestal subtraction input file." << endl;
        TFile *pedestal_file = TFile::Open(argv[3]);

        TTree* pedestal_tree = (TTree*)pedestal_file->Get("pedestal_tree");

        int kpix_pedestal, channel_pedestal, bucket_pedestal;
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
            auto index = keyhash(kpix_pedestal, channel_pedestal, bucket_pedestal);
            std::pair<double,double> pedestal_entry = std::make_pair(pedestal_median, pedestal_MAD);
            pedestals.emplace(index, pedestal_entry);
        }



        gROOT->GetListOfFiles()->Remove(calibration_file);
        calibration_file->Close();
        gROOT->GetListOfFiles()->Remove(pedestal_file);
        pedestal_file->Close();


		
	
	}
	
	
	//cout << "DEBUG: 2" << endl;
	//////////////////////////////////////////
	// Open Data file
	//////////////////////////////////////////
	
	if ( ! dataRead.open(argv[1])  ) {
		cout << "Error opening data file " << argv[1] << endl;
		return(1);
	}
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
	

	vector<double>* vec_corr_charge[n_kpix] = {nullptr}; //delete does not work if I do not initialize all vectors
	std::map<int, double> common_modes_median[n_kpix];
    while ( dataRead.next(&event)  &&  event.eventNumber() <= maxAcquisitions) // event read to check for filled channels and kpix to reduce number of empty histograms.
	{
		
		acqCount++;		
        acqProcessed++;


        for (auto x=0; x < event.count(); x++)
        {

            //// Get sample
            auto sample  = event.sample(x);
            auto kpix    = sample->getKpixAddress();
            auto channel = sample->getKpixChannel();
            auto bucket  = sample->getKpixBucket();
            auto value   = sample->getSampleValue();
            auto type    = sample->getSampleType();
            auto index = keyhash(kpix, channel, bucket);
            if (bucket >= n_buckets) continue;

            if ( type == KpixSample::Data )
            {
                //cout << kpix << endl;
                kpixFound_[kpix]          = true;
                channelFound_[kpix][channel] = true;
                    if (calibration_check == 1 && bucket == 0)
                    {
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

                            double corrected_charge_value_median = double(value - pedestals.at(index).first)/calibs.at(index).first;
                            vec_corr_charge[kpix]->push_back(corrected_charge_value_median);
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
                //cout << "Common modes median of EventNumber " << event.eventNumber()  << " kpix " << k  << " entry " << common_modes_median[k].at(event.eventNumber()) << endl;
                vec_corr_charge[k] = nullptr;

            }
        }
	}
	maxAcquisitions = acqCount;
    dataRead.close();

	//////////////////////////////////////////
	// Data read for all events for detailed look
	/////////////////////////////////////////SoNCut/
	dataRead.open(argv[1]); //open file again to start from the beginning
	std::vector<double>* corrected_charge_vec[n_kpix][n_channels][n_buckets] = {nullptr};
    while ( dataRead.next(&event)  &&  event.eventNumber() <= maxAcquisitions) //preread to determine noise value of each channel in each KPiX.
	{
        for (auto x=0; x < event.count(); x++)
        {
            //// Get sample
            auto sample  = event.sample(x);
            auto kpix    = sample->getKpixAddress();
            auto channel = sample->getKpixChannel();
            auto bucket  = sample->getKpixBucket();
            auto value   = sample->getSampleValue();
            auto type    = sample->getSampleType();
            auto index = keyhash(kpix, channel, bucket);
            if (bucket >= n_buckets) continue;
			if ( type == KpixSample::Data ) // If event is of type KPiX data
			{
				 if (bucket == 0)
				 {
                    if (pedestals.at(index).second != 0 && calibs.at(index).first != 0  && calibs.at(index).second > 0.85 ) //ensuring we ignore 0 MAD channels
					{
                        double charge_ped_corrected = value - pedestals.at(index).first;
                        double corrected_charge_value_median = charge_ped_corrected/calibs.at(index).first;
                        double true_charge = corrected_charge_value_median - common_modes_median[kpix].at(event.eventNumber());
						if (corrected_charge_vec[kpix][channel][bucket] == NULL)
						{
							corrected_charge_vec[kpix][channel][bucket] = new std::vector<double>;
							if (corrected_charge_vec[kpix][channel][bucket] == NULL)
							{
                                std::cerr << "Memory allocation error " << endl;
								exit(-1); // probably best to bail out
							}
						} 
                        corrected_charge_vec[kpix][channel][bucket]->push_back(true_charge);
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
//    marcel_file << "Sensor, Time (BCC), Noise (fC)" << endl;

    //cout << "DEBUG 1" << endl;
	 // END OF PREREAD

    // Open root file

    tmp_.str("");
    tmp_ << argv[1] << ".tree_noise.root";
    auto outRoot = tmp_.str();

    double noise_num;
    uint kpix_num, channel_num, bucket_num;
    TFile *rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
    rFile->cd(); // move into root folder base

    TTree* noise_tree = new TTree("noise_tree", "A ROOT Tree");
    noise_tree->Branch("kpix_num", &kpix_num, "kpix_num/I");
    noise_tree->Branch("channel_num", &channel_num, "channel_num/I");
    noise_tree->Branch("bucket_num", &bucket_num, "bucket_num/I");
    noise_tree->Branch("noise_num", &noise_num, "noise_num/D");



	 // BEGIN OF NOISE CALCULATION

    for (kpix_num = 0; kpix_num < n_kpix; ++kpix_num)
	{
        if (kpixFound_[kpix_num])
		{
            for (channel_num = 0; channel_num < n_channels; ++channel_num)
			{
                if (channelFound_[kpix_num][channel_num])
				{
//				  for (bucket_num = 0; bucket_num < n_buckets; ++bucket_num)
//				    {
					  auto index = keyhash(kpix_num,channel_num, 0);
				      if (pedestals.at(index).second != 0 && calibs.at(index).first != 0  && calibs.at(index).second > 0.85 )
					{
					  noise_num = 1.4826*MAD(corrected_charge_vec[kpix_num][channel_num][0]);
					  noise_tree->Fill();
					}
//				    }
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

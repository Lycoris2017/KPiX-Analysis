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

//	const unsigned int n_buckets = 4;
    const unsigned int n_kpix = 32;//24;
//	const unsigned int n_blocks = 32;
    const unsigned int n_channels = 1024;
//	const unsigned int n_BCC = 8192;
//	const unsigned int n_strips = 1840;

    // cycles to skip in front:
    long int                   skip_cycles_front;
    FILE*                  f_skipped_cycles;
    string                 outtxt;

    string                 calState;
    uint                   lastPct;
    uint                   currPct;
//	bool                   kpixFound[n_kpix] = {false}; // variable that gives true if a kpix at the index n (0<=n<32) was found
//	bool                   channelFound[n_kpix][n_channels] = {false};
//	bool                   bucketFound[n_kpix][n_channels][n_buckets] = {false};
    uint                   x;
    uint                   value;
    uint                   kpix;
    uint                   channel;
    uint                   bucket;
    uint                    time;
    double                  tstamp;
//	uint 					subCount;
//	double 					bunchClk;
    double 					y;
    string                 serial;
    KpixSample::SampleType type;



    // Stringstream initialization for histogram naming
    stringstream           tmp;
    stringstream           tmp_units;

    // Stringstream initialization for folder naming

    stringstream			FolderName;

    uint                   acqCount = 0; // acquisitionCount
    uint                   acqProcessed;
    string                 outRoot;
    TFile					*rFile;

    XmlVariables           config;

    // Calibration slope, is filled when
    //std::bitset<18> index;
    //std::unordered_map<std::bitset<18>, std::pair<double, double>> calibs;
    //std::unordered_map<std::bitset<18>, std::pair<double, double>> pedestals;

    uint index;
    std::unordered_map<uint, std::pair<double, double>> calibs;
    std::unordered_map<uint, std::pair<double, double>> pedestals;

//    std::unordered_map<uint, std::pair<uint, double>> KPIX_temp;

    double                  pearsson_cut = 0.8;

    //int						pedestal_check = 0;

    int maxAcquisitions = 50000;

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



    if ( ! dataRead.open(argv[1])  ) {
        cout << "Error opening data file " << argv[1] << endl;
        return(1);
    }
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
    tmp << argv[1] << "_" << pedestalname << ".noise.root";
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
    printf(ANSI_COLOR_YELLOW "(Warning) Maximum number of cycles that are being considered is set to %i" ANSI_COLOR_RESET "\n", maxAcquisitions);
    cout << "\rReading File: 0 %" << flush;  // Printout of progress bar
    //goodTimes       	= 0;

    // Open root file
    rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
    rFile->cd(); // move into root folder base





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
                time  = sample->getSampleTime();
                channel = sample->getKpixChannel();
                bucket  = sample->getKpixBucket();
                value   = sample->getSampleValue();
                type    = sample->getSampleType();
                int eventNum = sample->getEventNum();
                index = keyhash(kpix, channel, bucket);



                //cout << type <<endl;
                //cout << "DEBUG 2" << endl;
//                if (type == KpixSample::Temperature)
//                {
//                    cout << " KPiX " << kpix << " Temperature " << value << " Time " << tstamp << " eventNum " << eventNum << endl;
//                }

                if ( type == KpixSample::Data )
                {
                    //cout << kpix << endl;
//					bucketFound[kpix][channel][bucket] = true;
                    //cout << "Found KPIX " << kpix << endl;

                    if (bucket == 0)
                    {
//                        cout << "Just testing something: Index " << index << " KPiX " << kpix << " Channel " << channel << " Bucket " << bucket << endl;
//                        cout << "Previous calib slope " << calib_slope[kpix][channel][bucket] << endl;
//                        cout << "New calib slope " << calibs.at(index).first << endl;
                            //cout << "2nd Pedestal MAD " << pedestals.at(index).second << " kpix " << kpix << " channel " << channel << " bucket " << bucket << endl;
                            if (pedestals.at(index).second != 0 && calibs.at(index).first != 0  && calibs.at(index).second > 0.85)
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
                                double corrected_charge_value_median = double(value - pedestals.at(index).first)/calibs.at(index).first;
                                vec_corr_charge[kpix]->push_back(corrected_charge_value_median);
                                //cout << "Pointer of vec_corr_charge " << vec_corr_charge[kpix] << endl;
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

    dataRead.close();
    double weight = 1.0/acqProcessed;

    std::vector<double> vec_noise, vec_charge;
    std::vector<int> vec_kpix, vec_channel, vec_strip, vec_bucket, vec_time;
    //////////////////////////////////////////
    // TTree
    //////////////////////////////////////////

    TTree* noiseTimeTree = new TTree("noise_time_tree", "noise_time_tree");
    noiseTimeTree->Branch("kpix", &vec_kpix);
    noiseTimeTree->Branch("channel", &vec_channel);
    noiseTimeTree->Branch("strip", &vec_strip);
    noiseTimeTree->Branch("bucket", &vec_bucket);
    noiseTimeTree->Branch("time", &vec_time);
    noiseTimeTree->Branch("noise", &vec_noise);
    noiseTimeTree->Branch("charge", &vec_charge);
    TH1F* temp_hist = new TH1F("temp_hist", "temp_hist; temperature; Nr. of Entries", 256,-62, 90 );


    //TH1F* mean_noise = new TH1F("mean_noise_left", "mean_noise; noise(fC); entries", 100, -0.05, 0.95);


    //////////////////////////////////////////
    // Data read for all events for detailed look
    //////////////////////////////////////////
    dataRead.open(argv[1]); //open file again to start from the beginning







    rFile->cd();
    std::unordered_map<uint, std::vector<double>> QTrue;
    std::unordered_map<uint, double> noise;
    while ( dataRead.next(&event)  &&  event.eventNumber() <= maxAcquisitions) //preread to determine noise value of each channel in each KPiX.
    {
        std::vector<double> time_ext;


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
            time  = sample->getSampleTime();
            int strip = 9999;
            if (kpix%2 == 0) strip = kpix2strip_left.at(channel);// if left kpix
            else strip  = kpix2strip_right.at(channel); // if right kpix
            index = keyhash(kpix, channel, bucket);
            uint index1 = keyhash(kpix, channel, bucket, time);
            //cout << "DEBUG: " << kpix << " " << channel << " " << value << " " << type << endl;
            //cout << "DEBUG 0" << endl;
            if (type == KpixSample::Temperature)
            {
//                cout << "Temperature data" << endl;
//                std::string test = "";
                uint temp = value%256; // supposedly only the first 8 bit are temperature values
//                uint graytemp = time%256;
//                test = bin(temp, test);
//                cout << test << endl;
//                std::string greytest = binarytoGray(test);
//                std::string greytstamp = "";
//                greytstamp = bin(graytemp, greytstamp);
//                cout << greytest << " : " << greytstamp << endl;
                double tmp=0.598*(255-temp)-62;
                temp_hist->Fill(tmp);

            }
            if ( type == KpixSample::Data ) // If event is of type KPiX data
            {

                    if (pedestals.at(index).second != 0 && calibs.at(index).first != 0  && calibs.at(index).second > 0.85 ) //ensuring we ignore 0 MAD+Slope channels
                    {
                        //cout << "DEBUG: " << kpix << " " << channel << " " << value << " " << type << endl;
                        //if (kpix == 0 && channel == 188) cout << "Just testing things: " << calibs.at(index).first << " " << calibs.at(index).second;
                        double charge_ped_corrected = value - pedestals.at(index).first;
                        //cout << "DEBUG 1" << endl;
                        double corrected_charge_value_median = charge_ped_corrected/calibs.at(index).first;
                        double true_charge = corrected_charge_value_median - common_modes_median[kpix].at(event.eventNumber());
                        vec_charge.push_back(true_charge);
                        vec_bucket.push_back(bucket);
                        vec_channel.push_back(channel);
                        vec_kpix.push_back(kpix);
                        vec_strip.push_back(strip);
                        vec_time.push_back(tstamp);
                        //cout << "DEBUG 2" << endl;
                        //cout << "Time standard: " << tstamp <<  " and time index recalculated: " << index1/pow(2,17) << endl;
                        //cout << "Channel standard: " << channel <<  " and channel index recalculated: " << (index1 % int(pow(2,17)))/pow(2,7) << endl;
//                        cout << "KPiX standard: " << kpix <<  " and kpix index recalculated: " << ((index1 % pow(2,17))%pow(2,7))/pow(2,2) << endl;
//                        cout << "Bucket standard: " << bucket <<  " and bucket index recalculated: " << ((index1 % pow(2,17))%pow(2,7))%pow(2,2) << endl;
//                        QTrue.at(index1).push_back(charge_CM_corrected);
//                        cout << "CM corrected charge is " << charge_CM_corrected << " for KPiX " << kpix << " and channel " << channel << endl;


                    //else if (pedestals.at(index).second == 0 && kpix == 0) cout << "1KPIX " << kpix << " Channel " << channel << endl;
                    }
            }
            noiseTimeTree->Fill();
            vec_charge.clear();
            vec_bucket.clear();
            vec_channel.clear();
            vec_kpix.clear();
            vec_strip.clear();
            vec_time.clear();
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

    //cout << "DEBUG 2" << endl;
    /////////////////////////////////////////////////////////
    //// BEGIN OF CLUSTER READ
    dataRead.open(argv[1]);

//	std::vector<clustr> all_clusters[n_kpix/2];

    cout << endl;
    cout << "Writing root plots to " << outRoot << endl;
    cout << endl;

    rFile->Write();
    gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
    rFile->Close();



    dataRead.close();
    return(0);
}

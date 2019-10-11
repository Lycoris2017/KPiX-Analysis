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
struct tree_cluster_input
{
	double CoG ;
	double Charge;
	double Significance2; //Adding noises in quadrature and later dividing chargesum by this noise sum
	double Sigma;
	int Size;
};

//////////////////////////////////////////
// Functions
//////////////////////////////////////////

Double_t test2(Double_t *x, Double_t *par)
{
    return TMath::Landau(x[0], par[2], par[1]);
}
//Double_t landau2(Double_t x, Double_t par)
//{
//    return TMath::Landau(x, 4.6, 1);
//}

//Double_t fitFunction(Double_t x, Double_t par)
//{
//    return landau1(x,par)+landau2(x,par);
//}

//////////////////////////////////////////
// Begin of main analysis
//////////////////////////////////////////
int main ( int argc, char **argv )
{
	TH1::SetDefaultSumw2();

	const unsigned int		n_buckets = 4;
    const unsigned int		n_kpix = 32;//24;
	const unsigned int		n_sensor = n_kpix/2;//24;
//	const unsigned int		n_blocks = 32;
	const unsigned int		n_channels = 1024;
	const unsigned int		n_BCC = 8192;
//	const unsigned int		n_strips = 1840;
	unsigned int			n_events = 0;

	uint                   lastPct;
	uint                   currPct;

	bool					sensorFound[n_sensor] = {false};
	stringstream			FolderName;

	TFile					*rFile;
	int eventNumber, sensor, eventTime;
	tree_cluster_input		tree_clusters;
	stringstream           tmp;
	stringstream           tmp_units;
	string                 outRoot;

    TH1F 					*cluster_position_y[n_sensor][4];
	//TH1F 					*cluster_position_x[n_sensor][3];
    TH1F 					*num_clusters[n_sensor][4]; // only 2 because it makes no sense to ask for the number of clusters in the seed case cause it will always be 1.
    TH1F 					*cluster_charge[n_sensor][4];
    TH1F 					*cluster_significance2[n_sensor][4];
    TH1F 					*cluster_size[n_sensor][4];
    TH1F 					*cluster_sigma[n_sensor][4];

    TH2F					*cluster_correlation[3][n_sensor][n_sensor-1];
    TH1F                	*cluster_offset_y[3][n_sensor][n_sensor-1];
//	TH1F					*cluster_offset_x[4][n_kpix/2][n_kpix/2-1];


	string InName = argv[1];

	size_t name_start  = InName.find("/Run") + 1;
	size_t name_length = InName.find(".Tree") - name_start;

	InName = InName.substr(name_start, name_length);

	cout << "Name of output file is " <<  InName << endl;
	tmp.str("");
    tmp << InName << ".cluster_analysed.root";
	outRoot = tmp.str();
	cout << "\rReading File: 0 %" << flush;
	rFile = new TFile(outRoot.c_str(),"recreate"); // produce root file
	rFile->cd(); // move into root folder base

	//////////////////////////////////////////
	// Read Data
	//////////////////////////////////////////
	TFile *clusterFile = TFile::Open(argv[1]);
	TTree* clusterTree = (TTree*)clusterFile->Get("clusterTree");

	clusterTree->SetBranchAddress("CoG", &tree_clusters.CoG);
	clusterTree->SetBranchAddress("Size", &tree_clusters.Size);
	clusterTree->SetBranchAddress("Significance2", &tree_clusters.Significance2);
	clusterTree->SetBranchAddress("Sigma", &tree_clusters.Sigma);
	clusterTree->SetBranchAddress("Charge", &tree_clusters.Charge);
	clusterTree->SetBranchAddress("eventNumber", &eventNumber);
    clusterTree->SetBranchAddress("eventTime", &eventTime);
	clusterTree->SetBranchAddress("sensor", &sensor);

	long int NEntries = clusterTree->GetEntries();
	for (long int i = 0; i < NEntries; ++i)
	{
		clusterTree->GetEntry(i);
//		cout << "Sensor " << sensor << endl;
		sensorFound[sensor]          = true;
		if (eventNumber >= n_events)
			n_events = eventNumber;
	}
    cout << "DEBUG 1" << endl;
	for (sensor = 0; sensor < n_sensor; sensor++) //looping through all possible kpix
	{
		if (sensorFound[sensor])
		{
			rFile->cd(); //producing subfolder for kpix same as above for the event subfolder structure
			FolderName.str("");
			FolderName << "Sensor_" << sensor;
			rFile->mkdir(FolderName.str().c_str());
			TDirectory *sensor_folder = rFile->GetDirectory(FolderName.str().c_str());
			sensor_folder->cd();
			for (int s = sensor+1; s < n_sensor; ++s)
			{
				if (sensorFound[s])
				{
					tmp.str("");
					tmp << "cluster_offset_all_y_sens" << sensor << "_to_sens" << s << "_b0";
					cluster_offset_y[0][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
					//			tmp.str("");
					//			tmp << "cluster_offset_all_x_sens" << sensor << "_to_sens" << s << "_b0";
					//			cluster_offset_x[0][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
					tmp.str("");
					tmp << "cluster_correlation_all_sens" << sensor << "_to_sens" << s << "_b0";
					tmp_units.str("");
					tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
					cluster_correlation[0][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

					tmp.str("");
					tmp << "cluster_offset_CUTS_y_sens" << sensor << "_to_sens" << s << "_b0";
					cluster_offset_y[1][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
					//			tmp.str("");
					//			tmp << "cluster_offset_CUTS_x_sens" << sensor << "_to_sens" << s << "_b0";
					//			cluster_offset_x[1][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
					tmp.str("");
					tmp << "cluster_correlation_CUTS_sens" << sensor << "_to_sens" << s << "_b0";
					tmp_units.str("");
					tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
					cluster_correlation[1][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);

                    tmp.str("");
                    tmp << "cluster_offset_CORR_y_sens" << sensor << "_to_sens" << s << "_b0";
                    cluster_offset_y[2][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 100, -10000, 10000);
                    //			tmp.str("");
                    //			tmp << "cluster_offset_CUTS_x_sens" << sensor << "_to_sens" << s << "_b0";
                    //			cluster_offset_x[1][sensor][s] = new TH1F(tmp.str().c_str(), "offset; #mum; #entries", 200, -4000, 4000);
                    tmp.str("");
                    tmp << "cluster_correlation_CORR_sens" << sensor << "_to_sens" << s << "_b0";
                    tmp_units.str("");
                    tmp_units << "Strip correlation; Sensor " << sensor <<  " | Position (#mum); Sensor " << s << " | Position (#mum)";
                    cluster_correlation[2][sensor][s] = new TH2F(tmp.str().c_str(), tmp_units.str().c_str(), 230,-46000.5, 45999.5, 230,-46000, 46000);
				}
			}


			tmp.str("");
			tmp << "cluster_position_y_all_sensor" << sensor << "_b0";
			cluster_position_y[sensor][0] = new TH1F(tmp.str().c_str(), "cluster position y1; #mum; #Entries", 1840,-46000, 46000);
			tmp.str("");
			tmp << "cluster_charge_all_sensor" << sensor << "_b0";
			cluster_charge[sensor][0] = new TH1F(tmp.str().c_str(), "cluster charge1; Charge (fC); #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance2_all_sensor" << sensor << "_b0";
			cluster_significance2[sensor][0] = new TH1F(tmp.str().c_str(), "cluster significance2_1; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_sigma_all_sensor" << sensor << "_b0";
			cluster_sigma[sensor][0] = new TH1F(tmp.str().c_str(), "cluster sigma1; #sigma; #Entries", 200,-0.5, 4.95);
			tmp.str("");
			tmp << "cluster_size_all_sensor" << sensor << "_b0";
			cluster_size[sensor][0] = new TH1F(tmp.str().c_str(), "cluster size1; Size; #Entries", 10,-0.5, 9.5);
			tmp.str("");
			tmp << "clusters_all_sensor" << sensor << "_b0";
			num_clusters[sensor][0] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; #Entries", 100,-0.5, 99.5);

			tmp.str("");
            tmp << "cluster_position_y_CUTS_sensor" << sensor << "_b0";
            cluster_position_y[sensor][1] = new TH1F(tmp.str().c_str(), "cluster position y; #mum; #Entries", 1840,-46000, 46000);
			tmp.str("");
			tmp << "cluster_charge_CUTS_sensor" << sensor << "_b0";
            cluster_charge[sensor][1] = new TH1F(tmp.str().c_str(), "cluster charge; Charge (fC); #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance2_CUTS_sensor" << sensor << "_b0";
            cluster_significance2[sensor][1] = new TH1F(tmp.str().c_str(), "cluster significance2; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_sigma_CUTS_sensor" << sensor << "_b0";
            cluster_sigma[sensor][1] = new TH1F(tmp.str().c_str(), "cluster sigma; #sigma; #Entries", 200,-0.5, 4.95);
			tmp.str("");
			tmp << "cluster_size_CUTS_sensor" << sensor << "_b0";
            cluster_size[sensor][1] = new TH1F(tmp.str().c_str(), "cluster size; Size; #Entries", 10,-0.5, 9.5);
			tmp.str("");
			tmp << "clusters_CUTS_sensor" << sensor << "_b0";
			num_clusters[sensor][1] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; #Entries", 100,-0.5, 99.5);

            tmp.str("");
            tmp << "cluster_position_y_CORR_sensor" << sensor << "_b0";
            cluster_position_y[sensor][2] = new TH1F(tmp.str().c_str(), "cluster position y; #mum; #Entries", 1840,-46000, 46000);
            tmp.str("");
            tmp << "cluster_charge_CORR_sensor" << sensor << "_b0";
            cluster_charge[sensor][2] = new TH1F(tmp.str().c_str(), "cluster charge; Charge (fC); #Entries", 200,-0.5, 49.5);
            tmp.str("");
            tmp << "cluster_significance2_CORR_sensor" << sensor << "_b0";
            cluster_significance2[sensor][2] = new TH1F(tmp.str().c_str(), "cluster significance2; S/N; #Entries", 200,-0.5, 49.5);
            tmp.str("");
            tmp << "cluster_sigma_CORR_sensor" << sensor << "_b0";
            cluster_sigma[sensor][2] = new TH1F(tmp.str().c_str(), "cluster sigma; #sigma; #Entries", 200,-0.5, 4.95);
            tmp.str("");
            tmp << "cluster_size_CORR_sensor" << sensor << "_b0";
            cluster_size[sensor][2] = new TH1F(tmp.str().c_str(), "cluster size; Size; #Entries", 10,-0.5, 9.5);
            tmp.str("");
            tmp << "clusters_CORR_sensor" << sensor << "_b0";
            num_clusters[sensor][2] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; #Entries", 100,-0.5, 99.5);

            tmp.str("");
            tmp << "cluster_position_y_CORR_20k_sensor" << sensor << "_b0";
            cluster_position_y[sensor][3] = new TH1F(tmp.str().c_str(), "cluster position y; #mum; #Entries", 1840,-46000, 46000);
            tmp.str("");
            tmp << "cluster_charge_CORR_20k_sensor" << sensor << "_b0";
            cluster_charge[sensor][3] = new TH1F(tmp.str().c_str(), "cluster charge; Charge (fC); #Entries", 200,-0.5, 49.5);
            tmp.str("");
            tmp << "cluster_significance2_CORR_20k_sensor" << sensor << "_b0";
            cluster_significance2[sensor][3] = new TH1F(tmp.str().c_str(), "cluster significance2; S/N; #Entries", 200,-0.5, 49.5);
            tmp.str("");
            tmp << "cluster_sigma_CORR_20k_sensor" << sensor << "_b0";
            cluster_sigma[sensor][3] = new TH1F(tmp.str().c_str(), "cluster sigma; #sigma; #Entries", 200,-0.5, 4.95);
            tmp.str("");
            tmp << "cluster_size_CORR_20k_sensor" << sensor << "_b0";
            cluster_size[sensor][3] = new TH1F(tmp.str().c_str(), "cluster size; Size; #Entries", 10,-0.5, 9.5);
            tmp.str("");
            tmp << "clusters_CORR_20k_sensor" << sensor << "_b0";
            num_clusters[sensor][3] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; #Entries", 100,-0.5, 99.5);
		}
	}
     cout << "DEBUG 2" << endl;

	uint event = 0;
	vector<tree_cluster_input> multi_clusters[n_sensor];
	vector<tree_cluster_input> multi_clusters2[n_sensor];

    vector<tree_cluster_input> multi_clustersZ[n_sensor][n_events];

    double correlation_peak[n_sensor][n_sensor-1];

	for (long int i = 0; i < NEntries; ++i)
	{
		clusterTree->GetEntry(i);
//        cout << event << " : " << eventNumber << endl;
        if (event != eventNumber) // if an event is over and we move to the next then first calculate the correlation distribution of the previous event and clear the vectors afterwards
		{
			event = eventNumber;
			for (uint sensor1 = 0; sensor1 < n_sensor; sensor1++)
			{
				for (uint sensor2 = sensor1+1; sensor2 < n_sensor; sensor2++)
				{
					for (auto const& s1 : multi_clusters[sensor1])
					{
						for (auto const& s2 : multi_clusters[sensor2])
						{
							double y1 = yParameterSensor(s1.CoG, sensor1);
							double y2 = yParameterSensor(s2.CoG, sensor2);
							double clstroffset_y  = y1 - y2;
//							cout << clstroffset_y << endl;
							cluster_offset_y[0][sensor1][sensor2]->Fill(clstroffset_y);
							cluster_correlation[0][sensor1][sensor2]->Fill(y1,y2);
						}
					}
					for (auto const& s1 : multi_clusters2[sensor1])
					{
						for (auto const& s2 : multi_clusters2[sensor2])
						{
							double y1 = yParameterSensor(s1.CoG, sensor1);
							double y2 = yParameterSensor(s2.CoG, sensor2);
							double clstroffset_y  = y1 - y2;
//							cout << clstroffset_y << endl;
							cluster_offset_y[1][sensor1][sensor2]->Fill(clstroffset_y);
							cluster_correlation[1][sensor1][sensor2]->Fill(y1,y2);
						}
					}
				}
			}
			for (auto& v:multi_clusters)
			{
				v.clear();
            }
			for (auto& v:multi_clusters2)
			{
				v.clear();
			}
		}

		cluster_position_y[sensor][0]->Fill(yParameterSensor(tree_clusters.CoG,sensor));
		cluster_charge[sensor][0]->Fill(tree_clusters.Charge);
		cluster_size[sensor][0]->Fill(tree_clusters.Size);
		cluster_significance2[sensor][0]->Fill(tree_clusters.Significance2);
		cluster_sigma[sensor][0]->Fill(tree_clusters.Sigma);
        multi_clusters[sensor].push_back(tree_clusters);

        if (tree_clusters.Significance2 >= 7.0 && tree_clusters.Charge >= 1.0 && tree_clusters.Size <= 2)
		{
			cluster_position_y[sensor][1]->Fill(yParameterSensor(tree_clusters.CoG,sensor));
			cluster_charge[sensor][1]->Fill(tree_clusters.Charge);
			cluster_size[sensor][1]->Fill(tree_clusters.Size);
			cluster_significance2[sensor][1]->Fill(tree_clusters.Significance2);
			cluster_sigma[sensor][1]->Fill(tree_clusters.Sigma);
			multi_clusters2[sensor].push_back(tree_clusters);
            multi_clustersZ[sensor][eventNumber].push_back(tree_clusters);
		}



		currPct = (uint)(((double)i / (double)NEntries) * 100.0);
		if ( currPct != lastPct )
		{
			cout << "\rReading File: " << currPct << " %      " << flush;
			lastPct = currPct;
		}
	}
    for (uint sensor1 = 0; sensor1 < n_sensor; sensor1++)
    {
        if (sensorFound[sensor1])
        {
            for (uint sensor2 = sensor1+1; sensor2 < n_sensor; sensor2++)
            {
                if (sensorFound[sensor2])
                {
                    int maxbin = cluster_offset_y[1][sensor1][sensor2]->GetMaximumBin();
                    correlation_peak[sensor1][sensor2] = cluster_offset_y[1][sensor1][sensor2]->GetXaxis()->GetBinCenter(maxbin);
//                    cout << "max bin position is " << correlation_peak[sensor1][sensor2] << " for s1 : s2 " << sensor1 << " : " << sensor2 << endl;
                }
            }
        }
    }
    for (auto& v:multi_clusters2)
    {
        v.clear();
    }
    event = 0;
    for (long int i = 0; i < NEntries; ++i)
    {
        clusterTree->GetEntry(i);
        if (event != eventNumber) // if an event is over and we move to the next then first calculate the correlation distribution of the previous event and clear the vectors afterwards
        {
            event = eventNumber;
            for (uint sensor1 = 0; sensor1 < n_sensor; sensor1++)
            {
                for (uint sensor2 = sensor1+1; sensor2 < n_sensor; sensor2++)
                {
                    for (auto const& s1 : multi_clusters2[sensor1])
                    {
                        bool input = true;
                        for (auto const& s2 : multi_clusters2[sensor2])
                        {

                            double y1 = yParameterSensor(s1.CoG, sensor1);
                            double y2 = yParameterSensor(s2.CoG, sensor2);
                            double clstroffset_y  = y1 - y2;
                            if ( clstroffset_y > correlation_peak[sensor1][sensor2] - 300 &&  clstroffset_y < correlation_peak[sensor1][sensor2] + 300 )
                            {
                                if (input)
                                {
                                    cluster_position_y[sensor1][2]->Fill(y1);
                                    cluster_charge[sensor1][2]->Fill(s1.Charge);
                                    cluster_size[sensor1][2]->Fill(s1.Size);
                                    cluster_significance2[sensor1][2]->Fill(s1.Significance2);
                                    cluster_sigma[sensor1][2]->Fill(s1.Sigma);
                                    if (fabs(y1) < 20000)
                                    {
                                        cluster_position_y[sensor1][3]->Fill(y1);
                                        cluster_charge[sensor1][3]->Fill(s1.Charge);
                                        cluster_size[sensor1][3]->Fill(s1.Size);
                                        cluster_significance2[sensor1][3]->Fill(s1.Significance2);
                                        cluster_sigma[sensor1][3]->Fill(s1.Sigma);
                                    }
                                    input = false;
                                }
                                cluster_offset_y[2][sensor1][sensor2]->Fill(clstroffset_y);
                                cluster_correlation[2][sensor1][sensor2]->Fill(y1,y2);

                            }
                        }
                    }
                }
            }
            for (auto& v:multi_clusters2)
            {
                v.clear();
            }
        }
        if (tree_clusters.Significance2 >= 7.0 && tree_clusters.Charge >= 1.0 && tree_clusters.Size <= 2)
        {
            multi_clusters2[sensor].push_back(tree_clusters);
        }

    }

//    for (int sensor1 = 0 ; sensor1 < n_sensor; ++sensor1)
//    {
//        if (sensorFound[sensor1])
//        {
//            for (int q = 0; q < 4; ++q)
//                cluster_charge[sensor1][q]->Fit("landau");
//        }
//    }

    TF1 *fitFcn2 = new TF1("fitFcn2", test2, 0, 15, 2);
    cluster_charge[0][3]->Fit("fitFcn2");

    cout << "DEBUG END" << endl;


	gROOT->GetListOfFiles()->Remove(clusterFile);
	clusterFile->Close();

	rFile->Write();
	gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
	rFile->Close();
	
	cout << endl;
	cout << "Writing root plots to " << outRoot << endl;

	
	return(0);
}

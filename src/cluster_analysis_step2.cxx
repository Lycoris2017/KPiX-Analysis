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

	TH1F 					*cluster_position_y[n_sensor][2];
	//TH1F 					*cluster_position_x[n_sensor][3];
	TH1F 					*num_clusters[n_sensor][2]; // only 2 because it makes no sense to ask for the number of clusters in the seed case cause it will always be 1.
	TH1F 					*cluster_charge[n_sensor][2];
	TH1F 					*cluster_significance2[n_sensor][2];
	TH1F 					*cluster_size[n_sensor][2];
	TH1F 					*cluster_sigma[n_sensor][2];

	TH2F					*cluster_correlation[4][n_sensor][n_sensor-1];
	TH1F					*cluster_offset_y[2][n_sensor][n_sensor-1];
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
			cluster_position_y[sensor][1] = new TH1F(tmp.str().c_str(), "cluster position y2; #mum; #Entries", 1840,-46000, 46000);
			tmp.str("");
			tmp << "cluster_charge_CUTS_sensor" << sensor << "_b0";
			cluster_charge[sensor][1] = new TH1F(tmp.str().c_str(), "cluster charge2; Charge (fC); #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_significance2_CUTS_sensor" << sensor << "_b0";
			cluster_significance2[sensor][1] = new TH1F(tmp.str().c_str(), "cluster significance2_2; S/N; #Entries", 200,-0.5, 49.5);
			tmp.str("");
			tmp << "cluster_sigma_CUTS_sensor" << sensor << "_b0";
			cluster_sigma[sensor][1] = new TH1F(tmp.str().c_str(), "cluster sigma2; #sigma; #Entries", 200,-0.5, 4.95);
			tmp.str("");
			tmp << "cluster_size_CUTS_sensor" << sensor << "_b0";
			cluster_size[sensor][1] = new TH1F(tmp.str().c_str(), "cluster size2; Size; #Entries", 10,-0.5, 9.5);
			tmp.str("");
			tmp << "clusters_CUTS_sensor" << sensor << "_b0";
			num_clusters[sensor][1] = new TH1F(tmp.str().c_str(), "Clusters; #Clusters; #Entries", 100,-0.5, 99.5);
		}
	}


	uint event = 0;
	vector<tree_cluster_input> multi_clusters[n_sensor];
	vector<tree_cluster_input> multi_clusters2[n_sensor];
	for (long int i = 0; i < NEntries; ++i)
	{
		clusterTree->GetEntry(i);

		if (event != eventNumber)
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

		if (tree_clusters.Significance2 >= 7.0 && tree_clusters.Charge >= 1.0 && tree_clusters.Size <= 2)
		{
			cluster_position_y[sensor][1]->Fill(yParameterSensor(tree_clusters.CoG,sensor));
			cluster_charge[sensor][1]->Fill(tree_clusters.Charge);
			cluster_size[sensor][1]->Fill(tree_clusters.Size);
			cluster_significance2[sensor][1]->Fill(tree_clusters.Significance2);
			cluster_sigma[sensor][1]->Fill(tree_clusters.Sigma);
			multi_clusters2[sensor].push_back(tree_clusters);
		}

		multi_clusters[sensor].push_back(tree_clusters);

		currPct = (uint)(((double)i / (double)NEntries) * 100.0);
		if ( currPct != lastPct )
		{
			cout << "\rReading File: " << currPct << " %      " << flush;
			lastPct = currPct;
		}
	}



	gROOT->GetListOfFiles()->Remove(clusterFile);
	clusterFile->Close();

	rFile->Write();
	gROOT->GetListOfFiles()->Remove(rFile); //delete cross links that make production of subfolder structure take forever
	rFile->Close();
	
	cout << endl;
	cout << "Writing root plots to " << outRoot << endl;

	
	return(0);
}

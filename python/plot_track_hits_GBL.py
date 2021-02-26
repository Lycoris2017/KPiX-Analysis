#!usr/bin/python

import numpy as np
import string
from matplotlib import cm
import ROOT
import argparse
import argcomplete
from operator import add
import sys
from decimal import Decimal
import re

ROOT.gROOT.SetBatch(True)


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

ROOT.TH1.SetDefaultSumw2()

parser = MyParser()
parser.add_argument('file_in', help='name of the input file')
args = parser.parse_args()

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print parser.print_help()
		sys.exit(1)
	print ''
	outfile = args.file_in[:-4]+'.root'
	outHistFile = ROOT.TFile.Open(outfile, "RECREATE")
	outHistFile.cd()

	hist_residual = ROOT.TH1F("residual", "residual; Residual (mm); Nr. of Entries", 101, -50e-3, 50e-3)
	hist_residual2 = ROOT.TH1F("residual2", "residual2; Residual (mm); Nr. of Entries", 101, -50e-3, 50e-3)

	hist_residual_timed = ROOT.TH1F("residual_timed", "residual_timed; Residual (mm); Nr. of Entries", 101, -50e-3, 50e-3)

	size_fitpos_correlation = ROOT.TH2F("size_fitpos", "size_fitpos; size; distance (mm); Number of Entries ", 5, 0.5, 5.5, 51, 0, 51e-3)
	size_fitpos_correlation_timed = ROOT.TH2F("size_fitpos_timed", "size_fitpos_timed; size; distance (mm); Number of Entries ", 5, 0.5, 5.5, 51, 0, 51e-3)

	hist_profile_thesis = []
	hist_profile_thesis_timed = []
	hist_profile_thesis_timed2 = []
	for j in range(0,4,1):
		histname = "sub_cell_entries_csize"+str(j)
		hist_profile_thesis.append(ROOT.TH1F(histname, "sub_cell_entries; distance (#mum); Number of Entries ", 51, 0, 51e-3))
		hist_profile_thesis_timed.append(ROOT.TH1F(histname+"_timed", "sub_cell_entries; distance (#mum); Number of Entries ", 51, 0, 51))
		#hist_profile_thesis_timed_bump.append(ROOT.TH1F(histname+"_timed_bump", "sub_cell_entries; distance (mm); Number of Entries ", 51, 0, 51e-3))

	#size_fitpos_correlation_bump_timed = ROOT.TH2F("size_fitpos_timed", "size_fitpos_timed; size; distance (mm); Number of Entries ", 5, 0.5, 5.5, 51, 0, 51e-3)

	hist_y_track_local = ROOT.TH1F("y_track", "y_track; y (mm); Nr. of Entries;", 1840, -46, 46)
	hist_hits_layer = ROOT.TH1F("hits_all", "hits; layer; Nr. of Entries;", 6, 9.5, 15.5)
	hist_hits_layer_3p = ROOT.TH1F("hits_3p", "hits; layer; Nr. of Entries;", 6, 9.5, 15.5)
	hist_hits_layer_5p = ROOT.TH1F("hits_5p", "hits; layer; Nr. of Entries;", 6, 9.5, 15.5)
	hist_total = ROOT.TH1F("hits_total", "hits; layer; Efficiency;", 6, 9.5, 15.5)

	hist_hits_layer_untimed = ROOT.TH1F("hits_all_untimed1", "hits; layer; Nr. of Entries;", 6, 9.5, 15.5)
	hist_hits_layer_untimed2 = ROOT.TH1F("hits_all_untimed2", "hits; layer; Nr. of Entries;", 6, 9.5, 15.5)
	hits_hits_on_lycoris = ROOT.TH1F("hits_on_lycoris", "hits; Nr. of Hits; Nr. of Entries;", 7, -0.5, 6.5)

#hist_pos_track = ROOT.TH1F("track_pos", "track_pos; x (mm); y (mm); Nr. of Entries;", 200, -3, 3, 200, -3, 3)

	hist_charge = []
	hist_charge_timed = []
	for j in range(0,4,1):
		if j == 0:
			histname = "charge_csize_all"
		else:
			histname = "charge_csize"+str(j)
		hist_charge.append(ROOT.TH1F(histname, "charge; Charge (fC); Nr. of Entries ",  151, 0, 15))
		hist_charge_timed.append(ROOT.TH1F(histname+"_timed", "charge; Charge (fC); Nr. of Entries ",  151, 0, 15))

	hist_charge_floating = ROOT.TH1F("charge_floating", "charge_floating; Charge (fC); Nr. of Entries ",  151, 0, 15)
	hist_charge_floating_timed = ROOT.TH1F("charge_floating_timed", "charge_floating; Charge (fC); Number of Entries ",  201, 0, 20)
	hist_charge_readout = ROOT.TH1F("charge_readout", "charge_readout; Charge (fC); Nr. of Entries ",  151, 0, 15)
	hist_charge_readout_timed = ROOT.TH1F("charge_readout_timed", "charge_readout; Charge (fC); Number of Entries ",  201, 0, 20)
	hist_csize = ROOT.TH1F("cluster_size_distribution", "cluster_size_distribution; Size (Nr. of Strips); Nr. of Entries ", 5, 0.5, 5.5)
	hist_SbyN = ROOT.TH1F("SbyN", "SbyN; Signal over Noise; Nr. of Entries ", 100, -4.5, 44.5)

	#pEff = ROOT.TEfficiency("eff","my efficiency;x;#epsilon", 6, 9.5, 15.5);


	x_l = {}
	y_l = {}
	z_l = {}
	t_l = {}
	SbyN_l = {}
	size_l = {}
	charge_l = {}
	ID_l = {}
	x_t = {}
	y_t = {}
	chi2_t = {}
	x_t_loc = {}
	y_t_loc = {}
	x_l_loc = {}
	y_l_loc = {}
	n_lycoris = 0
	n_timepix = 0
	t_lycoris = None
	nr_tracks = 0
	nr_tracks_3p = 0
	nr_tracks_5p = 0
	nr_tracks_untimed = 0
	nr_tracks_all= 0
	with open(args.file_in) as inFile:
		for line in inFile:
			if "run" in line or "Track" in line:
				nr_tracks_all += 1
				if (n_lycoris > 0):
					nr_tracks_untimed += 1
				if (n_lycoris > 0 and 20 in t_l.keys()):
					t_diff = t_lycoris - t_l[20]
					if abs(t_diff) < 500:
						nr_tracks += 1
						for i in range(10,16,1):
							hist_total.Fill(i)
						if (n_lycoris > 2):
							nr_tracks_3p += 1
						if (n_lycoris > 4):
							nr_tracks_5p += 1
						for i in t_l.iterkeys():
							if i >= 10 and i < 20:
								hist_hits_layer.Fill(i)
								if (n_lycoris > 2):
									hist_hits_layer_3p.Fill(i)
								if (n_lycoris > 4):
									hist_hits_layer_5p.Fill(i)
						for i in range(10,16,1):
							if i in y_l_loc.iterkeys():
								size_fitpos_correlation_timed.Fill(size_l[i],y_t_loc[i]%0.05)
								temp_y = y_t_loc[i]%0.05
								temp_y = temp_y*1000.
								hist_profile_thesis_timed[0].Fill(temp_y)
								hist_charge_timed[0].Fill(charge_l[i])
								if (temp_y < 12.5 or temp_y > 37.5):
									hist_charge_readout_timed.Fill(charge_l[i])
								else:
									hist_charge_floating_timed.Fill(charge_l[i])
								if (size_l[i] > 3):
									hist_charge_timed[3].Fill(charge_l[i])
									hist_profile_thesis_timed[3].Fill(temp_y)
								else:
									#print size_l[i]
									hist_charge_timed[size_l[i]].Fill(charge_l[i])
									hist_profile_thesis_timed[size_l[i]].Fill(temp_y)

						if (10 in y_t_loc.iterkeys()):
							hist_residual_timed.Fill(y_t_loc[10] - y_l_loc[10])
				hits_hits_on_lycoris.Fill(n_lycoris)
				x_l = {}
				y_l = {}
				z_l = {}
				t_l = {}
				SbyN_l = {}
				size_l = {}
				charge_l = {}
				ID_l = {}
				x_t = {}
				y_t = {}
				chi2_t = {}
				x_t_loc = {}
				y_t_loc = {}
				x_l_loc = {}
				y_l_loc = {}
				n_lycoris = 0
				n_timepix = 0
				t_lycoris = None
				continue
			else:
				#fields = line.split( )
				fields =  re.split(r',| |[(|)]|\[|\]|\r\n', line) ## splitting fields using regular expression
				fields = filter(None, fields)
				#print fields
				layer = int(fields[0])
				p = int(fields[1])
				x_layer = float(fields[2])
				y_layer = float(fields[3])
				z_layer = float(fields[4])
				t_layer = float(fields[5])
				#hist_hits_layers.Fill(layer)
				if layer >= 10 and layer < 20:
					SbyN_layer = float(fields[6])
					size_layer = int(fields[7])
					charge_layer = float(fields[8])
					ID_layer = float(fields[9])
					x_loc_track = float(fields[10])
					y_loc_track = float(fields[11])
					loc2Sens00 = float(fields[12])
					loc2Sens01 = float(fields[13])
					resid = float(fields[14])
					#print fields[15]
					x_loc_meas = float(fields[16])
					y_loc_meas = float(fields[17])
					hist_hits_layer_untimed.Fill(layer)
					hist_hits_layer_untimed2.Fill(layer)
					hist_charge[0].Fill(charge_layer)
					hist_csize.Fill(size_layer)
					hist_SbyN.Fill(SbyN_layer)
					subcell_y = y_loc_meas%0.05
					if (subcell_y < 0.0125) or (subcell_y > 0.0375):
						hist_charge_readout.Fill(charge_layer)
					elif (subcell_y >= 0.0125) and (subcell_y <= 0.0375):
						hist_charge_floating.Fill(charge_layer)

#					print fields
#					print fields[14]
#					break
					#chi2_track = float(fields[11])
					size_fitpos_correlation.Fill(size_layer,y_loc_track%0.05)
					if layer == 10:
						hist_residual.Fill(resid)
						hist_residual2.Fill(y_loc_track-y_loc_meas)
					hist_profile_thesis[0].Fill(y_loc_track%0.05)
					if (size_layer > 3):
						hist_charge[3].Fill(charge_layer)
						hist_profile_thesis[3].Fill(y_loc_track%0.05)
					else:
						hist_charge[size_layer].Fill(charge_layer)
						hist_profile_thesis[size_layer].Fill(y_loc_track%0.05)
					hist_y_track_local.Fill(y_loc_track)
					t_lycoris = t_layer
					n_lycoris += 1
				elif layer == 20:
					SbyN_layer = 0
					size_layer = 0
					charge_layer = 0
					ID_layer = 0
					x_loc_track = float(fields[7])
					y_loc_track = float(fields[8])
					loc2Sens00 = float(fields[9])
					loc2Sens01 = float(fields[10])
					resid = float(fields[11])
					x_loc_meas = float(fields[12])
					y_loc_meas = float(fields[13])
					#chi2_track = float(fields[8])
					n_timepix+=1
				else:
					SbyN_layer = 0
					size_layer = 0
					charge_layer = 0
					ID_layer = 0
					x_track = float(fields[8])
					y_track = float(fields[9])
					#chi2_track = float(fields[10])
				x_l[layer] = x_layer
				y_l[layer] = y_layer
				z_l[layer] = z_layer
				t_l[layer] = t_layer
				SbyN_l[layer] = SbyN_layer
				size_l[layer] = size_layer
				charge_l[layer] = charge_layer
				ID_l[layer] = ID_layer
				x_t_loc[layer] = x_loc_track
				y_t_loc[layer] = y_loc_track
				x_l_loc[layer] = x_loc_meas
				y_l_loc[layer] = y_loc_meas
				#x_t[layer] = x_track
				#y_t[layer] = y_track
				#chi2_t[layer] = chi2_track


#				timeDiff = timeK-timeT
#				if hitsK > 0 and hitsT > 0:
#					hist_timeDiff.Fill(timeDiff)
#					timeK_cycles = timeK % 80E6
#					hist_XY_tpx.Fill(x,y)
#					#print timeK_cycles
#					hist_KTime_v_TTime.Fill(timeK%80E6, timeT%80E6)
#					if abs(timeDiff) < 500:
#					hist_KPiX_hits_select1.Fill(hitsK)
#					for i in fields[15:]:
#						if int(i) >= 10 and int(i) < 20:
#							hist_hits_layer.Fill(int(i))
#	NrOfTracks = hist_KPiX_hits_select1.GetEntries()
	if (ROOT.TEfficiency.CheckConsistency(hist_hits_layer, hist_total)):
		hist_efficiency = ROOT.TEfficiency(hist_hits_layer, hist_total)
		for i in range(1,7,1):
			print "Efficiency ", hist_efficiency.GetEfficiency(i)
			print "Efficiency err up ",hist_efficiency.GetEfficiencyErrorUp(i)
			print "Efficiency err down ",hist_efficiency.GetEfficiencyErrorLow(i)

		hist_efficiency.Write()
	if (nr_tracks > 0):
		hist_hits_layer.Scale(1./nr_tracks)
	if (nr_tracks_3p > 0):
		hist_hits_layer_3p.Scale(1./nr_tracks_3p)
	if (nr_tracks_5p > 0):
		hist_hits_layer_5p.Scale(1./nr_tracks_5p)
	if (nr_tracks_untimed > 0):
		hist_hits_layer_untimed.Scale(1./nr_tracks_untimed)
	if (nr_tracks_all > 0):
		hist_hits_layer_untimed2.Scale(1./nr_tracks_all)

	#hist_hits_layer_3p.Scale(1./nr_tracks_3p)
	#hist_hits_layer_5p.Scale(1./nr_tracks_5p)
	outHistFile.Write()

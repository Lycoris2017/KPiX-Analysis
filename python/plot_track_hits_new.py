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

	corr_yt_yk = []
	corr_yt_yk_matched = []
	for j in range(10,16,1):
		histname = "correlation_track_layer"+str(j)
		corr_yt_yk.append(ROOT.TH2F(histname, "correlation; y_track(mm); y_kpix(mm); Nr. of Entries ", 920, -46, 46, 920, -46, 46))
		histname = "correlation_track_layer_time_matched_"+str(j)
		corr_yt_yk_matched.append(ROOT.TH2F(histname, "correlation; y_track(mm); y_kpix(mm); Nr. of Entries ", 920, -46, 46, 920, -46, 46))

	hist_hits_layer_3p = ROOT.TH1F("KPiX_efficiency_3p", "KPiX_efficiency_3p; Layer Number; Nr. of Entries/Nr. of Tracks ", 6, 9.5, 15.5)
	hist_hits_layer_5p = ROOT.TH1F("KPiX_efficiency_5p", "KPiX_efficiency_5p; Layer Number; Nr. of Entries/Nr. of Tracks ", 6, 9.5, 15.5)

	hist_hits_layer = ROOT.TH1F("KPiX_efficiency_all", "KPiX_efficiency_all; Layer Number; Nr. of Entries/Nr. of Tracks ", 6, 9.5, 15.5)
	hist_chi2_track = ROOT.TH1F("track_chi2", "track_chi2; chi2; Nr. of Entries", 200, 0, 50)

	hist_meancharge_track = ROOT.TH1F("track_mean_charge", "track_mean_charge; Charge (fC); Nr. of Entries", 200, 0, 20)
	hist_meanSbyN_track = ROOT.TH1F("track_mean_SbyN", "track_mean_SbyN; SbyN; Nr. of Entries", 200, 0, 20)

	hist_pos_track = ROOT.TH2F("track_pos", "track_pos; x (mm); y (mm); Nr. of Entries;", 200, -30, 30, 200, -30, 30)
	hist_pos_track_timed = ROOT.TH2F("track_pos_timed", "track_pos_timed; x (mm); y (mm); Nr. of Entries;", 200, -30, 30, 200, -30, 30)

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
	n_lycoris = 0
	n_timepix = 0
	t_lycoris = None
	nr_tracks = 0
	nr_tracks_3p = 0
	nr_tracks_5p = 0
	with open(args.file_in) as inFile:
		for line in inFile:
			if "run" in line or "Track" in line:
				if (5 in t_l.keys()):
					hist_pos_track.Fill(x_t[5], y_t[5])
					if (20 in t_l.keys()):
						hist_pos_track_timed.Fill(x_t[5], y_t[5])
				if (n_lycoris > 0 and 20 in t_l.keys()):
					t_diff = t_lycoris - t_l[20]
					if abs(t_diff) < 500:
						nr_tracks += 1

						for i in t_l.iterkeys():
							if i >= 10 and i < 20:
								hist_hits_layer.Fill(i)
								corr_yt_yk_matched[i%10].Fill(y_l[i], y_t[i])
						if (n_lycoris >= 3):
							nr_tracks_3p += 1
							for i in t_l.iterkeys():
								if i >= 10 and i < 20:
									hist_hits_layer_3p.Fill(i)
						if (n_lycoris >= 5):
							nr_tracks_5p += 1
							for i in t_l.iterkeys():
								if i >= 10 and i < 20:
									hist_hits_layer_5p.Fill(i)
				#hist_chi2_track.Fill(chi2_track)
				charge_sum = 0.
				SbyN_sum = 0.
				if len(charge_l ) > 0:
					for key, value in charge_l.iteritems():
						charge_sum += value
						SbyN_sum += SbyN_l[key]
					hist_meancharge_track.Fill(charge_sum/len(charge_l))
					hist_meanSbyN_track.Fill(SbyN_sum/len(charge_l))
				for  key,value in chi2_t.iteritems():
					hist_chi2_track.Fill(value)
					break


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
				n_lycoris = 0
				n_timepix = 0
				t_lycoris = None
				continue
			else:
				#fields = line.split( )
				fields =  re.split(r',| |[(|)]|\[|\]|\r\n', line) ## splitting fields using regular expression
				fields = filter(None, fields)
				layer = int(fields[0])
				x_layer = float(fields[1])
				y_layer = float(fields[2])
				z_layer = float(fields[3])
				t_layer = float(fields[4])
				if layer >= 10 and layer < 20:
					SbyN_layer = float(fields[5])
					size_layer = float(fields[6])
					charge_layer = float(fields[7])
					ID_layer = float(fields[8])
					x_track = float(fields[9])
					y_track = float(fields[10])
					#chi2_track = float(fields[11])
					t_lycoris = t_layer
					n_lycoris += 1
					corr_yt_yk[layer%10].Fill(y_track, y_layer)
				elif layer == 20:
					SbyN_layer = 0
					size_layer = 0
					charge_layer = 0
					ID_layer = 0
					x_track = float(fields[6])
					y_track = float(fields[7])
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
				x_t[layer] = x_track
				y_t[layer] = y_track

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
	hist_hits_layer.Scale(1./nr_tracks)
	hist_hits_layer_3p.Scale(1./nr_tracks_3p)
	hist_hits_layer_5p.Scale(1./nr_tracks_5p)
	outHistFile.Write()

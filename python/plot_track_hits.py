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
if len(sys.argv) < 2:
	print parser.print_help()
	sys.exit(1)
print ''
outfile = args.file_in[:-4]+'.root'
outHistFile = ROOT.TFile.Open(outfile, "RECREATE")
outHistFile.cd()

hist_x = ROOT.TH1F("x(0)", "x(0); x(mm); Nr. of Entries ", 1840, -46, 46)
hist_y = ROOT.TH1F("y(0)", "y(0); y(mm); Nr. of Entries ", 1840, -46, 46)
hist_z1 = ROOT.TH1F("z(0)", "z(0); z(mm); Nr. of Entries ", 101, 50, 50)
hist_z2 = ROOT.TH1F("z(-1)", "z(-1); z(mm); Nr. of Entries ", 101, 50, 50)

hist_dx = ROOT.TH1F("dx", "dx; dx(rad); Nr. of Entries ", 1001, -0.1, 0.1)
hist_dy = ROOT.TH1F("dy", "y; dy(rad); Nr. ofEntries ", 1001, -0.02, 0.02)

hist_SbyN = ROOT.TH1F("SbyN", "SbyN; SbyN; Nr. of Entries ", 1000, 0, 100)

hist_HitsOnTrack = ROOT.TH1F("HitsOnTrack", "HitOnTrack; Nr. of Hits; Nr. of Entries ", 14, -0.5, 13.5)

hist_x_v_tracks = ROOT.TH2F("x_v_tracks", "x_v_tracks; x (mm); Nr. of tracks; Nr. of Entries",  1840,-46,46, 9,0.5,9.5)

hist_KPiX_hits = ROOT.TH1F("KPiXHits", "KPiXHits; Nr. of Hits; Nr. of Entries ", 7, -0.5, 6.5)
hist_Timepix_hits = ROOT.TH1F("TimepixHits", "TimepixHits; Nr. of Hits; Nr. of Entries ", 7, -0.5, 6.5)

hist_timeDiff = ROOT.TH1F("timeDiff",
                       "timeDiff; time (ns); #Entries ",
                       2001,
                       #-50e-2, 50e-2
                       -1000, 1000)
hist_KPiX_hits_select1 = ROOT.TH1F("KPiXHits_select", "KPiXHits_select; Nr. of Hits; Nr. of Entries ", 7, -0.5, 6.5)


hist_KTime_v_TTime = ROOT.TH2D("KTime_v_TTime", "KTime_v_TTime; time_kpix (us); time_timepix (us); Nr. of Entries",  1000,0,10000, 10000,0,10000)

hist_hits_layer = ROOT.TH1F("KPiX_efficiency", "KPiX_efficiency; Layer Number; Nr. of Entries/Nr. of Tracks ", 6, 9.5, 15.5)
hist_XY_all = ROOT.TH2F("x_v_y_all", "x_v_y_all; x (mm); y(mm); Nr. of Entries",  184,-46,46, 184,-46,46)
hist_XY_tpx = ROOT.TH2F("x_v_y_tpx", "x_v_y_tpx; x (mm); y(mm); Nr. of Entries",  184,-46,46, 184,-46,46)
#
#x0=array('f', [0])
#y0=array('f', [0])
#z1=array('f', [0])
#z2=array('f', [0])
#dx=array('f', [0])
#dy=array('f', [0])
#z1=array('f', [0])
#z2=array('f', [0])
#ccharge=array('f', [0])
#significance=array('f', [0])
#noise=array('f', [0])
#csize=array('i', [0])
#yPos=array('f', [0])
#onTrack=array('i', [0])

#corrClusterTree = ROOT.TTree("segments", "segments");
#corrClusterTree.Branch("x0", x0, 'x0/I');
#corrClusterTree.Branch("y0", ccharge, 'charge/F');
#corrClusterTree.Branch("significance", significance, 'significance/F');
#corrClusterTree.Branch("noise", noise, 'noise/F');
#corrClusterTree.Branch("csize", csize, 'csize/I');
#corrClusterTree.Branch("yPos", yPos, 'yPos/F');
#corrClusterTree.Branch("onTrack", onTrack, 'onTrack/I');


count = 0
matchCounter = 0
tracks = 0
_x = []
with open(args.file_in) as inFile:
	for line in inFile:
		if "run" in line or "Hits" in line:
			for i in _x:
				hist_x_v_tracks.Fill( i,len(_x))
			_x = []
			continue
		#print line
		fields = line.split( )
		hits = int(fields[0])
		x = float(fields[1])
		_x.append(x)
		y = float(fields[2])
		z1 = float(fields[3])
		z2 = float(fields[4])
		dx = float(fields[5])
		dy = float(fields[6])
		SbyN = float(fields[14])
		timeT = float(fields[13])
		hitsT = int(fields[12])
		timeK = float(fields[11])
		hitsK = float(fields[10])

		hist_x.Fill(x)
		hist_y.Fill(y)
		hist_z1.Fill(z1)
		hist_z2.Fill(z2)
		hist_dx.Fill(dx)
		hist_dy.Fill(dy)
		hist_HitsOnTrack.Fill(hits)
		hist_SbyN.Fill(SbyN)
		hist_KPiX_hits.Fill(hitsK)
		hist_Timepix_hits.Fill(hitsT)
		#hist_timepixTime.Fill(timeT)
		#hist_kpixTime.Fill(timeK)
		hist_XY_all.Fill(x,y)

		timeDiff = timeK-timeT
		if hitsK > 0 and hitsT > 0:
			hist_timeDiff.Fill(timeDiff)
			timeK_cycles = timeK % 80E6
			hist_XY_tpx.Fill(x,y)
			#print timeK_cycles
			hist_KTime_v_TTime.Fill(timeK%80E6, timeT%80E6)
			if abs(timeDiff) < 500:
				hist_KPiX_hits_select1.Fill(hitsK)
				for i in fields[15:]:
					if int(i) >= 10 and int(i) < 20:
						hist_hits_layer.Fill(int(i))

		if hits >= 7:
			matchCounter +=1
NrOfTracks = hist_KPiX_hits_select1.GetEntries()
hist_hits_layer.Scale(1./NrOfTracks)
#c1 = ROOT.TCanvas( 'test', 'Test', 1600, 900 )
#c1.cd()
print matchCounter
outHistFile.Write()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

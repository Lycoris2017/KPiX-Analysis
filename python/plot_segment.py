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

hist_HitsOnTrack = ROOT.TH1F("HitsOnTrack", "HitOnTrack; Nr. of Hits; Nr. of Entries ", 13, -0.5, 12.5)

hist_xColl = ROOT.TH1F("x(coll)", "x(coll); x(mm); #Entries ", 1840, -46, 46)
hist_yColl = ROOT.TH1F("y(coll)", "y(coll); y(mm); #Entries ", 1840, -46, 46)

hist_yColl_v_dy = ROOT.TH2F("yColl_v_dy", "yColl_v_dy; yColl (mm); dy (rad); Nr. of Entries", 1840,-46,46, 201, -0.02, 0.02)

hist_x_v_tracks = ROOT.TH2F("x_v_tracks", "x_v_tracks; x (mm); Nr. of tracks; Nr. of Entries",  1840,-46,46, 9,0.5,9.5)


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
		SbyN = float(fields[9])
		x_coll = float(fields[10])
		y_coll = float(fields[11])
		z_coll = float(fields[12])
		#print fields

		hist_x.Fill(x)
		hist_y.Fill(y)
		hist_z1.Fill(z1)
		hist_z2.Fill(z2)
		hist_dx.Fill(dx)
		hist_dy.Fill(dy)
		hist_HitsOnTrack.Fill(hits)
		hist_SbyN.Fill(SbyN)
		hist_xColl.Fill(x_coll)
		hist_yColl.Fill(y_coll)

		hist_yColl_v_dy.Fill(y_coll,dy)
		if hits >= 7:
			matchCounter +=1


#c1 = ROOT.TCanvas( 'test', 'Test', 1600, 900 )
#c1.cd()
print matchCounter
outHistFile.Write()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

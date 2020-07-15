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

hist_x = ROOT.TH1F("x", "x; x(mm); #Entries ", 1840, -46, 46)
hist_y = ROOT.TH1F("y", "y; y(mm); #Entries ", 1840, -46, 46)
hist_z1 = ROOT.TH1F("z1", "z1; z(mm); #Entries ", 101, 50, 50)
hist_z2 = ROOT.TH1F("z2", "z2; z(mm); #Entries ", 101, 50, 50)

hist_dx = ROOT.TH1F("dx", "dx; dx(rad); #Entries ", 1001, -0.1, 0.1)
hist_dy = ROOT.TH1F("dy", "y; dy(rad); #Entries ", 1001, -0.1, 0.1)

hist_SbyN = ROOT.TH1F("SbyN", "SbyN; SbyN; #Entries ", 1000, 0, 100)

hist_HitsOnTrack = ROOT.TH1F("HitsOnTrack", "HitOnTrack; Nr. of Hits; #Entries ", 7, 0, 7)

count = 0

with open(args.file_in) as inFile:
	for line in inFile:
		if "run" in line or "Hits" in line:
			continue
		#print line
		fields = line.split( )
		hits = int(fields[0])
		x = int(fields[1])
		y = int(fields[2])
		z1 = int(fields[3])
		z2 = int(fields[4])
		dx = int(fields[5])
		dy = int(fields[6])

		hist_x.Fill(x)
		hist_y.Fill(y)
		hist_z1.Fill(z1)
		hist_z2.Fill(z2)
		hist_dx.Fill(dx)
		hist_dy.Fill(dy)
		hist_HitsOnTrack.Fill(hits)


#c1 = ROOT.TCanvas( 'test', 'Test', 1600, 900 )
#c1.cd()

outHistFile.Write()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

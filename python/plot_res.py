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
parser.add_argument('file_out', help='name of the output file')
args = parser.parse_args()
if len(sys.argv) < 2:
	print parser.print_help()
	sys.exit(1)
print ''

outHistFile = ROOT.TFile.Open(args.file_out, "RECREATE")
outHistFile.cd()
residual_hist = []

residual_hist.append(ROOT.TH1F("residual_l10", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist.append(ROOT.TH1F("residual_l11", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist.append(ROOT.TH1F("residual_l12", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist.append(ROOT.TH1F("residual_l13", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist.append(ROOT.TH1F("residual_l14", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist.append(ROOT.TH1F("residual_l15", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))

residual_hist2 = []

residual_hist2.append(ROOT.TH1F("residual2_l10", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist2.append(ROOT.TH1F("residual2_l11", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist2.append(ROOT.TH1F("residual2_l12", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist2.append(ROOT.TH1F("residual2_l13", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist2.append(ROOT.TH1F("residual2_l14", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist2.append(ROOT.TH1F("residual2_l15", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))

count = 0

with open(args.file_in) as inFile:
	line = inFile.readline()
	if (count == 0):
		line = inFile.readline()
	count += 1
	while line:
		#print line
		splitline = line.split( )
		if (abs(float(splitline[5])) > 0.00005):
			if int(splitline[1]) is 10:
				residual_hist[0].Fill(float(splitline[5]))
			if int(splitline[1]) is 11:
				residual_hist[1].Fill(float(splitline[5]))
			if int(splitline[1]) is 12:
				residual_hist[2].Fill(float(splitline[5]))
			if int(splitline[1]) is 13:
				residual_hist[3].Fill(float(splitline[5]))
			if int(splitline[1]) is 14:
				residual_hist[4].Fill(float(splitline[5]))
			if int(splitline[1]) is 15:
				residual_hist[5].Fill(float(splitline[5]))
		if (abs(float(splitline[5])) > 0.00005):
			if int(splitline[1]) is 10:
				residual_hist2[0].Fill(float(splitline[18]))
			if int(splitline[1]) is 11:
				residual_hist2[1].Fill(float(splitline[18]))
			if int(splitline[1]) is 12:
				residual_hist2[2].Fill(float(splitline[18]))
			if int(splitline[1]) is 13:
				residual_hist2[3].Fill(float(splitline[18]))
			if int(splitline[1]) is 14:
				residual_hist2[4].Fill(float(splitline[18]))
			if int(splitline[1]) is 15:
				residual_hist2[5].Fill(float(splitline[18]))
		line = inFile.readline()

#c1 = ROOT.TCanvas( 'test', 'Test', 1600, 900 )
#c1.cd()
for i in residual_hist:
	i.Fit("gaus")
for i in residual_hist2:
	i.Fit("gaus")
outHistFile.Write()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

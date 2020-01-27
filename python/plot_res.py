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

outHistFile = ROOT.TFile.Open("res_results.root", "RECREATE")
outHistFile.cd()
residual_hist = []

residual_hist.append(ROOT.TH1F("residual_s0", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist.append(ROOT.TH1F("residual_s1", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))
residual_hist.append(ROOT.TH1F("residual_s2", "residual; distance (mm); #Entries ", 101, -50e-3, 50e-3))



with open(args.file_in) as inFile:
	line = inFile.readline()
	while line:
#		print line
		splitline = line.split( )
		if int(splitline[1]) is 13:
			residual_hist[0].Fill(float(splitline[4]))
		if int(splitline[1]) is 14:
			residual_hist[1].Fill(float(splitline[4]))
		if int(splitline[1]) is 15:
			residual_hist[2].Fill(float(splitline[4]))
		line = inFile.readline()

#c1 = ROOT.TCanvas( 'test', 'Test', 1600, 900 )
#c1.cd()
for i in residual_hist:
	i.Fit("gaus")
outHistFile.Write()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

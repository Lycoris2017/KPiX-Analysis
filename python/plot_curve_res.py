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

outHistFile = ROOT.TFile.Open("momentum_res_results.root", "RECREATE")
outHistFile.cd()
residual_hist = []

residual_hist.append(ROOT.TH1F("momentum_residual", "residual; dQbyP; #Entries ", 50, -50e-3, 50e-3))




with open(args.file_in) as inFile:
	line = inFile.readline()
	while line:
#		print line
		splitline = line.split( )
		residual_hist[0].Fill(float(splitline[2]))
		line = inFile.readline()

#c1 = ROOT.TCanvas( 'test', 'Test', 1600, 900 )
#c1.cd()
for i in residual_hist:
	i.Fit("gaus", "R", "", -0.01, 0.010)
outHistFile.Write()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

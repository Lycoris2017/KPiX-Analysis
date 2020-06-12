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
residual_hist = ROOT.TH1F("mom_residual", "mom_residual; momentum (GeV/c); #Entries ", 101, -50e-2, 50e-2)
loccov_hist = ROOT.TH1F("loccov", "loccov; A.U.; #Entries ", 101, -50e-3, 50e-3)

count = 0

with open(args.file_in) as inFile:
	line = inFile.readline()
	if (count == 0):
		line = inFile.readline()
	count += 1
	while line:
		print line
		splitline = line.split( )
		#if (abs(float(splitline[5])) > 0.00005):
		residual_hist.Fill(float(splitline[1])+float(splitline[2]))
		loccov_hist.Fill(float(splitline[3]))
		line = inFile.readline()


residual_hist.Fit("gaus")
outHistFile.Write()

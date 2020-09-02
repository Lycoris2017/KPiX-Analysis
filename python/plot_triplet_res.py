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
if len(sys.argv) < 2:
	print parser.print_help()
	sys.exit(1)
print ''
outfile = args.file_in[:-4]+'.root'
outHistFile = ROOT.TFile.Open(outfile, "RECREATE")
outHistFile.cd()


residual_x_total = ROOT.TH1F("residual_x_total",
                       "residual; distance (mm); #Entries ",
                       101,
                       #-50e-2, 50e-2
                       -1.0, 1.0)

residual_x_0 = ROOT.TH1F("residual_x_0",
                       "residual; distance (mm); #Entries ",
                       101,
                       #-50e-2, 50e-2
                       -1.0, 1.0)

residual_x_1 = ROOT.TH1F("residual_x_1",
                        "residual; distance (mm); #Entries ",
                        101,
                        #-50e-2, 50e-2
                        -1.0, 1.0)


residual_y_total = ROOT.TH1F("residual_y_total",
                        "residual; distance (mm); #Entries ",
                        101,
                        #-50e-2, 50e-2
                        -0.1, 0.1)

residual_y_0 = ROOT.TH1F("residual_y_0",
                        "residual; distance (mm); #Entries ",
                        101,
                        #-50e-2, 50e-2
                        -0.1, 0.1)

residual_y_1 = ROOT.TH1F("residual_y_1",
                        "residual; distance (mm); #Entries ",
                        101,
                        #-50e-2, 50e-2
                        -0.1, 0.1)


with open(args.file_in) as inFile:
	for line in inFile:
		fields =  re.split('=| |\r\n', line) ## splitting fields using regular expression
		fields = filter(None, fields) ## getting rid of empty fields

		if ('CassetteNum' in line):
			print fields
			continue
		dx = float(fields[3]) #residual x
		dy = float(fields[4])
		cx = float(fields[9]) #gbl correction x
		cy = float(fields[10])
		residual_x_total.Fill(dx-cx)
		residual_y_total.Fill(dy-cy)
		if (float(fields[11]) < 0):
			residual_x_0.Fill(dx-cx)
			residual_y_0.Fill(dy-cy)
		elif (float(fields[11]) > 0):
			residual_x_1.Fill(dx-cx)
			residual_y_1.Fill(dy-cy)
residual_x_total.Fit("gaus")
residual_x_0.Fit("gaus")
residual_x_1.Fit("gaus")

residual_y_total.Fit("gaus")
residual_y_0.Fit("gaus")
residual_y_1.Fit("gaus")

outHistFile.Write()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

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

hits_on_track= ROOT.TH1F("hits_on_track", "track hits; #hits; #Entries ", 6, +0.5, 6.5)



with open(args.file_in) as inFile:
	line = inFile.readline()
	line = inFile.readline()
	while line:
#		print line
		splitline = line.split( )
		print splitline[1]
		hits_on_track.Fill(int(splitline[1]))
		line = inFile.readline()

#c1 = ROOT.TCanvas( 'test', 'Test', 1600, 900 )
#c1.cd()
outHistFile.Write()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

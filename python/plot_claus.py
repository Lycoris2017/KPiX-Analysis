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
charge_hist = ROOT.TH1F("Charge", "Cluster Charge; Charge (fC); Number of Entries ", 200, 0, 20)
charge_hist2 = ROOT.TH1F("Charge_w_pos_cut", "Cluster Charge; Charge (fC); Number of Entries ", 200, 0, 20)
significance_hist = ROOT.TH1F("Significance", "Cluster Significance; S/N; Number of Entries ", 500, 0, 50)
significance_hist2 = ROOT.TH1F("Significance_w_pos_cut", "Cluster Significance; S/N; Number of Entries ", 500, 0, 50)
y_hist = ROOT.TH1F("y", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46)
x_hist = ROOT.TH1F("x", "Cluster x; position (mm); Number of Entries ", 1840, -46, 46)
corr_hist = ROOT.TH1F("Correlations", "Correlations; num of correlated planes; Number of Hits ", 10, 0, 9)
correlation_s0_s2 = ROOT.TH2F("correlation", "correlation y; Sensor 0 (mm); Sensor 2 (mm); Number of Hits", 1840, -46, 46, 1840, -46, 46)

trackhits_hist = ROOT.TH1F("Track hits", "Track hits; #Hits on track; Number of Entries ", 10, -0.5, 9.5)

significance_hists =  []
significance_hists.append(ROOT.TH1F("Significance_s0", "Cluster Significance; S/N; Number of Entries ", 200, 0, 100))
significance_hists.append(ROOT.TH1F("Significance_s1", "Cluster Significance; S/N; Number of Entries ", 200, 0, 100))
significance_hists.append(ROOT.TH1F("Significance_s2", "Cluster Significance; S/N; Number of Entries ", 200, 0, 100))

charge_hists = []
charge_hists.append(ROOT.TH1F("Charge_s0", "Cluster Charge; Charge (fC); Entries ", 200, 0, 20))
charge_hists.append(ROOT.TH1F("Charge_s1", "Cluster Charge; Charge (fC); Entries ", 200, 0, 20))
charge_hists.append(ROOT.TH1F("Charge_s2", "Cluster Charge; Charge (fC); Entries ", 200, 0, 20))

y_hists = []
y_hists.append(ROOT.TH1F("y_s0", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46))
y_hists.append(ROOT.TH1F("y_s1", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46))
y_hists.append(ROOT.TH1F("y_s2", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46))

x_hists = []
x_hists.append(ROOT.TH1F("x_s0", "Cluster x; position (mm); Number of Entries ", 1840, -1, 1))
x_hists.append(ROOT.TH1F("x_s1", "Cluster x; position (mm); Number of Entries ", 1840, -1, 1))
x_hists.append(ROOT.TH1F("x_s2", "Cluster x; position (mm); Number of Entries ", 1840, -1, 1))

size_hists = []
size_hists.append(ROOT.TH1F("size_s0", "size; #strips; Number of Entries ", 10, -0.5, 9.5))
size_hists.append(ROOT.TH1F("size_s1", "size; #strips; Number of Entries ", 10, -0.5, 9.5))
size_hists.append(ROOT.TH1F("size_s2", "size; #strips; Number of Entries ", 10, -0.5, 9.5))

all_hists = [charge_hist, charge_hist2, significance_hist, significance_hist2, y_hist, x_hist, corr_hist, trackhits_hist, charge_hists, y_hists[0:3], x_hists[0:3], size_hists[0:3], significance_hists ]


print all_hists
plane = []
x = []
y = []
z = []
corr = []
sig = []
size = []
charge = []

y_s0 = []
y_s2 = []

with open(args.file_in) as inFile:
	line = inFile.readline()
	hits_on_track = 0

	while line:
		if "run" in line:
			if (len(y_s0) is not 0  and len(y_s2) is not 0):
				for q in y_s0:
					for h in y_s2:
						correlation_s0_s2.Fill(float(q),float(h))
				y_s0 = []
				y_s2 = []
			trackhits_hist.Fill(hits_on_track)
			hits_on_track=0
			line = inFile.readline()
		else:
			a = line.split( )
			plane.append(a[0])
			x.append(a[1])
			y.append(a[2])
			z.append(a[3])
			corr.append(a[4])
			sig.append(a[5])
			size.append(a[6])
			charge.append(a[7])
			if "True" in a[8]:
				hits_on_track += 1
			if "10" in a[0]:
				y_s0.append(a[2])
			elif "12" in a[0]:
				y_s2.append(a[2])
			line = inFile.readline()

#c1 = ROOT.TCanvas( 'test', 'Test', 1600, 900 )
#c1.cd()

for q in xrange(1,len(charge)):
	corr_hist.Fill(int(corr[q]))
	if (int(corr[q]) >= 1):
		if ("10" in plane[q]):
			charge_hists[0].Fill(float(charge[q]))
			significance_hists[0].Fill(float(sig[q]))
			y_hists[0].Fill(float(y[q]))
			x_hists[0].Fill(float(x[q]))
			size_hists[0].Fill(int(size[q]))
		elif ("11" in plane[q]):
			charge_hists[1].Fill(float(charge[q]))
			significance_hists[1].Fill(float(sig[q]))
			y_hists[1].Fill(float(y[q]))
			x_hists[1].Fill(float(x[q]))
			size_hists[1].Fill(int(size[q]))
		elif ("12" in plane[q]):
			charge_hists[2].Fill(float(charge[q]))
			significance_hists[2].Fill(float(sig[q]))
			y_hists[2].Fill(float(y[q]))
			x_hists[2].Fill(float(x[q]))
			size_hists[2].Fill(int(size[q]))
		charge_hist.Fill(float(charge[q]))
		significance_hist.Fill(float(sig[q]))
		y_hist.Fill(float(y[q]))
		x_hist.Fill(float(x[q]))
		if (float(y[q]) > - 35 and float(y[q]) < -15):
			charge_hist2.Fill(float(charge[q]))
			significance_hist2.Fill(float(sig[q]))

for i in all_hists:
	if isinstance(i, list):
		for j in i:
			j.Scale(1.0/j.GetEntries())
	else:
		i.Scale(1.0/i.GetEntries())

outHistFile.Write()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

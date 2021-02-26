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
outfile = args.file_in[:-4]+'_redone.root'
outHistFile = ROOT.TFile.Open(outfile, "RECREATE")
outHistFile.cd()
residual_hist = []

residual_hist.append(ROOT.TH1F("residual_l10", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))
residual_hist.append(ROOT.TH1F("residual_l11", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))
residual_hist.append(ROOT.TH1F("residual_l12", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))
residual_hist.append(ROOT.TH1F("residual_l13", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))
residual_hist.append(ROOT.TH1F("residual_l14", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))
residual_hist.append(ROOT.TH1F("residual_l15", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))

residual_hist2 = []

residual_hist2.append(ROOT.TH1F("residual2_l10", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))
residual_hist2.append(ROOT.TH1F("residual2_l11", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))
residual_hist2.append(ROOT.TH1F("residual2_l12", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))
residual_hist2.append(ROOT.TH1F("residual2_l13", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))
residual_hist2.append(ROOT.TH1F("residual2_l14", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))
residual_hist2.append(ROOT.TH1F("residual2_l15", "residual; Distance (#mum); No. of Entries ", 101, -50, 50))



fitpos_mod_readout_pitch = ROOT.TH1F("fitpos_mod_readout_pitch", "fitpos_mod_readout_pitch; distance (mm); Entries ", 50, 0, 50e-3)
fitpos_mod_readout_pitch_bump = ROOT.TH1F("fitpos_mod_readout_pitch_bump", "fitpos_mod_readout_pitch_bump; distance (mm); Entries ", 50, 0, 50e-3)

size_fitpos_correlation = ROOT.TH2F("size_fitpos", "size_fitpos; size; distance (mm); Number of Entries ", 5, 0.5, 5.5, 51, 0, 51e-3)
size_fitpos_correlation_bump = ROOT.TH2F("size_fitpos_bump", "size_fitpos_bump; size; distance (mm); Number of Entries ", 5, 0.5, 5.5, 51, 0, 51e-3)

charge_all = ROOT.TH1F("charge_all", "charge_all; Charge (fC); Entries ", 151, 0, 15)
charge_readout = ROOT.TH1F("charge_readout", "charge_readout; Charge (fC); Entries ", 151, 0, 15)
charge_floating = ROOT.TH1F("charge_floating", "charge_floating; Charge (fC); Entries ", 151, 0, 15)
charge_transition = ROOT.TH1F("charge_transition", "charge_transition; Charge (fC); Entries ", 151, 0, 15)

charge_size1 = ROOT.TH1F("charge_size1", "charge_size1; Charge (fC); Entries ", 151, 0, 15)
charge_size2 = ROOT.TH1F("charge_size2", "charge_size2; Charge (fC); Entries ", 151, 0, 15)
charge_size3p = ROOT.TH1F("charge_size3+", "charge_size3+; Charge (fC); Entries ", 151, 0, 15)

count = 0

with open(args.file_in) as inFile:
	for line in inFile:
		if "run" in line or "label" in line:
			continue
		#print line
		splitline = line.split( )
		layer = int(splitline[1])
		res1 = float(splitline[5])*1000
		res2 = float(splitline[20])*1000
		csize = int(splitline[4])
		charge = float(splitline[19])
		posMod50 = float(splitline[21])
		#if (abs(res1) > 0.00005):
		if layer is 10:
			residual_hist[0].Fill(res1)
		if layer is 11:
			residual_hist[1].Fill(res1)
		if layer is 12:
			residual_hist[2].Fill(res1)
		if layer is 13:
			residual_hist[3].Fill(res1)
		if layer is 14:
			residual_hist[4].Fill(res1)
		if layer is 15:
			residual_hist[5].Fill(res1)
		#if (abs(res1) > 0.00005):
		if layer is 10:
			residual_hist2[0].Fill(res2)
		if layer is 11:
			residual_hist2[1].Fill(res2)
		if layer is 12:
			residual_hist2[2].Fill(res2)
		if layer is 13:
			residual_hist2[3].Fill(res2)
		if layer is 14:
			residual_hist2[4].Fill(res2)
		if layer is 15:
			residual_hist2[5].Fill(res2)
		fitpos_mod_readout_pitch.Fill(posMod50)
		size_fitpos_correlation.Fill(csize,posMod50)

		if (charge < 1.5):
			size_fitpos_correlation_bump.Fill(csize,posMod50)
			fitpos_mod_readout_pitch_bump.Fill(posMod50)
		charge_all.Fill(charge)
		if (posMod50 > 0.0125 and posMod50 < 0.0375):
			charge_floating.Fill(charge)
		else:
			charge_readout.Fill(charge)
		if (posMod50 > 0.07 and posMod50 < 0.017) or (posMod50 > 0.032 and posMod50 < 0.042):
			charge_transition.Fill(charge)
		if (csize == 1):
			charge_size1.Fill(charge)
		elif (csize == 2):
			charge_size2.Fill(charge)
		else:
			charge_size3p.Fill(charge)

#c1 = ROOT.TCanvas( 'test', 'Test', 1600, 900 )
#c1.cd()
for i in residual_hist:
	fitGaus =  i.Fit("gaus", "RE", "", -10, 10)
for i in residual_hist2:
	fitGaus =  i.Fit("gaus", "RE", "", -12, 12)
outHistFile.Write()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

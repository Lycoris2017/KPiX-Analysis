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
from array import array
import re

ROOT.gROOT.SetBatch(True)


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)



parser = MyParser()
parser.add_argument(
    '-r', '--root',
    dest='root_in',
    help='name of the input files'
)
parser.add_argument(
    '-g', '--gbl',
    dest='GBL_in',
    help='name of the input files'
)
parser.add_argument(
    '-o', '--output',
    dest='file_out',
    help='specifies the name of the output file'
)
args = parser.parse_args()
if len(sys.argv) < 2:
    print parser.print_help()
    sys.exit(1)
print ''

ROOT.TH1.SetDefaultSumw2()

#//////////////////////////////////////////
#// Read GBL file
#//////////////////////////////////////////
print 'Reading GBL File'

fitpos = {}
with open(args.GBL_in) as GBLFile:
    for line in GBLFile:
        if ('label' in line) or ('run' in line):
            continue
        #print line
        fields =  re.split(r',| |[(|)]|\[|\]|\r\n', line) ## splitting fields using regular expression
        fields = filter(None, fields)
        if ('New' not in fields):
            if (int(fields[0]) >= 10)  and (int(fields[0]) <= 15):
                #print fields[9]
                #print fields[11]
                fitpos[int(fields[9])] = float(fields[11])


#print IDonTrack
#//////////////////////////////////////////
#// Read ROOT Tree file
#//////////////////////////////////////////

outHistFile = ROOT.TFile.Open(args.file_out, "RECREATE")
outHistFile.cd()

layer=array('i', [0])
charge=array('f', [0])
signi=array('f', [0])
pos=array('f', [0])
ID=array('i', [0])

stripOnTrackTree = ROOT.TTree("stripOnTrack", "stripOnTrackTree");
stripOnTrackTree.Branch("layer", layer, 'layer/I')
stripOnTrackTree.Branch("charge", charge, 'charge/F')
stripOnTrackTree.Branch("signi", signi, 'signi/F')
stripOnTrackTree.Branch("pos", pos, 'pos/F')
stripOnTrackTree.Branch("ID", ID, 'ID/I')

eta_hists = []
eta_hists.append(ROOT.TH1F("Eta_total", "Eta_total; Eta; Number of Entries", 30, -0.5,1.5))
eta_hists.append(ROOT.TH1F("Eta_1", "Eta_1; Eta; Number of Entries", 30, -0.5,1.5))
eta_hists.append(ROOT.TH1F("Eta_2", "Eta_2; Eta; Number of Entries", 30, -0.5,1.5))
eta_hists.append(ROOT.TH1F("Eta_3+", "Eta_3+; Eta; Number of Entries", 30, -0.5,1.5))

eta_hists_v_pos = ROOT.TH2F("eta_pos", "eta_pos; #eta; sub-cell position (mm)", 30,-0.5,1.5, 51, 0, 51e-3)

skewed_hists = []
skewed_hists.append(ROOT.TH1F("skewed_total", "skewed_total; Eta; Number of Entries", 30, -0.5,1.5))
skewed_hists.append(ROOT.TH1F("skewed_1", "skewed_1; Eta; Number of Entries", 30, -0.5,1.5))
skewed_hists.append(ROOT.TH1F("skewed_2", "skewed_2; Eta; Number of Entries", 30, -0.5,1.5))
skewed_hists.append(ROOT.TH1F("skewed_3+", "skewed_3+; Eta; Number of Entries", 30, -0.5,1.5))


print 'Reading ROOT File'

f = ROOT.TFile(args.root_in)
stripTree = f.Get("vector_sub_cluster")

print 'Performing analysis'
charge_right = 0
charge_left = 0
count = 1
for entry in stripTree:
    for i in xrange(len(entry.sensor)):
        if fitpos.get(entry.ID[i], -99) != -99:
            layer[0] = entry.sensor[i]
            charge[0] = entry.charge[i]
            signi[0] = entry.signi[i]
            pos[0] = entry.pos[i]/1000
            if (ID[0] == entry.ID[i]):
                count = count + 1
                if (pos[0]-fitpos.get(entry.ID[i]) > 0):
                    charge_right = charge_right + charge[0]
                else:
                    charge_left = charge_left + charge[0]
            else:
                if (pos[0]-fitpos.get(entry.ID[i]) > 0):
                    charge_right = charge_right + charge[0]
                else:
                    charge_left = charge_left + charge[0]
                eta = charge_right/(charge_left + charge_right)
                eta_hists[0].Fill(eta)
                eta_hists_v_pos.Fill(eta, fitpos.get(entry.ID[i])%50e-3)
                #print "Subcell pos ", fitpos.get(entry.ID[i])
                #print "Eta ", eta
                if (count == 1):
                    eta_hists[1].Fill(eta)
                elif (count == 2):
                    eta_hists[2].Fill(eta)
                else:
                    eta_hists[3].Fill(eta)
                charge_left = 0
                charge_right = 0
                count = 1
            ID[0] = entry.ID[i]
            stripOnTrackTree.Fill()


#for j in eta_hists:
#    j.Scale(1.0/j.GetEntries())




#for event in f.vectr_sub_cluster:
#    print event.sensor

#outHistFile.cd()



outHistFile.Write()
outHistFile.Close()
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

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
parser.add_argument(
    '-t', '--timepix',
    dest='timepix',
    default=False,
    action='store_true',
    help='only uses hits that are time matched with TPX3 when called'
)
args = parser.parse_args()
if len(sys.argv) < 2:
    print(parser.print_help())
    sys.exit(1)
print('')

ROOT.TH1.SetDefaultSumw2()

#//////////////////////////////////////////
#// Read GBL file
#//////////////////////////////////////////
print('Reading GBL File')

fitpos = {}
t_layer = {}
pos_layer = {}
ID_layer = {}
layerNum = []
with open(args.GBL_in) as GBLFile:
    for line in GBLFile:
        if ('label' in line) or ('run' in line) or "Track" in line:
            if args.timepix:
                if len(t_layer) > 1 and 20 in layerNum: #need to have entries in the track and must have timepix entry
                    for i in range(10,16,1):
                        if i in layerNum:
                            t_lycoris = t_layer[i] #first lycoris layer that exists will be time recorded.
                            break
                    t_diff = t_lycoris - t_layer[20] #Time matching between timepix and lycoris
                    if abs(t_diff) < 500:
                        for i in range(10,16,1):
                            if i in layerNum:
                                fitpos[ID_layer[i]] = pos_layer[i]
            else:
                for i in range(10,16,1):
                    if i in layerNum:
                        fitpos[ID_layer[i]] = pos_layer[i]
            t_layer = {}
            pos_layer = {}
            ID_layer = {}
            layerNum = []
            continue
        #print line
        fields =  re.split(r',| |[(|)]|\[|\]|\r\n', line) ## splitting fields using regular expression
        fields = list(filter(None, fields))
        if ('Track' not in fields):
            layer = int(fields[0])
            time = float(fields[5])
            ID = float(fields[9])
            if (layer >= 10) and (layer <= 15):
                y_local = float(fields[11])
            else:
                y_local = -10000
            pos_layer[layer] = y_local
            ID_layer[layer] = ID
            t_layer[layer] = time;
            layerNum.append(layer)


#print IDonTrack
#//////////////////////////////////////////
#// Read ROOT Tree file
#//////////////////////////////////////////

print('Reading ROOT File')

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

eta_hists.append(ROOT.TH1F("Eta_rdt", "Eta_rdt; Eta; Number of Entries", 30, -0.5,1.5))
eta_hists.append(ROOT.TH1F("Eta_flt", "Eta_flt; Eta; Number of Entries", 30, -0.5,1.5))


eta_hists_v_pos = ROOT.TH2F("eta_pos", "eta_pos; #eta; sub-cell position (mm)", 30,-0.5,1.5, 51, 0, 51e-3)

skewed_hists = []
skewed_hists.append(ROOT.TH1F("skewed_total", "skewed_total; Eta; Number of Entries", 30, -0.5,1.5))
skewed_hists.append(ROOT.TH1F("skewed_1", "skewed_1; Eta; Number of Entries", 30, -0.5,1.5))
skewed_hists.append(ROOT.TH1F("skewed_2", "skewed_2; Eta; Number of Entries", 30, -0.5,1.5))
skewed_hists.append(ROOT.TH1F("skewed_3+", "skewed_3+; Eta; Number of Entries", 30, -0.5,1.5))

hist_profile = []
for j in range(0,4,1):
    histname = "sub_cell_entries_csize"+str(j)
    hist_profile.append(ROOT.TH1F(histname+"_timed", "sub_cell_entries; distance (#mum); Number of Entries ", 50, 0, 50))




f = ROOT.TFile(args.root_in)
stripTree = f.Get("vectorClusterStrip")

print('Performing analysis')
charge_right = 0
charge_left = 0
charge_right_self = 0
charge_left_self = 0
count = 1
for entry in stripTree:
    for i in range(len(entry.sensor)):
        if fitpos.get(entry.ID[i], -99) != -99:
            layer[0] = entry.sensor[i]
            charge[0] = entry.charge[i]
            signi[0] = entry.signi[i]
            pos[0] = entry.pos[i]
            #print "Position ", pos[0]
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
                #print "charge_right ", charge_right
                #print "charge_left ", charge_left
                eta_hists[0].Fill(eta)
                eta_hists_v_pos.Fill(eta, fitpos.get(entry.ID[i])%50)
                #print "Subcell pos ", fitpos.get(entry.ID[i])
                #print "Eta ", eta
                if (count == 1):
                    eta_hists[1].Fill(eta)
                    hist_profile[1].Fill((fitpos.get(entry.ID[i])%50e-3)*1000)
                    hist_profile[0].Fill((fitpos.get(entry.ID[i])%50e-3)*1000)
                elif (count == 2):
                    eta_hists[2].Fill(eta)
                    if eta > 0 and eta < 1:
                        hist_profile[2].Fill((fitpos.get(entry.ID[i])%50e-3)*1000)
                        hist_profile[0].Fill((fitpos.get(entry.ID[i])%50e-3)*1000)
                else:
                    eta_hists[3].Fill(eta)
                if ((fitpos.get(entry.ID[i])%50e-3)*1000 >= 12.5) and ((fitpos.get(entry.ID[i])%50e-3)*1000 < 37.5):
                    eta_hists[5].Fill(eta)
                else:
                    eta_hists[4].Fill(eta)
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

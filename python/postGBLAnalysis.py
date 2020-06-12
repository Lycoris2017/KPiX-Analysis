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

test = {}
with open(args.GBL_in) as GBLFile:
    line = GBLFile.readline()
    line = GBLFile.readline()
    while line:
        if "run" in line:
            line = GBLFile.readline()
        else:
            column = line.split( )
            test[int(column[8])] = column[9]
            line = GBLFile.readline()

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


print 'Reading ROOT File'

f = ROOT.TFile(args.root_in)
stripTree = f.Get("vector_sub_cluster")
for entry in stripTree:
    for i in xrange(len(entry.sensor)):
        if test.get(entry.ID[i], 'False') == 'True':
            layer[0] = entry.sensor[i]
            charge[0] = entry.charge[i]
            signi[0] = entry.signi[i]
            pos[0] = entry.pos[i]
            ID[0] = entry.ID[i]
            stripOnTrackTree.Fill()







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

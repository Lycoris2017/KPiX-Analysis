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

ROOT.TH1.SetDefaultSumw2()

parser = MyParser()
parser.add_argument('file_in', help='name of the input file')
args = parser.parse_args()
if len(sys.argv) < 2:
	print parser.print_help()
	sys.exit(1)
print ''

maxEvents = 100000
outfile = args.file_in[:-4]+'.root'
outHistFile = ROOT.TFile.Open(outfile, "RECREATE")
outHistFile.cd()
charge_hist = ROOT.TH1F("Charge", "Cluster Charge; Charge (fC); Number of Entries ", 200, 0, 20)
charge_size1 = ROOT.TH1F("Charge_on_Track_csize1", "Cluster Charge; Charge (fC); Number of Entries ", 200, 0, 20)
charge_size2 = ROOT.TH1F("Charge_on_Track_csize2", "Cluster Charge; Charge (fC); Number of Entries ", 200, 0, 20)
charge_size3 = ROOT.TH1F("Charge_on_Track_csize3+", "Cluster Charge; Charge (fC); Number of Entries ", 200, 0, 20)
charge_size1_3 = ROOT.TH1F("Charge_on_Track_csize1_3+", "Cluster Charge; Charge (fC); Number of Entries ", 200, 0, 20)
significance_hist = ROOT.TH1F("Significance", "Cluster Significance; S/N; Number of Entries ", 500, 0, 50)
y_hist = ROOT.TH1F("y", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46)
x_hist = ROOT.TH1F("x", "Cluster x; position (mm); Number of Entries ", 1840, -46, 46)
corr_hist = ROOT.TH1F("Correlations", "Correlations; num of correlated planes; Number of Hits ", 10, 0, 9)
correlation_l10_l11 = ROOT.TH2F("correlation", "correlation y; Layer 10 (mm); Layer 11 (mm); Number of Hits", 1840, -46, 46, 1840, -46, 46)



size_hist = ROOT.TH1F("cluster_size", "Cluster_Size; size; Number of Entries",  10, -0.5, 9.5)

trackhits_hist = ROOT.TH1F("Track_hits", "Track_hits; #Hits on track; Number of Entries ", 25, -0.5, 24.5)
HitsPerTrack_hist = ROOT.TH1F("hits_per_track", "hits_per_track; Number of hits per track; Number of Entries ", 7, -0.5, 6.5)
#trackhits_per_track_hist = ROOT.TH1F("Track_hits_per_track", "Track_hits_per_track; #Hits on track; Number of Entries ", 7, -0.5, 6.5)
trackNum = ROOT.TH1F("Nr._of_Tracks", "Nr._of_Tracks; #Nr. of Tracks; Number of Entries ", 10, -0.5, 9.5)

significance_hists =  []
significance_hists.append(ROOT.TH1F("Significance_l10", "Cluster Significance; S/N; Number of Entries ", 200, 0, 100))
significance_hists.append(ROOT.TH1F("Significance_l11", "Cluster Significance; S/N; Number of Entries ", 200, 0, 100))
significance_hists.append(ROOT.TH1F("Significance_l12", "Cluster Significance; S/N; Number of Entries ", 200, 0, 100))
significance_hists.append(ROOT.TH1F("Significance_l13", "Cluster Significance; S/N; Number of Entries ", 200, 0, 100))
significance_hists.append(ROOT.TH1F("Significance_l14", "Cluster Significance; S/N; Number of Entries ", 200, 0, 100))
significance_hists.append(ROOT.TH1F("Significance_l15", "Cluster Significance; S/N; Number of Entries ", 200, 0, 100))

charge_hists = []
charge_hists.append(ROOT.TH1F("Charge_l10", "Cluster Charge; Charge (fC); Entries ", 200, 0, 20))
charge_hists.append(ROOT.TH1F("Charge_l11", "Cluster Charge; Charge (fC); Entries ", 200, 0, 20))
charge_hists.append(ROOT.TH1F("Charge_l12", "Cluster Charge; Charge (fC); Entries ", 200, 0, 20))
charge_hists.append(ROOT.TH1F("Charge_l13", "Cluster Charge; Charge (fC); Entries ", 200, 0, 20))
charge_hists.append(ROOT.TH1F("Charge_l14", "Cluster Charge; Charge (fC); Entries ", 200, 0, 20))
charge_hists.append(ROOT.TH1F("Charge_l15", "Cluster Charge; Charge (fC); Entries ", 200, 0, 20))


y_hists = []
y_hists.append(ROOT.TH1F("y_l10", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46))
y_hists.append(ROOT.TH1F("y_l11", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46))
y_hists.append(ROOT.TH1F("y_l12", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46))
y_hists.append(ROOT.TH1F("y_l13", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46))
y_hists.append(ROOT.TH1F("y_l14", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46))
y_hists.append(ROOT.TH1F("y_l15", "Cluster y; position (mm); Number of Entries ", 1840, -46, 46))

x_hists = []
x_hists.append(ROOT.TH1F("x_l10", "Cluster x; position (mm); Number of Entries ", 1840, -1, 1))
x_hists.append(ROOT.TH1F("x_l11", "Cluster x; position (mm); Number of Entries ", 1840, -1, 1))
x_hists.append(ROOT.TH1F("x_l12", "Cluster x; position (mm); Number of Entries ", 1840, -1, 1))
x_hists.append(ROOT.TH1F("x_l13", "Cluster x; position (mm); Number of Entries ", 1840, -1, 1))
x_hists.append(ROOT.TH1F("x_l14", "Cluster x; position (mm); Number of Entries ", 1840, -1, 1))
x_hists.append(ROOT.TH1F("x_l15", "Cluster x; position (mm); Number of Entries ", 1840, -1, 1))

size_hists = []
size_hists.append(ROOT.TH1F("size_l10", "size; #strips; Number of Entries ", 10, -0.5, 9.5))
size_hists.append(ROOT.TH1F("size_l11", "size; #strips; Number of Entries ", 10, -0.5, 9.5))
size_hists.append(ROOT.TH1F("size_l12", "size; #strips; Number of Entries ", 10, -0.5, 9.5))
size_hists.append(ROOT.TH1F("size_l13", "size; #strips; Number of Entries ", 10, -0.5, 9.5))
size_hists.append(ROOT.TH1F("size_l14", "size; #strips; Number of Entries ", 10, -0.5, 9.5))
size_hists.append(ROOT.TH1F("size_l15", "size; #strips; Number of Entries ", 10, -0.5, 9.5))


all_hists = [charge_hist, significance_hist, y_hist, x_hist, corr_hist, charge_hists, y_hists[0:3], x_hists[0:3], size_hists[0:3], significance_hists, size_hist]

plane = []
x = []
y = []
z = []
corr = []
sig = []
size = []
charge = []
HitOnTrack = []
y_l10 = []
y_l11 = []
counter = 0
layer=array('i', [0])
ccharge=array('f', [0])
significance=array('f', [0])
noise=array('f', [0])
csize=array('i', [0])
yPos=array('f', [0])
onTrack=array('i', [0])

corrClusterTree = ROOT.TTree("corr_clusters", "corr_clusters");
corrClusterTree.Branch("layer", layer, 'layer/I');
corrClusterTree.Branch("charge", ccharge, 'charge/F');
corrClusterTree.Branch("significance", significance, 'significance/F');
corrClusterTree.Branch("noise", noise, 'noise/F');
corrClusterTree.Branch("csize", csize, 'csize/I');
corrClusterTree.Branch("yPos", yPos, 'yPos/F');
corrClusterTree.Branch("onTrack", onTrack, 'onTrack/I');
prev_run = '0'
track0Counter = 0
with open(args.file_in) as inFile:
    line = inFile.readline()
    hits_on_track = 0

    while line and counter < maxEvents:
        a = line.split( )
        if "run" in line:
            counter+=1
            if (len(y_l10) is not 0  and len(y_l11) is not 0):
                for q in y_l10:
                    for h in y_l11:
                        correlation_l10_l11.Fill(float(q),float(h))
                y_l10 = []
                y_l11 = []
            trackhits_hist.Fill(hits_on_track)
            if (counter > 1):
                trackNum.Fill(0)#int(tracks))
            hits_on_track=0
            if ('0' in tracks):
                track0Counter += 1
                print prev_run
            prev_run = a[2]
            line = inFile.readline()
        else:
            plane.append(a[0])
            tracks = a[10]
            x.append(a[1])
            y.append(a[2])
            z.append(a[3])
            corr.append(a[4])
            sig.append(a[5])
            size.append(a[6])
            charge.append(a[7])
            if "True" in a[9]:
                hits_on_track += 1
                HitOnTrack.append(1)
            else:
               HitOnTrack.append(0)
            if "10" in a[0]:
                y_l10.append(a[2])
            elif "11" in a[0]:
                y_l11.append(a[2])
            line = inFile.readline()

#c1 = ROOT.TCanvas( 'test', 'Test', 1600, 900 )
#c1.cd()

for q in xrange(1,len(charge)):
    corr_hist.Fill(int(corr[q]))
    if (int(corr[q]) >= 1):
        layer[0]=int(plane[q])
        ccharge[0]=float(charge[q])
        noise[0]=float(charge[q])/float(sig[q])
        significance[0]=float(sig[q])
        yPos[0]=float(y[q])
        csize[0]=int(size[q])
        onTrack[0]=HitOnTrack[q]
        corrClusterTree.Fill()
    if (HitOnTrack[q] == 1):
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
        elif ("13" in plane[q]):
            charge_hists[3].Fill(float(charge[q]))
            significance_hists[3].Fill(float(sig[q]))
            y_hists[3].Fill(float(y[q]))
            x_hists[3].Fill(float(x[q]))
            size_hists[3].Fill(int(size[q]))
        elif ("14" in plane[q]):
            charge_hists[4].Fill(float(charge[q]))
            significance_hists[4].Fill(float(sig[q]))
            y_hists[4].Fill(float(y[q]))
            x_hists[4].Fill(float(x[q]))
            size_hists[4].Fill(int(size[q]))
        elif ("15" in plane[q]):
            charge_hists[5].Fill(float(charge[q]))
            significance_hists[5].Fill(float(sig[q]))
            y_hists[5].Fill(float(y[q]))
            x_hists[5].Fill(float(x[q]))
            size_hists[5].Fill(int(size[q]))
        size_hist.Fill(int(size[q]))
        charge_hist.Fill(float(charge[q]))
        if (int(size[q]) == 1):
           charge_size1.Fill(float(charge[q]))
           charge_size1_3.Fill(float(charge[q]))
        elif (int(size[q]) == 2):
           charge_size2.Fill(float(charge[q]))
        elif (int(size[q]) >= 3):
           charge_size3.Fill(float(charge[q]))
           charge_size1_3.Fill(float(charge[q]))
        significance_hist.Fill(float(sig[q]))
        y_hist.Fill(float(y[q]))
        x_hist.Fill(float(x[q]))
for i in all_hists:
    if isinstance(i, list):
        for j in i:
            j.Scale(1.0/j.GetEntries())
    else:
        i.Scale(1.0/i.GetEntries())

outHistFile.Write()
outHistFile.Close()
#print track0Counter
#charge_hist.Draw("hist e")
#c1.Modified()
#c1.Update()
#c1.SaveAs("test.png")
#c1.Close()

#!usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import string
import ROOT
import argparse
import argcomplete
from operator import add
import sys


left_strip = None
right_strip = None

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def loopdir(keys):  # loop through all subdirectories of the root file and look for the vs channel histograms
	global left_strip, right_strip ##declare that the beforementioned variables are global variables
	for key_object in keys:
		if ('TDirectory' in key_object.GetClassName()):
			loopdir(key_object.ReadObj().GetListOfKeys())
		else:
			#print key_object.GetName()
			if 'Left_Strip_entries_k_26' in key_object.GetName() and '_total' in key_object.GetName():
				left_strip = key_object
			if 'Right_Strip_entries_k_28' in key_object.GetName() and '_total' in key_object.GetName():
				right_strip = key_object


parser = MyParser()
parser.add_argument('file_in', help='name of the input file')
args = parser.parse_args()

if len(sys.argv) < 2:
	print parser.print_help()
	sys.exit(1)


rootfile = ROOT.TFile(args.file_in)
key_root = rootfile.GetListOfKeys()
file_base_name = args.file_in[args.file_in.find('/20')+1:args.file_in.rfind('.bin')]
loopdir(key_root)
if left_strip:
    left_strip_hist = left_strip.ReadObj()
if right_strip:
    right_strip_hist = right_strip.ReadObj()


Input_left = []
Input_right = []
C = []

for chan in xrange(920):
	if left_strip:
            Input_left.append(left_strip_hist.GetBinContent(chan+1))
	if right_strip:
            Input_right.append(right_strip_hist.GetBinContent(chan+1))
	
if Input_left:
    C = [Input_left]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    axes = plt.axes() 
    axes.get_yaxis().set_visible(False)
    axes.set_xlim(0,920)
    pcm = ax.pcolormesh(C, cmap='viridis')
    plt.colorbar(pcm)
    filename = '/home/lycoris-dev/Documents/'+file_base_name+'_strip_map_left.png'
    print "File is saved in " +  filename
    plt.savefig(filename, dpi = 300)

if Input_right:
    C = [Input_right]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    axes = plt.axes() 
    axes.get_yaxis().set_visible(False)
    axes.set_xlim(0,920)
    pcm = ax.pcolormesh(C, cmap='viridis')
    plt.colorbar(pcm)
    filename = '/home/lycoris-dev/Documents/'+file_base_name+'_strip_map_right.png'
    print "File is saved in " +  filename
    plt.savefig(filename, dpi = 300)

plt.close()

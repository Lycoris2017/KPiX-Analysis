#!usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import string
import ROOT
import argparse
import argcomplete
from operator import add
import sys


left_strips = []
right_strips = []

left_strip = None
left_strip_b0 = None
left_strip_b1 = None
left_strip_b2 = None
left_strip_b3 = None
right_strip = None
right_strip_b0 = None
right_strip_b1 = None
right_strip_b2 = None
right_strip_b3 = None

timed_left_strip = None
timed_right_strip = None

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def loopdir(keys):  # loop through all subdirectories of the root file and look for the vs channel histograms
	global left_strip, left_strip_b0, left_strip_b1, left_strip_b2, left_strip_b3
	global right_strip, right_strip_b0, right_strip_b1, right_strip_b2, right_strip_b3 ##declare that the beforementioned variables are global variables
	for key_object in keys:
		if ('TDirectory' in key_object.GetClassName()):
			loopdir(key_object.ReadObj().GetListOfKeys())
		else:
			if 'Left_Strip' in key_object.GetName() and 'k'+args.kpix in key_object.GetName():
				print key_object.GetName()
				left_strips.append(key_object)
			if 'Right_Strip' in key_object.GetName()  and 'k'+args.kpix in key_object.GetName():
				print key_object.GetName()
				right_strips.append(key_object)
			#if 'Left_Strip_entries_k_'+args.kpix in key_object.GetName() and '_total' in key_object.GetName():
				#left_strips.append(key_object)
			#if 'Left_Strip_entries_k_'+args.kpix in key_object.GetName() and '_b0' in key_object.GetName():
				#left_strips.append(key_object)
			#if 'Left_Strip_entries_k_'+args.kpix in key_object.GetName() and '_b1' in key_object.GetName():
				#left_strips.append(key_object)
			#if 'Left_Strip_entries_k_'+args.kpix in key_object.GetName() and '_b2' in key_object.GetName():
				#left_strips.append(key_object)
			#if 'Left_Strip_entries_k_'+args.kpix in key_object.GetName() and '_b3' in key_object.GetName():
				#left_strips.append(key_object)
			#if 'Right_Strip_entries_k_'+args.kpix in key_object.GetName() and '_total' in key_object.GetName():
				#right_strip = key_object
			#if 'Right_Strip_entries_k_'+args.kpix in key_object.GetName() and '_b0' in key_object.GetName():
				#right_strip_b0 = key_object
			#if 'Right_Strip_entries_k_'+args.kpix in key_object.GetName() and '_b1' in key_object.GetName():
				#right_strip_b1 = key_object
			#if 'Right_Strip_entries_k_'+args.kpix in key_object.GetName() and '_b2' in key_object.GetName():
				#right_strip_b2 = key_object
			#if 'Right_Strip_entries_k_'+args.kpix in key_object.GetName() and '_b3' in key_object.GetName():
				#right_strip_b3 = key_object
			#if 'timed_left_strip_entries_k_'+args.kpix  in key_object.GetName() and '_total' in key_object.GetName():
				#left_strips.append(key_object)
			#if 'timed_right_strip_entries_k_'+args.kpix  in key_object.GetName() and '_total' in key_object.GetName():
				#timed_right_strip = key_object





parser = MyParser()
parser.add_argument('file_in', help='name of the input file')
parser.add_argument('-k', dest='kpix', help='which kpix to check, eg: 1 or 1 2')
args = parser.parse_args()

if len(sys.argv) < 2:
	print parser.print_help()
	sys.exit(1)


rootfile = ROOT.TFile(args.file_in)
key_root = rootfile.GetListOfKeys()
file_base_name = args.file_in[args.file_in.find('/20')+1:args.file_in.rfind('.bin')]
loopdir(key_root)


left_strip_hists = []
right_strip_hists = []
for i in left_strips:
	left_strip_hists.append(i.ReadObj())

for i in right_strips:
	right_strip_hists.append(i.ReadObj())



C = []
timed_Input_left = []


 

for i in left_strip_hists:
	inputs = []
	C = []
	for chan in xrange(920):
		inputs.append(i.GetBinContent(chan+1))
	maximum = max(inputs)
	C = [inputs]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	axes = plt.axes() 
	axes.get_yaxis().set_visible(False)
	axes.set_xlim(0,920)
	pcm = ax.pcolormesh(C, cmap='viridis', vmin=0, vmax=maximum)
	plt.colorbar(pcm)
	filename = '/home/lycoris-dev/Documents/'+file_base_name+'_'+i.GetName()+'.png'
	print "File is saved in " +  filename
	plt.savefig(filename, dpi = 300)
	
for i in right_strip_hists:
	inputs = []
	C = []
	for chan in xrange(920):
		inputs.append(i.GetBinContent(chan+1))
	maximum = max(inputs)
	C = [inputs]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	axes = plt.axes() 
	axes.get_yaxis().set_visible(False)
	axes.set_xlim(0,920)
	pcm = ax.pcolormesh(C, cmap='viridis', vmin=0, vmax=maximum)
	plt.colorbar(pcm)
	filename = '/home/lycoris-dev/Documents/'+file_base_name+'_'+i.GetName()+'.png'
	print "File is saved in " +  filename
	plt.savefig(filename, dpi = 300)


plt.close()

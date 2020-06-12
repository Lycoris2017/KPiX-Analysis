#!usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import string
import ROOT
import argparse
import argcomplete
from operator import add
from pyroot_functions import loopdir_new
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



parser = argparse.ArgumentParser()
parser.add_argument(
	'file_in',
	nargs='+',
	help='name of the input file'
)
parser.add_argument(
	'-n', '--name',
	dest='name',
	default=['everything'],
	nargs='*',
	help='used to specify the name of the plot which should be used'
)
parser.add_argument(
	'-o', '--output',
	dest='output_name',
	help='specifies the name and type of the output file (e.g. test.png, comparison.root etc...'
)
args = parser.parse_args()


root_file_list = []
filename_list = []

for root_file in args.file_in:
	root_file_list.append(ROOT.TFile(root_file))
	if '/Run_' in root_file:
		filename_list.append(root_file[root_file.find('/Run_')+1:root_file.rfind('.dat')+1])
	elif '/Calibration_' in root_file:
		filename_list.append(root_file[root_file.find('/Calibration_')+1:root_file.rfind('.dat')+1])

print filename_list
object_list = []
for x in root_file_list:
	key_root = x.GetListOfKeys()
	object_list= object_list + (loopdir_new(key_root, args.name))



C = []
timed_Input_left = []


 

for i in object_list:
	histogram  = i.ReadObj()
	inputs = []
	C = []
	for chan in xrange(1840):
		inputs.append(histogram.GetBinContent(chan+1))
	maximum = max(inputs)
	print inputs
	C = [inputs]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	axes = plt.axes() 
	axes.get_yaxis().set_visible(False)
	axes.set_xlim(0,1840)
	pcm = ax.pcolormesh(C, cmap='viridis', vmin=0, vmax=maximum)
	plt.colorbar(pcm)
	filename = '/home/lycoris-dev/Documents/'+args.output_name
	print "File is saved in " +  filename
	plt.savefig(filename, dpi = 300)
	
#for i in right_strip_hists:
#	inputs = []
#	C = []
#	for chan in xrange(920):
#		inputs.append(i.GetBinContent(chan+1))
#	maximum = max(inputs)
#	C = [inputs]
#	fig = plt.figure()
#	ax = fig.add_subplot(111)
#	axes = plt.axes()
#	axes.get_yaxis().set_visible(False)
#	axes.set_xlim(0,920)
#	pcm = ax.pcolormesh(C, cmap='viridis', vmin=0, vmax=maximum)
#	plt.colorbar(pcm)
#	filename = '/home/lycoris-dev/Documents/'+file_base_name+'_'+i.GetName()+'.png'
#	print "File is saved in " +  filename
#	plt.savefig(filename, dpi = 300)


plt.close()

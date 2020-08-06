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
import pyroot_functions as pf

ROOT.gROOT.SetBatch(True)

hist_list= []  #Global list of all chosen histograms (only works because its a list, with single variables they need to be declared as global in the function
graph_list = []





def arrangeStats(hists, statBoxW, statBoxH, name):
	i=0
	for h in hists:
		statBox = h.GetListOfFunctions().FindObject("stats")
		statBox.SetName('statBox' + str(i))
		
		listOfLines = statBox.GetListOfLines()

		tconst = statBox.GetLineWith(h.GetName());
		listOfLines.Remove(tconst);
		
		myt = ROOT.TLatex(0,0,name[len(hists)-i-1]);
		myt.SetTextFont(42);
		myt.SetTextColor(h.GetLineColor());
		listOfLines.AddFirst(myt)
		
		print name[i]
		#statBox.SetTitle(str(name[i]))
		statBox.SetX1NDC(args.legendloc[0] - statBoxW)
		statBox.SetY1NDC(args.legendloc[1] - i*(statBoxH+0.01) - statBoxH)
		statBox.SetX2NDC(args.legendloc[0])
		statBox.SetY2NDC(args.legendloc[1] - i*(statBoxH+0.01))
		statBox.SetTextColor(h.GetLineColor()) 
		#statBox.SetBorderColor(h.GetLineColor())
		#statBox.SetBorderSize(2)
		statBox.Draw("same")
		i+=1


#parser = argparse.ArgumentParser() #Command line argument parser.
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
	'-d', '--draw',
	dest='draw_option',
	default='',
	help='specify the drawing option as given by the root draw option, needs to be given as a single string (e.g. hist same or hist same multifile'
)

parser.add_argument(
	'-c', '--conditions',
	dest='conditions',
	default='',
	help='specify the root ttree conditions such as only ploy with y > 0 etc.'
)

parser.add_argument(
	'-v', '--variables',
	dest='variables',
	nargs='*',
	help='specify which ttree variable should be plotted'
)

parser.add_argument(
	'-o', '--output',
	dest='output_name',
	help='specifies the name and type of the output file (e.g. test.png, comparison.root etc...'
)
parser.add_argument(
	'--xrange',
	dest='xaxisrange',
	nargs='*',
	type=float,
	help='set a xrange for the plot to used with xmin xmax as the two arguments | type=float'
)
parser.add_argument(
	'--yrange',
	dest='yaxisrange',
	nargs='*',
	type=float,
	help='set a yrange for the plot to used with ymin ymax as the two arguments | type=float'
)
parser.add_argument(
	'--zrange',
	dest='zaxisrange',
	default=[9999],
	nargs='*',
	type=float,
	help='set a zrange for the plot to used with ymin ymax as the two arguments | type=float'
)
parser.add_argument(
	'--legend',
	dest='legend',
	default=None,
	nargs='*',
	help='list of names to be used as legend titles instead of the default filename+histogram name'
)
parser.add_argument(
	'--ylog',
	dest='ylog',
	default=False,
	action='store_true',
	help='if given as an option, set y axis to logarithmic. Remember to set the yrange to start above 0!'
)
parser.add_argument(
	'--zlog',
	dest='zlog',
	help='if given as an option, set z axis to logarithmic.'
)
parser.add_argument(
	'--color',
	dest='color',
	default=0,
	help='list of colors to be used'
)
#parser.add_argument('--color', dest='color', default=[590, 591, 593, 596, 600, 602, 604, 880, 860, 632, 840, 614], nargs='*', help='list of colors to be used')
parser.add_argument(
	'--xtitle',
	dest='xtitle',
	help='choose the name of the x axis title'
)
parser.add_argument(
	'--ytitle',
	dest='ytitle',
	help='choose the name of the y axis title'
)
parser.add_argument(
	'-b', '--binRange',
	dest='bin_range',
	nargs='*',
	help='set the number of bins'
)
parser.add_argument(
	'--nobox',
	dest='nobox',
	default=False,
	action='store_true',
	help='suppresses tstatbox from being added. Does not work with stack/nostack option'
)
parser.add_argument(
	'--order',
	dest='order',
	nargs='+',
	type=int,
	help='choose the order of plotting with same (to ensure no histograms overlap)'
)
args = parser.parse_args()

if (len(args.variables) == 0 or len(args.variables) > 2):
	print "Please enter 1 or 2 variables not more and not less"
	sys.exit()

mystyle = ROOT.TStyle("mystyle", "My Style")
mystyle, myMarker, myMarkerSize = pf.myROOTStyle(args.nobox)
myMarkerStyle = [myMarker, myMarkerSize]

mystyle.cd()
ROOT.gROOT.ForceStyle()
ROOT.gStyle.ls()
# finish setting root style

#legend_location = [0.65,0.65,0.98,0.85] # x_left, y_bottom, x_right, y_top
legend_location = [0.15,0.65,0.35,0.85] # x_left, y_bottom, x_right, y_top



root_file_list = []
filename_list = []
##-----------------
##loop through all given files and add them to the list. then loop through the keys for every file..
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
	#object_list = object_list + (loopdir_new(key_root, args.name))
	object_list.append((pf.loopdir_new(key_root, args.name)[0]))
	
folder_loc = '/scratch/plots/testbeam202003/'
##-----------------	
##general output
#print args.color

print 'Looking for objects'
print '----------------------'
print 'Name contains ', args.name
print 'Number of objects found is: ', len(object_list)
print object_list	
if (args.ylog):
	print 'Setting y axis to log, only works if the range was specified to start at y_min > 0'
##------------------
##start of the plotting.
if args.conditions:
	conditionsName = args.conditions.replace("<", "L")
	conditionsName = args.conditions.replace(">", "G")
if (args.output_name):
	if args.conditions:
		outName  = args.output_name+"_Cond_"+conditionsName
	else:
		outName  = args.output_name
else:
	stringStart = args.file_in[0].find('/testbeam')
	outName = 'Histogram'#'/scratch/plots/testbeam202003/'+args.file_in[0][stringStart:-5]
	for i in args.variables:
		if args.conditions:
			outName = outName+"_"+i+"_Cond_"+conditionsName
		else:
			outName = outName+"_"+i
print "Output name ", outName
bin_range = args.bin_range
variables = args.variables
draw_option = args.draw_option
yMax = None
yMin = None
if (len(object_list) is not 0):
	myHistograms = pf.plot_tree(object_list, args.conditions, variables ,bin_range, outName, )
else:
	print 'There are NO valid histograms/graphs in the current selection'
	print ''
	print ''
print myHistograms
yRangeScale = 1.5

if (args.xtitle):
	for i in myHistograms:
		i.GetXaxis().SetTitle(args.xtitle)
if (args.ytitle):
	for i in myHistograms:
		i.GetYaxis().SetTitle(args.ytitle)
if (args.yaxisrange):
	for i in myHistograms:
		i.GetYaxis().SetRangeUser(args.yaxisrange[0], args.yaxisrange[1])
else:
	for i in myHistograms:
		if (yMax is None):
			yMax = i.GetMaximum()
		elif(yMax < yRangeScale*i.GetMaximum()):
			yMax =  yRangeScale*i.GetMaximum()
		if (args.ylog):
			if (yMin is None):
				yMin = 1.0/(yRangeScale*i.GetEntries())
			elif (yMin > 1.0/(yRangeScale*i.GetEntries())):
				yMin = 1.0/(yRangeScale*i.GetEntries())
		else:
			if (yMin is None):
				yMin = (1./yRangeScale)*i.GetMinimum()
			elif (yMin >  (1./yRangeScale)*i.GetMinimum()):
				yMin =  (1.0/yRangeScale)*i.GetMinimum()
	print "Minimum y value ",yMin
	print "Maximum y value ",yMax
	if len(args.variables) == 1:
		for i in myHistograms:
			i.GetYaxis().SetRangeUser(yMin, yMax)
if (args.xaxisrange):
	for i in myHistograms:
		i.GetXaxis().SetRangeUser(args.xaxisrange[0], args.xaxisrange[1])

if (args.xaxisrange):
	xaxisrange = args.xaxisrange
pf.drawSame(myHistograms, draw_option, args.legend, myMarkerStyle, outName, args.ylog, args.order)
for x in root_file_list:
	ROOT.gROOT.GetListOfFiles().Remove(x)

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
from pyroot_functions import loopdir_new

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
	'-o', '--output',
	dest='output_name',
	help='specifies the name and type of the output file (e.g. test.png, comparison.root etc...'
)
parser.add_argument(
	'--xrange',
	dest='xaxisrange',
	default=[9999],
	nargs='*',
	type=float,
	help='set a xrange for the plot to used with xmin xmax as the two arguments | type=float'
)
parser.add_argument(
	'--yrange',
	dest='yaxisrange',
	default=[9999],
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
	nargs='*',
	help='list of names to be used as legend titles instead of the default filename+histogram name'
)
parser.add_argument(
	'--ylog',
	dest='ylog',
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
	default=[ 861, 1, 418,  810, 402,  908, 435, 880,60, 632, 840, 614],
	nargs='*',
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
args = parser.parse_args()


mystyle = ROOT.TStyle("mystyle", "My Style")


mystyle.SetPaintTextFormat("5.3f");


#set the background color to white
mystyle.SetFillColor(0)
mystyle.SetFrameFillColor(0)
mystyle.SetCanvasColor(0)
mystyle.SetPadColor(0)
mystyle.SetTitleFillColor(0)
mystyle.SetStatColor(0)

#dont put a colored frame around the plots
mystyle.SetFrameBorderMode(0)
mystyle.SetCanvasBorderMode(0)
mystyle.SetPadBorderMode(0)
mystyle.SetLegendBorderSize(0)
#
##use the primary color palette
##mystyle.SetPalette(1,0)
#
##set the default line color for a histogram to be black
#mystyle.SetHistLineColor(1)
#

#
##make the axis labels black
#mystyle.SetLabelColor(1,"xyz")
#
##set the default title color to be black
#mystyle.SetTitleColor(1)
mystyle.SetOptTitle(0)
#
##set the margins
mystyle.SetPadBottomMargin(0.16)
mystyle.SetPadTopMargin(0.05)
mystyle.SetPadRightMargin(0.15)
mystyle.SetPadLeftMargin(0.16)
#
##set axis label and title text sizes
mystyle.SetLabelFont(62,"xyz")
mystyle.SetLabelSize(0.04,"xyz")
mystyle.SetLabelOffset(0.003,"yz")
mystyle.SetLabelOffset(0.00,"x")
mystyle.SetTitleFont(62,"xyz")
mystyle.SetTitleSize(0.056,"xyz")
mystyle.SetTitleOffset(1.1,"y")
mystyle.SetTitleOffset(0.8,"z")
mystyle.SetTitleOffset(0.75,"x")
mystyle.SetStatFont(62)
mystyle.SetStatFontSize(0.05)


ROOT.TGaxis.SetMaxDigits(4)

#mystyle.SetTitleBorderSize(0)
#mystyle.SetStatBorderSize(0)
#mystyle.SetTextFont(42)

##set legend text size etc.
mystyle.SetLegendTextSize(0.04)
#
##set line widths
mystyle.SetFrameLineWidth(2)
mystyle.SetFuncWidth(2)
##set the default line color for a fit function to be red
mystyle.SetFuncColor(2)
mystyle.SetHistLineWidth(2)
#
##set the number of divisions to show
#mystyle.SetNdivisions(506, "xy")
#
##turn off xy grids
#mystyle.SetPadGridX(0)
#mystyle.SetPadGridY(0)
#
##set the tick mark style
#mystyle.SetPadTickX(1)
#mystyle.SetPadTickY(1)
#
##turn off stats
#mystyle.SetOptStat(0) ##removes stat box

#mystyle.SetOptStat(args.statBox)
#mystyle.SetOptFit(args.fitBox)
mystyle.SetOptStat(1111)
mystyle.SetOptFit(111)
#mystyle.SetOptStat(0000001) #only name
#
##marker settings
mystyle.SetMarkerStyle(8)
mystyle.SetMarkerSize(0.7)
mystyle.SetLineWidth(2)

#done
#mystyle.cd()
#ROOT.gROOT.ForceStyle()
#ROOT.gStyle.ls()




#legend_location = [0.65,0.65,0.98,0.85] # x_left, y_bottom, x_right, y_top
legend_location = [0.15,0.65,0.35,0.85] # x_left, y_bottom, x_right, y_top

##-----------------	
##produce empty root file and filename lists.

# finish setting root style
mystyle.cd()
ROOT.gROOT.ForceStyle()
ROOT.gStyle.ls()


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
	object_list= object_list + (loopdir_new(key_root, args.name))
	
folder_loc = '/home/lycoris-dev/Documents/testbeam201907/'
##-----------------	
##general output
#print args.color

print 'Looking for objects'
print '----------------------'
print 'Name contains ', args.name
print 'Number of objects found is: ', len(object_list)
print object_list	
if (args.ylog and args.yaxisrange[0] is 0):
	print 'Setting y axis to log, only works if the range was specified to start at y_min > 0'
##------------------
##start of the plotting.

if ((len(hist_list) > len(args.color) or len(graph_list) > len(args.color))) and ("same" in args.draw_option):
	print 'You do not have enough colors ', len(args.color), 'for the number of histograms you have ', len(hist_list)
	sys.exit(1)
if (len(object_list) is not 0):
	print 'Need to implement actual plotting here'
else:
	print 'There are NO valid histograms/graphs in the current selection'
	print ''
	print ''
for x in root_file_list:
	ROOT.gROOT.GetListOfFiles().Remove(x)

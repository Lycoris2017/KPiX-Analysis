#!usr/bin/python

import numpy as np
import string
from matplotlib import cm
import ROOT 
from ROOT import TColor
import argparse
import argcomplete
from operator import add
import sys
from decimal import Decimal
from array import array

ROOT.gROOT.SetBatch(True)

hist_list= []  #Global list of all chosen histograms (only works because its a list, with single variables they need to be declared as global in the function
graph_list = []



class MyParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)


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



def loopdir(keys):  # loop through all subdirectories of the root file and add all plots+histograms that contain the same string as given by the name, channel and bucket arguments to the global list of all chosen histograms hist_list
    for key_object in keys:
		if ('TDirectory' in key_object.GetClassName()):
			loopdir(key_object.ReadObj().GetListOfKeys())
		else:
			if ('TGraph' in key_object.GetClassName()):
				if(args.name == 'everything'):
					#print key_object.GetName()
					graph_list.append(key_object)
				else:
					if args.exact:
						for name in args.name:
							if (name == key_object.GetName()):
								graph_list.append(key_object)
					else:
						for channel in args.channel:
							for strip in args.strip:
								for bucket in args.bucket:
									for kpix in args.kpix:
										if (bucket == 4):
											if ((all(name in key_object.GetName() for name in args.name) and all(refuse not in key_object.GetName() for refuse in args.refuse)) \
											and ('_c'+str(channel)+'_' in key_object.GetName() or channel == '9999') and ('_s'+str(strip)+'_' in key_object.GetName() or strip == '9999') and ('b' not in key_object.GetName()) and ('k'+str(kpix) in key_object.GetName() or 'k_'+str(kpix) in key_object.GetName() or kpix == 9999)):
												#print 'Histogram found: ', key_object.GetName()
												graph_list.append(key_object)
										else:
											if ((all(name in key_object.GetName() for name in args.name) and all(refuse not in key_object.GetName() for refuse in args.refuse))  \
											and ('_c'+str(channel)+'_' in key_object.GetName() or channel == '9999') and ('_s'+str(strip)+'_' in key_object.GetName() or strip == '9999') and ('b'+str(bucket) in key_object.GetName() or bucket == 9999) and ('k'+str(kpix) in key_object.GetName() or 'k_'+str(kpix) in key_object.GetName() or kpix == 9999)):
												#print 'Histogram found: ', key_object.GetName()
												graph_list.append(key_object)
			else:
				if(args.name == 'everything' and key_object.ReadObj().GetEntries() != 0):
					#print key_object.GetName()
					hist_list.append(key_object)
				else:
					if args.exact:
						for name in args.name:
							if (name == key_object.GetName()):
								hist_list.append(key_object)
					else:
						for channel in args.channel:
							for strip in args.strip:
								for bucket in args.bucket:
									for kpix in args.kpix:
										if (bucket == 4):
											if ((all(name in key_object.GetName() for name in args.name) and all(refuse not in key_object.GetName() for refuse in args.refuse)) \
											and ('_c'+str(channel)+'_' in key_object.GetName() or channel == '9999') and ('_s'+str(strip)+'_' in key_object.GetName() or strip == '9999') and ('b' not in key_object.GetName()) and ('k'+str(kpix) in key_object.GetName() or 'k_'+str(kpix) in key_object.GetName() or kpix == 9999)):
												#print 'Histogram found: ', key_object.GetName()
												hist_list.append(key_object)
										else:
											if ((all(name in key_object.GetName() for name in args.name) and all(refuse not in key_object.GetName() for refuse in args.refuse)) \
											and ('_c'+str(channel)+'_' in key_object.GetName() or channel == '9999') and ('_s'+str(strip)+'_' in key_object.GetName() or strip == '9999') and ('b'+str(bucket) in key_object.GetName() or bucket == 9999) and ('k'+str(kpix) in key_object.GetName() or 'k_'+str(kpix) in key_object.GetName() or kpix == 9999)):
												#print 'Histogram found: ', key_object.GetName()
												hist_list.append(key_object)

mystyle = ROOT.TStyle("mystyle", "My Style")

#set the background color to white
mystyle.SetFillColor(10)
mystyle.SetFrameFillColor(10)
mystyle.SetCanvasColor(10)
mystyle.SetPadColor(10)
mystyle.SetTitleFillColor(0)
mystyle.SetStatColor(10)

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
##set the default line color for a fit function to be red
#mystyle.SetFuncColor(2)
#
##make the axis labels black
#mystyle.SetLabelColor(1,"xyz")
#
##set the default title color to be black
#mystyle.SetTitleColor(1)
mystyle.SetOptTitle(0)
#
##set the margins
mystyle.SetPadBottomMargin(0.15)
mystyle.SetPadTopMargin(0.04)
mystyle.SetPadRightMargin(0.18)
mystyle.SetPadLeftMargin(0.12)
#
##set axis label and title text sizes
mystyle.SetLabelFont(42,"xyz")
mystyle.SetLabelSize(0.08,"xy")
mystyle.SetLabelSize(0.06,"z")
mystyle.SetLabelOffset(0.003,"yz")
mystyle.SetLabelOffset(0.00,"x")
mystyle.SetTitleFont(42,"xyz")
mystyle.SetTitleSize(0.08,"xyz")
mystyle.SetTitleOffset(0.75,"yz")
mystyle.SetTitleOffset(0.75,"x")
mystyle.SetStatFont(42)
mystyle.SetStatFontSize(0.03)

#mystyle.SetTitleBorderSize(0)
#mystyle.SetStatBorderSize(0)
#mystyle.SetTextFont(42)

##set legend text size etc.
mystyle.SetLegendTextSize(0.03)
#
##set line widths
mystyle.SetFrameLineWidth(2)
mystyle.SetFuncWidth(2)
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
mystyle.SetOptStat(1001111)
#mystyle.SetOptFit(111)
#
##marker settings
mystyle.SetMarkerStyle(20)
mystyle.SetMarkerSize(0.7)
mystyle.SetLineWidth(1)

# Colour Style

NRGBs = 5
ncontours = 30

stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
red   = [ 0.00, 0.00, 0.87, 1.00, 0.51 ]
green = [ 0.00, 0.81, 1.00, 0.20, 0.00 ]
blue  = [ 0.51, 1.00, 0.12, 0.00, 0.00 ]
stopsArray = array('d', stops)
redArray = array('d', red)
greenArray = array('d', green)
blueArray = array('d', blue)

TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, ncontours)
mystyle.SetNumberContours(ncontours)


#done
mystyle.cd()
ROOT.gROOT.ForceStyle()
ROOT.gStyle.ls()







#parser = argparse.ArgumentParser() #Command line argument parser.
parser = MyParser()
parser.add_argument('file_in', nargs='+', help='name of the input file')
parser.add_argument('-n', '--name', dest='name', default=['everything'], nargs='*',  help='used to specify the name of the plot which should be used')
parser.add_argument('-c', '--channel', dest='channel', default=['9999'], nargs='*', help='used to specify the channel of the plot which should be used')
parser.add_argument('-s', '--strip', dest='strip', default=['9999'], nargs='*', help='used to specify the strip of the plot which should be used')
parser.add_argument('-b', '--bucket', dest='bucket', default=[9999], nargs='*', type=int, help='used to specify the bucket of the plot which should be used | type=int')
parser.add_argument('-k', '--kpix', dest='kpix', default=[9999], nargs='*', type=int, help='used to specify the bucket of the plot which should be used | type=int')
parser.add_argument('-d', '--draw', dest='draw_option', default='', help='specify the drawing option as given by the root draw option, needs to be given as a single string (e.g. hist same or hist same multifile')
parser.add_argument('-o', '--output', dest='output_name', help='specifies the name and type of the output file (e.g. test.png, comparison.root etc...')
parser.add_argument('--refuse', dest='refuse', default= ['nothing'], nargs='*', help='add string that should be exluded from histogram search')
parser.add_argument('-r', '--rebin', dest='rebin', default=1, type = int, help='add number to rebin the histograms | type=int')
parser.add_argument('--name2', dest='name2', default=['everything'], nargs='*',  help='used to specify the name of the plot which should be used')
parser.add_argument('--refuse2', dest='refuse2', default= ['nothing'], nargs='*', help='add string that should be exluded from histogram search')
parser.add_argument('--exact', dest='exact', default=False, help='if set to True, only histograms with the exact name will be used')
parser.add_argument('--xrange', dest='xaxisrange', default=[9999], nargs='*', type=float, help='set a xrange for the plot to used with xmin xmax as the two arguments | type=float')
parser.add_argument('--yrange', dest='yaxisrange', default=[9999], nargs='*', type=float, help='set a yrange for the plot to used with ymin ymax as the two arguments | type=float')
parser.add_argument('--zrange', dest='zaxisrange', default=[9999], nargs='*', type=float, help='set a yrange for the plot to used with ymin ymax as the two arguments | type=float')
parser.add_argument('--legend', dest='legend', nargs='*', help='list of names to be used as legend titles instead of the default filename+histogram name')
parser.add_argument('--ylog', dest='ylog', help='if given as an option, set y axis to logarithmic. Remember to set the yrange to start above 0!')
parser.add_argument('--zlog', dest='zlog', help='if given as an option, set z axis to logarithmic.')
parser.add_argument('--color', dest='color', default=[60, 1, 418,  810, 402,  908, 435, 880, 860, 632, 840, 614], nargs='*', help='list of colors to be used')
#parser.add_argument('--color', dest='color', default=[590, 591, 593, 596, 600, 602, 604, 880, 860, 632, 840, 614], nargs='*', help='list of colors to be used')
parser.add_argument('--xtitle', dest='xtitle', help='choose the name of the x axis title')
parser.add_argument('--ytitle', dest='ytitle', help='choose the name of the y axis title')
parser.add_argument('--ztitle', dest='ztitle', help='choose the name of the z axis title')
parser.add_argument('--order', dest='order', nargs='+', type=int,  help='choose the order of plotting with same (to ensure no histograms overlap)')
parser.add_argument('-q', '--newdaq', dest='newdaq', default=True, help='give as a command when using files from the new daq to ensure filename check etc. are correct')
parser.add_argument('-l', dest='legendloc', nargs='+', type=float, default = [0.98, 0.99], help='first argument is the left x position of the legend box and second argument is the upper y position of the legend box')
parser.add_argument('-f', dest='folder', default='summer', help='tb is testbeam folder elab is elab folder. default is elab folder.')
parser.add_argument('--aratio', dest='aratio', nargs='+', type=float,  default=[1200,900], help='aspect ratio of the output file')
args = parser.parse_args()
if len(sys.argv) < 2:
	print parser.print_help()
	sys.exit(1)
print ''


if ('everything' in args.name2):
	args.name2 = args.name

#print args.refuse

teststring = 'ab_cds'

if any(name in teststring for name in args.name) and all(refuse not in teststring for refuse in args.refuse):
	print 'accepted'
else:
	print teststring

#legend_location = [0.65,0.65,0.98,0.85] # x_left, y_bottom, x_right, y_top
legend_location = [0.15,0.65,0.35,0.85] # x_left, y_bottom, x_right, y_top

##-----------------
##produce empty root file and filename lists.

root_file_list = []
filename_list = []
##-----------------
##loop through all given files and add them to the list. then loop through the keys for every file..
for root_file in args.file_in:
	root_file_list.append(ROOT.TFile(root_file))

	#if (args.newdaq is True):
		#print "HUH"
	if '/Run_' in root_file:
		filename_list.append(root_file[root_file.find('/Run_')+1:root_file.rfind('.dat')+1])
	elif '/Calibration_' in root_file:
		filename_list.append(root_file[root_file.find('/Calibration_')+1:root_file.rfind('.dat')+1])
	#else:
		#filename_list.append(root_file[root_file.find('/20')+1:root_file.rfind('.external')])
print filename_list
for x in root_file_list:
	key_root = x.GetListOfKeys()
	loopdir(key_root)
if ('elab' in args.folder):
	folder_loc = '/home/lycoris-admin/Documents/elab201904/'
elif ('tb' in args.folder):
	folder_loc = '/home/lycoris-admin/Documents/testbeamPCMAG201904_05/'
elif ('summer' in args.folder):
	folder_loc = '/home/lycoris-admin/Documents/summer_student/'



print 'Looking for histograms/graphs'
print '----------------------'
print 'Name contains ', args.name
print 'Channel [9999 = everything] ',args.channel
print 'Strip [9999 = everything] ',args.strip
print 'Bucket [9999 = everything; 4 = only total] ',args.bucket
print 'KPiX [9999 = everything] ',args.kpix
print 'Refusing the following ', args.refuse
if (args.name2 is not args.name):
	print 'Name also contains ', args.name2


print 'Number of histograms found is: ', len(hist_list)
print hist_list
print 'Number of graphs found is: ', len(graph_list)
print graph_list
if (args.ylog and args.yaxisrange[0] is 0):
	print 'Setting y axis to log, only works if the range was specified to start at y_min > 0'
##------------------
##start of the plotting.

if ((len(hist_list) > len(args.color) or len(graph_list) > len(args.color))) and ("same" in args.draw_option):
	print 'You do not have enough colors', len(args.color), 'for the number of histograms you have', len(hist_list)
	sys.exit(1)
if (len(hist_list) is not 0):
	##------------------
	##initialize a canvas, a stack histogram and further variables.
	drawing_option = args.draw_option
	c1 = ROOT.TCanvas( args.output_name, 'Test', args.aratio[0], args.aratio[1] )
	c1.cd()
	#c1.SetFillColor(0)

	statBoxW = 0.15
	statBoxH = 0.105
	legend = ROOT.TLegend(args.legendloc[0], args.legendloc[1], args.legendloc[0]+statBoxW, args.legendloc[1]+statBoxH)
	hist_comp = ROOT.TH2F("test", "test", 32, -0.5, 31.5, 32, -0.5, 31.5)
	counter = 1
	x_title = None
	y_title = None
	z_title = None
	x_low = None
	x_high = None
	y_low = None
	y_high = None
	z_low = None
	z_high = None
	new_hist_list = []
	legendname = []
	z1 = []
	if args.order:
		for i in args.order:
			new_hist_list.append(hist_list[i])
	else:
		new_hist_list = hist_list
	for histogram in new_hist_list:
		obj = histogram.ReadObj()
		print 'Number of total entries = ', '%.2E' % Decimal(obj.GetEntries())
		xbins = obj.GetNbinsX()
		ybins = obj.GetNbinsY()
		c = 0
		if counter is 1:
			for  i in xrange(xbins):
				for j in xrange(ybins):
					z1.append(obj.GetBinContent(i,j))
					c+=1
		else:
			for  i in xrange(xbins):
				for j in xrange(ybins):
					z1[c]=z1[c]-obj.GetBinContent(i,j)
					hist_comp.SetBinContent(i,j,z1[c])
					c+=1
		counter+=1
	print z1
	hist_comp.Draw(drawing_option)
	if (args.xtitle):
		x_title = args.xtitle
	else:
		x_title = x_axis.GetTitle()
	if (args.ytitle):
		y_title = args.ytitle
	else:
		y_title = y_axis.GetTitle()
	if (args.ztitle):
		z_title = args.ztitle
	else:
		z_title = z_axis.GetTitle()

	xaxis = hist_comp.GetXaxis()
	xaxis.SetTitle(x_title)

	yaxis = hist_comp.GetYaxis()
	if 9999 not in args.yaxisrange:
		yaxis.SetRangeUser(y_low, y_high)
		print 'test'
	yaxis.SetTitle(y_title)

	zaxis = hist_comp.GetZaxis()
	if 9999 not in args.zaxisrange:
		z_low = args.zaxisrange[0]
		z_high = args.zaxisrange[1]
		zaxis.SetRangeUser(z_low, z_high)
		print 'test'
	zaxis.SetTitle(z_title)


	c1.Update()
	ROOT.TGaxis.SetMaxDigits(3)
	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetOptFit(0)
	c1.Modified()
	c1.Update()
	run_name = filename_list[0][:-1]
	if (args.output_name):
		outname = folder_loc+run_name+'_'+args.output_name
		print 'Creating '+outname
		#c1.SaveAs(outname+'.svg')
		c1.SaveAs(outname+'.png')
	else:
		outname = folder_loc+run_name+'_'+graph.GetName()
		print 'Creating '+outname+'.pvg'
		#c1.SaveAs(outname+'.svg')
		c1.SaveAs(outname+'.png')
	c1.Close()

	    ##------------------
		##loop through the histograms, get all parameters
#		obj = histogram.ReadObj()
#		print 'Number of total entries = ', '%.2E' % Decimal(obj.GetEntries())
#		x_axis = obj.GetXaxis()
#		y_axis = obj.GetYaxis()
#		##------------------
#		##rebin the histogram
#		if (args.rebin is not 1):
#			    obj.Rebin(args.rebin)
#		##------------------
#		##adjust the xrange
#		if 9999 in args.xaxisrange:
#			if (x_low is None):
#				x_low = obj.FindFirstBinAbove(0)-10
#			elif (x_low > obj.FindFirstBinAbove(0)-10):
#				x_low = obj.FindFirstBinAbove(0)-10
#			if (x_high is None):
#				x_high = obj.FindLastBinAbove(0)+10
#			elif (x_high < obj.FindLastBinAbove(0)+10):
#				x_high = obj.FindLastBinAbove(0)+10
#			if (x_high > obj.GetNbinsX()):  #avoids overflow bin
#			    x_high = obj.GetNbinsX()
#			if (x_low <= 0): #avoids underflow bin
#			    x_low = 1
#			x_axis.SetRange(x_low, x_high)
#		else:
#			x_low = args.xaxisrange[0]
#			x_high = args.xaxisrange[1]
#			x_axis.SetRangeUser(x_low, x_high)
#			print 'Number of normed Entries in range = ', obj.Integral()
#			print 'Number of unweighted Entries in range = ', '%.2E' % Decimal(obj.Integral() * obj.GetEntries())

#		if 9999 not in args.yaxisrange:
#			y_low = args.yaxisrange[0]
#			y_high = args.yaxisrange[1]
#			hist_comp.SetMaximum(y_high)
#			hist_comp.SetMinimum(y_low)
#		#print x_low, x_high
#		#x_axis.SetRangeUser(x_low, x_high)
#		obj.SetLineColor(args.color[counter-1])
#		obj.SetMarkerColor(args.color[counter-1]+3)
#		obj.SetFillColor(args.color[counter-1])

#		##------------------
#		##draw histograms into the same canvas (equivalent to option same)
#		hist_comp.Add(obj, "sames")


#		##------------------
#		##adjust legend and the x and y title name if chosen
#		if (not args.legend):
#			if len(filename_list) > 1:
#				legend.AddEntry(obj, filename_list[counter-1]+'_'+histogram.GetName())
#			else:
#				legend.AddEntry(obj, '_'+histogram.GetName())
#		else:
#			legend.AddEntry(obj, args.legend[counter-1])
#			legendname.append(args.legend[counter-1])
#		counter +=1
#		if (args.xtitle):
#			x_title = args.xtitle
#		else:
#			x_title = x_axis.GetTitle()
#		if (args.ytitle):
#			y_title = args.ytitle
#		else:
#			y_title = y_axis.GetTitle()



    ##------------------
	##set y axis to log

#	if args.ylog:
#		c1.SetLogy()
#		ROOT.gPad.SetLogy()
#	if args.zlog:
#		c1.SetLogz()
#		ROOT.gPad.SetLogz()
#	##------------------
#	##draw histogram + components and save the file

#

#

else:
	print 'There are NO valid histograms/graphs in the current selection'
	print ''
	print ''
for x in root_file_list:
	ROOT.gROOT.GetListOfFiles().Remove(x)



#raw_input('Press Enter to look at the next histogram')






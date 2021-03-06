#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import sys
import numpy as np
from decimal import Decimal
from array import array

def loopdir_new(keys, all_names, refuse_names, title=False):  # loop through all subdirectories of the root file and add all fitting object to a list
	object_list = []
	for key_object in keys:
		if ('TDirectory' in key_object.GetClassName()):
			object_list= object_list+loopdir_new(key_object.ReadObj().GetListOfKeys(), all_names, refuse_names)
		else:
			if ('everything' in all_names):
				object_list.append(key_object)
			else:
				if (title): #if titlesearch is true, will search wtihin histogram title instead of object name
					obj = key_object.ReadObj()
					if (None in refuse_names):
						if (all(name in obj.GetTitle() for name in all_names)):
							object_list.append(key_object)
							print(obj.GetTitle())
					elif (all(name in obj.GetTitle() for name in all_names) and  all(name not in obj.GetTitle() for name in refuse_names)):
						object_list.append(key_object)
						print(obj.GetTitle())
						key_object.Print()
				else:
					if (None in refuse_names):
						if (all(name in key_object.GetName() for name in all_names)):
							object_list.append(key_object)
							print(key_object.GetName())
					elif (all(name in key_object.GetName() for name in all_names) and  all(name not in key_object.GetName() for name in refuse_names)):
						object_list.append(key_object)
						print(key_object.GetName())
						key_object.Print()
	return object_list

def tree_to_hist(tree, conditions, variables, bin_range, name, keyCounter, norm=False):
	print("Number of variables ", len(variables))
	print("Binrange ", bin_range)
	if (len(variables) == 1):
		histName = name#+"_"+str(keyCounter)
		histlabels = histName+";"+variables[0]+";Entries"
		myh = ROOT.TH1F(histName, histlabels, int(bin_range[0]), float(bin_range[1]), float(bin_range[2]))
		draw = variables[0]+">>"+histName
		print("Histogram labels ", histlabels)
		print(draw)
	elif (len(variables) == 2):
		histName = name#+"_"+str(keyCounter)
		histlabels = histName+";"+variables[1]+";"+variables[0]
		myh = ROOT.TH2F(histName, histlabels, int(bin_range[0]), float(bin_range[1]), float(bin_range[2]), int(bin_range[3]), float(bin_range[4]), float(bin_range[5]))
		draw = variables[0]+":"+variables[1]+">>"+histName
		print("Histogram labels ", histlabels)
		print(draw)
	else:
		print("Someone or something fucked up in the variables")
		sys.exit()
	print("Conditions ", conditions)
	tree.Draw(draw, conditions)
	if norm:
		scaling = 1.0/myh.GetEntries()
		myh.Scale(scaling)
	return myh

def plot_tree(keys, conditions, variables, bin_range, name):
	print("Number of objects = ", len(keys))
	myHistograms = []
	for keyCounter, key_object in enumerate(keys):
		if ('TTree' in key_object.GetClassName()):
			inputTree = key_object.ReadObj()
			print("Tree Branches ", inputTree.Print())
			#outhist = ROOT.TH1F("myh", "myh;layers;entries", 6,10,16)
			myh = tree_to_hist(inputTree, conditions, variables, bin_range, name, keyCounter)
			myHistograms.append(myh)
	return myHistograms

def drawSame(hists, drawOption, legendName,  MarkerStyle, name, ylog, legendLoc, order=None):
	c1 = ROOT.TCanvas( 'Canvas', 'Canvas', 2000, 1500 )
	c1.cd()
	statBoxW = 0.1
	nColumns = 1.
	statBoxH = 0.05*len(hists)/nColumns
	if order:
		statBoxH = 0.05*len(order)/nColumns
	if ("L" in drawOption):
		lines = True
	else:
		lines = False
	drawOption = drawOption+"PLC PMC" #to use palette colors
	if legendName is not None:
		print("Legendnames ", legendName)
		legend = ROOT.TLegend(legendLoc[0], legendLoc[1]-statBoxH, legendLoc[0]+statBoxW, legendLoc[1])
		legend.SetNColumns(int(nColumns))
		legend.SetFillStyle(0)
	new_legendlist = []
	new_hist_list = []

	legendname = []
	if order:
		for i in order:
			new_hist_list.append(hists[i])
			#new_legendlist.append(legendName[i])
	else:
		new_hist_list = hists
		#new_legendlist = legendName
	new_legendlist = legendName  
	for counter, h in enumerate(new_hist_list):
		#print 'Number of total entries = ', '%.2E' % Decimal(h.GetEntries())
		print('Mean value = ', '%.3E' % Decimal(h.GetMean()))
		print('RMS = ', '%.3E' % Decimal(h.GetRMS()))
		print('Mean Error value = ', '%.3E' % Decimal(h.GetMeanError()))
		print('RMS Error = ', '%.3E' % Decimal(h.GetRMSError()))
		q = np.array([0.])
		prob = np.array([0.5])
		median = h.GetQuantiles(1, q, prob)
		print("Median value = ", q)
		print("")
		if("Graph" not in h.GetName()):
			print('Integral =', '%.3E' % Decimal(h.Integral()))
			print('Number of total entries = ', '%.3E' % Decimal(h.GetEntries()))
		x_axis = h.GetXaxis()
		y_axis = h.GetYaxis()
		if new_legendlist is not None:
			legEntry = legend.AddEntry(h, new_legendlist[counter])
		h.SetMarkerStyle(MarkerStyle[0][counter%6])
		h.SetMarkerSize(MarkerStyle[1][counter%6])
		h.SetLineWidth(4)
		fit=h.GetFunction("gaus")
		if fit: #if the fitpointer is not NULL
			fit.SetLineWidth(5)
		h.Draw(drawOption)
		if (lines):
			print("Drawing with connecting lines")
			h.Draw("l hist same")

		if legendName is not None:
			legend.Draw()
	if ylog:
		c1.SetLogy()
	c1.Modified()
	c1.Update()
	saveFile(c1, [""], 0, '/home/kraemeru/plots/2021/', str(name))


def drawGraph(hists, drawOption, legendName,  MarkerStyle, name, ylog, legendLoc='right', order=None):
	c1 = ROOT.TCanvas( 'Canvas', 'Canvas', 2000, 1500 )
	c1.cd()
	statBoxW = 0.1
	statBoxH = 0.05*len(hists)
	if ("L" in drawOption):
		lines = True
	else:
		lines = False
	if len(hists) > 1:
		drawOption = drawOption+" PLC PMC" #to use palette colors
	if legendName is not None:
		print("Legendnames ", legendName)
		legend = ROOT.TLegend(legendLoc[0], legendLoc[1]-statBoxH, legendLoc[0]+statBoxW, legendLoc[1])
		legend.SetFillStyle(0)
	new_legendlist = []
	new_hist_list = []
	multiGraph = ROOT.TMultiGraph()
	legendname = []
	yOld = [None, None]
	xOld = [None, None]
	if order:
		for i in order:
			new_hist_list.append(hists[i])
			new_legendlist.append(legendName[i])
	else:
		new_hist_list = hists
		new_legendlist = legendName
	cols = ROOT.TColor.GetPalette()
	for counter, h in enumerate(new_hist_list):
		#print 'Number of total entries = ', '%.2E' % Decimal(h.GetEntries())
		multiGraph.GetYaxis().SetTitle(h.GetYaxis().GetTitle())
		multiGraph.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
		print('Mean value = ', '%.3E' % Decimal(h.GetMean()))
		print('RMS = ', '%.3E' % Decimal(h.GetRMS()))
		xMin = ROOT.TMath.MinElement(h.GetN(),h.GetX())/1.01
		yMin = ROOT.TMath.MinElement(h.GetN(),h.GetY())/1.01
		xMax = ROOT.TMath.MaxElement(h.GetN(),h.GetX())*1.01
		yMax = ROOT.TMath.MaxElement(h.GetN(),h.GetY())*1.01
		if (None not in yOld):
			yMin = min(yOld[0],yMin)
			xMin = min(xOld[0],xMin)
			yMax = max(yOld[1],yMax)
			xMax = max(xOld[1],xMax)
		if("Graph" not in h.GetName()):
			print('Integral =', '%.3E' % Decimal(h.Integral()))
			print('Number of total entries = ', '%.3E' % Decimal(h.GetEntries()))
		x_axis = h.GetXaxis()
		y_axis = h.GetYaxis()
		if new_legendlist is not None:
			legEntry = legend.AddEntry(h, new_legendlist[counter])
		h.SetMarkerStyle(MarkerStyle[0][counter%6])
		h.SetMarkerSize(MarkerStyle[1][counter%6])
		multiGraph.Add(h)
		yOld = [yMin, yMax]
		xOld = [xMin, xMax]
        
	if ylog:
		c1.SetLogy()
	print("HEEEEEEEEEEEEEY ", multiGraph.GetListOfGraphs().GetEntries())
	multiGraph.Draw(drawOption)
	c1.Update()
	c1.Modified()
	multiGraph.GetXaxis().SetLimits(xMin,xMax)
	multiGraph.SetMaximum(yMax)
	multiGraph.SetMinimum(yMin)
	for gr in multiGraph:
		markerColor = gr.GetMarkerColor()
		gr.GetFunction('pol1').SetLineColor(markerColor)
		gr.GetFunction('pol1').SetLineStyle(9)
		gr.GetFunction('pol1').SetLineWidth(3)
	if legendName is not None:
		legend.Draw()

	c1.Modified()
	c1.Update()
	saveFile(c1, [""], 0, '/home/kraemeru/plots/2021/', str(name))


def saveFile(c1, filename_list, counter, folder_loc, outName):
	if (len(filename_list) == 1):
		run_name = filename_list[0][:-1]
	else:
		run_name = filename_list[counter][:-1]
	#c1.SaveAs(folder_loc+"svg/"+run_name+outName+'.svg')
	#c1.SaveAs(folder_loc+"eps/"+run_name+outName+'.eps')
	c1.SaveAs(folder_loc+"png/"+run_name+outName+'.png')
	print(folder_loc)
	#c1.SaveAs(folder_loc+"C/"+run_name+outName+'.C')

def myROOTStyle(markerScale, sbox, fbox, bSize, legendLoc=[0.4,0.9], nobox=True):
	mystyle = ROOT.TStyle("mystyle", "My Style")


	mystyle.SetPaintTextFormat("5.2f");
	## My adapted viridis (removed the yellow)
	red   = [ 26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255.]#, 246./255.]
	green = [  9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255.]#, 222./255.]
	blue  = [ 30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.]#,   0./255.]
	#stops = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625 ,0.75 ,0.875,1.0]
	stops = [0.0, 0.15, 0.3, 0.45, 0.6, 0.75 ,0.9 ,1.0]
	stopsArray = array('d', stops)
	redArray = array('d', red)
	greenArray = array('d', green)
	blueArray = array('d', blue)
	Idx = ROOT.TColor.CreateGradientColorTable(8, stopsArray, redArray, greenArray, blueArray, 255);

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
	mystyle.SetLegendFillColor(0)
	##set legend text size etc.
	mystyle.SetLegendTextSize(0.04)
	#
	##use the primary color palette
	mystyle.SetPalette(Idx)#112)
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
	mystyle.SetPadBottomMargin(0.11)
	#mystyle.SetPadBottomMargin(0.06)
	mystyle.SetPadTopMargin(0.05)
	mystyle.SetPadRightMargin(0.15)
	mystyle.SetPadLeftMargin(0.13)
	#mystyle.SetPadLeftMargin(0.06)
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
	mystyle.SetStatFormat("5.3f")
	


	ROOT.TGaxis.SetMaxDigits(4)

	#mystyle.SetTitleBorderSize(0)
	mystyle.SetStatBorderSize(0)
	#mystyle.SetTextFont(42)


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
	mystyle.SetOptStat(sbox) ##removes stat box
	mystyle.SetOptFit(fbox)
	#mystyle.SetOptStat(1110)
	#mystyle.SetStatY(0.9)
	#mystyle.SetStatX(0.8)
	mystyle.SetStatY(legendLoc[1])
	mystyle.SetStatX(legendLoc[0])
	mystyle.SetStatH(bSize[1])
	mystyle.SetStatW(bSize[0])#was 0.3
	#mystyle.SetOptFit(0)
	#mystyle.SetOptStat(0000001) #only name
	#
	##marker settings
	mystyle.SetMarkerStyle(8)
	mystyle.SetMarkerSize(2.5)
	mystyle.SetLineWidth(3)

	#default_colors=[["#08306b", "#8c2d04",	"#08519c", "#d94801",	"#2171b5", "#f16913",	"#4292c6", "#fd8d3c","#6baed6",	"#9ecae1",	"#c6dbef",	"#deebf7",	"#f7fbff"],
	#["#08306b",	"#08519c", 	"#2171b5", 	"#4292c6", "#6baed6",	"#9ecae1",	"#c6dbef",	"#deebf7",	"#f7fbff"],
	#["#08306b", 	"#08519c", 	"#2171b5", 	"#4292c6", "#fd8d3c", "#8c2d04","#d94801","#f16913","#6baed6",	"#9ecae1",	"#c6dbef",	"#deebf7",	"#f7fbff"],
	#["#0d0887", "#5302a3","#8b0aa5","#b83289","#db5c68","#f48849","#febd2a","#f0f921"]]
	if (nobox):
		print("Not printing stats box")
		mystyle.SetOptStat(0)
		mystyle.SetOptFit(0)
	myMarker = [20, 21, 22, 23, 33, 29] #FullCircle #FullSquare #UpTriangle #DownTriangle #Diamond #Star
	myMarkerSize = [2.2, 2.0, 2.5, 2.5, 3.2, 3.2]
	myMarkerSize = [x*markerScale for x in myMarkerSize] #Scaling up the sizes
	return mystyle, myMarker, myMarkerSize

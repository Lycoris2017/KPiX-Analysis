#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import sys
from decimal import Decimal

def loopdir_new(keys, all_names):  # loop through all subdirectories of the root file and add all fitting object to a list
	object_list = []
	for key_object in keys:
		if ('TDirectory' in key_object.GetClassName()):
			object_list= object_list+loopdir_new(key_object.ReadObj().GetListOfKeys(), all_names)
		else:
			if ('everything' in all_names):
				object_list.append(key_object)
			elif (all(name in key_object.GetName() for name in all_names)):
				object_list.append(key_object)
				key_object.Print()
	return object_list

def tree_to_hist(tree, conditions, variables, bin_range, name, norm=True):
	print "Number of variables ", len(variables)
	print "Binrange ", bin_range
	if (len(variables) == 1):
		histlabels = name+";"+variables[0]+";Nr. of Entries"
		myh = ROOT.TH1F(name, histlabels, int(bin_range[0]), float(bin_range[1]), float(bin_range[2]))
		draw = variables[0]+">>"+name
		print draw
	elif (len(variables) == 2):
		histlabels = name+";"+variables[1]+";"+variables[0]
		myh = ROOT.TH2F(name, histlabels, int(bin_range[0]), float(bin_range[1]), float(bin_range[2]), int(bin_range[3]), float(bin_range[4]), float(bin_range[5]))
		draw = variables[0]+":"+variables[1]+">>"+name
		print "Histogram labels ", histlabels
		print draw
	else:
		print "Someone or something fucked up in the variables"
		sys.exit()
	print "Conditions ", conditions
	tree.Draw(draw, conditions)
	if norm:
		scaling = 1.0/myh.GetEntries()
		myh.Scale(scaling)
	return myh

def plot_tree(keys, conditions, variables, bin_range, drawOption, name, rangeX, rangeY, colors, legendName, order, ylog):
	print "Number of objects = ", len(keys)
	myHistograms = []
	for keyCounter, key_object in enumerate(keys):
		if ('TTree' in key_object.GetClassName()):
			inputTree = key_object.ReadObj()
			print "Tree Branches ", inputTree.Print()
			#outhist = ROOT.TH1F("myh", "myh;layers;entries", 6,10,16)
			c1 = ROOT.TCanvas( str(keyCounter), str(keyCounter), 2000, 1500 )
			c1.cd()
			myh = tree_to_hist(inputTree, conditions, variables, bin_range, name)
			myHistograms.append(myh)
	drawSame(myHistograms, drawOption, rangeX, rangeY, colors, legendName, order, ylog)
#	for histCounter, hist in enumerate(myHistograms):
#		if ("same" in drawOption):
#		myh.Draw(drawOption)
#		myh.GetXaxis().SetRangeUser(rangeX[0], rangeX[1])
#		if (9999 not in rangeY):
#			myh.GetYaxis().SetRangeUser(rangeY[0], rangeY[1])
#		c1.Update()
#		outname = name+".png"
#		c1.SaveAs(outname)

def drawSame(hists, drawOption, rangeX, rangeY, colors, legendName, order=None, ylog=False):
	c1 = ROOT.TCanvas( 'Canvas', 'Canvas', 2000, 1500 )
	c1.cd()
	statBoxW = 0.1
	statBoxH = 0.05*len(hists)
	print "Legendnames ", legendName
	legend = ROOT.TLegend(0.7, 0.93-statBoxH, 0.7+statBoxW, 0.93)
	new_hist_list = []
	new_legendlist = []
	legendname = []
	if order:
		for i in order:
			new_hist_list.append(hists[i])
			new_legendlist.append(legendName[i])
	else:
		new_hist_list = hists
		new_legendlist = legendName

	for counter, h in enumerate(new_hist_list):
		print 'Number of total entries = ', '%.2E' % Decimal(h.GetEntries())
		print 'Mean value = ', '%.2E' % Decimal(h.GetMean())
		print 'RMS = ', '%.2E' % Decimal(h.GetRMS())
		x_axis = h.GetXaxis()
		y_axis = h.GetYaxis()
		h.SetLineColor(ROOT.TColor.GetColor(colors[counter]))
		h.SetMarkerColor(ROOT.TColor.GetColor(colors[counter]))

		legEntry = legend.AddEntry(h, new_legendlist[counter])
		legEntry.SetTextColor(ROOT.TColor.GetColor(colors[counter]))

		h.Draw(drawOption)
		legend.Draw()
	if ylog:
		c1.SetLogy()
	c1.Modified()
	c1.Update()
	saveFile(c1, ["test"], 0, '/home/lycoris-dev/Documents/testbeam202003/', 'test')

def saveFile(c1, filename_list, counter, folder_loc, outName):
	if (len(filename_list) == 1):
		run_name = filename_list[0][:-1]
	else:
		run_name = filename_list[counter][:-1]
	saveName = folder_loc+run_name+'_'+outName
	print 'Creating '+saveName
	c1.SaveAs(saveName+'.svg')
	c1.SaveAs(saveName+'.eps')
	c1.SaveAs(saveName+'.png')
	c1.SaveAs(saveName+'.C')

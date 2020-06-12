#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import sys

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

def plot_tree(keys, conditions, variables, bin_range, drawOption, name, rangeX, rangeY):
	for key_object in keys:
		if ('TTree' in key_object.GetClassName()):
			inputTree = key_object.ReadObj()
			print "Tree Branches ", inputTree.Print()
			#outhist = ROOT.TH1F("myh", "myh;layers;entries", 6,10,16)
			c1 = ROOT.TCanvas( 'Test', 'Test', 1600, 1200 )

			c1.cd()
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
			print "Drawing Option ", drawOption
			print "Conditions ", conditions
			inputTree.Draw(draw, conditions)
			myh.Draw(drawOption)
			myh.GetXaxis().SetRangeUser(rangeX[0], rangeX[1])
			if (9999 not in rangeY):
				myh.GetYaxis().SetRangeUser(rangeY[0], rangeY[1])
			c1.Update()

			outname = name+".png"
			c1.SaveAs(outname)

#			corr_clusters->Draw("noise>>h(100,0,2)", "onTrack == 1", "h e")

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

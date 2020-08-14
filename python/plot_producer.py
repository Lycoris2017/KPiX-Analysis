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
from pyroot_functions import saveFile
from pyroot_functions import myROOTStyle
from array import array

ROOT.gROOT.SetBatch(True)

hist_list= []  #Global list of all chosen histograms (only works because its a list, with single variables they need to be declared as global in the function
graph_list = []




#def saveFile(c1, filename_list, inCounter, counter, folder_loc, outName):
#	if (inCounter == 1):
#		run_name = filename_list[0][:-1]
#	else:
#		run_name = filename_list[counter][:-1]
#	saveName = folder_loc+run_name+'_'+outName
#	print 'Creating '+saveName
#	c1.SaveAs(saveName+'.svg')
#	c1.SaveAs(saveName+'.eps')
#	c1.SaveAs(saveName+'.png')
#	c1.SaveAs(saveName+'.C')

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


def hist_plotter():
	if ('stack' in args.draw_option or 'nostack' in args.draw_option):
		##------------------
		##initialize a canvas, a stack histogram and further variables.
		if ('nostack' not in args.draw_option):
			f = args.draw_option.find('stack')
			fend = f+5 ## End of stack
			l = len(args.draw_option)
			drawing_option = args.draw_option[0:f]+args.draw_option[fend:l]+"PLC PMC" ##cutting out 'stack' from the string of drawing options
			print "Using draw options: ", drawing_option
		else:
			drawing_option = args.draw_option+"PLC PMC" #exchange the same with a NOSTACK as I am using THStack
			print "Using draw options: ", drawing_option
		c1 = ROOT.TCanvas( args.output_name, 'Test', args.aratio[0], args.aratio[1] )
		c1.cd()
		#c1.SetFillColor(0)
		statBoxW = 0.15
		statBoxH = 0.105
		legend = ROOT.TLegend(args.legendloc[0], args.legendloc[1], args.legendloc[0]+statBoxW, args.legendloc[1]+statBoxH)
		hist_comp = ROOT.THStack()
		counter = 0
		x_title = None
		y_title = None
		x_low = None
		x_high = None
		y_low = None
		y_high = None
		new_hist_list = []
		legendname = []
		if args.order:
			for i in args.order:
				new_hist_list.append(hist_list[i]) 
		else:
			new_hist_list = hist_list
		for histogram in new_hist_list:
			##------------------
			##loop through the histograms, get all parameters
			obj = histogram.ReadObj()
			print 'Number of total entries = ', '%.2E' % Decimal(obj.GetEntries())
			if args.norm:
				norm = 1.0/obj.GetEntries()
				obj.Scale(norm)
			x_axis = obj.GetXaxis()
			y_axis = obj.GetYaxis()
			##------------------
			##rebin the histogram
			if (args.rebin is not 1):
					obj.Rebin(args.rebin)
			##------------------
			##adjust the xrange
			if 9999 in args.xaxisrange:
				if (x_low is None):
					x_low = obj.FindFirstBinAbove(0)-10
				elif (x_low > obj.FindFirstBinAbove(0)-10):
					x_low = obj.FindFirstBinAbove(0)-10
				if (x_high is None):
					x_high = obj.FindLastBinAbove(0)+10
				elif (x_high < obj.FindLastBinAbove(0)+10):
					x_high = obj.FindLastBinAbove(0)+10
				if (x_high > obj.GetNbinsX()):  #avoids overflow bin
					x_high = obj.GetNbinsX()
				if (x_low <= 0): #avoids underflow bin
					x_low = 1
				x_axis.SetRange(x_low, x_high)
			else:
				x_low = args.xaxisrange[0]
				x_high = args.xaxisrange[1]
				x_axis.SetRangeUser(x_low, x_high)
				print 'Number of normed Entries in range = ', obj.Integral()
				print 'Number of unweighted Entries in range = ', '%.2E' % Decimal(obj.Integral() * obj.GetEntries())

			if 9999 not in args.yaxisrange:
				y_low = args.yaxisrange[0]
				y_high = args.yaxisrange[1]
				hist_comp.SetMaximum(y_high)
				hist_comp.SetMinimum(y_low)
			#print x_low, x_high
			#x_axis.SetRangeUser(x_low, x_high)
			#obj.SetLineColor(ROOT.TColor.GetColor(args.color[counter]))
			#obj.SetMarkerColor(ROOT.TColor.GetColor(args.color[counter]))
			obj.SetMarkerStyle(myMarker[counter])
			obj.SetMarkerSize(myMarkerSize[counter])
			#if (args.fill):
			#	obj.SetFillColor(ROOT.TColor.GetColor(args.color[counter]))
			#obj.GetFunction('gaus').SetLineColor(ROOT.TColor.GetColor(args.color[counter]))

			##------------------
			##draw histograms into the same canvas (equivalent to option same)
			hist_comp.Add(obj, "sames")
			
			
			##------------------
			##adjust legend and the x and y title name if chosen
			if (not args.legend):
				if len(filename_list) > 1:
					legend.AddEntry(obj, filename_list[counter]+'_'+histogram.GetName())
					legendname.append(filename_list[counter]+'_'+histogram.GetName())
				else:
					legend.AddEntry(obj, histogram.GetName())
					legendname.append(histogram.GetName())
			else:
				legend.AddEntry(obj, args.legend[counter])
				legendname.append(args.legend[counter])
			if (args.xtitle):
				x_title = args.xtitle
			else:
				x_title = x_axis.GetTitle()
			if (args.ytitle):
				y_title = args.ytitle
			else:
				y_title = y_axis.GetTitle()
			counter +=1
#			ROOT.TGaxis.SetMaxDigits(3)
			
		##------------------
		##set y axis to log
		
		if args.ylog:
			c1.SetLogy()
			ROOT.gPad.SetLogy()
		if args.zlog:
			c1.SetLogz()
			ROOT.gPad.SetLogz()
		##------------------
		##draw histogram + components and save the file
	
		hist_comp.Draw(drawing_option)
		c1.Update()
	
		xaxis = hist_comp.GetXaxis()
		xaxis.SetTitle(x_title)
		print x_low
		print x_high
		if 9999 in args.xaxisrange:
			xaxis.SetRange(x_low, x_high)
		else:
			xaxis.SetRangeUser(x_low, x_high)
		yaxis = hist_comp.GetYaxis()
		if 9999 not in args.yaxisrange:
			yaxis.SetRangeUser(y_low, y_high) 
			print 'test'
		yaxis.SetTitle(y_title)
		histListForStats = []
		for j in range(len(new_hist_list)-1, -1, -1): 
			histListForStats.append(hist_comp.GetStack().At(j) )
		print "Make sure the stat box is not set to stats 0 as otherwise it gives an error"
		arrangeStats(histListForStats, statBoxW, statBoxH, legendname)
		ROOT.gStyle.SetOptStat(0)
		ROOT.gStyle.SetOptFit(0)
		c1.Modified()
		c1.Update()
		if (args.output_name):
			canvasName = args.output_name
		else:
			canvasName = histogram.GetName()
		saveFile(c1, filename_list, counter-1, folder_loc, canvasName)
		c1.Close()
	elif ('same' in args.draw_option):
		##------------------
		##initialize a canvas, a stack histogram and further variables.
		drawing_option = args.draw_option+"PLC PMC"
		c1 = ROOT.TCanvas( args.output_name, 'Canvas', args.aratio[0], args.aratio[1] )
		c1.cd()
		#c1.SetGrid(1)
		#c1.SetFillColor(0)
		statBoxW = 0.1
		statBoxH = 0.05*len(hist_list)
		print args.legendloc[0]
		legend = ROOT.TLegend(0.7, 0.93-statBoxH, 0.7+statBoxW, 0.93)
		#counter = 0
		x_title = None
		y_title = None
		x_low = None
		x_high = None
		yMin = None
		yMax = None
		new_hist_list = []
		new_legendlist = []
		legendname = []
		if args.order:
			for i in args.order:
				new_hist_list.append(hist_list[i]) 
				if args.legend:
					new_legendlist.append(args.legend[i])
		else:
			new_hist_list = hist_list
			if args.legend:
				new_legendlist = args.legend
		for counter, histogram in enumerate(new_hist_list):
			##------------------
			##loop through the histograms, get all parameters
			obj = histogram.ReadObj()
			if args.norm:
				norm = 1.0/obj.GetEntries()
				obj.Scale(norm)
			print 'Number of total entries = ', '%.2E' % Decimal(obj.GetEntries())
			print 'Mean value = ', '%.2E' % Decimal(obj.GetMean())
			print 'RMS = ', '%.2E' % Decimal(obj.GetRMS())
			x_axis = obj.GetXaxis()
			y_axis = obj.GetYaxis()
			##------------------
			##rebin the histogram
			if (args.rebin is not 1):
					obj.Rebin(args.rebin)
			##------------------
			##adjust the xrange
			if 9999 in args.xaxisrange:
				if (x_low is None):
					x_low = obj.FindFirstBinAbove(0)-10
				elif (x_low > obj.FindFirstBinAbove(0)-10):
					x_low = obj.FindFirstBinAbove(0)-10
				if (x_high is None):
					x_high = obj.FindLastBinAbove(0)+10
				elif (x_high < obj.FindLastBinAbove(0)+10):
					x_high = obj.FindLastBinAbove(0)+10
				if (x_high > obj.GetNbinsX()):  #avoids overflow bin
					x_high = obj.GetNbinsX()
				if (x_low <= 0): #avoids underflow bin
					x_low = 1
				x_axis.SetRange(x_low, x_high)
			else:
				x_low = args.xaxisrange[0]
				x_high = args.xaxisrange[1]
				x_axis.SetRangeUser(x_low, x_high)
				print 'Number of normed Entries in range = ', obj.Integral()
				print 'Number of unweighted Entries in range = ', '%.2E' % Decimal(obj.Integral() * obj.GetEntries())

			if args.yaxisrange:
				yMin = args.yaxisrange[0]
				yMax = args.yaxisrange[1]
			else:
				if (yMax is None):
					yMax = obj.GetMaximum()
				elif(yMax < yRangeScale*obj.GetMaximum()):
					yMax =  yRangeScale*obj.GetMaximum()
				if (args.ylog):
					if (yMin is None):
						yMin = 1.0/(yRangeScale*obj.GetEntries())
					elif (yMin > 1.0/(yRangeScale*obj.GetEntries())):
						yMin = 1.0/(yRangeScale*obj.GetEntries())
				else:
					if (yMin is None):
						yMin = (1./yRangeScale)*obj.GetMinimum()
					elif (yMin >  (1./yRangeScale)*obj.GetMinimum()):
						yMin =  (1.0/yRangeScale)*obj.GetMinimum()
			#print x_low, x_high
			#x_axis.SetRangeUser(x_low, x_high)
			#obj.SetLineColor(ROOT.TColor.GetColor(default_colors[args.color][counter]))
			#obj.SetMarkerColor(ROOT.TColor.GetColor(default_colors[args.color][counter]))
			obj.SetMarkerStyle(myMarker[counter])
			obj.SetMarkerSize(myMarkerSize[counter])
			#if (args.fill):
			#	obj.SetFillColor(ROOT.TColor.GetColor(default_colors[args.color][counter]))
			
			if (args.xtitle):
				x_title = args.xtitle
			else:
				x_title = x_axis.GetTitle()
			if (args.ytitle):
				y_title = args.ytitle
			else:
				y_title = y_axis.GetTitle()
			
			xaxis = obj.GetXaxis()
			xaxis.SetTitle(x_title)
			
			print x_low
			print x_high
			if 9999 in args.xaxisrange:
				xaxis.SetRange(x_low, x_high)
			else:
				xaxis.SetRangeUser(x_low, x_high)
			yaxis = obj.GetYaxis()
			yaxis.SetRangeUser(yMin, yMax)
			yaxis.SetTitle(y_title)
#			ROOT.TGaxis.SetMaxDigits(3)
			##------------------
			##draw histograms into the same canvas (equivalent to option same)
			
			
			
			##------------------
			##adjust legend and the x and y title name if chosen
			if (not args.legend):
				if len(filename_list) > 1:
					legEntry = legend.AddEntry(obj, filename_list[counter]+'_'+histogram.GetName())

				else:
					legEntry = legend.AddEntry(obj, '_'+histogram.GetName())
			else:
				legEntry = legend.AddEntry(obj, new_legendlist[counter])

				legendname.append(new_legendlist[counter])
			#counter +=1
			
			obj.Draw(drawing_option)
			legend.Draw()
			legEntry.SetTextColor(obj.GetLineColor())
			
		##------------------
		##set y axis to log
		
		if args.ylog:
			c1.SetLogy()
			ROOT.gPad.SetLogy()
		if args.zlog:
			c1.SetLogz()
			ROOT.gPad.SetLogz()
		##------------------
		##draw histogram + components and save the file
		c1.Update()
	
		
		
		ROOT.gStyle.SetOptStat(0)
		ROOT.gStyle.SetOptFit(0)
		
		c1.Modified()
		c1.Update()
		if (args.output_name):
			canvasName = args.output_name
		else:
			canvasName = histogram.GetName()
		saveFile(c1, filename_list, counter, folder_loc, canvasName)
		c1.Close()
	else:
		counter = 0
		file_counter = 0
		drawing_option = args.draw_option+"PLC PMC"
		for histogram in hist_list:
			##------------------
			##loop through the histograms, get all parameters and adjust the xrange
			#ROOT.gROOT.SetBatch(1)
			
			x_title = None
			y_title = None
			x_low = None
			x_high = None
			y_low = None
			y_high = None
			c1 = ROOT.TCanvas( 'test', 'Test', args.aratio[0], args.aratio[1] ) #
			obj = histogram.ReadObj()
			if args.norm:
				norm = 1.0/obj.GetEntries()
				obj.Scale(norm)
			print 'Number of total entries = ', '%.2E' % Decimal(obj.GetEntries())
			x_axis = obj.GetXaxis()
			y_axis = obj.GetYaxis()
			z_axis = obj.GetZaxis()


#			ROOT.TGaxis.SetMaxDigits(3)
			c1.cd()
			#print obj.GetEntries()
			##------------------
			##rebin the histogram
			if (args.rebin is not 1):
				obj.Rebin(args.rebin)
			
			if 9999 in args.xaxisrange:
				x_high = obj.FindLastBinAbove(0)+10
				if (x_high > obj.GetNbinsX()):  #avoids overflow bin
					x_high = obj.GetNbinsX()
				x_low = obj.FindFirstBinAbove(0)-10
				if (x_low <= 0): #avoids underflow bin
					x_low = 1
				x_axis.SetRange(x_low, x_high)
			else:
				x_low = args.xaxisrange[0]
				x_high = args.xaxisrange[1]
				x_axis.SetRangeUser(x_low, x_high)
				print 'Number of normed Entries in range = ', obj.Integral()
				print 'Number of unweighted Entries in range = ', '%.2E' % Decimal(obj.Integral() * obj.GetEntries())
			if  args.yaxisrange:
				y_low = args.yaxisrange[0]
				y_high = args.yaxisrange[1]
				y_axis.SetRangeUser(y_low, y_high)
#			if  args.zaxisrange:
#				z_low = args.zaxisrange[0]
#				z_high = args.zaxisrange[1]
#				z_axis.SetRangeUser(z_low, z_high)

			##obj.SetLineColor(ROOT.TColor.GetColor(default_colors[args.color][counter]))
			#obj.SetMarkerColor(ROOT.TColor.GetColor(default_colors[args.color][counter]))
			obj.SetMarkerStyle(myMarker[counter])
			obj.SetMarkerSize(myMarkerSize[counter])
			#if (args.fill):
				#obj.SetFillColor(ROOT.TColor.GetColor(default_colors[args.color][counter]))
			
			##------------------
			##set y axis to log
			if args.ylog:
				c1.SetLogy()
				ROOT.gPad.SetLogy()
			if args.zlog:
				c1.SetLogz()
				ROOT.gPad.SetLogz()
			##------------------
			##set axis titles if chosen through options.
			if (args.xtitle):
				x_axis.SetTitle(args.xtitle)
			if (args.ytitle):
				y_axis.SetTitle(args.ytitle)

			if (args.fit):
				print args.fit
				f1 = ROOT.TF1("f1", str(args.fit[0]), float(args.fit[1]), float(args.fit[2]))
				obj.Fit("f1", str(args.fit[0]), "", float(args.fit[1]), float(args.fit[2]) )
			obj.Draw(drawing_option)

#			if args.upperXaxis:
#				print "adding an upper x axis"
#				print x_low, y_low, x_high, y_high
#				axis = ROOT.TGaxis(x_low, y_high, x_high, y_high, x_low, x_high, 016, "-");
#				axis.SetLineColor(4);
#				axis.SetLabelFont(62);
#				axis.SetLabelSize(0.04);
#				axis.SetLabelColor(4);
#				axis.SetTitleFont(62);
#				axis.SetTitleSize(0.056);
#				axis.SetTitleOffset(1.3);
#				axis.SetTitleColor(4);
#				axis.SetTitle("electrons");
#				axis.ChangeLabel(-1, -1, -1, -1, -1, -1, "99864");
#				axis.ChangeLabel(1, -1, -1, -1, -1, -1, "0");
#				axis.ChangeLabel(2, -1, -1, -1, -1, -1, "12483");
#				axis.ChangeLabel(3, -1, -1, -1, -1, -1, "24966");
#				axis.ChangeLabel(4, -1, -1, -1, -1, -1, "37449");
#				axis.ChangeLabel(5, -1, -1, -1, -1, -1, "49932");
#				axis.ChangeLabel(6, -1, -1, -1, -1, -1, "62415");
#				axis.ChangeLabel(7, -1, -1, -1, -1, -1, "74898");
#				axis.ChangeLabel(8, -1, -1, -1, -1, -1, "87381");
#				axis.Draw()
#				c1.SetTopMargin(0.15)
			c1.Update()
			if args.nofit:
				obj.GetFunction('gaus').SetRange(0.0, 0.0, 0.01, 0.01) #IMPORTANT: workaround to remove fit from the plot
			if (args.output_name):
				canvasName = args.output_name
			else:
				canvasName = histogram.GetName()
			saveFile(c1, filename_list, file_counter, folder_loc, canvasName)
			c1.Close()
			counter= counter+1
			file_counter = counter/(len(hist_list)/len(filename_list))

			#for i in xrange(obj.FindFirstBinAbove(0),obj.FindLastBinAbove(0),1):
			#'	print obj.GetBinContent(i)							
	
def graph_plotter():
	if ('same' in args.draw_option):
		##------------------
		##initialize a canvas, a stack histogram and further variables.
		drawing_option = args.draw_option.replace('same', 'A')+"PLC PMC"
		if (args.output_name):
			c1 = ROOT.TCanvas( args.output_name, 'Test', args.aratio[0], args.aratio[1] )
		else:
			c1 = ROOT.TCanvas( 'test', 'Test', args.aratio[0], args.aratio[1] ) #
		c1.cd()
		c1.SetFillColor(0)
		statBoxW = 0.15
		statBoxH = 0.28
		legend = ROOT.TLegend(args.legendloc[0], args.legendloc[1], args.legendloc[0]+statBoxW, args.legendloc[1]+statBoxH)
		multi_graph = ROOT.TMultiGraph()
		counter = 1
		x_title = None
		y_title = None
		x_low = None
		x_high = None
		y_low = None
		y_high = None
		
		new_graph_list = []
		if args.order:
			for i in args.order:
				new_graph_list.append(graph_list[i]) 
		else:
			new_graph_list = graph_list
		
		for graph in new_graph_list:
			##------------------
			##loop through the histograms, get all parameters
			obj = graph.ReadObj()
			#print 'Number of entries =', obj.GetEntries()
			x_axis = obj.GetXaxis()
			y_axis = obj.GetYaxis()
			
			#chi_andreas_test = obj.GetFunction('pol1').GetChisquare()
			#print "CHI ANDREAS TEST = ", chi_andreas_test
			if args.nofit:
				obj.GetFunction('pol1').SetRange(0.0, 0.0, 0.01, 0.01) #IMPORTANT: workaround to remove fit from the plot
			obj.GetFunction('pol1').SetLineColor(797)

			##------------------
			##adjust the xrange
			#obj.SetLineColor(ROOT.TColor.GetColor(args.color[counter]))
			#obj.SetMarkerColor(ROOT.TColor.GetColor(args.color[counter]))
			obj.SetMarkerStyle(myMarker[counter])
			obj.SetMarkerSize(myMarkerSize[counter])
			#if (args.fill):
			#	obj.SetFillColor(ROOT.TColor.GetColor(args.color[counter]))
			
			##------------------
			##draw histograms into the same canvas (equivalent to option same)
			if (not args.legend):
				if len(filename_list) > 1:
					print "A"
					legend.AddEntry(obj, filename_list[counter-1]+'_'+graph.GetName())
				else:
					print "B"
					legend.AddEntry(obj, '_'+graph.GetName())
			else:
				print "C ", args.legend[counter-1]
				legend.AddEntry(obj, args.legend[counter-1])
				obj.SetName(args.legend[counter-1])
				obj.SetTitle(args.legend[counter-1])

			multi_graph.Add(obj,"PLA")
			

			
			##------------------
			##adjust legend and the x and y title name if chosen
			
			counter +=1
			if (args.xtitle):
				x_title = args.xtitle
			else:
				x_title = x_axis.GetTitle()
			if (args.ytitle):
				y_title = args.ytitle
			else:
				y_title = y_axis.GetTitle()
		##------------------
		##set y axis to logcolor
		if args.ylog:
			c1.SetLogy()
			ROOT.gPad.SetLogy()
		##------------------
		##draw histogram + components and save the file
		multi_graph.Draw(drawing_option)
		#legend.Draw()
		#c1.Update()
		c1.BuildLegend(0.2,0.7,0.3,0.90)
		xaxis = multi_graph.GetXaxis()
		xaxis.SetTitle(x_title)
		if 9999 not in args.xaxisrange:
			xaxis.SetRangeUser(args.xaxisrange[0], args.xaxisrange[1])
		yaxis = multi_graph.GetYaxis()
		if 9999 not in args.yaxisrange:
			yaxis.SetRangeUser(args.yaxisrange[0], args.yaxisrange[1]) 
			print 'test'
		yaxis.SetTitle(y_title)
		#if (args.ylog is True):
		if (args.output_name):
			canvasName = args.output_name
		else:
			canvasName = histogram.GetName()
		saveFile(c1, filename_list, counter, folder_loc, canvasName)
		c1.Close()
	else:
		counter = 0
		for graph in graph_list:
			##------------------
			##loop through the histograms, get all parameters and adjust the xrange
			#ROOT.gROOT.SetBatch(1)
			c1 = ROOT.TCanvas( 'test', 'Test', args.aratio[0], args.aratio[1] ) #
			obj = graph.ReadObj()
			#print 'Number of entries =', obj.GetEntries()
			x_axis = obj.GetXaxis()
			y_axis = obj.GetYaxis()
			c1.cd()
			#print obj.GetEntries()
			if 9999 not in args.xaxisrange:
				x_axis.SetRange(args.xaxisrange[0], args.xaxisrange[1])
			obj.SetLineColor(4) #Blue
			##------------------
			##set y axis to log
			if args.ylog:
				c1.SetLogy()
				ROOT.gPad.SetLogy()
			##------------------
			##set axis titles if chosen through options.
			if (args.xtitle):
				x_axis.SetTitle(args.xtitle)
			if (args.ytitle):
				y_axis.SetTitle(args.ytitle)
			obj.Draw(args.draw_option)
			#ROOT.gROOT.SetBatch(0)
			#raw_input('Press Enter to look at the next histogram')
			if (args.output_name):
				canvasName = args.output_name
			else:
				canvasName = graph.GetName()
			c1.SaveAs(canvasName+'.png')
			c1.Close()
			counter= counter+1
			#for i in xrange(obj.FindFirstBinAbove(0),obj.FindLastBinAbove(0),1):
			#'	print obj.GetBinContent(i)


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
	'-c','--channel',
	dest='channel',
	default=['9999'],
	nargs='*',
	help='used to specify the channel of the plot which should be used'
)

parser.add_argument(
	'-s','--strip',
	dest='strip',
	default=['9999'],
	nargs='*',
	help='used to specify the strip of the plot which should be used'
)
parser.add_argument(
	'-b', '--bucket',
	dest='bucket',
	default=[9999],
	nargs='*',
	type=int,
	help='used to specify the bucket of the plot which should be used | type=int'
)
parser.add_argument(
	'-k', '--kpix',
	dest='kpix',
	default=[9999],
	nargs='*',
	type=int,
	help='used to specify the bucket of the plot which should be used | type=int'
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
	'--refuse',
	dest='refuse',
	default= ['nothing'],
	nargs='*',
	help='add string that should be exluded from histogram search'
)
parser.add_argument(
	'-r', '--rebin',
	dest='rebin',
	default=1,
	type = int,
	help='add number to rebin the histograms | type=int'
)
parser.add_argument(
	'--name2',
	dest='name2',
	default=['everything'],
	nargs='*',
	help='used to specify the name of the plot which should be used'
)
parser.add_argument(
	'--refuse2',
	dest='refuse2',
	default= ['nothing'],
	nargs='*',
	help='add string that should be exluded from histogram search'
)
parser.add_argument(
	'--exact',
	dest='exact',
	default=False,
	help='if set to True, only histograms with the exact name will be used'
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
	'--order',
	dest='order',
	nargs='+',
	type=int,
	help='choose the order of plotting with same (to ensure no histograms overlap)'
)
parser.add_argument(
	'-q', '--olddaq',
	dest='olddaq',
	help='give as a command when using files from the new daq to ensure filename check etc. are correct'
)
parser.add_argument(
	'--other',
	dest='other',
	help='general plotting with non automatic file name generation'
)
parser.add_argument(
	'--nobox',
	dest='nobox',
	action='store_true',
	help='suppresses tstatbox from being added. Does not work with stack/nostack option'
)
parser.add_argument(
	'-l',
	dest='legendloc',
	nargs='+',
	type=float,
	default = [0.98, 0.99],
	help='first argument is the left x position of the legend box and second argument is the upper y position of the legend box'
)
parser.add_argument(
	'--folder',
	dest='folder',
	default='tb',
	help='tb is testbeam folder elab is elab folder. default is tb folder. there is also a thesis folder for uwes thesis plots'
)
parser.add_argument(
	'--aratio',
	dest='aratio',
	nargs='+',
	type=float,
	default=[2000,1500], help='aspect ratio of the output file'
)
parser.add_argument(
	'-f', '--fill',
	dest='fill',
	action='store_true',
	help='set whether to fill the area beneath the histogram with color'
)
parser.add_argument(
	'--upperaxis',
	dest='upperXaxis',
	action='store_true',
	help='needs to be adjusted manually all the time. adds a second x axis to the top of the plot'
)
parser.add_argument(
	'--fit',
	nargs='+',
	dest='fit',
	help='given if one wishes to add a fit to the histogram before plotting. First argument is type of function, second is left range, third is right range'
)
parser.add_argument(
	'--nofit',
	dest='nofit',
	action='store_true',
	help='uses a workaround to remove the plot from TGRAPH'
)
parser.add_argument(
	'--norm',
	dest='norm',
	action='store_true',
	help='Whether to normalize the histograms relativ to their entries.'
)
args = parser.parse_args()

# setting root style
if __name__ == '__main__':
	mystyle = ROOT.TStyle("mystyle", "My Style")
	mystyle, myMarker, myMarkerSize = myROOTStyle(args.nobox)

	mystyle.cd()
	ROOT.gROOT.ForceStyle()
	ROOT.gStyle.ls()
	# finish setting root style

	#legend_location = [0.65,0.65,0.98,0.85] # x_left, y_bottom, x_right, y_top
	legend_location = [0.15,0.65,0.35,0.85] # x_left, y_bottom, x_right, y_top

	if args.nofit:
		print ''
		print 'IMPORTANT: Please take note that NOFIT is just a very dirty hack'
		print ''
	##-----------------
	##produce empty root file and filename lists.




	root_file_list = []
	filename_list = []
	##-----------------
	##loop through all given files and add them to the list. then loop through the keys for every file..
	for root_file in args.file_in:
		root_file_list.append(ROOT.TFile(root_file))
	
		if (args.olddaq):
			filename_list.append(root_file[root_file.find('/2019_')+1:root_file.rfind('.bin')+1])
		elif (args.other):
			if 'for-lycoris' in root_file:
				nameStart = root_file.rfind('/run-') + 9
				nameEnd = root_file.rfind('root')
				print nameStart
				print nameEnd
				print 'Output name is: ', root_file[nameStart:nameEnd]
				filename_list.append(root_file[nameStart:nameEnd])
			else:
				filename_list.append(root_file[0:root_file.rfind('root')])
		else:
			if '/Run_' in root_file:
				filename_list.append(root_file[root_file.find('/Run_')+1:root_file.rfind('.dat')+1])
			elif '/Calibration_' in root_file:
				filename_list.append(root_file[root_file.find('/Calibration_')+1:root_file.rfind('.dat')+1])
	print filename_list
	#object_list = []
	for x in root_file_list:
		key_root = x.GetListOfKeys()
		loopdir(key_root)
		#object_list= object_list + (loopdir_new(key_root, args.name))
	if ('elab' in args.folder):
		folder_loc = '/scratch/plots/elab201904/'
	elif ('tb' in args.folder):
		folder_loc = '/scratch/plots/testbeam202003/'
	elif ('summer' in args.folder):
		folder_loc = '/scratch/plots/humidity/'
	elif ('thesis' in args.folder):
		folder_loc = '/scratch/plots/thesis/'
	##-----------------
	##general output
	#print args.color

	#print 'DEEEEEEBUUUUUUUG: ', object_list

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
	yRangeScale = 1.5
	if (len(hist_list) is not 0):
		hist_plotter()
	elif (len(graph_list) is not 0):
		graph_plotter()
	else:
		print 'There are NO valid histograms/graphs in the current selection'
		print ''
		print ''
	for x in root_file_list:
		ROOT.gROOT.GetListOfFiles().Remove(x)



#raw_input('Press Enter to look at the next histogram')




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


def hist_plotter():
	if ('stack' in args.draw_option or 'nostack' in args.draw_option):
		##------------------
		##initialize a canvas, a stack histogram and further variables.
		if ('nostack' not in args.draw_option):
			f = args.draw_option.find('stack')
			fend = f+5 ## End of stack
			l = len(args.draw_option)
			drawing_option = args.draw_option[0:f]+args.draw_option[fend:l] ##cutting out 'stack' from the string of drawing options
			print "Using draw options: ", drawing_option
		else:
			drawing_option = args.draw_option #exchange the same with a NOSTACK as I am using THStack
			print "Using draw options: ", drawing_option
		c1 = ROOT.TCanvas( args.output_name, 'Test', args.aratio[0], args.aratio[1] )
		c1.cd()
		#c1.SetFillColor(0)
		statBoxW = 0.15
		statBoxH = 0.105
		legend = ROOT.TLegend(args.legendloc[0], args.legendloc[1], args.legendloc[0]+statBoxW, args.legendloc[1]+statBoxH)
		hist_comp = ROOT.THStack()
		counter = 1
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
			obj.SetLineColor(args.color[counter-1])
			obj.SetMarkerColor(args.color[counter-1])
			if (args.fill):
				obj.SetFillColor(args.color[counter-1])


			##------------------
			##draw histograms into the same canvas (equivalent to option same)
			hist_comp.Add(obj, "sames")
			
			
			##------------------
			##adjust legend and the x and y title name if chosen
			if (not args.legend):
				if len(filename_list) > 1:
					legend.AddEntry(obj, filename_list[counter-1]+'_'+histogram.GetName())
				else:
					legend.AddEntry(obj, '_'+histogram.GetName())
			else:
				legend.AddEntry(obj, args.legend[counter-1])
				legendname.append(args.legend[counter-1])
			counter +=1
			if (args.xtitle):
				x_title = args.xtitle
			else:
				x_title = x_axis.GetTitle()
			if (args.ytitle):
				y_title = args.ytitle
			else:
				y_title = y_axis.GetTitle()
			
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
	elif ('same' in args.draw_option):
		##------------------
		##initialize a canvas, a stack histogram and further variables.
		drawing_option = args.draw_option
		c1 = ROOT.TCanvas( args.output_name, 'Test', args.aratio[0], args.aratio[1] )
		c1.cd()
		#c1.SetFillColor(0)
		statBoxW = 0.1
		statBoxH = 0.07
		print args.legendloc[0]
		legend = ROOT.TLegend(0.8, 0.85, 0.8+statBoxW, 0.85+statBoxH)
		counter = 1
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
			obj.SetLineColor(args.color[counter-1])
			obj.SetMarkerColor(args.color[counter-1])
			if (args.fill):
				obj.SetFillColor(args.color[counter-1])
			
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
			if 9999 not in args.yaxisrange:
				yaxis.SetRangeUser(y_low, y_high) 
				print 'test'
			yaxis.SetTitle(y_title)
#			ROOT.TGaxis.SetMaxDigits(3)
			##------------------
			##draw histograms into the same canvas (equivalent to option same)
			
			
			
			##------------------
			##adjust legend and the x and y title name if chosen
			if (not args.legend):
				if len(filename_list) > 1:
					legend.AddEntry(obj, filename_list[counter-1]+'_'+histogram.GetName())
				else:
					legend.AddEntry(obj, '_'+histogram.GetName())
			else:
				legend.AddEntry(obj, args.legend[counter-1])
				legendname.append(args.legend[counter-1])
			counter +=1
			
			obj.Draw(drawing_option)
			legend.Draw()
			
			
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
	else:
		counter = 0
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
			if 9999 in args.yaxisrange:
				if (y_low is None):
					y_low = obj.FindFirstBinAbove(0,2)-10
				elif (y_low > obj.FindFirstBinAbove(0,2)-10):
					y_low = obj.FindFirstBinAbove(0,2)-10
				if (y_high is None):
					y_high = obj.FindLastBinAbove(0,2)+10
				elif (y_high < obj.FindLastBinAbove(0,2)+10):
					y_high = obj.FindLastBinAbove(0,2)+10
				if (y_high > obj.GetNbinsY()):  #avoids overflow bin
					y_high = obj.GetNbinsY()
				if (y_low <= 0): #avoids underflow bin
					y_low = 1
				y_axis.SetRange(y_low, y_high)
			else:
				y_low = args.yaxisrange[0]
				y_high = args.yaxisrange[1]
				y_axis.SetRangeUser(y_low, y_high)
			if 9999 not in args.zaxisrange:
				z_low = args.zaxisrange[0]
				z_high = args.zaxisrange[1]
				z_axis.SetRangeUser(z_low, z_high)
			obj.SetLineColor(args.color[0]) #Blue
			obj.SetMarkerColor(args.color[0]) #Blue
			if (args.fill):
				obj.SetFillColor(args.color[0])
			
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
			obj.Draw(args.draw_option)
			run_name = filename_list[0][:-1]			
			if (args.output_name):
				outname = folder_loc+run_name+'_'+args.output_name
				print 'Creating '+outname
				#c1.SaveAs(outname+'.svg')
				c1.SaveAs(outname+'.png')
			else:
				outname = folder_loc+run_name+'_'+histogram.GetName()
				print 'Creating '+outname+'.pvg'
				#c1.SaveAs(outname+'.svg')
				c1.SaveAs(outname+'.png')
			c1.Close()
			counter= counter+1
			#for i in xrange(obj.FindFirstBinAbove(0),obj.FindLastBinAbove(0),1):
			#'	print obj.GetBinContent(i)							
	
def graph_plotter():
	if ('same' in args.draw_option):
		##------------------
		##initialize a canvas, a stack histogram and further variables.
		drawing_option = args.draw_option.replace('same', 'A') 
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
			#obj.GetFunction('pol1').SetRange(0.0, 0.0, 0.01, 0.01)
			obj.GetFunction('pol1').SetLineColor(797)
			
			##------------------
			##adjust the xrange
			obj.SetLineColor(args.color[counter-1])
			obj.SetMarkerColor(args.color[counter-1])
			if (args.fill):
				obj.SetFillColor(args.color[counter-1])
			
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
				print "C"
				legend.AddEntry(obj, args.legend[counter-1])
			
			
			multi_graph.Add(obj,"PL")
			
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
		legend.Draw()
		#c1.BuildLegend()
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
		run_name = filename_list[0][:-1]
		if (args.output_name):
			outname = folder_loc+run_name+'_'+args.output_name
			print 'Creating '+outname
			#c1.SaveAs(outname+'.svg')
			c1.SaveAs(outname+'.png')
		else:
			search_name=''
			for q in args.name:
				search_name = search_name + '_' + q
			print search_name
			outname = folder_loc+run_name+'graphs_w'+search_name
			print 'Creating '+outname+'.png'
			#c1.SaveAs(outname+'.svg')
			c1.SaveAs(outname+'.png')
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
			run_name = filename_list[0][:-1]
			if (args.output_name):
				outname = folder_loc+run_name+'_'+args.output_name
				print 'Creating '+outname
				#c1.SaveAs(outname+'.svg')
				c1.SaveAs(outname+'.png')
			else:
				search_name=''
				for q in args.name:
					search_name = search_name + '_' + q
				print search_name
				outname = folder_loc+run_name+'_'+graph.GetName()
				print 'Creating '+outname+'.pvg'
				#c1.SaveAs(outname+'.svg')
				c1.SaveAs(outname+'.png')
			c1.Close()
			counter= counter+1
			#for i in xrange(obj.FindFirstBinAbove(0),obj.FindLastBinAbove(0),1):
			#'	print obj.GetBinContent(i)


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
##mystyle.SetPadBottomMargin(0.18)
##mystyle.SetPadTopMargin(0.08)
mystyle.SetPadRightMargin(0.14)
mystyle.SetPadLeftMargin(0.14)
#
##set axis label and title text sizes
mystyle.SetLabelFont(42,"xyz")
mystyle.SetLabelSize(0.055,"xyz")
mystyle.SetLabelOffset(0.003,"yz")
mystyle.SetLabelOffset(0.00,"x")
mystyle.SetTitleFont(42,"xyz")
mystyle.SetTitleSize(0.06,"xyz")
mystyle.SetTitleOffset(1.1,"yz")
mystyle.SetTitleOffset(0.75,"x")
mystyle.SetStatFont(42)
mystyle.SetStatFontSize(0.03)


ROOT.TGaxis.SetMaxDigits(4)

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
parser.add_argument('--zrange', dest='zaxisrange', default=[9999], nargs='*', type=float, help='set a zrange for the plot to used with ymin ymax as the two arguments | type=float')
parser.add_argument('--legend', dest='legend', nargs='*', help='list of names to be used as legend titles instead of the default filename+histogram name')
parser.add_argument('--ylog', dest='ylog', help='if given as an option, set y axis to logarithmic. Remember to set the yrange to start above 0!')
parser.add_argument('--zlog', dest='zlog', help='if given as an option, set z axis to logarithmic.')
parser.add_argument('--color', dest='color', default=[ 861, 1, 418,  810, 402,  908, 435, 880,60, 632, 840, 614], nargs='*', help='list of colors to be used')
#parser.add_argument('--color', dest='color', default=[590, 591, 593, 596, 600, 602, 604, 880, 860, 632, 840, 614], nargs='*', help='list of colors to be used')
parser.add_argument('--xtitle', dest='xtitle', help='choose the name of the x axis title')
parser.add_argument('--ytitle', dest='ytitle', help='choose the name of the y axis title')
parser.add_argument('--order', dest='order', nargs='+', type=int,  help='choose the order of plotting with same (to ensure no histograms overlap)')
parser.add_argument('-q', '--olddaq', dest='olddaq', help='give as a command when using files from the new daq to ensure filename check etc. are correct')
parser.add_argument('-l', dest='legendloc', nargs='+', type=float, default = [0.98, 0.99], help='first argument is the left x position of the legend box and second argument is the upper y position of the legend box')
parser.add_argument('--folder', dest='folder', default='tb', help='tb is testbeam folder elab is elab folder. default is elab folder.')
parser.add_argument('--aratio', dest='aratio', nargs='+', type=float,  default=[1200,900], help='aspect ratio of the output file')
parser.add_argument('-f', '--fill', dest='fill', action='store_true', help='set whether to fill the area beneath the histogram with color')
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
	
	if (args.olddaq):
		filename_list.append(root_file[root_file.find('/2019_')+1:root_file.rfind('.bin')+1])
	else:
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
	folder_loc = '/home/lycoris-dev/Documents/elab201904/'
elif ('tb' in args.folder):
	folder_loc = '/home/lycoris-dev/Documents/testbeam201907/'
elif ('summer' in args.folder):
	folder_loc = '/home/lycoris-dev/Documents/humidity/'
##-----------------	
##general output
#print args.color


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




#!usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import string
import ROOT
import argparse
import argcomplete
from operator import add
import sys


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


class MakeStrip:
    def __init__(self, file_in):
        if not file_in: 
            exit("Your file_in is empty!")
              

        self.rootfile = ROOT.TFile(file_in)
       
        self.file_base_name = args.file_in[args.file_in.find('/20')+1:args.file_in.rfind('.bin')]

        self.left_strip = None
        self.right_strip = None

        self.kpix = []
        
        return

    def __del__(self):
        return
    
    def __loophist(self, hist_obj):
        if not self.kpix or len(self.kpix)!=2 :
            print "__loophist: your kpix input is empty!"
            return
        
	if ('Left_Strip_entries_k_'+str(self.kpix[0])) in hist_obj.GetName() and '_total' in hist_obj.GetName():
	    self.left_strip = hist_obj
            print self.left_strip.GetName()
            
        if ('Right_Strip_entries_k_'+str(self.kpix[1])) in hist_obj.GetName() and '_total' in hist_obj.GetName():
	    self.right_strip = hist_obj
            print self.right_strip.GetName()
            
        return
    
    def __loopdir(self, keys):
        for key_object in keys:
            if ('TDirectory' in key_object.GetClassName()):
	        self.__loopdir(key_object.ReadObj().GetListOfKeys())
	    else:
		#print key_object.GetName()
                self.__loophist(key_object)
	return

    def run(self, kpix):
        if len(kpix)!=2 :
            exit("StripMaker: your kpix input wrong!")

        self.kpix = kpix
        print self.kpix, len(self.kpix)
        
        key_root = self.rootfile.GetListOfKeys()
        
        self.__loopdir(key_root)

        if self.left_strip:
            left_strip_hist = self.left_strip.ReadObj()
        if self.right_strip:
	    right_strip_hist = self.right_strip.ReadObj()


        Input_left = []
        Input_right = []
        C = []

        
        for chan in xrange(920):
            if self.left_strip:
                Input_left.append(left_strip_hist.GetBinContent(chan+1))
	    if self.right_strip:
                Input_right.append(right_strip_hist.GetBinContent(chan+1))

        if Input_left:
            C = [Input_left]
            fig = plt.figure()
            ax = fig.add_subplot(111)
            axes = plt.axes() 
            axes.get_yaxis().set_visible(False)
            axes.set_xlim(0,920)
            pcm = ax.pcolormesh(C, cmap='viridis')
            plt.colorbar(pcm)
            filename = '/home/lycoris-dev/Documents/'+self.file_base_name+'_strip_map_left.png'
            print "File is saved in " +  filename
            plt.savefig(filename, dpi = 300)
            
        if Input_right:
            C = [Input_right]
            fig = plt.figure()
            ax = fig.add_subplot(111)
            axes = plt.axes() 
            axes.get_yaxis().set_visible(False)
            axes.set_xlim(0,920)
            pcm = ax.pcolormesh(C, cmap='viridis')
            plt.colorbar(pcm)
            filename = '/home/lycoris-dev/Documents/'+self.file_base_name+'_strip_map_right.png'
            print "File is saved in " +  filename
            plt.savefig(filename, dpi = 300)

        plt.close()

        ROOT.gROOT.GetListOfFiles().Remove(self.rootfile)
        self.rootfile.Close()


        return
        
######### Main - BEG #########
   
parser = MyParser()
parser.add_argument('-f', dest='file_in', help='input root file')
parser.add_argument('-k1', dest='lkpix', nargs='+', default=[], type=int, help='which left kpix to plot, eg: 26 or 26 20')
parser.add_argument('-k2', dest='rkpix', nargs='+', default=[], type=int, help='which right kpix to plot, eg: 28 or 28 22')

if len(sys.argv) < 2 :
	print parser.print_help()
	sys.exit(1)

args = parser.parse_args()
stripMaker = MakeStrip(args.file_in)
        
for k1 in args.lkpix:
    print " left kpix:  ", k1
for k2 in args.rkpix:
    print " right kpix: ", k2

#if len(args.lkpix)!=len(args.rkpix):

    
stripMaker.run([26,0])

sys.exit(1)
######### Main - END #########
##-------------------------------------##

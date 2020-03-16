#!/usr/bin/env python
'''
2018 Summer
 - Originally by Uwe Kraemer <uwe.kraemer@desy.de>
 - Update     by Mengqing <mengqing.wu@desy.de>
 - Target:    For making kpix channel disable mapping, 
          based on calibration /and pedestal files.
 - Copyright reserved for DESY
 - Latest Modification: 
   * 16/03/2020, MW: adapt to new calibration histogram nomeclature
'''

import numpy as np
import os 
import ROOT
import argparse
import argcomplete
import sys


class MyParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

class ChannelDisabler:
	def __init__(self, file_in=[], kpix='0', ecal=False):

		'''
		Achtung! the file_in can only contain 1 of each type of root file, 
		i.e. it can have only 1 x calibration.root, 1 x pedestal.root, like this.
		'''     
		if not file_in:
			exit("Your file_in is empty")
		self.ecal = ecal  # Default = False
		self.slope_name =  'slope_vs_channel_k'
		self.RMS_name = 'RMSfc_vs_channel_k'
		self.acquire_name = 'Channel_entries_k'
	
		self.kpix = kpix
		if not kpix:
			print "Warning: no kpix chosen, BE SURE your input root file has 1 kpix running!"
		else:
			self.slope_name += ''+str(kpix)+'_b0'
			self.RMS_name += ''+str(kpix)+'_b0'
			self.acquire_name += str(kpix)

		print self.slope_name
		print self.RMS_name 
		print self.acquire_name 
			
		self.root_file_list=[]
		for root_file in file_in:
			self.root_file_list.append(ROOT.TFile(root_file)) # open the file and add address to an array
		
		self.RMS_hist = None #TKey
		self.slope_hist = None #TKey
		self.acquire_hist = None #TKey
		
		mapping_file_name = './include/kpix_left_and_right.h'
		self.mapping_file = open(mapping_file_name, "r")
	
		if not kpix:
			self.filename_map = './data/disable.txt' # output file
		else:
			self.filename_map = './data/disabler_k'+str(kpix)+'.txt' 

	def ChangeKpix(self, newkpix):
		if not newkpix:
			return

		self.slope_name.replace( str(self.kpix), str(newkpix) )
		self.RMS_name.replace( str(self.kpix), str(newkpix) ) 
		self.acquire_name.replace( str(self.kpix), str(newkpix) )
		self.filename_map = './data/disabler_k'+str(newkpix)+'.txt' 

		self.kpix = newkpix
		
	def __del__(self):
		self.mapping_file.close()
		
	def __loophist(self, hist_obj):
		#print "test3",  self.acquire_name
		if self.slope_name == hist_obj.GetName():
			self.slope_hist = hist_obj
		elif self.RMS_name == hist_obj.GetName():
			self.RMS_hist = hist_obj
		elif self.acquire_name in hist_obj.GetName() and '_total' in hist_obj.GetName() and 'timed' not in hist_obj.GetName() and 'no' not in hist_obj.GetName():# 
			self.acquire_hist = hist_obj
		return
	
	def __loopdir(self, keys):
		# loop through all subdirectories of the root file and look for the vs channel histograms
		## declare that the beforementioned variables are global variables
		#global RMS_hist, slope_hist, acquire_hist

		for key_object in keys:
			if ('TDirectory' in key_object.GetClassName()):
				self.__loopdir(key_object.ReadObj().GetListOfKeys())
			else:
			#print key_object.GetName()
				self.__loophist(key_object)
		return        

	def run(self):
		
		for x in self.root_file_list:
			key_root = x.GetListOfKeys()
			self.__loopdir(key_root)

		noisy_channels = []
		noisy_channels_acquire = []
		dead_channels_slope = []
		dead_channels_RMS = []
		dc_channels = []
		if self.slope_hist:
			print self.slope_hist
			slope_obj = self.slope_hist.ReadObj()
		if self.RMS_hist:
			print self.RMS_hist
			RMS_obj = self.RMS_hist.ReadObj()
		if self.acquire_hist:
			print self.acquire_hist
			acquire_obj = self.acquire_hist.ReadObj()
				
		for chan in xrange(1024): #level1: loop over histos
			if self.slope_hist: #level2
				if (slope_obj.GetBinContent(chan+1) <= 1.0):
					#GetBinContent(0) is the underflow bin, here the bin counting starts at 1 therefore chan+1
					dead_channels_slope.append(chan)
			if self.RMS_hist: #level2
				if (RMS_obj.GetBinContent(chan+1) <= 0.05):
					dead_channels_RMS.append(chan)
				if (RMS_obj.GetBinContent(chan+1) >= 4):
					noisy_channels.append(chan)
			if self.acquire_hist: #level2
				if (acquire_obj.GetBinContent(chan+1) >= 0.4):
					noisy_channels_acquire.append(chan)

		if (not self.ecal):
			for line in self.mapping_file: #level1: check channel map for disconnected ones
				if 'm1.insert' in line:
					line=line[22:-4]
					line_split=line.split(',')
					if (line_split[1] == '9999') and (int(line_split[0]) not in dc_channels):
						dc_channels.append(int(line_split[0]))	

				
		kpix = np.chararray(1024)		#create mapper kpix
		#filename_map = './data/disable.txt'	#choose filename
	
		print 'Dead channels from calibration slope = ', dead_channels_slope
		print 'Dead channels from calibration RMS = ', dead_channels_RMS
		print 'Very high noise channels = ', noisy_channels
		print 'Channels not connected to strips based on mapping = ', dc_channels
		print 'High noise channels in acquire running = ', noisy_channels_acquire
		print ''


		for chan in xrange(1024):		#create the character mapping of the channels
			
			if (chan in dc_channels) or (chan in dead_channels_slope) or (chan in dead_channels_RMS) or (chan in noisy_channels) or (chan in noisy_channels_acquire):		#A for active and D for deactivated channels
				kpix[chan] = 'D'
			else :
				kpix[chan] = 'A'

		for i in xrange(1024):
			if ((i%32 == 0) and (i is not 0)):
				print ''
				sys.stdout.write(kpix[i])
			else:
				sys.stdout.write(kpix[i])
			if (i == 7):
				sys.stdout.write(' ')
			elif (i > 7 and (i-7)%8 == 0 and i is not 0 and i%32 is not 0):
				sys.stdout.write(' ')
					
		disChanCount=0
		for chan in kpix:
			if chan == 'D':
				disChanCount+=1

		kpix = np.reshape(kpix, (32, 32))	#reshape it for fitting in xml file

		print '\n\n'
		np.savetxt(self.filename_map, kpix, fmt= ['%.1c' ,'%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%2.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%2.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%2.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c'], delimiter = '')
		print 'Saving deactivation mapping in: ', self.filename_map
		
		##save it in filename with a format that fits in the xml file, it adds the space after 8 characters so it's easier to see (but a pain to write it like this, must be another way) 
		
		#mapping_file.close() #close tracker_to_kpix
		
		
		for x in self.root_file_list:
			ROOT.gROOT.GetListOfFiles().Remove(x)
			x.Close()

		print '\n How many channels disabled? ', disChanCount, '\n'
		return

######### Parser #########
parser = MyParser()
#parser = argparse.ArgumentParser()
#parser.add_argument('file_in', nargs='+', help='name of the input files')
parser.add_argument('--calib', dest='calib', nargs='?', help='input calibration file')
parser.add_argument('--acquire', dest='acquire', nargs='?', help='input pedestal file')
parser.add_argument('-k', dest='kpix', nargs='+', default=[], type=str, help='which kpix to check, eg: 1 or 1 2')
parser.add_argument('--ecal', dest='ecal', type=bool, nargs='?', default=False, help='ecal sensor used here')

if len(sys.argv) < 2:
	print parser.print_help()
	sys.exit(1)

args = parser.parse_args()
######### Parser #########

######### Main: #########

file_in=[]
if args.calib:
	file_in.append(args.calib)
if args.acquire:
	file_in.append(args.acquire)

print args.kpix
if args.kpix:
	for k in args.kpix:
		disabler = ChannelDisabler(file_in, k, args.ecal)
		disabler.run()
else:
	disabler = ChannelDisabler(file_in, 0, args.ecal)
	disabler.run()
		





#!usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import string
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.ticker as mticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import ROOT
from ROOT import TFile, gDirectory, TDirectory
import argparse
import argcomplete
from operator import add

hist_list= [] #Global list of all chosen histograms (only works because its a list, with single variables they need to be declared as global in the function

def numberlist(count,num_of_x):
    listofnums = [count] * int(num_of_x)
    return listofnums

def ECAL_map(column, row, num_of_x, ymin, ymax):
	for count in xrange(39):
		if ( count == 0 ):
			column.extend(np.arange(ymax, ymin, -1))
			num_of_x = ymax-ymin
			row.extend(numberlist(count*np.cos(np.pi/3), num_of_x))
		if (( count > 0 ) and (count <= 7)):
			ymax = ymax + 1.5
			ymin = ymin - 1.5
			column.extend(np.arange(ymax, ymin, -1))
			num_of_x = ymax-ymin
			row.extend(numberlist(count*np.cos(np.pi/3), num_of_x))
		if ((count >= 8) and (count <= 11)):
			ymax = ymax + 0.5
			ymin = ymin - 0.5
			column.extend(np.arange(ymax, ymin, -1))
			num_of_x = ymax-ymin
			row.extend(numberlist(count*np.cos(np.pi/3), num_of_x))
		if ((count >= 12) and (count < 28)):
			if (count%2 == 1):
				ymax = ymax - 0.5
				ymin = ymin + 0.5
				num_of_x = ymax-ymin
			else:
				ymax = ymax + 0.5
				ymin = ymin - 0.5
				num_of_x = ymax-ymin
			if ((count >= 15) and (count <= 23)):
				if ((count == 15) or (count == 23)):
					column.extend(np.arange(ymax,ymin+num_of_x/2+0.5, -1))
					column.extend(np.arange(ymax-num_of_x/2+0.5, ymin+num_of_x/2-0.5, -0.5))
					column.extend(np.arange(ymax-num_of_x/2-0.5,ymin, -1))
					num_of_x = ymax-ymin + 1
					#print num_of_x
					row.extend(numberlist(count*np.cos(np.pi/3), num_of_x))
				elif ((count == 16) or (count == 22)):
					column.extend(np.arange(ymax,ymin+num_of_x/2+3, -1))
					column.extend(np.arange(ymax-num_of_x/2+3, ymin+num_of_x/2-3, -0.5))
					column.extend(np.arange(ymax-num_of_x/2-3,ymin, -1))
					num_of_x = ymax-ymin + 6
					#print num_of_x
					row.extend(numberlist(count*np.cos(np.pi/3), num_of_x))
				elif ((count == 17) or (count == 19) or (count == 21)):
					column.extend(np.arange(ymax,ymin+num_of_x/2+3.5, -1))
					column.extend(np.arange(ymax-num_of_x/2+3.5, ymin+num_of_x/2-3.5, -0.5))
					column.extend(np.arange(ymax-num_of_x/2-3.5,ymin, -1))
					num_of_x = ymax-ymin + 7
					#print num_of_x
					row.extend(numberlist(count*np.cos(np.pi/3), num_of_x))
				elif ((count == 18) or (count == 20)):
					column.extend(np.arange(ymax,ymin+num_of_x/2+4, -1))
					column.extend(np.arange(ymax-num_of_x/2+4, ymin+num_of_x/2-4, -0.5))
					column.extend(np.arange(ymax-num_of_x/2-4,ymin, -1))
					num_of_x = ymax-ymin + 8
					#print num_of_x
					row.extend(numberlist(count*np.cos(np.pi/3), num_of_x))
			else:
				column.extend(np.arange(ymax, ymin, -1))
				row.extend(numberlist(count*np.cos(np.pi/3), num_of_x))
		if ((count >= 28) and (count <= 31)):
			ymax = ymax - 0.5
			ymin = ymin + 0.5
			num_of_x = ymax-ymin
			column.extend(np.arange(ymax, ymin, -1))
			row.extend(numberlist(count*np.cos(np.pi/3), num_of_x))
		if (count >= 32):
			ymax = ymax - 1.5
			ymin = ymin + 1.5
			num_of_x = ymax-ymin
			column.extend(np.arange(ymax, ymin, -1))
			row.extend(numberlist(count*np.cos(np.pi/3), num_of_x))
	

def map_text(ax, channel_kpix, mapping, row, column):
	pixel = 0
	for i in channel_kpix:
		if (mapping_sensor_to_kpix[pixel] < 10):
			ax.text(row[i]-0.10,column[i], mapping[pixel], fontsize=3)
		elif (mapping_sensor_to_kpix[pixel] < 100):
			ax.text(row[i]-0.13,column[i], mapping[pixel], fontsize=3)
		elif (mapping_sensor_to_kpix[pixel] < 1000):
			ax.text(row[i]-0.18,column[i], mapping[pixel], fontsize=3)
		elif (mapping_sensor_to_kpix[pixel] >= 1000):
			ax.text(row[i]-0.24,column[i], mapping[pixel], fontsize=3)
		pixel = pixel +1
	return ax


def map_plot_ecal(channel_kpix, mapping, sensor_row, sensor_column, color, list_of_channel_values, bucket_name, bucket_nr, color_values, type_name, color_type):
		fig = plt.figure()
		ax = fig.add_subplot(111) 
		map_text(ax, channel_kpix, mapping, row, column)
		plt.xlim(-1, 24)
		plt.ylim(-24, 16)
		plt.axis('off')
		if ('log' in color_type):
			plsctr = plt.scatter(sensor_row, sensor_column, c=list_of_channel_values[bucket_nr], s = 100, marker='H', cmap='viridis', norm=mpl.colors.SymLogNorm(linthresh = color_values[1], vmin=color_values[0], vmax=color_values[-1]))
		elif ('log' not in color_type):
			plsctr = plt.scatter(sensor_row, sensor_column,  c=list_of_channel_values[bucket_nr], s = 100, marker='H', cmap='viridis', vmin=color_values[0], vmax=color_values[-1])
		if args.file_out:
			if ('mean' in type_name):
				Filename=args.file_out  + '_mean_' + args.kpix + '_' +  bucket_name + '_ecal_map.png'
			elif ('rms' in type_name):
				Filename = args.file_out +'_RMS_' +  args.kpix + '_' +  bucket_name + '_ecal_map.png'
			elif ('slope' in type_name):
				Filename = args.file_out +'_slope_' +  args.kpix + '_' +  bucket_name +'_ecal_map.png'
			elif ('entry' in type_name):
				Filename = args.file_out +'_entries_' +  args.kpix + '_' +  bucket_name + '_ecal_map.png'
		elif ('mean' in type_name):
			Filename = args.file_in +'_mean_' +  args.kpix + '_' +  bucket_name + '_ecal_map.png'
		elif ('rms' in type_name):
			Filename = args.file_in +'_RMS_' +  args.kpix + '_' +  bucket_name + '_ecal_map.png'
		elif ('slope' in type_name):
			Filename = args.file_in +'_slope_' +  args.kpix + '_' +  bucket_name + '_ecal_map.png'
		elif ('entry' in type_name):
			Filename = args.file_in +'_entries_' +  args.kpix + '_' +  bucket_name + '_ecal_map.png'
		plt.colorbar(plsctr)
		print "File is saved in " +  Filename
		plt.savefig(Filename, dpi = 300)
		plt.close()	
	
def log_tick_formatter(val, pos=None):
    return "{:.2e}".format(10**val)
    
    
def loopdir(keys):  # loop through all subdirectories of the root file and add all histograms that fulfill the below given criteria to a list
	for key_object in keys:
		if ('TDirectory' in key_object.GetClassName()):
			loopdir(key_object.ReadObj().GetListOfKeys())
		else:
			if (('hist_s' in key_object.GetName()) and ('_b' in key_object.GetName()) and ('_k'+args.kpix in key_object.GetName()) and ('time' not in key_object.GetName())):
				hist_list.append(key_object)


parser = argparse.ArgumentParser()
parser.add_argument('file_in', help='name of the input file')                
parser.add_argument('-k', '--kpix', dest='kpix', type = str, help='Choose the KPiX connection slot on the DAQ board')
#parser.add_argument('-f', '--file_in', dest='file_in', help="Input filename")
parser.add_argument('-o', '--file_out', dest='file_out', help="Output filename")
parser.add_argument('--calib', dest='calib', help='True if calibration file')
parser.add_argument('-t', '--timed', dest='timed', help='Choose whether or not to take time stamped plots')
args = parser.parse_args()


channel_kpix = range(1024)
column = []
row = []

ymax = 0
ymin = -6

#file_map_kpix_to_sensor = '/afs/desy.de/user/k/kraemeru/public/gui_sort.txt'
#file_map_sensor_to_kpix = '/afs/desy.de/user/k/kraemeru/public/gui_hama.txt'
file_map_kpix_to_sensor_ecal = 'python/gui_sort.txt'
file_map_sensor_to_kpix_ecal = 'python/gui_hama.txt'

bad = []
noise = []
bad_fit = []

mapping_kpix_to_sensor = []
mapping_sensor_to_kpix = []




adcmean_kpix = np.zeros((4,1024))
adcmean_sensor = np.zeros((4,1024))

adcRMS_kpix = np.zeros((4,1024))
adcRMS_sensor = np.zeros((4,1024))

adchannel_kpix = np.zeros((4,1024))

channel_slope_kpix = np.zeros((4,1024))
channel_slope_sensor = np.zeros((4,1024))

num_of_entries_kpix = np.zeros((4,1024))
num_of_entries_sensor = np.zeros((4,1024))


if ((not args.file_in) or ('root' not in args.file_in)):
	print "Please specify a root file as a first argument"
	quit()

	
root_file = ROOT.TFile(args.file_in)
bucket = ['b0', 'b1', 'b2', 'b3']


key_root = root_file.GetListOfKeys()
loopdir(key_root)

# Read the file and get the mean histogram value and the corresponding channel number
total = 0
maximum_value = 0
for histogram in hist_list:
	obj = histogram.ReadObj()
	num_of_entries_kpix
	str1_start = obj.GetName().find('_c') + 2
	str1_end = obj.GetName().find('_c') +6
	str2_start = obj.GetName().find('_b')+2
	str2_end = obj.GetName().find('_b')+3
	channel_nr = int(obj.GetName()[str1_start:str1_end])
	bucket_nr = int(obj.GetName()[str2_start:str2_end])
	entries = 0
	for c in xrange(8192):
		entries = entries + obj.GetBinContent(c)
	num_of_entries_kpix[bucket_nr][channel_nr] = entries
	if (entries > maximum_value):
		maximum_value = entries
	total = total + obj.GetEntries()
print "Registered events = ", total

#print num_of_entries_kpix

#mapping
with open(file_map_kpix_to_sensor_ecal) as mapfile:
	for line in mapfile:
		linesplit = str(line).split( )
		mapping_kpix_to_sensor.append(int(linesplit[0])-1)
with open(file_map_sensor_to_kpix_ecal) as mapfile:
	for line in mapfile:
		linesplit = str(line).split( )
		mapping_sensor_to_kpix.append(int(linesplit[1]))
	
			


for i in channel_kpix:
	for b in xrange(len(bucket)):
		adcmean_sensor[b][i] = adcmean_kpix[b][mapping_sensor_to_kpix[i]]
		adcRMS_sensor[b][i] = adcRMS_kpix[b][mapping_sensor_to_kpix[i]]
		num_of_entries_sensor[b][i] = num_of_entries_kpix[b][mapping_sensor_to_kpix[i]]
		if (args.calib == True):
			channel_slope_sensor[b][i] = channel_slope_kpix[b][mapping_sensor_to_kpix[i]]
	

#print num_of_entries_kpix[0][499]
#print num_of_entries_sensor[0][0]
pixel_count = []
x = 0
num_of_x = 0
## PRODUCING AN EXTERNAL FUNCTION FOR MAPPING and subsequent plotting
ECAL_map(column, row, num_of_x, ymin, ymax)  ## PRoducing the scatter map for the ECAL sensor, row and column coordinates
color_entry = [0.0, maximum_value] 
bucket_nr = 0
for b in bucket: #entries map
	map_plot_ecal(channel_kpix, mapping_sensor_to_kpix, row, column, color_entry, num_of_entries_sensor, b, bucket_nr, color_entry, 'entry', 'norm')
	bucket_nr = bucket_nr + 1


	

#!usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import string
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.ticker as mticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import argparse
import re
import itertools
import pandas as pd
from operator import add

hist_list= [] #Global list of all chosen histograms (only works because its a list, with single variables they need to be declared as global in the function

def numberlist(count,num_of_x):
    listofnums = [count] * int(num_of_x)
    return listofnums

def ECAL_map(column, row, num_of_x, ymin, ymax):
	for count in range(39):
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
	

def map_text(ax, channel_kpix, mapping, row, column): # writes the kpix channel numbers into the plot
	pixel = 0
	for i in channel_kpix:
		if (mapping_sensor_to_kpix[pixel] < 10): #make sure the text is centered that is why we shift it further left the larger the text is.
			ax.text(row[i]-0.08,column[i], mapping[pixel], fontsize=3, color='white')
		elif (mapping_sensor_to_kpix[pixel] < 100):
			ax.text(row[i]-0.16,column[i], mapping[pixel], fontsize=3, color='white')
		elif (mapping_sensor_to_kpix[pixel] < 1000):
			ax.text(row[i]-0.22,column[i], mapping[pixel], fontsize=3, color='white')
		elif (mapping_sensor_to_kpix[pixel] >= 1000):
			ax.text(row[i]-0.30,column[i], mapping[pixel], fontsize=3, color='white')
		pixel = pixel +1
	return ax


def map_plot_ecal(channel_kpix, mapping, sensor_row, sensor_column, color, list_of_channel_values, color_values, type_name, color_type): #produces the colored plot
		fig = plt.figure() #generate a new figure
		ax = fig.add_subplot(111)  #generate a new plot

		map_text(ax, channel_kpix, mapping, row, column)
		plt.xlim(-1, 20) #set x axis limit
		plt.ylim(-20, 15) #set y axis limit
		if args.no_axis:
			plt.axis('off') # do not show x and y axis
		if ('log' in color_type):
			plsctr = plt.scatter(sensor_row, sensor_column, c=list_of_channel_values, s = 70, marker='H', cmap='viridis', norm=mpl.colors.SymLogNorm(linthresh = color_values[1], vmin=color_values[0], vmax=color_values[-1]))
		elif ('log' not in color_type):
			plsctr = plt.scatter(sensor_row, sensor_column,  c=list_of_channel_values, s = 70, marker='H', cmap='viridis', vmin=color_values[0], vmax=color_values[-1])
		if args.file_out:
			Filename=args.file_out + '_ecal_map.png'
		elif not args.file_in:
			Filename = 'ecal_map.png'
		else:
			Filename = args.file_in + '_ecal_map.png'
		cbar = plt.colorbar(plsctr)
		cbar.set_label(args.zaxistitle, rotation=270)
		print("File is saved in " +  Filename)
		plt.savefig(Filename, dpi = 300)
		plt.close()	
	
def log_tick_formatter(val, pos=None):
    return "{:.2e}".format(10**val)
    


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', dest='file_in', help='name of the input file')
parser.add_argument('-o', '--file_out', dest='file_out', help="Output filename")
parser.add_argument('-m', '--max', dest='maximum_value', default = 20, help="maximum value on the colored z axis")
parser.add_argument('-z', '--z_axis_title', dest='zaxistitle', default = "Z", help="title of the z axis")
parser.add_argument('-n', '--no_axis', dest='no_axis', action='store_true', help="called when x and y axis should not be drawn")

args = parser.parse_args()


if __name__ == '__main__':
    channel_kpix = range(1024)
    column = []
    row = []

    ymax = 0
    ymin = -6

    #file_map_kpix_to_sensor = '/afs/desy.de/user/k/kraemeru/public/gui_sort.txt'
    #file_map_sensor_to_kpix = '/afs/desy.de/user/k/kraemeru/public/gui_hama.txt'
    file_map_kpix_to_sensor_ecal = 'python/gui_sort.txt' #specify location of the mapping gui_sort adjust to your file path
    file_map_sensor_to_kpix_ecal = 'python/gui_hama.txt' #specify location of the mapping gui_hama adjust to your file path


    mapping_kpix_to_sensor = []
    mapping_sensor_to_kpix = []

    #mapping
    with open(file_map_kpix_to_sensor_ecal) as mapfile:
        for line in mapfile:
            linesplit = str(line).split( )
            mapping_kpix_to_sensor.append(int(linesplit[0])-1)
    with open(file_map_sensor_to_kpix_ecal) as mapfile:
        for line in mapfile:
            linesplit = str(line).split( )
            mapping_sensor_to_kpix.append(int(linesplit[1]))


    entries_kpix = np.zeros(1024)
    entries_sensor = np.zeros(1024)


    if not args.file_in: #if no input file was given
        print("Please specify a text file as argument using random normal distributions as test")
        entries_kpix = np.random.normal(10,6,1024)
    else: #if an input file was given
        file = pd.read_csv(args.file_in, header=None, nrows=1024)  # some of the files have a weird entry in row 4096 that I do not care about that is why I tell him to only read 2014 rows
        file_entries = file.to_numpy()
        for i in file_entries:
            entries_kpix[int(i[0])] = i[1] # I assign to the channel specified by the first column the value in the second column. This is done so in case one of your files does not have channel 0 in row 0 it does not screw up.



    for i in range(len(entries_kpix)):
        entries_sensor[i] = entries_kpix[mapping_sensor_to_kpix[i]]
    x = 0
    num_of_x = 0
    ECAL_map(column, row, num_of_x, ymin, ymax)  ## PRoducing the scatter map for the ECAL sensor, row and column coordinates
    color_entry = [0.0, args.maximum_value]
    map_plot_ecal(channel_kpix, mapping_sensor_to_kpix, row, column, color_entry, entries_sensor, color_entry, 'entry', 'norm')



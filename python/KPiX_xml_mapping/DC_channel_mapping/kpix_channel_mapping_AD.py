import numpy as np
import os 

filename_tracker_to_kpix = "tracker_to_kpix_left.txt"
f = open(filename_tracker_to_kpix, "r")#open the tracker_to_kpix (left/right) file 

#filename_noisy_channels = "channel_list_noise.txt"
#c = open(filename_noisy_channels, "r") #open noisy channel file

channels = {}
noisy_channels = []

for line in f:				#convert the file and numbers into a dictionairy
	test = line.split() 	#grab single line as list
	test = map(int, test) 	#convert string of numbers to int
	if len(test) == 0:
		continue
	key = test[0]			#first number kpix_channel
	val = test[1]			#second number strip_number
	channels[int(key)] = val	
	
"""
for line in c:
	line = line.split()
	line = map(int,line)
	if len(line) == 0:
		continue
	noisy_channels.append(line[0])
	
print(noisy_channels)
"""	

dc_channels = []


for key,value in channels.items():	#store disconnected channels with strip value of 9999 in list dc_channels
	if value == 9999:
		dc_channels.append(key)
		
kpix = np.chararray(1024)		#create mapper kpix
filename_map = 'map_DC_first128_xml.txt'		#choose filename
	
for char in range(1024):		#create the character mapping of the channels
	if char in dc_channels :		#A for active and D for deactivated channels 
		kpix[char] = 'D'
	else :
		kpix[char] = 'A'
print len(dc_channels)

for char in range(128):
	kpix[char] = 'D'
	
kpix = np.reshape(kpix, (32, 32))	#reshape it for fitting in xml file

np.savetxt(filename_map, kpix, fmt= ['%.1c' ,'%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%2.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%2.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%2.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c', '%.1c' , '%.1c'], delimiter = '')

#save it in filename with a format that fits in the xml file, it adds the space after 8 characters so it's easier to see (but a pain to write it like this, must be another way) 

f.close() #close tracker_to_kpix






		





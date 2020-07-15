#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import argcomplete
import sys
import re

parser = argparse.ArgumentParser()
parser.add_argument(
	'file_in',
	nargs='+',
	help='name of the input file'
)
args = parser.parse_args()
if len(sys.argv) < 2:
	print parser.print_help()
	sys.exit(1)
print ''

def fieldassign(fields):
	layer = int(fields[0])
	x = float(fields[1])
	y = float(fields[2])
	z = float(fields[3])
	SoN = float(fields[4])
	cSize = int(fields[5])
	charge = float(fields[6])
	ID = int(fields[7])
	return layer,x,y,z,SoN,cSize,charge,ID



if __name__ == '__main__':
	if len(args.file_in) != 2:
		print "Please enter two file arguments for files to be compared"
		sys.exit(1)
	file1 = open(args.file_in[0])
	file2 = open(args.file_in[1])

	outFile = open("file_compare.txt", 'w')

	eventNumber = 0
	maxEvent = 0
	eventContent1 = {}
	for line in file1:
		fields =  re.split('=| |\n', line) ## splitting fields using regular expression
		fields = filter(None, fields) ## getting rid of empty fields
		if len(fields) == 3:
			eventNumber = int(fields[1])
			if (eventNumber >= maxEvent):
				maxEvent = eventNumber
			eventContent1[eventNumber] = []
			continue
		eventContent1[eventNumber].append(line)

	eventNumber = 0
	eventContent2 = {}
	for line in file2:
		fields =  re.split('=| |\n', line) ## splitting fields using regular expression
		fields = filter(None, fields) ## getting rid of empty fields
		if len(fields) == 3:
			eventNumber = int(fields[1])
			if (eventNumber >= maxEvent):
				maxEvent = eventNumber
			eventContent2[eventNumber] = []
			continue
		eventContent2[eventNumber].append(line)
	print maxEvent

	file1No = 0
	file2No = 0
	fileBothNo = 0
	IDadded = []
	for i in range(maxEvent):
		if (i not in eventContent1) and (i not in eventContent2):
			#print "Event ", i, " was not found in either sample"
			fileBothNo += 1
		if (i not in eventContent1) and (i in eventContent2):
			#print "Event ", i, " was not found in sample 1"
			file1No += 1
		if (i in eventContent1) and (i not in eventContent2):
			#print "Event ", i, " was not found in sample 2"
			file2No += 1
		if (i in eventContent1) and (i in eventContent2):
			outFile.write('0 %i 0.0 \n' %i)
			for line2 in eventContent2[i]:
				fields2 =  re.split('=| |\n', line2) ## splitting fields using regular expression
				fields2 = filter(None, fields2) ## getting rid of empty fields
				layer2,x2,y2,z2,SoN2,cSize2,charge2,ID2 = fieldassign(fields2)
				for line1 in eventContent1[i]:
					fields1 =  re.split('=| |\n', line1) ## splitting fields using regular expression
					fields1 = filter(None, fields1) ## getting rid of empty fields
					layer1,x1,y1,z1,SoN1,cSize1,charge1,ID1 = fieldassign(fields1)
					if (layer1 == layer2) and (abs(y1-y2) < 1.0) and (ID2 not in IDadded):
						outFile.write('%i %f %f %f %f %i %f %i \n' %(layer2,x2,y2,z2,SoN2,cSize2,charge2,ID2))
						IDadded.append(ID2)



	print "Did not find a corresponding event in file 1 ", file1No, " times"
	print "Did not find a corresponding event in file 2 ", file2No, " times"
	print "Did not find a corresponding event in both files ", fileBothNo, " times"

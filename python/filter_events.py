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

if __name__ == '__main__':
	if len(args.file_in) != 2:
		print "Please enter two file arguments for files to be compared"
		sys.exit(1)
	file1 = open(args.file_in[0])
	file2 = open(args.file_in[1])
	outName = args.file_in[1][:-4]+'_filter.txt'
	outFile = open(outName, 'w')

	eventsDieter = []
	eventsUwe = []

	for line in file1:
		fields =  re.split('=| |\n|\r\n', line) ## splitting fields using regular expression
		fields = filter(None, fields) ## getting rid of empty fields
		eventsDieter.append(int(fields[0]))
	writeFlag = True
	for line in file2:
		fields =  re.split('=| |\n|\r\n', line) ## splitting fields using regular expression
		fields = filter(None, fields) ## getting rid of empty fields
		#print len(fields)
		if (len(fields) == 4):
			eventNumber = int(fields[2])
			print "Currently at event: ", eventNumber
			if eventNumber in eventsDieter:
				writeFlag = True
			else:
				writeFlag = False
		if writeFlag == True:
			outFile.write(line)


	print 'Writing file to ', outName

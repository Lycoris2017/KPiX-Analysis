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

	outFile = open("missing_dieter.txt", 'w')

	eventsDieter = []
	eventsUwe = []

	for line in file1:
		fields =  re.split('=| |\n', line) ## splitting fields using regular expression
		fields = filter(None, fields) ## getting rid of empty fields
		eventsDieter.append(int(fields[0]))
	for line in file2:
		fields =  re.split('=| |\n', line) ## splitting fields using regular expression
		fields = filter(None, fields) ## getting rid of empty fields
		eventsUwe.append(int(fields[0]))

	#print eventsDieter
	#print eventsUwe
	count = 0
	for i in eventsDieter:
		if i not in eventsUwe:
			outFile.write('%i \r\n' %i)
			count += 1
	print count

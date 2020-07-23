#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import ROOT

def yParameterSensor(strip, sensor):
	return (50*(strip-920))*1E-3

def sensor2layer(sensor):
	if (sensor == 0):
		return 10
	if (sensor == 1):
		return 11
	if (sensor == 2):
		return 12
	if (sensor == 3):
		return 15
	if (sensor == 4):
		return 14
	if (sensor == 5):
		return 13



if __name__ == '__main__':
	inFile = '/scratch2/for-lycoris/survey_full.txt'
	outFile = '/scratch2/for-lycoris/run-895/run895_dieter.dat'
	fout = open(outFile, 'w')
	with open(inFile) as input:
		line = input.readline()
		lineNum = 0
		while line:
			#newline = [x.strip(' ') for x in line]
			#print line
			#print newline
			lineNum += 1
			fields =  re.split('=| |\n', line)
			newfields = filter(None, fields)
			Num=0
			fout.write('0 %s 0.0 \r\n' %newfields[0] )
			for Num, x in enumerate(newfields):
				if (Num%4 == 3 and Num >= 7): #Starting at field 7 every fourth field is the y position
					y=yParameterSensor(float(x), int((Num-7)/4))
#					print y
#					print (Num-8)/2
					#print x
					#print Num
					layer = sensor2layer(int((Num-7)/4))
					charge = float(newfields[Num+1])
					#print charge
					fout.write("%i 0. %f 0. %f %i %f %i \r\n" % (layer, y, 15, 2, charge, 0))
				fout.write(' ')
			line = input.readline()


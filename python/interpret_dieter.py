#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import ROOT

def yParameterSensor(strip, sensor):
#	if (sensor == 0 or sensor == 5 or sensor == 4): # kpix side showing towards beam movement beam side  KPIX >  < Beam
#		return (-50*(-strip+920))*1E-3
#	else:  # kpix side in direction of beam movement KPIX < < BEAM
	return (-50*(strip-920))*1E-3

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
	inFile = '/scratch/data/survey_xtra.txt'
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
			fields = re.split('=| |\n', line)
			newfields = filter(None, fields)
			Num=0

			for x in newfields:


				if (Num%2 == 1 and Num >= 8):
					y=yParameterSensor(float(x), int((Num-8)/2))
#					print y
#					print (Num-8)/2
#					print x
					layer = sensor2layer(int((Num-8)/2))
					fout.write("%i 0. %f 0. %f %i %f %i\n" % (layer, y, 15, 2, 3, 0))
				fout.write(' ')
				Num+=1
			fout.write('0 %i 0.0 \n' %lineNum)
			line = input.readline()


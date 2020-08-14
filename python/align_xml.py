#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import argparse
import xml.etree.ElementTree as ET
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument(
	'-x',
	'--xml',
	dest='xml',
	help='gear.xml file'
)
parser.add_argument(
	'-p',
	'--pede',
	dest='millepede',
	help='millepede.res file'
)
parser.add_argument(
	'-a',
	'--align',
	dest='align',
	help='which parameters to align'
)
args = parser.parse_args()



if __name__ == '__main__':
	xmlFile = args.xml
	dx = {}
	dy = {}
	dz = {}
	drot = {}
	#setting some parameters to 0 because they do not appear in the millepede.res file
	dy[12]=0
	dy[13]=0

	dz[10]=0
	dz[11]=0
	dz[12]=0
	dz[13]=0
	dz[14]=0
	dz[15]=0

	drot[10]=0
	drot[11]=0
	drot[12]=0
	drot[13]=0
	drot[14]=0
	drot[15]=0
	#opening the millepede file
	with open(args.millepede) as milleFile:
		for line in milleFile:
			if 'Parameter' in line:
				continue
			fields = line.split( )
			layer = int(fields[0])

			if '-1.0000' not in fields[2] and float(fields[1]) != 0.: #if the change is 0 or the parameter is static then fields[4] is empty so I just move on to the next line since there is nothing to be done

				if layer%10 == 1:
					if 2*float(fields[4]) >= abs(float(fields[1])): #if the error on the parameter is similar to the actual value (similar in this case being a factor 2 smaller) then do not perform a change
						dx[int(layer/10)] = 0
					else:
						dx[int(layer/10)] = float(fields[1]) #otherwise read in the value
				elif layer%10 == 2:
					if 2*float(fields[4]) >= abs(float(fields[1])):
						dy[int(layer/10)] = 0
					else:
						dy[int(layer/10)] = float(fields[1])
				elif layer%10 == 3:
					if 2*float(fields[4]) >= abs(float(fields[1])):
						dz[int(layer/10)] = 0
					else:
						dz[int(layer/10)] = float(fields[1])
				elif layer%10 == 6:
					if 2*float(fields[4]) >= abs(float(fields[1])):
						drot[int(layer/10)] = 0
					else:
						drot[int(layer/10)] = float(fields[1])*(180./1000 /np.pi) #millepede.res is in mrad and gear file is in degrees
			else: # since I need every dictionary to have an entry I just set them to 0 if no changes are to be made.
				if layer%10 == 1:
					dx[int(layer/10)] = 0
				if layer%10 == 2:
					dy[int(layer/10)] = 0
				if layer%10 == 3:
					dz[int(layer/10)] = 0
				if layer%10 == 6:
					drot[int(layer/10)] = 0

	print(dx,dy,dz)

	outFile = args.xml[:-4]+"_new.xml"
	fout = open(outFile, 'w')
	tree = ET.parse(xmlFile)
	root = tree.getroot()
	print(root.tag)
	for level1 in root: #start of a massive amount of for loops simply to loop down to the ladder level in the xml
		if ('detectors' in level1.tag):
			for level2 in level1:
				for level3 in level2:
					print(level3.tag)
					if ('layers' in level3.tag):
						for level4 in level3:
							for level5 in level4:
								if 'ladder' in level5.tag:
									#print(level5.tag)
									#print(level5.attrib)
									layer = int(level5.get('ID'))
									x = float(level5.get('positionX'))
									y = float(level5.get('positionY'))
									z = float(level5.get('positionZ'))
									rot = float(level5.get('rotationXY'))
									if 'x' in args.align:
										x = x-dx[layer]
										level5.set('positionX', str(x)) #setting the new value to the corrected value
									if 'y' in args.align:
										print(layer)
										y = y-dy[layer]
										level5.set('positionY', str(y))
									if 'z' in args.align:
										z = z-dz[layer]
										level5.set('positionZ', str(z))
									if 'rot' in args.align:
										rot = rot-drot[layer]
										level5.set('rotationXY', str(rot))
	tree.write(outFile)
	print(outFile)
#	with open(inFile) as input:
#		line = input.readline()
#		lineNum = 0
#		while line:
#			#newline = [x.strip(' ') for x in line]
#			#print line
#			#print newline
#			lineNum += 1
#			fields =  re.split('=| |\n', line)
#			newfields = filter(None, fields)
#			Num=0
#			fout.write('0 %s 0.0 \r\n' %newfields[0] )
#			for Num, x in enumerate(newfields):
#				if (Num%4 == 3 and Num >= 7): #Starting at field 7 every fourth field is the y position
#					y=yParameterSensor(float(x), int((Num-7)/4))
##					print y
##					print (Num-8)/2
#					#print x
#					#print Num
#					layer = sensor2layer(int((Num-7)/4))
#					charge = float(newfields[Num+1])
#					#print charge
#					fout.write("%i 0. %f 0. %f %i %f %i \r\n" % (layer, y, 15, 2, charge, 0))
#				fout.write(' ')
#			line = input.readline()


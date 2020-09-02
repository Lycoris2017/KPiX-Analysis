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
					if 1.1*float(fields[4]) >= abs(float(fields[1])): #if the error on the parameter is similar to the actual value (similar in this case being less than a factor 1.1 smaller) then do not perform a change
						dx[int(layer/10)] = 0
					else:
						dx[int(layer/10)] = float(fields[1]) #otherwise read in the value
				elif layer%10 == 2:
					if 1.1*float(fields[4]) >= abs(float(fields[1])):
						dy[int(layer/10)] = 0
					else:
						dy[int(layer/10)] = float(fields[1])
				elif layer%10 == 3:
					if 1.1*float(fields[4]) >= abs(float(fields[1])):
						dz[int(layer/10)] = 0
					else:
						dz[int(layer/10)] = float(fields[1])
				elif layer%10 == 6:
					if 1.1*float(fields[4]) >= abs(float(fields[1])):
						drot[int(layer/10)] = 0
					else:
						drot[int(layer/10)] = float(fields[1])*(180. /np.pi) #millepede.res is in rad and gear file is in degrees
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
	if (all(value == 0 for value in dx.values()) and all(value == 0 for value in dy.values()) and all(value == 0 for value in dz.values()) and all(value == 0 for value in drot.values())):
		print("No corrections being applied!")
		sys.exit(2)
	ending = args.xml.find("_step")
	if (ending == -1):
		outFile = args.xml[:-4]+"_step1.xml"
	else:
		number = int(args.xml[ending+5:-4])
		number = number+1
		outFile = args.xml[:ending]+"_step"+str(number)+".xml"
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


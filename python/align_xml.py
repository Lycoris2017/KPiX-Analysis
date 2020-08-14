#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import argparse
import xml.etree.ElementTree as ET


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
	with open(args.millepede) as milleFile:
		for line in milleFile:
			if 'Parameter' in line:
				continue
			fields = line.split( )
			layer = int(fields[0])
			if layer%10 == 1:
				dx[int(layer/10)] = float(fields[1])
			elif layer%10 == 2:
				dy[int(layer/10)] = float(fields[1])
			elif layer%10 == 3:
				dz[int(layer/10)] = float(fields[1])
			elif layer%10 == 6:
				drot[int(layer/10)] = float(fields[1])
	dy[12]=0
	dy[13]=0
	print(dx,dy,dz)

	outFile = args.xml[:-4]+"_new.xml"
	fout = open(outFile, 'w')
	tree = ET.parse(xmlFile)
	root = tree.getroot()
	print(root.tag)
	for level1 in root:
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
										level5.set('positionX', str(x))
									if 'y' in args.align:
										print(layer)
										y = y-dy[layer]
										level5.set('positionY', str(y))
									if 'z' in args.align:
										z = z-dz[layer]
										level5.set('positionZ', z)
									if 'rot' in args.align:
										rot = rot-drot[layer]
										level5.set('rotationXY', rot)
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


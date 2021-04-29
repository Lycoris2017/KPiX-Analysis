#!usr/bin/python

import numpy as np
import string
import argparse
import argcomplete
from operator import add
import sys
import re
import ROOT
import queue

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("-c",
                    "--corry",
                    dest='corry',
                    help="corry textwriter file"
                    )
parser.add_argument("-k",
                    "--kpix",
                    dest='kpix',
                    help="kpix cluster file with trigger timestamps .txt"
                    )
args = parser.parse_args()

outfile = open('test_time.dat', 'w')




outHistFile = ROOT.TFile.Open('test_time.root', "RECREATE")
outHistFile.cd()


syncedTime = ROOT.TH1F("matched_time",
                       "synced_events_v_time; time (s); No. of Entries ",
                       750,
                       #-50e-2, 50e-2
                       -0.5, 2999.5)


eventTimeKpix = ROOT.TH1F("eventTimeKpix",
                       "eventTimeKpix; time (s); No. of Entries ",
                       750,
                       #-50e-2, 50e-2
                       -0.5, 2999.5)

eventTimeTpx = ROOT.TH1F("eventTimeTpx3",
                       "eventTimeTpx3; time (s); No. of Entries ",
                       750,
                       #-50e-2, 50e-2
                       -0.5, 2999.5)
MimosaFrameTime = ROOT.TH1F("MimosaFrameTime",
                            "MimosaFrameTime; time (us); No. of Entries",
                            1000,
                            0, 1000)
minTimeDiff = ROOT.TH1F("minTimeDiff",
                            "minTimeDiff; time (us); No. of Entries",
                            200,
                            -1000, 1000)


eventCounter = 0
timeStart = queue.Queue()
timeEnd = queue.Queue()
for line in iter(open(args.corry)):
	if ('===' in line) or ('Start:' in line) or ('End:' in line):
		fields =  re.split(r',| |[(|)]|\[|\]|\r\n|\n', line) ## splitting fields using regular expression
		fields = list(filter(None, fields)) ## getting rid of empty fields
		if ('===' in fields):
			maxEvt = int(fields[1])
		elif ('Start:' in fields):
			timeStart.put(float(fields[1]))
			eventTimeTpx.Fill(float(fields[1])/1E9)
		elif 'End:' in fields:
			timeEnd.put(float(fields[1]))
	else:
		continue
print("Found ", maxEvt, " Events in TPX3 file")
eventCounter = 0
time = 0
syncedEvents = 0
test = 0
for i in range(len(timeStart)):
        timeFrame = abs(timeStart[i]-timeEnd[i])
        MimosaFrameTime.Fill(timeFrame/1000.)
                                

with open(args.kpix) as inFile:
	for line in inFile:
		fields =  re.split(r',| |[(|)]|\[|\]|\r\n|\n', line) ## splitting fields using regular expression
		fields = list(filter(None, fields)) ## getting rid of empty fields
		#print fields
#		if (syncedEvents == 1 and test == 1):
#			print fields
		if (eventCounter > maxEvt):
			print ("End of TPX3 data. Stopping")
			break
		if ('runtime' in fields):
			continue
		if time != float(fields[7]):
			if (test):
				outfile.write('6566 '+str(eventCounter)+' 0.0 \r\n')
			time = float(fields[7])
			eventTimeKpix.Fill(time/1E9)
			tDiff = float("Inf")
			
			while(timeStart != queue.Empty && timeEnd != queue.Empty
			
			for i in range(len(timeStart[eventCounter:])):
				if abs(tDiff) > abs(time-timeStart[eventCounter+i]):
					tDiff = time-timeStart[eventCounter+i]
				if time >= timeStart[eventCounter+i] and time <= timeEnd[eventCounter+i]:
					eventCounter = eventCounter + i
					sys.stdout.write('\r' + str(eventCounter)+'/'+str(maxEvt)+" | KpiX Event Number "+fields[0])
					sys.stdout.flush()
					syncedEvents += 1
					syncedTime.Fill(time/1E9)
					test = 1
					break
				else: test = 0
			minTimeDiff.Fill(tDiff/1000.)
		if (test):
			for i in fields:
				outfile.write(i+' ')
			outfile.write(str(eventCounter)+' ')
			outfile.write('\r\n')
print('\n', syncedEvents)
outHistFile.Write()

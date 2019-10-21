#!/usr/bin/env python

import ROOT

histobjs=[]
slope_digits=5

def loopdir(keys):
    
    for key_object in keys:
        if ('TDirectory' in key_object.GetClassName()):
            loopdir( key_object.ReadObj().GetListOfKeys() )
        else:
            if ("slope_vs_channel_k" in key_object.GetName()
                and "_b0" in key_object.GetName() ):
                histobjs.append(key_object)
                #print key_object.GetName()
        #endif
    #endfor
    return 

def PrintSlopes():
    if not histobjs:
        return

    print "kpix,channel,bucket,slope"
    for  histobj in histobjs:
        hist = histobj.ReadObj()
        
        kpix=hist.GetName().strip("slope_vs_channel_k")
        kpix=kpix.rstrip('0')
        kpix=kpix.strip("_b")
        bucket = 0
        for channel in range(0, 1024):
            slope = hist.GetBinContent(channel+1)
            slope = round(slope, slope_digits)
            print "{},{},{},{}".format(kpix,channel,bucket,slope)
        
def main():
    #print "Hello world"
    #fname="/opt/data/calibration/Calibration_20190724_182714.dat.ymlcalib.root"
    fname="/opt/data/calibration/Calibration_20190724_182714.dat.ymlcalib.root"
    rfile=ROOT.TFile(fname)
    keys = rfile.GetListOfKeys()

    loopdir(keys)
    PrintSlopes()

    
if __name__ == '__main__':
    main()

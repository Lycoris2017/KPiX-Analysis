#!/bin/sh

#g++ -o StandardPlotter StandardPlotter.cxx `root-config --cflags --glibs`
g++ -o StandardPlotter StandardPlotter.cxx  `root-config --libs --glibs --cflags --ldflags --auxlibs --auxcflags` -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore


# Variables
CC=g++ -std=c++11

# exe sources:
SRC=$(PWD)/src
INC=$(PWD)/include
binaries=$(patsubst $(SRC)/%.cxx,$(BIN)/%,$(wildcard $(SRC)/*.cxx))

# external sources (kpix etc):
EXT_DIR=$(PWD)/extern
EXT_OBJ=$(patsubst $(EXT_DIR)%.cpp,$(OBJ)/%.o,$(wildcard $(EXT_DIR)/*.cpp))

#LDIR=$(PWD)/lib
#KP_CFLAGS=-I/home/mengqing/kpix/kpix -I/home/mengqing/kpix/generic
#KP_LFLAGS=-L$(LDIR) -lkpix 
#KP_RLIB=-Wl,-rpath=$(LDIR)

KP_RLIB=-Wl,-rpath=/usr/local/lib/kpix/
KP_CFLAGS=-I/usr/local/include/kpix
KP_LFLAGS=-L/usr/local/lib/kpix/ -lkpix

CFLAGS=-Wall $(shell xml2-config --cflags) $(shell root-config --cflags) $(KP_CFLAGS) -I$(INC) -I$(PWD)/extern
LFLAGS=$(shell xml2-config --libs) $(shell root-config --glibs) $(KP_LFLAGS)

# dir:
OBJ=$(PWD)/.obj
BIN=$(PWD)/bin

#default:
all: dir map LycorisDict.cxx $(EXT_OBJ) $(binaries)
ana: dir map analysis


LycorisDict.cxx: $(INC)/ntupleMaker.h $(INC)/LinkDef.h
	rootcint -f $@ -c $(CFLAGS) -p $^

$(OBJ)/%.o: $(EXT_DIR)/%.cpp $(EXT_DIR)/%.h 
	$(CC) -c $(CFLAGS) $(KP_RLIB) -o $@ $<

$(BIN)/%: $(SRC)/%.cxx $(EXT_OBJ)
	$(CC) -o $@ $< $(CFLAGS) $(KP_RLIB) $(LFLAGS) 

#analysis: $(SRC)/analysis.cxx
#	$(CC) -o $(BIN)/analysis $(SRC)/analysis.cxx $(CFLAGS) $(KP_RLIB) $(LFLAGS) 

# executable directory
dir:
	test -d $(BIN) || mkdir $(BIN)
	test -d $(OBJ) || mkdir $(OBJ)

map:
	python python/stripMapMaker.py data/tracker_to_kpix_left.dat data/tracker_to_kpix_right.dat include/kpix_left_and_right.h

clean:
	rm -r ./.obj
	rm ./bin/*

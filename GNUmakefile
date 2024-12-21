#############################################
#                                           #
############  makefile  #####################
#                                           #
#############################################
# Author: Yuzhi Che 2022.10
# Modified: Jiaxuan Wang 2023.4 
#
#
##set up control#######
TOP:= $(shell pwd)
DIR_OBJ:= $(TOP)/obj
DIR_BIN:= $(TOP)/bin
DIR_SRC:= $(TOP)/src
DIR_INCLUDE:= $(TOP)/include


SRC = $(wildcard $(DIR_SRC)/*.cxx) 	
OBJ = $(patsubst $(DIR_SRC)/%.cxx,$(DIR_OBJ)/%.o,$(SRC)) 
#ALL:
#	@echo $(SRC)
#	@echo $(OBJ)

##ROOT#######################################
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMinuit
ROOTGLIBS     = $(shell root-config --glibs)

CPPLIBS = $(ROOTLIBS) $(ROOTGLIBS) #-ltbb
DIR_INCLUDE+=$(ROOTCFLAGS)

# #set up compilers#
CXX = g++
CPPFLAGS = -g -Wall -I$(DIR_INCLUDE)

####### Make Execuatbles#######
all: pedestal calib mip hl_ratio auto_decode hl_decode extract # decode track_fit

# decode: main/decode.cxx $(OBJ) $(DIR_BIN)
#	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

hl_decode: main/hl_decode.cxx $(OBJ) $(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

auto_decode: main/auto_decode.cxx $(OBJ) $(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

pedestal: main/pedestal.cxx $(OBJ) $(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

# track_fit: main/track_fit.cxx $(OBJ) $(DIR_BIN)
# 	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS) -o $(DIR_BIN)/$(notdir $@)

mip: main/mip.cxx $(OBJ) $(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

calib: main/calib.cxx $(OBJ) $(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

hl_ratio: main/hl_ratio.cxx $(OBJ) $(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)

extract: main/extract.cxx $(OBJ) $(DIR_BIN)
	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notdir $@)


# event_display: event_display.cxx $(OBJ) #$(DIR_BIN)
# 	$(CXX) $(CPPFLAGS) $< $(OBJ) $(CPPLIBS)  -o $(DIR_BIN)/$(notidr $@)
#################################
$(DIR_OBJ)/%.o:$(DIR_SRC)/%.cxx  
	$(CXX) $(CPPFLAGS)  -c $(DIR_SRC)/$(notdir $<)  -o $(DIR_OBJ)/$(notdir $@)

#	$(AR) r $(DIR_SRC)/$(notdir $@) $(DIR_OBJ)/$(notdir $@)
#	$(RM) $(DIR_OBJ)/$(notdir $@)
#@echo 

.PHONY:clean
clean: 
	rm -f $(DIR_OBJ)/*.o rm -f $(DIR_BIN)/*

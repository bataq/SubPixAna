TARGET = SubPixAna
SRCS = $(TARGET).cpp
OBJS = $(TARGET).o

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
BOOSTOPTFLAG = -lboost_program_options

CXXFLAGS = $(ROOTCFLAGS) -Wall -fPIC
CXXLIBS = $(ROOTLIBS)
CC = g++

$(TARGET): $(OBJS)
	$(CC) $(CXXLIBS) $(OBJS) -o $@ $(BOOSTOPTFLAG)
.cpp.o:
	$(CC) $(CXXFLAGS) -c $<
clean:
	rm -rf $(TARGET)
	rm -rf $(OBJS)

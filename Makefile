
########################################################################
#   Preamble.
########################################################################

ifeq ($(OPTFLAGS),)
	OPTFLAGS = -g
endif

CXX = g++
JAVAC = javac
CPPFLAGS = -I/user1/faculty/stanchen/pub/boost
CXXFLAGS = -Wall
LDFLAGS = $(OPTFLAGS)
LDLIBS = -lm

CXXFLAGS += $(OPTFLAGS)

#   GNU make's default rule uses $(CC) for linking
LINK.o = $(CXX) $(LDFLAGS) $(TARGET_ARCH)


########################################################################
#   Rules.
########################################################################

all : lab1 lab1_dtw Lab1

clean:
	rm -f lab1 lab1_dtw *.o Lab1.class FrontEnd.class Lab1_DTW.class

util.o : util.C util.H

front_end.o : front_end.C front_end.H util.H

lab1.o : lab1.C util.H front_end.H

lab1_dtw.o : lab1_dtw.C util.H

lab1 : lab1.o util.o front_end.o

lab1_dtw : lab1_dtw.o util.o

Lab1 :
	$(JAVAC) Lab1.java FrontEnd.java

Lab1_DTW :
	$(JAVAC) Lab1_DTW.java FrontEnd.java


########################################################################
#   
########################################################################



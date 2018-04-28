
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

FRONT_END_OBJ = /user1/faculty/stanchen/e6870/lab2data/front_end.o


########################################################################
#   Rules.
########################################################################

all : lab2_vit lab2_train lab2_fb Lab2Vit Lab2Train Lab2Fb

clean:
	rm -f lab2_vit lab2_train lab2_fb *.o \
	    Lab2Vit.class Lab2Train.class Lab2Fb.class GmmStats.class

util.o : util.C util.H

gmm_util.o : gmm_util.C gmm_util.H util.H

main.o : main.C

lab2_vit.o : lab2_vit.C lab2_vit.H util.H front_end.H

lab2_train.o : lab2_train.C util.H front_end.H gmm_util.H

lab2_fb.o : lab2_fb.C util.H front_end.H gmm_util.H

lab2_vit : lab2_vit.o util.o $(FRONT_END_OBJ) main.o

lab2_train : lab2_train.o util.o $(FRONT_END_OBJ) gmm_util.o main.o

lab2_fb : lab2_fb.o util.o $(FRONT_END_OBJ) gmm_util.o main.o

Lab2Vit :
	$(JAVAC) Lab2Vit.java

Lab2Train :
	$(JAVAC) Lab2Train.java GmmStats.java

Lab2Fb :
	$(JAVAC) Lab2Fb.java GmmStats.java


########################################################################
#   
########################################################################



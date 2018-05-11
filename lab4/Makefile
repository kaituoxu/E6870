
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

all : lab4_vit Lab4Vit

clean:
	rm -f lab4_vit front_end.o lab4_vit.o util.o Lab4Vit.class

util.o : util.C util.H

front_end.o : front_end.C front_end.H util.H

main.o : main.C

lab4_vit.o : lab4_vit.C lab4_vit.H util.H front_end.H

lab4_vit : lab4_vit.o util.o front_end.o gmm_set_attila.o main.o libattila.so /usr/lib64/libpthread.so

Lab4Vit :
	$(JAVAC) Lab4Vit.java


########################################################################
#   
########################################################################



OUTEXE=$(ODIR)/nui
INCLUDEDIRS =. ../mfl /usr/local/include/ 
LIBDIRS=
LIBS=termcap
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
ARCHIVES = ../mfl/$(MFLLIB) /usr/local/lib/ncl/libncl.a /usr/local/lib/libedit.a
include ../Makefile.inc
-include $(DEPS)

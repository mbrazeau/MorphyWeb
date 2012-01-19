CC=gcc
CPP=g++
AR=ar
CFLAGS=$(patsubst %,-I%,$(INCLUDEDIRS)) -Wall
ODIR=obj
MFLLIB=$(ODIR)/libmfl.a

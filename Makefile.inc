INCLUDEDIRS =. ../mfl
CC=gcc
AR=ar
CFLAGS=$(patsubst %,-I%,$(INCLUDEDIRS)) -Wall
ODIR=obj
MFLLIB=$(ODIR)/libmfl.a
TUIEXE=$(ODIR)/morphy



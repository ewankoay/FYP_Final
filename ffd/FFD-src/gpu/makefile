#######################################################
# LINUX OPERATING SYSTEMS
#
# Makefile to compile FFD gpu code
# 
# This was written by Cary Faulkner based on a
# previous makefile written by Michael Wetter
#
# Michael Wetter (MWetter@lbl.gov) October 24, 2012
#######################################################
SHELL = /bin/sh
ARCH = $(shell getconf LONG_BIT)
CL_DIR=$(CURC_CUDA_INC)
# Directory where executable will be copied to
BINDIR = ../../Library/linux$(ARCH)

#######################################################
## Compilation flags
CC = gcc

#Note that Dymola use 32bit compiler, so generated executable only support 32bit loaded library
CC_FLAGS_32 = -Wall -lm -m32 -std=c89 -pedantic -msse2 -mfpmath=sse
CC_FLAGS_64 = -Wall -lm -m64 -std=c89 -pedantic -msse2 -mfpmath=sse

SRCS = advection.c boundary.c chen_zero_equ_model.c \
       data_writer.c diffusion.c ffd.c ffd_data_reader.c geometry.c initialization.c \
       interpolation.c parameter_reader.c projection.c sci_reader.c solver.c solver_gs.c \
       timing.c utility.c opencl_main.c

OBJS = advection.o boundary.o chen_zero_equ_model.o \
       data_writer.o diffusion.o ffd.o ffd_data_reader.o geometry.o initialization.o \
       interpolation.o parameter_reader.o projection.o sci_reader.o solver.o solver_gs.o \
       timing.o utility.o opencl_main.o

LIBS = -lpthread -lc -lm -lOpenCL

# Note that -fPIC is recommended on Linux according to the Modelica specification

all: clean
	$(CC) -o ffd_gpu $(SRCS) $(LIBS) -I$(CL_DIR)
	rm -f $(OBJS)

clean:
	rm -f $(OBJS) 

# To enable RootMakefile, add fellow empty targets
doc:
cleandoc:

#######################################################
# LINUX/MAC OPERATING SYSTEMS
#######################################################
# Makefile to build FFD executable using GCC
# Wei Tian (Wei.Tian@Schneider-Electric.com) April 13, 2018
# Michael Wetter (MWetter@lbl.gov) October 24, 2012
#######################################################
SHELL = /bin/sh
ARCH = $(shell getconf LONG_BIT)

## Name of FFD executable
LIB_C_MAC = FFD_C
LIB_C_LINUX = FFD_C
LIB_O_MAC = FFD_OPENCL
LIB_O_LINUX = FFD_OPENCL
LIB_O_DETECT_LINUX = DEVICE_DETECT

# Directory where executable will be copied to
BINDIR_C_MAC = ./Build_C
BINDIR_C_LINUX = ./Build_C
BINDIR_O_MAC = ./Build_OPENCL
BINDIR_O_LINUX = ./Build_OPENCL
BINDIR_O_DETECT = ./Build_OPENCL_Detect

#######################################################
## Compilation flags
## -w: suppress warnings
CC = gcc
CC_FLAGS_32 = -Wall -lm -m32 -std=c99 -pedantic -msse2 -mfpmath=sse -w
CC_FLAGS_64 = -Wall -lm -m64 -std=c99 -pedantic -msse2 -mfpmath=sse -w

SRCS_C = advection.c boundary.c chen_zero_equ_model.c \
       data_writer.c diffusion.c ffd.c ffd_data_reader.c geometry.c initialization.c \
       interpolation.c parameter_reader.c projection.c sci_reader.c solver.c solver_gs.c \
       timing.c utility.c main.c

OBJS_C = advection.o boundary.o chen_zero_equ_model.o \
       data_writer.o diffusion.o ffd.o ffd_data_reader.o geometry.o initialization.o \
       interpolation.o parameter_reader.o projection.o sci_reader.o solver.o solver_gs.o \
       timing.o utility.o main.o

SRCS_O = advection.c boundary.c chen_zero_equ_model.c \
      data_writer.c diffusion.c ffd.c ffd_data_reader.c geometry.c initialization.c \
      interpolation.c parameter_reader.c projection.c sci_reader.c solver.c solver_gs.c \
      timing.c utility.c opencl_solver.c

OBJS_O = advection.o boundary.o chen_zero_equ_model.o \
      data_writer.o diffusion.o ffd.o ffd_data_reader.o geometry.o initialization.o \
      interpolation.o parameter_reader.o projection.o sci_reader.o solver.o solver_gs.o \
      timing.o utility.o opencl_solver.o

SRCS_O_R = $(SRCS_O) opencl_deploy_main.c

OBJS_O_R = $(OBJS_O) opencl_deploy_main.o

SRCS_O_T = $(SRCS_O) opencl_dev_main.c

OBJS_O_T = $(OBJS_O) opencl_dev_main.o

SRCS_DETECT_O = platform_device_detect.c

OBJS_DETECT_O = platform_device_detect.o

## Build C executable on Mac
FFD_C_MAC:
	## framework OpenCL GL and GLUT are needed for compiling the visualization codes
	## delete and recreate the exisiting Directory
	rm -f $(OBJS_C) $(BINDIR_C_MAC)$(LIB_C_MAC)
	if [ -d $(BINDIR_C_MAC) ]; then rm -rf $(BINDIR_C_MAC); fi
	mkdir -p $(BINDIR_C_MAC)
	## start compiling
	$(CC) $(CC_FLAGS_64) -c $(SRCS_C)
	$(CC) -m64 -o $(LIB_C_MAC) $(OBJS_C)
	rm -f $(OBJS_C)
	mv $(LIB_C_MAC) $(BINDIR_C_MAC)
	@echo "==== library generated in $(BINDIR_C_MAC)"

## Compile C executable on Cygwin on Windows
FFD_C_LINUX:
	## Due to a bug in the Cygwin, the linker should be put at end of the command
	## delete and recreate the exisiting Directory
	rm -f $(OBJS_C) $(BINDIR_C_LINUX)$(LIB_C_LINUX)
	if [ -d $(BINDIR_C_LINUX) ]; then rm -rf $(BINDIR_C_LINUX); fi
	mkdir -p $(BINDIR_C_LINUX)
	## start compiling
	$(CC) $(CC_FLAGS_64) -c $(SRCS_C)
	$(CC) -m64 -o $(LIB_C_LINUX) $(OBJS_C) -lm
	rm -f $(OBJS_C)
	mv $(LIB_C_LINUX) $(BINDIR_C_LINUX)
	@echo "==== library generated in $(BINDIR_C_LINUX)"

## Build OpenCL executable on Mac
FFD_OPENCL_MAC:
	rm -f $(OBJS_O_R) $(BINDIR_O_MAC)$(LIB_O_MAC)
	if [ -d $(BINDIR_O_MAC) ]; then rm -rf $(BINDIR_O_MAC); fi
	mkdir -p $(BINDIR_O_MAC)
	$(CC) $(CC_FLAGS_64) -c $(SRCS_O_R)
	$(CC) -m64 -framework OpenCL -o $(LIB_O_MAC) $(OBJS_O_R)
	rm -f $(OBJS_O_R)
	mv $(LIB_O_MAC) $(BINDIR_O_MAC)
	@echo "==== library generated in $(BINDIR_O_MAC)"

## Compile OpenCL executable on Cygwin on Windows
FFD_OPENCL_LINUX:
	rm -f $(OBJS_O_R) $(BINDIR_O_LINUX)$(LIB_O_LINUX)
	if [ -d $(BINDIR_O_LINUX) ]; then rm -rf $(BINDIR_O_LINUX); fi
	mkdir -p $(BINDIR_O_LINUX)
	$(CC) $(CC_FLAGS_64) -c $(SRCS_O_R)
	$(CC) -m64 -o $(LIB_O_LINUX) $(OBJS_O_R) -lm -lOpenCL
	rm -f $(OBJS_O_R)
	mv $(LIB_O_LINUX) $(BINDIR_O_LINUX)
	@echo "==== library generated in $(BINDIR_O_LINUX)"

## Compile an OpenCL program to detect the platform and devices
DEVICE_DETECT_LINUX:
	rm -f $(OBJS_DETECT_O) $(BINDIR_O_DETECT)$(LIB_O_DETECT_LINUX)
	if [ -d $(BINDIR_O_DETECT) ]; then rm -rf $(BINDIR_O_DETECT); fi
	mkdir -p $(BINDIR_O_DETECT)
	$(CC) $(CC_FLAGS_64) -c $(SRCS_DETECT_O)
	$(CC) -m64 -o $(LIB_O_DETECT_LINUX) $(OBJS_DETECT_O) -lOpenCL
	rm -f $(OBJS_DETECT_O)
	mv $(LIB_O_DETECT_LINUX) $(BINDIR_O_DETECT)
	@echo "==== library generated in $(BINDIR_O_DETECT)"

## Compile an OpenCL program to detect the platform and devices
DEVICE_DETECT_MAC:
	rm -f $(OBJS_DETECT_O) $(BINDIR_O_DETECT)$(LIB_O_DETECT_LINUX)
	if [ -d $(BINDIR_O_DETECT) ]; then rm -rf $(BINDIR_O_DETECT); fi
	mkdir -p $(BINDIR_O_DETECT)
	$(CC) $(CC_FLAGS_64) -c $(SRCS_DETECT_O)
	$(CC) -m64 -framework OpenCL -o $(LIB_O_DETECT_LINUX) $(OBJS_DETECT_O)
	rm -f $(OBJS_DETECT_O)
	mv $(LIB_O_DETECT_LINUX) $(BINDIR_O_DETECT)
	@echo "==== library generated in $(BINDIR_O_DETECT)"

## Clean build folders
CLEAN:
	if [ -d $(BINDIR_C_MAC) ]; then rm -rf $(BINDIR_C_MAC); fi
	if [ -d $(BINDIR_C_LINUX) ]; then rm -rf $(BINDIR_C_LINUX); fi
	if [ -d $(BINDIR_O_MAC) ]; then rm -rf $(BINDIR_O_MAC); fi
	if [ -d $(BINDIR_O_LINUX) ]; then rm -rf $(BINDIR_O_LINUX); fi
	if [ -d $(BINDIR_O_DETECT) ]; then rm -rf $(BINDIR_O_DETECT); fi

## Compile prorgrams to detect platform and device information
DEVICE_DETECT:
ifeq ($(shell uname -s),Darwin)
	$(info mac os system is detected and will build mac exe)
	$(MAKE) DEVICE_DETECT_MAC
else ifeq ($(shell uname -s),Linux)
	$(info linux system is detected and will build linux exe)
	$(MAKE) DEVICE_DETECT_LINUX
else ifeq ($(shell uname -o),Cygwin)
	$(info cygwin system is detected and will build cygwin exe)
	$(MAKE) DEVICE_DETECT_LINUX
else ifeq ($(OS),Windows_NT)
	$(info this makefile would not work in Windows. Will exit)
endif

## Compile FFD executables from C and OpenCL source codes
FFD:
ifeq ($(shell uname -s),Darwin)
	$(info mac os system is detected and will build mac exe)
	$(MAKE) FFD_OPENCL_MAC
	$(MAKE) FFD_C_MAC
else ifeq ($(shell uname -s),Linux)
	$(info linux system is detected and will build linux exe)
	$(MAKE) FFD_OPENCL_LINUX
	$(MAKE) FFD_C_LINUX
else ifeq ($(shell uname -o),Cygwin)
	$(info cygwin system is detected and will build cygwin exe)
	$(MAKE) FFD_OPENCL_LINUX
	$(MAKE) FFD_C_LINUX
else ifeq ($(OS),Windows_NT)
 	$(info this makefile would not work in Windows. Will exit)
endif

## Compile programs to detect platform and device information and FFD executables
ALL: DEVICE_DETECT FFD

##########################################################################
## Compile OpenCL executable on Cygwin on Windows (test)
FFD_OPENCL_LINUX_TEST:
	rm -f $(OBJS_O_T) $(BINDIR_O_LINUX)$(LIB_O_LINUX)
	if [ -d $(BINDIR_O_LINUX) ]; then rm -rf $(BINDIR_O_LINUX); fi
	mkdir -p $(BINDIR_O_LINUX)
	$(CC) $(CC_FLAGS_64) -c $(SRCS_O_T)
	$(CC) -m64 -o $(LIB_O_LINUX) $(OBJS_O_T) -lOpenCL
	rm -f $(OBJS_O_T)
	mv $(LIB_O_LINUX) $(BINDIR_O_LINUX)
	@echo "==== library generated in $(BINDIR_O_LINUX)"

## Compile an OpenCL program to detect the platform and devices (test)
FFD_OPENCL_MAC_TEST:
	rm -f $(OBJS_DETECT_O_T) $(BINDIR_O_DETECT)$(LIB_O_DETECT_LINUX)
	if [ -d $(BINDIR_O_DETECT) ]; then rm -rf $(BINDIR_O_DETECT); fi
	mkdir -p $(BINDIR_O_DETECT)
	$(CC) $(CC_FLAGS_64) -c $(SRCS_DETECT_O_T)
	$(CC) -m64 -framework OpenCL -o $(LIB_O_DETECT_LINUX) $(OBJS_DETECT_O_T)
	rm -f $(OBJS_DETECT_O_T)
	mv $(LIB_O_DETECT_LINUX) $(BINDIR_O_DETECT)
	@echo "==== library generated in $(BINDIR_O_DETECT)"

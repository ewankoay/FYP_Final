///////////////////////////////////////////////////////////////////////////////
///
/// \file   main.c
///
/// \brief  main entrance of parallel FFD program
///
/// \author Wei Tian
///         University of Miami, Schneider Electric
///         w.tian@umiami.edu, Wei.Tian@Schneider-Electric.com
///         Thomas Sevilla
///         University of Miami
///         t.sevilla@umiami.edu
///
/// \date   7/07/2017
///
///\ This codes have been extensively validated and for the results please refer
///   to the journal paper:
///   Tian, Wei, Thomas Alonso Sevilla, and Wangda Zuo. "A systematic
///   evaluation of accelerating indoor airflow simulations using cross-platform parallel computing."
///   Journal of Building Performance Simulation 10.3 (2017): 243-255.
///
///\  All RIGHTS RESERVED.
///////////////////////////////////////////////////////////////////////////////
#ifndef _OPENCL_MAIN_H
#define _OPENCL_MAIN_H

#include "timing.h"
#include "data_writer.h"
#include "initialization.h"
#include "geometry.h"
#include "utility.h"
#include "solver.h"
#include "kernel.h"
#include "data_structure.h"
#include "sci_reader.h"

#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#ifdef __APPLE__ 
  #include <OpenCL/opencl.h>
#elif __linux__
  #include <CL/opencl.h>	
#else
  #include <CL\cl.h>
#endif

#ifdef _MSC_VER
  #define _CRT_SECURE_NO_WARNINGS
#endif

#ifdef _MSC_VER
  #include <windows.h>
#elif defined __GNUC__
  #ifdef __WIN64 || __WIN32
    #include <windows.h>
  #elif __APPLE__ || __linux__
    #include <unistd.h>
  #else
    #include <unistd.h>
  #endif
#else
  #include <unistd.h>
#endif

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

///////////////////////////////////////////////////////////////////////////////
/// Main routine of FFD
///
///\para coupled simulation Integer to identify the simulation type
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int ffd_prep(int cosimulation);

///////////////////////////////////////////////////////////////////////////////
/// Allocate memory for variables
///
///\param para Pointer to FFD parameters
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
int allocate_memory_opencl(PARA_DATA *para);

#endif

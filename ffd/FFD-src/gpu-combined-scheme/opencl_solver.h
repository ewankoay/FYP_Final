///////////////////////////////////////////////////////////////////////////////
///
/// \file   opencl_solver.c
///
/// \brief  main entrance of parallel FFD program
///
/// \author Wei Tian
///         University of Miami, Schneider Electric
///         w.tian@umiami.edu, Wei.Tian@Schneider-Electric.com
///         Thomas Sevilla
///         University of Miami
///         t.sevilla@umiami.edu
///////////////////////////////////////////////////////////////////////////////
#ifndef _OPENCL_SOLVER_H
#define _OPENCL_SOLVER_H

#include "data_structure.h"
#include "data_writer.h"
#include "geometry.h"
#include "initialization.h"
#include "sci_reader.h"
#include "solver.h"
#include "timing.h"
#include "utility.h"


#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#elif defined __linux__
#include <CL/cl.h>
#else
#include <CL/cl.h>
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

void mypause(int x);

int ffd_prep(int cosimulation);

int allocate_memory_opencl(PARA_DATA *para);

int create_opencl_platform(cl_platform_id *my_platform);

int create_opencl_device(cl_platform_id *platforms, cl_device_id *my_devices);

int create_opencl_context(cl_context *context, cl_device_id *my_devices);

int create_opencl_cq(cl_command_queue *commandQueue, cl_context context,
                     cl_device_id *my_devices);

int ffd_opencl_solve(cl_context context, cl_device_id device,
                     cl_command_queue commandQueue);

#endif

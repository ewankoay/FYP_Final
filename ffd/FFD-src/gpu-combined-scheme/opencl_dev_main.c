/// \file   main entrance for opencl development
///
/// \author Wei Tian, Wei.Tian@SE.com
///
///
/// \date   3/14/2019
///
///
///////////////////////////////////////////////////////////////////////////////
#include "opencl_solver.h"
///////////////////////////////////////////////////////////////////////////////
/// Main routine of OPENCL FFD
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////

int main() {
  ffd_log("Start running Fast Fluid Dynamics simulations provided by Schneider "
          "Electric; Report to Wei.Tian@Schneider-Electric.com for problems, "
          "bugs and questions.",
          FFD_NEW);
  cl_platform_id platforms[1];
  cl_device_id devices[1];
  cl_context contexts[1];
  cl_command_queue commandQueues[1];
  cl_device_id device;
  cl_context context;
  cl_command_queue commandQueue;

  // create platform
  if (create_opencl_platform(platforms) != 0) {
    ffd_log("main(): cannot create platform", FFD_ERROR);
    return 1;
  }
  // create device
  if (create_opencl_device(platforms, devices) != 0) {
    ffd_log("main(): cannot create device", FFD_ERROR);
    return 1;
  }
  // create context
  if (create_opencl_context(contexts, devices) != 0) {
    ffd_log("main(): cannot create context", FFD_ERROR);
    return 1;
  }
  // create command queue
  if (create_opencl_cq(commandQueues, contexts[0], devices) != 0) {
    ffd_log("main(): cannot create command queue", FFD_ERROR);
    return 1;
  }
  // choose platform,contexts,device,and command queues
  device = devices[0];
  context = contexts[0];
  commandQueue = commandQueues[0];
  // run FFD simulation
  ffd_opencl_solve(context, device, commandQueue);

  return 0;
}

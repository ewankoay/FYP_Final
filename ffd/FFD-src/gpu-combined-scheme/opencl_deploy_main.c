///////////////////////////////////////////////////////////////////////////////
///
/// \file   opencl_iso_main.c
///
/// \brief  main entrance of parallel FFD program with platform and device
/// number
///  as the input argument for deployment purpose
///
/// \author Wei Tian, Wei.Tian@SE.com
///////////////////////////////////////////////////////////////////////////////

#include "opencl_solver.h"
#define PRINT_OUT 0
///////////////////////////////////////////////////////////////////////////////
/// Main routine of OPENCL FFD
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
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

  // create platform and device
  if (configure_platform_device(argc, argv, devices, contexts, commandQueues) !=
      0) {
    ffd_log("main(): cannot create platform", FFD_ERROR);
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
///////////////////////////////////////////////////////////////////////////////////////
// configure the platform and device
//////////////////////////////////////////////////////////////////////////////////////
int configure_platform_device(int argc, char *argv[], cl_device_id *my_devices,
                              cl_context *my_contexts,
                              cl_command_queue *my_commandQueues) {
  int platform_device[3] = {-1, -1, -1};
  /*Translate and validate the targeted platform and device from the
   * command-line argument*/
  parse_argument(argc, argv, platform_device);

  /***************************CHECK PLATFORM AND
   * DEVICES*****************************************/
  cl_uint numPlatforms;  // the NO. of platforms
  cl_platform_id *platforms = NULL;
  cl_platform_id platform = NULL;  // the chosen platform
  int platform_id = 0;
  sprintf(msg, "main(): start getting platform information");
  ffd_log(msg, FFD_NORMAL);
  // get number of platforms and store it into the platforms variable
  cl_int status = clGetPlatformIDs(0, NULL, &numPlatforms);
  platforms = (cl_platform_id *)malloc(numPlatforms * sizeof(cl_platform_id));
  status = clGetPlatformIDs(numPlatforms, platforms, NULL);
  if (status != CL_SUCCESS) {
    sprintf(msg, "Error: Getting platforms!");
    ffd_log(msg, FFD_ERROR);
    return FAILURE;
  }

  /*For clarity, choose the first available platform. */
  if (PRINT_OUT) {
    printf("NUMBER OF PLATFORMS: %d\n", numPlatforms);
  }

  /*Compare the user specified platform and available platforms */
  int USE_PLATFORM = 0;
  if (platform_device[0] + 1 > numPlatforms) {
    sprintf(msg, "chosen platform is larger than the available platforms");
    ffd_log(msg, FFD_ERROR);
    exit(1);
  }
  else {
    USE_PLATFORM = platform_device[0];
    // print out the selected platform information
    char *prof_info = NULL;
    size_t prof_len;
    status = clGetPlatformInfo(platforms[USE_PLATFORM], CL_PLATFORM_NAME, NULL,
                               NULL, &prof_len);
    prof_info = (char *)malloc(sizeof(char) * prof_len);
    status = clGetPlatformInfo(platforms[USE_PLATFORM], CL_PLATFORM_NAME,
                               prof_len, prof_info, NULL);
    if (PRINT_OUT) {
      printf("%d: %s\n", USE_PLATFORM, prof_info);
    }
    free(prof_info);
    // choose the platform according to the device to be used
    platform = platforms[USE_PLATFORM];
    free(platforms);
  }

  sprintf(msg, "main(): finishing checking avaiblable platforms");
  ffd_log(msg, FFD_NORMAL);

  /*Step 2:Query the platform and choose the first GPU device if has
   * one.Otherwise use the CPU as device.*/
  cl_uint numDevices = 0;
  cl_device_id devices[1] = {NULL};
  cl_device_id *devices_cpu = NULL;
  cl_device_id *devices_gpu = NULL;

  // Determine which platform to be used
  int USE_DEVICE = 0;
  if (platform_device[1] == -1 && platform_device[2] == -1) {
    sprintf(msg, "selected neighter GPU or CPU, not supported");
    ffd_log(msg, FFD_ERROR);
    exit(1);
  }
  else if (platform_device[1] != -1 && platform_device[2] != -1) {
    sprintf(msg, "selected both GPU and CPU, not supported");
    ffd_log(msg, FFD_ERROR);
    exit(1);
  }
  else if (platform_device[1] == -1) {
    USE_DEVICE = 0;  // No CPU is specified and use GPU
  }
  else if (platform_device[2] == -1) {
    USE_DEVICE = 1;  // No GPU is specified and use CPU
  }
  else {
    sprintf(msg, "selected device is not available");
    ffd_log(msg, FFD_ERROR);
    exit(1);
  }

  //=========================================//
  char *device_info = NULL;
  size_t device_len;
  int device_id = 0;
  int gpu_id = 0;
  if (!USE_DEVICE)  // GPU available.
  {
    // get the number of GPU device
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
    // check if selected device is available or not
    if (platform_device[2] + 1 > numDevices) {
      sprintf(msg, "chosen GPU is out of the available platforms");
      ffd_log(msg, FFD_ERROR);
      exit(1);
    }
    else {
      device_id = platform_device[2];
      // allocation memory for device
      devices_gpu = (cl_device_id *)malloc(numDevices * sizeof(cl_device_id));
      // get values for device
      status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices,
                              devices_gpu, NULL);
      // output device information
      status = clGetDeviceInfo(devices_gpu[device_id], CL_DEVICE_NAME, NULL,
                               NULL, &device_len);
      device_info = (char *)malloc(sizeof(char) * device_len);
      status = clGetDeviceInfo(devices_gpu[device_id], CL_DEVICE_NAME,
                               device_len, device_info, NULL);
      if (PRINT_OUT) {
        printf("%s\n", device_info);
      }
      free(device_info);
      // choose device
      devices[0] = devices_gpu[device_id];
    }
  }
  else {
    // get the number of CPU device
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &numDevices);
    // check if selected device is available or not
    if (platform_device[1] + 1 > numDevices) {
      sprintf(msg, "chosen GPU is out of the available platforms");
      ffd_log(msg, FFD_ERROR);
      exit(1);
    }
    else {
      device_id = platform_device[1];
      // allocation memory for device
      devices_cpu = (cl_device_id *)malloc(numDevices * sizeof(cl_device_id));
      // get values for device
      status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, numDevices,
                              devices_cpu, NULL);
      // print out CPU information
      status = clGetDeviceInfo(devices_cpu[device_id], CL_DEVICE_NAME, NULL,
                               NULL, &device_len);
      device_info = (char *)malloc(sizeof(char) * device_len);
      status = clGetDeviceInfo(devices_cpu[device_id], CL_DEVICE_NAME,
                               device_len, device_info, NULL);
      if (PRINT_OUT) {
        printf("%s\n", device_info);
      }
      free(device_info);
      // choose CPU device
      devices[0] = devices_cpu[device_id];
    }
  }

  sprintf(msg, "main(): finishing chosing the avaiblable platform");
  ffd_log(msg, FFD_NORMAL);

  // create context
  cl_context context = clCreateContext(NULL, 1, devices, NULL, NULL, NULL);
  sprintf(msg, "main(): finishing creating context for the platform");
  ffd_log(msg, FFD_NORMAL);

  // create command queue
  cl_command_queue commandQueue =
      clCreateCommandQueue(context, devices[0], NULL, NULL);
  sprintf(msg, "main(): finishing creating queue associated with the context");
  ffd_log(msg, FFD_NORMAL);

  // store the device, context, and command queues
  my_devices[0] = devices[0];
  my_contexts[0] = context;
  my_commandQueues[0] = commandQueue;

  return 0;
}

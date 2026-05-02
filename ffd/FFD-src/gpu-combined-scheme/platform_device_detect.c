///////////////////////////////////////////////////////////////////////////////
///
/// \file   platform_device_detect.c
///
/// \brief  a small program to detect available platform and device
///
/// \author Wei Tian
///         Schneider Electric
///         Wei.Tian@Schneider-Electric.com
///
/// \date   4/20/2018
///
///////////////////////////////////////////////////////////////////////////////

#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#ifdef __APPLE__ || __linux__
  #include <OpenCL/opencl.h>
#else
  #include <CL/cl.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#define WRITE_TO_JSON 1 // write the platform and device info in JSON to stdout
static char msg[100000]; // define the char vairable to store the info in JSON

#ifdef _MSC_VER
  #define _CRT_SECURE_NO_WARNINGS
#endif

int main () {
  //Step1: Getting platforms and choose an available one
  // variables for platform
  cl_uint numPlatforms;  //the NO. of platforms
  cl_platform_id* platforms = NULL;
  cl_platform_id platform = NULL;  //the chosen platform
  int platform_id = 0;
  // variabl;es for device
  cl_device_id   devices[1] = { NULL };
  cl_device_id   *devices_cpu = NULL;
  cl_device_id   *devices_gpu = NULL;
  cl_uint numGPU = 0;
  cl_uint numCPU = 0;
  int device_id = 0;
  // get number of platforms and store it into the platforms variable
  cl_int  status = clGetPlatformIDs(0, NULL, &numPlatforms);
  platforms = (cl_platform_id*)malloc(numPlatforms * sizeof(cl_platform_id));
  status = clGetPlatformIDs(numPlatforms, platforms, NULL);
  if (status != CL_SUCCESS){
    printf("Error: Getting platforms!\n");
    return 1;
  }
  // print out the number of available platoforms
  if(!WRITE_TO_JSON) printf("NUMBER OF PLATFORMS: %d\n", numPlatforms);
  // If there is at least one platform avaible, print out the available platform and associated devices
  if (numPlatforms > 0){
    // Loop all the platform
    for (platform_id = 0; platform_id<numPlatforms; platform_id++){
      // print out platform info
      char *prof_info= NULL;
      size_t prof_len;
      status = clGetPlatformInfo(platforms[platform_id], CL_PLATFORM_NAME, NULL,NULL, &prof_len);
      prof_info = (char*)malloc(sizeof(char)*prof_len);
      status = clGetPlatformInfo(platforms[platform_id], CL_PLATFORM_NAME, prof_len, prof_info, NULL);
      if(!WRITE_TO_JSON) printf("%d: %s\n", platform_id, prof_info);
      // write the platform information to JSON format
      if (platform_id == 0) {
        printf ("{\"%d\":{\"label\":\"%s\",\n",platform_id,prof_info);
      }
      else {
        printf ("\"%d\":{\"label\":\"%s\",\n",platform_id,prof_info);
      }
      free(prof_info);
      // get the number of GPU device
      status = clGetDeviceIDs(platforms[platform_id], CL_DEVICE_TYPE_GPU, 0, NULL, &numGPU);
      // write number of GPU in JSON formats
      printf ("\"num_gpus\":%d,\n",numGPU);
      // write empty vector if number of gpu is 0
      if (numGPU==0) {
        printf("\"gpus\":[],\n");
      }
      // get the number of CPU device
      status = clGetDeviceIDs(platforms[platform_id], CL_DEVICE_TYPE_CPU, 0, NULL, &numCPU);
      // allocation memory for GPU device
      if (numGPU > 0) devices_gpu = (cl_device_id*)malloc(numGPU * sizeof(cl_device_id));
      // get values for GPU device
      status = clGetDeviceIDs(platforms[platform_id], CL_DEVICE_TYPE_GPU, numGPU, devices_gpu, NULL);
      // allocation memory for CPU device
      if (numCPU > 0) devices_cpu = (cl_device_id*)malloc(numCPU * sizeof(cl_device_id));
      // get values for CPU device
      status = clGetDeviceIDs(platforms[platform_id], CL_DEVICE_TYPE_CPU, numCPU, devices_cpu, NULL);
      // output GPU device information
      if(!WRITE_TO_JSON) printf("\tNUMBER OF GPUs: %d\n", numGPU);
      for (device_id = 0; device_id < numGPU; device_id++) {
        char *gpu_device_info = NULL;
        size_t gpu_device_len;
        status = clGetDeviceInfo(devices_gpu[device_id], CL_DEVICE_NAME, NULL, NULL, &gpu_device_len);
        gpu_device_info = (char*)malloc(sizeof(char)*gpu_device_len);
        status = clGetDeviceInfo(devices_gpu[device_id], CL_DEVICE_NAME, gpu_device_len, gpu_device_info, NULL);
        if(!WRITE_TO_JSON) printf("\t\t%d: %s\n", device_id, gpu_device_info);
        // write name of GPU in JSON formats
        if (device_id == 0 && numGPU>1){
          printf ("\"gpus\":[{\"dev_index\":%d,\"label\":\"%s\"},\n",device_id,gpu_device_info);
        }
        else if (device_id == 0 && numGPU==1){
          printf ("\"gpus\":[{\"dev_index\":%d,\"label\":\"%s\"}],\n",device_id,gpu_device_info);
        }
        else if (device_id == numGPU-1){
          printf ("{\"dev_index\":%d,\"label\":\"%s\"}],\n",device_id,gpu_device_info);
        }
        else {
          printf ("{\"dev_index\":%d,\"label\":\"%s\"}",device_id,gpu_device_info);
        }
      }//end of for loop
      // output GPU device information
      if(!WRITE_TO_JSON) printf("\tNUMBER OF CPUs: %d\n", numCPU);
      // write number of CPU in JSON formats
      printf ("\"num_cpus\":%d,\n",numCPU);
      // write empty vector if number of cpu is 0
      if (numCPU==0) {
        printf("\"cpus\":[]");
      }
      for (device_id = 0; device_id < numCPU; device_id++) {
        char *cpu_device_info = NULL;
        size_t cpu_device_len;
        status = clGetDeviceInfo(devices_cpu[device_id], CL_DEVICE_NAME, NULL, NULL, &cpu_device_len);
        cpu_device_info = (char*)malloc(sizeof(char)*cpu_device_len);
        status = clGetDeviceInfo(devices_cpu[device_id], CL_DEVICE_NAME, cpu_device_len, cpu_device_info, NULL);
        if(!WRITE_TO_JSON) printf("\t\t%d: %s\n", device_id, cpu_device_info);
        // write name of GPU in JSON formats
        if (device_id == 0 && numCPU>1){
          printf ("\"cpus\":[{\"dev_index\":%d,\"label\":\"%s\"},\n",device_id,cpu_device_info);
        }
        else if (device_id == 0 && numCPU==1){
          printf ("\"cpus\":[{\"dev_index\":%d,\"label\":\"%s\"}]",device_id,cpu_device_info);
        }
        else if (device_id == numCPU-1){
          printf ("{\"dev_index\":%d,\"label\":\"%s\"}]",device_id,cpu_device_info);
        }
        else {
          printf ("{\"dev_index\":%d,\"label\":\"\"%s\"},\n",device_id,cpu_device_info);
        }
      }//end of for loop
      // write the platform information (end) to JSON format
      if (platform_id != numPlatforms-1) {
        printf ("},\n");
      }
      else {
        printf ( "}\n");
      }
    } //end of for loop

  } //end of if loop
  printf ( "}\n"); //end of JSON
  // write the JSON format to a file
  return 0;
}

/*
 * OCLManager.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: pfandedd
 */

#include <iostream>
#include <sstream>

#include "OCLManager.hpp"

#include <sgpp/base/exception/operation_exception.hpp>

namespace SGPP {
namespace base {

OCLManager::OCLManager(base::OpenCLConfigurationParameters parameters) :
    parameters(parameters) {
  verbose = parameters.getAsBoolean("OCL_MANAGER_VERBOSE");
  cl_int err;

  //upper limit for number of devices of a single platform

  // determine number of available OpenCL platforms
  err = clGetPlatformIDs(0, nullptr, &num_platforms);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get number of OpenCL platforms. Error Code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  if (verbose) {
    std::cout << "OCL Info: " << num_platforms << " OpenCL Platforms have been found" << std::endl;
  }

  // get available platforms
  platform_ids = new cl_platform_id[num_platforms];
  err = clGetPlatformIDs(num_platforms, platform_ids, nullptr);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get Platform ID. Error Code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  for (cl_uint ui = 0; ui < num_platforms; ui++) {
    char vendor_name[128] = { 0 };
    err = clGetPlatformInfo(platform_ids[ui], CL_PLATFORM_VENDOR, 128 * sizeof(char), vendor_name, nullptr);
    if (CL_SUCCESS != err) {
      std::stringstream errorString;
      errorString << "OCL Error: Can't get platform vendor!" << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    } else {
      if (vendor_name != nullptr && verbose) {
        std::cout << "OCL Info: Platform " << ui << " vendor name: " << vendor_name << std::endl;
      }
    }

    char platform_name[128] = { 0 };
    err = clGetPlatformInfo(platform_ids[ui], CL_PLATFORM_NAME, 128 * sizeof(char), platform_name, nullptr);
    if (CL_SUCCESS != err) {
      std::stringstream errorString;
      errorString << "OCL Error: Can't get platform name!" << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    } else {
      if (platform_name != nullptr && verbose) {
        std::cout << "OCL Info: Platform " << ui << " name: " << platform_name << std::endl;
      }

      if (parameters["PLATFORM"].compare(platform_name) == 0) {
        platform_id = platform_ids[ui];
        if (verbose) {
          std::cout << "platform selected" << std::endl;
        }
      }
    }
  }

  if (verbose) {
    std::cout << std::endl;
  }

  if (parameters["PLATFORM"].compare("first") == 0) {
    if (verbose) {
      std::cout << "using first platform" << std::endl;
    }
    platform_id = platform_ids[0];
  }

  // Find out how many devices there are
  if (parameters["DEVICE_TYPE"] == "CL_DEVICE_TYPE_CPU") {
    if (verbose) {
      std::cout << "OCL Info: looking for CPU device" << std::endl;
    }
    // get the number of devices
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 0, nullptr, &num_devices);
    device_ids = new cl_device_id[num_devices];
    // get the device ids
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, num_devices, device_ids, nullptr);
  } else if (parameters["DEVICE_TYPE"] == "CL_DEVICE_TYPE_GPU") {
    if (verbose) {
      std::cout << "OCL Info: looking for GPU device" << std::endl;
    }
    // get the number of devices
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 0, nullptr, &num_devices);
    device_ids = new cl_device_id[num_devices];
    // get the device ids
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, num_devices, device_ids, nullptr);
  } else if (parameters["DEVICE_TYPE"] == "CL_DEVICE_TYPE_ACCELERATOR") {
    if (verbose) {
      std::cout << "OCL Info: looking for device of accelerator type" << std::endl;
    }
    // get the number of devices
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, 0, nullptr, &num_devices);
    device_ids = new cl_device_id[num_devices];
    // get the device ids
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, num_devices, device_ids, nullptr);
  } else if (parameters["DEVICE_TYPE"] == "CL_DEVICE_TYPE_ALL") {
    if (verbose) {
      std::cout << "OCL Info: looking for device of all available devices" << std::endl;
    }
    // get the number of devices
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, 0, nullptr, &num_devices);
    device_ids = new cl_device_id[num_devices];
    // get the device ids
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, num_devices, device_ids, nullptr);
  } else {
    throw SGPP::base::operation_exception(
        "OCL Error: No device found or unknown type specified (supported are: CPU, GPU, accelerator and all)");
  }

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Unable to get Device IDs. Error Code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  if (verbose) {
    std::cout << "OCL Info: " << num_devices << " OpenCL devices have been found!" << std::endl;
  }

  cl_uint maxDevices = (cl_uint) parameters.getAsUnsigned("MAX_DEVICES");

  if (parameters["SELECT_SPECIFIC_DEVICE"].compare("DISABLED") != 0) {
    if (maxDevices != 1) {
      std::stringstream errorString;
      errorString
          << "OCL Error: Cannot select a specific device if more than one device is used, MAX_DEVICES be set incorrectly"
          << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
    size_t selectedDevice = parameters.getAsUnsigned("SELECT_SPECIFIC_DEVICE");
    if (selectedDevice > num_devices) {
      std::stringstream errorString;
      errorString << "OCL Error: Illegal value set for \"SELECT_SPECIFIC_DEVICE\"" << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
    device_ids[0] = device_ids[selectedDevice];

    if (verbose) {
      std::cout << "OCL Info: select device number " << selectedDevice << std::endl;
    }
  }

  if (maxDevices != 0 && maxDevices < num_devices) {
    num_devices = maxDevices;
  }

  if (verbose) {
    std::cout << "OCL Info: using " << num_devices << " device/s" << std::endl;
  }

  // allocate arrays
  command_queue = new cl_command_queue[num_devices];

  // Create OpenCL context
  context = clCreateContext(0, num_devices, device_ids, nullptr, nullptr, &err);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to create OpenCL context! Error Code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  // Creating the command queues
  for (size_t i = 0; i < num_devices; i++) {
    char buffer[128];
    err = clGetDeviceInfo(device_ids[i], CL_DEVICE_NAME, 128 * sizeof(char), &buffer, nullptr);
    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to read the device name for device: " << i << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }

    if (verbose) {
      std::cout << "OCL Info: device name: " << buffer << std::endl;
    }

    command_queue[i] = clCreateCommandQueue(context, device_ids[i],
    CL_QUEUE_PROFILING_ENABLE, &err);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create command queue! Error Code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
  }

  if (verbose) {
    std::cout << "OCL Info: Successfully initialized OpenCL (local workgroup size: "
        << parameters.getAsUnsigned("LOCAL_SIZE") << ")" << std::endl << std::endl;
  }
}

/**
 * @brief buildKernel builds the program that is represented by @a program_src and creates @a num_devices kernel objects
 * that are stored into the array @a kernel (must be already allocated with at least @a num_devices )
 *
 * @param program_src the source of the program to compile
 * @param kernel_name name of the kernel function (in program_src) to create the kernel for
 * @param context OpenCL context
 * @param num_devices number of OpenCL devices
 * @param device_ids array with device ids, necessary for displaying build info
 * @param kernel already allocated array: the resulting kernels are put into this array, one for each device (=> at least num_devices entries)
 * @return
 */
void OCLManager::buildKernel(const std::string& program_src, const char* kernel_name, cl_context context,
    size_t num_devices, cl_device_id* device_ids, cl_kernel* kernel) {
  cl_int err;

  // setting the program
  const char* kernel_src = program_src.c_str();
  cl_program program = clCreateProgramWithSource(context, 1, &kernel_src,
  NULL, &err);

  if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: Failed to create program! Error Code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
  }

  std::string build_opts;
  if (parameters.getAsBoolean("ENABLE_OPTIMIZATIONS")) {
    //TODO: user should be able to change
    build_opts = parameters["OPTIMIZATION_FLAGS"]; // -O5  -cl-mad-enable -cl-denorms-are-zero -cl-no-signed-zeros -cl-unsafe-math-optimizations -cl-finite-math-only -cl-fast-relaxed-math
  } else {
    build_opts = "-cl-opt-disable"; // -g
  }

  //TODO: check multi device support
  // compiling the program
  err = clBuildProgram(program, 0, NULL, build_opts.c_str(), NULL, NULL);

  //if (err != CL_SUCCESS) {


  if (parameters.getAsBoolean("SHOW_BUILD_LOG")) {
    size_t len;
    char buffer[4096];

    // get the build log
    clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
    if (verbose) {
      std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
    }
  }

    if (err != CL_SUCCESS) {
    std::stringstream errorString;
    errorString << "OCL Error: OpenCL build error. Error code: " << err << std::endl;
    throw SGPP::base::operation_exception(errorString.str());
    }



  // creating the kernel
  for (size_t i = 0; i < num_devices; i++) {
    kernel[i] = clCreateKernel(program, kernel_name, &err);

    if (err != CL_SUCCESS) {
      std::stringstream errorString;
      errorString << "OCL Error: Failed to create kernel! Error code: " << err << std::endl;
      throw SGPP::base::operation_exception(errorString.str());
    }
  }

  if (program) {
    clReleaseProgram(program);
  }
}

}
}
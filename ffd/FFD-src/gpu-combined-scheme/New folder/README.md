# Compiling on Ubuntu 18.04 LTS

`OpenCL` headers and library are the prerequisite. Before compiling the codes by hitting `make FFD_OPENCL_LINUX`, one should make sure everything is all set.

## Install dependencies

`sudo apt updatesudo`
`apt install ocl-icd-opencl-dev`

### Debugging and further information

`which gcc` to figure out the version of gcc.

`echo | gcc -E -Wp,-v -` and find if `CL/cl.h` is in the searching path. See [this post ](https://stackoverflow.com/questions/17939930/finding-out-what-the-gcc-include-path-is) for detail.

`ldconfig -v 2>/dev/null | grep -v ^$'\t'` for dependent library search path. See [this post](https://stackoverflow.com/questions/9922949/how-to-print-the-ldlinker-search-path) for detail.

## Build

`build.sh` will build both the OPEN_CL device detection application as well as
the FFD engine application


## Other debugging issues
`cat /proc/driver/nvidia/version`

```
NVRM version: NVIDIA UNIX x86_64 Kernel Module  390.48  Thu Mar 22 00:42:57 PDT 2018
GCC version:  gcc version 7.3.0 (Ubuntu 7.3.0-16ubuntu3)
```

Run `nvidia-smi`
If that fails, check
`dmesg`

```
NVRM: API mismatch: the client has the version 390.116, but
                 NVRM: this kernel module has the version 390.48.  Please
                 NVRM: make sure that this kernel module and all NVIDIA driver
                 NVRM: components have the same version.
```
A mismatch like this can likely be solved with a reboot

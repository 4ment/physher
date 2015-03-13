# Binaries #

Binaries are only available for intel-based Mac [physher-mac.tar.gz](https://drive.google.com/folderview?id=0B1fZePvFRVeBemNjUE4zZWdPVVU&usp=sharing). These were compiled with Mac OSX 10.9 with backward compatibility for 10.6 and above.

_physher_ is available in 3 flavours: SSE, Pthread, and SSE+Pthread. For bootstrapping, discrete and local clock algorithms it is strongly advised to use one of the 2 Pthread versions.

# Compile source code #

## Source ##

Download [physher-src.tar.gz](https://drive.google.com/folderview?id=0B1fZePvFRVeBemNjUE4zZWdPVVU&usp=sharing)

## Requirements ##

  * C compiler: only [GCC](http://gcc.gnu.org) and Clang are supported
  * [CMake](http://www.cmake.org/cmake/resources/software.html)

## Optional ##
  * [OpenMP](http://openmp.org) library: for multithreading or
  * Pthreads library: for multithreading.

As far as I know GCC supports OpenMP since version 4.4. The current version of LLVM/clang does not support OpenMP but a newer version can be downloaded (see below).

## Note for Mac users ##

There are at least 3 ways to get a compiler on a Mac:
  * GCC is usually installed with Xcode. Unfortunately, the default compiler, LLVM (clang) compiler version 5.0 does not support OpenMP. However,  with this compiler one can compile a newer version of GCC or clang that supports OpenMP but other libraries need to be compiled before compiling it.

  * Package managers like [HomeBrew](http://brew.sh) and [MacPorts](http://www.macports.org) can be used to install GCC and probably other compilers.

  * Download a precompiled compiler from a website like [HPC Mac OSX](http://hpc.sourceforge.net). That's probably the easiest way.

## Installation ##

Once CMake and a compiler are installed download [physher-src.tar.gz](https://drive.google.com/folderview?id=0B1fZePvFRVeBemNjUE4zZWdPVVU&usp=sharing)

  * Linux and Mac:
```
tar -zxvf physher-src.tar.gz
cd physher-src
mkdir Release
cd Release
cmake ..
make -j5
sudo make install
```

The last step is optional. It simply installs _physher_, _simultron_, _bootstrap_... in /usr/local/bin

_cmake_ can take arguments in order to enable vectorisation and parallelisation. Some computers are not able to use SSE, Pthreads, and/or OpenMP

To use SSE and OpenMP together (RECOMMENDED)
```
cmake -DUSE_OPENMP_SUPPORT=1 -DUSE_SSE_SUPPORT=1 ..
```


To use SSE and Pthreads together (RECOMMENDED IF USING DEFAULT COMPILER ON MAC)
```
cmake -DUSE_PTHREAD_SUPPORT=1 -DUSE_SSE_SUPPORT=1 ..
```

To use a different compiler use something like
```
cmake -DCMAKE_C_COMPILER=gcc4.8 ..
```

IMPORTANT NOTE: Physher is written in C and for some reason _cmake_ can complain about not finding OpenMP even though it is present. Before giving up on OpenMP uncomment ENABLE\_LANGUAGE(CXX) in both CMakeLists.txt files in the source and phyc directories.


Alternatively there are make files in the source directory.
```
tar -zxvf physher-src.tar.gz
cd physher-src
make -f makefile.sse
```



  * Windows
I have not tried to compile on Windows but it should be doable with Cygwin.
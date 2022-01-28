# IPMPI - Improved Profiler for MPI
## How to use
***
To use the profiler, first of all change the compiler in Makefile to either mpic++ (for MPICH) or mpiicpc (for Intel MPI) and then simply execute the below instructions:

```bash
$ bash setup.sh
$ source $HOME/.ipmpi
```

_.ipmpi_ - will be present under your home directory

_.setup.sh_ - will be present under "ipmpi" directory and is used for the initial setup.

**Note:** After sourcing the file, you can simply execute your programs using MPI and can get profiled data by default under _$HOME/ipmpi_profiles_ or can change the profiling directory by exporting like this:
```bash
export PROFILE_DIR=path/to/profiling/directory
```

## About the code:
***
1. Intercepting methods for MPI collectives are defined in _ipmpi.cpp_ which uses a number of algorithms for each collective which are defined in _helper.cpp_ file.
2. The method declarations are present inside _include/helper.h_.
3. Intercepting methods include - MPI_Bcast, MPI_Scatter, MPI_Gather, MPI_Alltoall, MPI_Allgather, MPI_Allreduce, MPI_Sendrecv, MPI_Reduce, MPI_Allreduce.


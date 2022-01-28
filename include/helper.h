#ifndef _helper_h
#define _helper_h
// Used to define the datatype for comm_matrix >> not used currently
#define MAT_TYPE long
// Include all the functions definitions here for making the code reusable and readable
#define MPIR_CVAR_ALLGATHER_LONG_MSG_SIZE 524288
#define MPIR_CVAR_ALLGATHER_SHORT_MSG_SIZE 81920
#define MPIR_CVAR_ALLTOALL_SHORT_MSG_SIZE 256
#define MPIR_CVAR_ALLTOALL_MEDIUM_MSG_SIZE 32768
#define MPIR_CVAR_ALLTOALL_THROTTLE 32

#define MPIR_CVAR_BCAST_SHORT_MSG_SIZE 12288
#define MPIR_CVAR_BCAST_MIN_PROCS 8
#define MPIR_CVAR_BCAST_LONG_MSG_SIZE 524288
#define MPIR_CVAR_GATHER_VSMALL_MSG_SIZE 1024
#define MPIR_CVAR_REDUCE_SHORT_MSG_SIZE 2048
// These next 2 definition is relevant for comm_matrix_count only
#define BIN_RESIZE_FACTOR 2
#define BIN_RESIZE_VALUE 0

// These are essential for MPI_Allreduce
#define MPIR_CVAR_ALLREDUCE_SHORT_MSG_SIZE 2048
#define HANDLE_GET_KIND(a) (((unsigned)(a)&HANDLE_KIND_MASK)>>HANDLE_KIND_SHIFT)
#define HANDLE_KIND_BUILTIN  0x1
/* Mask assumes that ints are at least 4 bytes */
#define HANDLE_KIND_MASK 0xc0000000
#define HANDLE_KIND_SHIFT 30

#include <vector>
#include <string>
#include <cstring>
#include <iostream>

#include <unistd.h>
#include <pwd.h>

#include "mpi.h"
#include <cmath>
#include<unordered_map>
#include <cstdlib>

/***********************Variable Declaration******************************/

extern int comm_size;
extern std::vector<std::vector<long>> comm_matrix_sum, comm_matrix_max, comm_matrix_min, p2p_mat;
// extern FILE *logger;
extern std::vector<std::vector<std::vector<int>>> comm_matrix_count;

extern std::unordered_map<std::string, int> log_calls;

// For computing #sends per collective
extern std::unordered_map<std::string, std::unordered_map<int, std::unordered_map<int, int>>> sends;
extern std::string parent_func;
extern int rank;
extern long BIN_SIZE;

/***********************Methods Declaration******************************/

void reset_comm_matrix(std::vector<std::vector<long>> &comm_matrix, int value);
void reset_count_matrix(std::vector<std::vector<std::vector<int>>> &count_matrix);

// for message count matrix
inline long getBinNumber(long number) {
  return number / BIN_SIZE;
}

inline void updateCommMatrix(long msg_size, int src, int dst) {
  
  /*****************Updating the comm_matrix******************/
  sends[parent_func][src][dst]++;
  
  if(comm_matrix_sum[src][dst] == -1)
    comm_matrix_sum[src][dst]++;
  if(comm_matrix_sum[dst][src] == -1)
    comm_matrix_sum[dst][src]++;
  
  // Updating the sum matrix
  comm_matrix_sum[src][dst] += msg_size;
  comm_matrix_sum[dst][src] += msg_size;
  
  // Updating the max matrix
  comm_matrix_max[src][dst] = std::max(comm_matrix_max[src][dst], msg_size);
  comm_matrix_max[dst][src] = std::max(comm_matrix_max[src][dst], msg_size);

  // Updating the min matrix
  comm_matrix_min[src][dst] = (comm_matrix_min[src][dst] != -1) ? std::min(comm_matrix_min[src][dst], msg_size) : msg_size;
  comm_matrix_min[dst][src] = (comm_matrix_min[dst][src] != -1) ? std::min(comm_matrix_min[dst][src], msg_size) : msg_size;
 
  // Updating the count matrix
  long bin = getBinNumber(msg_size);
  long sender_size = comm_matrix_count[src][dst].size();
  long receiver_size = comm_matrix_count[dst][src].size();
  if(bin >= sender_size)
    comm_matrix_count[src][dst].resize(BIN_RESIZE_FACTOR*bin, BIN_RESIZE_VALUE);
  if(bin >= receiver_size)
    comm_matrix_count[dst][src].resize(BIN_RESIZE_FACTOR*bin, BIN_RESIZE_VALUE);
  comm_matrix_count[src][dst][bin]++;
  comm_matrix_count[dst][src][bin]++;
}

void populateCommMatrix(std::vector<std::vector<long>> &comm_matrix, std::string& s, std::string msg);
void populateCountMatrix(std::vector<std::vector<std::vector<int>>> &count_matrix, std::string& s, std::string msg);
void logging(std::string str, int count, int type_size);
void logging_dst(std::string str, int count, int type_size, std::string type, int dst);
std::string getFunc(char *temp);

/***************************Profiling Methods****************************/

/*---------Bcast-----------*/
void profile_Bcast_intra_binomial(int root, int type_size, int count);
void profile_Scatter_for_bcast(int root, int nbytes);
void profile_Bcast_intra_scatter_recursive_doubling_allgather(int root, int type_size, int count);
void profile_Bcast_intra_scatter_ring_allgather(int root, int type_size, int count);

/*---------Allgather-----------*/
void profile_Allgather_intra_recursive_doubling(int &recvcount, int &type_size);
void profile_Allgather_intra_brucks(int &recvcount, int &type_size);
void profile_Allgather_intra_ring(int &recvcount, int &type_size);

/*------------Allreduce------------*/
void profile_Allreduce_intra_recursive_doubling(int count, int type_size);
void profile_Allreduce_intra_reduce_scatter_allgather(int count, int type_size);

/***************Scatter*****************/
void profile_Scatter_binomial(int root, int sendtype_size, int sendcount, int recvtype_size, int recvcount);

/***************Gather*****************/
void profile_Gather_binomial(int root, MPI_Datatype sendtype, int sendcount, MPI_Datatype recvtype, int recvcount);

/***************Reduce*****************/
void profile_Reduce_intra_reduce_scatter_gather(MPI_Datatype datatype, int root, int count);
void profile_Reduce_intra_binomial(MPI_Datatype datatype, int root, int count, MPI_Op op);

/***************Alltoall****************/
void profile_Alltoall_intra_pairwise_sendrecv_replace(MPI_Datatype sendtype, MPI_Datatype recvtype, int sendcount, int recvcount);
void profile_Alltoall_intra_brucks(MPI_Datatype sendtype, MPI_Datatype recvtype, int sendcount, int recvcount);
void profile_Alltoall_intra_scattered(MPI_Datatype sendtype, MPI_Datatype recvtype, int sendcount, int recvcount);
void profile_Alltoall_intra_pairwise(MPI_Datatype sendtype, MPI_Datatype recvtype, int sendcount, int recvcount);

/**************Barrier***********/
void profile_Barrier();
#endif

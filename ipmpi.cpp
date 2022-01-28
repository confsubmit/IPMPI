#include <stdio.h>
#include <execinfo.h>

#include <limits.h>  
#include "include/helper.h"

// Used for logging
#define LOG_MPI_CALLS

int rank = -1;
const char *homedir;
static FILE *logger;
enum IpmpiEvent { Error, Compute, Recv, Send, Bcast, Scatter, Gather };
// const char *homedir;

int comm_size = 0;
long BIN_SIZE=102400;
std::vector<std::vector<long>> comm_matrix_sum, comm_matrix_max, comm_matrix_min, p2p_mat;
std::vector<std::vector<std::vector<int>>> comm_matrix_count;

std::unordered_map<std::string, int> log_calls;
std::unordered_map<std::string, int> log_calls_w_cs;

// For computing #sends per collective
std::unordered_map<std::string, std::unordered_map<int, std::unordered_map<int, int>>> sends;

std::string parent_func="";
// static std::vector<std::vector<vector<long>>> comm_matrix_sum;
std::string getFunc(char *temp) {
  std::string s = temp;
  std::string res="";
  int i;
  bool flag = false;
  for(i=0;i<s.size();i++) {
    if(flag and s[i] != '+')
      res += s[i]; 
    if(s[i] == '(')
      flag = true;
    if(s[i] == '+')
      flag = false;
  }
  return res;
}
void logging(std::string str, int count, int type_size) {
  void *array[10];
  char **strings;
  int size;
  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);
  if (strings != NULL)
  {
    #ifdef DEBUG
      std::string func_name1 = getFunc(strings[1]);
      std::string func_name2 = getFunc(strings[2]);
      log_calls_w_cs[func_name1 + " " + func_name2]++;
      fprintf(logger,"%s @ %s %s %d %s %d\n", str.c_str(), strings[2], " Count: ",count, " TypeSize: ", type_size);  
    #endif
  }
  free (strings);
}
void logging_dst(std::string str, int count, int type_size, std::string type, int dst) {
void *array[10];
  char **strings;
  int size;

  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);
  if (strings != NULL)
  {
#ifdef DEBUG
  std::string func_name1 = getFunc(strings[1]);
  std::string func_name2 = getFunc(strings[2]);
  log_calls_w_cs[func_name1 + " " + func_name2]++;
  fprintf(logger,"%s @ %s %s %d %s %d %s %d\n", str.c_str(), strings[2]," Count: ",count, " TypeSize: ", type_size, type.c_str(), dst);
#endif
  }
  free (strings);
}
long initBinSize(char * ptr) {
    long defaultBinSize=102400, binSize=0;
    if(ptr == NULL)
        return defaultBinSize;
    int digit=0;
    for(int idx=0; idx<std::strlen(ptr);idx++) {
        if(std::isdigit(ptr[idx])) {
            digit=ptr[idx]-'0';
            binSize = 10*binSize+digit;
        }
        else
            return defaultBinSize;
    }
    // std::cout << "BIN_SIZE: "<<binSize << "\n";
    return binSize > 0 ? binSize : defaultBinSize;
}
int MPI_Init(int *argc, char ***argv) {
  parent_func = __FUNCTION__;
  // std::cout << "In Init\n";
#ifdef LOG_MPI_CALLS
  log_calls[__FUNCTION__]++;
#endif
  // Default BIN_SIZE for count matrix
  BIN_SIZE = initBinSize(getenv("BIN_SIZE"));
  std::cout<<"BIN_SIZE: "<<BIN_SIZE<<"\n";
  /********Get HOME directory*******/
  struct passwd *pw = getpwuid(getuid());
  homedir = pw->pw_dir;
  // Get the PROFILE_DIR to store the profile files
  char *PROFILE_DIR = getenv("PROFILE_DIR");

  int result = PMPI_Init(argc, argv);
  PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
  PMPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  char logName[100];

  if (PROFILE_DIR != NULL) {
    // std::cout<<"PROFILE_DIR: "<<PROFILE_DIR<<"\n";
    sprintf(logName, "%s/ipmpi-%d.log", PROFILE_DIR, rank);
  } else {
    sprintf(logName, "%s/ipmpi_profiles/ipmpi-%d.log", homedir, rank);
  }
  logger= fopen(logName, "w+");
  // std::cout << "In Init\n";

  // Reset the communication matrix with -1 implying no communication initially
  int val = -1;
  reset_comm_matrix(comm_matrix_sum, val);
  reset_comm_matrix(comm_matrix_max, val);
  reset_comm_matrix(comm_matrix_min, val);
  reset_comm_matrix(p2p_mat, val);
  reset_count_matrix(comm_matrix_count);
  // std::cout << "In end of MPI Init\n";
  return result;
}
int MPI_Init_thread(int *argc, char ***argv, int required, int *provided) {
  parent_func = __FUNCTION__;
  // std::cout << "In Init\n";
#ifdef LOG_MPI_CALLS
  log_calls[__FUNCTION__]++;
#endif
#ifdef START_FUNC
  std::cout <<"In Init thread\n";
#endif
  // Default BIN_SIZE for count matrix
  BIN_SIZE = initBinSize(getenv("BIN_SIZE"));
  std::cout<<"BIN_SIZE: "<<BIN_SIZE<<"\n";
  /********Get HOME directory*******/
  struct passwd *pw = getpwuid(getuid());
  homedir = pw->pw_dir;
  // Get the PROFILE_DIR to store the profile files
  char *PROFILE_DIR = getenv("PROFILE_DIR");
  
  int result = PMPI_Init_thread(argc, argv, required, provided);
  PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
  PMPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  char logName[100];
  
  if (PROFILE_DIR != NULL) {
    // std::cout<<"PROFILE_DIR: "<<PROFILE_DIR<<"\n";
    sprintf(logName, "%s/ipmpi-%d.log", PROFILE_DIR, rank);
  } else {
    sprintf(logName, "%s/ipmpi_profiles/ipmpi-%d.log", homedir, rank);
  }

  logger= fopen(logName, "w+");

  // Reset the communication matrix with -1 implying no communication initially
  int val = -1;
  reset_comm_matrix(comm_matrix_sum, val);
  reset_comm_matrix(comm_matrix_max, val);
  reset_comm_matrix(comm_matrix_min, val);
  reset_comm_matrix(p2p_mat, val);
  reset_count_matrix(comm_matrix_count);
  return result;
}
int MPI_Finalize() {
  parent_func = __FUNCTION__;
#ifdef LOG_MPI_CALLS
  log_calls[__FUNCTION__]++;
#endif
#ifdef START_FUNC
  std::cout <<"In Finalize\n";
#endif
  // Logging comm_matrix:
  std::string s, st;
  populateCommMatrix(comm_matrix_sum, s, "Communication Sum Matrix");
  fprintf(logger,"%s", s.c_str());
  populateCommMatrix(comm_matrix_max, s, "Communication Max Matrix");
  fprintf(logger,"%s", s.c_str());
  populateCommMatrix(comm_matrix_min, s, "Communication Min Matrix");
  fprintf(logger,"%s", s.c_str());
  populateCommMatrix(p2p_mat, s, "P2P Matrix");
  fprintf(logger,"%s", s.c_str());

  populateCountMatrix(comm_matrix_count, st, "Count Matrix");
  
  // fprintf(logger,"%s", s.c_str());
  fprintf(logger,"%s\n", st.c_str());
  
  // Logging MPI calls
  for(auto &x:log_calls) {
    fprintf(logger, "%s %d\n", x.first.c_str(), x.second);
  }
  // // Logging MPI calls with callsites
  // for(auto &x:log_calls_w_cs) {
  //   fprintf(logger, "%s %d\n", x.first.c_str(), x.second);
  // }
  
  // // Logging #sends per MPI_Collectives
  // for(auto &x:sends) {
  //   fprintf(logger, "%s\n", x.first.c_str());
  //   for(int i=0;i<comm_size;i++) {
  //     for(int j=0;j<comm_size;j++) {
  //       fprintf(logger, "%d\t", x.second[i][j]);
  //     }
  //     fprintf(logger, "%s", "\n");
  //   }
  //   fprintf(logger, "%s", "\n");
  // }

  fclose(logger);
  fprintf(stderr, "[%d] wrapping up\n", rank);
  // std::cout << "In Finalize\n";
  return PMPI_Finalize();
}

int MPI_Send(const void *buffer, int count, MPI_Datatype datatype, int dest,
             int tag, MPI_Comm comm) {
  parent_func = __FUNCTION__;
#ifdef LOG_MPI_CALLS
  log_calls[__FUNCTION__]++;
#endif

  int size;
  PMPI_Type_size(datatype, &size); /* Compute size */  
  
#ifdef START_FUNC
  std::cout << "In Send "<< (size*count) <<"\n";
#endif
  // Creating a p2p matrix
  p2p_mat[rank][dest] += p2p_mat[rank][dest] == -1 ? 1+ size * count : size * count;
  p2p_mat[dest][rank] += p2p_mat[dest][rank] == -1 ? 1+ size * count : size * count;

  updateCommMatrix(size * count, rank, dest);

  int result = PMPI_Send(buffer, count, datatype, dest, tag, comm);
#ifdef DEST_LOG
  logging_dst("MPI_Send", count, size, " Dst: ", dest);
#endif
  
  return result;
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status) {
  parent_func = __FUNCTION__;
#ifdef LOG_MPI_CALLS
  log_calls[__FUNCTION__]++;
#endif  
  int actual_count, size;
  int result = PMPI_Recv(buf, count, datatype, source, tag, comm, status);
  PMPI_Type_size(datatype, &size);                 /* Compute size */
  PMPI_Get_count(status, datatype, &actual_count); /* Compute count */
#ifdef DEST_LOG
  logging_dst("MPI_Recv", count, size, " Src: ", source);
#endif
  return result;
}

int MPI_Isend(const void *buffer, int count, MPI_Datatype datatype, int dest,
             int tag, MPI_Comm comm, MPI_Request *request) {
  parent_func = __FUNCTION__;
#ifdef LOG_MPI_CALLS
  log_calls[__FUNCTION__]++;
#endif
#ifdef START_FUNC
  std::cout <<"In " << parent_func << "\n";
#endif
  // std::cout << "In Send\n";
  int size;
  PMPI_Type_size(datatype, &size); /* Compute size */  

  // Creating a p2p matrix
  p2p_mat[rank][dest] += p2p_mat[rank][dest] == -1 ? 1+ size * count : size * count;
  p2p_mat[dest][rank] += p2p_mat[dest][rank] == -1 ? 1+ size * count : size * count;
  
  updateCommMatrix(size * count, rank, dest);

  int result = PMPI_Isend(buffer, count, datatype, dest, tag, comm, request);
#ifdef DEST_LOG
  logging_dst("MPI_Isend", count, size, " Dst: ", dest);
#endif
  return result;
}
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Request * request) {
  parent_func = __FUNCTION__;
#ifdef LOG_MPI_CALLS
  log_calls[__FUNCTION__]++;
#endif  
  int actual_count, size;
  int result = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
  PMPI_Type_size(datatype, &size);                 /* Compute size */
  // PMPI_Get_count(status, datatype, &actual_count); /* Compute count */
#ifdef DEST_LOG
  logging_dst("MPI_Irecv", count, size, " Src: ", source);
#endif
  return result;
}

int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, 
  void* recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
  parent_func = __FUNCTION__;
#ifdef LOG_MPI_CALLS
    log_calls[__FUNCTION__]++;
#endif  
  #ifdef START_FUNC
    std::cout <<"In " << parent_func << "\n";
  #endif
    int sendtype_size;
    PMPI_Type_size(sendtype, &sendtype_size); /* Compute size */
    
    // Creating a p2p matrix
    p2p_mat[rank][dest] += p2p_mat[rank][dest] == -1 ? 1+ sendcount*sendtype_size : sendcount*sendtype_size;
    p2p_mat[dest][rank] += p2p_mat[dest][rank] == -1 ? 1+ sendcount*sendtype_size : sendcount*sendtype_size;

    updateCommMatrix(sendcount*sendtype_size, rank, dest);
    // fprintf(logger,"%s", "MPI_SendRecv\n");
    // logging("MPI_SendRecv", sendcount, sendtype_size);
#ifdef DEST_LOG
    logging_dst("MPI_Sendrecv", sendcount, sendtype_size, " Dst: ", dest);
#endif
    int result = PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm,status);
    // fprintf(logger,"%d %d\n", rank, status->MPI_SOURCE);
    return result;
}

int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
              MPI_Comm comm) {
  parent_func = __FUNCTION__;
  #ifdef LOG_MPI_CALLS
    log_calls[__FUNCTION__]++;
  #endif
  #ifdef START_FUNC
    std::cout <<"In " << parent_func << "\n";
  #endif
  MPI_Aint nbytes = 0;
  int type_size;
  PMPI_Type_size(datatype, &type_size); /* Compute size */
  nbytes = type_size * count;
  /* Profiling Starts */
  if((nbytes > 0) && (comm_size > 1)) {
    if((nbytes < MPIR_CVAR_BCAST_SHORT_MSG_SIZE) || ((comm_size < MPIR_CVAR_BCAST_MIN_PROCS))) {
      // Condition: nbytes < 12KB || comm_size < 8
      // fprintf(logger, "%s", "MPI_Bcast_intra_binomial\n");
#ifdef DEST_LOG
      logging_dst("MPI_Bcast_intra_binomial", count, type_size, " Root: ", root);
#endif
      profile_Bcast_intra_binomial(root, type_size, count);
    } else {
      if((nbytes < MPIR_CVAR_BCAST_LONG_MSG_SIZE) && !(comm_size & (comm_size - 1))){
        // Condition: nbytes < 512KB && comm_size is pof2
        // fprintf(logger, "%s", "MPI_Bcast_intra_scatter_recursive_doubling_allgather\n");
        // logging("MPI_Bcast_intra_scatter_recursive_doubling_allgather", count, type_size);
#ifdef DEST_LOG
        logging_dst("MPI_Bcast_intra_scatter_recursive_doubling_allgather", count, type_size, " Root: ", root);
#endif
        profile_Bcast_intra_scatter_recursive_doubling_allgather(root, type_size, count);
      } else {
        // Implement profile_Bcast_intra_scatter_ring_allgather
#ifdef DEST_LOG
        logging_dst("MPI_Bcast_intra_scatter_ring_allgather", count, type_size, " Root: ", root);
#endif
        // logging("MPI_Bcast_intra_scatter_ring_allgather", count, type_size);
        // fprintf(logger, "%s", "MPI_Bcast_intra_scatter_ring_allgather\n");
        profile_Bcast_intra_scatter_ring_allgather(root, type_size, count);
      }
    } 
  }
  
  int result = PMPI_Bcast(buffer, count, datatype, root, comm);
  // fprintf(stderr, "[%d] bcast %d %d\n", rank, count * size, root);
  // fprintf(logger,"%d %d %d %d\n", rank, Bcast, count * size, root);
  return result;
}
/* Some problem needs to be fixed */
int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
                MPI_Comm comm) {
  parent_func = __FUNCTION__;
  #ifdef LOG_MPI_CALLS
    log_calls[__FUNCTION__]++;
  #endif
  #ifdef START_FUNC
    std::cout <<"In " << parent_func << "\n";
  #endif
  int sendtype_size, recvtype_size;
  int result = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                            recvtype, root, comm);
  PMPI_Type_size(sendtype, &sendtype_size); /* Compute size */
  PMPI_Type_size(recvtype, &recvtype_size); /* Compute size */
  profile_Scatter_binomial(root, sendtype_size, sendcount, recvtype_size, recvcount);
  // logging("MPI_Scatter_Binomial", sendcount, sendtype_size);
#ifdef DEST_LOG
  logging_dst("MPI_Scatter_Binomial", sendcount, sendtype_size, " Root: ", root);
#endif
  // fprintf(stderr, "[%d] scatter %d %d\n", rank, sendcount * size, root);
  // fprintf(logger,"%d %d %d %d\n", rank, Scatter, sendcount * size, root);

  return result;
}

int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
               MPI_Comm comm) {
  parent_func = __FUNCTION__;
  #ifdef LOG_MPI_CALLS
    log_calls[__FUNCTION__]++;
  #endif
  #ifdef START_FUNC
    std::cout <<"In " << parent_func << "\n";
  #endif
  profile_Gather_binomial(root, sendtype, sendcount, recvtype, recvcount);
  int result = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                           recvtype, root, comm);
  int sendtype_size;
  PMPI_Type_size(sendtype, &sendtype_size);
  // logging("MPI_Gather_Binomial", sendcount, sendtype_size);
#ifdef DEST_LOG
  logging_dst("MPI_Gather_Binomial", sendcount, sendtype_size, " Root: ", root);
#endif
  // fprintf(logger,"MPI_Gather_Binomial\n");
  return result;
}


// Change function body from here
int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm) {
  parent_func = __FUNCTION__;
  #ifdef LOG_MPI_CALLS
    log_calls[__FUNCTION__]++;
  #endif
  #ifdef START_FUNC
    std::cout<<"in reduce\n";
  #endif
  int type_size;
  int result = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  
  PMPI_Type_size(datatype, &type_size); /* Compute size */
  
  int pof2 = (int)pow(2, (int)log2(comm_size));
  if ((count * type_size > MPIR_CVAR_REDUCE_SHORT_MSG_SIZE) &&
      (HANDLE_GET_KIND(op) == HANDLE_KIND_BUILTIN) && (count >= pof2)) {
      /* do a reduce-scatter followed by gather to root. */
      // Condition: count*type_size > 2KB and Built-in Op and count >=pof2
      profile_Reduce_intra_reduce_scatter_gather(datatype, root, count);
      // fprintf(logger,"MPI_reduce_scatter_gather\n");
      // logging("MPI_reduce_scatter_gather", count, type_size);
#ifdef DEST_LOG
      logging_dst("MPI_reduce_scatter_gather", count, type_size, " Root: ", root);
#endif
  } else {
      /* use a binomial tree algorithm */
      profile_Reduce_intra_binomial(datatype, root, count, op);
      // fprintf(logger,"MPI_reduce_intra_binomial\n");
#ifdef DEST_LOG
      logging_dst("MPI_reduce_intra_binomial", count, type_size, " Root: ", root);
#endif
      // logging("MPI_reduce_intra_binomial", count, type_size);
  }

  // fprintf(logger,"MPI_Reduce\n");
  return result;
}

int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  parent_func = __FUNCTION__;
  #ifdef LOG_MPI_CALLS
    log_calls[__FUNCTION__]++;
  #endif
  #ifdef START_FUNC
    std::cout <<"In " << parent_func << "\n";
  #endif
  MPI_Aint tot_bytes;
  // papi_log_compute();
  int type_size;
  PMPI_Type_size(recvtype, &type_size); /* Compute size */
  // Implement our algorithm for profiling the code:
  tot_bytes = (MPI_Aint) recvcount *comm_size * type_size;
  if((tot_bytes < MPIR_CVAR_ALLGATHER_LONG_MSG_SIZE) && !(comm_size & (comm_size - 1))) {
    // Condition: tot_bytes < 512KB and np is pof2
    // Intra_recursive_doubling
    profile_Allgather_intra_recursive_doubling(recvcount, type_size);
    // fprintf(logger,"%s", "Allgather_RD\n");
    logging("Allgather_RD", sendcount, type_size);
  } else if (tot_bytes < MPIR_CVAR_ALLGATHER_SHORT_MSG_SIZE) {
    // Condition: tot_bytes < 80KB
    // INTRA_BRUCKS
    profile_Allgather_intra_brucks(recvcount, type_size);
    // fprintf(logger,"%s", "Allgather_Brucks\n");
    logging("Allgather_Brucks", sendcount, type_size);
  } else {
    // INTRA_RING
    profile_Allgather_intra_ring(recvcount, type_size);
    // fprintf(logger, "%s", "Allgather_Ring\n");
    logging("Allgather_Ring", sendcount, type_size);
  }
  int result = PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                           recvtype, comm);
  
  // fprintf(logger, "MPI_Allgather\n");
  return result;
}

int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, MPI_Comm comm) {
  parent_func = __FUNCTION__;
  #ifdef LOG_MPI_CALLS
    log_calls[__FUNCTION__]++;
  #endif
  #ifdef START_FUNC
    std::cout <<"In " << parent_func << "\n";
  #endif
  int type_size;
  int pof2 = (int)pow(2, (int)log2(comm_size));
  PMPI_Type_size(datatype, &type_size); /* Compute size */
  MPI_Aint nbytes = 0;
  nbytes = type_size * count;
  if ((nbytes <= MPIR_CVAR_ALLREDUCE_SHORT_MSG_SIZE) ||
        (HANDLE_GET_KIND(op) != HANDLE_KIND_BUILTIN) || (count < pof2)) {
    // nbytes <= 2KB || count < pof2
    profile_Allreduce_intra_recursive_doubling(count, type_size);
    // fprintf(logger,"%s", "Allreduce_intra_recursive_doubling\n");
    logging("Allreduce_intra_recursive_doubling", count, type_size);
    // fprintf(logger,"%s %d %s", "Allreduce_intra_recursive_doubling ", count , "\n");
  } else {
    profile_Allreduce_intra_reduce_scatter_allgather(count, type_size);
    logging("Allreduce_intra_reduce_scatter_allgather", count, type_size);
    // fprintf(logger,"%s %d %d %s", "Allreduce_intra_reduce_scatter_allgather - Count: ", count , " ", type_size, "\n");
  }
  
  int result = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  
  // fprintf(logger, "MPI_Allreduce\n");
  
  return result;
}
int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
               int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
  parent_func = __FUNCTION__;
  #ifdef LOG_MPI_CALLS
    log_calls[__FUNCTION__]++;
  #endif
  #ifdef START_FUNC
    std::cout <<"In " << parent_func << "\n";
  #endif
  int sendtype_size, nbytes;
  PMPI_Type_size(sendtype, &sendtype_size); /* Compute size */
  nbytes = sendtype_size * sendcount;
  if (sendbuf == MPI_IN_PLACE) {
    profile_Alltoall_intra_pairwise_sendrecv_replace(sendtype, recvtype, sendcount, recvcount);
    fprintf(logger,"%s", "Alltoall_intra_pairwise_sendrecv_replace\n");
  } else if ((nbytes <= MPIR_CVAR_ALLTOALL_SHORT_MSG_SIZE) && (comm_size >= 8)) {
    profile_Alltoall_intra_brucks(sendtype, recvtype, sendcount, recvcount);
    fprintf(logger,"%s", "Alltoall_intra_brucks\n");
  } else if (nbytes <= MPIR_CVAR_ALLTOALL_MEDIUM_MSG_SIZE) {
    profile_Alltoall_intra_scattered(sendtype, recvtype, sendcount, recvcount);
    fprintf(logger,"%s", "Alltoall_intra_scattered\n");
  } else {
    profile_Alltoall_intra_pairwise(sendtype, recvtype, sendcount, recvcount);
    fprintf(logger,"%s", "Alltoall_intra_pairwise\n");
  }
  int result = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                           recvtype, comm);
  
  // fprintf(logger, "MPI_Alltoall\n");
  return result;
}

int MPI_Barrier(MPI_Comm comm) {
  parent_func = __FUNCTION__;
  #ifdef LOG_MPI_CALLS
    log_calls[__FUNCTION__]++;
  #endif
  // Profile Barrier
  // profile_Barrier();
  int result = PMPI_Barrier(comm);

  return result;  
}

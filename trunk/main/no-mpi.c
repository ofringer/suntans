/*
 * File: no-mpi.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains dummy mpi definitions for use when the mpich libraries
 * are not employed and is meant for use with the serial version of suntans.
 * It is automatically compiled when MPIHOME is not set in Makefile.in.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include<stdlib.h>
#include<sys/time.h>
#include "suntans.h"
#include "string.h"
#include "no-mpi.h"

#define evaltime(timeval_time) (double)timeval_time.tv_sec \
                             + (double)timeval_time.tv_usec*1e-6 

void MPI_Init(int *argc, char ***argv) {}

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *comm_out ) {
  *comm_out = comm;
  return 0;
}

void MPI_Comm_free(MPI_Comm *comm) {}

void MPI_Finalize(void) {
  exit(0);
}

void MPI_Comm_size(MPI_Comm comm, int *size) {
  *size = 1;
}

void MPI_Comm_rank(MPI_Comm comm, int *rank) {
  *rank = 0;
}

int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, 
	     int tag, MPI_Comm comm ) {}
  
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, 
	     int tag, MPI_Comm comm, MPI_Status *status ) {}

int MPI_Bcast (void *buffer, int count, MPI_Datatype datatype, int root, 
	       MPI_Comm comm ) {}

void MPI_Barrier(MPI_Comm comm) {} 

int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
	      MPI_Comm comm, MPI_Request *request ) {}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, 
	      MPI_Comm comm, MPI_Request *request) {}

int MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]) {}

int MPI_Reduce (void *sendbuf, void *recvbuf, int count, 
		MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm ) {
  memcpy(recvbuf,sendbuf,datatype);
}

int MPI_Gather (void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
		void *recvbuf, int recvcount, MPI_Datatype recvtype, 
		int root, MPI_Comm comm ) {
  memcpy(recvbuf,sendbuf,sendtype);
}

double MPI_Wtime(void) {
  struct timeval timeval_time;
  gettimeofday(&timeval_time,NULL);
  return evaltime(timeval_time);
}

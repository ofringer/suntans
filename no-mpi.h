/*
 * File: no-mpi.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for no-mpi.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _mpi_h
#define _mpi_h

#define MPI_DOUBLE 8
#define MPI_INT 4
#define MPI_COMM_WORLD 0
#define MPI_SUM 3
#define MPI_STATUS_IGNORE NULL
#define MPI_MIN 0
#define MPI_MAX 0
#define MPI_IN_PLACE 0

typedef int MPI_Comm;
typedef int MPI_Status;
typedef int MPI_Request;
typedef int MPI_Datatype;
typedef int MPI_Op;

void MPI_Init(int *argc, char ***argv);
int MPI_Comm_dup(MPI_Comm comm,MPI_Comm *comm_out );
void MPI_Comm_free(MPI_Comm *comm);
void MPI_Finalize(void);
void MPI_Comm_size(MPI_Comm comm, int *size);
void MPI_Comm_rank(MPI_Comm comm, int *rank);
int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, 
	     int tag, MPI_Comm comm );
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, 
	     int tag, MPI_Comm comm, MPI_Status *status );
int MPI_Bcast (void *buffer, int count, MPI_Datatype datatype, int root, 
	       MPI_Comm comm );
void MPI_Barrier(MPI_Comm comm);
int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
	      MPI_Comm comm, MPI_Request *request );
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, 
	      MPI_Comm comm, MPI_Request *request);
int MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
int MPI_Reduce (void *sendbuf, void *recvbuf, int count, 
		MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm );
int MPI_Gather (void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
		void *recvbuf, int recvcount, MPI_Datatype recvtype, 
		int root, MPI_Comm comm );
double MPI_Wtime(void);

#endif

/*
 * spreadout_alltoall.cpp
 *
 *  Created on: Sep 2, 2022
 *      Author: kokofan
 */

#include "radix_r_bruck.h"

void spreadout_alltoall(char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm) {

	int rank, nprocs;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &nprocs);

	int typesize;
	MPI_Type_size(sendtype, &typesize);
	int unit_size = sendcount * typesize;

	MPI_Request* req = (MPI_Request*)malloc(2*nprocs*sizeof(MPI_Request));
	MPI_Status* stat = (MPI_Status*)malloc(2*nprocs*sizeof(MPI_Status));
	for (int i = 0; i < nprocs; i++) {
		int src = (rank + i) % nprocs; // avoid always to reach first master node

		if (rank == 1)
			std::cout << "send " << rank << " " << i << " " << src << " " << std::endl;
		MPI_Irecv(&recvbuf[src*recvcount*typesize], recvcount*typesize, MPI_CHAR, src, 0, comm, &req[i]);
	}

	for (int i = 0; i < nprocs; i++) {
		int dst = (rank - i + nprocs) % nprocs;

		if (rank == 1)
			std::cout << "receive " << rank << " " << i << " " << dst << " " << std::endl;
		MPI_Isend(&sendbuf[dst*sendcount*typesize], sendcount*typesize, MPI_CHAR, dst, 0, comm, &req[i+nprocs]);
	}

	MPI_Waitall(2*nprocs, req, stat);
	free(req);
	free(stat);
}



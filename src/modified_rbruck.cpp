/*
 * modified_r_bruck.cpp
 *
 *  Created on: Sep 1, 2022
 *      Author: kokofan
 */

#include "radix_r_bruck.h"

void uniform_modified_radix_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm) {

	int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int typesize;
    MPI_Type_size(sendtype, &typesize);

    int unit_size = sendcount * typesize;
    int w = ceil(log(nprocs) / log(r)); // calculate the number of digits when using r-representation
	int nlpow = pow(r, w-1);
	int d = (pow(r, w) - nprocs) / nlpow; // calculate the number of highest digits

	// convert rank to base r representation
    int* rank_r_reps = (int*) malloc(nprocs * w * sizeof(int));
	for (int i = 0; i < nprocs; i++) {
		std::vector<int> r_rep = convert10tob(w, i, r);
		std::memcpy(&rank_r_reps[i*w], r_rep.data(), w*sizeof(int));
	}

	for (int i = 0; i < nprocs; i++) {
		int index = (2*rank-i+nprocs)%nprocs;
		memcpy(recvbuf+(index*unit_size), sendbuf+(i*unit_size), unit_size);
	}

	int sent_blocks[nlpow];
	int di = 0;
	int ci = 0;

	int comm_steps = (r - 1)*w - d;
	char* temp_buffer = (char*)malloc(nlpow * unit_size); // temporary buffer

    for (int x = 0; x < w; x++) {
    	int ze = (x == w - 1)? r - d: r;
    	for (int z = 1; z < ze; z++) {
    		// get the sent data-blocks
    		// copy blocks which need to be sent at this step
    		di = 0;
    		ci = 0;
    		for (int i = 0; i < nprocs; i++) {
    			if (rank_r_reps[i*w + x] == z){
    				int id = (i + rank) % nprocs;
    				memcpy(&temp_buffer[unit_size*ci++], &recvbuf[id*unit_size], unit_size);
    				sent_blocks[di++] = id;
    			}
    		}

    		// send and receive
    		int distance = z * myPow(r, x); // pow(1, 51) = 51, int d = pow(1, 51); // 50
    		int recv_proc = (rank + distance) % nprocs; // receive data from rank - 2^step process
    		int send_proc = (rank - distance + nprocs) % nprocs; // send data from rank + 2^k process
    		long long comm_size = di * unit_size;
    		MPI_Sendrecv(temp_buffer, comm_size, MPI_CHAR, send_proc, 0, sendbuf, comm_size, MPI_CHAR, recv_proc, 0, comm, MPI_STATUS_IGNORE);

    		// replace with received data
    		for (int i = 0; i < di; i++)
    		{
    			long long offset = sent_blocks[i] * unit_size;
    			memcpy(recvbuf+offset, sendbuf+(i*unit_size), unit_size);
    		}
    	}
    }

	free(rank_r_reps);
	free(temp_buffer);
}



/*
 * basic_rbruck.cpp
 *
 *  Created on: Sep 1, 2022
 *      Author: kokofan
 */

#include "radix_r_bruck.h"

void uniform_radix_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm) {

    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int typesize;
    MPI_Type_size(sendtype, &typesize);

    int unit_size = sendcount * typesize;
    int w = ceil(log(nprocs) / log(r)); // calculate the number of digits when using r-representation
	int nlpow = myPow(r, w-1);
	int d = (myPow(r, w) - nprocs) / nlpow; // calculate the number of highest digits

    // local rotation
    std::memcpy(recvbuf, sendbuf, nprocs*unit_size);
    std::memcpy(&sendbuf[(nprocs - rank)*unit_size], recvbuf, rank*unit_size);
    std::memcpy(sendbuf, &recvbuf[rank*unit_size], (nprocs-rank)*unit_size);

    // convert rank to base r representation
    int* rank_r_reps = (int*) malloc(nprocs * w * sizeof(int));
	for (int i = 0; i < nprocs; i++) {
		std::vector<int> r_rep = convert10tob(w, i, r);
		std::memcpy(&rank_r_reps[i*w], r_rep.data(), w*sizeof(int));
	}

	int sent_blocks[nlpow];
	int sent_blocks_comp[nlpow];
	int di = 0;
	int ci = 0;

	int comm_steps = (r - 1)*w - d;
	char* temp_buffer = (char*)malloc(nlpow * unit_size); // temporary buffer

	// communication steps = (r - 1)w - d
    for (int x = 0; x < w; x++) {
    	int ze = (x == w - 1)? r - d: r;
    	for (int z = 1; z < ze; z++) {
    		// get the sent data-blocks
    		// copy blocks which need to be sent at this step
    		di = 0;
    		ci = 0;
    		for (int i = 0; i < nprocs; i++) {
    			if (rank_r_reps[i*w + x] == z) {
    				sent_blocks[di++] = i;
    				memcpy(&temp_buffer[unit_size*ci++], &sendbuf[i*unit_size], unit_size);
    			}
    		}

    		// send and receive
    		int distance = z * myPow(r, x); // pow(1, 51) = 51, int d = pow(1, 51); // 50
    		int recv_proc = (rank - distance + nprocs) % nprocs; // receive data from rank - 2^step process
    		int send_proc = (rank + distance) % nprocs; // send data from rank + 2^k process
    		long long comm_size = di * unit_size;
    		MPI_Sendrecv(temp_buffer, comm_size, MPI_CHAR, send_proc, 0, recvbuf, comm_size, MPI_CHAR, recv_proc, 0, comm, MPI_STATUS_IGNORE);

    		// replace with received data
    		for (int i = 0; i < di; i++) {
    			long long offset = sent_blocks[i] * unit_size;
    			memcpy(sendbuf+offset, recvbuf+(i*unit_size), unit_size);
    		}
    	}
    }

    free(rank_r_reps);
    free(temp_buffer);

    // local rotation
	for (int i = 0; i < nprocs; i++) {
		int index = (rank - i + nprocs) % nprocs;
		memcpy(&recvbuf[index*unit_size], &sendbuf[i*unit_size], unit_size);
	}
}


void optimized_radix_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm) {

    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int typesize;
    MPI_Type_size(sendtype, &typesize);

    int unit_size = sendcount * typesize;
    int w = ceil(log(nprocs) / log(r)); // calculate the number of digits when using r-representation
	int nlpow = myPow(r, w-1);
	int d = (myPow(r, w) - nprocs) / nlpow; // calculate the number of highest digits

    // local rotation
    std::memcpy(recvbuf, sendbuf, nprocs*unit_size);
    std::memcpy(&sendbuf[(nprocs - rank)*unit_size], recvbuf, rank*unit_size);
    std::memcpy(sendbuf, &recvbuf[rank*unit_size], (nprocs-rank)*unit_size);

	int sent_blocks[nlpow];
	int di = 0, ci = 0;

	int comm_steps = (r - 1)*w - d;
	char* temp_buffer = (char*)malloc(nlpow * unit_size); // temporary buffer

	// communication steps = (r - 1)w - d
	int spoint = 1, distance = 1, next_distance = r;
    for (int x = 0; x < w; x++) {
    	for (int z = 1; z < r; z++) {
    		di = 0; ci = 0;
    		spoint = z * distance;
    		if (spoint > nprocs - 1) {break;}

    		// get the sent data-blocks
    		for (int i = spoint; i < nprocs; i += next_distance) {
    			for (int j = i; j < (i+distance); j++) {
    				if (j > nprocs - 1 ) { break; }
    				sent_blocks[di++] = j;
    				memcpy(&temp_buffer[unit_size*ci++], &sendbuf[j*unit_size], unit_size);
    			}
    		}

    		// send and receive
    		int recv_proc = (rank - spoint + nprocs) % nprocs; // receive data from rank - 2^step process
    		int send_proc = (rank + spoint) % nprocs; // send data from rank + 2^k process
    		long long comm_size = di * unit_size;
    		MPI_Sendrecv(temp_buffer, comm_size, MPI_CHAR, send_proc, 0, recvbuf, comm_size, MPI_CHAR, recv_proc, 0, comm, MPI_STATUS_IGNORE);

    		// replace with received data
    		for (int i = 0; i < di; i++) {
    			long long offset = sent_blocks[i] * unit_size;
    			memcpy(sendbuf+offset, recvbuf+(i*unit_size), unit_size);
    		}
    	}
		distance *= r;
		next_distance *= r;
    }
    free(temp_buffer);

    // local rotation
	for (int i = 0; i < nprocs; i++) {
		int index = (rank - i + nprocs) % nprocs;
		memcpy(&recvbuf[index*unit_size], &sendbuf[i*unit_size], unit_size);
	}
}


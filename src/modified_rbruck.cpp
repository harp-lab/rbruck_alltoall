/*
 * modified_r_bruck.cpp
 *
 *  Created on: Sep 1, 2022
 *      Author: kokofan
 */

#include "radix_r_bruck.h"

std::vector<double> iteTimes;

void uniform_modified_radix_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype,
		char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm) {

	int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int typesize;
    MPI_Type_size(sendtype, &typesize);

    int unit_size = sendcount * typesize;
    int w = ceil(log(nprocs) / log(r)); // calculate the number of digits when using r-representation
	int nlpow = pow(r, w-1);

	for (int i = 0; i < nprocs; i++) {
		int index = (2*rank-i+nprocs)%nprocs;
		memcpy(recvbuf+(index*unit_size), sendbuf+(i*unit_size), unit_size);
	}

	int sent_blocks[nlpow];
	int di = 0;
	int ci = 0;

	char* temp_buffer = (char*)malloc(nlpow * unit_size); // temporary buffer
	int spoint = 1, distance = 1, next_distance = r;
    for (int x = 0; x < w; x++) {
    	for (int z = 1; z < r; z++) {

    		double st = MPI_Wtime();
    		// get the sent data-blocks
    		// copy blocks which need to be sent at this step
    		spoint = z * distance;
    		if (spoint > nprocs - 1) {break;}
    		di = 0; ci = 0;
			for (int i = spoint; i < nprocs; i += next_distance) {
				for (int j = i; j < (i+distance); j++) {
					if (j > nprocs - 1 ) { break; }
					int id = (j + rank) % nprocs;
					sent_blocks[di++] = id;
					memcpy(&temp_buffer[unit_size*ci++], &recvbuf[id*unit_size], unit_size);
				}
			}

    		// send and receive
    		int recv_proc = (rank + spoint) % nprocs; // receive data from rank - 2^step process
    		int send_proc = (rank - spoint + nprocs) % nprocs; // send data from rank + 2^k process
    		long long comm_size = di * unit_size;
    		MPI_Sendrecv(temp_buffer, comm_size, MPI_CHAR, send_proc, 0, sendbuf, comm_size, MPI_CHAR, recv_proc, 0, comm, MPI_STATUS_IGNORE);

    		// replace with received data
    		for (int i = 0; i < di; i++) {
    			long long offset = sent_blocks[i] * unit_size;
    			memcpy(recvbuf+offset, sendbuf+(i*unit_size), unit_size);
    		}

    		double et = MPI_Wtime();
    		double curtime = et - st;

    		iteTimes.push_back(curtime);
    	}
		distance *= r;
		next_distance *= r;
    }

	free(temp_buffer);
}



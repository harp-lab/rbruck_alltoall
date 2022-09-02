/*
 * group_rbruck.cpp
 *
 *  Created on: Sep 1, 2022
 *      Author: kokofan
 */

#include "radix_r_bruck.h"

void uniform_isplit_r_bruck(int n, int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm) {

	int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int typesize;
    MPI_Type_size(sendtype, &typesize);

    int unit_size = sendcount * typesize;

	int ngroup = nprocs / n; // number of groups
	int sw = ceil(log(n) / log(r)); // required digits for intra-Bruck
	int gw = ceil(log(ngroup) / log(r)); // required digits for inter-Bruck

	int nlpow = pow(r, sw-1); // the largest power of r that smaller than n
	int sd = (pow(r, sw) - n) / nlpow;

	int glpow = pow(r, gw-1); // the largest power of r that smaller than ngroup
	int gd = (pow(r, gw) - ngroup) / glpow;

	int maxNum = (n > ngroup)? n: ngroup;
	int maxW = (gw > sw)? gw: sw;

	// convert rank to base r representation
	int* rank_r_reps = (int*) malloc(maxNum * maxW * sizeof(int));
	for (int i = 0; i < maxNum; i++) {
		std::vector<int> r_rep = convert10tob(maxW, i, r);
		std::memcpy(&rank_r_reps[i*maxW], r_rep.data(), maxW*sizeof(int));
	}

	int grank = rank % n; // rank of each process in a group
	int gid = rank / n; // group id

	// Initial rotation phase for intra-Bruck
	for (int i = 0; i < ngroup; i++) {
		int gsp = i*n;
		for (int j = 0; j < n; j++) {
			int index = gsp + (2 * grank - j + n) % n;
			memcpy(recvbuf+(index*unit_size), sendbuf+((i*n+j)*unit_size), unit_size);
		}
	}

	int max_sd = glpow * n; // max send data block count

	int sent_blocks[max_sd];
	int di = 0, ci = 0;

	char* temp_buffer = (char*)malloc(max_sd * unit_size); // temporary buffer

	// Intra-Bruck
    for (int x = 0; x < sw; x++) {
    	int ze = (x == sw - 1)? r - sd: r;
    	for (int z = 1; z < ze; z++) {
    		di = 0; ci = 0;
    		for (int i = 0; i < nprocs; i++) { // get the sent data-blocks
    			int gi = i % n;
    			int gn = i / n;
    			if (rank_r_reps[gi*maxW + x] == z){
    				int id = gn*n + (gi + grank) % n; // copy blocks which need to be sent at this step
    				memcpy(&temp_buffer[unit_size*ci++], &recvbuf[id*unit_size], unit_size);
    				sent_blocks[di++] = id;
    			}
    		}

    		// send and receive
    		int distance = z * myPow(r, x); // pow(1, 51) = 51, int d = pow(1, 51); // 50
    		int recv_proc = gid*n + (grank + distance) % n; // receive data from rank - 2^step process
    		int send_proc = gid*n + (grank - distance + n) % n; // send data from rank + 2^k process

    		long long comm_size = di * unit_size;
    		MPI_Sendrecv(temp_buffer, comm_size, MPI_CHAR, send_proc, 0, sendbuf, comm_size, MPI_CHAR, recv_proc, 0, comm, MPI_STATUS_IGNORE);


    		// replace with received data
    		for (int i = 0; i < di; i++) {
    			long long offset = sent_blocks[i] * unit_size;
    			memcpy(recvbuf+offset, sendbuf+(i*unit_size), unit_size);
    		}
    	}
    }

    unit_size = n * sendcount * typesize;
	// Initial rotation phase for inter-Bruck
	for (int i = 0; i < ngroup; i++) {
		int index = (2 * gid - i + ngroup) % ngroup;
		memcpy(sendbuf+(index*unit_size), recvbuf+(i*unit_size), unit_size);
	}

	// Inter-Bruck
    for (int x = 0; x < gw; x++) {
    	int ze = (x == gw - 1)? r - gd: r;
    	for (int z = 1; z < ze; z++) {
    		di = 0; ci = 0;
    		for (int i = 0; i < ngroup; i++) {
    			if (rank_r_reps[i*maxW + x] == z){
    				int id = (i + gid) % ngroup;
    				memcpy(&temp_buffer[unit_size*ci++], &sendbuf[id*unit_size], unit_size);
    				sent_blocks[di++] = id;
    			}
    		}

    		int distance = z * myPow(r, x) * n;
    		int recv_proc = (gid*n + (grank + distance)) % nprocs; // receive data from rank - 2^step process
    		int send_proc = (gid*n + (grank - distance + nprocs)) % nprocs; // send data from rank + 2^k process

    		long long comm_size = di * unit_size;
    		MPI_Sendrecv(temp_buffer, comm_size, MPI_CHAR, send_proc, 0, recvbuf, comm_size, MPI_CHAR, recv_proc, 0, comm, MPI_STATUS_IGNORE);

    		for (int i = 0; i < di; i++) {
    			long long offset = sent_blocks[i] * unit_size;
    			memcpy(sendbuf+offset, recvbuf+(i*unit_size), unit_size);
    		}
    	}
    }
    memcpy(recvbuf, sendbuf, nprocs*sendcount*typesize);

	free(rank_r_reps);
	free(temp_buffer);
}


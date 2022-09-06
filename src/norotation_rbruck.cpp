/*
 * norotation_rbruck.cpp
 *
 *  Created on: Sep 1, 2022
 *      Author: kokofan
 */

#include "radix_r_bruck.h"

void uniform_norot_radix_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm) {

    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int typesize;
    MPI_Type_size(sendtype, &typesize);

    int unit_size = sendcount * typesize;
    int w = ceil(log(nprocs) / log(r)); // calculate the number of digits when using r-representation
	int nlpow = myPow(r, w-1);
	int d = (myPow(r, w) - nprocs) / nlpow; // calculate the number of highest digits

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
    		int tsd = 0;
    		for (int i = spoint; i < nprocs; i += next_distance) {
    			int dis = (i + distance < nprocs)? distance: (nprocs - i);
    			tsd += dis;
//    			for (int j = i; j < (i+distance); j++) {
//    				if (j > nprocs - 1 ) { break; }
//    				int id = (j + rank) % nprocs;
//    				sent_blocks[di++] = id;
//    				memcpy(&temp_buffer[unit_size*ci++], &sendbuf[id*unit_size], unit_size);
//    			}
    		}
//
//    		// send and receive
//    		int recv_proc = (rank - spoint + nprocs) % nprocs; // receive data from rank - 2^step process
//    		int send_proc = (rank + spoint) % nprocs; // send data from rank + 2^k process
//    		long long comm_size = di * unit_size;
//    		MPI_Sendrecv(temp_buffer, comm_size, MPI_CHAR, send_proc, 0, recvbuf, comm_size, MPI_CHAR, recv_proc, 0, comm, MPI_STATUS_IGNORE);
//
//    		// replace with received data
//    		for (int i = 0; i < di; i++) {
//    			long long offset = sent_blocks[i] * unit_size;
//    			memcpy(sendbuf+offset, recvbuf+(i*unit_size), unit_size);
//    		}
    	}
		distance *= r;
		next_distance *= r;
    }
    free(temp_buffer);


//	for (int i = 0; i < nprocs; i++) {
//		if (rank == 3) {
//			long long a;
//			memcpy(&a, &sendbuf[i*unit_size], unit_size);
//			std::cout << a << std::endl;
//		}
//	}

}

void uniform_norotation_radix_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm) {

    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int typesize;
    MPI_Type_size(sendtype, &typesize);

    int unit_size = sendcount * typesize;
    int w = ceil(log(nprocs) / log(r)); // calculate the number of digits when using r-representation
	int nlpow = myPow(r, w-1);
	int d = (myPow(r, w) - nprocs) / nlpow; // calculate the number of highest digits

	// convert rank to base r representation
    int* rank_r_reps = (int*) malloc(nprocs * w * sizeof(int));
	for (int i = 0; i < nprocs; i++) {
		std::vector<int> r_rep = convert10tob(w, i, r);
		std::memcpy(&rank_r_reps[i*w], r_rep.data(), w*sizeof(int));
	}

	// create local index array after rotation
	memcpy(recvbuf+rank*unit_size, sendbuf+rank*unit_size, unit_size);
	int rotate_array[nprocs];
    for (int i = 0; i < nprocs; i++)
    	rotate_array[i] = (2*rank-i+nprocs)%nprocs;

    char* stemp_buffer = (char*)malloc(nlpow * unit_size); // temporary buffer
    char* rtemp_buffer = (char*)malloc(nlpow * unit_size); // temporary buffer

	int sent_blocks[nlpow];
	int di = 0;
	int ci = 0;

    // communication steps = (r - 1)w - d
	for (int x = 0; x < w; x++) {
		int ze = (x == w - 1)? r - d: r;
		for (int z = 1; z < ze; z++) {

			// get the sent data-blocks
			// copy blocks which need to be sent at this step
			di = 0;
			ci = 0;
			for (int i = 0; i < nprocs; i++) {
				if (rank_r_reps[i*w + x] == z){
					int sbs =(i + rank) % nprocs;
					sent_blocks[di++] = sbs;
					memcpy(&stemp_buffer[unit_size*ci++], &sendbuf[rotate_array[sbs]*unit_size], unit_size);
				}
			}

			int distance = z * myPow(r, x);
			int recv_proc = (rank + distance) % nprocs; // receive data from rank - 2^step process
			int send_proc = (rank - distance + nprocs) % nprocs; // send data from rank + 2^k process
			long long comm_size = di * unit_size;
			MPI_Sendrecv(stemp_buffer, comm_size, MPI_CHAR, send_proc, 0, rtemp_buffer, comm_size, MPI_CHAR, recv_proc, 0, comm, MPI_STATUS_IGNORE);

			for (int i = 0; i < di; i++)
			{
				long long offset = rotate_array[sent_blocks[i]] * unit_size;
				memcpy(recvbuf+(sent_blocks[i]*unit_size), rtemp_buffer+(i*unit_size), unit_size);
				memcpy(sendbuf+offset, rtemp_buffer+(i*unit_size), unit_size);
			}
		}
	}

	free(rank_r_reps);
	free(stemp_buffer);
	free(rtemp_buffer);
}


void uniform_norotation_radix_r_bruck_dt(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm) {

	int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int typesize;
    MPI_Type_size(sendtype, &typesize);

    int unit_size = sendcount * typesize;
    int w = ceil(log(nprocs) / log(r)); // calculate the number of digits when using r-representation
	int nlpow = myPow(r, w-1);
	int d = (myPow(r, w) - nprocs) / nlpow; // calculate the number of highest digits

	// convert rank to base r representation
    int* rank_r_reps = (int*) malloc(nprocs * w * sizeof(int));
	for (int i = 0; i < nprocs; i++) {
		std::vector<int> r_rep = convert10tob(w, i, r);
		std::memcpy(&rank_r_reps[i*w], r_rep.data(), w*sizeof(int));
	}

	// create local index array after rotation
	memcpy(recvbuf+rank*unit_size, sendbuf+rank*unit_size, unit_size);
	int rotate_array[nprocs];
    for (int i = 0; i < nprocs; i++)
    	rotate_array[i] = (2*rank-i+nprocs)%nprocs;

    char* rtemp_buffer = (char*)malloc(nprocs * unit_size); // temporary buffer

	int sent_blocks[nlpow];
	int send_displs[nlpow];
	int di = 0;
	int ci = 0;

    // communication steps = (r - 1)w - d
	for (int x = 0; x < w; x++) {
		int ze = (x == w - 1)? r - d: r;
		for (int z = 1; z < ze; z++) {

			// get the sent data-blocks
			// copy blocks which need to be sent at this step
			di = 0;
			ci = 0;
			for (int i = 0; i < nprocs; i++) {
				if (rank_r_reps[i*w + x] == z){
					int sbs = (i + rank) % nprocs;
					sent_blocks[di] = sbs;
					send_displs[di] = rotate_array[sbs]*unit_size;
					di++;
				}
			}


			MPI_Datatype send_type;
			MPI_Type_create_indexed_block(di, unit_size, send_displs, MPI_CHAR, &send_type);
			MPI_Type_commit(&send_type);

			int distance = z * myPow(r, x);
			int recv_proc = (rank + distance) % nprocs; // receive data from rank - 2^step process
			int send_proc = (rank - distance + nprocs) % nprocs; // send data from rank + 2^k process
			MPI_Sendrecv(sendbuf, 1, send_type, send_proc, 0, rtemp_buffer, 1, send_type, recv_proc, 0, comm, MPI_STATUS_IGNORE);

			MPI_Type_free(&send_type);


			for (int i = 0; i < di; i++)
			{
				long long offset = rotate_array[sent_blocks[i]] * unit_size;
				memcpy(sendbuf+offset, rtemp_buffer+offset, unit_size);
				memcpy(recvbuf+(sent_blocks[i]*unit_size), rtemp_buffer+offset, unit_size);
			}
		}
	}

	free(rank_r_reps);
	free(rtemp_buffer);

}


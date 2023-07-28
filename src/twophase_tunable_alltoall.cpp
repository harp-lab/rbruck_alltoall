/*
 * twophase_bruck.cpp
 *
 *      Author: kokofan
 */

#include "radix_r_bruck.h"

void twophase_tuable_radix_alltoallv(int r, char *sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, char *recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) {

	int rank, nprocs;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &nprocs);

	if (r > nprocs) r = nprocs;
	if (r < 2) r = 2;

	int typesize;
	MPI_Type_size(sendtype, &typesize);

	// 1. Find max send count
	int local_max_count = 0;
	for (int i = 0; i < nprocs; i++) {
		if (sendcounts[i] > local_max_count)
			local_max_count = sendcounts[i];
	}
	int max_send_count = 0;
	MPI_Allreduce(&local_max_count, &max_send_count, 1, MPI_INT, MPI_MAX, comm);


	int w = ceil(log(nprocs) / log(r)); // calculate the number of digits when using r-representation
	int nlpow = myPow(r, w-1);
	int d = (myPow(r, w) - nprocs) / nlpow; // calculate the number of highest digits



//	if (rank == 0) {
//		int inx = 0;
//		for (int i = 0; i < nprocs; i++) {
//			std::cout << i << ", " << sendcounts[i] << std::endl;
//			for (int j = 0; j < sendcounts[i]; j++) {
//				long a = 0;
//				memcpy(&a, &sendbuf[inx*typesize], typesize);
//				std::cout << a << std::endl;
//				inx++;
//			}
//		}
//	}

	// 2. create local index array after rotation
//	int rotate_index_array[nprocs];
//	for (int i = 0; i < nprocs; i++)
//		rotate_index_array[i] = (2*rank-i+nprocs)%nprocs;

	// 3. exchange data with log(P) steps
//	int max_send_elements = (nprocs+1)/2;
	char* extra_buffer = (char*) malloc(max_send_count*typesize*nprocs);
	char* temp_send_buffer = (char*) malloc(max_send_count*typesize*nlpow);
	char* temp_recv_buffer = (char*) malloc(max_send_count*typesize*nlpow);
	int pos_status[nprocs];
	memset(pos_status, 0, nprocs*sizeof(int));
	memcpy(&recvbuf[rdispls[rank]*typesize], &sendbuf[sdispls[rank]*typesize], recvcounts[rank]*typesize);


	if (rank == 1) {
		int inx = 0;
		for (int i = 0; i < nprocs; i++) {
			std::cout << i << ", " << sendcounts[i] << std::endl;
			for (int j = 0; j < sendcounts[i]; j++) {
				long a = 0;
				memcpy(&a, &sendbuf[inx*typesize], typesize);
				std::cout << a << std::endl;
				inx++;
			}
		}
	}

	int di = 0, ci = 0;
	int spoint = 1, distance = 1, next_distance = r;
	int sent_blocks[nlpow];

	for (int x = 0; x < w; x++) {
		int ze = (x == w - 1)? r - d: r;
		for (int z = 1; z < ze; z++) {

    		di = 0; ci = 0;
    		spoint = z * distance;

//    		// get the sent data-blocks
//    		for (int i = spoint; i < nprocs; i += next_distance) {
//    			for (int j = i; j < (i+distance); j++) {
//    				if (j > nprocs - 1 ) { break; }
//    				sent_blocks[di++] = j;
//    				memcpy(&temp_buffer[unit_size*ci++], &sendbuf[j*unit_size], unit_size);
//    			}
//    		}

		}
	}

//	for (int k = 1; k < nprocs; k <<= 1) {
//		// 1) find which data blocks to send
//		int send_indexes[max_send_elements];
//		int sendb_num = 0;
//		for (int i = k; i < nprocs; i++) {
//			if (i & k)
//				send_indexes[sendb_num++] = (rank+i)%nprocs;
//		}
//
//		// 2) prepare metadata and send buffer
//		int metadata_send[sendb_num];
//		int sendCount = 0;
//		int offset = 0;
//		for (int i = 0; i < sendb_num; i++) {
//			int send_index = rotate_index_array[send_indexes[i]];
//			metadata_send[i] = sendcounts[send_index];
//			if (pos_status[send_index] == 0)
//				memcpy(&temp_send_buffer[offset], &sendbuf[sdispls[send_index]*typesize], sendcounts[send_index]*typesize);
//			else
//				memcpy(&temp_send_buffer[offset], &extra_buffer[send_indexes[i]*max_send_count*typesize], sendcounts[send_index]*typesize);
//			offset += sendcounts[send_index]*typesize;
//		}
//
//		// 3) exchange metadata
//		int sendrank = (rank - k + nprocs) % nprocs;
//		int recvrank = (rank + k) % nprocs;
//		int metadata_recv[sendb_num];
//		MPI_Sendrecv(metadata_send, sendb_num, MPI_INT, sendrank, 0, metadata_recv, sendb_num, MPI_INT, recvrank, 0, comm, MPI_STATUS_IGNORE);
//
//		for(int i = 0; i < sendb_num; i++)
//			sendCount += metadata_recv[i];
//
//		// 4) exchange data
//		MPI_Sendrecv(temp_send_buffer, offset, MPI_CHAR, sendrank, 1, temp_recv_buffer, sendCount*typesize, MPI_CHAR, recvrank, 1, comm, MPI_STATUS_IGNORE);
//
//		// 5) replace
//		offset = 0;
//		for (int i = 0; i < sendb_num; i++) {
//			int send_index = rotate_index_array[send_indexes[i]];
//
//			if ((send_indexes[i] - rank + nprocs) % nprocs < (k << 1))
//				memcpy(&recvbuf[rdispls[send_indexes[i]]*typesize], &temp_recv_buffer[offset], metadata_recv[i]*typesize);
//			else
//				memcpy(&extra_buffer[send_indexes[i]*max_send_count*typesize], &temp_recv_buffer[offset], metadata_recv[i]*typesize);
//
//			offset += metadata_recv[i]*typesize;
//			pos_status[send_index] = 1;
//			sendcounts[send_index] = metadata_recv[i];
//		}
//	}
	free(temp_send_buffer);
	free(temp_recv_buffer);
	free(extra_buffer);
}




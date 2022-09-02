/*
 * r_radix_bruck.cpp
 *
 *  Created on: Jul 09, 2022
 *      Author: kokofan
 */

#include "radix_r_bruck.h"

std::vector<int> convert10tob(int w, int N, int b)
{
	std::vector<int> v(w);
	int i = 0;
	while(N) {
	  v[i++] = (N % b);
	  N /= b;
	}
//	std::reverse(v.begin(), v.end());
	return v;
}

int myPow(int x, unsigned int p) {
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = myPow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}

void uniform_radix_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm)
{
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
	int di = 0;
	int ci = 0;

	int comm_steps = (r - 1)*w - d;
	char* temp_buffer = (char*)malloc(nlpow * unit_size); // temporary buffer

	// communication steps = (r - 1)w - d
//	double pre_time = 0, comm_time = 0, replace_time = 0;
    for (int x = 0; x < w; x++) {
    	int ze = (x == w - 1)? r - d: r;
    	for (int z = 1; z < ze; z++) {
    		// get the sent data-blocks
    		// copy blocks which need to be sent at this step
    		di = 0;
    		ci = 0;
    		for (int i = 0; i < nprocs; i++) {
    			if (rank_r_reps[i*w + x] == z){
    				sent_blocks[di++] = i;
    				memcpy(&temp_buffer[unit_size*ci++], &sendbuf[i*unit_size], unit_size);

    			}
    		}

    		// send and receive
//    		s = MPI_Wtime();
    		int distance = z * myPow(r, x); // pow(1, 51) = 51, int d = pow(1, 51); // 50
    		int recv_proc = (rank - distance + nprocs) % nprocs; // receive data from rank - 2^step process
    		int send_proc = (rank + distance) % nprocs; // send data from rank + 2^k process
    		long long comm_size = di * unit_size;
    		MPI_Sendrecv(temp_buffer, comm_size, MPI_CHAR, send_proc, 0, recvbuf, comm_size, MPI_CHAR, recv_proc, 0, comm, MPI_STATUS_IGNORE);
//    		e = MPI_Wtime();
//    		comm_time += e - s;

//    		s = MPI_Wtime();
    		// replace with received data
    		for (int i = 0; i < di; i++)
    		{
    			long long offset = sent_blocks[i] * unit_size;
    			memcpy(sendbuf+offset, recvbuf+(i*unit_size), unit_size);
    		}
//    		e = MPI_Wtime();
//    		replace_time += e - s;
    	}
    }

    free(rank_r_reps);
    free(temp_buffer);

    // local rotation
//    s = MPI_Wtime();
	for (int i = 0; i < nprocs; i++) {
		int index = (rank - i + nprocs) % nprocs;
		memcpy(&recvbuf[index*unit_size], &sendbuf[i*unit_size], unit_size);
	}
//	e = MPI_Wtime();
//	double second_time = e - s;

//  double te = MPI_Wtime();
//	double total_time = te - ts;

//	timelist[it][0] = total_time;
//	timelist[it][1] = first_time;
//	timelist[it][2] = conv_time;
//	timelist[it][3] = pre_time;
//	timelist[it][4] = comm_time;
//	timelist[it][5] = replace_time;
//	timelist[it][6] = second_time;

}


void uniform_norot_radix_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm)
{
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

	int sent_blocks[nlpow];
	int di = 0;
	int ci = 0;

	int comm_steps = (r - 1)*w - d;
	char* temp_buffer = (char*)malloc(nlpow * unit_size); // temporary buffer

	int distance = 0;
	// communication steps = (r - 1)w - d
    for (int x = 0; x < w; x++) { // 0-w
    	int ze = (x == w - 1)? r - d: r;
    	for (int z = 1; z < ze; z++) { //1-ze
    		// get the sent data-blocks
    		// copy blocks which need to be sent at this step
    		di = 0;
    		ci = 0;
    		for (int i = 0; i < nprocs; i++) {
    			if (rank_r_reps[i*w + x] == z){
    				int id = ((i + rank) % nprocs) - distance;
    				if (rank == 0) {
						long long a;
						memcpy(&a, &sendbuf[id*unit_size], unit_size);
						std::cout << x << " " << z << " " << distance << " " << a << std::endl;
    				}
    				memcpy(&temp_buffer[unit_size*ci++], &sendbuf[id*unit_size], unit_size);
       				sent_blocks[di++] = id;
    			}
    		}


    		// send and receive

    		distance = z * myPow(r, x); // pow(1, 51) = 51, int d = pow(1, 51); // 50
    		int recv_proc = (rank - distance + nprocs) % nprocs; // receive data from rank - 2^step process
    		int send_proc = (rank + distance) % nprocs; // send data from rank + 2^k process
    		long long comm_size = di * unit_size;
    		MPI_Sendrecv(temp_buffer, comm_size, MPI_CHAR, send_proc, 0, recvbuf, comm_size, MPI_CHAR, recv_proc, 0, comm, MPI_STATUS_IGNORE);

    		// replace with received data
    		for (int i = 0; i < di; i++) {
    			int id = (i - distance + di) % di;
//    			if (rank == 3) {
//    				std::cout << i << " " << id << " " << sent_blocks[id] << std::endl;
//    			}
    			long long offset = sent_blocks[id] * unit_size;
    			memcpy(sendbuf+offset, recvbuf+(i*unit_size), unit_size);
    		}
    	}
    }
//
//	for (int i = 0; i < nprocs; i++) {
//
//		if (rank == 3) {
//			long long a;
//			memcpy(&a, &sendbuf[i*unit_size], unit_size);
//			std::cout << a << std::endl;
//		}
//	}

    free(rank_r_reps);
    free(temp_buffer);

    // local rotation
//    s = MPI_Wtime();
//	for (int i = 0; i < nprocs; i++) {
//		int index = (rank - i + nprocs) % nprocs;
//		memcpy(&recvbuf[index*unit_size], &sendbuf[i*unit_size], unit_size);
//	}

}


void uniform_modified_radix_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm)
{
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




void uniform_modified_inverse_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm)
{
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

    for (int x = w-1; x > -1; x--) {
    	int ze = (x == w - 1)? r - d: r;
    	for (int z = ze-1; z > 0; z--) {
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


void uniform_norotation_radix_r_bruck(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm)
{
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

//	if (rank == 3)
//	{
//		for (int a = 0; a < sendcount*nprocs; a++) {
//			long long value;
//			memcpy(&value, &recvbuf[a*sizeof(long long)], sizeof(long long));
//			std::cout << a << ", " << value << std::endl;
//		}
//	}

	free(rank_r_reps);
	free(stemp_buffer);
	free(rtemp_buffer);
}


void uniform_norotation_radix_r_bruck_dt(int r, char *sendbuf, int sendcount, MPI_Datatype sendtype, char *recvbuf, int recvcount, MPI_Datatype recvtype,  MPI_Comm comm)
{
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

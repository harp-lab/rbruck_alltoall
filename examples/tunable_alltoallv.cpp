#include "radix_r_bruck.h"

static int rank, nprocs;
static void run_tunable_alltoallv(int loopcount, int nprocs, std::vector<int> bases, int warmup);

int main(int argc, char **argv) {
    // MPI Initial
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
        std::cout << "ERROR: MPI_Init error\n" << std::endl;
    if (MPI_Comm_size(MPI_COMM_WORLD, &nprocs) != MPI_SUCCESS)
    	std::cout << "ERROR: MPI_Comm_size error\n" << std::endl;
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    	std::cout << "ERROR: MPI_Comm_rank error\n" << std::endl;

    if (argc < 3) {
    	std::cout << "Usage: mpirun -n <nprocs> " << argv[0] << " <loop-count> <base-list> " << std::endl;
    	return -1;
    }

    int loopCount = atoi(argv[1]);
    std::vector<int> bases;
    for (int i = 2; i < argc; i++)
    	bases.push_back(atoi(argv[i]));

    // for warm up only
//    run_tunable_alltoallv(1, nprocs, bases, 1);

    // run code
    run_tunable_alltoallv(loopCount, nprocs, bases, 0);

	MPI_Finalize();
    return 0;
}


static void run_tunable_alltoallv(int loopcount, int nprocs, std::vector<int> bases, int warmup) {

	int basecount = bases.size();

	for (int n = 4; n <= 4; n = n*2) {

		int sendcounts[nprocs]; // the size of data each process send to other process
		memset(sendcounts, 0, nprocs*sizeof(int));
		int sdispls[nprocs];
		int soffset = 0;

		// Uniform random distribution
		srand(time(NULL));
		for (int i=0; i < nprocs; i++) {
			int random = rand() % 100;
			sendcounts[i] = (n * random) / 100;
		}

		// Random shuffling the sendcounts array
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::shuffle(&sendcounts[0], &sendcounts[nprocs], std::default_random_engine(seed));


		// Initial send offset array
		for (int i = 0; i < nprocs; ++i) {
			sdispls[i] = soffset;
			soffset += sendcounts[i];
		}

		// Initial receive counts and offset array
		int recvcounts[nprocs];
		MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
		int rdispls[nprocs];
		int roffset = 0;
		for (int i = 0; i < nprocs; ++i) {
			rdispls[i] = roffset;
			roffset += recvcounts[i];
		}

		// Initial send buffer
		long long* send_buffer = new long long[soffset];
		long long* recv_buffer = new long long[roffset];

		int index = 0;
		for (int i = 0; i < nprocs; i++) {
			for (int j = 0; j < sendcounts[i]; j++)
				send_buffer[index++] = i + rank * 10;
		}

		int scounts[nprocs]; // a copy of sendcounts for each iteration

		MPI_Barrier(MPI_COMM_WORLD);


//		long long* sendbuf = new long long[soffset];


		// two-phase algorithm
		for (int i = 0; i < basecount; i++){
			for (int it = 0; it < loopcount; it++) {

				double st = MPI_Wtime();
				twophase_tuable_radix_alltoallv(bases[i], (char*)send_buffer, sendcounts, sdispls, MPI_UNSIGNED_LONG_LONG, (char*)recv_buffer, recvcounts, rdispls, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
				double et = MPI_Wtime();
				double total_time = et - st;

				// check correctness
				for (int i=0; i < roffset; i++) {
					if ( (recv_buffer[i] % 10) != (rank % 10) )
						std::cout << "TwoPhase EROOR: " << rank << " " << i << " " << recv_buffer[i] << std::endl;
				}

				if (warmup == 0) {
					double max_time = 0;
					MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
					if (total_time == max_time)
						std::cout << "[TwoPhase] " << nprocs << " " << n << " " <<  max_time << std::endl;
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);


		delete[] send_buffer;
		delete[] recv_buffer;
//		delete[] sendbuf;

	}
}

/*
 * pingpongcheck.cpp
 *
 *  Created on: Jul 28, 2022
 *      Author: kokofan
 */

#include "../rbrucks.h"

static int rank, nprocs;

void run(int ite_count, int warmup);

int main(int argc, char **argv) {
    // MPI Initial
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
        std::cout << "ERROR: MPI_Init error\n" << std::endl;
    if (MPI_Comm_size(MPI_COMM_WORLD, &nprocs) != MPI_SUCCESS)
    	std::cout << "ERROR: MPI_Comm_size error\n" << std::endl;
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    	std::cout << "ERROR: MPI_Comm_rank error\n" << std::endl;

    run(30, 1); // warm-up only
    MPI_Barrier(MPI_COMM_WORLD);

    run(100, 0);

	MPI_Finalize();
    return 0;
}

void run(int ite_count, int warmup) {

    for (int count = 4; count <= 2048; count *= 2) {
		int send_data[count];
		std::fill_n(send_data, count, rank);
		int receve_data[count];
		std::fill_n(receve_data, count, -1);

    	for (int i = 0; i < ite_count; i++) {
			double Start = MPI_Wtime();
			if (rank == (nprocs-1))
				MPI_Recv(&receve_data, count, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (rank == 0)
				MPI_Send(&send_data, count, MPI_INT, (nprocs-1), 0, MPI_COMM_WORLD);
			double end = MPI_Wtime();
			double time = end - Start;

			if (warmup == 0) {
				double max_time = 0;
				MPI_Allreduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				if (time == max_time)
					std::cout << nprocs << " " << count << " " << max_time << std::endl;
			}
    	}
    }
}

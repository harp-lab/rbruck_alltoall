/*
 * generate_ncomm_ndb.cpp
 *
 *  Created on: Sep 1, 2022
 *      Author: kokofan
 */

#include "radix_r_bruck.h"

int main(int argc, char **argv) {

    if (argc < 2) {
    	std::cout << "Usage: mpirun -n <nprocs> " << argv[0] << " <baselist>" << std::endl;
    	return -1;
    }

    int nprocs = atoi(argv[1]);

	int ebr = ceil(sqrt(nprocs));
	std::cout << "Estimated Best R of " << nprocs << " is " << ebr << std::endl;

    for ( int r = 2; r < nprocs; r++) {

    	int total_dc = 0;
    	std::vector<int> the_sd_pstep;

    	calculate_commsteps_and_datablock_counts(nprocs, r, the_sd_pstep);

//		std::cout << nprocs << " " << r << " " <<the_sd_pstep.size() << " [";
		for (int s = 0; s < the_sd_pstep.size(); s++) {
//			std::cout << the_sd_pstep[s] << " ";
			total_dc += the_sd_pstep[s];
		}
//		std::cout << "] " << total_dc << std::endl;

    	std::cout << r << " " << the_sd_pstep.size() << " " << total_dc << std::endl;
    }

    return 0;
}


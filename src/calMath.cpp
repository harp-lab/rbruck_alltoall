/*
 * calMath.cpp
 *
 *  Created on: Sep 1, 2022
 *      Author: kokofan
 */

#include "radix_r_bruck.h"

void calculate_commsteps_and_datablock_counts(int nprocs, int r, std::vector<int>& the_sd_pstep) {

	int w = ceil(log(nprocs) / log(r)); // calculate the number of digits when using r-representation
	int nlpow = myPow(r, w-1);
	int d = (myPow(r, w) - nprocs) / nlpow; // calculate the number of highest digits

	for (int x = 0; x < w; x++) {
		int ze = (x == w - 1)? r - d: r;
		for (int z = 1; z < ze; z++) {

			int xhpow = myPow(r, x+1);
			int xpow = myPow(r, x);
			int div = floor(nprocs / xhpow );
			int re = nprocs % xhpow;
			int t  = re - z * xpow;

			int dc = div * xpow; // number of sent data-blocks per step
			if (t > 0) {
				if (t / xpow > 0) { dc += xpow; }
				else { dc += t % xpow; }
			}
			the_sd_pstep.push_back(dc);
		}
	}
}

/*
 * pingpong.cpp
 *
 *  Created on: Jun 10, 2016
 *      Author: Robert Sawko
 */

#include "../mui/mui.h"

int main() {
    using namespace mui;

    uniface1d interface("mpi://ping/ifs");
    int state = 0;
    for (int t = 0; t < 100; ++t) {
        state++;
        interface.push("data", 0, state);
        interface.commit(t);
        state = interface.fetch("data", 0, t, sampler_exact1d<int>(),
                                chrono_sampler_exact1d());

    }
    printf("Final ping state: %d\n", state);

    return 0;
}

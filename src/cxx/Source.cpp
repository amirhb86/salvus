//
// Created by Michael Afanasiev on 2016-02-04.
//

#include <iostream>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <vector>
#include "Source.h"

std::vector<Source *> Source::factory(Options options) {

    std::vector<Source*> sources;
    for (auto i = 0; i < options.NumberSources(); i++)
    try {
        if (options.SourceType() == "ricker") {
            sources.push_back(new Ricker(options, i));
        } else {
            throw std::runtime_error("Runtime error: Source type " + options.SourceType() + " not supported.");
        }

    } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        MPI::COMM_WORLD.Abort(-1);
    }

    return sources;
}

Ricker::Ricker(Options options, int number) {

    SetLocationX(options.SourceLocationX()[number]);
    SetLocationY(options.SourceLocationY()[number]);
    SetLocationZ(options.SourceLocationZ()[number]);

    mTimeDelay = options.SourceRickerTimeDelay()[number];
    mAmplitude = options.SourceRickerAmplitude()[number];
    mCenterFreq = options.SourceRickerCenterFreq()[number];
}

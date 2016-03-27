//
// Created by Michael Afanasiev on 2016-03-27.
//

#define CATCH_CONFIG_RUNNER
#include "catch.h"
#include <petsc.h>
#include <Element/Element.h>
#include <Element/HyperCube/Quad/Acoustic.h>

int main(int argc, char *argv[]) {

    PetscInitialize(&argc, &argv, NULL, NULL);

    argc = 1;
//    argv = [];
    int result = Catch::Session().run(argc, argv);

    PetscFinalize();

    return result;
}

TEST_CASE("Test whether simple stuff works.", "[element]") {

    Options options;
    options.setOptions();
    std::cout << options.ElementShape() << std::endl;

}


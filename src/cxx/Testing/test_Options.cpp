#include <iostream>
#include <salvus.h>
#include "catch.h"


TEST_CASE("Unit test options", "[options]") {

  SECTION("Homogeneous Dirichlet") {
    PetscOptionsClear(NULL);
    const char *arg[] = {
        "salvus_test",
        "--testing", "true",
        "--homogeneous-dirichlet", "x0,x1",
        NULL};

    /* Fake setting via command line. */
    char **argv = const_cast<char **> (arg);
    int argc = sizeof(arg) / sizeof(const char *) - 1;
    PetscOptionsInsert(NULL, &argc, &argv, NULL);

    std::unique_ptr<Options> options(new Options);
    options->setOptions();
    std::vector<std::string> bnds {"x0", "x1"};
    REQUIRE(options->HomogeneousDirichlet() == bnds);

  }

}

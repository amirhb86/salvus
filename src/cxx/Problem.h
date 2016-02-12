//
// Created by Michael Afanasiev on 2016-01-27.
//

#ifndef SALVUS_SOLVER_H
#define SALVUS_SOLVER_H


#include <iosfwd>
#include <string>

class Problem {

public:

    static Problem *factory(std::string solver_type);

    virtual void initialize() = 0;
    virtual void solve() = 0;


};

class TimeDomain : public Problem {

public:

    virtual void solve();
    virtual void initialize();

};

class FrequencyDomain : public Problem {

public:

    virtual void solve();
    virtual void initialize();

};


#endif //SALVUS_SOLVE

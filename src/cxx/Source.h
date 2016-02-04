//
// Created by Michael Afanasiev on 2016-02-04.
//

#ifndef SALVUS_SOURCE_H
#define SALVUS_SOURCE_H


#include "Options.h"

class Source {

    double mLocationX;
    double mLocationY;
    double mLocationZ;

public:

    static std::vector<Source*> factory(Options options);

    inline void SetLocationX(double location_x) { mLocationX = location_x; }
    inline void SetLocationY(double location_y) { mLocationY = location_y; }
    inline void SetLocationZ(double location_z) { mLocationZ = location_z; }

    inline double LocationX() { return mLocationX; }
    inline double LocationY() { return mLocationY; }
    inline double LocationZ() { return mLocationZ; }

};

class Ricker: public Source {

    double mAmplitude;
    double mTimeDelay;
    double mCenterFreq;

public:

    Ricker(Options options, int number);

};


#endif //SALVUS_SOURCE_H

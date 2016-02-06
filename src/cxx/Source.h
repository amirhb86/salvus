//
// Created by Michael Afanasiev on 2016-02-04.
//

#ifndef SALVUS_SOURCE_H
#define SALVUS_SOURCE_H


#include "Options.h"

class Source {

    double mPhysicalLocationX;
    double mPhysicalLocationY;
    double mPhysicalLocationZ;

    double mReferenceLocationEps;
    double mReferenceLocationYname;
    double mReferenceLocationEta;

public:

    static std::vector<Source*> factory(Options options);

    inline void SetPhysicalLocationX(double location_x) { mPhysicalLocationX = location_x; }
    inline void SetPhysicalLocationY(double location_y) { mPhysicalLocationY = location_y; }
    inline void SetPhysicalLocationZ(double location_z) { mPhysicalLocationZ = location_z; }

    inline void setReferenceLocationEps(double location_eps) { mReferenceLocationEps = location_eps; }
    inline void setReferenceLocationYname(double location_yname) { mReferenceLocationYname = location_yname; }
    inline void setReferenceLocationEta(double location_eta) { mReferenceLocationEta = location_eta; }

    inline double PhysicalLocationX() { return mPhysicalLocationX; }
    inline double PhysicalLocationY() { return mPhysicalLocationY; }
    inline double PhysicalLocationZ() { return mPhysicalLocationZ; }

    inline double ReferenceLocationEps() { return mReferenceLocationEps; }
    inline double ReferenceLocationYname() { return mReferenceLocationYname; }
    inline double ReferenceLocationEta() { return mReferenceLocationEta; }

    virtual double fire(const double &time) = 0;

};

class Ricker: public Source {

    double mAmplitude;
    double mTimeDelay;
    double mCenterFreq;

public:

    Ricker(Options options, int number);
    double fire(const double &time);

};


#endif //SALVUS_SOURCE_H

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

    /**
     * Set the source x coordinate in physical (model) space.
     */
    inline void SetPhysicalLocationX(double location_x) { mPhysicalLocationX = location_x; }

    /**
     * Set the source y coordinate in physical (model) space.
     */
    inline void SetPhysicalLocationY(double location_y) { mPhysicalLocationY = location_y; }

    /**
     * Set the source z coordinate in physical (model) space.
     */
    inline void SetPhysicalLocationZ(double location_z) { mPhysicalLocationZ = location_z; }

    /**
     * Set the source epsilon coordinate on the reference element.
     */
    inline void setReferenceLocationEps(double location_eps) { mReferenceLocationEps = location_eps; }

    /**
     * Set the source MAYBE WE SHOULD CHANGE COORDINATE NAMES coordinate on the reference element.
     */
    inline void setReferenceLocationYname(double location_yname) { mReferenceLocationYname = location_yname; }

    /**
     * Set the source eta coordinate on the reference element.
     */
    inline void setReferenceLocationEta(double location_eta) { mReferenceLocationEta = location_eta; }

    inline double PhysicalLocationX() { return mPhysicalLocationX; }
    inline double PhysicalLocationY() { return mPhysicalLocationY; }
    inline double PhysicalLocationZ() { return mPhysicalLocationZ; }

    inline double ReferenceLocationEps() { return mReferenceLocationEps; }
    inline double ReferenceLocationYname() { return mReferenceLocationYname; }
    inline double ReferenceLocationEta() { return mReferenceLocationEta; }

    /**
     * Returns a value for the force, given a certain time. This needs to be implemented by each derived class. For
     * instance, a Ricker source will need to implement the source time characeristics of a ricker source time
     * function.
     */
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

//
// Created by Michael Afanasiev on 2016-01-29.
//

#ifndef SALVUS_OPTIONS_H
#define SALVUS_OPTIONS_H

#include <iosfwd>
#include <string>
#include <vector>
#include "petscsys.h"

class Options {

    // Integer options.
    PetscInt mNumberSources;
    PetscInt mPolynomialOrder;

    // String options.
    std::string mMeshType;
    std::string mExodusMeshFile;
    std::string mElementShape;
    std::string mPhysicsSystem;
    std::string mExodusModelFile;
    std::string mSourceType;
    // Determines the output name of the movie.
    std::string mOutputMovieFile;

    std::vector<double> mSourceLocationX;
    std::vector<double> mSourceLocationY;
    std::vector<double> mSourceLocationZ;
    std::vector<double> mSourceRickerAmplitude;
    std::vector<double> mSourceRickerCenterFreq;
    std::vector<double> mSourceRickerTimeDelay;

public:

    void setOptions();

    // Integer geters
    inline PetscInt PolynomialOrder() { return mPolynomialOrder; }
    inline PetscInt NumberSources() { return mNumberSources; }

    // String geters
    inline std::string PhysicsSystem() { return mPhysicsSystem; }
    inline std::string ExodusMeshFile() { return mExodusMeshFile; }
    inline std::string MeshType() { return mMeshType; }
    inline std::string ElementShape() { return mElementShape; }
    inline std::string ExodusModelFile() { return mExodusModelFile; }
    inline std::string SourceType() { return mSourceType; }
    inline std::string OutputMovieFile() { return mOutputMovieFile; }

    // Vector getattrs.
    inline std::vector<double> SourceLocationX() { return mSourceLocationX; }
    inline std::vector<double> SourceLocationY() { return mSourceLocationY; }
    inline std::vector<double> SourceLocationZ() { return mSourceLocationZ; }
    inline std::vector<double> SourceRickerAmplitude() { return mSourceRickerAmplitude; }
    inline std::vector<double> SourceRickerCenterFreq() { return mSourceRickerCenterFreq; }
    inline std::vector<double> SourceRickerTimeDelay() { return mSourceRickerTimeDelay; }


};


#endif //SALVUS_OPTIONS_H_H

//
// Created by Michael Afanasiev on 2016-01-29.
//

#ifndef SALVUS_OPTIONS_H
#define SALVUS_OPTIONS_H

#include <iosfwd>
#include <string>
#include <vector>
#include <petscsys.h>

class Options {

    // Integer options.
    PetscInt mDimension;
    PetscInt mNumberSources;
    PetscInt mPolynomialOrder;

    // Double options.
    // Simulation duration
    double mDuration;
    // Time step.
    double mTimeStep;

    // String options.
    std::string mMeshType;
    std::string mExodusMeshFile;
    std::string mElementShape;
    std::string mPhysicsSystem;
    std::string mExodusModelFile;
    std::string mSourceType;
    std::string mProblemType;
    std::string mTimeStepType;

    std::vector<std::string> mDirchletBoundaryNames;
    
    // Determines the output name of the movie.
    std::string mOutputMovieFile;
    // Save movie?, and how often.
    PetscBool mSaveMovie;
    PetscInt mSaveFrameEvery;

    // Run initial condition test with exact solution
    PetscBool mTestIC;
    // exact solution parameters
    PetscReal mCenter_x;
    PetscReal mCenter_z;
    PetscReal mSquareSide_L;

    std::vector<double> mSourceLocationX;
    std::vector<double> mSourceLocationY;
    std::vector<double> mSourceLocationZ;
    std::vector<double> mSourceRickerAmplitude;
    std::vector<double> mSourceRickerCenterFreq;
    std::vector<double> mSourceRickerTimeDelay;

public:

    PetscErrorCode setOptions();

    // Bool getters
    inline PetscBool SaveMovie() { return mSaveMovie; }
    inline PetscBool TestIC() { return mTestIC; }    
    
    // Integer getters
    inline PetscInt PolynomialOrder() { return mPolynomialOrder; }
    inline PetscInt NumberSources() { return mNumberSources; }
    inline PetscInt SaveFrameEvery() { return mSaveFrameEvery; }
    
    // Double getters
    inline double Duration() { return mDuration; }
    inline double TimeStep() { return mTimeStep; }
    inline PetscReal IC_Center_x() { return mCenter_x; }
    inline PetscReal IC_Center_z() { return mCenter_z; }
    inline PetscReal IC_SquareSide_L() { return mSquareSide_L; }

    // String getters
    inline std::string PhysicsSystem() { return mPhysicsSystem; }
    inline std::string ExodusMeshFile() { return mExodusMeshFile; }
    inline std::string MeshType() { return mMeshType; }
    inline std::string ElementShape() { return mElementShape; }
    inline std::string ExodusModelFile() { return mExodusModelFile; }
    inline std::string SourceType() { return mSourceType; }
    inline std::string OutputMovieFile() { return mOutputMovieFile; }
    inline std::string ProblemType() { return mProblemType; }

    inline std::vector<std::string> DirichletBoundaries() { return mDirchletBoundaryNames; }
    
    // Vector getters
    inline std::vector<double> SourceLocationX() { return mSourceLocationX; }
    inline std::vector<double> SourceLocationY() { return mSourceLocationY; }
    inline std::vector<double> SourceLocationZ() { return mSourceLocationZ; }
    inline std::vector<double> SourceRickerAmplitude() { return mSourceRickerAmplitude; }
    inline std::vector<double> SourceRickerCenterFreq() { return mSourceRickerCenterFreq; }
    inline std::vector<double> SourceRickerTimeDelay() { return mSourceRickerTimeDelay; }
    
};


#endif //SALVUS_OPTIONS_H_H

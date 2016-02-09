//
// Created by Michael Afanasiev on 2016-01-29.
//

#ifndef SALVUS_ELEMENT_H
#define SALVUS_ELEMENT_H


#include <petscsys.h>
#include <petscvec.h>
#include <petscdmtypes.h>
#include <vector>
#include <Eigen/Dense>
#include "../Options.h"
#include "../Model/ExodusModel.h"
#include "../Source.h"
#include "../Mesh.h"

class Element {

private:

    bool mContainsSource;
    std::string mElementShape;
    std::string mPhysicsSystem;
    PetscInt mNumberVertex;
    PetscInt mLocalElementNumber;

protected:

    double mTime;
    PetscInt mPolynomialOrder;
    PetscInt mNumberDimensions;
    PetscInt mNumberDofVertex;
    PetscInt mNumberDofEdge;
    PetscInt mNumberDofFace;
    PetscInt mNumberDofVolume;
    PetscInt mNumberIntegrationPoints;

    Eigen::VectorXi mClosureMapping;
    std::vector<Source*> mSources;

    DM mDistributedMesh;
    PetscSection mMeshSection;

public:

    static Eigen::VectorXd GllPointsForOrder(const int order);
    static Eigen::VectorXd GllIntegrationWeightForOrder(const int order);
    static Eigen::VectorXi ClosureMapping(const int order, const int dimension);

    // Utility methods.
    static Element *factory(Options options);
    virtual Element* clone() const = 0;
    void printInfoToScreen() const;

    // Generic methods.
    void registerMesh(DM &distributed_mesh, PetscSection &section);

    // Pure virtual methods.
    virtual void attachSource(std::vector<Source*> sources) = 0;
    virtual void readOperators() = 0;
    virtual void registerFieldVectors(Mesh *mesh) = 0;
    virtual void attachVertexCoordinates() = 0;
    virtual void attachIntegrationPoints() = 0;
    virtual void constructStiffnessMatrix(Mesh *mesh) = 0;
    virtual void interpolateMaterialProperties(ExodusModel &model) = 0;

    // Integer setattr.
    inline void SetLocalElementNumber(PetscInt number) { mLocalElementNumber = number; }
    inline void SetContainsSource(const bool contains_source) { mContainsSource = contains_source; }
    inline void SetElementShape(const std::string element_shape) { mElementShape = element_shape; }
    inline void SetPhysicsSystem(const std::string physics_system) { mPhysicsSystem = physics_system; }
    inline void SetPolynomialOrder(const int polynomial_order) { mPolynomialOrder = polynomial_order; }
    inline void SetNumberVertex(const int number_vertex) { mNumberVertex = number_vertex; }
    inline void SetNumberDimensions(const int number_dimensions) { mNumberDimensions = number_dimensions; }
    inline void SetCurrentTime(const double time) { mTime = time; }

    // Integer geters.
    inline bool ContainsSource() const { return mContainsSource; }
    inline PetscInt NumberVertex() const { return mNumberVertex; }
    inline PetscInt NumberDimensions() const { return mNumberDimensions; }
    inline PetscInt NumberDofVertex() const { return mNumberDofVertex; }
    inline PetscInt NumberDofEdge() const { return mNumberDofEdge; }
    inline PetscInt NumberDofFace() const { return mNumberDofFace; }
    inline PetscInt NumberDofVolume() const { return mNumberDofVolume; }
    inline PetscInt LocalElementNumber() const { return mLocalElementNumber; }
    inline PetscInt NumberIntegrationPoints() const { return mNumberIntegrationPoints; }


};


#endif //SALVUS_ELEMENT_H

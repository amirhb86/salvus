//
// Created by Michael Afanasiev on 2016-01-27.
//

#ifndef SALVUS_MESH_H
#define SALVUS_MESH_H


#include <iosfwd>
#include <string>
#include <petscdmtypes.h>
#include <vector>
#include <petscistypes.h>
#include <petscvec.h>
#include <Eigen/Dense>
#include "Options.h"

struct vec_struct {

    bool check_in, check_out;
    std::string name;
    Vec field_globals, field_locals;

};

class Mesh {

    PetscInt mNumberElementsLocal;
    PetscInt mNumberDimensions;

    std::string mExodusFileName;

    DM mDistributedMesh;
    PetscSection mMeshSection;

protected:

    Vec mMassMatrix;
    std::map<std::string,vec_struct> mFields;
    std::vector<Vec> mFieldVectorGlobals;
    std::vector<Vec> mFieldVectorLocals;
    std::vector<bool> mFieldVectorCheckin;
    std::vector<bool> mFieldVectorCheckout;
    std::vector<std::string> mFieldVectorNames;
    std::vector<std::string> mCheckedOutFields;

public:

    static Mesh *factory(Options options);

    void read(Options options);
    void setupGlobalDof(PetscInt number_dof_vertex, PetscInt number_dof_edge,
                        PetscInt number_dof_face, PetscInt number_dof_volume,
                        PetscInt number_dimensions);

    void registerFieldVector(const std::string &name, const bool &check_out, const bool &check_in);

    void checkOutFields();
    void checkInFieldsEnd();
    void checkInFieldsBegin();

    void getFieldOnElement(Eigen::VectorXd &field, Eigen::VectorXi &closure, const std::string &name,
                           const int &element_number);
    void setFieldOnElement(const Eigen::VectorXd &field, Eigen::VectorXi &closure, const std::string &name,
                           const int &element_number);


    // Integer getattr.
    inline PetscInt NumberElementsLocal() { return mNumberElementsLocal; }

    virtual void advanceField() = 0;

    // Distributed mesh getattr.
    inline DM &DistributedMesh() { return mDistributedMesh; }
    inline PetscSection &MeshSection() { return mMeshSection; }

};

class ScalarNewmark : public Mesh {

    Vec acceleration_;

public:

    virtual void advanceField();

};

#endif //SALVUS_MESH_H

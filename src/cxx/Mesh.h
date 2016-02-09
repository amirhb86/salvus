//
// Created by Michael Afanasiev on 2016-01-27.
//

#ifndef SALVUS_MESH_H
#define SALVUS_MESH_H

#include <map>
#include <iosfwd>
#include <string>
#include <petscdmtypes.h>
#include <vector>
#include <petscistypes.h>
#include <petscvec.h>
#include <Eigen/Dense>
#include "Options.h"
#include "Utilities.h"

struct vec_struct {

    bool check_in, check_out;
    std::string name;
    Vec field_globals, field_locals;

};

class Mesh {

    int mNumberElementsLocal;
    int mNumberDimensions;

    std::string mExodusFileName;

    DM mDistributedMesh;
    PetscSection mMeshSection;

protected:

    Vec mMassMatrix;
    std::map<int,vec_struct> mFields;
    std::vector<int> mFieldVectorCheckin;
    std::vector<int> mFieldVectorCheckout;

public:

    static Mesh *factory(Options options);

    void read(Options options);
    void setupGlobalDof(PetscInt number_dof_vertex, PetscInt number_dof_edge,
                        PetscInt number_dof_face, PetscInt number_dof_volume,
                        PetscInt number_dimensions);

    void registerFieldVectors(const int &num, const bool &check_out, const bool &check_in,
                              const std::string &name);

    void checkOutFields();
    void checkInFieldsEnd();
    void checkInFieldsBegin();

    Eigen::VectorXd getFieldOnElement(const int &field_num, const int &element_number,
                                            const Eigen::VectorXi &closure);
    void setFieldOnElement(const int &field_num, const int &element_number,
                           const Eigen::VectorXi &closure, const Eigen::VectorXd &field);


    // Integer getattr.
    inline PetscInt NumberElementsLocal() { return mNumberElementsLocal; }

    virtual void advanceField() = 0;
    virtual void registerFields() = 0;

    // Distributed mesh getattr.
    inline DM &DistributedMesh() { return mDistributedMesh; }
    inline PetscSection &MeshSection() { return mMeshSection; }

};

class ScalarNewmark : public Mesh {

    Vec acceleration_;

public:

    virtual void registerFields() {

        registerFieldVectors((int) AcousticFields::acceleration_, false, false, "acceleration_");
        registerFieldVectors((int) AcousticFields::acceleration, false, false, "acceleration");
        registerFieldVectors((int) AcousticFields::displacement, true, false, "displacement");
        registerFieldVectors((int) AcousticFields::velocity, false, false, "velocity");
        registerFieldVectors((int) AcousticFields::force, false, true, "force");

    }

    virtual void advanceField();

};

#endif //SALVUS_MESH_H

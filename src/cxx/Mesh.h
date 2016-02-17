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

    std::string name;
    Vec glb, loc;

};

class Mesh {

    int mNumberElementsLocal;
    int mNumberDimensions;
    int mTime;

    std::string mExodusFileName;

    DM mDistributedMesh;
    PetscSection mMeshSection;

protected:

    Vec mMassMatrix;
    std::map<std::string,vec_struct> mFields;

    PetscViewer mViewer;

public:

    static Mesh *factory(Options options);

    void read(Options options);
    void setupGlobalDof(PetscInt number_dof_vertex, PetscInt number_dof_edge,
                        PetscInt number_dof_face, PetscInt number_dof_volume,
                        PetscInt number_dimensions);

    /**
     * Registers both the global (across parallel partition) and local (on a single parallel partition) vectors for a
     * given name. Vector information is stored in an `std::map` under mFields, with `name` as the key.
     * @param name Name of field to register.
     */
    void registerFieldVectors(const std::string &name);

    /**
     * Begins (and ends) the gloabl -> local MPI sends for a given field name. By the time this function returns, you
     * can be confident that the MPI sends have been completed (i.e. it is blocking), and will be available to each
     * element.
     * @param name Name of field to checkout.
     */
    void checkOutField(const std::string &name);

    /**
     * Begins the local -> global MPI sends for a given field name. The send is performed with an implied sum, i.e.
     * the value of field on coincident GLL points are properly summed together. Note that this function is
     * non-blocking!! This MUST be paired with an equivalent call to checkInFieldEnd.
     * @param name Name of field to checkin.
     */
    void checkInFieldBegin(const std::string &name);

    /**
     * Ends the local -> global MPI send for a given field name. This function should come after an equivalent
     * checkInFieldBegin. This method is blocking, so you can be confident that when it returns the desired field has
     * been scattered and summed into the global (parallel) degrees of freedom.
     * @param name Name of field to checkin.
     */
    void checkInFieldEnd(const std::string &name);

    /**
     * Returns an ordered vector of a field (i.e. x-displacement) on an element, via a call to DMPlexVecGetClosure.
     * Note that a vector containing the closure mapping must also be passed in -- this should change in the future.
     * @param [in] name Name of field.
     * @param [in] element_number Element number (on the local processor) for which the field is requested.
     * @param [in] closure A vector of size mNumberIntegrationPoints, which specifies the mapping from the Plex element
     * closure to the desired gll point ordering.
     * @ return The ordered field on an element.
     */
    Eigen::VectorXd getFieldOnElement(const std::string &name, const int &element_number,
                                      const Eigen::VectorXi &closure);

    /**
     * Sums a field on an element into the degrees of freedom owned by the local processor, via a call to
     * DMPlexVecSetClosure.
     * @param [in] field_name Name of field.
     * @param [in] element_number Element number (on the local processor) for which the field is requested.
     * @param [in] closure A vector of size mNumberIntegrationPoints, which specifies the mapping from the Plex element
     * closure to the desired gll point ordering.
     * @param [in] field The element-ordered field (i.e. x-displacement) to sum into the mesh.
     */
    void setFieldOnElement(const std::string &name, const int &element_number,
                           const Eigen::VectorXi &closure, const Eigen::VectorXd &field);


    // Integer getattr.
    inline PetscInt NumberElementsLocal() { return mNumberElementsLocal; }

    void setUpMovie(const std::string &movie_filename);
    void saveFrame();
    void finalizeMovie();
    void zeroFields(const std::string &name);
    virtual void advanceField() = 0;
    virtual void registerFields() = 0;
    virtual void applyInverseMassMatrix() = 0;

    // Distributed mesh getattr.
    inline DM &DistributedMesh() { return mDistributedMesh; }
    inline PetscSection &MeshSection() { return mMeshSection; }

};

class ScalarNewmark : public Mesh {

public:

    virtual void registerFields() {
        registerFieldVectors("acceleration_");
        registerFieldVectors("acceleration");
        registerFieldVectors("displacement");
        registerFieldVectors("velocity");
        registerFieldVectors("force");
        registerFieldVectors("mass_matrix");
    }

    virtual void advanceField();
    virtual void applyInverseMassMatrix();

};

class ElasticNewmark2D: public Mesh {

public:

    virtual void registerFields() {
        registerFieldVectors("acceleration_x_");
        registerFieldVectors("acceleration_z_");
        registerFieldVectors("acceleration_x");
        registerFieldVectors("acceleration_z");
        registerFieldVectors("displacement_x");
        registerFieldVectors("displacement_z");
        registerFieldVectors("velocity_x");
        registerFieldVectors("velocity_z");
        registerFieldVectors("force_x");
        registerFieldVectors("force_z");
        registerFieldVectors("mass_matrix");
    }

    virtual void advanceField();
    virtual void applyInverseMassMatrix();


};

#endif //SALVUS_MESH_H

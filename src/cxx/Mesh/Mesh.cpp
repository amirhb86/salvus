//
// Created by Michael Afanasiev on 2016-01-27.
//

#include "Mesh.h"
#include "ScalarNewmark2D.h"
#include "ElasticNewmark2D.h"
#include "../Utilities/Utilities.h"

void exodusError(const int retval, std::string func_name) {

    if (retval) {
        throw std::runtime_error("Error in exodus function: " + func_name);
    }

}

extern "C" {
#include <exodusII.h>
}

Mesh *Mesh::factory(Options options) {

    std::string mesh_type(options.MeshType());
    try {
        if (mesh_type == "newmark") {
            return new ScalarNewmark2D;
        } else if (mesh_type == "newmark_2d_elastic"){
            return new ElasticNewmark2D;
        } else {
            throw std::runtime_error("Runtime Error: Mesh type " + mesh_type + " not supported");
        }
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI::COMM_WORLD.Abort(-1);
        return nullptr;
    };

}

#undef __FUNCT__
#define __FUNCT__ "readBoundaryNames"
int Mesh::readBoundaryNames(Options options) {
    
    PetscFunctionBegin;
    std::vector<std::string> boundary_names;
    if (MPI::COMM_WORLD.Get_rank() == 0) {
     
        int num_boundaries = -1;
        float exodus_version;
        int io_ws = 0;
        int comp_ws = 8;
        try {
            int exoid = ex_open(options.ExodusMeshFile().c_str() , EX_READ, &comp_ws, &io_ws, &exodus_version);
            if (exoid < 0) { throw std::runtime_error("Error opening exodus model file."); }
            float fdum = 0.0;
            char* cdum = 0;
            exodusError(ex_inquire(exoid,EX_INQ_SIDE_SETS,&num_boundaries,&fdum,cdum),"ex_inq");
            for(int b=1;b<=num_boundaries;b++) {
                char name[1024];
                int ret=-1;
                // retrieve name from id={1,2,3..}
                // note that blank spaces given in Trelis get an underscore.
                exodusError(ex_get_name(exoid,EX_SIDE_SET,b,name),"ex_get_name");
                std::string namestr(name);
                boundary_names.push_back(name);
            }
        } catch (std::exception &e) {
            std::cerr << "ERROR(rank=" << MPI::COMM_WORLD.Get_rank() << ")!: " << e.what();
            MPI_Abort(PETSC_COMM_WORLD, -1);
        }
    } // rank==0        
        
    boundary_names = utilities::broadcastStringVecFromFroot(boundary_names);
    // Build mapping boundary name -> label value (id). `id` starts counting at 1.
    for(int i=0;i<boundary_names.size();i++) {
        mBoundaryIds[boundary_names[i]] = i+1;
    }        
    MPI::COMM_WORLD.Barrier();
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "setupBoundaries"
int Mesh::setupBoundaries(Options options) {
    PetscFunctionBegin;
    
    int ierr = 0;
    
    readBoundaryNames(options);

    int num_sideset_ids=-1;
    
    // PETSc gives all side sets the label name "Face Sets" Each face
    // is then given value, which is the order found in the exodus
    // file; Ex: two boundaries "absorbing_boundaries", "free_surface"
    // faces on absorbing boundary get label value "1", and faces on
    // free surface get value "2".
    DMLabel label;
    ierr = DMPlexGetLabel(mDistributedMesh,"Face Sets",&label);CHKERRQ(ierr);
    // if no side sets are present, proceed with no boundaries (naturally a free surface)
    if(label == NULL) { return 0; }
    ierr = DMLabelGetNumValues(label, &num_sideset_ids);CHKERRQ(ierr);
    if(num_sideset_ids == 0) {return 0;}
    
    IS idIS;
    const PetscInt* ids;
    ierr = DMLabelGetValueIS(label, &idIS);CHKERRQ(ierr);
    ierr = ISGetIndices(idIS, &ids);CHKERRQ(ierr);
    
    for(int i=0;i<num_sideset_ids;i++) {
        
        IS              pointIS;
        const PetscInt *faces;
        PetscInt        numFaces, p;
        
        ierr = DMLabelGetStratumSize(label, ids[i], &numFaces);CHKERRQ(ierr);
        ierr = DMLabelGetStratumIS(label, ids[i], &pointIS);CHKERRQ(ierr);
        ierr = ISGetIndices(pointIS, &faces);CHKERRQ(ierr);

        for(int p=0;p<numFaces;p++) {

            PetscInt support_size;
            const PetscInt* support;
            int face = faces[p];
            DMPlexGetSupportSize(mDistributedMesh,face,&support_size);
            DMPlexGetSupport(mDistributedMesh,face,&support);
            for(int elem=0;elem<support_size;elem++) {                
                mBoundaryElementFaces[ids[i]][support[elem]].push_back(face);
            }
            
        }
        ISRestoreIndices(pointIS, &faces);
        ISDestroy(&pointIS);
    }
    ierr = ISRestoreIndices(idIS, &ids);CHKERRQ(ierr);
    ISDestroy(&idIS);
    
    // visualize structure
    // for(int rank=0;rank<MPI::COMM_WORLD.Get_size();rank++) {
    //     if(rank==MPI::COMM_WORLD.Get_rank()) {
    //         printf("proc=%d\n",MPI::COMM_WORLD.Get_rank());
    //         for(auto& id_key : mBoundaryElementFaces) {
    //             int i=id_key.first;
    //             for(auto& elem : id_key.second) {
    //                 printf("%d:(%d:{",i,elem.first);
    //                 for(int face=0;face<elem.second.size();face++) {                
    //                     printf("%d,",elem.second[face]);
    //                 }
    //                 printf("})\n");
    //             }
    //         }
    //     }
    //     MPI::COMM_WORLD.Barrier();
    // }
    
    // Set some important information.
    
    PetscFunctionReturn(0);
}

void Mesh::read(Options options) {

    // Class variables.
    mDistributedMesh = NULL;
    mExodusFileName = options.ExodusMeshFile();

    // Function variables.
    DM dm = NULL;
    PetscInt partition_overlap = 0;
    PetscBool interpolate_edges = PETSC_TRUE;

    // Read exodus file.
    DMPlexCreateExodusFromFile(PETSC_COMM_WORLD, mExodusFileName.c_str(), interpolate_edges, &dm);

    // May be a race condition on distribute.
    MPI::COMM_WORLD.Barrier();
    DMPlexDistribute(dm, partition_overlap, NULL, &mDistributedMesh);

    // We don't need the serial mesh anymore.
    if (mDistributedMesh) { DMDestroy(&dm); }

    DMGetDimension(mDistributedMesh, &mNumberDimensions);
    DMPlexGetDepthStratum(mDistributedMesh, mNumberDimensions, NULL, &mNumberElementsLocal);

    setupBoundaries(options);
}

void Mesh::setupGlobalDof(int number_dof_vertex, int number_dof_edge, int number_dof_face,
                          int number_dof_volume, int number_dimensions) {

    // Ensure that the mesh and the elements are the same dimension.
    assert(number_dimensions == mNumberDimensions);

    // Only define 1 field here because we're taking care of multiple fields manually.
    int number_fields = 1;
    int number_components = 1;
    int number_dof_per_element[mNumberDimensions + 1];

    // Number of dof on vertex, edge, face, volume.
    number_dof_per_element[0] = number_dof_vertex;
    number_dof_per_element[1] = number_dof_edge;
    number_dof_per_element[2] = number_dof_face;
    if (mNumberDimensions == 3) { number_dof_per_element[3] = number_dof_volume; }

    // Setup the global and local (distributed) degrees of freedom.
    DMPlexCreateSection(mDistributedMesh, mNumberDimensions, number_fields, &number_components,
                        number_dof_per_element, 0, NULL, NULL, NULL, NULL, &mMeshSection);
    DMSetDefaultSection(mDistributedMesh, mMeshSection);

    // supposed to improve performance of closure calls
    DMPlexCreateClosureIndex(mDistributedMesh, mMeshSection);    
}

void Mesh::registerFieldVectors(const std::string &name) {

    Vec field_vector_local;
    Vec field_vector_global;
    DMCreateLocalVector(mDistributedMesh, &field_vector_local);
    DMCreateGlobalVector(mDistributedMesh, &field_vector_global);

    double zero = 0.0;
    VecSet(field_vector_local, zero);
    VecSet(field_vector_global, zero);
    PetscObjectSetName((PetscObject) field_vector_global, name.c_str());

    vec_struct registrar;
    registrar.name = name;
    registrar.loc = field_vector_local;
    registrar.glb = field_vector_global;
    mFields[name] = registrar;

}

void Mesh::checkOutField(const std::string &name) {

    assert(mFields.find(name) != mFields.end());

    // Begin the MPI broadcast global -> local.
    DMGlobalToLocalBegin(mDistributedMesh, mFields[name].glb, INSERT_VALUES,
                         mFields[name].loc);

    // End the MPI broadcast global -> local.
    DMGlobalToLocalEnd(mDistributedMesh, mFields[name].glb, INSERT_VALUES,
                       mFields[name].loc);

}

Eigen::VectorXd Mesh::getFieldOnFace(const std::string &name, const int &face_number,
                                     const Eigen::VectorXi &closure) {

    PetscScalar *val = NULL;
    Eigen::VectorXd field(closure.size());
    PetscInt num_nodes_face = -1;
    DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                        face_number, &num_nodes_face, &val);    
    for (auto j = 0; j < closure.size(); j++) { field(closure(j)) = val[j]; }
    DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                            face_number, NULL, &val);
    return field;
}

void Mesh::setFieldFromFace(const std::string &name, const int face_number,
                               const Eigen::VectorXi &face_closure, const Eigen::VectorXd &field) {

    Eigen::VectorXd val(face_closure.size());
    // map "our" nodal ordering back onto PETSC ordering
    for (auto j = 0; j < face_closure.size(); j++) { val(j) = field(face_closure(j)); }
    DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                        face_number, val.data(), INSERT_VALUES);
}

void Mesh::addFieldFromFace(const std::string &name, const int face_number,
                            const Eigen::VectorXi &face_closure, const Eigen::VectorXd &field) {

    Eigen::VectorXd val(face_closure.size());
    // map "our" nodal ordering back onto PETSC ordering
    for (auto j = 0; j < face_closure.size(); j++) { val(j) = field(face_closure(j)); }
    DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                        face_number, val.data(), ADD_VALUES);
}

Eigen::VectorXd Mesh::getFieldOnElement(const std::string &name, const int &element_number,
                                        const Eigen::VectorXi &closure) {

    PetscScalar *val = NULL;
    Eigen::VectorXd field(closure.size());
    DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                        element_number, NULL, &val);
    for (auto j = 0; j < closure.size(); j++) { field(closure(j)) = val[j]; }
    DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                            element_number, NULL, &val);
    return field;

}

void Mesh::setFieldFromElement(const std::string &name, const int element_number,
                               const Eigen::VectorXi &closure, const Eigen::VectorXd &field) {

    Eigen::VectorXd val(closure.size());
    // map "our" nodal ordering back onto PETSC ordering
    for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
    DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                        element_number, val.data(), INSERT_VALUES);
}

void Mesh::addFieldFromElement(const std::string &name, const int element_number,
                               const Eigen::VectorXi &closure, const Eigen::VectorXd &field) {

    Eigen::VectorXd val(closure.size());
    // map "our" nodal ordering back onto PETSC ordering
    for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
    DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                        element_number, val.data(), ADD_VALUES);
}

void Mesh::assembleLocalFieldToGlobal(const std::string &name) {

    assembleLocalFieldToGlobalBegin(name);
    assembleLocalFieldToGlobalEnd(name);   
}

void Mesh::assembleLocalFieldToGlobalBegin(const std::string &name) {

    // Make sure the field exists in our dictionary.
    assert(mFields.find(name) != mFields.end());

    // Begin MPI broadcast local -> global.
    DMLocalToGlobalBegin(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
}

void Mesh::assembleLocalFieldToGlobalEnd(const std::string &name) {

    // Make sure the field exists in our dictionary.
    assert(mFields.find(name) != mFields.end());

    // Begin MPI broadcast local -> global.
    DMLocalToGlobalEnd(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
}

void Mesh::setLocalFieldToGlobal(const std::string &name) {
    
    // Make sure the field exists in our dictionary.
    assert(mFields.find(name) != mFields.end());

    // Do "communication". `INSERT_VALUE` will result in no communication
    DMLocalToGlobalBegin(mDistributedMesh, mFields[name].loc, INSERT_VALUES, mFields[name].glb);
    DMLocalToGlobalEnd(mDistributedMesh, mFields[name].loc, INSERT_VALUES, mFields[name].glb);

}

void Mesh::checkInFieldBegin(const std::string &name) {

    // Make sure the field exists in our dictionary.
    assert(mFields.find(name) != mFields.end());

    // Begin MPI broadcast local -> global.
    DMLocalToGlobalBegin(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
}


void Mesh::checkInFieldEnd(const std::string &name) {

    // Make sure the field exists in our dictionary.
    assert(mFields.find(name) != mFields.end());

    // Begin MPI broadcast local -> global.
    DMLocalToGlobalEnd(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);

}

void Mesh::zeroFields(const std::string &name) {
    double zero = 0.0;
    VecSet(mFields[name].loc, zero);
    VecSet(mFields[name].glb, zero);
}

void Mesh::setUpMovie(const std::string &movie_filename) {
    mViewer = nullptr;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, movie_filename.c_str(), FILE_MODE_WRITE, &mViewer);
    PetscViewerHDF5PushGroup(mViewer, "/");
    DMView(mDistributedMesh, mViewer);
}

void Mesh::saveFrame(std::string name,PetscInt timestep) {

    DMSetOutputSequenceNumber(mDistributedMesh, timestep, timestep);
    VecView(mFields[name].glb, mViewer);
}

void Mesh::finalizeMovie() {

    PetscViewerHDF5PopGroup(mViewer);
    PetscViewerDestroy(&mViewer);

}


#include <petscviewerhdf5.h>
#include <stdexcept>

#include <mpi.h>
#include <fstream>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Utilities/Utilities.h>
#include <Utilities/Logging.h>
#include <petscdmplex.h>
#include <Utilities/Types.h>

extern "C" {
#include <exodusII.h>
}

using std::unique_ptr;

Mesh::Mesh(const std::unique_ptr<Options> &options) {
  mExodusFileName = options->MeshFile();
  mDistributedMesh = NULL;
  mMeshSection = NULL;
  mNumDim = 0;
}

std::unique_ptr<Mesh> Mesh::Factory(const std::unique_ptr<Options> &options) {
  return std::unique_ptr<Mesh> (new Mesh(options));
}

void Mesh::read(int dim, int numCells, int numVerts, int numVertsPerElem,
                int* cells, double* vertex_coords) {

  // Class variables.
  mDistributedMesh = NULL;
  DM dm = NULL;
  PetscBool interpolate_edges = PETSC_TRUE;
  
  DMPlexCreateFromCellList(PETSC_COMM_WORLD, dim, numCells, numVerts, numVertsPerElem,
                           interpolate_edges, cells, dim, vertex_coords,&dm);
  
  DMPlexDistribute(dm, 0, NULL, &mDistributedMesh);

  // We don't need the serial mesh anymore.
  if (mDistributedMesh) { DMDestroy(&dm); }
    // mDistributedMesh == NULL when only 1 proc is used.
  else { mDistributedMesh = dm; }

  DMGetDimension(mDistributedMesh, &mNumDim);
  DMPlexGetDepthStratum(mDistributedMesh, mNumDim, NULL, &mNumberElementsLocal);
}

void Mesh::read() {

  // Class variables.
  mDistributedMesh = NULL;

  // check if file exists
  PetscInt rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    std::ifstream f(mExodusFileName.c_str());
    if (!f.good()) {
      throw std::runtime_error("Mesh file not found! You requested '" + mExodusFileName + "'.");
    }
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  // Function variables.
  DM dm = NULL;
  PetscBool interpolate_edges = PETSC_TRUE;

  /* Read exodus file. */
  DMPlexCreateExodusFromFile(PETSC_COMM_WORLD, mExodusFileName.c_str(), interpolate_edges, &dm);

  /* Distribute mesh. */
  DMPlexDistribute(dm, 0, NULL, &mDistributedMesh);

  /* We don't need the serial mesh anymore if we're in parallel. */
  if (mDistributedMesh) { DMDestroy(&dm); }
  else { mDistributedMesh = dm; }

  DMGetDimension(mDistributedMesh, &mNumDim);
  // Get number of elements (duh.)
  DMPlexGetDepthStratum(mDistributedMesh, mNumDim, NULL, &mNumberElementsLocal);

}

#undef __FUNCT__
#define __FUNCT__ "FixOrientation"
PetscErrorCode FixOrientation(DM dm, PetscSection s, int* user_Nc, int user_Nf, int user_k, int user_dim)
{
  PetscInt       f, o, i, j, k, c, d;
  DMLabel        depthLabel;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetLabel(dm,"depth",&depthLabel);CHKERRQ(ierr);
  for (f = 0; f < user_Nf; f++) {
    PetscSectionSym sym;

    if (user_k < 3) continue; /* No symmetries needed for order < 3, because no cell, facet, edge or vertex has more than one node */
    ierr = PetscSectionSymCreateLabel(PetscObjectComm((PetscObject)s),depthLabel,&sym);CHKERRQ(ierr);

    for (d = 0; d <= user_dim; d++) {
      /* Edge. */
      if (d == 1) {
        PetscInt        numDof  = user_k - 1;
        PetscInt        numComp = user_Nc[f];
        PetscInt        minOrnt = -2;
        PetscInt        maxOrnt = 2;
        PetscInt        **perms;

        ierr = PetscCalloc1(maxOrnt - minOrnt,&perms);CHKERRQ(ierr);
        for (o = minOrnt; o < maxOrnt; o++) {
          PetscInt *perm;

          if (o == -1 || !o) { /* identity */
            perms[o - minOrnt] = NULL;
          } else {
            ierr = PetscMalloc1(numDof * numComp, &perm);CHKERRQ(ierr);
            for (i = numDof - 1, k = 0; i >= 0; i--) {
              for (j = 0; j < numComp; j++, k++) perm[k] = i * numComp + j;
            }
            perms[o - minOrnt] = perm;
          }
        }
        ierr = PetscSectionSymLabelSetStratum(sym,d,numDof*numComp,minOrnt,maxOrnt,PETSC_OWN_POINTER,(const PetscInt **) perms,NULL);CHKERRQ(ierr);
        /* Surface */
      } else if (d == 2) {
        PetscInt        perEdge = user_k - 1;
        PetscInt        numDof  = perEdge * perEdge;
        PetscInt        numComp = user_Nc[f];
        PetscInt        minOrnt = -4;
        PetscInt        maxOrnt = 4;
        PetscInt        **perms;

        ierr = PetscCalloc1(maxOrnt-minOrnt,&perms);CHKERRQ(ierr);
        for (o = minOrnt; o < maxOrnt; o++) {
          PetscInt *perm;

          if (!o) continue; /* identity */
          ierr = PetscMalloc1(numDof * numComp, &perm);CHKERRQ(ierr);
          /* We want to perm[k] to list which *localArray* position the *sectionArray* position k should go to for the given orientation*/
          switch (o) {
          case 0:
            break; /* identity */
          case -4: /* flip along (-1,-1)--( 1, 1), which swaps edges 0 and 3 and edges 1 and 2.  This swaps the i and j variables */
            for (i = 0, k = 0; i < perEdge; i++) {
              for (j = 0; j < perEdge; j++, k++) {
                for (c = 0; c < numComp; c++) {
                  perm[k * numComp + c] = (perEdge * j + i) * numComp + c; // old
                }
              }
            }
            break;
          case -3: /* flip along (-1, 0)--( 1, 0), which swaps edges 0 and 2.  This reverses the i variable */
            for (i = 0, k = 0; i < perEdge; i++) {
              for (j = 0; j < perEdge; j++, k++) {
                for (c = 0; c < numComp; c++) {
                  perm[k * numComp + c] = (perEdge * (perEdge - 1 - i) + j) * numComp + c;
                }
              }
            }
            break;
          case -2: /* flip along ( 1,-1)--(-1, 1), which swaps edges 0 and 1 and edges 2 and 3.  This swaps the i and j variables and reverse both */
            for (i = 0, k = 0; i < perEdge; i++) {
              for (j = 0; j < perEdge; j++, k++) {
                for (c = 0; c < numComp; c++) {
                  perm[k * numComp + c] = (perEdge * (perEdge - 1 - j) + (perEdge - 1 - i)) * numComp + c;
                }
              }
            }
            break;
          case -1: /* flip along ( 0,-1)--( 0, 1), which swaps edges 3 and 1.  This reverses the j variable */
            for (i = 0, k = 0; i < perEdge; i++) {
              for (j = 0; j < perEdge; j++, k++) {
                for (c = 0; c < numComp; c++) {
                  perm[k * numComp + c] = (perEdge * i + (perEdge - 1 - j)) * numComp + c;
                }
              }
            }
            break;
          case  1: /* rotate section edge 1 to local edge 0.  This swaps the i and j variables and then reverses the j variable */
            for (i = 0, k = 0; i < perEdge; i++) {
              for (j = 0; j < perEdge; j++, k++) {
                for (c = 0; c < numComp; c++) {
                  perm[k * numComp + c] = (perEdge * (perEdge - 1 - j) + i) * numComp + c;
                }
              }
            }
            break;
          case  2: /* rotate section edge 2 to local edge 0.  This reverse both i and j variables */
            for (i = 0, k = 0; i < perEdge; i++) {
              for (j = 0; j < perEdge; j++, k++) {
                for (c = 0; c < numComp; c++) {
                  perm[k * numComp + c] = (perEdge * (perEdge - 1 - i) + (perEdge - 1 - j)) * numComp + c;
                }
              }
            }
            break;
          case  3: /* rotate section edge 3 to local edge 0.  This swaps the i and j variables and then reverses the i variable */
            for (i = 0, k = 0; i < perEdge; i++) {
              for (j = 0; j < perEdge; j++, k++) {
                for (c = 0; c < numComp; c++) {
                  perm[k * numComp + c] = (perEdge * j + (perEdge - 1 - i)) * numComp + c;
                }
              }
            }
            break;
          default:
            break;
          }
          perms[o - minOrnt] = perm;
        }
        ierr = PetscSectionSymLabelSetStratum(sym,d,numDof*numComp,minOrnt,maxOrnt,PETSC_OWN_POINTER,(const PetscInt **) perms,NULL);CHKERRQ(ierr);
      }
    }
    ierr = PetscSectionSetFieldSym(s,f,sym);CHKERRQ(ierr);
    ierr = PetscSectionSymDestroy(&sym);CHKERRQ(ierr);
  }
  ierr = PetscSectionViewFromOptions(s,NULL,"-section_with_sym_view");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void Mesh::setupTopology(const unique_ptr<ExodusModel> &model,
                         const unique_ptr<Options> &options) {

  /* Ensure mesh was read. */
  if (!mNumDim) { throw std::runtime_error("Mesh appears to have zero dimensions. Have you called read()?"); }

  /* Find all the mesh boundaries. */
  DMLabel label; DMGetLabel(mDistributedMesh, "Face Sets", &label);
  PetscInt boundary_size = 0;
  // label==NULL when there are no side sets, so there are no labeled boundaries
  if (label) {
    DMLabelGetNumValues(label, &boundary_size);
    IS idIS; DMLabelGetValueIS(label, &idIS);
    const PetscInt *ids; ISGetIndices(idIS, &ids);
    for (PetscInt i = 0; i < boundary_size; i++) {
      PetscInt numFaces; DMLabelGetStratumSize(label, ids[i], &numFaces);
      IS pointIs; DMLabelGetStratumIS(label, ids[i], &pointIs);
      const PetscInt *faces; ISGetIndices(pointIs, &faces);
      /* Tuple describing boundary set i for face j. */
      for (PetscInt j = 0; j < numFaces; j++) {
          mBndPts.insert(std::tie(i, faces[j]));
        /* We need to walk down the graph and apply boundaries to all points. */
        PetscInt num_closure, *val_closure = NULL;
        DMPlexGetTransitiveClosure(mDistributedMesh, faces[j], PETSC_TRUE,
                                   &num_closure, &val_closure);
        for (PetscInt k = 0; k < 2 * num_closure; k += 2) {
          mBndPts.insert(std::tie(i, val_closure[k]));
        }
        DMPlexRestoreTransitiveClosure(mDistributedMesh, faces[j], PETSC_TRUE,
                                       &num_closure, &val_closure);
        }
    }
  }

  /* Walk through the mesh and extract element types. */
  for (PetscInt i = 0; i < mNumberElementsLocal; i++) {

    /* Get and save the type of element i. */
    RealMat vtx = getElementCoordinateClosure(i); RealVec ctr(vtx.cols());
    for (PetscInt j = 0; j < mNumDim; j++) { ctr(j) = vtx.col(j).mean(); }
    std::string type = model->getElementType(ctr); mPointFields[i].insert(type);

    /* Add the type of element i to all mesh points connected via the Hasse graph. */
    PetscInt num_pts; const PetscInt *pts = NULL;
    DMPlexGetCone(mDistributedMesh, i, &pts); DMPlexGetConeSize(mDistributedMesh, i, &num_pts);
    /* for all d-1 mesh points attached to this element (i.e., edges/) */
    for (PetscInt j = 0; j < num_pts; j++) {
      /* insert this element type... */
      mPointFields[pts[j]].insert(type);
      /* for all mesh boundaries... */
      for (PetscInt k = 0; k < boundary_size; k++) {
        /* if this particular mesh point is on boundary set k... */
        if (std::find(mBndPts.begin(), mBndPts.end(), std::tie(k,pts[j]))
            != mBndPts.end()) {
          /* if boundary set k is labeled as homogeneous dirichlet... */
          auto hd = options->HomogeneousDirichlet();
          if (std::find(hd.begin(), hd.end(), model->SideSetName(k)) != hd.end()) {
            mPointFields[pts[j]].insert("boundary_homo_dirichlet");
            /* If we're on a boundary, it's important to the entire graph. */
            PetscInt num_closure; PetscInt *pts_closure = NULL;
            DMPlexGetTransitiveClosure(mDistributedMesh, pts[j], PETSC_TRUE, &num_closure,
                                       &pts_closure);
            /* Remember closure includes orientations (k += 2) */
            for (PetscInt l = 0; l < 2 * num_closure; l += 2) {
              mPointFields[pts_closure[l]].insert("boundary_homo_dirichlet");
            }
            DMPlexRestoreTransitiveClosure(mDistributedMesh, pts[j], PETSC_TRUE, &num_closure,
                                           &pts_closure);
          }
          /* default to free surface... insert nothing. */
        }
      }
    }

    /* Finally, add the type to the global fields. */
    mMeshFields.insert(type);

  }

}

void Mesh::setupGlobalDof(unique_ptr<Element> const &element,
                          unique_ptr<Options> const &options) {

  /* Mesh topology must be set up first. */
  if (!mMeshFields.size()) {
    throw std::runtime_error("Before setting up the degrees of freedom, you must attach the topology by "
                                 "calling setupTopology.");
  }

  /* Use the information extracted above to inform the global DOF layout. */
  PetscInt num_bc = 0;
  PetscInt poly_order = options->PolynomialOrder();
  /* Note here we are treating each field component as a separate vector (for now).
   * mMeshFields.size(); */
  PetscInt num_fields = 1;
  PetscInt *num_comps; PetscMalloc1(num_fields, &num_comps);
  {
    /* Note here that we are allocating each field component as a separate vector.
     * This may change in the future.
     * for (auto &f: mMeshFields) { num_comps[i++] = numFieldPerPhysics(f); } */
    for(int i=0; i<num_fields; i++) { num_comps[i] = 1; }
  }
  PetscInt *num_dof; PetscMalloc1(num_fields * (mNumDim + 1), &num_dof);

  /* This now just allocated the necessary memory per element. */
  for (PetscInt f = 0; f < num_fields; f++) {
    /* For each element point... (vertex, edge, face, volume) */
    for (PetscInt d = 0; d < (mNumDim + 1); d++) {

      PetscInt index = f * (mNumDim + 1) + d;
      switch (d) {

        /* 0-dimensional points (vertices). */
        case(0):

          num_dof[index] = element->NumDofVtx() * num_comps[f];
          break;

        /* 1-dimensional points (edges). */
        case(1):

          num_dof[index] = element->NumDofEdg() * num_comps[f];
          break;

        /* 2-dimensional points (faces). */
        case(2):

          num_dof[index] = element->NumDofFac() * num_comps[f];
          break;

        /* 3-dimensional points (cells). */
        case(3):

          num_dof[index] = element->NumDofVol() * num_comps[f];
          break;

        default:
          throw std::runtime_error("Unrecognized dimension! Dim: " + std::to_string(d));

      }

    }
  }

  
  /* Allocate the section. */
  DMPlexCreateSection(mDistributedMesh, mNumDim, num_fields, num_comps, num_dof,
                      num_bc, NULL, NULL, NULL, NULL, &mMeshSection);
  
  /* Attach some meta-information to the section. */
//  {
//    PetscInt i = 0;
//    for (auto &f: mMeshFields) { PetscSectionSetFieldName(mMeshSection, i++, f.c_str()); }
//  }

  /* Attach the section to our DM */
  DMSetDefaultSection(mDistributedMesh, mMeshSection);

  /* Only create a spectral basis if it makes sense. */
  if ((this->baseElementType() == "quad") ||
      (this->baseElementType() == "hex")  ||
      (this->baseElementType() == "tri")) {
    // Fixes edges and surfaces with rotated orientation given by neighboring element
    FixOrientation(mDistributedMesh, mMeshSection, num_comps, num_fields, poly_order, mNumDim);
    if ((this->baseElementType() == "quad") ||
        (this->baseElementType() == "hex")) {        
      // Sets up tensorized ordering for quads and hexes    
      DMPlexCreateSpectralClosurePermutation(mDistributedMesh, NULL);
    }
  }

  PetscFree(num_dof); 
  PetscFree(num_comps);
}

std::vector<PetscInt> Mesh::EdgeNumbers(const PetscInt elm) {

  std::vector<PetscInt> edges;
  /* Step up one level in the graph (reduce dim). */
  PetscInt csize; DMPlexGetConeSize(mDistributedMesh, elm, &csize);
  const PetscInt *cone; DMPlexGetCone(mDistributedMesh, elm, &cone);
  for (PetscInt i = 0; i < csize; i++) { edges.push_back(cone[i]); }
  return edges;

}

PetscInt Mesh::GetNeighbouringElement(const PetscInt interface, const PetscInt this_elm) const {

  PetscInt neighbour = -1;
  PetscInt ssize; DMPlexGetSupportSize(mDistributedMesh, interface, &ssize);
  const PetscInt *supp; DMPlexGetSupport(mDistributedMesh, interface, &supp);
  for (PetscInt i = 0; i < ssize; i++) {
    if (supp[i] != this_elm) { neighbour = supp[i]; }
  }
  if (neighbour == -1) { throw std::runtime_error(
        "Element " + std::to_string(this_elm) + " has no neighbour on face "
                   + std::to_string(interface)); }
  return neighbour;

}

std::vector<std::tuple<PetscInt,std::vector<std::string>>> Mesh::CouplingFields(
    const PetscInt elm) {

  std::vector<std::tuple<PetscInt,std::vector<std::string>>> couple;
  /* Step one level up in graph (reduce dim). */
  PetscInt csize; DMPlexGetConeSize(mDistributedMesh, elm, &csize);
  const PetscInt *cone; DMPlexGetCone(mDistributedMesh, elm, &cone);
  /*---------------------------------------------------------------------------
                    Get physics of neighbouring element.
   *--------------------------------------------------------------------------*/
  for (PetscInt i = 0; i < csize; i++) {
    /* Step through the support of each point (increase dim to element). */
    PetscInt ssize; DMPlexGetSupportSize(mDistributedMesh, cone[i], &ssize);
    const PetscInt *supp; DMPlexGetSupport(mDistributedMesh, cone[i], &supp);
    for (PetscInt j = 0; j < ssize; j++) {
      /* Get fields on this element, but skip parent element. */
      std::vector<std::string> fields;
      if (supp[j] == elm) { continue; }
      for (auto &f: mPointFields[supp[j]]) { fields.push_back(f); }
      couple.push_back(std::make_tuple(cone[i], fields));
    }
  }
  return couple;
}

std::vector<std::string> Mesh::TotalCouplingFields(const PetscInt elm) {

  std::set<std::string> coupling_fields;
  PetscInt num_pts; PetscInt *pts = NULL;
  DMPlexGetTransitiveClosure(mDistributedMesh, elm, PETSC_TRUE, &num_pts, &pts);
  for (PetscInt i = 0; i < 2*num_pts; i += 2) {
    for (auto &f: mPointFields[pts[i]]) {
      if (mPointFields[elm].find(f) == mPointFields[elm].end()) {
        coupling_fields.insert(f);
      }
    }
  }
  DMPlexRestoreTransitiveClosure(mDistributedMesh, elm, PETSC_TRUE, &num_pts, &pts);
  return std::vector<std::string> (coupling_fields.begin(), coupling_fields.end());

}

Eigen::MatrixXd Mesh::getElementCoordinateClosure(PetscInt elem_num) {

  Vec coord;
  DMGetCoordinatesLocal(mDistributedMesh, &coord);
  PetscSection coord_section;
  DMGetCoordinateSection(mDistributedMesh, &coord_section);

  PetscInt coord_buf_size;
  PetscReal *coord_buf = NULL;
  DMPlexVecGetClosure(mDistributedMesh, coord_section, coord, elem_num, &coord_buf_size, &coord_buf);

  int num_vtx = coord_buf_size / mNumDim;
  Eigen::MatrixXd vtx(num_vtx, mNumDim);
  for (int i = 0; i < num_vtx; i++) {
    vtx(i, 0) = coord_buf[mNumDim * i + 0];
    vtx(i, 1) = coord_buf[mNumDim * i + 1];
    if (mNumDim == 3) {
      vtx(i, 2) = coord_buf[mNumDim * i + 2];
    }
  }

  DMPlexVecRestoreClosure(mDistributedMesh, coord_section, coord, elem_num, &coord_buf_size, &coord_buf);
  return vtx;

}

PetscInt Mesh::numFieldPerPhysics(std::string physics) {
  PetscInt num;
    if (physics == "fluid") { num = 1; }
    else if (physics == "2delastic") { num = 2; }
    else if (physics == "3delastic") { num = 3; }
    else {
      throw std::runtime_error("Derived type " + physics + " is not known.");
    }
  return num;
}

std::string Mesh::baseElementType() {
  std::string type;
  RealMat vtx = getElementCoordinateClosure(0);
  if (vtx.rows() == 3) {
    type = "tri";
  } else if (vtx.rows() == 4 && mNumDim == 2) {
    type = "quad";
  } else if (vtx.rows() == 4 && mNumDim == 3) {
    type = "tet";
  } else if (vtx.rows() == 8 && mNumDim == 3) {
    type = "hex";
  } else {
    ERROR() << "Element type not detected: vtx.rows()=" << vtx.rows() << " and dim = " << mNumDim;
  }
  return type;
}


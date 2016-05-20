#include <mpi.h>

#include <petscdm.h>
#include <petscdmplex.h>
#include "Hexahedra.h"

// enables std::cout << vector;
template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  out << "[";
  size_t last = v.size() - 1;
  for(size_t i = 0; i < v.size(); ++i) {
    out << v[i];
    if (i != last) 
      out << ", ";
  }
  out << "]";
  return out;
}


Hexahedra::Hexahedra(Options options) {

  // Basic properties.
  mNumDim = 3;
  mNumVtx = 8;
  mPlyOrd = options.PolynomialOrder();

  mElmCtr.resize(3, 1);
  mVtxCrd.resize(3, mNumVtx);

  // Gll points.
  mNumDofVtx = 1;
  mNumDofEdg = mPlyOrd - 1;
  mNumDofFac = (mPlyOrd - 1) * (mPlyOrd - 1);
  mNumDofVol = (mPlyOrd - 1) * (mPlyOrd - 1) * (mPlyOrd - 1);

  // Integration points.
  mIntCoordR = Hexahedra::GllPointsForOrder(options.PolynomialOrder());
  mIntCoordS = Hexahedra::GllPointsForOrder(options.PolynomialOrder());
  mIntCoordT = Hexahedra::GllPointsForOrder(options.PolynomialOrder());
  mIntWeightsR = Hexahedra::GllIntegrationWeightForOrder(options.PolynomialOrder());
  mIntWeightsS = Hexahedra::GllIntegrationWeightForOrder(options.PolynomialOrder());
  mIntWeightsT = Hexahedra::GllIntegrationWeightForOrder(options.PolynomialOrder());

  // Save number of integration points.
  mNumIntPtsR = mIntCoordR.size();
  mNumIntPtsS = mIntCoordS.size();
  mNumIntPtsT = mIntCoordT.size();
  mNumIntPnt = mNumIntPtsR * mNumIntPtsS * mNumIntPtsT;

  // setup evaluated derivatives of test functions
  setupGradientOperator();
}

std::vector<int> surface_orient(int surface_orientation, std::vector<int> input_map) {
  if(surface_orientation == 0) {
    return input_map;
  }
  else if(surface_orientation == -4) {
    std::reverse(input_map.begin(),input_map.end());
    return input_map;
  }
  else {
    std::cerr <<
      "ERROR: We do not currently know how to deal with this surface orientation: "
              << surface_orientation << "\n";
    // exit(1);
  }
}

std::vector<int> edge_orient(std::vector<int> edge_orientation, std::vector<std::vector<int>> input_map) {

  std::vector<int> final_map(input_map.size()*2);
  for(int iedge=0; iedge < input_map.size(); iedge++) {
    if(edge_orientation[iedge] == 0 || -2) {
      final_map[iedge*2 + 0] = input_map[iedge][0];
      final_map[iedge*2 + 1] = input_map[iedge][1];
    }
    else if(edge_orientation[iedge] == -1) {
      final_map[iedge*2 + 0] = input_map[iedge][1];
      final_map[iedge*2 + 1] = input_map[iedge][0];
    }
  }
  return final_map;
}

void Hexahedra::BuildClosureMapping(DM &distributed_mesh) {
  mClsMap = ClosureMapping(mPlyOrd, mNumDim, distributed_mesh);
}

// void BuildRST(int order, int numdim, Mesh* mesh) {
//   if( order == 3) {
//     int* points = NULL;
//     int numPoints;
//     DMPlexGetTransitiveClosure(mesh->DistributedMesh(),mElmNum,PETSC_TRUE,&numPoints,&points);
//     std::vector<int> edge_orientations(12);
//     std::vector<int> surface_orientations(6);
//     int surface_counter = 0;
//     int edge_counter = 0;
//     for(int i=0;i<numPoints;i++) {
//       // surfaces
//       if(i>=1 && i<7) {
//         // note: points are grouped by 2s (point_id, surface_orientation) x numPoints
//         surface_orientations[surface_counter] = points[2*i+1];
//         surface_counter++;        
//       }
//       else if(i>=7 && i<19) {
//         // note: points are grouped by 2s (point_id, surface_orientation) x numPoints
//         edge_orientations[edge_counter] = points[2*i+1];
//         edge_counter++;
//         if(points[2*i+1] != 0 && points[2*i+1] != -1 && points[2*i+1] != -2) {
//           std::cerr <<
//             "ERROR: We do not currently know how to deal with this edge orientation: " << points[i] << "\n";
//           exit(1);
//         }
//       }
//     }
    
//   }
// }

// given precomputed point set
std::vector<int> getVertsFromPoint(int point, int numVerts, int numPoints, int* points);

// computes closure every time
std::vector<int> getVertsFromPoint(int point, int numVerts, DM &distributed_mesh);

void edgeHandler(std::vector<std::vector<int>> canonical_edges,
                 std::vector<std::vector<int>> edge_verts,
                 std::vector<int> o,
                 std::vector<std::vector<int>> std_layout,
                 int &dof_counter,
                 VectorXi &element_dof_map
                 );

void surfHandler(std::vector<int> canonical_verts, // the vertices as we expect them
                 std::vector<int> surf_verts, // vertices given from actual surface
                 int o, // mesh given orientation of this surface
                 std::vector<int> std_layout, // expected degrees of freedom layout for this surface
                 int &dof_counter, // current dof count
                 VectorXi &element_dof_map // final dof layout (output)
                 );

VectorXi internalMappingOrder3(int element, DM &distributed_mesh);

VectorXi Hexahedra::ClosureMapping(int order, int numdim, DM &distributed_mesh) {
  if( order == 3) {
    auto petsc_mapping = internalMappingOrder3(mElmNum,distributed_mesh);
    
    /**
       |  DOFs are ordered: {cell, faces, edges, vertices}
       |  The faces are ordered:
       |   * Bottom
       |   * Top   
       |   * Front 
       |   * Back  
       |   * Right 
       |   * Left  
    */
    
    VectorXi mapping_linear(64);
    for(int i=0;i<64;i++) { mapping_linear[i] = i; }
    
    std::vector<int> interior_mapping =
      {21, 22, 25, 26,
       37, 38, 41, 42};
    std::vector<int> bottom_surface_mapping = {5,9,10,6};
        
    std::vector<int> top_surface_mapping = {53,54,58,57};
    
    std::vector<int> front_surface_mapping = {17,18,34,33};
    
    std::vector<int> back_surface_mapping = {30,29,45,46};

    std::vector<int> right_surface_mapping = {23,27,43,39};

    std::vector<int> left_surface_mapping = {20,36,40,24};
    
    // edges    
    std::vector<int> bottom_edge_mapping = {4,8,13,14,11,7,2,1};
    
    std::vector<int> top_edge_mapping = {49,50,55,59,62,61,56,52};
    
    std::vector<int> front_edge_mapping = {19,35,32,16};
    
    std::vector<int> back_edge_mapping = {28,44,47,31};

    // right and left surface-edges are already covered
    std::vector<int> vertex_mapping = {0,12,15,3,48,51,63,60};

    std::vector<int> full_mapping;
    // concatenate them all together!
    full_mapping.insert(full_mapping.end(),
                        interior_mapping.begin(),
                        interior_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        bottom_surface_mapping.begin(),
                        bottom_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        top_surface_mapping.begin(),
                        top_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        front_surface_mapping.begin(),
                        front_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        back_surface_mapping.begin(),
                        back_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        right_surface_mapping.begin(),
                        right_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        left_surface_mapping.begin(),
                        left_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        bottom_edge_mapping.begin(),
                        bottom_edge_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        top_edge_mapping.begin(),
                        top_edge_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        front_edge_mapping.begin(),
                        front_edge_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        back_edge_mapping.begin(),
                        back_edge_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        vertex_mapping.begin(),
                        vertex_mapping.end());
    Map<VectorXi> mapping(full_mapping.data(),full_mapping.size());    
    
    // MatrixXi viewMap(64,2);
    // for(int i=0;i<64;i++) { viewMap(i,0) = i; }
    // viewMap.col(1) = petsc_mapping;    
    // // std::cout << "Hello!\n";
    // std::cout << "viewMap:\n" << viewMap << "\n";
    // return mapping_linear;
    VectorXi new_mapping(64);
    for(int i=0;i<64;i++) {
      new_mapping[i] = full_mapping[petsc_mapping[i]];
    }
    return new_mapping;
    // return petsc_mapping;
    // return mapping;
    // return mapping_linear;
  }
  
  else {
    std::cerr << "ERROR(ClosureMapping): Order for hexahedra not supported\n";
    VectorXi wrong_mapping(1);
    return wrong_mapping;
    
  }
}

VectorXd Hexahedra::GllPointsForOrder(const int order) {
  VectorXd gll_points(order + 1);
  if (order == 1) {
    gll_coordinates_order1_square(gll_points.data());
  } else if (order == 2) {
    gll_coordinates_order2_square(gll_points.data());
  } else if (order == 3) {
    gll_coordinates_order3_square(gll_points.data());
  } else if (order == 4) {
    gll_coordinates_order4_square(gll_points.data());
  } else if (order == 5) {
    gll_coordinates_order5_square(gll_points.data());
  } else if (order == 6) {
    gll_coordinates_order6_square(gll_points.data());
  } else if (order == 7) {
    gll_coordinates_order7_square(gll_points.data());
  } else if (order == 8) {
    gll_coordinates_order8_square(gll_points.data());
  } else if (order == 9) {
    gll_coordinates_order9_square(gll_points.data());
  } else if (order == 10) {
    gll_coordinates_order10_square(gll_points.data());
  }
  return gll_points;
}

VectorXd Hexahedra::GllIntegrationWeightForOrder(const int order) {
  VectorXd integration_weights(order + 1);
  if (order == 1) {
    gll_weights_order1_square(integration_weights.data());
  } else if (order == 2) {
    gll_weights_order2_square(integration_weights.data());
  } else if (order == 3) {
    gll_weights_order3_square(integration_weights.data());
  } else if (order == 4) {
    gll_weights_order4_square(integration_weights.data());
  } else if (order == 5) {
    gll_weights_order5_square(integration_weights.data());
  } else if (order == 6) {
    gll_weights_order6_square(integration_weights.data());
  } else if (order == 7) {
    gll_weights_order7_square(integration_weights.data());
  } else if (order == 8) {
    gll_weights_order8_square(integration_weights.data());
  } else if (order == 9) {
    gll_weights_order9_square(integration_weights.data());
  } else if (order == 10) {
    gll_weights_order10_square(integration_weights.data());
  }
  return integration_weights;
}

void Hexahedra::setupGradientOperator() {

  double s = mIntCoordS[0];
  mGrd.resize(mNumIntPtsR,mNumIntPtsS);
  Eigen::MatrixXd test(mNumIntPtsS, mNumIntPtsR);
  for (int i=0; i<mNumIntPtsR; i++) {
    double r = mIntCoordR[i];
    if (mPlyOrd == 1) {
      interpolate_eps_derivative_order1_square(s, test.data());
    } else if (mPlyOrd == 2) {
      interpolate_eps_derivative_order2_square(r, s, test.data());
    } else if (mPlyOrd == 3) {
      interpolate_eps_derivative_order3_square(r, s, test.data());
    } else if (mPlyOrd == 4) {
      interpolate_eps_derivative_order4_square(r, s, test.data());
    } else if (mPlyOrd == 5) {
      interpolate_eps_derivative_order5_square(r, s, test.data());
    } else if (mPlyOrd == 6) {
      interpolate_eps_derivative_order6_square(r, s, test.data());
    } else if (mPlyOrd == 7) {
      interpolate_eps_derivative_order7_square(r, s, test.data());
    } else if (mPlyOrd == 8) {
      interpolate_eps_derivative_order8_square(r, s, test.data());
    } else if (mPlyOrd == 9) {
      interpolate_eps_derivative_order9_square(r, s, test.data());
    } else if (mPlyOrd == 10) {
      interpolate_eps_derivative_order10_square(r, s, test.data());
    }
    mGrd.row(i) = test.col(0);
  }
  
  
}

Eigen::Map<const VectorXd> Hexahedra::rVectorStride(const VectorXd &function,
                                        const int &s_index, const int &t_index) {
  return Map<const VectorXd>(function.data() + s_index * mNumIntPtsS + t_index*mNumIntPtsS*mNumIntPtsR,
                             mNumIntPtsR);
}

Eigen::Map<const VectorXd, 0, InnerStride<>> Hexahedra::sVectorStride(const VectorXd &function,
                                                          const int &r_index,const int &t_index) {
  return Map<const VectorXd, 0, InnerStride<>>(
                                               function.data() + r_index + t_index*mNumIntPtsR*mNumIntPtsS, mNumIntPtsS,
                                               InnerStride<>(mNumIntPtsS));
}

Eigen::Map<const VectorXd, 0, InnerStride<>> Hexahedra::tVectorStride(const VectorXd &function,
                                                          const int &r_index,const int &s_index) {
  return Map<const VectorXd, 0, InnerStride<>>(
                                               function.data() + r_index + mNumIntPtsR*s_index, mNumIntPtsR,
                                               InnerStride<>(mNumIntPtsR*mNumIntPtsS));
}

VectorXd Hexahedra::__attachMaterialProperties(ExodusModel *model,std::string parameter_name) {

  Eigen::VectorXd material_at_vertices(mNumVtx);

  for (auto i = 0; i < mNumVtx; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(
                                                                           mElmCtr, parameter_name, i);
  }
  return material_at_vertices;
  
}

void Hexahedra::attachVertexCoordinates(DM &distributed_mesh) {

  // needs building somewhere
  BuildClosureMapping(distributed_mesh);
  
  Vec coordinates_local;
  PetscInt coordinate_buffer_size;
  PetscSection coordinate_section;
  PetscReal *coordinates_buffer = NULL;

  DMGetCoordinatesLocal(distributed_mesh, &coordinates_local);
  DMGetCoordinateSection(distributed_mesh, &coordinate_section);
  DMPlexVecGetClosure(distributed_mesh, coordinate_section, coordinates_local, mElmNum,
                      &coordinate_buffer_size, &coordinates_buffer);
  std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer + coordinate_buffer_size);
  DMPlexVecRestoreClosure(distributed_mesh, coordinate_section, coordinates_local, mElmNum,
                          &coordinate_buffer_size, &coordinates_buffer);

  for (int i = 0; i < mNumVtx; i++) {
    mVtxCrd(0, i) = coordinates_element[mNumDim * i + 0];
    mVtxCrd(1, i) = coordinates_element[mNumDim * i + 1];
    mVtxCrd(2, i) = coordinates_element[mNumDim * i + 2];
  }

  // Save element center
  mElmCtr <<
    mVtxCrd.row(0).mean(),
    mVtxCrd.row(1).mean(),
    mVtxCrd.row(2).mean();
  
}

/* Derivation for reference to "real" element mapping
/*  reference hex, with r,s,t=[-1,1]x[-1,1]x[-1,1]
 *  
 *                                   (v7)                    (v6)
 *                                     /--------------d------/+
 *                                  /--|                  /-- |
 *                              /---   |              /---    |
 *                           /--       |       f   /--        |
 *                        /--          |         --           |
 *                 (v4) +---------------b--------+ (v5)       |
 *                      |              |         |            |
 *                      |              |         |            |
 *                      |              |       g |            |
 *    ^                 |              |         |            |
 *    | (t)             |              |         |            |
 *    |                 |            --+-------- |----c-------+ (v2)
 *    |                 |        ---/ (v1)       |         /--
 *    |                 |     --/              e |     /---
 *    |         (s)     | ---/                   |  /--
 *    |         /-      +/--------------a--------+--
 *    |     /---       (v0)                    (v3)
 *    | /---
 *    +-----------------> (r)
 * 
 *  
 * // note (r-1)/2.0 to bring range from [-1,1] to [0,1]
 * a = v0 + (v3-v0)*(r+1.0)/2.0
 * b = v4 + (v5-v4)*(r+1.0)/2.0
 * c = v1 + (v2-v1)*(r+1.0)/2.0
 * d = v7 + (v6-v7)*(r+1.0)/2.0
 * 
 * e = a + (c-a)*(s+1.0)/2.0
 * f = b + (d-b)*(s+1.0)/2.0
 * 
 * g = e + (f-e)*(t+1.0)/2.0
 * using sympy to expand `g` as a function of only vertices v0..7 and r,s,t gives the formula below
 */
PetscReal referenceToElementMapping(PetscReal v0, PetscReal v1, PetscReal v2, PetscReal v3, PetscReal v4, PetscReal v5, PetscReal v6, PetscReal v7, PetscReal r, PetscReal s, PetscReal t);

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> Hexahedra::buildNodalPoints() {
  assert(mNumIntPnt == mNumIntPtsR * mNumIntPtsS * mNumIntPtsT);
  
  VectorXd nodalPoints_x(mNumIntPnt);
  VectorXd nodalPoints_y(mNumIntPnt);
  VectorXd nodalPoints_z(mNumIntPnt);

  // assumes right-hand rule vertex layout (same as PETSc)
  double v0x = mVtxCrd(0,0);
  double v1x = mVtxCrd(0,1);
  double v2x = mVtxCrd(0,2);
  double v3x = mVtxCrd(0,3);
  double v4x = mVtxCrd(0,4);
  double v5x = mVtxCrd(0,5);
  double v6x = mVtxCrd(0,6);
  double v7x = mVtxCrd(0,7);

  double v0y = mVtxCrd(1,0);
  double v1y = mVtxCrd(1,1);
  double v2y = mVtxCrd(1,2);
  double v3y = mVtxCrd(1,3);
  double v4y = mVtxCrd(1,4);
  double v5y = mVtxCrd(1,5);
  double v6y = mVtxCrd(1,6);
  double v7y = mVtxCrd(1,7);

  double v0z = mVtxCrd(2,0);
  double v1z = mVtxCrd(2,1);
  double v2z = mVtxCrd(2,2);
  double v3z = mVtxCrd(2,3);
  double v4z = mVtxCrd(2,4);
  double v5z = mVtxCrd(2,5);
  double v6z = mVtxCrd(2,6);
  double v7z = mVtxCrd(2,7);

  
  
  int idx = 0;
  for (auto k = 0; k < mNumIntPtsT; k++) {
    for (auto j = 0; j < mNumIntPtsS; j++) {
      for (auto i = 0; i < mNumIntPtsR; i++) {
        double r = mIntCoordR(i);
        double s = mIntCoordS(j);
        double t = mIntCoordT(k);

        // this represents `g` above (for each x,y,z) (copy paste from hexahedra_element_metrics.py)
        nodalPoints_x(idx) = referenceToElementMapping(v0x,v1x,v2x,v3x,v4x,v5x,v6x,v7x,r,s,t);
        nodalPoints_y(idx) = referenceToElementMapping(v0y,v1y,v2y,v3y,v4y,v5y,v6y,v7y,r,s,t);
        nodalPoints_z(idx) = referenceToElementMapping(v0z,v1z,v2z,v3z,v4z,v5z,v6z,v7z,r,s,t);

        // if(mElmNum == 0) {
        //   printf("(%f,%f,%f)->(%f,%f,%f)\n",r,s,t,nodalPoints_x[idx],nodalPoints_y[idx],nodalPoints_z[idx]);
        // }
        
        idx++;
      }
    }
  }
  return std::make_tuple(nodalPoints_x,nodalPoints_y,nodalPoints_z);
}

std::tuple<Matrix3d, PetscReal> Hexahedra::inverseJacobianAtPoint(PetscReal r,
                                                                  PetscReal s,
                                                                  PetscReal t) {

  // assumes right-hand rule vertex layout (same as PETSc)
  double v0x = mVtxCrd(0,0);
  double v1x = mVtxCrd(0,1);
  double v2x = mVtxCrd(0,2);
  double v3x = mVtxCrd(0,3);
  double v4x = mVtxCrd(0,4);
  double v5x = mVtxCrd(0,5);
  double v6x = mVtxCrd(0,6);
  double v7x = mVtxCrd(0,7);

  double v0y = mVtxCrd(1,0);
  double v1y = mVtxCrd(1,1);
  double v2y = mVtxCrd(1,2);
  double v3y = mVtxCrd(1,3);
  double v4y = mVtxCrd(1,4);
  double v5y = mVtxCrd(1,5);
  double v6y = mVtxCrd(1,6);
  double v7y = mVtxCrd(1,7);

  double v0z = mVtxCrd(2,0);
  double v1z = mVtxCrd(2,1);
  double v2z = mVtxCrd(2,2);
  double v3z = mVtxCrd(2,3);
  double v4z = mVtxCrd(2,4);
  double v5z = mVtxCrd(2,5);
  double v6z = mVtxCrd(2,6);
  double v7z = mVtxCrd(2,7);

  PetscReal dxdr=-0.5*v0x + 0.5*v3x + 0.5*(s + 1.0)*(0.5*v0x - 0.5*v1x + 0.5*v2x - 0.5*v3x) + 0.5*(t + 1.0)*(0.5*v0x - 0.5*v3x - 0.5*v4x + 0.5*v5x - 0.5*(s + 1.0)*(0.5*v0x - 0.5*v1x + 0.5*v2x - 0.5*v3x) + 0.5*(s + 1.0)*(0.5*v4x - 0.5*v5x + 0.5*v6x - 0.5*v7x));
PetscReal dydr=-0.5*v0y + 0.5*v3y + 0.5*(s + 1.0)*(0.5*v0y - 0.5*v1y + 0.5*v2y - 0.5*v3y) + 0.5*(t + 1.0)*(0.5*v0y - 0.5*v3y - 0.5*v4y + 0.5*v5y - 0.5*(s + 1.0)*(0.5*v0y - 0.5*v1y + 0.5*v2y - 0.5*v3y) + 0.5*(s + 1.0)*(0.5*v4y - 0.5*v5y + 0.5*v6y - 0.5*v7y));
PetscReal dzdr=-0.5*v0z + 0.5*v3z + 0.5*(s + 1.0)*(0.5*v0z - 0.5*v1z + 0.5*v2z - 0.5*v3z) + 0.5*(t + 1.0)*(0.5*v0z - 0.5*v3z - 0.5*v4z + 0.5*v5z - 0.5*(s + 1.0)*(0.5*v0z - 0.5*v1z + 0.5*v2z - 0.5*v3z) + 0.5*(s + 1.0)*(0.5*v4z - 0.5*v5z + 0.5*v6z - 0.5*v7z));
PetscReal dxds=-0.5*v0x + 0.5*v1x - 0.25*(r + 1.0)*(-v0x + v3x) + 0.25*(r + 1.0)*(-v1x + v2x) + 0.5*(t + 1.0)*(0.5*v0x - 0.5*v1x - 0.5*v4x + 0.5*v7x + 0.25*(r + 1.0)*(-v0x + v3x) - 0.25*(r + 1.0)*(-v1x + v2x) - 0.25*(r + 1.0)*(-v4x + v5x) + 0.25*(r + 1.0)*(v6x - v7x));
PetscReal dyds=-0.5*v0y + 0.5*v1y - 0.25*(r + 1.0)*(-v0y + v3y) + 0.25*(r + 1.0)*(-v1y + v2y) + 0.5*(t + 1.0)*(0.5*v0y - 0.5*v1y - 0.5*v4y + 0.5*v7y + 0.25*(r + 1.0)*(-v0y + v3y) - 0.25*(r + 1.0)*(-v1y + v2y) - 0.25*(r + 1.0)*(-v4y + v5y) + 0.25*(r + 1.0)*(v6y - v7y));
PetscReal dzds=-0.5*v0z + 0.5*v1z - 0.25*(r + 1.0)*(-v0z + v3z) + 0.25*(r + 1.0)*(-v1z + v2z) + 0.5*(t + 1.0)*(0.5*v0z - 0.5*v1z - 0.5*v4z + 0.5*v7z + 0.25*(r + 1.0)*(-v0z + v3z) - 0.25*(r + 1.0)*(-v1z + v2z) - 0.25*(r + 1.0)*(-v4z + v5z) + 0.25*(r + 1.0)*(v6z - v7z));
PetscReal dxdt=-0.5*v0x + 0.5*v4x - 0.25*(r + 1.0)*(-v0x + v3x) + 0.25*(r + 1.0)*(-v4x + v5x) - 0.25*(s + 1.0)*(-v0x + v1x - 0.5*(r + 1.0)*(-v0x + v3x) + 0.5*(r + 1.0)*(-v1x + v2x)) + 0.25*(s + 1.0)*(-v4x + v7x - 0.5*(r + 1.0)*(-v4x + v5x) + 0.5*(r + 1.0)*(v6x - v7x));
PetscReal dydt=-0.5*v0y + 0.5*v4y - 0.25*(r + 1.0)*(-v0y + v3y) + 0.25*(r + 1.0)*(-v4y + v5y) - 0.25*(s + 1.0)*(-v0y + v1y - 0.5*(r + 1.0)*(-v0y + v3y) + 0.5*(r + 1.0)*(-v1y + v2y)) + 0.25*(s + 1.0)*(-v4y + v7y - 0.5*(r + 1.0)*(-v4y + v5y) + 0.5*(r + 1.0)*(v6y - v7y));
PetscReal dzdt=-0.5*v0z + 0.5*v4z - 0.25*(r + 1.0)*(-v0z + v3z) + 0.25*(r + 1.0)*(-v4z + v5z) - 0.25*(s + 1.0)*(-v0z + v1z - 0.5*(r + 1.0)*(-v0z + v3z) + 0.5*(r + 1.0)*(-v1z + v2z)) + 0.25*(s + 1.0)*(-v4z + v7z - 0.5*(r + 1.0)*(-v4z + v5z) + 0.5*(r + 1.0)*(v6z - v7z));
    
  Matrix3d J;
  J <<
    dxdr, dydr, dzdr,
    dxds, dyds, dzds,
    dxdt, dydt, dzdt;

  PetscReal detJ = J.determinant();

  // detJ * invJ_ = [ dyds*dzdt - dydt*dzds, -dydr*dzdt + dydt*dzdr,  dydr*dzds - dyds*dzdr]
  //                [-dxds*dzdt + dxdt*dzds,  dxdr*dzdt - dxdt*dzdr, -dxdr*dzds + dxds*dzdr]
  //                [ dxds*dydt - dxdt*dyds, -dxdr*dydt + dxdt*dydr,  dxdr*dyds - dxds*dydr]  
  
  // need to test if this quasi-analytical formula is faster than simply J.inverse();
  
  double rx = (1 / detJ) * (dyds*dzdt - dydt*dzds);
  double sx = (1 / detJ) * (-dydr*dzdt + dydt*dzdr);
  double tx = (1 / detJ) * (dydr*dzds - dyds*dzdr);

  double ry = (1 / detJ) * (-dxds*dzdt + dxdt*dzds);
  double sy = (1 / detJ) * ( dxdr*dzdt - dxdt*dzdr);
  double ty = (1 / detJ) * (-dxdr*dzds + dxds*dzdr);
  
  double rz = (1 / detJ) * (dxds*dydt - dxdt*dyds);
  double sz = (1 / detJ) * (-dxdr*dydt + dxdt*dydr);
  double tz = (1 / detJ) * (dxdr*dyds - dxds*dydr);

  Matrix3d inverseJacobian;
  inverseJacobian <<
    rx, sx, tx,
    ry, sy, ty,
    rz, sz, tz;

  return std::make_tuple(inverseJacobian, detJ);
}

Vector3d Hexahedra::inverseCoordinateTransform(const double &x_real, const double &y_real, const double &z_real) {

  // assumes right-hand rule vertex layout (same as PETSc)
  double v0x = mVtxCrd(0,0);
  double v1x = mVtxCrd(0,1);
  double v2x = mVtxCrd(0,2);
  double v3x = mVtxCrd(0,3);
  double v4x = mVtxCrd(0,4);
  double v5x = mVtxCrd(0,5);
  double v6x = mVtxCrd(0,6);
  double v7x = mVtxCrd(0,7);

  double v0y = mVtxCrd(1,0);
  double v1y = mVtxCrd(1,1);
  double v2y = mVtxCrd(1,2);
  double v3y = mVtxCrd(1,3);
  double v4y = mVtxCrd(1,4);
  double v5y = mVtxCrd(1,5);
  double v6y = mVtxCrd(1,6);
  double v7y = mVtxCrd(1,7);

  double v0z = mVtxCrd(2,0);
  double v1z = mVtxCrd(2,1);
  double v2z = mVtxCrd(2,2);
  double v3z = mVtxCrd(2,3);
  double v4z = mVtxCrd(2,4);
  double v5z = mVtxCrd(2,5);
  double v6z = mVtxCrd(2,6);
  double v7z = mVtxCrd(2,7);

  // Using Newton iterations
  // https://en.wikipedia.org/wiki/Newton%27s_method#Nonlinear_systems_of_equations
  // J_F(xn)(x_{n+1} - x_n) = -F(x_n)
  // Solve for x_{n+1}
  // where J_F(x_n) is jacobian:
  // https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
  double tol = 1e-6;
  int num_iter = 0;
  Vector3d solution{0.0, 0.0, 0.0}; // initial guess at (0.0,0.0)
  while (true) {

    double r = solution(0);
    double s = solution(1);
    double t = solution(2);

    Matrix3d jacobian_inverse_t;
    PetscReal detJ;

    // mapping from reference hex [-1,1]x[-1,1]x[-1,1] to *this* element
    PetscReal Tx = referenceToElementMapping(v0x,v1x,v2x,v3x,v4x,v5x,v6x,v7x,r,s,t);
    PetscReal Ty = referenceToElementMapping(v0y,v1y,v2y,v3y,v4y,v5y,v6y,v7y,r,s,t);
    PetscReal Tz = referenceToElementMapping(v0z,v1z,v2z,v3z,v4z,v5z,v6z,v7z,r,s,t);
     
    Vector3d objective_function{x_real - Tx, y_real - Ty, z_real - Tz};
    
    std::tie(jacobian_inverse_t,detJ) = inverseJacobianAtPoint(r,s,t);
    
    if ((objective_function.array().abs() < tol).all()) {
      return solution;
    } else {
      solution += (jacobian_inverse_t.transpose() * objective_function);
    }
    if (num_iter > 10) {
      std::cerr << "inverseCoordinateTransform: TOO MANY ITERATIONS!\n";
      exit(1);
    }
    num_iter++;

  }

}

Eigen::VectorXd Hexahedra::interpolateAtPoint(double r, double s, double t) {
  /**  
   *  reference hex, with r,s,t=[-1,1]x[-1,1]x[-1,1]
   * 
   *                                      (v7)                    (v6)
   *                                        /--------------d------/+
   *                                     /--|                  /-- |
   *                                 /---   |              /---    |
   *                              /--       |       f   /--        |
   *                           /--          |         --           |
   *                    (v4) +---------------b--------+ (v5)       |        \gamma
   *                         |              |         |            |          ^
   *                         |              |         |            |          |
   *                         |              |       g |            |          |
   *    ^                    |              |         |            |          |
   *    | (t)                |              |         |            |          |
   *    |                    |            --+-------- |----c-------+ (v2)     |
   *    |                    |        ---/ (v1)       |         /--          ---
   *    |                    |     --/              e |     /---               /> \beta
   *    |         (s)        | ---/                   |  /--                 /-
   *    |         /-         +/--------------a--------+--                  --
   *    |     /---          (v0)                    (v3)
   *    | /---               |--------------> \alpha
   *    +-----------------> (r)
   * 
   * alpha = (r+1)/2.0
   * beta = (s+1)/2.0
   * gamma = (t+1)/2.0
   * 
   * a = alpha*v3 + (1-alpha)*v0
   * b = alpha*v5 + (1-alpha)*v4
   * c = alpha*v2 + (1-alpha)*v1
   * d = alpha*v6 + (1-alpha)*v7
   * 
   * e = beta*c + (1-beta)*a
   * f = beta*d + (1-beta)*b
   * 
   * g = gamma*f + (1-gamma)*e
   
   * group vf by terms v0..v7                                       
   * = [v0,...,v7].dot([-0.125*r*s*t + 0.125*r*s + ... - 0.125*s - 0.125*t + 0.125,
   *                    ...,
   *                    -0.125*r*s*t - 0.125*r*s - ... + 0.125*s + 0.125*t + 0.125]
   */
  
  Eigen::VectorXd interpolator(8);
  interpolator <<
    -0.125*r*s*t + 0.125*r*s + 0.125*r*t - 0.125*r + 0.125*s*t - 0.125*s - 0.125*t + 0.125,
    0.125*r*s*t - 0.125*r*s + 0.125*r*t - 0.125*r - 0.125*s*t + 0.125*s - 0.125*t + 0.125,
    -0.125*r*s*t + 0.125*r*s - 0.125*r*t + 0.125*r - 0.125*s*t + 0.125*s - 0.125*t + 0.125,
    0.125*r*s*t - 0.125*r*s - 0.125*r*t + 0.125*r + 0.125*s*t - 0.125*s - 0.125*t + 0.125,
    0.125*r*s*t + 0.125*r*s - 0.125*r*t - 0.125*r - 0.125*s*t - 0.125*s + 0.125*t + 0.125,
    -0.125*r*s*t - 0.125*r*s + 0.125*r*t + 0.125*r - 0.125*s*t - 0.125*s + 0.125*t + 0.125,
    0.125*r*s*t + 0.125*r*s + 0.125*r*t + 0.125*r + 0.125*s*t + 0.125*s + 0.125*t + 0.125,
    -0.125*r*s*t - 0.125*r*s - 0.125*r*t - 0.125*r + 0.125*s*t + 0.125*s + 0.125*t + 0.125;
  return interpolator;
}

void Hexahedra::attachSource(std::vector<std::shared_ptr<Source>> sources) {
  printf("TODO: attachSource");
  exit(1);
}

void Hexahedra::prepareStiffness() {}

void Hexahedra::assembleElementMassMatrix(Mesh *mesh) {}

MatrixXd Hexahedra::computeSourceTerm(double time) {

  printf("TODO: computeSourceTerm");
  exit(1);
  MatrixXd ret;
  return ret;
}

MatrixXd Hexahedra::computeStiffnessTerm(const MatrixXd &displacement) {

  printf("TODO: computeStiffnessTerm");
  exit(1);
  MatrixXd ret;
  return ret;
  
}

void Hexahedra::attachMaterialProperties(ExodusModel *model) {
  printf("TODO: attachMaterialProperties");
  exit(1);
}


void Hexahedra::applyDirichletBoundaries(Mesh *mesh, Options &options, const std::string &fieldname) {

  if (! mBndElm) return;

  double value = 0;
  auto dirchlet_boundary_names = options.DirichletBoundaries();
  for (auto &bndry: dirchlet_boundary_names) {
    auto faceids = mBnd[bndry];
    for (auto &faceid: faceids) {
      auto field = mesh->getFieldOnFace(fieldname, faceid);
      field = 0 * field.array() + value;
      mesh->setFieldFromFace(fieldname, faceid, field);
    }
  }
}


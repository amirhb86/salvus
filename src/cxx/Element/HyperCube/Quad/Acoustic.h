#pragma once

#include <petscvec.h>
#include "Utilities/Options.h"
#include "Element/HyperCube/Quad.h"
#include "Model/ExodusModel.h"
#include "Utilities/Utilities.h"
#include "Mesh/Mesh.h"

using namespace Eigen;

class AcousticQuad: public Quad {

  Vector4d mMaterialVelocityAtVertices;
  MatrixXd mElementStrain;

 public:

  AcousticQuad(Options options);

  AcousticQuad *clone() const { return new AcousticQuad(*this); }

  /**
   * Empty as its not needed for quadrilaterals (build the stiffness matrix on the fly).
   */
  void prepareStiffness() { };

  void assembleElementMassMatrix(Mesh *mesh);
  void interpolateMaterialProperties(ExodusModel *model);
  void setInitialCondition(Mesh *mesh, VectorXd &pts_x, VectorXd &pts_z,
                           double L, double x0, double z0);

  /***************************************************************************
   *                              TIME LOOP.
   ***************************************************************************/
  void computeSurfaceTerm();
  MatrixXd computeSourceTerm(double time);
  MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);
  std::vector<std::string> PullElementalFields() const { return {"u"}; }
  std::vector<std::string> PushElementalFields() const { return {"a"}; }

  /***************************************************************************
   *                            TESTING HELPERS
   ***************************************************************************/

  /**
   * Computes the exact solution (for the case of the first
   * Eigenfunctions with Dirchelet boundaries) at a given physical point.
   * @param [in] pts_x Physical GLL locations (x).
   * @param [in] ptx_z Physical GLL locations (z).
   * @param [in] L Length of one side of mesh.
   * @param [in] x0 Center of displacement (x).
   * @param [in] z0 Center of displacement (z).
   * @param [in] time Simulation time.
   * @returns A vector containing the exact solution at GLL points.
   */
  VectorXd exactSolution(Eigen::VectorXd &pts_x, Eigen::VectorXd &pts_z,
                                double L, double x0, double z0, double time);

  /**
   * Setup initial conditions for tests.
   * @param [in] mesh The mesh.
   * @param [in] options The options class.
   */
  void setupTest(Mesh *mesh, Options options);

  /**
   * Check exact solution against current displacement
   * @param [in] mesh The mesh.
   * @param [in] displacement SEM (approximate) displacement.
   * @param [in] time Simulation time.
   */
  double checkTest(Mesh *mesh, Options options,
                   const MatrixXd &displacement, double time);


};


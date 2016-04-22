#pragma once

#include <petscvec.h>
#include "Utilities/Options.h"
#include "Element/HyperCube/Hexahedra.h"
#include "Model/ExodusModel.h"
#include "Utilities/Utilities.h"
#include "Mesh/Mesh.h"

using namespace Eigen;

class AcousticHex: public Hexahedra {

  VectorXd mMaterialVelocityAtVertices;
  MatrixXd mElementStrain;

 public:

  AcousticHex(Options options);

  /**
   * Empty as its not needed for quadrilaterals (build the stiffness matrix on the fly).
   */
  void prepareStiffness() { };

  void assembleElementMassMatrix(Mesh *mesh);
  void attachMaterialProperties(ExodusModel *model);
  void setInitialCondition(Mesh *mesh, VectorXd &pts_x, VectorXd &pts_y, VectorXd &pts_z,
                           double L, double x0, double y0, double z0);

  /***************************************************************************
   *                              TIME LOOP.
   ***************************************************************************/
  void computeSurfaceTerm();
  MatrixXd computeSourceTerm(double time);
  MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);
  std::vector<std::string> PullElementalFields() const { return {"u"}; }
  std::vector<std::string> PushElementalFields() const { return {"a"}; }

  void attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers);
  /***************************************************************************
   *                            TESTING HELPERS
   ***************************************************************************/

  /**
   * Computes the exact solution (for the case of the first
   * Eigenfunctions with Dirchelet boundaries) at a given physical point.
   * @param [in] pts_x Physical GLL locations (x).
   * @param [in] ptx_z Physical GLL locations (y).
   * @param [in] ptx_z Physical GLL locations (z).
   * @param [in] L Length of one side of mesh.
   * @param [in] x0 Center of displacement (x).
   * @param [in] y0 Center of displacement (y).
   * @param [in] z0 Center of displacement (z).
   * @param [in] time Simulation time.
   * @returns A vector containing the exact solution at GLL points.
   */
  VectorXd exactSolution(Eigen::VectorXd &pts_x, Eigen::VectorXd &pts_y, Eigen::VectorXd &pts_z,
                         double L, double x0, double y0, double z0, double time);

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

  MatrixXd interpolateFieldAtPoint(const VectorXd &pnt);
  void recordField(const MatrixXd &u);
  
  std::shared_ptr<Element> clone() const { return std::shared_ptr<Element> (new AcousticHex(*this)); }
  

};


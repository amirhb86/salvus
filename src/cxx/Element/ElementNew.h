#pragma once

#include <Eigen/Dense>
#include <Model/ExodusModel.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>
#include <mpi.h>
#include <petsc.h>
#include <Receiver/Receiver.h>
#include <memory>

class ElementNew {

 public:

  virtual std::shared_ptr<ElementNew> clone() const = 0;
//  virtual ~ElementNew() {};
//  static std::shared_ptr<ElementNew> Factory(Options options) {};
  virtual void attachMaterialPropertiesNew(ExodusModel *model) = 0;
  virtual void attachVertexCoordinates(DM &distributed_mesh) = 0;
  virtual void attachSource(std::vector<std::shared_ptr<Source>> sources) = 0;
  virtual void attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers) = 0;
  virtual void prepareStiffness() = 0;
  virtual void assembleElementMassMatrix(Mesh *mesh) = 0;
  virtual void recordField(const Eigen::Ref<const Eigen::MatrixXd>& field) = 0;
  virtual Eigen::MatrixXd computeSourceTerm(const double time) = 0;
  virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::Ref<const Eigen::MatrixXd>& u) = 0;
  virtual std::vector<std::string> PullElementalFields() const = 0;
  virtual std::vector<std::string> PushElementalFields() const = 0;
  virtual Eigen::MatrixXd interpolateFieldAtPoint(const Eigen::Ref<const Eigen::VectorXd>& pnt) = 0;
  virtual void setBoundaryConditions(Mesh *mesh) = 0;
  virtual void setupTest(Mesh *mesh, Options options) = 0;
  virtual double checkTest(Mesh *mesh, Options options, const Eigen::Ref<const Eigen::MatrixXd>& u, double time) = 0;


  virtual inline bool BndElm() const = 0;
  virtual inline void SetNum(const int num) = 0;
  virtual inline int Num() const = 0;
  virtual inline int NumDim() const = 0;
  virtual inline int NumDofVol() const = 0;
  virtual inline int NumDofFac() const = 0;
  virtual inline int NumDofEdg() const = 0;
  virtual inline int NumDofVtx() const = 0;
  virtual inline int NumIntPnt() const = 0;
  virtual inline Eigen::VectorXi ClsMap() const = 0;

  virtual void applyDirichletBoundaries(Mesh *mesh, Options options, const std::string &fieldname) = 0;

};
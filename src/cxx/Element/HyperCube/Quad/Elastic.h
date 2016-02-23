//
// Created by Michael Afanasiev on 2016-02-17.
//

#ifndef SALVUS_ELASTIC_H
#define SALVUS_ELASTIC_H


#include <Eigen/Dense>
#include "../Quad.h"

class Elastic : public Quad {

    Eigen::Vector4d mC11AtVertices;
    Eigen::Vector4d mC13AtVertices;
    Eigen::Vector4d mC15AtVertices;
    Eigen::Vector4d mC31AtVertices;
    Eigen::Vector4d mC33AtVertices;
    Eigen::Vector4d mC35AtVertices;
    Eigen::Vector4d mC51AtVertices;
    Eigen::Vector4d mC53AtVertices;
    Eigen::Vector4d mC55AtVertices;
    Eigen::Vector4d mRhoAtVertices;

    Eigen::Matrix<double, 2, Eigen::Dynamic> mElementStrain;

public:

    Elastic(Options options);
    virtual Elastic *clone() const { return new Elastic(*this); }

    virtual void checkInField(Mesh *mesh);
    virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);
    virtual Eigen::MatrixXd computeSourceTerm();
    virtual void computeSurfaceTerm();
    virtual void assembleElementMassMatrix(Mesh *mesh);
    virtual void interpolateMaterialProperties(ExodusModel *model);
    virtual void setInitialCondition(Mesh *mesh, Eigen::VectorXd& pts_x, Eigen::VectorXd& ptz_z);


};


#endif //SALVUS_ELASTIC_H

//
// Created by Michael Afanasiev on 2016-02-17.
//

#ifndef SALVUS_ELASTIC_H
#define SALVUS_ELASTIC_H


#include <Eigen/Dense>
#include "../Quad.h"

class Elastic : public Quad {

    Eigen::Vector4d mC11AtVertices;
    Eigen::Vector4d mC12AtVertices;
    Eigen::Vector4d mC13AtVertices;
    Eigen::Vector4d mC21AtVertices;
    Eigen::Vector4d mC22AtVertices;
    Eigen::Vector4d mC23AtVertices;
    Eigen::Vector4d mC31AtVertices;
    Eigen::Vector4d mC32AtVertices;
    Eigen::Vector4d mC33AtVertices;

    Eigen::Matrix<double, 2, Eigen::Dynamic> mElementStrain;

public:

    Elastic(Options options);
    virtual Elastic *clone() const { return new Elastic(*this); }

    virtual void checkInField(Mesh *mesh);
    virtual Eigen::MatrixXd checkOutField(Mesh *mesh, const std::string name);
    virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);
    virtual void computeSourceTerm();
    virtual void computeSurfaceTerm();
    virtual void assembleMassMatrix();
    virtual void interpolateMaterialProperties(ExodusModel *model);

};


#endif //SALVUS_ELASTIC_H

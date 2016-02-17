//
// Created by Michael Afanasiev on 2016-02-17.
//

#include "Elastic.h"

Elastic::Elastic(Options options): Quad(options) {

    mMassMatrix.setZero(mNumberIntegrationPoints);
    mElementStrain.setZero(2, mNumberIntegrationPoints*2);

}

void Elastic::checkInField(Mesh *mesh) {

}

Eigen::MatrixXd Elastic::computeSourceTerm() {

}

void Elastic::computeSurfaceTerm() {

}

void Elastic::assembleElementMassMatrix(Mesh *mesh) {

}

Eigen::MatrixXd Elastic::computeStiffnessTerm(const Eigen::MatrixXd &displacement) {
    Eigen::MatrixXd test;
    return test;
}

void Elastic::interpolateMaterialProperties(ExodusModel *model) {

}

void Elastic::setInitialCondition(Mesh *mesh, Eigen::VectorXd &pts_x, Eigen::VectorXd &ptz_z) {

}

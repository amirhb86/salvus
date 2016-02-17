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

Eigen::MatrixXd Elastic::checkOutField(Mesh *mesh, const std::string name) {
    return Quad::checkOutField(mesh, name);
}

void Elastic::computeSourceTerm() {

}

void Elastic::computeSurfaceTerm() {

}

void Elastic::assembleMassMatrix() {

}

Eigen::MatrixXd Elastic::computeStiffnessTerm(const Eigen::MatrixXd &displacement) {
    Eigen::MatrixXd test;
    return test;
}

void Elastic::interpolateMaterialProperties(ExodusModel *model) {

}

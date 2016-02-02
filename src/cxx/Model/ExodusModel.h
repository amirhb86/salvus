//
// Created by Michael Afanasiev on 2016-02-01.
//

#ifndef SALVUS_EXODUSMODEL_H
#define SALVUS_EXODUSMODEL_H

#include <vector>
#include "../Options.h"
extern "C" {
#include "../kdtree.h"
#include "exodusII.h"
};

class ExodusModel {

    int mExodusId;
    int mNumberDimension;
    int mNumberVertices;
    int mNumberElements;
    int mNumberElementBlocks;
    int mNumberNodeSets;
    int mNumberSideSets;
    int mNumberNodalVariables;

    char mExodusTitle[MAX_LINE_LENGTH+1];

    float mExodusVersion;

    kdtree *mKdTree;

    std::string mExodusFileName;

    std::vector<int> mKdTreeData;
    std::vector<double> mNodalX;
    std::vector<double> mNodalY;
    std::vector<double> mNodalZ;
    std::vector<double> mNodalVariables;

    std::vector<std::string> mNodalVariableNames;

    void getInitialization();
    void readCoordinates();
    void readNodalVariables();
    void createKdTree();
    void exodusError(const int retval, std::string func_name);

public:

    ExodusModel(Options options);
    void initializeParallel();
    PetscReal getMaterialParameterAtPoint(const std::vector<double> point,
                                          const std::string parameter_name);

};


#endif //SALVUS_EXODUSMODEL_H

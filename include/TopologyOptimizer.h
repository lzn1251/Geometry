#ifndef GEOMETRY_TOPOLOGY_OPTIMIZER_H
#define GEOMETRY_TOPOLOGY_OPTIMIZER_H 

#include "Basic.h"

class TopologyOptimizer {
public:
    struct Parameter {
        int nelx = 60;            // number of elements in x direction
        int nely = 20;            // number of elements in y direction
        double volfrac = 0.4;     // volume constraint
        double rmin = 1.5;        // filter radius
        int maxIter = 100;        // maximum number of iterations
        double penal = 3.0;       // SIMP penaly factor
    };

    TopologyOptimizer(const Parameter& param);

    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& optimize();

    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& getDensity() const;

private:
    void FEA(Eigen::VectorXd& U);

    void filterSensitivity(const Eigen::VectorXd& dc, Eigen::VectorXd& dc_filtered);

    void ocUpdate(const Eigen::VectorXd& dc, double& change);

    Parameter m_Params;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_Density;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_OldDensity;
    Eigen::VectorXd m_Dc;                       // sensitivity
    Eigen::SparseMatrix<double> m_K;            // stiffness matrix
    std::vector<Eigen::Triplet<double>> m_Triplets;
};

#endif
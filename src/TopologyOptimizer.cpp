
#include "TopologyOptimizer.h"

TopologyOptimizer::TopologyOptimizer(const Parameter& params)
    : m_Params(params)
{
    m_Density = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Ones(m_Params.nely, m_Params.nelx);
    m_OldDensity = m_Density;
    m_Dc.resize(m_Params.nely * m_Params.nelx);
}

const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& TopologyOptimizer::getDensity() const
{
    return m_Density;
}

void TopologyOptimizer::FEA(Eigen::VectorXd& U) 
{
    int nelx = m_Params.nelx;
    int nely = m_Params.nely;
    int ndof = 2 * (nelx + 1) * (nely + 1);

    // init
    m_Triplets.clear();
    m_K = Eigen::SparseMatrix<double>(ndof, ndof);

    // element stiffness matrix
    Eigen::Matrix<double, 8, 8> KE;
    double E = 1.0;
    double nu = 0.3;
    double k = E / (1.0 - nu * nu);

    KE <<
       k * (1.0 / 3.0), k * (1.0 / 6.0), 0.0, k * (-1.0 / 6.0), k * (-1.0 / 3.0), k * (-1.0 / 6.0), 0.0, k * (1.0 / 6.0),
       k * (1.0 / 6.0), k * (1.0 / 3.0), k * (1.0 / 6.0), 0.0, k * (-1.0 / 6.0), k * (-1.0 / 3.0), k * (-1.0 / 6.0), 0.0,
       0.0, k * (1.0 / 6.0), k * (1.0 / 3.0), k * (1.0 / 6.0), 0.0, k * (-1.0 / 6.0), k * (-1.0 / 3.0), k * (-1.0 / 6.0),
       k * (-1.0 / 6.0), 0.0, k * (1.0 / 6.0), k * (1.0 / 3.0), k * (-1.0 / 6.0), 0.0, k * (-1.0 / 6.0), k * (-1.0 / 3.0),
       k * (-1.0 / 3.0), k * (-1.0 / 6.0), 0.0, k * (1.0 / 6.0), k * (1.0 / 3.0), k * (1.0 / 6.0), 0.0, k * (-1.0 / 6.0), 
       k * (-1.0 / 6.0), k * (-1.0 / 3.0), k * (-1.0 / 6.0), 0.0, k * (1.0 / 6.0), k * (1.0 / 3.0), k * (1.0 / 6.0), 0.0,
       0.0, k * (-1.0 / 6.0), k * (-1.0 / 3.0), k * (-1.0 / 6.0), 0.0, k * (1.0 / 6.0), k * (1.0 / 3.0), k * (1.0 / 6.0),
       k * (1.0 / 6.0), 0.0, k * (-1.0 / 6.0), k * (-1.0 / 3.0), k * (-1.0 / 6.0), 0.0, k * (1.0 / 6.0), k * (1.0 / 3.0);

    for (int ely = 0; ely < nely; ++ely) {
        for (int elx = 0; elx < nelx; ++elx) { 
            int n1 = (nely + 1 - ely) * (nelx + 1) + elx;
            int n2 = (nely + 1 - ely) * (nelx + 1) + elx + 1;
            int n3 = (nely - ely) * (nelx + 1) + elx + 1;
            int n4 = (nely - ely) * (nelx + 1) + elx;

            std::vector<int> edof = {
                2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1, 2 * n3, 2 * n3 + 1, 2 * n4, 2 * n4 + 1
            };

            // SIMP
            double rho = m_Density(ely, elx);
            double penal = m_Params.penal;
            double Ee = 1e-3 + rho * rho * rho * (1.0 - 1e-3);

            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < 8; ++j) {
                    m_Triplets.emplace_back(edof[i], edof[j], Ee * KE(i, j));
                }
            }
        }
    }

    m_K.setFromSortedTriplets(m_Triplets.begin(), m_Triplets.end());

    Eigen::VectorXd F = Eigen::VectorXd::Zero(ndof);
    F(2 * (nelx + 1) * (nely + 1) - 1) = -1.0;       // last dof

    std::vector<bool> fixed(ndof, false);
    for (int i = 0; i <= nely; ++i) {
        fixed[2 * i * (nelx + 1)] = true;
        fixed[2 * i * (nelx + 1) + 1] = true;
    }

    // KU = F
    Eigen::SparseMatrix<double> K_reduced;
    Eigen::VectorXd F_reduced;
    std::vector<int> free_dofs;
    for (int i = 0; i < ndof; ++i) {
        if (!fixed[i]) free_dofs.push_back(i);
    }

    std::vector<Eigen::Triplet<double>> reduced_triplets;
    for (const auto& t : m_Triplets) {
        if (!fixed[t.row()] && !fixed[t.col()]) {
            int newRow = std::lower_bound(free_dofs.begin(), free_dofs.end(), t.row()) - free_dofs.begin();
            int newCol = std::lower_bound(free_dofs.begin(), free_dofs.end(), t.col()) - free_dofs.begin();
            reduced_triplets.emplace_back(newRow, newCol, t.value());
        }
    }

    K_reduced.resize(free_dofs.size(), free_dofs.size());
    K_reduced.setFromTriplets(reduced_triplets.begin(), reduced_triplets.end());

    F_reduced.resize(free_dofs.size());
    for (size_t i = 0; i < free_dofs.size(); ++i) {
        F_reduced[i] = F(free_dofs[i]);
    }

    // solve
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(K_reduced);
    Eigen::VectorXd U_reduced = solver.solve(F_reduced);

    U = Eigen::VectorXd::Zero(ndof);
    for (size_t i = 0; i < free_dofs.size(); ++i) {
        U(free_dofs[i]) = U_reduced(i);
    }
}

void TopologyOptimizer::filterSensitivity(const Eigen::VectorXd& dc, Eigen::VectorXd& dc_filtered) {
    int nelx = m_Params.nelx;
    int nely = m_Params.nely;
    dc_filtered.setZero();

    for (int i = 0; i < nelx; ++i) {
        for (int j = 0; j < nely; ++j) {
            double sum = 0.0;
            for (int k = std::max(i - (int)m_Params.rmin, 0); k <= std::min(i + (int)m_Params.rmin, nely - 1); ++k) {
                for (int l = std::max(j - (int)m_Params.rmin, 0); l <= std::min(j + (int)m_Params.rmin, nelx - 1); ++l) {
                    double dist = std::sqrt((i - k) * (i - k) + (j - l) * (j - l));
                    if (dist <= m_Params.rmin) {
                        double weight = m_Params.rmin - dist;
                        sum += weight;
                        dc_filtered(i * nelx + j) += weight * dc(k * nelx + l);
                    }
                }
            }
            dc_filtered(i * nelx + j) /= sum;
        }
    }
}

void TopologyOptimizer::ocUpdate(const Eigen::VectorXd& dc, double& change) {
    int nelx = m_Params.nelx;
    int nely = m_Params.nely;
    double l1 = 0.0, l2 = 1e9, lmid;
    double move = 0.2;

    // 拉格朗日乘子
    while (l2 - l1 > 1e-6) {
        lmid = (l1 + l2) / 2;
        double vol = 0.0;
        for (int i = 0; i < nely; ++i) {
            for (int j = 0; j < nelx; ++j) {
                double xnew = std::max(0.0, std::max(m_Density(i, j) - move,
                          std::min(1.0, std::min(m_Density(i, j) + move, 
                          m_Density(i, j) * std::sqrt(-dc(i * nelx + j) / lmid)))));
                vol += xnew;
            }
        }
        if (vol > m_Params.volfrac * nelx * nely) {
            l1 = lmid;
        } else {
            l2 = lmid;
        }
    }

    change = 0.0;
    for (int i = 0; i < nely; i++) {
        for (int j = 0; j < nelx; j++) {
            double xnew = std::max(0.0, std::max(m_Density(i, j) - move, 
                std::min(1.0, std::min(m_Density(i, j) + move, 
                m_Density(i, j) * std::sqrt(-dc(i * nelx + j) / l2)))));
            change = std::max(change, std::abs(xnew - m_Density(i, j)));
            m_OldDensity(i, j) = m_Density(i, j);
            m_Density(i, j) = xnew;
        }
    }
    
}


const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& TopologyOptimizer::optimize() {
    Eigen::VectorXd U;
    double change = 1.0;

    for (int iter = 0; iter < m_Params.maxIter && change > 0.01; ++iter) {
        FEA(U);

        m_Dc.setZero();
        int nelx = m_Params.nelx;
        int nely = m_Params.nely;
        for (int ely = 0; ely < nely; ++ely) {
            for (int elx = 0; elx < nelx; ++elx) {
                int n1 = (nely + 1 - ely) * (nelx + 1) + elx;
                int n2 = (nely + 1 - ely) * (nelx + 1) + elx + 1;
                int n3 = (nely - ely) * (nelx + 1) + elx + 1;
                int n4 = (nely - ely) * (nelx + 1) + elx;
                std::vector<int> edof = {
                    2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1, 2 * n3, 2 * n3 + 1, 2 * n4, 2 * n4 + 1
                };

                Eigen::Matrix<double, 8, 1> Ue;
                for (int i = 0; i < 8; ++i) {
                    Ue(i) = U(edof[i]);
                }

                double rho = m_Density(ely, elx);
                double penal = m_Params.penal;
                double dc_val = -penal * rho * (rho - 1e-3) * (rho - 1e-3) * Ue.transpose() *
                        (Eigen::Matrix<double, 8, 8>::Identity() * (1.0 - 1e-3)) * Ue;
                m_Dc(ely * nelx + elx) = dc_val;
            }
        }

        Eigen::VectorXd dc_filtered;
        filterSensitivity(m_Dc, dc_filtered);

        ocUpdate(dc_filtered, change);

        double vol = m_Density.sum() / (nelx * nely);
        std::cout << "Iter " << iter << ", Volume: " << vol << ", Change = " << change << std::endl;
    }

    return m_Density;
}
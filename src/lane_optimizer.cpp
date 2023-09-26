// Copyright 2023 watson.wang

#include "lane_optimizer.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <iostream>

LaneOptimizer::LaneOptimizer(P_weight weight, A_setups constrains)
    : _weight(weight), _constrains(constrains)
{
    qp_solver_.settings()->setMaxIteration(200000);
    qp_solver_.settings()->setAdaptiveRhoInterval(0);  // 0 means automatic
    qp_solver_.settings()->setRelativeTolerance(1.0e-4);  // def: 1.0e-4
    qp_solver_.settings()->setAbsoluteTolerance(1.0e-8);  // def: 1.0e-4
    //qp_solver_.settings()->setVerbosity(false);
    qp_solver_.settings()->setWarmStart(true);
}

LaneOptimizer::OptimizationResult LaneOptimizer::optimize(const OptimizationData & data)
{
    const auto sample = data.sample;
    const Eigen::Index N = sample.size();
    const double l_end = data.l_end;
    const double dl_end = data.dl_end;
    const double ddl_end = data.ddl_end;

    const Eigen::Index IDX_L = 0;
    const Eigen::Index IDX_DL = N;
    const Eigen::Index IDX_DDL = 2 * N;
    const Eigen::Index l_variables = 3 * N;
    const Eigen::Index l_constraints = _constrains.border * N + _constrains.heading * N +
       _constrains.kappa * N + _constrains.jerk * (N - 1) + _constrains.start * 3 + 
       _constrains.ddl2dl * (N - 1) + _constrains.dl2l * (N - 1);


    // the matrix size depends on constraint numbers.
    Eigen::SparseMatrix<double> A(l_constraints, l_variables);
    Eigen::VectorXd lower_bound(l_constraints);
    Eigen::VectorXd upper_bound(l_constraints);

    // Object Variables
    Eigen::SparseMatrix<double> P(l_variables, l_variables);
    Eigen::VectorXd q = Eigen::VectorXd::Zero(l_variables);

    // Object Function
    for (Eigen::Index i = 0; i < N; ++i) {
        double s2 = sample[i].ref_ds_ * sample[i].ref_ds_;
        P.coeffRef(IDX_L + i, IDX_L + i) += 2*_weight.l;
        P.coeffRef(IDX_L + i, IDX_L + i) += 2 * _weight.l;
        P.coeffRef(IDX_DL + i, IDX_DL + i) += 2 * _weight.dl;
        P.coeffRef(IDX_DDL + i, IDX_DDL + i) += 2 * (_weight.ddl + _weight.dddl / s2);
        if (i < (N - 1)) {
            P.coeffRef(IDX_DDL + i, IDX_DDL + i + 1) += 2 * (-2 * _weight.ddl / s2);
        }
        if ((i > 0) && (i < (N - 1))) {
            P.coeffRef(IDX_DDL + i, IDX_DDL + i + 1) += 2 * (_weight.dddl / s2);
        }
    }

    for (Eigen::Index i = 0; i < N; ++i) {
        P.coeffRef(IDX_L + i, IDX_L + i) += 2 * _weight.l_ref;
        q(IDX_L + i) += -2.0 * _weight.l_ref * sample.at(i).ref_l_;
    }


    P.coeffRef(N - 1, N - 1) += 2 * _weight.l_end;
    P.coeffRef(N * 2 - 1, N * 2 - 1) += 2 * _weight.dl_end;
    P.coeffRef(N * 3 - 1, N * 3 - 1) += 2 * _weight.ddl_end;

    q(N - 1) += -2.0 * l_end * _weight.l_end;
    q(N * 2 - 1) += -2.0 * dl_end * _weight.dl_end;
    q(N * 3 - 1) += -2.0 * ddl_end * _weight.ddl_end;


    // Constraint
    int constr_idx = 0;

    if (_constrains.border) {
       for (int i = 0; i < N; ++i, ++constr_idx) {
           A.insert(constr_idx, IDX_L + i) = 1.0;
           upper_bound(constr_idx) = sample.at(i).max_l_;
           lower_bound(constr_idx) = sample.at(i).min_l_;
       }
    }

    if (_constrains.heading) {
       for (int i = 0; i < N; ++i, ++constr_idx) {
          A.insert(constr_idx, IDX_DL + i) = 1.0;
          upper_bound(constr_idx) = 2.0;
          lower_bound(constr_idx) = -2.0;
       }
    }

    if (_constrains.kappa) {
       for (int i = 0; i < N; ++i, ++constr_idx) {
          A.insert(constr_idx, IDX_DDL + i) = 1.0;
          upper_bound(constr_idx) = data.acc_kappa + sample.at(i).kappa_;
          lower_bound(constr_idx) = -data.acc_kappa + sample.at(i).kappa_;
       }
    }


    if (_constrains.jerk) {
       for (int i = 0; i < N - 1; ++i, ++constr_idx) {
          A.insert(constr_idx, IDX_DDL + i) = -1.0;
          A.insert(constr_idx, IDX_DDL + i + 1) = 1.0;
          upper_bound(constr_idx) = sample.at(i).ref_ds_ * sample.at(i).jerk_;
          lower_bound(constr_idx) = -sample.at(i).ref_ds_ * sample.at(i).jerk_;
       }
    }

    if (_constrains.start) {
       A.insert(constr_idx, 0) = 1.0;
       upper_bound(constr_idx) = data.l_start;
       lower_bound(constr_idx) = data.l_start;
       constr_idx++;

       A.insert(constr_idx, N) = 1.0;
       upper_bound(constr_idx) = data.dl_start;
       lower_bound(constr_idx) = data.dl_start;
       constr_idx++;

       A.insert(constr_idx, 2 * N) = 1.0;
       upper_bound(constr_idx) = data.ddl_start;
       lower_bound(constr_idx) = data.ddl_start;
       constr_idx++;
    }

    if (_constrains.ddl2dl) {
       for (int i = 0; i < N - 1; ++i, ++constr_idx) {
          A.insert(constr_idx, IDX_DL + i) = -1.0;
          A.insert(constr_idx, IDX_DL + i + 1) = 1.0;
          A.insert(constr_idx, IDX_DDL + i) = -0.5 * sample[i].ref_ds_;
          A.insert(constr_idx, IDX_DDL + i + 1) = -0.5 * sample[i].ref_ds_;
          upper_bound(constr_idx) = 0.0;
          lower_bound(constr_idx) = 0.0;
       }
    }

    if (_constrains.dl2l) {
       for (int i = 0; i < N - 1; ++i, ++constr_idx) {
          A.insert(constr_idx, IDX_L + i) = -1.0;
          A.insert(constr_idx, IDX_L + i + 1) = 1.0;
          A.insert(constr_idx, IDX_DL + i) = -sample[i].ref_ds_;
          A.insert(constr_idx, IDX_DDL + i) = -0.5 * sample[i].ref_ds_ * sample[i].ref_ds_;
          A.insert(constr_idx, IDX_DDL + i + 1) = -(1.0 / 6.0) * sample[i].ref_ds_ * sample[i].ref_ds_ * sample[i].ref_ds_;
          upper_bound(constr_idx) = 0;
          lower_bound(constr_idx) = 0;
       }
    }

    qp_solver_.data()->setNumberOfVariables((int)l_variables);
    qp_solver_.data()->setNumberOfConstraints((int)l_constraints);

    auto qp_solver_run = [&]()-> bool {

        if (!qp_solver_.data()->setHessianMatrix(P)) return false;
        if (!qp_solver_.data()->setGradient(q)) return false;
        if (!qp_solver_.data()->setLinearConstraintsMatrix(A)) return false;
        if (!qp_solver_.data()->setLowerBound(lower_bound)) return false;
        if (!qp_solver_.data()->setUpperBound(upper_bound)) return false;

        // instantiate the solver
        if (!qp_solver_.initSolver()) return false;

        // execute optimization
        if (OsqpEigen::ErrorExitFlag::NoError != qp_solver_.solveProblem()) return false;

        return true;
    };

    auto sts = qp_solver_run();
    Eigen::VectorXd QPSolution;
    OptimizationResult optimized_result;

    if (sts) {
        QPSolution = qp_solver_.getSolution();
        //const auto has_nan =
        //  std::any_of(optval.begin(), optval.end(), [](const auto v) { return std::isnan(v); });
        //if (has_nan) std::cerr << "optimization failed : result contains NaN values\n";
        // 
        //std::cout << QPSolution << std::endl;

        std::vector<double> opt_time(N);
        std::vector<double> opt_pos(N);
        std::vector<double> opt_vel(N);
        std::vector<double> opt_acc(N);
        std::vector<double> opt_jerk(N);

        for (int i = 0; i < N; ++i) {
            opt_time.at(i) = sample.at(i).ref_ds_;
            opt_pos.at(i) = QPSolution(IDX_L + i);
            opt_vel.at(i) = (std::max)(QPSolution(IDX_DL + i), 0.0);
            opt_acc.at(i) = QPSolution(IDX_DDL + i);
        }
        opt_vel.back() = 0.0;

        optimized_result.t = opt_time;
        optimized_result.s = opt_pos;
        optimized_result.v = opt_vel;
        optimized_result.a = opt_acc;
        optimized_result.j = opt_jerk;
    }
    else {
        std::cerr << "optimization failed " << std::endl;
    }


    return optimized_result;
}

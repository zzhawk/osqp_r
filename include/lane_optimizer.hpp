// Copyright 2023 watson.wang

#ifndef LAINE_OPTIMIZER_HPP_
#define LAINE_OPTIMIZER_HPP_

#include "s_boundary.hpp"
#include "OsqpEigen/OsqpEigen.h"

#include <vector>

class LaneOptimizer
{
public:
   struct P_weight {
      double l;
      double dl;
      double ddl;
      double dddl;
      double l_end;
      double dl_end;
      double ddl_end;
      double l_ref;
   };

   struct A_setups {
      bool border;
      bool heading;
      bool kappa;
      bool jerk;
      bool start;
      bool ddl2dl;
      bool dl2l;
   };

    struct OptimizationData
    {
        double l_start;
        double dl_start;
        double ddl_start;
        double l_end;
        double dl_end;
        double ddl_end;
        double acc_kappa;
        SamplePoints sample;
    };

  struct OptimizationResult
  {
    std::vector<double> t;
    std::vector<double> s;
    std::vector<double> v;
    std::vector<double> a;
    std::vector<double> j;
  };

  explicit LaneOptimizer(P_weight weight, A_setups constrains);

  OptimizationResult optimize(const OptimizationData & data);

private:
    // Parameter
   P_weight _weight{};
   A_setups _constrains{};

    // QPSolver
    OsqpEigen::Solver qp_solver_;
};

#endif  // LAINE_OPTIMIZER_HPP_

#ifndef SAMCODE__LINEARSOLVERS__TRILINOS_SOLVER_HPP
#define SAMCODE__LINEARSOLVERS__TRILINOS_SOLVER_HPP

#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

namespace Teuchos
{
    class ParameterList;
}


enum TrilinosSolverType {BiCGStab,
                         CG,
                         GMRES,
                         BiCGStab_ML,
                         CG_ML,
                         GMRES_ML,
                         DeterminedByInit};

class TrilinosSolver
{
public:
    TrilinosSolver(const TrilinosSolverType solver_type);
    ~TrilinosSolver(){}

    void init();
    void solve(int           size,
               const int*    ia,
               const int*    ja,
               const double* sa,
               const double* b,
               double*       x) const;
private:
    boost::multi_array<int, 1>    options_;
    boost::multi_array<double, 1> parameters_;

    boost::shared_ptr<Teuchos::ParameterList> ml_parameters_;

    bool multilevel_preconditioner_;
    bool static_profile_;
    bool a_matrix_do_copy_;
    bool b_vector_do_copy_;
    bool x_result_do_copy_;
    bool unknown_solver_;
};

#endif

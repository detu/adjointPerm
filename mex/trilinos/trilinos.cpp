#include "trilinos.hpp"

#include <cassert>
#include <iostream>

#include <boost/scoped_ptr.hpp>

#include <trilinos/Epetra_SerialComm.h>
#include <trilinos/Epetra_Map.h>
#include <trilinos/Epetra_CrsMatrix.h>
#include <trilinos/Epetra_Vector.h>
#include <trilinos/AztecOO.h>
#include <trilinos/Teuchos_ParameterList.hpp>
#include <trilinos/ml_defs.h>
#include <trilinos/ml_config.h>
#include <trilinos/ml_MultiLevelPreconditioner.h>


namespace
{
    std::string id_CG       = "CG";
    std::string id_GMRES    = "GMRES";
    std::string id_BiCGStab = "BiCGStab";
    std::string id_TFQMR    = "TFQMR";
}

TrilinosSolver::TrilinosSolver(const TrilinosSolverType solver_type)
    :
    options_         (boost::extents[AZ_OPTIONS_SIZE]),
    parameters_      (boost::extents[AZ_PARAMS_SIZE]),
    ml_parameters_   (new Teuchos::ParameterList),
    multilevel_preconditioner_(false),
    static_profile_  (false), a_matrix_do_copy_(false),
    b_vector_do_copy_(false), x_result_do_copy_(false),
    unknown_solver_  (false)
{
    AZ_defaults(&options_[0], &parameters_[0]);
    switch(solver_type) {
        case BiCGStab:
            options_[AZ_solver]         = AZ_bicgstab;
            multilevel_preconditioner_  = false;
            break;

        case CG:
            options_[AZ_solver]         = AZ_cg;
            multilevel_preconditioner_  = false;
            break;

        case GMRES:
            options_[AZ_solver]         = AZ_gmres;
            multilevel_preconditioner_  = false;
            break;

        case BiCGStab_ML:
            options_[AZ_solver]         = AZ_bicgstab;
            multilevel_preconditioner_  = true;
            break;

        case CG_ML:
            options_[AZ_solver]         = AZ_cg;
            multilevel_preconditioner_  = true;
            break;

        case GMRES_ML:
            options_[AZ_solver]         = AZ_gmres;
            multilevel_preconditioner_  = true;
            break;

        case DeterminedByInit:
            break;

        default:
            ;//("No recognizable solver type given.");
    }
}

void TrilinosSolver::init() {
    assert(options_.size() == AZ_OPTIONS_SIZE);
    assert(parameters_.size() == AZ_PARAMS_SIZE);


    if (false)//(estimate_condition_number)
    {
        if (options_[AZ_solver] == AZ_bicgstab)
        {
            std::cout << "Warning: Can not estimate condition number for BiCGStab\n";
        }
        else if (options_[AZ_solver] == AZ_cg)
        {
            options_[AZ_solver] = AZ_cg_condnum;
        }
        else
        {
            assert(options_[AZ_solver] == AZ_gmres);
            options_[AZ_solver] = AZ_gmres_condnum;
        }
    }
    int solver_output = 2;
    if (solver_output == 0)
    {
        options_[AZ_output]      = AZ_none;
    }
    else if (solver_output == 1)
    {
        options_[AZ_output]      = AZ_last;
    }
    else if (solver_output == 2)
    {
        options_[AZ_output]      = AZ_warnings;
    }
    else
    {
        options_[AZ_output]      = 1;
    }

    options_[AZ_scaling]         = AZ_none;
    options_[AZ_precond]         = AZ_none;
    options_[AZ_conv]            = AZ_r0;
    options_[AZ_pre_calc]        = AZ_calc;
    options_[AZ_max_iter]        = 5000;
    options_[AZ_kspace]          = 30;
    options_[AZ_reorder]         = 1;
    options_[AZ_keep_info]       = 0;
    options_[AZ_orthog]          = AZ_classic;
    options_[AZ_aux_vec]         = AZ_resid;
    //options_[AZ_subdomain_solve] =
    //options_[AZ_graph_fill]      =
    //options_[AZ_poly_ord]        =
    //options_[AZ_overlap]         =
    //options_[AZ_type_overlap]    =

    parameters_[AZ_tol]          = 1e-8;
    //parameters_[AZ_drop]         =
    //parameters_[AZ_ilut_fill]    =
    //parameters_[AZ_omega]        =
    //parameters_[AZ_weights]      =

    if (multilevel_preconditioner_) {
        ML_Epetra::SetDefaults("SA", *ml_parameters_);

        ml_parameters_->set("output", 0);
        ml_parameters_->set("cycle applications", 1);
        ml_parameters_->set("max levels", 5);
        ml_parameters_->set("PDE equations", 1);
        ml_parameters_->set("increasing or decreasing", "increasing");
        ml_parameters_->set("aggregation: threshold", 0.06);
        ml_parameters_->set("aggregation: type", "MIS");
        ml_parameters_->set("aggregation: damping factor", 1.33);
        ml_parameters_->set("aggregation: symmetrize", true);
        ml_parameters_->set("smoother: type", "MLS");
        ml_parameters_->set("smoother: pre or post", "both");
        //ml_parameters_->set("smoother: damping factor", 0.67);
        ml_parameters_->set("smoother: sweeps", 4);
        //ml_parameters_->set("smoother: MLS alpha", 30.0);
        ml_parameters_->set("smoother: MLS polynomial order", 3);
        ml_parameters_->set("coarse: max size", 16);
        ml_parameters_->set("coarse: type", "Amesos-KLU");
        ml_parameters_->set("coarse: sweeps", 1);
        //ml_parameters_->set("coarse: damping factor", 0.67);
    }
}

void TrilinosSolver::solve(int           size,
                           const int*    ia,
                           const int*    ja,
                           const double* sa,
                           const double* b,
                           double*       x) const
{
    Epetra_SerialComm Comm;
    Epetra_Map        row_map(size, 0, Comm);
    Epetra_Map        col_map(size, 0, Comm);
    Epetra_DataAccess cv = a_matrix_do_copy_ ? Copy : View;

    boost::multi_array<int, 1> sizes(boost::extents[size]);

    for (int i = 0; i < size; ++i)
    {
        sizes[i] = ia[i+1] - ia[i];
    }

    Epetra_CrsMatrix matrix(cv, row_map, col_map, &sizes[0], static_profile_);

    for (int i = 0; i < size; ++i)
    {
        matrix.InsertGlobalValues(i, sizes[i], &const_cast<double&>(sa[ia[i]]), &const_cast<int&>(ja[ia[i]]));
    }
    matrix.FillComplete();

    Epetra_Vector vector(Epetra_DataAccess(b_vector_do_copy_ ? Copy : View),
                         Epetra_Map(size, 0, Epetra_SerialComm()), &const_cast<double&>(b[0]));
    Epetra_Vector result(Epetra_DataAccess(x_result_do_copy_ ? Copy : View),
                         Epetra_Map(size, 0, Epetra_SerialComm()), &x[0]);

    AztecOO solver(&matrix, &result, &vector);

    // Be careful: the order of setting options and preconditioner is not commutative!!!!
    boost::multi_array<int,    1> options_copy   (options_   );
    boost::multi_array<double, 1> parameters_copy(parameters_);

    solver.SetAllAztecOptions(&options_copy   [0]);
    solver.SetAllAztecParams (&parameters_copy[0]);

    boost::scoped_ptr<ML_Epetra::MultiLevelPreconditioner> ml_prec;
    Teuchos::ParameterList ml_parameters_copy(*ml_parameters_);
    if (multilevel_preconditioner_)
    {
        ml_prec.reset(new ML_Epetra::MultiLevelPreconditioner(matrix,
                                                              ml_parameters_copy,
                                                              true));
        solver.SetPrecOperator(ml_prec.get());
    }

    solver.Iterate(options_copy[AZ_max_iter], parameters_copy[AZ_tol]);

    if (x_result_do_copy_) {
        for (int i = 0; i < size; ++i) {
            x[i] = result[i];
        }
    }
}

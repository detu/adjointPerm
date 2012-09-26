#ifdef MATLAB_MEX_FILE
#include <matlab/mex.h>

#include <string>
#include <vector>

#include "trilinos.hpp"




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    if (nrhs != 2) {
	// FAILURE - wrong number of arguments
	nlhs = 0;
	mexErrMsgTxt("matlab_solve requires 3 arguments.\nCall with 'matlab_solve(A, b, parameterfile)'" );

    }

    int      m        = mxGetM (prhs[0]);
    int      n        = mxGetN (prhs[0]);
    mwIndex* array_ia = mxGetJc(prhs[0]);
    mwIndex* array_ja = mxGetIr(prhs[0]);
    const double*  array_da = mxGetPr(prhs[0]);
    const double*  b        = mxGetPr(prhs[1]);


    TrilinosSolver solver(CG_ML);
    solver.init();
    const std::vector<int> ia(&array_ia[0], &array_ia[m+1]);
    int N = array_ia[m];
    const std::vector<int> ja(&array_ja[0], &array_ja[N]);

    std::vector<double> x(m, 0.0);
    solver.solve(*(&m), &ia[0], &ja[0], &array_da[0], &b[0], &x[0]);

    nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    double* ret = mxGetPr(plhs[0]);
    for (int i = 0; i < m; ++i) {
	ret[i] = x[i];
    }
}
#endif

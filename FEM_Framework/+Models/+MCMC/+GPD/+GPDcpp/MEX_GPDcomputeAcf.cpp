/*
 * MEX_MET_FNTcomputeAcf.cpp
 *
 *  Created on: Aug 14, 2012
 *      Author: kaisero
 */
#include <iostream>
#include <sstream>
#include <exception>
#include "mex.h"
#include "gsl/gsl_matrix.h"
#include "GSL_Helper.h"
#include "MCMCMet.h"
#include "MCMCMetGPD.h"
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	GSL_Helper helper;
	stringstream ss;

	if (nrhs != 2) {
		ss.clear();
		ss << "2 inputs required, got " << nrhs << ", try again";
		mexErrMsgTxt(ss.str().c_str());
	}
	if (nlhs != 1) {
		ss.clear();
		ss << "one output required. got " << nlhs << ", try again";
		mexErrMsgTxt(ss.str().c_str());
	}
	double *o_gamma, *o_resid;

	// mex matrizen
	const mxArray *mxGamma, *mxResid;

	/*  create a pointer to the input matrix */
	mxGamma = prhs[0];
	mxResid = prhs[1];


	/* get dimensions */
	int K = mxGetM(mxGamma);
	int T = mxGetN(mxGamma);

	if ((mxGetM(mxResid) != K) || (mxGetN(mxResid) != T) ) {
		ss.clear();
		ss << "got size(gamma) => (" << K << "," << T << ")" << endl;
		ss << "bad size of resid, expected (" << K << "," << T << "), got ("<< mxGetM(mxResid) << "," << mxGetN(mxResid)  << "), try again";
		mexErrMsgTxt(ss.str().c_str());
	}


	/* get input matrix */
	o_gamma = mxGetPr(mxGamma);
	o_resid = mxGetPr(mxResid);

	try {
		// change gamma from column to row  representation (transpose gamma)
		gsl_matrix_view gammaTv = gsl_matrix_view_array(o_gamma, T, K);
		gsl_matrix *gamma = helper.get_managed_gsl_matrix(K, T);
		helper.matrixTranspose(&gammaTv.matrix, gamma);

		// change prev_resid from column to row  representation (transpose prev_resid)
		gsl_matrix_view residTv = gsl_matrix_view_array(o_resid, T, K);
		gsl_matrix *resid = helper.get_managed_gsl_matrix(K, T);
		helper.matrixTranspose(&residTv.matrix, resid);

		double acf;

		MCMCMet *met = new MCMCMetGPD();

		plhs[0] = mxCreateDoubleScalar(met->computeAcf(gamma, resid));

	} catch (exception& e) {
		mexErrMsgTxt(e.what());
	}
}

/*
 * MEX_GEV_FNTestimate.cpp
 *
 *  Created on: Aug 14, 2012
 *      Author: igdalov, changed to adpativ by kaisero
 */
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <exception>
#include "mex.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "GSL_Helper.h"
#include "MCMCMet.h"
#include "MCMCMetGPD.h"
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	GSL_Helper helper;
	stringstream ss;

	if (nrhs != 17) {
		ss.clear();
		ss << "17 inputs required, got " << nrhs << ", try again";
		mexErrMsgTxt(ss.str().c_str());
	}
	if (nlhs != 1) {
		ss.clear();
		ss << "one outputs required. got " << nlhs << ", try again";
		mexErrMsgTxt(ss.str().c_str());
	}
	double *o_xt, *o_ut_ksi, *o_ut_si, *o_gamma, *o_prev_resid, *o_startpar;

	// mex matrizen
	const mxArray *mxGamma, *mxUtwo_ksi, *mxUtwo_si, *mxXt, *mxResiPrev, *mxStartpar;

	/*  create a pointer to the input matrix */
	mxGamma = prhs[0];
	mxResiPrev = prhs[1];
	mxXt = prhs[2];
	mxUtwo_ksi = prhs[3];
	mxUtwo_si = prhs[4];
	mxStartpar = prhs[5];

	/* get the cfg stuff TODO: change to struct */
	int samples = mxGetScalar(prhs[6]);
	int proposMaxTimes = mxGetScalar(prhs[7]);
	double noise = mxGetScalar(prhs[8]);
	double noise_more = mxGetScalar(prhs[9]);
	double noise_less = mxGetScalar(prhs[10]);
	double beta = mxGetScalar(prhs[11]);
	double beta_factor = mxGetScalar(prhs[12]);
	int seed = mxGetScalar(prhs[13]);
	bool gaussFLAG = mxIsLogicalScalarTrue(prhs[14]);
	double lassoFLAG = mxGetScalar(prhs[15]);
	double ridgeFLAG = mxGetScalar(prhs[16]);

	/* get dimensions */
	int K = mxGetN(mxGamma);
	int tBinsN = mxGetM(mxGamma);
    
	int dimUt_ksi = mxGetN(mxUtwo_ksi);
	int dimUt_si  = mxGetN(mxUtwo_si);
	int len_param = dimUt_ksi + dimUt_si;

	if (mxGetM(mxUtwo_ksi) != tBinsN) {
		ss.clear();
		ss << "got size(gamma) => (" << K << "," << tBinsN << ")," <<" size(ut_ksi) => (" << dimUt_ksi << "," << mxGetM(mxUtwo_ksi)  << ")" << endl;
		ss << "bad length of ut_ksi, length(ut_ksi) != length(gamma), expected " << tBinsN << ", got " << mxGetM(mxUtwo_ksi) << ", try again";
		mexErrMsgTxt(ss.str().c_str());
	}
	if (mxGetM(mxUtwo_si) != tBinsN) {
		ss.clear();
		ss << "got size(gamma) => (" << K << "," << tBinsN << ")," <<" size(ut_si) => (" << dimUt_si << "," << mxGetM(mxUtwo_si)  << ")" << endl;
		ss << "bad length of ut_si, length(ut_si) != length(gamma), expected " << tBinsN << ", got " << mxGetM(mxUtwo_si) << ", try again";
		mexErrMsgTxt(ss.str().c_str());
	}
	if ((mxGetM(mxXt) != tBinsN) || (mxGetN(mxXt) != 1)) {
		ss.clear();
		ss << "got size(gamma) => (" << K << "," << tBinsN << ")," << endl;
		ss << "bad size of xt, expected (1," << tBinsN << "), got (" << mxGetN(mxXt) << "," << mxGetM(mxXt) << "), try again";
		mexErrMsgTxt(ss.str().c_str());
	}
	if ((mxGetN(mxStartpar) != (len_param)) || (mxGetM(mxStartpar) != K)) {
		ss.clear();
		ss << "got size(gamma) => (" << K << "," << tBinsN << ")," <<" size(ut_all) => (" << len_param << "," << mxGetM(mxUtwo_ksi)  << ")" << endl;
		ss << "bad size of startPar, expected (" << K << "," << len_param << ") " << ", got (" << mxGetM(mxStartpar) << "," << mxGetN(mxStartpar) << "), try again";
		mexErrMsgTxt(ss.str().c_str());
	}
	int RT = mxGetM(mxResiPrev);
    
    
  
    
	if ((RT != 0) && ((RT != tBinsN) || (mxGetN(mxResiPrev) != K))) {
		ss.clear();
		ss << "got size(gamma) => (" << K << "," << tBinsN << ")," <<" size(ut_all) => (" << len_param << "," << mxGetM(mxUtwo_ksi)  << ")" << endl;
		ss << "bad size of prev_resid, expected ("<< K << "," << tBinsN << ") or []" << ", got (" << mxGetN(mxResiPrev) << "," << RT << "), try again";
		mexErrMsgTxt(ss.str().c_str());
	}

	/* get input matrix */
	o_gamma = mxGetPr(mxGamma);
	o_prev_resid = mxGetPr(mxResiPrev);
	o_xt = mxGetPr(mxXt);
	o_ut_ksi = mxGetPr(mxUtwo_ksi);
	o_ut_si  = mxGetPr(mxUtwo_si);
	o_startpar = mxGetPr(mxStartpar);

	gsl_matrix* prev_resid = 0;

	try {
        // gamma
		gsl_matrix_view gamma_v = gsl_matrix_view_array(o_gamma, K, tBinsN);
        gsl_matrix* Gamma = helper.get_managed_gsl_matrix(tBinsN, K);
        helper.matrixTranspose(&gamma_v.matrix, Gamma);

		// prev resid may be empty
		if (RT != 0) {
			gsl_matrix_view prev_resid_v = gsl_matrix_view_array(o_prev_resid, K, tBinsN);
            prev_resid = helper.get_managed_gsl_matrix(tBinsN, K);
            helper.matrixTranspose(&prev_resid_v.matrix, prev_resid);
		}

		// xt is just a vector
		gsl_vector_view xt_v = gsl_vector_view_array(o_xt, tBinsN);

		// ut 
        gsl_matrix_view utT_v_ksi = gsl_matrix_view_array(o_ut_ksi, dimUt_ksi, tBinsN);
        gsl_matrix* utT_ksi = helper.get_managed_gsl_matrix(tBinsN, dimUt_ksi);
        helper.matrixTranspose(&utT_v_ksi.matrix, utT_ksi);
        
        gsl_matrix_view utT_v_si = gsl_matrix_view_array(o_ut_si, dimUt_si, tBinsN);
        gsl_matrix* utT_si = helper.get_managed_gsl_matrix(tBinsN, dimUt_si);
        helper.matrixTranspose(&utT_v_si.matrix, utT_si);
        
		// startpar 
		gsl_matrix_view startPar_v = gsl_matrix_view_array(o_startpar, len_param, K);
        gsl_matrix* startPar = helper.get_managed_gsl_matrix(K, len_param);
        helper.matrixTranspose(&startPar_v.matrix, startPar);

		// alloc memory for result parameters and result resid
		gsl_matrix* resultPAR = helper.get_managed_gsl_matrix(K, len_param);
		gsl_matrix* resultRESID = helper.get_managed_gsl_matrix(tBinsN, K);

//   testing block
//  				cout << " gamma  " << endl;
//  				helper.printMatrix(Gamma);
//  				cout  << endl;
// 
// 				if (prev_resid!= 0){
// 					cout << " prev_resid  " << endl;
// 					helper.printMatrix(prev_resid);
// 				}else{
// 					cout << "prev resid empty" << endl;
// 				}
// 				cout << " xt vector " << endl;
// 				helper.printVector(&xt_v.vector);
// 				cout <<  endl;
// 				cout << " ut_ksi  " << endl;
// 				helper.printMatrix(utT_ksi);
// 				cout  << endl;
// 				cout << " ut_si " << endl;
// 				helper.printMatrix(utT_si);
// 				cout <<  endl;
//     			cout << " startpar right way  " << endl;
// 				helper.printMatrix(startPar);
// 				cout << endl;
// 				cout << "samples " << samples << endl;
// 				cout << "proposMaxTimes " << proposMaxTimes << endl;
// 				cout << "noise " << noise << endl;
// 				cout << "noise_more " << noise_more << endl;
// 				cout << "noise_less " << noise_less << endl;
// 				cout << "beta " << beta << endl;
// 				cout << "beta_factor " << beta_factor << endl;
// 				cout << "seed " << seed << endl;

		double resultNRG;

		MCMCMet *met = new MCMCMetGPD(samples, proposMaxTimes, noise, noise_more, noise_less, beta, beta_factor, seed, gaussFLAG, lassoFLAG, ridgeFLAG);

		bool valid = met->estimate(Gamma, prev_resid, &xt_v.vector, utT_ksi, utT_si, startPar, resultPAR, resultRESID, &resultNRG);
   
		if (valid) {            
			//---------------------------------------------------------------------
			// mxCreateDoubleMatrix initialize the memory with 0
			// this is a performanze issue
			// to optimize this http://www.mathworks.com/matlabcentral/newsreader/view_thread/108658

			/*  set the output pointer to the output matrix */
			plhs[0] = mxCreateDoubleMatrix(resultPAR->size1, resultPAR->size2, mxREAL);
			//plhs[1] = mxCreateDoubleMatrix(resultRESID->size1, resultRESID->size2, mxREAL);
			//plhs[2] = mxCreateDoubleScalar(1);

			/*  create a C pointer to a copy of the output matrix */
			double *p_parameter = mxGetPr(plhs[0]);
			//double *p_resid = mxGetPr(plhs[1]);
			//double *p_resultNRG = mxGetPr(plhs[2]);

			// bring the output in the col representation back (for matlab)
			gsl_matrix_view resultPARv = gsl_matrix_view_array(p_parameter, resultPAR->size2, resultPAR->size1);
			//gsl_matrix_view resultRESIDv = gsl_matrix_view_array(p_resid, resultRESID->size2, resultRESID->size1);
			helper.matrixTranspose(resultPAR, &resultPARv.matrix);
			//helper.matrixTranspose(resultRESID, &resultRESIDv.matrix);
			//*p_resultNRG = resultNRG;

		} else {
			// in case of no result we should return empty parameters and empty resid
			//*p_parameter = [];//*p_resid = [];//*p_resultNRG = inf;
			/*  set the output pointer to the output matrix */
			plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
			//plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
			//plhs[2] = mxCreateDoubleScalar(1);

			/*  create a C pointer to a copy of the output matrix */
			//double *p_resultNRG = mxGetPr(plhs[2]);
			//*p_resultNRG = mxGetInf();

		}
	} catch (exception& e) {
		mexErrMsgTxt(e.what());
	}
}

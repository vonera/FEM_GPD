/*
 * MEX_GEV_FNT_resid.cpp
 *
 *  Created on: Aug 14, 2012
 *      Author: kaisero
 */
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <exception>
#include "mex.h"
#include <math.h>
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "GSL_Helper.h"
using namespace std;


//-----------------------------------------------------------------------------
void computeResid(gsl_matrix* param, gsl_vector* xt, gsl_matrix* utT_ksi, gsl_matrix* utT_si, gsl_matrix* result) {

	GSL_Helper helper;
	double dimUt_ksi = utT_ksi->size2;
	double dimUt_si  = utT_si ->size2;
	gsl_matrix_set_all(result, 3e+10);

	for (unsigned int t = 0; t < result->size1; t++) {

		double xtt = gsl_vector_get(xt, t);
		gsl_vector_view ut_tt_ksi = gsl_matrix_row(utT_ksi, t);
		gsl_vector_view ut_tt_si  = gsl_matrix_row(utT_si,  t);

		for (unsigned int i = 0; i < result->size2; i++) {

			// TODO: optimize this, it should be done only ones instead of tBinsN times
			gsl_vector_view iparam = gsl_matrix_row(param, i);

			// estimate ksi(ut) si(ut)
			gsl_vector_view tmp = gsl_vector_subvector(&iparam.vector, 0, dimUt_ksi);
			double ksi = helper.dotProd(&tmp.vector, &ut_tt_ksi.vector);

			tmp = gsl_vector_subvector(&iparam.vector, dimUt_ksi, dimUt_si);
			double si = helper.dotProd(&tmp.vector, &ut_tt_si.vector);

			if (ksi != 0) {
				double zt = 1 + ksi * xtt / si;
				if (zt > 0) {
					gsl_matrix_set(result, t, i, log(si) + (1 + 1 / ksi) * log(zt));
				}
			} else{
				gsl_matrix_set(result, t, i, log(si) + xtt / si);
			}
		}
	}
}
//-----------------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    GSL_Helper helper;
    stringstream ss;

    if (nrhs != 4) {
        ss.clear();
        ss << "4 inputs required, got " << nrhs << ", try again";
        mexErrMsgTxt(ss.str().c_str());
    }
    if (nlhs != 1) {
        ss.clear();
        ss << "one output required. got " << nlhs << ", try again";
        mexErrMsgTxt(ss.str().c_str());
    }
    double *o_param, *o_xt, *o_ut_ksi, *o_ut_si;

    // mex matrizen
    const mxArray *mxParam, *mxXt, *mxUt_ksi, *mxUt_si;

    /*  create a pointer to the input matrix */
    mxParam = prhs[0];
    mxXt = prhs[1];
    mxUt_ksi = prhs[2];
    mxUt_si  = prhs[3];
 
    /* get dimensions */
    int K = mxGetM(mxParam);
    int tBinsN = mxGetM(mxXt);
    int dimUt_ksi = mxGetN(mxUt_ksi);
    int dimUt_si = mxGetN(mxUt_si);
     int len_param = dimUt_ksi + dimUt_si;

    /* get input matrix */
    o_param = mxGetPr(mxParam);
    o_xt = mxGetPr(mxXt);
    o_ut_ksi = mxGetPr(mxUt_ksi);
    o_ut_si  = mxGetPr(mxUt_si);


    try {
        // change gamma from column to row  representation (transpose gamma)
        gsl_matrix_view paramTv = gsl_matrix_view_array(o_param,len_param, K);
        gsl_matrix *param = helper.get_managed_gsl_matrix(K, len_param);
        helper.matrixTranspose(&paramTv.matrix, param);

        // xt is just a vector
        gsl_vector_view xtv = gsl_vector_view_array(o_xt, tBinsN);

        // ut
        gsl_matrix_view utT_v_ksi = gsl_matrix_view_array(o_ut_ksi, dimUt_ksi, tBinsN);
        gsl_matrix* utT_ksi = helper.get_managed_gsl_matrix(tBinsN, dimUt_ksi);
        helper.matrixTranspose(&utT_v_ksi.matrix, utT_ksi);
        
        gsl_matrix_view utT_v_si = gsl_matrix_view_array(o_ut_si, dimUt_si, tBinsN);
        gsl_matrix* utT_si = helper.get_managed_gsl_matrix(tBinsN, dimUt_si);
        helper.matrixTranspose(&utT_v_si.matrix, utT_si);
       

        // alloc memory for result parameters and result resid
        gsl_matrix* resultRESID = helper.get_managed_gsl_matrix(tBinsN, K);


        computeResid(param, &xtv.vector, utT_ksi, utT_si, resultRESID);

        /*  set the output pointer to the output matrix */
        plhs[0] = mxCreateDoubleMatrix(resultRESID->size1, resultRESID->size2, mxREAL);

        /*  create a C pointer to a copy of the output matrix */
        double *p_resid = mxGetPr(plhs[0]);

        // bring the output in the col representation back (for matlab)

        gsl_matrix_view resultRESIDv = gsl_matrix_view_array(p_resid, resultRESID->size2, resultRESID->size1);
        helper.matrixTranspose(resultRESID, &resultRESIDv.matrix);
        
        
        helper.free_managed(param);
        helper.free_managed(utT_ksi);
        helper.free_managed(utT_si);
        helper.free_managed(resultRESID);
        
    } catch (exception& e) {
        mexErrMsgTxt(e.what());
    }
}



/*
 * MCMCMetGevFNT.cpp
 *
 *  Created on: Feb 28, 2013
 *      Author: kaisero
 *      This version contains the fully non stationary GPD MCMC parameter estimation
 */
#include <math.h>
#include <list>
#include <iostream>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include "GSL_Helper.h"
#include "MCMCMetGPD.h"

using namespace std;

//-----------------------------------------------------------------------------
void MCMCMetGPD::computeResid(const gsl_matrix* Param, const gsl_vector* xt, const gsl_matrix* UtT_ksi, const gsl_matrix* UtT_si, gsl_matrix* Result) {

	GSL_Helper helper;
	double dimUt_ksi = UtT_ksi->size2;
	double dimUt_si  = UtT_si ->size2;
	gsl_matrix_set_all(Result, 3e+10);

	for (unsigned int t = 0; t < Result->size1; t++) {

		double xtt = gsl_vector_get(xt, t);
		gsl_vector_const_view utt_ksi_v = gsl_matrix_const_row(UtT_ksi, t);
		gsl_vector_const_view utt_si_v  = gsl_matrix_const_row(UtT_si,  t);

		for (unsigned int i = 0; i < Result->size2; i++) {

			// TODO: optimize this, it should be done only ones instead of T times
			gsl_vector_const_view iparam_v = gsl_matrix_const_row(Param, i);

			// estimate ksi(ut) si(ut) mu(ut)
			gsl_vector_const_view tmpksi_v = gsl_vector_const_subvector(&iparam_v.vector, 0, dimUt_ksi);
			double ksi = helper.dotProd(&tmpksi_v.vector, &utt_ksi_v.vector);

			gsl_vector_const_view tmpsi_v = gsl_vector_const_subvector(&iparam_v.vector, dimUt_ksi, dimUt_si);
			double si = helper.dotProd(&tmpsi_v.vector, &utt_si_v.vector);

			if (ksi != 0) {
				double zt = 1 + ksi * xtt / si;
				if (zt > 0) {
					gsl_matrix_set(Result, t, i, log(si)  + (1 + 1 / ksi) * log( zt ));
				}
			} else {
				gsl_matrix_set(Result, t, i, log(si) + xtt / si);
			}
		}
	}
}


//-----------------------------------------------------------------------------
void MCMCMetGPD::simulate(gsl_matrix *Param, gsl_matrix *Gamma, gsl_matrix *Start, gsl_matrix* Ut_ksi, gsl_matrix* Ut_si, gsl_matrix* Result) {
	// TODO implement
}


//-----------------------------------------------------------------------------
/**
 * utT is ut Transposed,
 */
bool MCMCMetGPD::proveConstraints(const gsl_matrix* param, const gsl_matrix* gamma, const gsl_vector* xt,  const gsl_matrix* utT_ksi, const gsl_matrix* utT_si) {

	GSL_Helper helper;
	bool result;
	double dimUt_ksi = utT_ksi->size2;
	double dimUt_si  = utT_si ->size2;

	// check constrains on ksi and si:
	if (!proveConstraints_ksi(param, utT_ksi) ) {
		return result = false;
	}
	if ( !proveConstraints_si(param, utT_si, dimUt_ksi) ) {
		return result = false;
	}

	// prepare for cheking contrant 1+ksi/si(xt-mu) > 0
	gsl_matrix* utTcluster_i_full_ksi = helper.get_managed_gsl_matrix(utT_ksi->size1, utT_ksi->size2);
	gsl_matrix* utTcluster_i_full_si  = helper.get_managed_gsl_matrix(utT_si ->size1, utT_si ->size2);
	gsl_vector* xtCluster_i_full = helper.get_managed_gsl_vector(xt->size);

	for (unsigned int i = 0; i < param->size1; i++) {

		gsl_vector_const_view igamma = gsl_matrix_const_column(gamma, i);
		double max_igamma = gsl_vector_max (&igamma.vector);
http://www.myswitzerland.com/de-ch/events/event-203120649.html
		if(max_igamma > 0){

			result = false;

			// get all parameters for the cluster i
			gsl_vector_const_view iparam = gsl_matrix_const_row(param, i);
			gsl_vector_const_view iparam_ksi = gsl_vector_const_subvector(&iparam.vector, 0, dimUt_ksi);
			gsl_vector_const_view iparam_si  = gsl_vector_const_subvector(&iparam.vector, dimUt_ksi, dimUt_si);

			int next_free = 0;
			for (unsigned int t = 0; t < gamma->size1; t++) {
				// here are both variants save, TODO check the performance of both.
				//if (gsl_vector_get(&gamma_row.vector , t)  > 0) {
				if (igamma.vector.data[t * igamma.vector.stride] > 0) {
					helper.copyMatrixRowToMatrixRow(utT_ksi, t, utTcluster_i_full_ksi, next_free);
					helper.copyMatrixRowToMatrixRow(utT_si,  t, utTcluster_i_full_si,  next_free);
					// TODO test performance of all three
					//xtG1_full->data[next_free] = xt->data[t];
					//xtG1_full->data[next_free*xtG1_full->stride] = xt->data[t*xt->stride];
					gsl_vector_set(xtCluster_i_full, next_free, gsl_vector_get(xt, t));
					next_free += 1;
				}
			}

			//TODO check performance of that, wouldn't it be faster to put everithing in fo(if()) ?
			// reduce utT_cluster_i_full_* and xtCluster_i_full to the actual size of next_free
			gsl_matrix_view utT_cluster_i_ksi = gsl_matrix_submatrix(utTcluster_i_full_ksi, 0, 0, next_free, utTcluster_i_full_ksi->size2);
			gsl_matrix_view utT_cluster_i_si  = gsl_matrix_submatrix(utTcluster_i_full_si,  0, 0, next_free, utTcluster_i_full_si ->size2);
			gsl_vector_view xt_zt  = gsl_vector_subvector(xtCluster_i_full, 0, next_free);

			// initialize ksi(ut), si(ut), mu(ut) for this cluster
			gsl_vector* ksi_ut = helper.get_managed_gsl_vector(next_free);
			gsl_vector* si_ut  = helper.get_managed_gsl_vector(next_free);

			// instead of  "xt_cluster_i * ut_cluster_i" "utT_cluster_i * xt_cluster_i"
			helper.matrixTimesVector(&utT_cluster_i_ksi.matrix, &iparam_ksi.vector, ksi_ut);
			helper.matrixTimesVector(&utT_cluster_i_si.matrix,  &iparam_si.vector,  si_ut);

			// compute xt_zt = (1+ksi/si*xt)
			helper.divideVector(si_ut, ksi_ut);
			helper.VectorDotTimesVector(ksi_ut, &xt_zt.vector);
			helper.addConstVector(&xt_zt.vector, 1);

			helper.free_managed(ksi_ut);
			helper.free_managed(si_ut);

			// must be xt_zt > 0
			if (gsl_vector_min(&xt_zt.vector) < 0) {
				break;
			}
		}

		result = true;
	}

	helper.free_managed(xtCluster_i_full);
	helper.free_managed(utTcluster_i_full_ksi);
	helper.free_managed(utTcluster_i_full_si);
	return result;
}



//-----------------------------------------------------------------------------
bool MCMCMetGPD::proveConstraints_ksi(const gsl_matrix* Param,  const gsl_matrix* UtT_ksi) {

	GSL_Helper helper;
	bool result;
	double dimUt_ksi = UtT_ksi->size2;

	for (unsigned int i = 0; i < Param->size1; i++) {

		result = false;
		// get all parameters for the cluster i
		gsl_vector_const_view iparam = gsl_matrix_const_row(Param, i);
		gsl_vector_const_view iparam_ksi = gsl_vector_const_subvector(&iparam.vector, 0, dimUt_ksi);
		gsl_vector* ksi_ut = helper.get_managed_gsl_vector(UtT_ksi->size1);

		//because ut is transposed we can do instead of  "iparam_ksiT * Ut_ksi" "UtT_ksi * iparam_ksiT"
		helper.matrixTimesVector(UtT_ksi, &iparam_ksi.vector, ksi_ut);


		//  check, must be:  -0.5 < ksi(ut) < 0.5
		if ( gsl_vector_max( ksi_ut ) > 0.5 || gsl_vector_min( ksi_ut) < -0.5 ) {
			helper.free_managed(ksi_ut);
			break;
		}
		helper.free_managed(ksi_ut);
		result = true;
	}


	return result;
}


//-----------------------------------------------------------------------------
bool MCMCMetGPD::proveConstraints_si(const gsl_matrix* Param,  const gsl_matrix* UtT_si, int dimUt_ksi) {

	GSL_Helper helper;
	bool result;
	double dimUt_si = UtT_si->size2;

	for (unsigned int i = 0; i < Param->size1; i++) {

		result = false;
		// get si(ut) parameters for the cluster i
		gsl_vector_const_view iparam_v = gsl_matrix_const_row(Param, i);
		gsl_vector_const_view iparam_si_v = gsl_vector_const_subvector(&iparam_v.vector, dimUt_ksi, dimUt_si);
		gsl_vector* si_ut = helper.get_managed_gsl_vector(UtT_si->size1);

		//because ut is transposed we can do instead of  "iparam_siT * ut_si" "UtT_si * iparam_siT"
		helper.matrixTimesVector(UtT_si, &iparam_si_v.vector, si_ut);

		//  check, must be: si(ut) > 0
		if ( gsl_vector_min( si_ut ) <= 0.0 ) {
			helper.free_managed(si_ut);
			break;
		}
		helper.free_managed(si_ut);
		result = true;
	}

	return result;
}


//-----------------------------------------------------------------------------
/**-----------------------------------------
 * If the number of accepted parameter is small (<= 2*length(param))
 * then use next ~ N(ilast,0.1^2/length(param) * ID) for each cluster i independently:
 */
void MCMCMetGPD::fillNextNormal(const gsl_matrix* Last, const double noise, gsl_matrix* Next) {

	GSL_Helper helper;
	unsigned int K = Last->size1;
	unsigned int len_param = Last->size2;

	gsl_vector *r_norm = helper.get_managed_gsl_vector(len_param);
	gsl_matrix *iCov   = helper.get_managed_gsl_matrix(len_param, len_param);

	gsl_matrix_set_identity(iCov);
	helper.scaleMatrix(iCov, noise * 0.01/len_param);

	for (unsigned int i = 0; i < K ; i++) {

		gsl_vector_const_view ilast_v = gsl_matrix_const_row(Last, i);
		helper.mvnrnd(r_seed, &ilast_v.vector, iCov, r_norm);

		gsl_matrix_view Next_row_v =  gsl_matrix_view_vector(r_norm, 1, len_param);
		helper.copyMatrixRowToMatrixRow(&Next_row_v.matrix, 0, Next, i);
	}

	helper.free_managed(r_norm);
	helper.free_managed(iCov);
}


//-----------------------------------------------------------------------------
/**-----------------------------------------
 * If the number of accepted parameter is big enough (>2*length(param))
 * use Gaussian Distribution N(mean, cov) for each cluster i independently:
 * next = (1-noise) * N(ilast, 5.6644/len_param * icovMatrix) + noise*N(ilast,0.1^2/length(param) * covID)
 * NOTE: here if the parameter noise -> 0 the acceptance is decreasing,
 *       and as noise -> 1 the acceptance is increasing.
 */
void MCMCMetGPD::fillNextGauss(gsl_matrix **CovCell, const gsl_matrix* Last, const double noise, gsl_matrix* Next) {

	GSL_Helper helper;
	unsigned int K = Last->size1;
	unsigned int len_param = Last->size2;

	gsl_vector *r_gaus = helper.get_managed_gsl_vector(len_param);
	gsl_vector *r_norm = helper.get_managed_gsl_vector(len_param);
	gsl_vector *inext  = helper.get_managed_gsl_vector(len_param);
	gsl_matrix *iCov   = helper.get_managed_gsl_matrix(len_param, len_param);
	gsl_matrix *id_Cov = helper.get_managed_gsl_matrix(len_param, len_param);

	gsl_matrix_set_identity(id_Cov);
	helper.scaleMatrix(id_Cov, 0.01/len_param);

	for (unsigned int i = 0; i < K ; i++) {
		helper.copyMatrix(CovCell[i], iCov);
		helper.scaleMatrix(iCov, 5.6644 / len_param);

		gsl_vector_const_view ilast_v = gsl_matrix_const_row(Last, i);
		helper.mvnrnd(r_seed, &ilast_v.vector, iCov, r_gaus);
		helper.mvnrnd(r_seed, &ilast_v.vector, id_Cov, r_norm);

		helper.av_plus_bu(1-noise, r_gaus, noise, r_norm, inext);
		gsl_matrix_view Next_row_v =  gsl_matrix_view_vector(inext, 1, len_param);
		helper.copyMatrixRowToMatrixRow(&Next_row_v.matrix, 0, Next, i);
	}

	helper.free_managed(r_gaus);
	helper.free_managed(r_norm);
	helper.free_managed(inext);
	helper.free_managed(iCov);
	helper.free_managed(id_Cov);
}


//-----------------------------------------------------------------------------
/**-----------------------------------------
 * This function propose the next parameter with Gaussian sampling involving the covariance
 * matrix of all accepted parameters. The implementation of the algorithm is according to
 * [1] Optimal Proposal Distributions and Adaptive MCMC, Jeffrey S. Rosenthal,
 * [2] Adaptive Optimal Scaling of MH Algo Using the Robbins-Monro Process, P.H. GARTHWAITE & Y.FAN
 */
void MCMCMetGPD::proposeNext( const gsl_matrix* Last, const gsl_matrix* Gamma, const gsl_vector* xt, const gsl_matrix* Ut_ksi, const gsl_matrix* Ut_si, double noise, gsl_matrix **CovCell, const unsigned int countAcc, gsl_matrix* Next ) {

	GSL_Helper helper;

	bool found = false;
	unsigned int len_param = Last->size2;

	for (int i = 0; i < proposMaxTimes; i++) {

		// propose next parameter
		if ( gaussFLAG && countAcc > 2*len_param ) {
			try{
				fillNextGauss(CovCell, Last, noise, Next);
			}
			catch(eOperationFaild& e){
				fillNextNormal(Last, noise, Next);
			}
		}else{
			fillNextNormal(Last, noise, Next);
		}

		// proveCostraints
		found = proveConstraints(Next, Gamma, xt, Ut_ksi, Ut_si);
		if (found) {
			break;
		}
	}
	if (!found) {
		helper.copyMatrix(Last, Next);
	}
}


//-----------------------------------------------------------------------------
double MCMCMetGPD::computeEnergy(const gsl_matrix* Param, const gsl_vector* xt, const gsl_matrix* UtT_ksi, const gsl_matrix* UtT_si, const gsl_matrix* Gamma, gsl_matrix* Result) {

	GSL_Helper helper;
	computeResid(Param, xt, UtT_ksi, UtT_si, Result);


	return helper.sumMatrixDotTimesMatix(Result, Gamma);
}

/*
 * MCMCMet.cpp
 *
 *  Created on: 23.08.2012
 *      Author: kaisero
 */

#include "MCMCMet.h"
#include <algorithm>
#include <math.h>
#include <iostream>
#include <list>
using namespace std;

//-----------------------------------------------------------------------------
MCMCMet::MCMCMet() :
		initialized(false), r_seed(0), samples(0), proposMaxTimes(0), initial_noise(0), noise_more(0), noise_less(0), initial_beta(0), beta_factor(0), gaussFLAG(0), lassoFLAG(0), ridgeFLAG(0) {
}


//-----------------------------------------------------------------------------
MCMCMet::MCMCMet(int samples, int proposMaxTimes, double noise, double noise_more, double noise_less, double beta, double beta_factor, unsigned long int seed, bool gaussFLAG, double lassoFLAG, double ridgeFLAG) :
		initialized(true), samples(samples), proposMaxTimes(proposMaxTimes), initial_noise(noise), noise_more(noise_more), noise_less(noise_less), initial_beta(
				beta), beta_factor(beta_factor), gaussFLAG(gaussFLAG), lassoFLAG(lassoFLAG), ridgeFLAG(ridgeFLAG) {
	const gsl_rng_type * T;
	T = gsl_rng_mt19937;
	r_seed = gsl_rng_alloc(T);
	gsl_rng_set(r_seed, seed);
//	cout.precision(15);
//	for (int i = 0; i < 100000; i++){
//		cout << gsl_rng_uniform(r) << endl;
//	}

}


//-----------------------------------------------------------------------------
MCMCMet::~MCMCMet() {

	if (initialized) {
		//free random
		gsl_rng_free(r_seed);
	}
}


//-----------------------------------------------------------------------------
bool MCMCMet::estimate(gsl_matrix* gamma, gsl_matrix* prev_resid, gsl_vector* xt, gsl_matrix* ut_ksi, gsl_matrix* ut_si, gsl_matrix* startpar, gsl_matrix* resultPAR, gsl_matrix* resultRESID, double *p_resultNRG) {

	GSL_Helper helper;
	double noise = initial_noise;
	double beta  = initial_beta;
	unsigned long int acceptet_times = 0;
	int countAcc = 0;
	int K = gamma->size2;

	// initialize the covMatrix for proposing the next parameter
	gsl_matrix *cov_temp   = helper.get_managed_gsl_matrix(startpar->size2 , startpar->size2);
	gsl_vector *mean_temp  = helper.get_managed_gsl_vector(startpar->size2);
	gsl_matrix_set_zero(cov_temp);
	gsl_vector_set_zero(mean_temp);

	gsl_matrix **covCell  = new gsl_matrix*[K];
	gsl_vector **meanCell = new gsl_vector*[K];
	for(int i = 0; i < K; i++){
		covCell[i]  = helper.get_managed_gsl_matrix(startpar->size2 , startpar->size2);
		meanCell[i] = helper.get_managed_gsl_vector(startpar->size2);
		helper.copyMatrix(cov_temp,  covCell[i]);
		helper.copyVector(mean_temp, meanCell[i]);
	}

	// TODO here we can check the constraints on the startpars and if they are not good we can simply return false.
	helper.copyMatrix(startpar, resultPAR);
	*p_resultNRG = computeEnergy(resultPAR, xt, ut_ksi, ut_si, gamma, resultRESID);
	*p_resultNRG = shrinkageEnergyUpdate(resultPAR, ut_ksi->size2, ut_si->size2, *p_resultNRG);

	bool result_valide = true;
	double initialNRG;

	if (prev_resid == 0) {
		initialNRG = *p_resultNRG;
	} else {
		double tmp = helper.sumMatrixDotTimesMatix(gamma, prev_resid);
		initialNRG = shrinkageEnergyUpdate(resultPAR,ut_ksi->size2, ut_si->size2, tmp);
	}

	// TODO check if >= here is better because now we accept startParameters that are better or EQUALLY good
	if (*p_resultNRG >= initialNRG) {

		result_valide = false;
		gsl_matrix* proposedPar = helper.get_managed_gsl_matrix(resultPAR->size1, resultPAR->size2);
		double proposedNRG;
		
		for (int sample = 1; sample < samples; sample++) {

			proposeNext(resultPAR, gamma, xt, ut_ksi, ut_si, noise, covCell, countAcc, proposedPar);
			proposedNRG = computeEnergy(proposedPar, xt, ut_ksi, ut_si, gamma, resultRESID);
            proposedNRG = shrinkageEnergyUpdate(proposedPar, ut_ksi->size2, ut_si->size2, proposedNRG);

			if (accept(beta, *p_resultNRG, proposedNRG)) {

				*p_resultNRG = proposedNRG;
				helper.copyMatrix(proposedPar,resultPAR);
				acceptet_times += 1;

				// update matrix and mean
				countAcc += 1;
				if(countAcc > 2){
					updateCovMean(meanCell, covCell, resultPAR, countAcc, K);
				}

				// TODO <=??
				if (*p_resultNRG < initialNRG) {
					result_valide = true;
					break;
				}
			}
			if (sample >= 50) {
				adapt(&noise, &beta, acceptet_times, sample);
			}
		}
		helper.free_managed(proposedPar);
	}

	free(covCell);
	free(meanCell);
	helper.free_managed(cov_temp);
	helper.free_managed(mean_temp);

	return result_valide;
}

//-----------------------------------------------------------------------------
double MCMCMet::computeAcf(gsl_matrix* Gamma, gsl_matrix* Resid){
	GSL_Helper helper;
	return helper.sumMatrixDotTimesMatix(Resid, Gamma);
}


//-----------------------------------------------------------------------------
void MCMCMet::adapt(double* noise, double* beta, int acceptet_times, int number_stepps) {

	if (number_stepps % 50 == 0) {

		double acceptance_rate = (double)acceptet_times / number_stepps;
		//cout << "acceptance_rate " << acceptance_rate << endl;

		if (acceptance_rate > 0.7) {

			*noise = *noise * noise_more;
			//cout << "reset noise more" << *noise << endl;
		} else if (acceptance_rate < 0.3) {

			*noise = *noise * noise_less;
			//cout << "reset noise less" << *noise << endl;
		}
	}
	if (number_stepps % 100 == 0) {

		*beta = *beta * beta_factor;
		//cout << "increase beta" << endl;
	}
}


//-----------------------------------------------------------------------------
/**
 * all comments copied from o.kaisers matlab implementation
 */
bool MCMCMet::accept(double beta, double acceptedNRG, double proposedNRG) {

	bool accepted = false;

	// this is enough "gsl_rng_uniform(r) <= min(1.0, exp(beta * (acceptedNRG - proposedNRG )))"
	// but we can shorten the whole thing by first checking the acceptedNRG >= proposedNRG  condition.
	if ( acceptedNRG >= proposedNRG  ||  gsl_rng_uniform(r_seed) <= min(1.0, exp(beta * (acceptedNRG - proposedNRG )))) {

		accepted = true;
	}

	return accepted;
}


//--------------------------------------------------------------------------------------
void MCMCMet::updateCovMean(gsl_vector **meanCell, gsl_matrix **CovCell, const gsl_matrix *ResultPAR, const int countAcc, const int K){

	GSL_Helper helper;
	gsl_vector * mean_temp = helper.get_managed_gsl_vector(ResultPAR->size2);

	for( int i = 0; i < K; i++){

		helper.copyVector(meanCell[i], mean_temp );
		gsl_vector_const_view itemp = gsl_matrix_const_row(ResultPAR, i);

		// update mean:  ((countAcc-1)*meanP + acceptedPar(i,:)') / countAcc;
		helper.av_plus_bu( 1.0*(countAcc-1.0) / countAcc ,  mean_temp, 1.0 / countAcc, &itemp.vector, meanCell[i] );

		//	update covMatrix
		//  cov_new = (countAcc-2)/(countAcc-1) * covM + meanP*meanP' - countAcc * (meanP_new*meanP_new') / (countAcc-1) + acceptedPar(i,:)'*acceptedPar(i,:) /(countAcc-1);
		gsl_matrix *vmTvmT = helper.get_managed_gsl_matrix(mean_temp->size,mean_temp->size);
		gsl_matrix *vm_newTvm_newT = helper.get_managed_gsl_matrix(mean_temp->size,mean_temp->size);
		gsl_matrix *v_accParamTv_accParamT = helper.get_managed_gsl_matrix(mean_temp->size,mean_temp->size);

		helper.VTimesVT(mean_temp, mean_temp, vmTvmT);
		helper.VTimesVT(meanCell[i], meanCell[i], vm_newTvm_newT);
		helper.VTimesVT(&itemp.vector, &itemp.vector, v_accParamTv_accParamT);

		helper.scaleMatrix( CovCell[i], ( countAcc-2.0 ) / ( countAcc-1.0 ));
		helper.addMatrix(vmTvmT, CovCell[i]);

		helper.scaleMatrix( vm_newTvm_newT, (-1.0*countAcc) / (countAcc -1.0) );
		helper.addMatrix(vm_newTvm_newT, CovCell[i]);

		helper.scaleMatrix( v_accParamTv_accParamT, 1.0 / (countAcc - 1.0) );
		helper.addMatrix(v_accParamTv_accParamT, CovCell[i]);

		helper.free_managed(vmTvmT);
		helper.free_managed(vm_newTvm_newT);
		helper.free_managed(v_accParamTv_accParamT);
	}

	helper.free_managed(mean_temp);

}




//-------------------------------------------------------------------------
/**
 * This function is used to update the energy wrt shrinakge: Lasso (sum |bi|) or Ridgde (sum (bi)^2) 
 * only on the coefficients of the covas,
 * Lasso/Ridge are included via Lagrangian with the lagrange parameter: lassoFLAG/ridgeFLAG.
 */
double MCMCMet::shrinkageEnergyUpdate(const gsl_matrix* Param, int dimUt_ksi, int dimUt_si, const double energy){

	GSL_Helper helper;
	double shrinkage_energy = 0;
		
	for (unsigned int i = 0; i < Param->size1; i++) {	
		gsl_vector_const_view iparam_v = gsl_matrix_const_row(Param, i);
		gsl_vector_const_view iparam_ksi_v = gsl_vector_const_subvector(&iparam_v.vector, 0, dimUt_ksi);
		gsl_vector_const_view iparam_si_v = gsl_vector_const_subvector(&iparam_v.vector, dimUt_ksi, dimUt_si);
		double sum_ksi = 0;
		double sum_si = 0;
			
		if (lassoFLAG != 0){
			for (int j = 1; j < dimUt_ksi; j++){
				sum_ksi = sum_ksi + fabs(gsl_vector_get(&iparam_ksi_v.vector, j));
			}
			for (int j = 1; j < dimUt_si; j++){
				sum_si = sum_si + fabs(gsl_vector_get(&iparam_si_v.vector, j));
			}
			shrinkage_energy = shrinkage_energy + lassoFLAG*(sum_ksi + sum_si);
	 	}else if(ridgeFLAG != 0){
			for (int j = 1; j < dimUt_ksi; j++){
				sum_ksi = sum_ksi + pow(gsl_vector_get(&iparam_ksi_v.vector, j), 2);
			}
			for (int j = 1; j < dimUt_si; j++){
				sum_si = sum_si + pow(gsl_vector_get(&iparam_si_v.vector, j), 2);
			}
			shrinkage_energy = shrinkage_energy + ridgeFLAG*(sum_ksi + sum_si);
		}
	}
	
	

	return energy + shrinkage_energy;
}









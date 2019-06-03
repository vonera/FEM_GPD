/*
 * MCMCMet.h
 *
 *  Created on: Feb 28, 2013
 *      Author: kaisero
 */

#ifndef MCMCMet_H_
#define MCMCMet_H_
#include "GSL_Helper.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_blas.h"

class MCMCMet {

private:

protected:
	bool initialized;
	gsl_rng *r_seed;
	int samples;
	int proposMaxTimes;
	double initial_noise;
	double noise_more;
	double noise_less;
	double initial_beta;
	double beta_factor;
	double lassoFLAG;
	double ridgeFLAG;
	bool gaussFLAG;

public:
	MCMCMet();
	MCMCMet(int samples, int proposMaxTimes, double noise, double noise_more, double noise_less, double beta, double beta_factor, unsigned long int seed, bool gaussFLAG, double lassoFLAG, double ridgeFLAG);

	bool estimate(gsl_matrix* gamma, gsl_matrix* prev_resid, gsl_vector* xt, gsl_matrix* ut_ksi, gsl_matrix* ut_si, gsl_matrix* startpar, gsl_matrix* resultPAR, gsl_matrix* resultRESID, double *resultNRG);

	double computeAcf(gsl_matrix* gamma, gsl_matrix* resid);

	virtual void computeResid(const gsl_matrix *param, const gsl_vector *xt, const gsl_matrix* ut_ksi, const gsl_matrix* ut_si, gsl_matrix* result) = 0;
	virtual void simulate(gsl_matrix *param, gsl_matrix *gamma, gsl_matrix *start, gsl_matrix* ut_ksi, gsl_matrix* ut_si, gsl_matrix* result) = 0;
	virtual ~MCMCMet();
//protected:
	void adapt(double *noise, double *beta, int acceptet_times, int number_stepps);
	bool accept(double beta, double acceptedNRG, double proposedNRG);
	void updateCovMean(gsl_vector **meanCell, gsl_matrix **covCell, const gsl_matrix *resultPAR, const int countAcc, const int K);
	virtual void proposeNext(const gsl_matrix* last, const gsl_matrix* gamma, const gsl_vector* xt, const gsl_matrix* ut_ksi, const gsl_matrix* ut_si, const double noise, gsl_matrix **paramCov, const unsigned int countAcc, gsl_matrix* next) = 0;
	virtual double computeEnergy(const gsl_matrix *param, const gsl_vector *xt, const gsl_matrix* ut_ksi, const gsl_matrix* ut_si, const gsl_matrix *gamma, gsl_matrix* resid) = 0;
	double shrinkageEnergyUpdate(const gsl_matrix* Param, int dimUt_ksi, int dimUt_si, double energy);
};


#endif /* MCMCMet_H_ */

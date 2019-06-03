/*
 * MCMCMetGPD.h
 *
 *  Created on: Feb 28, 2013
 *      Author: kaisero
 */

#ifndef MCMCMETGPD_H_
#define MCMCMETGPD_H_

#include "MCMCMet.h"

class MCMCMetGPD: public MCMCMet {

public:
	MCMCMetGPD() :
		MCMCMet() {
	}
	MCMCMetGPD(int samples, int proposMaxTimes, double noise, double noise_more, double noise_less, double beta, double beta_factor, unsigned long int seed, bool gaussFLAG, double lassoFLAG, double ridgeFLAG) :
		MCMCMet(samples, proposMaxTimes, noise, noise_more, noise_less, beta, beta_factor, seed, gaussFLAG, lassoFLAG, ridgeFLAG) {
	}
	virtual ~MCMCMetGPD(){};

	void computeResid(const gsl_matrix* Param, const gsl_vector* xt, const gsl_matrix* UtT_ksi, const gsl_matrix* UtT_si,  gsl_matrix* Result);
	void simulate(gsl_matrix *Param, gsl_matrix *Gamma, gsl_matrix *Start, gsl_matrix* Ut_ksi, gsl_matrix* Ut_si, gsl_matrix* Result);
//protected:
	bool proveConstraints(const gsl_matrix* Param, const gsl_matrix* Gamma, const gsl_vector* xt,  const gsl_matrix* UtT_ksi, const gsl_matrix* UtT_si);
	bool proveConstraints_ksi(const gsl_matrix* Param,  const gsl_matrix* UtT_ksi);
	bool proveConstraints_si(const gsl_matrix* Param,  const gsl_matrix* UtT_si, int dimUt_ksi);
	void proposeNext(const gsl_matrix* Last, const gsl_matrix* Gamma, const gsl_vector* xt, const gsl_matrix* Ut_ksi, const gsl_matrix* Ut_si, const double noise, gsl_matrix **CovCell, const unsigned int countAcc, gsl_matrix* Next );
	double computeEnergy(const gsl_matrix* Param, const gsl_vector* xt, const gsl_matrix* UtT_ksi, const gsl_matrix* UtT_si,  const gsl_matrix* Gamma, gsl_matrix* Resid);
	void fillNextNormal(const gsl_matrix* Last, const double noise, gsl_matrix* Next);
	void fillNextGauss(gsl_matrix **CovCell, const gsl_matrix* Last, const double noise, gsl_matrix* Next);

};

#endif /* MCMCMETGPD_H_ */

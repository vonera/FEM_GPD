/*
 * GLSHelper.h
 *
 *  Created on: 24.08.2012
 *      Author: igdalovd
 */

#ifndef GSLHELPER_H_
#define GSLHELPER_H_
#include <exception>
#include <list>
#include <sstream>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_rng.h"


class eGSLException: public std::exception {
private:
	const char* text;
	int errnum;

public:

	eGSLException(const char* text, int errnum) {
		this->text = text;
		this->errnum = errnum;
	}

	virtual const char* what() const throw () {
		std::string message(text);
        std::stringstream ss;
        
		message.append(" : errornum : ");
        ss << errnum;
        message.append(ss.str());
        
        //message.append(static_cast<std::ostringstream*>(&(std::ostringstream() << errnum))->str());

		return message.c_str();
	}
};
class eAllocationError: public eGSLException {
public:
	eAllocationError(const char* text, int errnum) :
			eGSLException(text, errnum) {
	}
};
class eOperationFaild: public eGSLException {
public:
	eOperationFaild(const char* text, int errnum) :
			eGSLException(text, errnum) {
	}
};


class GSL_Helper {
private:
	std::list<gsl_matrix*> my_gsl_matrices;
	std::list<gsl_vector*> my_gsl_vectors;
public:
	//-----------------------------------------------------------------------------
	GSL_Helper();
	virtual ~GSL_Helper();

	//-----------------------------------------------------------------------------
	void printMatrix(const gsl_matrix* toprint) throw ();

	void printVector(const gsl_vector* toprint) throw ();

	void printVector(const gsl_vector *toprint, bool transpose) throw ();

	void switch_pointers(void **p1, void **p2) throw ();

	//-----------------------------------------------------------------------------
	gsl_vector * get_managed_gsl_vector(size_t) throw (eAllocationError);

	gsl_matrix * get_managed_gsl_matrix(size_t, size_t) throw (eAllocationError);

	gsl_matrix* get_managed_gsl_matrix(gsl_matrix* to_copy_from) throw (eAllocationError, eOperationFaild);

	gsl_matrix* get_managed_gsl_matrix_zeros(size_t, size_t) throw (eAllocationError);

	void free_managed(gsl_matrix *M) throw ();

	void free_managed(gsl_vector *v) throw ();

	void unmanage(gsl_matrix *M) throw ();

	void unmanage(gsl_vector *v) throw ();

	//-----------------------------------------------------------------------------
	double sumMatrixDotTimesMatix(const gsl_matrix *A, const gsl_matrix *B) throw (eAllocationError, eOperationFaild);

	double dotProd(const gsl_vector *v, const gsl_vector *u) throw (eOperationFaild);

	void matrixTimesVector(const gsl_matrix *A, const gsl_vector *v, gsl_vector* result) throw (eOperationFaild);

	void av_plus_bu(double alpha, const gsl_vector *v, double beta, const gsl_vector *u, gsl_vector* result) throw (eOperationFaild);

	void aM_plus_bA(double alpha, const gsl_matrix *M, double beta, const gsl_matrix *A, gsl_matrix *result) throw (eAllocationError, eOperationFaild);

	void matrixTranspose(const gsl_matrix *A, gsl_matrix *result) throw (eOperationFaild);

	void copyMatrix(const gsl_matrix * from, gsl_matrix * to) throw (eOperationFaild);

	void copyVector(const gsl_vector * from, gsl_vector * to) throw (eOperationFaild);

	void copyMatrixRowToMatrixRow(const gsl_matrix *src, int srow, gsl_matrix *dest, int drow) throw (eOperationFaild);

	void copySubMatrixToMatrix(const gsl_matrix *src, gsl_matrix* dest, int ifrom, int jfrom);

	void copySubMatrixToSubMatrix(const gsl_matrix *src, gsl_matrix* dest, int si, int sj, int di, int dj, int nrow, int ncol) throw (eOperationFaild);

	void copySubVectorToVector1(const gsl_vector *src, gsl_vector* dest, int from, int to, int num);

	void copySubVectorToVector2(const gsl_vector *src, gsl_vector* dest, int from, int to, int num) throw (eOperationFaild);

	void copySubVectorToVector3(const gsl_vector *src, gsl_vector* dest, int from, int to, int num) throw (eOperationFaild);

	void copySubVectorToVector4(const gsl_vector *src, gsl_vector* dest, int from);

	void addVector(const gsl_vector * thisOne, gsl_vector * toThat) throw (eOperationFaild);

	void addMatrix(const gsl_matrix * thisOne, gsl_matrix * toThat) throw (eOperationFaild);

	void scaleVector(gsl_vector * toScale, double by) throw (eOperationFaild);

	void scaleMatrix(gsl_matrix * toScale, double by) throw (eOperationFaild);

	void addConstVector(gsl_vector * vec, double constant) throw (eOperationFaild);

	void subVector(const gsl_vector * that, gsl_vector * fromThat) throw (eOperationFaild);

	void divideVector(const gsl_vector *that, gsl_vector *byThat) throw (eOperationFaild);

	void VectorDotTimesVector(const gsl_vector *that, gsl_vector *withThat) throw (eOperationFaild);

	void MatrixTimesMatrix(const gsl_matrix* A, const gsl_matrix* M, gsl_matrix* result) throw (eOperationFaild);

	void VTimesVT(const gsl_vector* v, const gsl_vector* u, gsl_matrix* result) throw (eOperationFaild);

	void mvnrnd(const gsl_rng *r, const gsl_vector *mean, const gsl_matrix *cov, gsl_vector *result) throw (eOperationFaild);
	//----------------------------------------------------------------------------------------------
};

#endif /* GSLHELPER_H_ */

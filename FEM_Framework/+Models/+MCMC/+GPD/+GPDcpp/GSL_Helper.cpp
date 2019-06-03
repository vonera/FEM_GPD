/*
 * GLSHelper.cpp
 *
 *  Created on: 24.08.2012
 *      Author: igdalovd
 */

#include "GSL_Helper.h"
#include <string.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>

#include <iostream>
using namespace std;

//-----------------------------------------------------------------------------
GSL_Helper::GSL_Helper() {
}

//-----------------------------------------------------------------------------
GSL_Helper::~GSL_Helper() {

	list<gsl_matrix*>::iterator itm;
	for (itm = my_gsl_matrices.begin(); itm != my_gsl_matrices.end(); itm++) {
		if (*itm != 0) {
			gsl_matrix_free(*itm);
		}
	}
	my_gsl_matrices.clear();

	list<gsl_vector*>::iterator itv;
	for (itv = my_gsl_vectors.begin(); itv != my_gsl_vectors.end(); itv++) {
		if (*itv != 0) {
			gsl_vector_free(*itv);
		}
	}
	my_gsl_vectors.clear();
}

//-----------------------------------------------------------------------------
void GSL_Helper::printMatrix(const gsl_matrix* toprint) throw () {
	cout << endl;
	for (unsigned int row = 0; row < toprint->size1; row++) {
		for (unsigned int col = 0; col < toprint->size2; col++) {
			cout << gsl_matrix_get(toprint, row, col) << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

//-----------------------------------------------------------------------------
void GSL_Helper::printVector(const gsl_vector* toprint) throw () {
	printVector(toprint, false);
}

//-----------------------------------------------------------------------------
void GSL_Helper::printVector(const gsl_vector* toprint, bool transpose) throw () {
	cout << endl;
	const char *terminal = "\t";
	if (!transpose) {
		terminal = "\n";
	}
	for (unsigned int i = 0; i < toprint->size; i++) {
		cout << gsl_vector_get(toprint, i) << terminal;
	}
	cout << endl;
}

//-----------------------------------------------------------------------------
void GSL_Helper::switch_pointers(void **p1, void **p2) throw () {
	void *ptr = *p1;
	*p1 = *p2;
	*p2 = ptr;
}

//-----------------------------------------------------------------------------
gsl_vector* GSL_Helper::get_managed_gsl_vector(size_t size) throw (eAllocationError) {

	gsl_vector *result;

	if ((result = gsl_vector_alloc(size)) == 0) {
		throw eAllocationError("Error in GSL_Helper::get_managed_gsl_vector : gsl_vector_alloc", 0);
	}
	my_gsl_vectors.push_back(result);

	return result;
}

//-----------------------------------------------------------------------------
gsl_matrix* GSL_Helper::get_managed_gsl_matrix(gsl_matrix* to_copy_from) throw (eAllocationError, eOperationFaild) {

	gsl_matrix* result;
	if ((result = gsl_matrix_alloc(to_copy_from->size1, to_copy_from->size2)) == 0) {
		throw eAllocationError("Error in GSL_Helper::get_managed_gsl_matrix : gsl_matrix_alloc", 0);
	}
	my_gsl_matrices.push_back(result);

	int err = gsl_matrix_memcpy(result, to_copy_from);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::get_managed_gsl_matrix : gsl_matrix_memcpy", err);
	}

	return result;
}

//-----------------------------------------------------------------------------
gsl_matrix* GSL_Helper::get_managed_gsl_matrix(size_t rows, size_t cols) throw (eAllocationError) {

	gsl_matrix* result;
	if ((result = gsl_matrix_alloc(rows, cols)) == 0) {
		throw eAllocationError("Error in GSL_Helper::get_managed_gsl_matrix : gsl_matrix_alloc", 0);
	}
	my_gsl_matrices.push_back(result);

	return result;
}

//-----------------------------------------------------------------------------
gsl_matrix* GSL_Helper::get_managed_gsl_matrix_zeros(size_t rows, size_t cols) throw (eAllocationError) {

	gsl_matrix* result;
	if ((result = gsl_matrix_calloc(rows, cols)) == 0) {
		throw eAllocationError("Error in GSL_Helper::get_managed_gsl_matrix_zeros : gsl_matrix_calloc", 0);
	}
	my_gsl_matrices.push_back(result);

	return result;
}


//-----------------------------------------------------------------------------
void GSL_Helper::free_managed(gsl_matrix* M) throw () {
	my_gsl_matrices.remove(M);
	if (M != 0) {
		gsl_matrix_free(M);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::free_managed(gsl_vector* v) throw () {
	my_gsl_vectors.remove(v);
	if (v != 0) {
		gsl_vector_free(v);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::unmanage(gsl_matrix* M) throw () {
	my_gsl_matrices.remove(M);
}

//-----------------------------------------------------------------------------
void GSL_Helper::unmanage(gsl_vector* v) throw () {
	my_gsl_vectors.remove(v);
}

//-----------------------------------------------------------------------------
double GSL_Helper::sumMatrixDotTimesMatix(const gsl_matrix* A, const gsl_matrix* B) throw (eAllocationError, eOperationFaild) {

	gsl_matrix *tmp = get_managed_gsl_matrix(A->size1, A->size2);

	int err = gsl_matrix_memcpy(tmp, A);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::sumMatrixDotTimesMatix : gsl_matrix_memcpy", err);
	}

	err = gsl_matrix_mul_elements(tmp, B);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::sumMatrixDotTimesMatix : gsl_matrix_mul_elements", err);
	}
	double result = 0;


	//slow
	for (unsigned int i = 0; i < tmp->size1; i++){
		for (unsigned int j = 0; j < tmp->size2; j++){
			result += gsl_matrix_get(tmp, i,j);
		}
	}
	//  fast
	//for (unsigned int i = 0; i < (tmp->size1 * tmp->size2); i++) {
	//	result += tmp->data[i];
	//}

	free_managed(tmp);

	return result;
}

//-----------------------------------------------------------------------------
double GSL_Helper::dotProd(const gsl_vector* v, const gsl_vector* u) throw (eOperationFaild) {

	double result;
	int err = gsl_blas_ddot(v, u, &result);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::dotProd : gsl_blas_ddot", err);
	}
	return result;
}

//-----------------------------------------------------------------------------
/**
 * result is completely overwritten
 * */
void GSL_Helper::matrixTimesVector(const gsl_matrix* A, const gsl_vector* v, gsl_vector* result) throw (eOperationFaild) {
	int err = gsl_blas_dgemv(CblasNoTrans, 1.0, A, v, 0.0, result);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::matrixTimesVector : gsl_blas_dgemv", err);
	}
}

//-----------------------------------------------------------------------------
/**
 * result is completely overwritten
 * */
void GSL_Helper::av_plus_bu(double alpha, const gsl_vector* v, double beta, const gsl_vector* u, gsl_vector* result) throw (eOperationFaild) {

	int err = gsl_vector_memcpy(result, u);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::av_plus_bu : gsl_vector_memcpy", err);
	}

	if (beta != 1.0) {
		err = gsl_vector_scale(result, beta);
		if (err) {
			throw eOperationFaild("Error in GSL_Helper::av_plus_bu : gsl_vector_scale", err);
		}
	}

	err = gsl_blas_daxpy(alpha, v, result);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::av_plus_bu : gsl_blas_daxpy", err);
	}
}

//-----------------------------------------------------------------------------
/**
 * result is completely overwritten
 * */
void GSL_Helper::aM_plus_bA(double alpha, const gsl_matrix* M, double beta, const gsl_matrix* A, gsl_matrix* result) throw (eAllocationError, eOperationFaild) {

	int err = gsl_matrix_memcpy(result, M);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::aM_plus_bA : gsl_matrix_memcpy", err);
	}
	err = gsl_matrix_scale(result, alpha);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::aM_plus_bA : gsl_matrix_scale", err);
	}
	if (beta == 1.0) {

		err = gsl_matrix_add(result, A);
		if (err) {
			throw eOperationFaild("Error in GSL_Helper::aM_plus_bA : gsl_matrix_add", err);
		}
	} else {
		gsl_matrix *tmp = get_managed_gsl_matrix(A->size1, A->size2);

		int err = gsl_matrix_memcpy(tmp, A);
		if (err) {
			throw eOperationFaild("Error in GSL_Helper::aM_plus_bA : gsl_matrix_memcpy", err);
		}

		err = gsl_matrix_scale(tmp, beta);
		if (err) {
			throw eOperationFaild("Error in GSL_Helper::aM_plus_bA : gsl_matrix_scale", err);
		}

		err = gsl_matrix_add(result, tmp);
		if (err) {
			throw eOperationFaild("Error in GSL_Helper::aM_plus_bA : gsl_matrix_add", err);
		}

		free_managed(tmp);
	}
}

//-----------------------------------------------------------------------------
/**
 * result is completely overwritten
 * */
void GSL_Helper::matrixTranspose(const gsl_matrix* A, gsl_matrix* result) throw (eOperationFaild) {

	int err = gsl_matrix_transpose_memcpy(result, A);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::matrixTranspose : gsl_matrix_transpose_memcpy", err);
	}

}

//-----------------------------------------------------------------------------
void GSL_Helper::copyMatrix(const gsl_matrix * from, gsl_matrix * to) throw (eOperationFaild) {
	int err = gsl_matrix_memcpy(to, from);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::copyMatrix : gsl_matrix_memcpy", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::copyVector(const gsl_vector * from, gsl_vector * to) throw (eOperationFaild) {
	int err = gsl_vector_memcpy(to, from);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::copyVector : gsl_vector_memcpy", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::copyMatrixRowToMatrixRow(const gsl_matrix *src, int srow, gsl_matrix *dest, int drow) throw (eOperationFaild) {
	gsl_vector_view dest_row = gsl_matrix_row(dest, drow);
	int err = gsl_matrix_get_row(&dest_row.vector, src, srow);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::copyMatrixRowToMatrixRow : gsl_matrix_get_row", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::copySubMatrixToMatrix(const gsl_matrix *src, gsl_matrix* dest, int ifrom, int jfrom){

	int N = src->size1;
	int M = src->size2;

	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			gsl_matrix_set(dest, ifrom + i, jfrom + j, gsl_matrix_get(src, i, j) );
		}
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::copySubMatrixToSubMatrix(const gsl_matrix *src, gsl_matrix* dest, int si, int sj, int di, int dj, int nrow, int ncol) throw (eOperationFaild) {
	gsl_matrix_const_view subSrc =  gsl_matrix_const_submatrix (src, si, sj, nrow, ncol);
	if (subSrc.matrix.data == 0){
		throw eOperationFaild("Error in GSL_Helper::copySubMatrixToSubMatrix : gsl_matrix_const_submatrix returned null", -1);
	}
	gsl_matrix_view subDest  = gsl_matrix_submatrix (dest, di, dj, nrow, ncol);
	if (subDest.matrix.data == 0){
		throw eOperationFaild("Error in GSL_Helper::copySubMatrixToSubMatrix : gsl_matrix_const_submatrix returned null", -1);
	}
	int err = gsl_matrix_memcpy(&subDest.matrix, &subSrc.matrix);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::copySubMatrixToSubMatrix : gsl_matrix_memcpy", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::copySubVectorToVector1(const gsl_vector *src, gsl_vector* dest, int from, int to, int num){
	memcpy(&dest->data[to], &src->data[from], sizeof(double)*num);
}

//-----------------------------------------------------------------------------
void GSL_Helper::copySubVectorToVector2(const gsl_vector *src, gsl_vector* dest, int from, int to, int num) throw (eOperationFaild) {
	if (dest->size < to + num || src->size < from + num){
		throw eOperationFaild("Error in GSL_Helper::copySubVectorToVector2 : bad size", -1);
	}
	memcpy(&dest->data[to], &src->data[from], sizeof(double)*num);
}

//-----------------------------------------------------------------------------
void GSL_Helper::copySubVectorToVector3(const gsl_vector *src, gsl_vector* dest, int from, int to, int num) throw (eOperationFaild) {
	gsl_vector_view d = gsl_vector_subvector(dest, to, num);
	gsl_vector_const_view s = gsl_vector_const_subvector(src, from, num);
	int err = gsl_vector_memcpy(&d.vector, &s.vector);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::copySubVectorToVector3 : gsl_vector_memcpy", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::copySubVectorToVector4(const gsl_vector *src, gsl_vector* dest, int from){

	for(unsigned int i = 0; i < src->size; i++){
		gsl_vector_set(dest, from + i, gsl_vector_get(src, i) );
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::addVector(const gsl_vector * thisOne, gsl_vector * toThat) throw (eOperationFaild) {

	int err = gsl_vector_add(toThat, thisOne);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::addVector : gsl_vector_add", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::scaleVector(gsl_vector * toScale, double by) throw (eOperationFaild) {

	int err = gsl_vector_scale(toScale, by);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::scaleVector : gsl_vector_scale", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::addMatrix(const gsl_matrix * thisOne, gsl_matrix * toThat) throw (eOperationFaild) {

	int err = gsl_matrix_add(toThat, thisOne);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::addMatrix : gsl_matrix_add", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::scaleMatrix(gsl_matrix * toScale, double by) throw (eOperationFaild) {

	int err = gsl_matrix_scale(toScale, by);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::scaleMatrix : gsl_matrix_scale", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::addConstVector(gsl_vector * vec, double constant) throw (eOperationFaild) {
	int err = gsl_vector_add_constant(vec, constant);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::addConstVector : gsl_vector_add_constant", err);
	}
}

//-----------------------------------------------------------------------------
// result in fromThat
void GSL_Helper::subVector(const gsl_vector * that, gsl_vector * fromThat) throw (eOperationFaild) {
	int err = gsl_vector_sub(fromThat, that);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::subVector : gsl_vector_sub", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::divideVector(const gsl_vector * that, gsl_vector * byThat) throw (eOperationFaild) {

	int err = gsl_vector_div(byThat, that);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::v_divide_u : gsl_vector_div", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::VectorDotTimesVector(const gsl_vector* that, gsl_vector* withThat) throw (eOperationFaild) {

	int err = gsl_vector_mul(withThat, that);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::v_divide_u : gsl_vector_div", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::MatrixTimesMatrix(const gsl_matrix* A, const gsl_matrix* M, gsl_matrix* result) throw (eOperationFaild) {

	int err = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, M, 0.0, result);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::matrixTimesMatrix : gsl_blas_dgemm", err);
	}
}

//-----------------------------------------------------------------------------
void GSL_Helper::VTimesVT(const gsl_vector* v, const gsl_vector* u, gsl_matrix* result) throw (eOperationFaild) {

	gsl_matrix_const_view A_view = gsl_matrix_const_view_array(v->data, v->size, 1);
	gsl_matrix_const_view B_view = gsl_matrix_const_view_array(u->data, 1, u->size);
	const gsl_matrix *A = &A_view.matrix;
	const gsl_matrix *B = &B_view.matrix;

	// A :=v, B := uT
	int err = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, result);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::matrixTimesMatrix : gsl_blas_dgemm", err);
	}
}

//-----------------------------------------------------------------------------
/** -----------------------------------------------------------------------------
 *
 * r is the seed for random number generator
 * mean is the mean value of Gaussian
 * cov is the covariance matrix of the Gaussian
 */
void GSL_Helper::mvnrnd(const gsl_rng *r, const gsl_vector *mean, const gsl_matrix *cov, gsl_vector *result) throw (eOperationFaild) {

	gsl_matrix* work_matrix;
	if ((work_matrix = gsl_matrix_alloc(mean->size, mean->size)) == 0) {
		throw eAllocationError("Error in GSL_Helper::mvnrnd(), gsl_matrix_alloc", 0);
	}

	int err = gsl_matrix_memcpy(work_matrix, cov);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::mvnrnd(), gsl_matrix_memcpy", err);
	}


	gsl_error_handler_t * pDefaultHandler = gsl_set_error_handler_off() ;

	// get cholesky decomposition
	err = gsl_linalg_cholesky_decomp(work_matrix);
	if (err) {
		throw eOperationFaild("Error in GSL_Helper::mvnrnd(), gsl_linalg_cholesky_decomp", err);
	}

	gsl_set_error_handler (pDefaultHandler);


	// get uniform vector of size mean->size
	 for(unsigned int k = 0; k < mean->size; k++){
	    gsl_vector_set( result, k, gsl_ran_ugaussian(r) );
	 }


	 // transform wrt cholesky
	 err = gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work_matrix, result);
	 if (err) {
	 	 throw eOperationFaild("Error in GSL_Helper::mvnrnd(), gsl_blas_dtrmv", err);
	 }

	 // add mean value
	 err = gsl_vector_add(result, mean);
	 if (err) {
		throw eOperationFaild("Error in GSL_Helper::mvnrnd(), gsl_vector_add", err);
	 }

	 gsl_matrix_free(work_matrix);
}

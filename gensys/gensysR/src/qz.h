/*
 * Copyright 2023 James M. Murray
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * If you have any questions about this software, contact James Murray
 * at james@murraylax.org
*/

/**
 * @file qz.h
 *
 * @brief 
 * This is a QZ decomposition for matrices defined by the Eigen package. 
 * It is a wrapper for the C function, LAPACKE_zgges(), in the lapacke package.
 * 
 * @author James Murray
 * Contact: james@murraylax.org
*/


#ifndef _INCL_QZ
#define _INCL_QZ

#include <complex>
#include <lapacke.h>
#include <Eigen/Dense>


/**
 * @brief Eigenvalue selection function passed to LAPACKE_zgges.
 *
 * Returns true (selects) eigenvalues with modulus >= 1 (unstable roots),
 * causing LAPACK to order them last in the generalized Schur decomposition.
 *
 * @param alpha Numerator of the generalized eigenvalue
 * @param beta  Denominator of the generalized eigenvalue
 * @return Non-zero if the eigenvalue should be placed in the unstable block
 */
extern "C" lapack_logical eigenvalue_threshold_function(const lapack_complex_double* alpha, const lapack_complex_double* beta);

/**
 * @brief Compute the QZ (generalized Schur) decomposition of the matrix pencil (A, B).
 *
 * Wraps LAPACKE_zgges to produce the decomposition \f$ A = Q S Z^H \f$, \f$ B = Q T Z^H \f$,
 * with eigenvalues ordered so that stable roots (modulus < 1) appear first.
 *
 * @param mcdQ      (output) Unitary matrix Q
 * @param mcdZ      (output) Unitary matrix Z
 * @param mcdS      (output) Upper quasi-triangular matrix S
 * @param mcdT      (output) Upper quasi-triangular matrix T
 * @param vdLambda  (output) Generalized eigenvalue moduli |alpha/beta|
 * @param mdA       (input)  Real matrix A
 * @param mdB       (input)  Real matrix B
 * @return Number of stable eigenvalues (modulus < 1)
 */
int qz(Eigen::MatrixXcd& mcdQ, Eigen::MatrixXcd& mcdZ, Eigen::MatrixXcd& mcdS, Eigen::MatrixXcd& mcdT, Eigen::VectorXd& vdLambda,
    const Eigen::MatrixXd& mdA, const Eigen::MatrixXd& mdB);

/**
 * @brief Check if the column space of B is a subset of the column space of A.
 *
 * Uses SVD of A; tests whether projections of B's columns onto the left null space of A are zero.
 *
 * @param A Matrix whose column space is the reference
 * @param B Matrix whose columns are tested
 * @return true if col(B) ⊆ col(A) within a default tolerance of 1e-10
 */
bool check_colspace(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B);

#endif
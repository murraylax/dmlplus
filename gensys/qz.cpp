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
 * @file qz.cpp
 *
 * @brief 
 * This is a QZ decomposition for matrices defined by the Eigen package. 
 * It is a wrapper for the C function, LAPACKE_zgges(), in the lapacke package.
 * 
 * @author James Murray
 * Contact: james@murraylax.org
*/

#include <iostream>

#include "qz.h"
#include "utils.h"

using namespace std;
using namespace Eigen;

/**
 * extern "C" lapack_logical eigenvalue_threshold_function(const lapack_complex_double* alpha, const lapack_complex_double* beta)
 * 
 * Threshold function for LAPACKE_zgges() to identify eigenvalues that are less than 1.0
 * 
 * @param alpha Numerator of an eigenvalue from the QZ decomposition (lambda = alpha / beta)
 * @param beta Numerator of an eigenvalue from the QZ decomposition (lambda = alpha / beta)
*/
extern "C" lapack_logical eigenvalue_threshold_function(const lapack_complex_double* alpha, const lapack_complex_double* beta) {
    // Convert to std::complex for easy handling
    std::complex<double> eigenvalue = std::complex<double>(*beta) / std::complex<double>(*alpha);

    // Check the magnitude of the eigenvalue
    if (std::abs(eigenvalue) > 1.0) {
        return 0;
    } else {
        return 1;
    }
}

/**
 * qz(Eigen::MatrixXcd& mcdQ, Eigen::MatrixXcd& mcdZ, Eigen::MatrixXcd& mcdS, Eigen::MatrixXcd& mcdT, 
 *           Eigen::VectorXcd& vcdLambda, const Eigen::MatrixXd& mdA, const Eigen::MatrixXd& mdB) 
 * 
 * Compute a complex QZ decomposition, with eigenvalues > 1.0 appearing at the bottom of the decomposition
 * 
 * For given matrices A and B, a QZ decomposition returns matrices Q, Z, S, and T, and eigenvalues lambda, such that...
 * 
 * Q' A Z = S
 * Q' B Z = T
 * Q'Q = Z'Z = I
 * S and T are upper triangular
 * And the eigenvalues in lambda are magnitude of the diagonal elements of T divided by the magnitude of the diagonal elements of S
 * 
 * @param mcdQ Output parameter, the complex matrix Q in the decomposition above
 * @param mcdZ Output parameter, the complex matrix Z in the decomposition above
 * @param mcdS Output parameter, the complex matrix S in the decomposition above
 * @param mcdT Output parameter, the complex matrix T in the decomposition above
 * @param vdLambda Output parameter, the vector lambda in the decomposition above
 * @param mdA Input parameter, the matrix A above
 * @param mdB Input parameter, the matrix B above
 * 
 * @return Returns an integer equal to the number of stable eigenvalues (i.e. lambda_i < 1.0)
*/
int qz(Eigen::MatrixXcd& mcdQ, Eigen::MatrixXcd& mcdZ, Eigen::MatrixXcd& mcdS, Eigen::MatrixXcd& mcdT, Eigen::VectorXd& vdLambda,
    const Eigen::MatrixXd& mdA, const Eigen::MatrixXd& mdB) {

    // Check matrix dimensions
    int n = mdA.rows();
    int n1 = mdA.cols();
    int n2 = mdB.rows();
    int n3 = mdB.cols();
    if(n != n1 || n != n1 || n != n2 || n != n3) {
        std::cerr << "Error in qz_lapack(): Matrices are not square or do not have the same dimension." << endl;
        std::exit(EXIT_FAILURE);
    }

    // Copy Eigen::MatrixXd to lapack_complex_double
    lapack_complex_double* lcA = new lapack_complex_double[mdA.size()];
    lapack_complex_double* lcB = new lapack_complex_double[mdB.size()];

    copy_matrix_to_lapack_complex(lcA, mdA);
    copy_matrix_to_lapack_complex(lcB, mdB);

    char jobvsl = 'V'; // Compute left Schur vectors
    char jobvsr = 'V'; // Compute right Schur vectors

    char sort = 'S'; // Sort the eigenvalues.
    
    // Outputs
    lapack_int sdim;
    lapack_complex_double* alpha = new lapack_complex_double[n];
    lapack_complex_double* beta = new lapack_complex_double[n];
    lapack_complex_double* vsl = new lapack_complex_double[n * n];
    lapack_complex_double* vsr = new lapack_complex_double[n * n];

    // Workspace
    lapack_complex_double* work = new lapack_complex_double[3 * n]; // Workspace
    double* rwork = new double[8 * n]; // Real workspace for complex routines

    LAPACKE_zgges(LAPACK_COL_MAJOR, jobvsl, jobvsr, sort, eigenvalue_threshold_function, n, lcA, n, lcB, n, &sdim, alpha, beta, vsl, n, vsr, n);

    // Copy results
    copy_lapack_complex_to_matrix(mcdQ, vsl);
    copy_lapack_complex_to_matrix(mcdZ, vsr);
    copy_lapack_complex_to_matrix(mcdS, lcA); // A is overwritten to be S
    copy_lapack_complex_to_matrix(mcdT, lcB); // B is overwritten to be T

    int nstable = 0;
    for(int i=0; i<n; i++) {
        vdLambda(i) = std::abs( std::complex<double>(beta[i]) / std::complex<double>(alpha[i]) ); 
        if(vdLambda(i)<=1.0) nstable++;
    }

    // Clean up
    delete[] lcA;
    delete[] lcB;
    delete[] alpha;
    delete[] beta;
    delete[] vsl;
    delete[] vsr;
    delete[] work;
    delete[] rwork;

    return nstable;
}

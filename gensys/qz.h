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


extern "C" lapack_logical eigenvalue_threshold_function(const lapack_complex_double* alpha, const lapack_complex_double* beta);

int qz(Eigen::MatrixXcd& mcdQ, Eigen::MatrixXcd& mcdZ, Eigen::MatrixXcd& mcdS, Eigen::MatrixXcd& mcdT, Eigen::VectorXd& vdLambda,
    const Eigen::MatrixXd& mdA, const Eigen::MatrixXd& mdB);

#endif
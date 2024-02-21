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
 * @file utils.h
 *
 * @brief 
 * This is a collection of commonly used utilities for copying 
 * matrices and vectors from one format to another, including outputting 
 * text or R code to bring into R.
 * 
 * @author James Murray
 * Contact: james@murraylax.org
*/


#ifndef _INCL_UTILS
#define _INCL_UTILS

#include <complex>
#include <lapacke.h>
#include <Eigen/Dense>


void copy_matrix_to_lapack_complex(lapack_complex_double* lapackMatrix, const Eigen::MatrixXd& matrix);
void copy_lapack_complex_to_matrix(Eigen::MatrixXcd& complex_matrix, const lapack_complex_double* lapack_complex);
void copy_lapack_complex_to_vector(Eigen::VectorXcd& complex_vector, const lapack_complex_double* lapack_complex);
std::string print_vector_to_rcode(const Eigen::VectorXcd& v, const std::string& vecname);
std::string print_vector_to_rcode(const Eigen::VectorXd& v, const std::string& vecname);
std::string print_matrix_to_rcode(const Eigen::MatrixXcd& m, const std::string& matname);
std::string print_matrix_to_rcode(const Eigen::MatrixXd& m, const std::string& matname);

#endif
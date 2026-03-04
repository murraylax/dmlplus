/*
 * Copyright 2024 James M. Murray
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <lapacke.h>
#include <Eigen/Dense>
#include <chrono>

/**
 * @brief Write an Eigen::MatrixXd to a CSV file (one row per matrix row, comma-separated).
 * @param mat      Matrix to write
 * @param filepath Output file path
 */
void write_eigen_csv(const Eigen::MatrixXd& mat, const std::string& filepath);

/**
 * @brief Write an Eigen::VectorXd to a file with a variable name header line.
 * @param vec      Vector to write
 * @param varname  Variable name written on the first line
 * @param filepath Output file path
 */
void write_eigen_csv(Eigen::VectorXd& vec, std::string& varname, std::string& filepath);

/**
 * @brief Write an Eigen::VectorXd to a CSV file (one value per line).
 * @param vec      Vector to write
 * @param filepath Output file path
 */
void write_eigen_csv(const Eigen::VectorXd& vec, const std::string& filepath);

/**
 * @brief Write an Eigen::MatrixXd to a CSV file with a header row of variable names.
 * @param mat      Matrix to write
 * @param varnames Column variable names for the header row
 * @param filepath Output file path
 */
void write_eigen_csv(Eigen::MatrixXd& mat, std::vector<std::string>& varnames, std::string& filepath);

/// Type alias for steady-clock time points used by start_timer() / stop_timer().
using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;

/**
 * @brief Start a wall-clock timer.
 * @return The current time point to be passed to stop_timer()
 */
std::chrono::time_point<std::chrono::steady_clock> start_timer();

/**
 * @brief Stop a wall-clock timer and print the elapsed time to stdout.
 * @param start_time Time point returned by start_timer()
 */
void stop_timer(const TimePoint& start_time);

/**
 * @brief Copy an Eigen::VectorXd into a pre-allocated gsl_vector.
 * @param gsl_output   Destination gsl_vector (must be pre-allocated with the same size)
 * @param eigen_input  Source Eigen vector
 */
void copy_vector_to_gsl(gsl_vector* gsl_output, const Eigen::VectorXd& eigen_input);

/**
 * @brief Copy a gsl_vector into an Eigen::VectorXd.
 * @param eigen_output Destination Eigen vector (resized to match)
 * @param gsl_input    Source gsl_vector
 */
void copy_gsl_to_vector(Eigen::VectorXd& eigen_output, const gsl_vector* gsl_input);

/**
 * @brief Copy an Eigen::MatrixXd into a pre-allocated gsl_matrix.
 * @param gsl_output   Destination gsl_matrix (must be pre-allocated with the same dimensions)
 * @param eigen_input  Source Eigen matrix
 */
void copy_matrix_to_gsl(gsl_matrix* gsl_output, const Eigen::MatrixXd& eigen_input);

/**
 * @brief Copy a gsl_matrix into an Eigen::MatrixXd.
 * @param eigen_output Destination Eigen matrix (resized to match)
 * @param gsl_input    Source gsl_matrix
 */
void copy_gsl_to_matrix(Eigen::MatrixXd& eigen_output, const gsl_matrix* gsl_input);

/**
 * @brief Copy an Eigen::MatrixXd into a LAPACK column-major complex double array.
 * @param lapackMatrix Destination array (must be pre-allocated, size rows*cols)
 * @param matrix       Source real Eigen matrix
 */
void copy_matrix_to_lapack_complex(lapack_complex_double* lapackMatrix, const Eigen::MatrixXd& matrix);

/**
 * @brief Copy a LAPACK column-major complex double array into an Eigen::MatrixXcd.
 * @param complex_matrix Destination Eigen complex matrix
 * @param lapack_complex Source LAPACK array
 */
void copy_lapack_complex_to_matrix(Eigen::MatrixXcd& complex_matrix, const lapack_complex_double* lapack_complex);

/**
 * @brief Copy a LAPACK complex double array into an Eigen::VectorXcd.
 * @param complex_vector Destination Eigen complex vector
 * @param lapack_complex Source LAPACK array
 */
void copy_lapack_complex_to_vector(Eigen::VectorXcd& complex_vector, const lapack_complex_double* lapack_complex);

/**
 * @brief Generate an R assignment statement for a complex vector.
 * @param v       Source complex vector
 * @param vecname R variable name to use in the output string
 * @return String of the form "vecname <- c(...)"
 */
std::string print_vector_to_rcode(const Eigen::VectorXcd& v, const std::string& vecname);

/**
 * @brief Generate an R assignment statement for a real vector.
 * @param v       Source real vector
 * @param vecname R variable name to use in the output string
 * @return String of the form "vecname <- c(...)"
 */
std::string print_vector_to_rcode(const Eigen::VectorXd& v, const std::string& vecname);

/**
 * @brief Generate an R assignment statement for a complex matrix.
 * @param m       Source complex matrix
 * @param matname R variable name to use in the output string
 * @return String of the form "matname <- matrix(c(...), nrow=..., ncol=...)"
 */
std::string print_matrix_to_rcode(const Eigen::MatrixXcd& m, const std::string& matname);

/**
 * @brief Generate an R assignment statement for a real matrix.
 * @param m       Source real matrix
 * @param matname R variable name to use in the output string
 * @return String of the form "matname <- matrix(c(...), nrow=..., ncol=...)"
 */
std::string print_matrix_to_rcode(const Eigen::MatrixXd& m, const std::string& matname);

#endif
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
 * at jmurray@uwlax.edu
*/

/**
 * @file dml_multiroot.h
 *
 * @brief 
 * This is a numeric multivariate root finder using vectors/matrices defined by the Eigen package. 
 * It is a wrapper for the C gsl_multiroot functions from the GNU Scientific Library (https://www.gnu.org/software/gsl/doc/html/multiroots.html)
 * 
 * @author James Murray
 * Contact: jmurray@uwlax.edu
*/

#ifndef _INCL_MULTIROOT
#define _INCL_MULTIROOT

#include <Eigen/Dense>
#include <gsl/gsl_multiroots.h>
#include <functional>

/**
 * @brief Parameters passed to the GSL-wrapped root-finding function.
 */
struct multiroot_function_data {
    std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void*)> func;    ///< Function \f$ f : \mathbb{R}^n \to \mathbb{R}^n \f$ whose zeros are sought
    std::function<Eigen::MatrixXd(const Eigen::VectorXd&, const void*)> jacfunc; ///< Jacobian of func (optional; unused if not provided)
    const void* params; ///< User-supplied parameter pointer passed through to func()
};

/**
 * @brief Construct a multiroot_function_data struct with a function, Jacobian, and parameters.
 *
 * @param func     Function whose zeros are sought
 * @param jacfunc  Jacobian of func
 * @param params   User-supplied parameter pointer
 * @return Populated multiroot_function_data
 */
multiroot_function_data setup_multiroot_function_data(
    std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
    std::function<Eigen::MatrixXd(const Eigen::VectorXd&, const void* data)> jacfunc,
    const void* params);

/**
 * @brief Construct a multiroot_function_data struct with a function and parameters (no Jacobian).
 *
 * @param func   Function whose zeros are sought
 * @param params User-supplied parameter pointer
 * @return Populated multiroot_function_data
 */
multiroot_function_data setup_multiroot_function_data(
    std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
    const void* params);

/**
 * @brief GSL-compatible wrapper that evaluates the user-defined function.
 *
 * @param gsl_x  Current point as a gsl_vector
 * @param data   Pointer to a multiroot_function_data struct
 * @param gsl_f  Output: function value at gsl_x
 * @return GSL_SUCCESS on success
 */
int multiroot_gsl_f(const gsl_vector* gsl_x, void* data, gsl_vector* gsl_f);

/**
 * @brief GSL-compatible wrapper for the Jacobian (currently a no-op placeholder).
 *
 * @param gsl_x Current point as a gsl_vector
 * @param data  Pointer to a multiroot_function_data struct
 * @param J     Output: Jacobian matrix
 * @return GSL_SUCCESS on success
 */
int multiroot_gsl_df(const gsl_vector* gsl_x, void* data, gsl_matrix* J);

/**
 * @brief GSL-compatible wrapper that evaluates both the function and its Jacobian.
 *
 * @param x    Current point as a gsl_vector
 * @param data Pointer to a multiroot_function_data struct
 * @param f    Output: function value at x
 * @param J    Output: Jacobian matrix at x
 * @return GSL_SUCCESS on success
 */
int multiroot_gsl_fdf(const gsl_vector* x, void* data, gsl_vector* f, gsl_matrix* J);

/**
 * @brief Find zeros of a multivariate function using GSL (verbose defaults to false).
 *
 * @param initial_guess Starting point for the search
 * @param func          Function \f$ f(x, \text{params}) \f$ whose zeros are sought
 * @param params        User-supplied parameter pointer
 * @return Vector at which f is zero
 */
Eigen::VectorXd dml_multiroot(const Eigen::VectorXd& initial_guess, 
                              std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
                              const void* params);

/**
 * @brief Find zeros of a multivariate function with optional verbose output.
 *
 * @param initial_guess Starting point for the search
 * @param func          Function \f$ f(x, \text{params}) \f$ whose zeros are sought
 * @param params        User-supplied parameter pointer
 * @param verbose       If true, print iteration progress to stdout
 * @return Vector at which f is zero
 */
Eigen::VectorXd dml_multiroot(const Eigen::VectorXd& initial_guess, 
                              std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
                              const void* params, bool verbose);

/**
 * @brief Find zeros of a multivariate function using a user-supplied Jacobian.
 *
 * @param initial_guess Starting point for the search
 * @param func          Function \f$ f(x, \text{params}) \f$ whose zeros are sought
 * @param jacfunc       Jacobian of func
 * @param params        User-supplied parameter pointer
 * @param verbose       If true, print iteration progress to stdout
 * @return Vector at which f is zero
 */
Eigen::VectorXd dml_multiroot(const Eigen::VectorXd& initial_guess, 
        std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
        std::function<Eigen::MatrixXd(const Eigen::VectorXd&, const void* data)> jacfunc,
        const void* params, bool verbose);

#endif
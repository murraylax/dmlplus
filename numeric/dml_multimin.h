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
 * @file dml_multimin.h
 *
 * @brief 
 * This is a numeric multivariate minimizer using vectors/matrices defined by the Eigen package. 
 * It is a wrapper for the C gsl_multimin functions from the GNU Scientific Library (https://www.gnu.org/software/gsl/doc/html/multimin.html)
 * 
 * @author James Murray
 * Contact: jmurray@uwlax.edu
*/

#ifndef _INCL_MULTIMIN
#define _INCL_MULTIMIN

#include <Eigen/Dense>
#include <gsl/gsl_multimin.h>
#include <functional>

/**
 * @brief Parameters passed to the GSL-wrapped objective function.
 */
struct multimin_function_data {
    std::function<double(const Eigen::VectorXd&, const void*)> func; ///< Objective function \f$ f : \mathbb{R}^n \to \mathbb{R} \f$ to minimize
    const void* params;          ///< User-supplied parameter pointer passed through to func()
    Eigen::VectorXd upper_bounds; ///< Element-wise upper bounds on the solution
    Eigen::VectorXd lower_bounds; ///< Element-wise lower bounds on the solution
};

/**
 * @brief Construct a multimin_function_data struct.
 *
 * @param func         Objective function to minimize
 * @param params       User-supplied parameter pointer
 * @param lower_bounds Element-wise lower bounds
 * @param upper_bounds Element-wise upper bounds
 * @return Populated multimin_function_data
 */
multimin_function_data setup_multimin_function_data(std::function<double(const Eigen::VectorXd&, const void* data)> func,
                                 const void* params,
                                 Eigen::VectorXd& lower_bounds,
                                 Eigen::VectorXd& upper_bounds);

/**
 * @brief GSL-compatible wrapper that calls the user-defined objective function.
 *
 * @param gsl_x Current point as a gsl_vector
 * @param data  Pointer to a multimin_function_data struct
 * @return Objective function value at gsl_x
 */
double multimin_gsl_f(const gsl_vector* gsl_x, void* data);

/**
 * @brief Find the minimum of a multivariate function using GSL.
 *
 * @param initial_guess Starting point for the search
 * @param func          Objective function \f$ f(x, \text{params}) \f$
 * @param params        User-supplied parameter pointer
 * @return Vector at the minimum
 */
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params);

/**
 * @brief Find the minimum of a multivariate function with box constraints.
 *
 * @param initial_guess Starting point for the search
 * @param lower_bounds  Element-wise lower bounds
 * @param upper_bounds  Element-wise upper bounds
 * @param func          Objective function \f$ f(x, \text{params}) \f$
 * @param params        User-supplied parameter pointer
 * @return Vector at the minimum
 */
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              const Eigen::VectorXd& lower_bounds,
                              const Eigen::VectorXd& upper_bounds,
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params);

/**
 * @brief Find the minimum of a multivariate function with optional verbose output.
 *
 * @param initial_guess Starting point for the search
 * @param func          Objective function \f$ f(x, \text{params}) \f$
 * @param params        User-supplied parameter pointer
 * @param verbose       If true, print iteration progress to stdout
 * @return Vector at the minimum
 */
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params, bool verbose);

/**
 * @brief Find the minimum of a multivariate function with box constraints and optional verbose output.
 *
 * @param initial_guess Starting point for the search
 * @param lower_bounds  Element-wise lower bounds
 * @param upper_bounds  Element-wise upper bounds
 * @param func          Objective function \f$ f(x, \text{params}) \f$
 * @param params        User-supplied parameter pointer
 * @param verbose       If true, print iteration progress to stdout
 * @return Vector at the minimum
 */
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              const Eigen::VectorXd& lower_bounds,
                              const Eigen::VectorXd& upper_bounds,
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params, 
                              bool verbose);

#endif
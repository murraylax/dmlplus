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

// The parameter to pass to the GSL wrapped function. It includes both the user-provided function and the user-provided parameters
struct FunctionData {
    // func is the function to be solved. This is a function from R^n to R^n for which to find the zeros
    std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void*)> func;

    // This a generic object of parameters to pass to func()
    const void* params; 
};

// GSL wrapper for the user-defined function to minimize
int multiroot_gsl_f(const gsl_vector* gsl_x, void* data, gsl_vector* gsl_f);
// GSL wrapper for the Jacobian of the function. It does not currently do anything
int multiroot_gsl_df(const gsl_vector* gsl_x, void* data, gsl_matrix* J);
// Combined GSL wrapper for the user-defined function and the Jacobian
int multiroot_gsl_fdf(const gsl_vector* x, void* data, gsl_vector* f, gsl_matrix* J);

// Here's the big function
Eigen::VectorXd dml_multiroot(const Eigen::VectorXd& initial_guess, 
                              std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
                              const void* params);

#endif
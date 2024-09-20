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
 * It is a wrapper for the C gsl_multimin functions from the GNU Scientific Library (https://www.gnu.org/software/gsl/doc/html/multiroots.html)
 * 
 * @author James Murray
 * Contact: jmurray@uwlax.edu
*/

#ifndef _INCL_MULTIMIN
#define _INCL_MULTIMIN

#include <Eigen/Dense>
#include <gsl/gsl_multimin.h>
#include <functional>

// The parameter to pass to the GSL wrapped function. It includes both the user-provided function and the user-provided parameters
struct multimin_function_data {
    // func is the function to be solved. This is a function from R^n to R^1 to minimize
    std::function<double(const Eigen::VectorXd&, const void*)> func;

    // This a generic object of parameters to pass to func()
    const void* params; 

    // Upper bounds 
    Eigen::VectorXd upper_bounds;
    // Lower bounds 
    Eigen::VectorXd lower_bounds;

};

// A function to set up a FunctionData struct 
multimin_function_data setup_multimin_function_data(std::function<double(const Eigen::VectorXd&, const void* data)> func,
                                 const void* params,
                                 Eigen::VectorXd& lower_bounds,
                                 Eigen::VectorXd& upper_bounds);

// GSL wrapper for the user-defined function to minimize
double multimin_gsl_f(const gsl_vector* gsl_x, void* data);

// Here are all the versions of the dml_multimin() function

// Multiple dimensional minimizer with no lower and upper bounds, and no value for verbose
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params);


// Multiple dimensional minimizer with upper and lower bounds given, but no value for verbose
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              const Eigen::VectorXd& lower_bounds,
                              const Eigen::VectorXd& upper_bounds,
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params);

// Multiple dimensional minimizer with no lower and upper bounds given
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params, bool verbose);

// Multiple dimensional minimizer with lower and upper bounds
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              const Eigen::VectorXd& lower_bounds,
                              const Eigen::VectorXd& upper_bounds,
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params, 
                              bool verbose);

#endif
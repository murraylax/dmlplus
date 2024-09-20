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
 * @file dml_multimin.cpp
 *
 * @brief 
 * This is a numeric multivariate minimizer using vectors/matrices defined by the Eigen package. 
 * It is a wrapper for the C gsl_multimin functions from the GNU Scientific Library (https://www.gnu.org/software/gsl/doc/html/multiroots.html)
 * 
 * @author James Murray
 * Contact: jmurray@uwlax.edu
*/


#include <iostream>
#include <dml_multimin.h>
#include <utils.h>

using namespace std;

multimin_function_data setup_multimin_function_data(std::function<double(const Eigen::VectorXd&, const void* data)> func,
                                 const void* params,
                                 const Eigen::VectorXd& lower_bounds,
                                 const Eigen::VectorXd& upper_bounds) {

    multimin_function_data funcData;
    funcData.func = func;
    funcData.params = params;
    funcData.lower_bounds = lower_bounds;
    funcData.upper_bounds = upper_bounds;

    return funcData;
}

double multimin_gsl_f(const gsl_vector* gsl_x, void* data) {
    int n = gsl_x->size;
    multimin_function_data* funcData = static_cast<multimin_function_data*>(data);
    Eigen::VectorXd current_x(n);
    copy_gsl_to_vector(current_x, gsl_x);

    double fpenalty_multiplier = 1e+7;

    // Penalty if beyond the bounds
    double fpenalty = 0.0;

    if(funcData->lower_bounds.size() > 0 && funcData->upper_bounds.size() > 0) {
        for(int i = 0; i < n; i++) {
            if(current_x(i) > funcData->upper_bounds(i)) {
                fpenalty += fpenalty_multiplier * pow(current_x(i) - funcData->upper_bounds(i), 2.0);
                current_x(i) = funcData->upper_bounds(i); // Set the current_x to upper bound
            } 
            if(current_x(i) < funcData->lower_bounds(i)) {
                fpenalty += fpenalty_multiplier * pow(current_x(i) - funcData->lower_bounds(i), 2.0);
                current_x(i) = funcData->lower_bounds(i); // Set the current_x to lower bound
            }     
        }
    }

    // Call the user's function 
    double result = funcData->func(current_x, funcData->params);

    // Impose penalty if necessary
    result += fpenalty;

    return result;
}

// Multiroot finder with no lower and upper bounds, and no value for verbose
Eigen::VectorXd dml_multiroot(const Eigen::VectorXd& initial_guess, 
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params) {

    Eigen::VectorXd lower_bounds(0);
    Eigen::VectorXd upper_bounds(0);
    bool verbose = false;

    return dml_multimin(initial_guess, lower_bounds, upper_bounds, func, params, verbose);
}

// Multiroot finder with upper and lower bounds given, but no value for verbose
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              const Eigen::VectorXd& lower_bounds,
                              const Eigen::VectorXd& upper_bounds,
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params) {
    return dml_multimin(initial_guess, lower_bounds, upper_bounds, func, params, false);
}




// Multiroot finder with no lower and upper bounds given
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params, bool verbose) {

    Eigen::VectorXd lower_bounds(0);
    Eigen::VectorXd upper_bounds(0);

    return dml_multimin(initial_guess, lower_bounds, upper_bounds, func, params, verbose);
}
                             
// Multi-dimensional minimizer with lower and upper bounds
Eigen::VectorXd dml_multimin(const Eigen::VectorXd& initial_guess, 
                              const Eigen::VectorXd& lower_bounds,
                              const Eigen::VectorXd& upper_bounds,
                              std::function<double(const Eigen::VectorXd&, const void* data)> func,
                              const void* params, 
                              bool verbose) {

    int max_iterations = 100000; // Maximum number of iterations
    double ftol = 1e-8; // 

    const size_t dimension = initial_guess.size();

    // Set up the function data
    multimin_function_data funcData = setup_multimin_function_data(func, params, lower_bounds, upper_bounds);

    // Set up GSL minimizer
    const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(T, dimension);

    // Set up GSL function
    gsl_multimin_function system;
    system.n = dimension;           // Number of variables
    system.f = multimin_gsl_f;               // Objective function
    system.params = &funcData;      // Pass FunctionData as params

    // Initial Guess
    gsl_vector* gsl_x = gsl_vector_alloc(dimension);
    copy_vector_to_gsl(gsl_x, initial_guess);

    // Initial step size
    gsl_vector* step_size = gsl_vector_alloc(dimension);
    gsl_vector_set_all(step_size, 0.1);  

    // Set up solver
    gsl_multimin_fminimizer_set(minimizer, &system, gsl_x, step_size);

    // This will be the return solution
    Eigen::VectorXd result(dimension);

    int status = GSL_CONTINUE;
    int iter = 0;
    double size;

    while(status==GSL_CONTINUE && iter < max_iterations) {
        iter++;

        if(verbose) {
            copy_gsl_to_vector(result, minimizer->x);
            cout << "Iteration: " << iter << "\nCurent value: " << result.transpose() << endl;
        }
        
        status = gsl_multimin_fminimizer_iterate(minimizer);
        if(status) {
            throw std::runtime_error("Minimizer iteration failed: " + std::string(gsl_strerror(status)));
        }

        // Test for convergence
        size = gsl_multimin_fminimizer_size(minimizer);
        status = gsl_multimin_test_size(size, ftol);
    }

    if(status != GSL_SUCCESS) {
        throw std::runtime_error("Minimizer did not converge: " + std::string(gsl_strerror(status)));
    }
    
    // Convert back to Eigen::VectorXd
    copy_gsl_to_vector(result, minimizer->x);

    // Clean up
    gsl_multimin_fminimizer_free(minimizer);
    gsl_vector_free(step_size);
    gsl_vector_free(gsl_x);

    return result;
}
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
 * @file dml_multiroot.cpp
 *
 * @brief 
 * This is a numeric multivariate root finder using vectors/matrices defined by the Eigen package. 
 * It is a wrapper for the C gsl_multiroot functions from the GNU Scientific Library (https://www.gnu.org/software/gsl/doc/html/multiroots.html)
 * 
 * @author James Murray
 * Contact: jmurray@uwlax.edu
*/

#include <iostream>
#include <utils.h>
#include "dml_multiroot.h"

using namespace std;

int multiroot_gsl_f(const gsl_vector* gsl_x, void* data, gsl_vector* gsl_f) {
    FunctionData* funcData = static_cast<FunctionData*>(data);
    Eigen::VectorXd x(gsl_x->size);
    copy_gsl_to_vector(x, gsl_x);

    // Call the user's function 
    Eigen::VectorXd result = funcData->func(x, funcData->params);
    // Copy the output to a gsl_vector for the output
    copy_vector_to_gsl(gsl_f, result);

    return GSL_SUCCESS;
}

// The Jacobian of the function, wrapped for GSL
int multiroot_gsl_df(const gsl_vector* gsl_x, void* data, gsl_matrix* J) {
    // Optional: You can implement the Jacobian here if available, or use GSL's numerical derivative solver.
    // Returning GSL_ENOMEM to signal we don't provide an analytical Jacobian here.
    return GSL_ENOMEM;
}

// The combined function and Jacobian for GSL
int multiroot_gsl_fdf(const gsl_vector* x, void* data, gsl_vector* f, gsl_matrix* J) {
    multiroot_gsl_f(x, data, f);
    multiroot_gsl_df(x, data, J);
    return GSL_SUCCESS;
}

Eigen::VectorXd dml_multiroot(const Eigen::VectorXd& initial_guess, 
                              std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
                              const void* params) {

    int max_iterations = 100000; // Maximum number of iterations
    double ftol = 1e-8; // 

    const size_t dimension = initial_guess.size();

    // Set up the function data
    FunctionData funcData;
    funcData.func = func;
    funcData.params = params;

    // Set up GSL multiroot solver
    const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(T, dimension);

    // GSL Multiroot System
    gsl_multiroot_function system = {&multiroot_gsl_f, dimension, &funcData};

    // Initial Guess
    gsl_vector* gsl_x = gsl_vector_alloc(dimension);
    copy_vector_to_gsl(gsl_x, initial_guess);

    // Set up solver
    gsl_multiroot_fsolver_set(solver, &system, gsl_x);

    Eigen::VectorXd result(dimension);

    int status = GSL_CONTINUE;
    int iter = 0;
    while(status==GSL_CONTINUE && iter < max_iterations) {
        iter++;

        status = gsl_multiroot_fsolver_iterate(solver);
        if(status) {
            throw std::runtime_error("Howdy! Solver failed: " + std::string(gsl_strerror(status)));
        }

        status = gsl_multiroot_test_residual(solver->f, ftol);
    }

    if(status != GSL_SUCCESS) {
        throw std::runtime_error("Solver did not converge: " + std::string(gsl_strerror(status)));
    }
    
    // Convert back to Eigen::VectorXd
    copy_gsl_to_vector(result, solver->x);

    // Clean up
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(gsl_x);

    return result;
}
                             
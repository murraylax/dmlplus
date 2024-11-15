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
#include <dml_multiroot.h>
#include <utils.h>

using namespace std;
 
multiroot_function_data setup_multiroot_function_data(
        std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
        std::function<Eigen::MatrixXd(const Eigen::VectorXd&, const void* data)> jacfunc,
        const void* params,
        const Eigen::VectorXd& lower_bounds,
        const Eigen::VectorXd& upper_bounds) {

    multiroot_function_data funcData;
    funcData.func = func;
    funcData.jacfunc = jacfunc;
    funcData.params = params;
    funcData.lower_bounds = lower_bounds;
    funcData.upper_bounds = upper_bounds;

    return funcData;
}

multiroot_function_data setup_multiroot_function_data(
        std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
        const void* params,
        const Eigen::VectorXd& lower_bounds,
        const Eigen::VectorXd& upper_bounds) {

    multiroot_function_data funcData;
    funcData.func = func;
    funcData.params = params;
    funcData.lower_bounds = lower_bounds;
    funcData.upper_bounds = upper_bounds;

    return funcData;
}

int multiroot_gsl_f(const gsl_vector* gsl_x, void* data, gsl_vector* gsl_f) {
    int n = gsl_x->size;
    multiroot_function_data* funcData = static_cast<multiroot_function_data*>(data);
    Eigen::VectorXd current_x(n);
    copy_gsl_to_vector(current_x, gsl_x);

    double fpenalty_multiplier = 1000.0;

    // Penalty if beyond the bounds
    Eigen::VectorXd penalty(n);
    penalty.setConstant(0.0);

    if(funcData->lower_bounds.size() > 0 && funcData->upper_bounds.size() > 0) {
        for(int i = 0; i < n; i++) {
            if(current_x(i) > funcData->upper_bounds(i)) {
                penalty(i) = fpenalty_multiplier * pow(current_x(i) - funcData->upper_bounds(i), 2.0);
                current_x(i) = funcData->upper_bounds(i); // Set the current_x to upper bound
            } 
            if(current_x(i) < funcData->lower_bounds(i)) {
                penalty(i) = fpenalty_multiplier * pow(current_x(i) - funcData->lower_bounds(i), 2.0);
                current_x(i) = funcData->lower_bounds(i); // Set the current_x to lower bound
            }     
        }
    }

    // Call the user's function 
    Eigen::VectorXd result = funcData->func(current_x, funcData->params);

    // Impose penalty if necessary
    for(int i = 0; i < n; i++) {
        if(penalty(i) > 0) {
            result(i) = result(i)*result(i) + penalty(i);
        }
    }

    // Copy the result to a gsl_vector for the output
    copy_vector_to_gsl(gsl_f, result);

    return GSL_SUCCESS;
}

// The Jacobian of the function, wrapped for GSL
int multiroot_gsl_df(const gsl_vector* gsl_x, void* data, gsl_matrix* gsl_J) {
    int n = gsl_x->size;
    multiroot_function_data* funcData = static_cast<multiroot_function_data*>(data);
 
    Eigen::VectorXd current_x(n);
    copy_gsl_to_vector(current_x, gsl_x);

    Eigen::MatrixXd jac = funcData->jacfunc(current_x, funcData->params);
    copy_matrix_to_gsl(gsl_J, jac);

    return GSL_SUCCESS;
}

// The combined function and Jacobian for GSL
int multiroot_gsl_fdf(const gsl_vector* x, void* data, gsl_vector* f, gsl_matrix* J) {
    multiroot_gsl_f(x, data, f);
    multiroot_gsl_df(x, data, J);
    return GSL_SUCCESS;
}

// Multiroot finder with no lower and upper bounds, and no value for verbose
Eigen::VectorXd dml_multiroot(const Eigen::VectorXd& initial_guess, 
                              std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
                              const void* params) {

    Eigen::VectorXd lower_bounds(0);
    Eigen::VectorXd upper_bounds(0);
    bool verbose = false;

    return dml_multiroot(initial_guess, lower_bounds, upper_bounds, func, params, verbose);
}

// Multiroot finder with upper and lower bounds given, but no value for verbose
Eigen::VectorXd dml_multiroot(const Eigen::VectorXd& initial_guess, 
                              const Eigen::VectorXd& lower_bounds,
                              const Eigen::VectorXd& upper_bounds,
                              std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
                              const void* params) {
    return dml_multiroot(initial_guess, lower_bounds, upper_bounds, func, params, false);
}




// Multiroot finder with no lower and upper bounds given
Eigen::VectorXd dml_multiroot(const Eigen::VectorXd& initial_guess, 
                              std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
                              const void* params, bool verbose) {

    Eigen::VectorXd lower_bounds(0);
    Eigen::VectorXd upper_bounds(0);

    return dml_multiroot(initial_guess, lower_bounds, upper_bounds, func, params, verbose);
}
        
// Multiroot finder with lower and upper bounds - No Jacobian
Eigen::VectorXd dml_multiroot(
        const Eigen::VectorXd& initial_guess, 
        const Eigen::VectorXd& lower_bounds,
        const Eigen::VectorXd& upper_bounds,
        std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
        const void* params, 
        bool verbose) {

    int max_iterations = 100000; // Maximum number of iterations
    double ftol = 1e-8; // 

    const size_t dimension = initial_guess.size();

    // Initial Guess
    gsl_vector* gsl_x = gsl_vector_alloc(dimension);
    copy_vector_to_gsl(gsl_x, initial_guess);

    // Set up the function data
    multiroot_function_data funcData = setup_multiroot_function_data(func, params, lower_bounds, upper_bounds);

    // Set up GSL multiroot solver
    const gsl_multiroot_fsolver_type* ftype = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver* fsolver;
    gsl_multiroot_function fsystem; 

    fsystem.f = &multiroot_gsl_f;
    fsystem.n = dimension;
    fsystem.params = &funcData;
    fsolver = gsl_multiroot_fsolver_alloc(ftype, dimension);
    gsl_multiroot_fsolver_set(fsolver, &fsystem, gsl_x);

    // This will be the return solution
    Eigen::VectorXd result(dimension);

    int status = GSL_CONTINUE;
    int iter = 0;
    while(status==GSL_CONTINUE && iter < max_iterations) {
        iter++;
        
        if(verbose) {
            copy_gsl_to_vector(result, fsolver->x);
            cout << "Iteration: " << iter << "\nCurent value: " << result.transpose() << endl;
        }

        // Iterate the solver
        status = gsl_multiroot_fsolver_iterate(fsolver);

        if(status) {
            throw std::runtime_error("Solver failed: " + std::string(gsl_strerror(status)));
        }

        // Check residual
        status = gsl_multiroot_test_residual(fsolver->f, ftol);
    }

    if(status != GSL_SUCCESS) {
        throw std::runtime_error("Solver did not converge: " + std::string(gsl_strerror(status)));
    }
    
    // Convert back to Eigen::VectorXd
    copy_gsl_to_vector(result, fsolver->x);

    // Clean up
    gsl_multiroot_fsolver_free(fsolver);
    gsl_vector_free(gsl_x);

    return result;
}

// Multiroot finder with Jacobian, with no lower and upper bounds given
Eigen::VectorXd dml_multiroot(const Eigen::VectorXd& initial_guess, 
        std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
        std::function<Eigen::MatrixXd(const Eigen::VectorXd&, const void* data)> jacfunc,
        const void* params, bool verbose) {

    Eigen::VectorXd lower_bounds(0);
    Eigen::VectorXd upper_bounds(0);

    return dml_multiroot(initial_guess, lower_bounds, upper_bounds, func, jacfunc, params, verbose);
}

// Multiroot finder with lower and upper bounds and Jacobian given
Eigen::VectorXd dml_multiroot(
        const Eigen::VectorXd& initial_guess, 
        const Eigen::VectorXd& lower_bounds,
        const Eigen::VectorXd& upper_bounds,
        std::function<Eigen::VectorXd(const Eigen::VectorXd&, const void* data)> func,
        std::function<Eigen::MatrixXd(const Eigen::VectorXd&, const void* data)> jacfunc,
        const void* params, 
        bool verbose) {

    int max_iterations = 100000; // Maximum number of iterations
    double ftol = 1e-8; // 

    const size_t dimension = initial_guess.size();

    bool jacobian_exists = static_cast<bool>(jacfunc);

    // Initial Guess
    gsl_vector* gsl_x = gsl_vector_alloc(dimension);
    copy_vector_to_gsl(gsl_x, initial_guess);

    // Set up the function data
    multiroot_function_data funcData = setup_multiroot_function_data(func, jacfunc, params, lower_bounds, upper_bounds);

    // Set up GSL multiroot solver
    const gsl_multiroot_fdfsolver_type* fdftype = gsl_multiroot_fdfsolver_hybridsj;
    gsl_multiroot_fdfsolver* fdfsolver;
    gsl_multiroot_function_fdf fdfsystem; 

    fdfsystem.f = &multiroot_gsl_f;
    fdfsystem.df = &multiroot_gsl_df;
    fdfsystem.fdf = &multiroot_gsl_fdf;
    fdfsystem.n = dimension;
    fdfsystem.params = &funcData;
    fdfsolver = gsl_multiroot_fdfsolver_alloc(fdftype, dimension);
    gsl_multiroot_fdfsolver_set(fdfsolver, &fdfsystem, gsl_x);

    // This will be the return solution
    Eigen::VectorXd result(dimension);

    int status = GSL_CONTINUE;
    int iter = 0;
    while(status==GSL_CONTINUE && iter < max_iterations) {
        iter++;
        
        if(verbose) {
            copy_gsl_to_vector(result, fdfsolver->x);
            cout << "Iteration: " << iter << "\nCurent value: " << result.transpose() << endl;
        }

        // Iterate the solver
        status = gsl_multiroot_fdfsolver_iterate(fdfsolver);

        if(status) {
            throw std::runtime_error("Solver failed: " + std::string(gsl_strerror(status)));
        }

        // Check residual
        status = gsl_multiroot_test_residual(fdfsolver->f, ftol);
    }

    if(status != GSL_SUCCESS) {
        throw std::runtime_error("Solver did not converge: " + std::string(gsl_strerror(status)));
    }
    
    // Convert back to Eigen::VectorXd
    copy_gsl_to_vector(result, fdfsolver->x);

    // Clean up
    gsl_multiroot_fdfsolver_free(fdfsolver);
    gsl_vector_free(gsl_x);

    return result;
}
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
 * @file gensys.cpp
 *
 * @brief
 * Use Sims (2002) method for solving a linear dynamic general equilibrium model
 *
 * Model has the following form: Gamma_0 x_t = C + Gamma_1 x_t-1 + Psi z_t + Pi eta_t
 *
 * where...
 *  x_t is a vector of variables
 *  z_t is a vector of exogenous shocks
 *  eta_t is a vector of endogenous expectation errors: eta_t = x_t - E_t-1 x_t
 *
 *  C is a vector of constants
 *  Gamma_0 and Gamma_1 are square matrices of coefficients, dimensions nvar x nvar
 *  Psi is a matrix, dimension nvar x nshocks
 *  Pi is a matrix, dimension nvar x nendo
 *
 * The solution is given as,
 *
 * x_t = D_sol + G_sol x_t-1 + M_sol z_t
 *
 * The method includes a return value for existence and uniqueness of solution:
 *   Return value > 0: Indeterminacy, there are `return_value` number of loose endogenous variables, ignore outputs for D_sol, G_sol, and M_sol
 *   Return value = 0: Unique solution is found. There are 0 loose endogenous variables.
 *   Return value = -1: No solution exists.
 *
 * @author James Murray
 * Contact: james@murraylax.org
*/

// [[Rcpp::import(Rcpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include "gensys.h"
#include <RcppEigen.h>

/**
 * checksysR(const Eigen::Map<Eigen::MatrixXd>& Gamma0, const Eigen::Map<Eigen::MatrixXd>& Gamma1, Eigen::Map<Eigen::MatrixXd>& Psi,
              const Eigen::Map<Eigen::MatrixXd>& Pi, Eigen::Map<Eigen::VectorXd>& C)
 *
 * Use Sims (2002) method to check for a unique solution for a linear dynamic general equilibrium model
 *
 * Model has the following form: Gamma_0 x_t = C + Gamma_1 x_t-1 + Psi z_t + Pi eta_t
 *
 * where...
 *  x_t is a vector of variables
 *  z_t is a vector of exogenous shocks
 *  eta_t is a vector of endogenous expectation errors: eta_t = x_t - E_t-1 x_t
 *
 *  C is a vector of constants
 *  Gamma_0 and Gamma_1 are square matrices of coefficients, dimensions nvar x nvar
 *  Psi is a matrix, dimension nvar x nshocks
 *  Pi is a matrix, dimension nvar x nendo
 *
 * The solution is given as,
 *
 * x_t = D_sol + G_sol x_t-1 + M_sol z_t
 *
 * The method includes a return value for existence and uniqueness of solution:
 *   Return value > 0: Indeterminacy, there are `return_value` number of loose endogenous variables, ignore outputs for D_sol, G_sol, and M_sol
 *   Return value = 0: Unique solution is found. There are 0 loose endogenous variables.
 *   Return value = -1: No solution exists.
 *
 * @param Gamma0 (input) The coefficient matrix Gamma_0
 * @param Gamma1 (input) The coefficient matrix Gamma_1
 * @param Psi (input) The coefficient matrix Psi
 * @param Pi (input) The coefficient matrix Pi
 * @param C (input) The constant vector C
 *
 * @return Integer equal to the number of loose parameters. =0 when there is a unique solution, =-1 for no solution, >0 for indeterminacy
*/
// [[Rcpp::export]]
int checksysR(const Eigen::Map<Eigen::MatrixXd>& Gamma0, const Eigen::Map<Eigen::MatrixXd>& Gamma1, const Eigen::Map<Eigen::MatrixXd>& Psi,
              const Eigen::Map<Eigen::MatrixXd>& Pi, const Eigen::Map<Eigen::VectorXd>& C) {

    return checksys(Gamma0, Gamma1, Psi, Pi, C);
}



/**
 * Rcpp::List gensysR(const Eigen::Map<Eigen::MatrixXd>& Gamma0, const Eigen::Map<Eigen::MatrixXd>& Gamma1,
                    const Eigen::Map<Eigen::MatrixXd>& Psi, const Eigen::Map<Eigen::MatrixXd>& Pi, const Eigen::Map<Eigen::VectorXd>& C)
 *
 * Use Sims (2002) method for solving a linear dynamic general equilibrium model
 *
 * Model has the following form: Gamma_0 x_t = C + Gamma_1 x_t-1 + Psi z_t + Pi eta_{t}
 *
 * where...
 *  x_t is a vector of variables
 *  z_t is a vector of exogenous shocks
 *  eta_t is a vector of endogenous expectation errors: eta_t = x_t - E_t-1 x_t
 *
 *  C is a vector of constants
 *  Gamma_0 and Gamma_1 are square matrices of coefficients, dimensions nvar x nvar
 *  Psi is a matrix, dimension nvar x nshocks
 *  Pi is a matrix, dimension nvar x nendo
 *
 * The solution is given as,
 *
 * x_t = D_sol + G_sol x_t-1 + M_sol z_t
 *
 * The method includes a return value for existence and uniqueness of solution:
 *   Return value > 0: Indeterminacy, there are `return_value` number of loose endogenous variables, ignore outputs for D_sol, G_sol, and M_sol
 *   Return value = 0: Unique solution is found. There are 0 loose endogenous variables.
 *   Return value = -1: No solution exists.
 *
 * @param Gamma0 The coefficient matrix Gamma_0
 * @param Gamma1 The coefficient matrix Gamma_1
 * @param Psi The coefficient matrix Psi
 * @param Pi The coefficient matrix Pi
 * @param C The constant vector C
 *
 * @return A list with the solution and an integer specifying whether the solution exists and is unique.
 *   The elements of the list are as follows:
 *      Gsol The solution matrix G_sol
 *      Msol The solution matrix M_sol
 *      Dsol The solution vector D_sol
 *      nloose The number of loose parameters. nloose=0 when there is a unique solution, nloose=-1 for no solution, nloose>0 for indeterminacy
*/
// [[Rcpp::export]]
Rcpp::List gensysR(const Eigen::Map<Eigen::MatrixXd>& Gamma0, const Eigen::Map<Eigen::MatrixXd>& Gamma1,
                    const Eigen::Map<Eigen::MatrixXd>& Psi, const Eigen::Map<Eigen::MatrixXd>& Pi, const Eigen::Map<Eigen::VectorXd>& C) {

    int n = Gamma0.rows();
    int k = Psi.cols();

    Eigen::MatrixXd Gsol(n,n);
    Eigen::MatrixXd Msol(n,k);
    Eigen::VectorXd Dsol(n);

    int nloose = gensys(Gsol, Msol, Dsol, Gamma0, Gamma1, Psi, Pi, C);

    return Rcpp::List::create(
        Rcpp::Named("Gsol") = Gsol,
        Rcpp::Named("Msol") = Msol,
        Rcpp::Named("Dsol") = Dsol,
        Rcpp::Named("nloose") = nloose
    );
}


/**
 * Rcpp::List gensysR_qzdetails(const Eigen::Map<Eigen::MatrixXd>& Gamma0, const Eigen::Map<Eigen::MatrixXd>& Gamma1,
                    const Eigen::Map<Eigen::MatrixXd>& Psi, const Eigen::Map<Eigen::MatrixXd>& Pi, const Eigen::Map<Eigen::VectorXd>& C)
 *
 * Use Sims (2002) method for solving a linear dynamic general equilibrium model.
 * Include in the return value the full details of the QZ decomposition
 *
 * Model has the following form: Gamma_0 x_t = C + Gamma_1 x_t-1 + Psi z_t + Pi eta_{t}
 *
 * where...
 *  x_t is a vector of variables
 *  z_t is a vector of exogenous shocks
 *  eta_t is a vector of endogenous expectation errors: eta_t = x_t - E_t-1 x_t
 *
 *  C is a vector of constants
 *  Gamma_0 and Gamma_1 are square matrices of coefficients, dimensions nvar x nvar
 *  Psi is a matrix, dimension nvar x nshocks
 *  Pi is a matrix, dimension nvar x nendo
 *
 * The solution is given as,
 *
 * x_t = D_sol + G_sol x_t-1 + M_sol z_t
 *
 * The method includes a return value for existence and uniqueness of solution:
 *   Return value > 0: Indeterminacy, there are `return_value` number of loose endogenous variables, ignore outputs for D_sol, G_sol, and M_sol
 *   Return value = 0: Unique solution is found. There are 0 loose endogenous variables.
 *   Return value = -1: No solution exists.
 *
 * Solution involves the QZ decomposition of Gamma0, Gamma1, where
 *   Q' Gamma_0 Z = S
 *   Q' Gamma_1 Z = T
 *   Q'Q = Q Q' = Z'Z = Z Z' = I
 *   Eigenvalues, lambda_i = abs(T_ii / S_ii)
 *   And eigenvalues are ordered where all those at the bottom are >= 1.0 in magnitude
 *
 * @param Gamma0 The coefficient matrix Gamma_0
 * @param Gamma1 The coefficient matrix Gamma_1
 * @param Psi The coefficient matrix Psi
 * @param Pi The coefficient matrix Pi
 * @param C The constant vector C
 *
 * @return A list with the solution and an integer specifying whether the solution exists and is unique.
 *   The elements of the list are as follows:
 *      Gsol The solution matrix G_sol
 *      Msol The solution matrix M_sol
 *      Dsol The solution vector D_sol
 *      S The QZ decomposition matrix S above
 *      T The QZ decomposition matrix T above
 *      Q The QZ decomposition matrix Q above
 *      Z The QZ decomposition matrix Z above
 *      Lambda The real magnitude of the eigenvalues
 *      nunstable The number of unstable eigenvalues
 *      nloose The number of loose parameters. nloose=0 when there is a unique solution, nloose=-1 for no solution, nloose>0 for indeterminacy
*/
// [[Rcpp::export]]
Rcpp::List gensysR_qzdetails(const Eigen::Map<Eigen::MatrixXd>& Gamma0, const Eigen::Map<Eigen::MatrixXd>& Gamma1,
                    const Eigen::Map<Eigen::MatrixXd>& Psi, const Eigen::Map<Eigen::MatrixXd>& Pi, const Eigen::Map<Eigen::VectorXd>& C) {

    int n = Gamma0.rows();
    int k = Psi.cols();

    Eigen::MatrixXd Gsol(n,n);
    Eigen::MatrixXd Msol(n,k);
    Eigen::VectorXd Dsol(n);
    Eigen::MatrixXcd Q(n,n);
    Eigen::MatrixXcd Z(n,n);
    Eigen::MatrixXcd S(n,n);
    Eigen::MatrixXcd T(n,n);
    Eigen::VectorXd Lambda(n);
    int nunstable;

    int nloose = gensys_qzdetails(Gsol, Msol, Dsol,
        S, T, Q, Z, Lambda, nunstable,
        Gamma0, Gamma1, Psi, Pi, C);

    return Rcpp::List::create(
        Rcpp::Named("Gsol") = Gsol,
        Rcpp::Named("Msol") = Msol,
        Rcpp::Named("Dsol") = Dsol,
        Rcpp::Named("S") = S,
        Rcpp::Named("T") = T,
        Rcpp::Named("Q") = Q,
        Rcpp::Named("Z") = Z,
        Rcpp::Named("Lambda") = Lambda,
        Rcpp::Named("nunstable") = nunstable,
        Rcpp::Named("nloose") = nloose
    );
}

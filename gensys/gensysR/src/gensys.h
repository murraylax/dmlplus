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
 * @file gensys.h
 *
 * @brief 
 * Use Sims (2002) method for solving a linear dynamic general equilibrium model
 * 
 * Model has the following form: \Gamma_0 x_t = C + \Gamma_1 x_t-1 + \Psi z_t + \Pi \eta_t
 * 
 * where...
 *  x_t is a vector of variables
 *  z_t is a vector of exogenous shocks
 *  \eta_t is a vector of endogenous expectation errors: \eta_t = x_t - E_t-1 x_t
 * 
 *  C is a vector of constants
 *  \Gamma_0 and \Gamma_1 are square matrices of coefficients, dimensions nvar x nvar
 *  \Psi is a matrix, dimension nvar x nshocks
 *  \Pi is a matrix, dimension nvar x nendo
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

#ifndef _INCL_GENSYS
#define _INCL_GENSYS

#include <Eigen/Dense>

/**
 * gensys(Eigen::MatrixXd& mdGsol, Eigen::MatrixXd& mdMsol, Eigen::MatrixXd& Dsol, const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC)
 * 
 * Use Sims (2002) method for solving a linear dynamic general equilibrium model
 * 
 * Model has the following form: \Gamma_0 x_t = C + \Gamma_1 x_t-1 + \Psi z_t + \Pi \eta_t
 * 
 * where...
 *  x_t is a vector of variables
 *  z_t is a vector of exogenous shocks
 *  \eta_t is a vector of endogenous expectation errors: \eta_t = x_t - E_t-1 x_t
 * 
 *  C is a vector of constants
 *  \Gamma_0 and \Gamma_1 are square matrices of coefficients, dimensions nvar x nvar
 *  \Psi is a matrix, dimension nvar x nshocks
 *  \Pi is a matrix, dimension nvar x nendo
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
 * @param mdGsol (output) Eigen::MatrixXd with the solution matrix G_{sol}
 * @param mdMsol (output) Eigen::MatrixXd with the solution matrix M_{sol}
 * @param vcDsol (output) Eigen::VectorXd with the solution vector D_{sol}
 * 
 * @param mdGamma0 (input) const Eigen::MatrixXd with the coefficient matrix \Gamma_0
 * @param mdGamma1 (input) const Eigen::MatrixXd with the coefficient matrix \Gamma_1
 * @param mdPsi (input) const Eigen::MatrixXd with the coefficient matrix \Psi
 * @param mdPi (input) const Eigen::MatrixXd with the coefficient matrix \Pi
 * @param vdC (input) const Eigen::VectorXd with the constant vector C
 * 
*/
int gensys(Eigen::MatrixXd& mdGsol, Eigen::MatrixXd& mdMsol, Eigen::VectorXd& vDsol, const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC);

int gensys_qzdetails(Eigen::MatrixXd& mdGsol, Eigen::MatrixXd& mdMsol, Eigen::VectorXd& vdDsol, 
     Eigen::MatrixXcd& mcdS, Eigen::MatrixXcd& mcdT, Eigen::MatrixXcd& mcdQ, Eigen::MatrixXcd& mcdZ, 
     Eigen::VectorXd& vdLambda, int& nunstable,
     const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, 
     const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC);
/**
 * checksys(Eigen::MatrixXd& mdGsol, Eigen::MatrixXd& mdMsol, Eigen::MatrixXd& Dsol, const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC)
 * 
 * Use Sims (2002) method for solving a linear dynamic general equilibrium model to determine if there is a unique solution
 * 
 * Model has the following form: \Gamma_0 x_t = C + \Gamma_1 x_t-1 + \Psi z_t + \Pi \eta_t
 * 
 * where...
 *  x_t is a vector of variables
 *  z_t is a vector of exogenous shocks
 *  \eta_t is a vector of endogenous expectation errors: \eta_t = x_t - E_t-1 x_t
 * 
 *  C is a vector of constants
 *  \Gamma_0 and \Gamma_1 are square matrices of coefficients, dimensions nvar x nvar
 *  \Psi is a matrix, dimension nvar x nshocks
 *  \Pi is a matrix, dimension nvar x nendo
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
 * @param mdGamma0 (input) const Eigen::MatrixXd with the coefficient matrix \Gamma_0
 * @param mdGamma1 (input) const Eigen::MatrixXd with the coefficient matrix \Gamma_1
 * @param mdPsi (input) const Eigen::MatrixXd with the coefficient matrix \Psi
 * @param mdPi (input) const Eigen::MatrixXd with the coefficient matrix \Pi
 * @param vdC (input) const Eigen::VectorXd with the constant vector C
 * 
*/
int checksys(const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC);


#endif
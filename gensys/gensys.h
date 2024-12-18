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
 * write_irf(const Eigen::MatrixXd& mdIRF, const std::vector<std::string>& varnames, const std::string& filepath)
 * 
 * Write the impulse response functions to a text file to be easily read in R
 * 
 * @param mdIRF Matrix of impulse responses, return varlue of gensys_irf()
 * @param varnames Vector of strings for the variable names
 * @param filepath The filepath for the output
 */
void write_irf(const Eigen::MatrixXd& mdIRF, const std::vector<std::string>& varnames, const std::string& filepath);

/**
 * Eigen::MatrixXd gensys_irf(const Eigen::MatrixXd& mdG, const Eigen::MatrixXd& mdM, double fshock, size_t shock_idx, size_t nirf)
 * 
 * Compute impulse responses for a given shock for a given solution of the dynamic system,  x_t = D + G x_t-1 + M z_t
 * 
 * @param mdG: Matrix G in the solution
 * @param mdM: Matrix M in the solution
 * @param fshock: Magnitude of the shock at time t=0
 * @param shock_idx: Index into z_t for the specific shock
 * @param nirf: Number of periods for the impulse response
 * 
 * @return Matrix size (nirf x nvar) for the impulse responses, where each row t is the response of x_t 
 */
Eigen::MatrixXd gensys_irf(const Eigen::MatrixXd& mdG, const Eigen::MatrixXd& mdM, double fshock, size_t shock_idx, size_t nirf);

/**
 * Eigen::MatrixXd gensys_irf(const Eigen::MatrixXd& mdG, const Eigen::MatrixXd& mdM, size_t nirf)
 * 
 * Compute impulse responses for a all shocks for a given solution of the dynamic system,  x_t = D + G x_t-1 + M z_t
 * Assumes a magnitude for the shock = 0.01
 * 
 * @param mdG: Matrix G in the solution
 * @param mdM: Matrix M in the solution
 * @param nirf: Number of periods for the impulse response
 * 
 * @return Matrix size ((nirf*nshocks) x (nvar+2)) for the impulse responses, where each row t is the response of x_t. The second-to-last column is the index of the shock, last column is the time period.
 */
Eigen::MatrixXd gensys_irf(const Eigen::MatrixXd& mdG, const Eigen::MatrixXd& mdM, size_t nirf);

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

/**
 * Write impulse response functions to a csv file, designed for easy import into R
 * 
 * @param mdIRF: Eigen::MatrixXd of impulse response functions
 * @param varnames: Vector of strings for the variable names - columns of mdIRF
 * @param shocknames: Vector of strings for the shock names 
 * @param desc: String for a description of the IRF - will be entered in the last column
 * @param csvfile: file to write, std::ofstream
 */
void write_irf_to_csvfile(const Eigen::MatrixXd& mdIRF, std::vector<std::string>& varnames, std::vector<std::string>& shocknames, std::string& desc, std::ofstream& csvfile);

/**
 * Write impulse response functions to a csv file, designed for easy import into R
 * 
 * No header row
 * 
 * @param mdIRF: Eigen::MatrixXd of impulse response functions
 * @param shocknames: Vector of strings for the shock names 
 * @param desc: String for a description of the IRF - will be entered in the last column
 * @param csvfile: file to write, std::ofstream
 */
void write_irf_to_csvfile(const Eigen::MatrixXd& mdIRF, std::vector<std::string>& shocknames, std::string& desc, std::ofstream& csvfile);

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
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
 * @brief Write impulse response functions to a CSV file with a header row, designed for easy import into R.
 *
 * @param mdIRF Eigen::MatrixXd of impulse response functions
 * @param varnames Vector of strings for the variable names — columns of mdIRF
 * @param shocknames Vector of strings for the shock names
 * @param desc String description of the IRF — written into the last column
 * @param csvfile Open output file stream to write to
 */
void write_irf_to_csvfile(const Eigen::MatrixXd& mdIRF, std::vector<std::string>& varnames, std::vector<std::string>& shocknames, std::string& desc, std::ofstream& csvfile);

/**
 * @brief Write impulse response functions to a CSV file with no header row, designed for easy import into R.
 *
 * @param mdIRF Eigen::MatrixXd of impulse response functions
 * @param shocknames Vector of strings for the shock names
 * @param desc String description of the IRF — written into the last column
 * @param csvfile Open output file stream to write to
 */
void write_irf_to_csvfile(const Eigen::MatrixXd& mdIRF, std::vector<std::string>& shocknames, std::string& desc, std::ofstream& csvfile);

/**
 * @brief Solve a linear dynamic general equilibrium model using the Sims (2002) method.
 *
 * Model form: \f$ \Gamma_0 x_t = C + \Gamma_1 x_{t-1} + \Psi z_t + \Pi \eta_t \f$
 *
 * where \f$ x_t \f$ is the vector of variables, \f$ z_t \f$ is the vector of exogenous shocks,
 * and \f$ \eta_t = x_t - E_{t-1} x_t \f$ is the vector of endogenous expectation errors.
 *
 * The solution is: \f$ x_t = D_{sol} + G_{sol} x_{t-1} + M_{sol} z_t \f$
 *
 * @param mdGsol  (output) Solution matrix \f$ G_{sol} \f$
 * @param mdMsol  (output) Solution matrix \f$ M_{sol} \f$
 * @param vDsol   (output) Solution vector \f$ D_{sol} \f$
 * @param mdGamma0 (input) Coefficient matrix \f$ \Gamma_0 \f$, size nvar x nvar
 * @param mdGamma1 (input) Coefficient matrix \f$ \Gamma_1 \f$, size nvar x nvar
 * @param mdPsi    (input) Coefficient matrix \f$ \Psi \f$, size nvar x nshocks
 * @param mdPi     (input) Coefficient matrix \f$ \Pi \f$, size nvar x nendo
 * @param vdC      (input) Constant vector \f$ C \f$, size nvar
 * @return 0 for a unique solution, > 0 for indeterminacy (number of loose endogenous variables), -1 for no solution
 */
int gensys(Eigen::MatrixXd& mdGsol, Eigen::MatrixXd& mdMsol, Eigen::VectorXd& vDsol, const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC);

/**
 * @brief Solve a linear DSGE model using Sims (2002) and also return the full QZ decomposition details.
 *
 * Same as gensys() but additionally returns the generalized Schur matrices, eigenvalues,
 * and the number of unstable roots for inspection.
 *
 * @param mdGsol    (output) Solution matrix \f$ G_{sol} \f$
 * @param mdMsol    (output) Solution matrix \f$ M_{sol} \f$
 * @param vdDsol    (output) Solution vector \f$ D_{sol} \f$
 * @param mcdS      (output) Upper quasi-triangular Schur matrix S from QZ decomposition
 * @param mcdT      (output) Upper quasi-triangular Schur matrix T from QZ decomposition
 * @param mcdQ      (output) Unitary matrix Q from QZ decomposition
 * @param mcdZ      (output) Unitary matrix Z from QZ decomposition
 * @param vdLambda  (output) Generalized eigenvalues |alpha/beta|
 * @param nunstable (output) Number of unstable (explosive) eigenvalues
 * @param mdGamma0  (input)  Coefficient matrix \f$ \Gamma_0 \f$, size nvar x nvar
 * @param mdGamma1  (input)  Coefficient matrix \f$ \Gamma_1 \f$, size nvar x nvar
 * @param mdPsi     (input)  Coefficient matrix \f$ \Psi \f$, size nvar x nshocks
 * @param mdPi      (input)  Coefficient matrix \f$ \Pi \f$, size nvar x nendo
 * @param vdC       (input)  Constant vector \f$ C \f$, size nvar
 * @return 0 for a unique solution, > 0 for indeterminacy (number of loose endogenous variables), -1 for no solution
 */
int gensys_qzdetails(Eigen::MatrixXd& mdGsol, Eigen::MatrixXd& mdMsol, Eigen::VectorXd& vdDsol, 
     Eigen::MatrixXcd& mcdS, Eigen::MatrixXcd& mcdT, Eigen::MatrixXcd& mcdQ, Eigen::MatrixXcd& mcdZ, 
     Eigen::VectorXd& vdLambda, int& nunstable,
     const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, 
     const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC);


/**
 * @brief
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
 * @return 0 for a unique solution, > 0 for indeterminacy (number of loose endogenous variables), -1 for no solution
 */
int checksys(const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC);

/**
 * @brief Check if the column space of B is a subset of the column space of A.
 *
 * Uses a singular value decomposition of A to determine its effective column space.
 * The left singular vectors corresponding to singular values below a numerical
 * tolerance span the left null space of A (the orthogonal complement of col(A)).
 * col(B) ⊆ col(A) if and only if U_2^H * B ≈ 0, where U_2 are those vectors.
 *
 * @param A Input matrix whose column space is tested against.
 * @param B Input matrix whose columns are tested for membership in col(A).
 * @return true if every column of B lies in the column space of A (within numerical tolerance), false otherwise.
 */
bool check_colspace(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B);

/**
 * @brief Check if the column space of B is a subset of the column space of A.
 *
 * Uses a singular value decomposition of A to determine its effective column space.
 * The left singular vectors corresponding to singular values below a numerical
 * tolerance span the left null space of A (the orthogonal complement of col(A)).
 * col(B) ⊆ col(A) if and only if U_2^H * B ≈ 0, where U_2 are those vectors.
 *
 * @param A Input matrix whose column space is tested against.
 * @param B Input matrix whose columns are tested for membership in col(A).
 * @param tol Numerical tolerance for determining the rank of A and the effective column space.
 * @return true if every column of B lies in the column space of A (within numerical tolerance), false otherwise.
 */
bool check_colspace(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, double tol);

#endif
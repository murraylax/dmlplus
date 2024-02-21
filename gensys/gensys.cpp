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

#include "gensys.h"
#include "qz.h"

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
 * @param mdGsol (output) Eigen::MatrixXd with the solution matrix G_sol
 * @param mdMsol (output) Eigen::MatrixXd with the solution matrix M_sol
 * @param vcDsol (output) Eigen::VectorXd with the solution vector D_sol
 * 
 * @param mdGamma0 (input) const Eigen::MatrixXd with the coefficient matrix \Gamma_0
 * @param mdGamma1 (input) const Eigen::MatrixXd with the coefficient matrix \Gamma_1
 * @param mdPsi (input) const Eigen::MatrixXd with the coefficient matrix \Psi
 * @param mdPi (input) const Eigen::MatrixXd with the coefficient matrix \Pi
 * @param vdC (input) const Eigen::VectorXd with the constant vector C
 * 
 * @return Integer equal to the number of loose parameters. =0 when there is a unique solution, =-1 for no solution, >0 for indeterminacy
 * 
*/
int gensys(Eigen::MatrixXd& mdGsol, Eigen::MatrixXd& mdMsol, Eigen::VectorXd& vdDsol, const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC) {
    int nvar = mdGamma0.rows();
    int nshocks = mdPsi.cols();
    int nendo = mdPi.cols();

    // Convert input matrices to complex matricies
    Eigen::MatrixXcd mcdGamma0 = mdGamma0.cast<std::complex<double>>();
    Eigen::MatrixXcd mcdGamma1 = mdGamma1.cast<std::complex<double>>();
    Eigen::MatrixXcd mcdPsi = mdPsi.cast<std::complex<double>>();
    Eigen::MatrixXcd mcdPi = mdPi.cast<std::complex<double>>();
    Eigen::VectorXcd vcdC = vdC.cast<std::complex<double>>();

    Eigen::MatrixXcd mcdS(nvar,nvar);
    Eigen::MatrixXcd mcdT(nvar,nvar);
    Eigen::MatrixXcd mcdQ(nvar,nvar);
    Eigen::MatrixXcd mcdZ(nvar,nvar);
    Eigen::VectorXd vdLambda(nvar);

    int nstable = qz(mcdQ, mcdZ, mcdS, mcdT, vdLambda, mdGamma0, mdGamma1);
    int nunstable = nvar - nstable;

    Eigen::MatrixXcd mcdQt = mcdQ.adjoint();
    Eigen::MatrixXcd mcdPsi_tilde = mcdQt * mcdPsi;
    Eigen::MatrixXcd mcdPi_tilde = mcdQt * mcdPi;
    Eigen::VectorXcd vcdC_tilde = mcdQt * vcdC;

    Eigen::MatrixXcd mcdS11 = mcdS.block(0,0,nstable,nstable);
    Eigen::MatrixXcd mcdS12 = mcdS.block(0,nstable,nstable,nunstable);
    Eigen::MatrixXcd mcdS22 = mcdS.block(nstable,nstable,nunstable,nunstable);
    Eigen::MatrixXcd mcdT11 = mcdT.block(0,0,nstable,nstable);
    Eigen::MatrixXcd mcdT22 = mcdT.block(nstable,nstable,nunstable,nunstable);
    Eigen::MatrixXcd mcdPsi1 = mcdPsi_tilde.block(0,0,nstable,nshocks);
    Eigen::MatrixXcd mcdPsi2 = mcdPsi_tilde.block(nstable,0,nunstable,nshocks);
    Eigen::MatrixXcd mcdPi1 = mcdPi_tilde.block(0,0,nstable,nendo);
    Eigen::MatrixXcd mcdPi2 = mcdPi_tilde.block(nstable,0,nunstable,nendo);
    Eigen::VectorXcd vcdC1 = vcdC_tilde.segment(0,nstable);
    Eigen::VectorXcd vcdC2 = vcdC_tilde.segment(nstable, nunstable);
    Eigen::MatrixXcd mcdZ1 = mcdZ.block(0,0,nvar,nstable);
    Eigen::MatrixXcd mcdZ2 = mcdZ.block(0,nstable,nvar,nunstable);

    // Check the rank of mdPi2. For a unique solution, it should be of rank nendo
    int nrank_pi2 = Eigen::FullPivLU<Eigen::MatrixXcd>(mcdPi2).rank();
    int nloose = 0;

    // Check if there is either indeterminacy or no solution, and if so, figure out which, and get out of this function.
    if(nrank_pi2 != nendo) { 
        // Make an augmented matrix [mcdPi2 mcdPsi2 vcdC2]
        Eigen::MatrixXcd mcdCombined(nunstable, nendo+nshocks+1);
        mcdCombined.block(0, 0, nunstable, nendo) = mcdPi2;
        mcdCombined.col(nendo) = vcdC2;
        mcdCombined.block(0, nendo+1, nunstable, nshocks) = mcdPsi2;
        int nrank_combined = Eigen::FullPivLU<Eigen::MatrixXcd>(mcdCombined).rank();
        if(nrank_combined > nrank_pi2) {
            nloose = -1;
        } else {
            nloose = nendo - nrank_pi2;
        }
        if(nloose<0) nloose = -1;
        return nloose; // Leave the function, nothing more to do here
    } 

    // G = Z_1 S_{11}^{-1} T_{11} Z_1' 
    Eigen::MatrixXcd mcdGsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve(mcdT11 * mcdZ1.adjoint());

    // M = Z_1 S_{11}^{-1} ( \tilde{\Psi}_1 - \tilde{\Pi_1} \tilde{\Pi}_2^{-1} \tilde{\Psi}_2)
    Eigen::MatrixXcd mcdMsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( mcdPsi1 - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(mcdPsi2) );

    // D = Z_1 * S_{11}^{-1} (  \tilde{C}_1 - \tilde{\Pi_1} \tilde{\Pi}_2^{-1} \tilde{C}_2)
    Eigen::VectorXcd vcdDsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( vcdC1 - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(vcdC2)) ;

    mdGsol = mcdGsol.real();
    mdMsol = mcdMsol.real();
    vdDsol = vcdDsol.real();

    return 0;
}

int gensys_qzdetails(Eigen::MatrixXd& mdGsol, Eigen::MatrixXd& mdMsol, Eigen::VectorXd& vdDsol, 
     Eigen::MatrixXcd& mcdS, Eigen::MatrixXcd& mcdT, Eigen::MatrixXcd& mcdQ, Eigen::MatrixXcd& mcdZ, 
     Eigen::VectorXd& vdLambda, int& nunstable,
     const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, 
     const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC) {

    int nvar = mdGamma0.rows();
    int nshocks = mdPsi.cols();
    int nendo = mdPi.cols();

    // Convert input matrices to complex matricies
    Eigen::MatrixXcd mcdGamma0 = mdGamma0.cast<std::complex<double>>();
    Eigen::MatrixXcd mcdGamma1 = mdGamma1.cast<std::complex<double>>();
    Eigen::MatrixXcd mcdPsi = mdPsi.cast<std::complex<double>>();
    Eigen::MatrixXcd mcdPi = mdPi.cast<std::complex<double>>();
    Eigen::VectorXcd vcdC = vdC.cast<std::complex<double>>();

    int nstable = qz(mcdQ, mcdZ, mcdS, mcdT, vdLambda, mdGamma0, mdGamma1);
    nunstable = nvar - nstable;

    Eigen::MatrixXcd mcdQt = mcdQ.adjoint();
    Eigen::MatrixXcd mcdPsi_tilde = mcdQt * mcdPsi;
    Eigen::MatrixXcd mcdPi_tilde = mcdQt * mcdPi;
    Eigen::VectorXcd vcdC_tilde = mcdQt * vcdC;

    Eigen::MatrixXcd mcdS11 = mcdS.block(0,0,nstable,nstable);
    Eigen::MatrixXcd mcdS12 = mcdS.block(0,nstable,nstable,nunstable);
    Eigen::MatrixXcd mcdS22 = mcdS.block(nstable,nstable,nunstable,nunstable);
    Eigen::MatrixXcd mcdT11 = mcdT.block(0,0,nstable,nstable);
    Eigen::MatrixXcd mcdT22 = mcdT.block(nstable,nstable,nunstable,nunstable);
    Eigen::MatrixXcd mcdPsi1 = mcdPsi_tilde.block(0,0,nstable,nshocks);
    Eigen::MatrixXcd mcdPsi2 = mcdPsi_tilde.block(nstable,0,nunstable,nshocks);
    Eigen::MatrixXcd mcdPi1 = mcdPi_tilde.block(0,0,nstable,nendo);
    Eigen::MatrixXcd mcdPi2 = mcdPi_tilde.block(nstable,0,nunstable,nendo);
    Eigen::VectorXcd vcdC1 = vcdC_tilde.segment(0,nstable);
    Eigen::VectorXcd vcdC2 = vcdC_tilde.segment(nstable, nunstable);
    Eigen::MatrixXcd mcdZ1 = mcdZ.block(0,0,nvar,nstable);
    Eigen::MatrixXcd mcdZ2 = mcdZ.block(0,nstable,nvar,nunstable);

    // Check the rank of mdPi2. For a unique solution, it should be of rank nendo
    int nrank_pi2 = Eigen::FullPivLU<Eigen::MatrixXcd>(mcdPi2).rank();
    int nloose = 0;

    // Check if there is either indeterminacy or no solution, and if so, figure out which, and get out of this function.
    if(nrank_pi2 != nendo) { 
        // Make an augmented matrix [mcdPi2 mcdPsi2 vcdC2]
        Eigen::MatrixXcd mcdCombined(nunstable, nendo+nshocks+1);
        mcdCombined.block(0, 0, nunstable, nendo) = mcdPi2;
        mcdCombined.col(nendo) = vcdC2;
        mcdCombined.block(0, nendo+1, nunstable, nshocks) = mcdPsi2;
        int nrank_combined = Eigen::FullPivLU<Eigen::MatrixXcd>(mcdCombined).rank();
        if(nrank_combined > nrank_pi2) {
            nloose = -1;
        } else {
            nloose = nendo - nrank_pi2;
        }
        if(nloose<0) nloose = -1;
        return nloose; // Leave the function, nothing more to do here
    } 

    // G = Z_1 S_{11}^{-1} T_{11} Z_1' 
    Eigen::MatrixXcd mcdGsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve(mcdT11 * mcdZ1.adjoint());

    // M = Z_1 S_{11}^{-1} ( \tilde{\Psi}_1 - \tilde{\Pi_1} \tilde{\Pi}_2^{-1} \tilde{\Psi}_2)
    Eigen::MatrixXcd mcdMsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( mcdPsi1 - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(mcdPsi2) );

    // D = Z_1 * S_{11}^{-1} (  \tilde{C}_1 - \tilde{\Pi_1} \tilde{\Pi}_2^{-1} \tilde{C}_2)
    Eigen::VectorXcd vcdDsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( vcdC1 - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(vcdC2)) ;

    mdGsol = mcdGsol.real();
    mdMsol = mcdMsol.real();
    vdDsol = vcdDsol.real();

    return nloose;
}


/**
 * checksys(Eigen::MatrixXd& mdGsol, Eigen::MatrixXd& mdMsol, Eigen::MatrixXd& Dsol, const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC)
 * 
 * Use Sims (2002) method for checking for a solution to a linear dynamic general equilibrium model 
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
 * The method returns a value for existence and uniqueness of solution:
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
 * @return Integer equal to the number of loose parameters. =0 when there is a unique solution, =-1 for no solution, >0 for indeterminacy
 * 
*/
int checksys(const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC) {
    int nvar = mdGamma0.rows();
    int nshocks = mdPsi.cols();
    int nendo = mdPi.cols();

    // Convert input matrices to complex matrices
    Eigen::MatrixXcd mcdGamma0 = mdGamma0.cast<std::complex<double>>();
    Eigen::MatrixXcd mcdGamma1 = mdGamma0.cast<std::complex<double>>();
    Eigen::MatrixXcd mcdPsi = mdPsi.cast<std::complex<double>>();
    Eigen::MatrixXcd mcdPi = mdPi.cast<std::complex<double>>();
    Eigen::VectorXcd vcdC = vdC.cast<std::complex<double>>();

    Eigen::MatrixXcd mcdS(nvar,nvar);
    Eigen::MatrixXcd mcdT(nvar,nvar);
    Eigen::MatrixXcd mcdQ(nvar,nvar);
    Eigen::MatrixXcd mcdZ(nvar,nvar);
    Eigen::VectorXd vdLambda(nvar);

    int nstable = qz(mcdQ, mcdZ, mcdS, mcdT, vdLambda, mdGamma0, mdGamma1);
    int nunstable = nvar - nstable;

    Eigen::MatrixXcd mcdQt = mcdQ.adjoint();
    Eigen::MatrixXcd mcdPsi_tilde = mcdQt * mcdPsi;
    Eigen::MatrixXcd mcdPi_tilde = mcdQt * mcdPi;
    Eigen::VectorXcd vcdC_tilde = mcdQt * vcdC;

    Eigen::MatrixXcd mcdPsi2 = mcdPsi_tilde.block(nstable,0,nvar-nstable,nshocks);
    Eigen::MatrixXcd mcdPi2 = mcdPi_tilde.block(nstable,0,nvar-nstable,nendo);
    Eigen::VectorXcd vcdC2 = vcdC_tilde.segment(nstable, nunstable);

    // Check the rank of mdPi2. For a unique solution, it should be of rank nendo
    int nrank_pi2 = Eigen::FullPivLU<Eigen::MatrixXcd>(mcdPi2).rank();
    int nloose = 0;

    // Check if there is either indeterminacy or no solution, and if so, figure out which, and get out of this function.
    if(nrank_pi2 != nendo) { 
        // Make an augmented matrix [mcdPi2 mcdPsi2 vcdC2]
        Eigen::MatrixXcd mcdCombined(nunstable, nendo+nshocks+1);
        mcdCombined.block(0, 0, nunstable, nendo) = mcdPi2;
        mcdCombined.col(nendo) = vcdC2;
        mcdCombined.block(0, nendo+1, nunstable, nshocks) = mcdPsi2;
        int nrank_combined = Eigen::FullPivLU<Eigen::MatrixXcd>(mcdCombined).rank();
        if(nrank_combined > nrank_pi2) {
            nloose = -1;
        } else {
            nloose = nendo - nrank_pi2;
        }
        if(nloose<0) nloose = -1;
    } 

    return nloose;
}

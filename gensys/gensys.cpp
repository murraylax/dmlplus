/*
 * Copyright 2026 James M. Murray
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
#include <fstream>
#include <iostream>

// Write impulse response functions to a text file, with a header row of variable names.
void write_irf(const Eigen::MatrixXd& mdIRF, const std::vector<std::string>& varnames, const std::string& filepath) {
    std::ofstream irffile(filepath);
    size_t nvar = mdIRF.cols() - 1;
    for(size_t v=0; v<nvar; v++) {
        irffile << varnames[v] << "  ";
    }
    irffile << "Time" << std::endl;

    irffile << mdIRF;
    irffile.close();

    return;
}

// Compute impulse responses for a single shock over nirf periods. Returns an (nirf x nvar+2) matrix;
// the second-to-last column is the shock index, the last column is the time period.
Eigen::MatrixXd gensys_irf(const Eigen::MatrixXd& mdG, const Eigen::MatrixXd& mdM, double fshock, size_t shock_idx, size_t nirf) {
    size_t nvar = mdG.rows();
    size_t nshocks = mdM.cols();
    Eigen::MatrixXd mdIRF(nirf, nvar+2);
    Eigen::VectorXd vdIRF0(nvar);
    Eigen::VectorXd vdIRF1(nvar);
    Eigen::VectorXd vdZ(nshocks);
    vdZ.setZero();
    vdZ(shock_idx) = fshock;

    // Time t=0
    vdIRF0 = mdM * vdZ;
    mdIRF.row(0).segment(0, nvar) = vdIRF0.transpose();
    mdIRF(0,nvar) = shock_idx; 
    mdIRF(0,nvar+1) = 0; // Time period
    
    // All other t
    for(int t=1; t<nirf; t++) {
        vdIRF1 = mdG * vdIRF0;
        mdIRF.row(t).segment(0, nvar) = vdIRF1.transpose();
        mdIRF(t,nvar) = shock_idx;
        mdIRF(t,nvar+1) = t; // Time period
        vdIRF0 = vdIRF1;
    }

    return mdIRF;
}

// Compute impulse responses for all shocks (magnitude 0.01) over nirf periods.
// Returns a stacked ((nirf*nshocks) x nvar+2) matrix; second-to-last column is shock index, last is time period.
Eigen::MatrixXd gensys_irf(const Eigen::MatrixXd& mdG, const Eigen::MatrixXd& mdM, size_t nirf) {
    size_t nvar = mdG.rows();
    size_t nshocks = mdM.cols();
    Eigen::MatrixXd mdIRF_s(nirf, nvar+2);
    Eigen::MatrixXd mdIRF_all(nirf*nshocks, nvar+2);
    double fshock = 0.01;

    for(size_t s=0; s<nshocks; s++) {
        mdIRF_s = gensys_irf(mdG, mdM, fshock, s, nirf);
        mdIRF_all.block(s*nirf, 0, nirf, nvar+2) = mdIRF_s;
    }

    return mdIRF_all;
}


// Write impulse response functions to a CSV file with a header row of variable names, shock names, and a description column.
void write_irf_to_csvfile(const Eigen::MatrixXd& mdIRF, std::vector<std::string>& varnames, std::vector<std::string>& shocknames, std::string& desc, std::ofstream& csvfile) {

    size_t nvar = mdIRF.cols() - 2;
    size_t nrow = mdIRF.rows();

    for(int i=0; i<nvar; i++) {
        csvfile << varnames[i] << ", ";
    }
    csvfile << "Shock, Time, Description\n";

    for(size_t r=0; r<nrow; r++) {
        for(size_t i=0; i<nvar; i++) {
            csvfile << mdIRF(r,i);
            if(i==(nvar-1)) {
                csvfile << ", \"" << shocknames[(size_t)(mdIRF(r,nvar))] << "\", " << (size_t)(mdIRF(r,nvar+1)) << ", " << "\"" << desc <<"\"\n";
            } else {
                csvfile << ", ";
            }
        }
    }

    return;
}

// Write impulse response functions to a CSV file with no header row.
void write_irf_to_csvfile(const Eigen::MatrixXd& mdIRF, std::vector<std::string>& shocknames, std::string& desc, std::ofstream& csvfile) {

    size_t nvar = mdIRF.cols() - 2;
    size_t nrow = mdIRF.rows();

    for(size_t r=0; r<nrow; r++) {
        for(size_t i=0; i<nvar; i++) {
            csvfile << mdIRF(r,i);
            if(i==(nvar-1)) {
                csvfile << ", \"" << shocknames[(size_t)(mdIRF(r,nvar))] << "\", " << (size_t)(mdIRF(r,nvar+1)) << ", " << "\"" << desc <<"\"\n";
            } else {
                csvfile << ", ";
            }
        }
    }

    return;
}

// Solve the linear DSGE model using the Sims (2002) QZ method. Returns 0 for a unique solution,
// -1 for no solution, or a positive integer for the number of loose endogenous variables (indeterminacy).
int gensys(Eigen::MatrixXd& mdGsol, Eigen::MatrixXd& mdMsol, Eigen::VectorXd& vdDsol, const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC) {
    double tol = 1e-10; // Tolerance for checking rank and indeterminacy
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

    // Partition Q^H by rows: bottom nunstable rows correspond to unstable eigenvalues
    Eigen::MatrixXcd mcdQt   = mcdQ.adjoint();
    Eigen::MatrixXcd mcdQ2   = mcdQt.block(nstable, 0, nunstable, nvar);
    Eigen::MatrixXcd mcdQ2Pi  = mcdQ2 * mcdPi;   // nunstable x nendo
    Eigen::MatrixXcd mcdQ2Psi = mcdQ2 * mcdPsi;  // nunstable x nshocks

    // SVD on Q2 * Pi
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(mcdQ2Pi, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXcd U = svd.matrixU();   // nunstable x nunstable
    Eigen::MatrixXcd V = svd.matrixV();   // nendo x nendo
    Eigen::VectorXd s = svd.singularValues();  // min(nunstable, nendo) singular values    
    int nrank = (s.array() > tol).count();

    if(nrank < nunstable) {
        bool indeterminacy = check_colspace(mcdQ2Pi, mcdQ2Psi);
        if (indeterminacy) {   
            return nunstable - nrank; // Indeterminacy, return number of loose endogenous variables
        } else {
            return -1; // No solution
        }
    }




    /* I think all this is wrong 
    Eigen::MatrixXcd mcdQt = mcdQ.adjoint();
    Eigen::MatrixXcd mcdPsi_tilde = mcdQt * mcdPsi;
    Eigen::MatrixXcd mcdPi_tilde = mcdQt * mcdPi;
    Eigen::VectorXcd vcdC_tilde = mcdQt * vcdC;

    Eigen::MatrixXcd mcdS11 = mcdS.block(0,0,nstable,nstable);
    Eigen::MatrixXcd mcdS12 = mcdS.block(0,nstable,nstable,nunstable);
    Eigen::MatrixXcd mcdS22 = mcdS.block(nstable,nstable,nunstable,nunstable);
    Eigen::MatrixXcd mcdT11 = mcdT.block(0,0,nstable,nstable);
    Eigen::MatrixXcd mcdT12 = mcdT.block(0,nstable,nstable,nunstable);
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

    // Forward solution for unstable block: phi_z = -Pi2^{-1} Psi2
    Eigen::MatrixXcd mcdPhi_z = -mcdPi2.colPivHouseholderQr().solve(mcdPsi2);

    // Forward solution for unstable block constant: phi_c = (S22 - T22)^{-1} C2
    Eigen::MatrixXcd mcdS22_minus_T22 = mcdS22 - mcdT22;
    Eigen::VectorXcd vcdPhi_c = mcdS22_minus_T22.colPivHouseholderQr().solve(vcdC2);

    // G = Z_1 S_{11}^{-1} (T_{11} Z_1' + T_{12} Z_2')
    Eigen::MatrixXcd mcdGsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve(mcdT11 * mcdZ1.adjoint() + mcdT12 * mcdZ2.adjoint());

    // M = Z_1 S_{11}^{-1} (Psi1 - S12 * phi_z - Pi1 * Pi2^{-1} * Psi2) + Z_2 * phi_z
    Eigen::MatrixXcd mcdMsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( 
        mcdPsi1 - mcdS12 * mcdPhi_z - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(mcdPsi2) ) 
        + mcdZ2 * mcdPhi_z;

    // D = Z_1 * S_{11}^{-1} (C1 - S12 * phi_c - Pi1 * Pi2^{-1} * C2) + Z_2 * phi_c
    Eigen::VectorXcd vcdDsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( 
        vcdC1 - mcdS12 * vcdPhi_c - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(vcdC2) ) 
        + mcdZ2 * vcdPhi_c;

    mdGsol = mcdGsol.real();
    mdMsol = mcdMsol.real();
    vdDsol = vcdDsol.real();

    // // G = Z_1 S_{11}^{-1} T_{11} Z_1' 
    // Eigen::MatrixXcd mcdGsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve(mcdT11 * mcdZ1.adjoint());

    // // M = Z_1 S_{11}^{-1} ( \tilde{\Psi}_1 - \tilde{\Pi_1} \tilde{\Pi}_2^{-1} \tilde{\Psi}_2)
    // Eigen::MatrixXcd mcdMsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( mcdPsi1 - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(mcdPsi2) );

    // // D = Z_1 * S_{11}^{-1} (  \tilde{C}_1 - \tilde{\Pi_1} \tilde{\Pi}_2^{-1} \tilde{C}_2)
    // Eigen::VectorXcd vcdDsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( vcdC1 - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(vcdC2)) ;

    // mdGsol = mcdGsol.real();
    // mdMsol = mcdMsol.real();
    // vdDsol = vcdDsol.real();

    return 0;
    */
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
    Eigen::MatrixXcd mcdT12 = mcdT.block(0,nstable,nstable,nunstable);
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

    // Forward solution for unstable block: phi_z = -Pi2^{-1} Psi2
    Eigen::MatrixXcd mcdPhi_z = -mcdPi2.colPivHouseholderQr().solve(mcdPsi2);

    // Forward solution for unstable block constant: phi_c = (S22 - T22)^{-1} C2
    Eigen::MatrixXcd mcdS22_minus_T22 = mcdS22 - mcdT22;
    Eigen::VectorXcd vcdPhi_c = mcdS22_minus_T22.colPivHouseholderQr().solve(vcdC2);

    // G = Z_1 S_{11}^{-1} (T_{11} Z_1' + T_{12} Z_2')
    Eigen::MatrixXcd mcdGsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve(mcdT11 * mcdZ1.adjoint() + mcdT12 * mcdZ2.adjoint());

    // M = Z_1 S_{11}^{-1} (Psi1 - S12 * phi_z - Pi1 * Pi2^{-1} * Psi2) + Z_2 * phi_z
    Eigen::MatrixXcd mcdMsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( 
        mcdPsi1 - mcdS12 * mcdPhi_z - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(mcdPsi2) ) 
        + mcdZ2 * mcdPhi_z;

    // D = Z_1 * S_{11}^{-1} (C1 - S12 * phi_c - Pi1 * Pi2^{-1} * C2) + Z_2 * phi_c
    Eigen::VectorXcd vcdDsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( 
        vcdC1 - mcdS12 * vcdPhi_c - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(vcdC2) ) 
        + mcdZ2 * vcdPhi_c;

    mdGsol = mcdGsol.real();
    mdMsol = mcdMsol.real();
    vdDsol = vcdDsol.real();

    // // G = Z_1 S_{11}^{-1} T_{11} Z_1' 
    // Eigen::MatrixXcd mcdGsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve(mcdT11 * mcdZ1.adjoint());

    // // M = Z_1 S_{11}^{-1} ( \tilde{\Psi}_1 - \tilde{\Pi_1} \tilde{\Pi}_2^{-1} \tilde{\Psi}_2)
    // Eigen::MatrixXcd mcdMsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( mcdPsi1 - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(mcdPsi2) );

    // // D = Z_1 * S_{11}^{-1} (  \tilde{C}_1 - \tilde{\Pi_1} \tilde{\Pi}_2^{-1} \tilde{C}_2)
    // Eigen::VectorXcd vcdDsol = mcdZ1 * mcdS11.colPivHouseholderQr().solve( vcdC1 - mcdPi1 * mcdPi2.colPivHouseholderQr().solve(vcdC2)) ;

    // mdGsol = mcdGsol.real();
    // mdMsol = mcdMsol.real();
    // vdDsol = vcdDsol.real();

    return nloose;
}


// Check existence and uniqueness of the DSGE solution without computing G, M, or D.
// Returns 0 (unique), -1 (no solution), or a positive integer (number of loose endogenous variables).
int checksys(const Eigen::MatrixXd& mdGamma0, const Eigen::MatrixXd& mdGamma1, const Eigen::MatrixXd& mdPsi, const Eigen::MatrixXd& mdPi, const Eigen::VectorXd& vdC) {
    int nvar = mdGamma0.rows();
    int nshocks = mdPsi.cols();
    int nendo = mdPi.cols();

    // Convert input matrices to complex matrices
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

// Check whether col(B) ⊆ col(A) using SVD of A; returns true if the projection of B onto the left null space of A is ~zero.
bool check_colspace(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, double tol) {

    // Compute SVD of A: A = U * S * V^H
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeFullU);
    Eigen::MatrixXcd U = svd.matrixU();
    Eigen::VectorXd singularValues = svd.singularValues();

    // Determine the rank of A: number of singular values above tolerance
    int rank = 0;
    for(int i = 0; i < singularValues.size(); i++) {
        if(singularValues(i) > tol) rank++;
    }

    // If A is full rank, col(B) is trivially a subset of col(A)
    if(rank == A.rows()) return true;

    // U_2: left singular vectors spanning the left null space of A (columns beyond `rank`)
    // These are orthogonal to col(A), so col(B) ⊆ col(A) iff U_2^H * B ≈ 0
    Eigen::MatrixXcd U_2 = U.rightCols(U.cols() - rank);

    return (U_2.adjoint() * B).norm() < tol;
}

// Overload with default tolerance of 1e-10.
bool check_colspace(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B) {
    double tol = 1e-10;
    return check_colspace(A, B, tol);
}
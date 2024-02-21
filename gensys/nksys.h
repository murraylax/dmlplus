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
 * @file nksys.h
 *
 * @brief 
 * This is a New Keynesian DSGE model with two types of agents, optimizing and rule-of-thumb
 * 
 * @author James Murray
 * Contact: james@murraylax.org
*/


#ifndef _INCL_NKSYS
#define _INCL_NKSYS

#include <Eigen/Dense>

// Indices into the vector for the state variables
enum var {
    i_cO,
    i_cR,
    i_c,
    i_NO,
    i_NR,
    i_N,
    i_bO,
    i_b,
    i_iO,
    i_i,
    i_kO,
    i_k,
    i_w,
    i_R,
    i_Rk,
    i_Pi,
    i_lambda,
    i_Q,
    i_Psi,
    i_pstar,
    i_y,
    i_ystar,
    i_pi,
    i_piA,
    i_piB,
    i_piF,
    i_tauO,
    i_tauR,
    i_tau,
    i_g,
    i_s,
    i_d,
    i_a,
    i_uI,
    i_uD,
    i_uC,
    i_Ea,
    i_EPi,
    i_Elambda,
    i_Ei,
    i_ERk,
    i_EQ,
    i_Epi,
    i_EpiF,
    i_lagcO,
    i_cum_g,
    i_cum_tauO,
    i_cum_tauR,
    i_cum_tau,
    i_cum_d,
    i_cum_y,
    nvar
};

// Indices for endogenous expectation errors
enum exp {
    i_a,
    i_Pi,
    i_lambda,
    i_i,
    i_Rk,
    i_Q,
    i_pi,
    i_piF,
    nexp
}

// Indices for exogenous shocks
enum shock {
    i_shock_a,
    i_shock_uD,
    i_shock_uI,
    i_shock_r,
    i_shock_uC,
    i_shock_fg,
    i_shock_fO,
    i_shock_fR,
    nshocks
}



#endif
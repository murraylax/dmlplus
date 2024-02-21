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
 * @file nksys.cpp
 *
 * @brief 
 * This is a New Keynesian DSGE model with two types of agents, optimizing and rule-of-thumb
 * 
 * @author James Murray
 * Contact: james@murraylax.org
*/

#include "nksys.h"
#include "gensys.h"

void nksys(Eigen::MatrixXd& mdGamma0, Eigen::MatrixXd& mdGamma1, Eigen::MatrixXd& mdPsi, Eigen::MatrixXd& mdPi) {
    int eq=0;

    // FOC labor and consumption - Equation D.1
    eq++;
    mdGamma0(eq, var::i_cO) = 1.0;
    mdGamma0(eq, var::i_w) = -1.0*(1.0 - parm_h / parm_gamma);
    mdGamma0(eq, var::i_NO) = parm_psi*(1.0 - parm_h / parm_gamma);
    mdGamma0(eq, var::i_a) = parm_h / parm_gamma;
    mdGamma0(eq, var::i_uD) = -1.0*(1.0 - parm_h / parm_gamma);
    

}



#endif
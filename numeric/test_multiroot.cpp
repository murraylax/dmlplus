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
 * @file multiroot.cpp
 *
 * @brief 
 * This is a numeric multivariate root finder using vectors/matrices defined by the Eigen package. 
 * It is a wrapper for the C gsl_multiroot functions from the GNU Scientific Library (https://www.gnu.org/software/gsl/doc/html/multiroots.html)
 * 
 * @author James Murray
 * Contact: jmurray@uwlax.edu
*/

#include <Eigen/Dense>
#include <iostream>
#include "dml_multiroot.h"

using namespace std;

Eigen::VectorXd myfunc(Eigen::VectorXd x, const void* params) {
    const Eigen::VectorXd* p = static_cast<const Eigen::VectorXd*>(params);

    Eigen::VectorXd result(x.size());
    int eq = -1;
    eq++;
    result(eq) = (*p)[0] * (1 - x[0]);
    eq++;
    result(eq) = (*p)[1] * (x[1] - x[0] * x[0]);


    return result;
}

int main() {
    Eigen::VectorXd params(2);
    params << 1, 10;
    Eigen::VectorXd x(2);
    x << 2, 2;

    Eigen::VectorXd sol = dml_multiroot(x, myfunc, static_cast<void*>(&params));
    cout << "Solution:\n" << sol << endl;

    return 0;
}
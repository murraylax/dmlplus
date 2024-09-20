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
 * @file test_multimin.cpp
 *
 * @brief 
 * Test out the numeric multidimensional minimizer dml_multimin
 * 
 * @author James Murray
 * Contact: jmurray@uwlax.edu
*/

#include <Eigen/Dense>
#include <iostream>
#include "dml_multimin.h"

using namespace std;

double myfunc(Eigen::VectorXd x, const void* params) {
    const Eigen::VectorXd* p = static_cast<const Eigen::VectorXd*>(params);

    Eigen::VectorXd result(x.size());
    int eq = -1;
    eq++;
    result(eq) = (*p)[0] * (1 - x[0]);
    eq++;
    result(eq) = (*p)[1] * (x[1] - x[0] * x[0]);

    double fmin = result.squaredNorm();

    return fmin;
}

int main() {
    Eigen::VectorXd params(2);
    params << 1, 10;
    Eigen::VectorXd x(2);
    x << 2, 2;

    Eigen::VectorXd lower_bounds(2);
    Eigen::VectorXd upper_bounds(2);
    lower_bounds << 0.5, 0.5;
    upper_bounds << 2.5, 2.5;


    Eigen::VectorXd sol = dml_multimin(x, lower_bounds, upper_bounds, myfunc, static_cast<void*>(&params), true);
    cout << "Solution:\n" << sol << endl;

    return 0;
}
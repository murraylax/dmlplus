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
 * @file utils.cpp
 *
 * @brief 
 * This is a collection of commonly used utilities for copying 
 * matrices and vectors from one format to another, including outputting 
 * text or R code to bring into R.
 * 
 * @author James Murray
 * Conact: james@murraylax.org
*/

#include "utils.h"

#include <sstream>
#include <iomanip>

using namespace std;
using namespace Eigen;

/**
 * copy_matrix_to_lapack_complex(lapack_complex_double* lapackMatrix, const Eigen::MatrixXd& matrix)
 * 
 * Copy the elements of an Eigen::MatrixXd to a lapack_complex_double* 
 * 
 * @param lapackMatrix (Output) The elements of this lapack_complex_double* matrix will be overwritten with elements of the `matrix` parameter
 * @param matrix This is a const Eigen::MatrixXd whose elements will be copies to the `lapackMatrix` parameter
*/
void copy_matrix_to_lapack_complex(lapack_complex_double* lapackMatrix, const Eigen::MatrixXd& matrix) {
    int rows = matrix.rows();
    int cols = matrix.cols();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Convert each real number in 'matrix' to a complex number and store it in the 'lapackMatrix'
            lapackMatrix[j * rows + i] = {matrix(i, j), 0.0};
        }
    }
}

/**
 * copy_lapack_complex_to_matrix(Eigen::MatrixXcd& complex_matrix, const lapack_complex_double* lapack_complex)
 * 
 * Copy the elements of a lapack_complex_double* to an Eigen::MatrixXd 
 * 
 * @param complex_matrix (Output) This is a Eigen::MatrixXcd whose elements will overwritten with elements of the `lapackMatrix` parameter
 * @param lapackMatrix This is a lapack_complex_double* array whose elements will be copies to the `complex_matrix` parameter
*/
void copy_lapack_complex_to_matrix(Eigen::MatrixXcd& complex_matrix, const lapack_complex_double* lapack_complex) {
    int rows = complex_matrix.rows();
    int cols = complex_matrix.cols();

    // Copy data from the lapack_complex_double array to the Eigen matrix
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Note: LAPACK uses column-major order, same as Eigen
            complex_matrix(i, j) = lapack_complex[j * rows + i];
        }
    }
}

/**
 * void copy_lapack_complex_to_vector(Eigen::VectorXcd& complex_vector, const lapack_complex_double* lapack_complex)
 * 
 * Copy the elements of a lapack_complex_double* to an Eigen::MatrixXd 
 * 
 * @param complex_vector (Output) This is a Eigen::VectorXcd whose elements will overwritten with elements of the `lapack_complex` parameter
 * @param lapack_complex This is a lapack_complex_double* array whose elements will be copies to the `complex_vector` parameter
*/
void copy_lapack_complex_to_vector(Eigen::VectorXcd& complex_vector, const lapack_complex_double* lapack_complex) {
    int n = complex_vector.size();

    for(int i = 0; i < n; i++) {
        complex_vector(i) = lapack_complex[i];
    }
}

/**
 * std::string print_vector_to_rcode(const Eigen::VectorXcd& v, const std::string& vecname)
 * 
 * Create R code to create a complex R matrix with the contents of a Eigen::VectorXcd
 * 
 * @param v A const Eigen::VectorXcd with the elements to be printed to R code
 * @param vecname A string with the name of the vector variable for the R code
*/
std::string print_vector_to_rcode(const Eigen::VectorXcd& v, const std::string& vecname) {
    int nrow = v.size();

    std::ostringstream oss;

    oss << vecname << " <- matrix(data = c(";

    for(int i=0; i<nrow; i++) {
        oss << v(i).real() << " + " << v(i).imag() << "i";
        if(i!=(nrow-1)) {
            oss << ", ";
        }          
    }
    
    oss << "), nrow = " << nrow << ", ncol = 1)" << endl;

    std::string ss = oss.str();
    return ss;
}

/**
 * @overload std::string print_vector_to_rcode(const Eigen::VectorXd& v, const std::string& vecname)
 * 
 * Create R code to create an R matrix with the contents of a real Eigen::VectorXd
 * 
 * @param v A const Eigen::VectorXd with the elements to be printed to R code
 * @param vecname A string with the name of the vector variable for the R code
*/
std::string print_vector_to_rcode(const Eigen::VectorXd& v, const std::string& vecname) {
    int nrow = v.size();

    std::ostringstream oss;

    oss << vecname << " <- matrix(data = c(";

    for(int i=0; i<nrow; i++) {
        oss << v(i);
        if(i!=(nrow-1)) {
            oss << ", ";
        }          
    }
    
    oss << "), nrow = " << nrow << ", ncol = 1)" << endl;

    std::string ss = oss.str();
    return ss;
}

/**
 * std::string print_matrix_to_rcode(const Eigen::MatrixXcd& m, const std::string& matname) 
 * 
 * Create R code to create an R matrix with the contents of a complex Eigen::MatrixXcd
 * 
 * @param m A const Eigen::MatrixXcd with the elements to be printed to R code
 * @param matname A string with the name of the matrix variable for the R code
*/
std::string print_matrix_to_rcode(const Eigen::MatrixXcd& m, const std::string& matname) {
    int nrow = m.rows();
    int ncol = m.cols();

    std::ostringstream oss;

    oss << matname << " <- matrix(data = c(";

    for(int j=0; j<ncol; j++) {
        for(int i=0; i<nrow; i++) {
            oss << m(i,j).real() << " + " << m(i,j).imag() << "i";
            if(j!=(ncol-1) || i!=(nrow-1)) {
                oss << ", ";
            }          
        }
    }
    oss << "), nrow = " << nrow << ", ncol = " << ncol << ")" << endl;

    std::string ss = oss.str();
    return ss;
}

/**
 * @overload std::string print_matrix_to_rcode(const Eigen::MatrixXd& m, const std::string& matname)
 * 
 * Create R code to create an R matrix with the contents of a Eigen::MatrixXd
 * 
 * @param m A const Eigen::MatrixXd with the elements to be printed to R code
 * @param matname A string with the name of the matrix variable for the R code
*/
std::string print_matrix_to_rcode(const Eigen::MatrixXd& m, const std::string& matname) {
    int nrow = m.rows();
    int ncol = m.cols();

    std::ostringstream oss;

    oss << matname << " <- matrix(data = c(";

    for(int j=0; j<ncol; j++) {
        for(int i=0; i<nrow; i++) {
            oss << m(i,j);
            if(j!=(ncol-1) || i!=(nrow-1)) {
                oss << ", ";
            }          
        }
    }
    oss << "), nrow = " << nrow << ", ncol = " << ncol << ")" << endl;

    std::string ss = oss.str();
    return ss;
}


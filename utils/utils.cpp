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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>

using namespace std;
using namespace Eigen;

/**
 * write_eigen_csv(Eigen::MatrixXd& mat, std::string& filepath)
 * 
 * Write an Eigen::MatrixXd to a csv file
 * 
 * @param mat Eigem::MatrixXd to write to file
 * @param filepath String that contains the path and filename to write to
 */
void write_eigen_csv(const Eigen::MatrixXd& mat, const std::string& filepath) {
    std::ofstream csvfile(filepath);
    size_t nvar = mat.cols();
    size_t nrow = mat.rows();

    for(int r=0; r<nrow; r++) {
        for(int i=0; i<nvar; i++) {
            csvfile << mat(r,i);
            if(i==(nvar-1)) {
                csvfile << "\n";
            } else {
                csvfile << ", ";
            }
        }
    }

    csvfile.close();
}

/**
 * write_eigen_csv(Eigen::VectorXd& vec, std::string& filepath)
 * 
 * Write an Eigen::VectorXd to a csv file
 * 
 * @param vec Eigen::VectorXd to write to file
 * @param filepath String that contains the path and filename to write to
 * 
 */
void write_eigen_csv(const Eigen::VectorXd& vec, const std::string& filepath) {
    std::ofstream csvfile(filepath);

    size_t nrow = vec.size();

    for(int r=0; r<nrow; r++) {
        csvfile << vec(r) << "\n";
    }

    csvfile.close();
}

/**
 * write_eigen_csv(Eigen::VectorXd& vec, std::string& varname, std::string& filepath)
 * 
 * Write an Eigen::VectorXd to a file, with a variable name at the top
 * 
 * @param vec Eigen::VectorXd to write to file
 * @param varname String that is the variable name for the vector
 * @param filepath String that contains the path and filename to write to
 */
void write_eigen_csv(Eigen::VectorXd& vec, std::string& varname, std::string& filepath) {
    std::ofstream csvfile(filepath);

    size_t nrow = vec.rows();

    csvfile << varname << "\n";
    for(int r=0; r<nrow; r++) {
        csvfile << vec(r) << "\n";
    }

    csvfile.close();
}

/**
 * write_eigen_csv(Eigen::MatrixXd& mat, std::string& filepath)
 * 
 * Write an Eigen::MatrixXd to a csv file, including variable names
 * 
 * @param mat Eigen::MatrixXd to write to file
 * @param varnames Vector of strings for the variable names, must have the same size as the number of columns of mat
 * @param filepath String that contains the path and filename to write to
 */
void write_eigen_csv(Eigen::MatrixXd& mat, std::vector<std::string>& varnames, std::string& filepath) {
    ofstream csvfile(filepath);
    size_t nvar = mat.cols();
    size_t nrow = mat.rows();
    for(int i=0; i<nvar; i++) {
        csvfile << "\"" << varnames[i] << "\"";
        if(i==(nvar=1)) {
            csvfile << "\n";
        } else {
            csvfile << ", ";
        }
    }
    for(int r=0; r<nrow; r++) {
        for(int i=0; i<nvar; i++) {
            csvfile << mat(r,i);
            if(i==(nvar-1)) {
                csvfile << "\n";
            } else {
                csvfile << ", ";
            }
        }
    }

    csvfile.close();
}

// Function to start the timer and return the start time point
std::chrono::time_point<std::chrono::steady_clock> start_timer() {
    return std::chrono::steady_clock::now();
}

// Function to stop the timer, calculate the elapsed time, and print it
void stop_timer(const std::chrono::time_point<std::chrono::steady_clock>& start_time) {
    auto end_time = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    int minutes = duration.count() / 60000;
    int seconds = (duration.count() % 60000) / 1000;
    int milliseconds = duration.count() % 1000;

    std::cout << "Elapsed time: " << minutes << " minute(s), "
              << seconds << " second(s), and "
              << milliseconds << " millisecond(s)." << std::endl;

    return;
}

/**
 * copy_vector_to_gsl(const Eigen::VectorXd& eigen_input, gsl_vector* gsl_output)
 * 
 * Copy the elements of an Eigen::VectorXd to a gsl_vector* 
 * 
 * @param gsl_output (Output) The elements of this gsl_vector* will be overwritten with elements of the `eigen_input` parameter
 * @param eigen_input This is a const Eigen::VectorXd whose elements will be copied to the `gsl_output` parameter
*/
void copy_vector_to_gsl(gsl_vector* gsl_output, const Eigen::VectorXd& eigen_input) {
    int eigen_size = eigen_input.size();
    int gsl_size = gsl_output->size;
    if(eigen_size != gsl_size) {
        throw std::invalid_argument("Size mismatch: Eigen::VectorXd and gsl_vector must have the same size.");
    }

    for(int i=0; i<eigen_input.size(); i++) {
        gsl_vector_set(gsl_output, i, eigen_input[i]);
    }
    return;
}

/**
 * copy_gsl_to_vector(Eigen::VectorXd& eigen_output, const gsl_vector* gsl_input)
 * 
 * Copy the elements of a gsl_vector* to an Eigen::VectorXd  
 * 
 * @param eigen_output (Output) The elements of this Eigen::VectorXd* vector will be overwritten with elements of the `gsl_input` parameter
 * @param gsl_input This is a const gsl_vector* whose elements will be copied to the `eigen_output` parameter
*/
void copy_gsl_to_vector(Eigen::VectorXd& eigen_output, const gsl_vector* gsl_input) {
    int eigen_size = eigen_output.size();
    int gsl_size = gsl_input->size;
    if(eigen_size != gsl_size) {
        throw std::invalid_argument("Size mismatch: Eigen::VectorXd and gsl_vector must have the same size.");
    }

    for(int i=0; i<eigen_output.size(); i++) {
        eigen_output[i] = gsl_vector_get(gsl_input, i);
    }
    return;
}

/**
 * copy_matrix_to_gsl(const Eigen::MatrixXd& eigen_input, gsl_matrix* gsl_output)
 * 
 * Copy the elements of an Eigen::MatrixXd to a gsl_matrix* 
 * 
 * @param gsl_output (Output) The elements of this gsl_matrix* will be overwritten with elements of the `eigen_input` parameter
 * @param eigen_input This is a const Eigen::MatrixXd whose elements will be copied to the `gsl_output` parameter
*/
void copy_matrix_to_gsl(gsl_matrix* gsl_output, const Eigen::MatrixXd& eigen_input) {
    int eigen_nrow = eigen_input.rows();
    int eigen_ncol = eigen_input.cols();
    int gsl_nrow = gsl_output->size1;
    int gsl_ncol = gsl_output->size2;
    if(eigen_nrow != gsl_nrow || eigen_ncol != gsl_ncol) {
        throw std::invalid_argument("Size mismatch: Eigen::MatrixXd and gsl_matrix* must have the same dimensions.");
    }

    for(int i = 0; i < eigen_nrow; i++) {
        for(int j = 0; j < eigen_ncol; j++) {
            gsl_matrix_set(gsl_output, i, j, eigen_input(i,j));
        }
    }
    return;
}

/**
 * copy_matrix_to_gsl(const Eigen::MatrixXd& eigen_input, gsl_matrix* gsl_output)
 * 
 * Copy the elements of an Eigen::MatrixXd to a gsl_matrix* 
 * 
 * @param eigen_output (Output) This is a Eigen::MatrixXd whose elements will be overwritten with elements of the `gsl_intput` parameter
 * @param gsl_input This is a const gsl_matrix* whose elements will be copied to the `eigen_output` parameter
*/
void copy_gsl_to_matrix(Eigen::MatrixXd& eigen_output, const gsl_matrix* gsl_input) {
    int eigen_nrow = eigen_output.rows();
    int eigen_ncol = eigen_output.cols();
    int gsl_nrow = gsl_input->size1;
    int gsl_ncol = gsl_input->size2;
    if(eigen_nrow != gsl_nrow || eigen_ncol != gsl_ncol) {
        throw std::invalid_argument("Size mismatch: Eigen::MatrixXd and gsl_matrix* must have the same dimensions.");
    }

    for(int i = 0; i < eigen_nrow; i++) {
        for(int j = 0; j < eigen_ncol; j++) {
            eigen_output(i,j) = gsl_matrix_get(gsl_input, i, j);
        }
    }
    return;
}

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


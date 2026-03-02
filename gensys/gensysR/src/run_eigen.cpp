#include "utils.h"
#include "qz.h"

#include <iostream>
#include <sstream>
#include <Eigen/Dense>
#include <iomanip>
#include <complex>
#include <lapacke.h>
#include <algorithm>

using namespace std;
using namespace Eigen;

void test_inverse() {
    MatrixXd mat(2,2);
    mat << 1, 2, 3, 4;
    cout << "Matrix mat: " << mat << endl;

    double fdet = mat.determinant();
    cout << "Determinant of mat: " << setprecision(2) << fdet << endl;

    MatrixXd invmat = mat.inverse();
    cout << "Inverse of max: " << setprecision(2) << invmat << endl;
}


int main() {
    int n = 3;
    MatrixXd mdA(n,n);
    MatrixXd mdB(n,n);

    MatrixXcd mcdS(n,n);
    MatrixXcd mcdT(n,n);
    MatrixXcd mcdQ(n,n);
    MatrixXcd mcdZ(n,n);
    VectorXd vdLambda(n);

    mdA << 1, 2, 13, 14, 15, 6, 7, 8, 9;
    mdB << 11, 22, 3, 4, 5, 6, 7, 18, 19;

    std::string sA = print_matrix_to_rcode(mdA, "A");
    std::string sB = print_matrix_to_rcode(mdB, "B");

    cout << sA << endl;
    cout << sB << endl;

    qz(mcdQ, mcdZ, mcdS, mcdT, vdLambda, mdA, mdB);

    std::string sQ = print_matrix_to_rcode(mcdQ, "Q");
    std::string sZ = print_matrix_to_rcode(mcdZ, "Z");
    std::string sS = print_matrix_to_rcode(mcdS, "S");
    std::string sT = print_matrix_to_rcode(mcdT, "T");
    std::string sLambda = print_vector_to_rcode(vdLambda, "lambda");
    
    cout << sQ << endl << sZ << endl << sS << endl << sT << endl << sLambda << endl;

    return 1;
}
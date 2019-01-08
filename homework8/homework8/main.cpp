//
//  main.cpp
//  homework8
//
//  Created by Jun Zhou on 11/1/18.
//  Copyright © 2018 jzhou. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include "newmat.h"
#include "newmatio.h"
#include <vector>

Matrix repeated_squaring(Matrix A, int exponent, int no_rows) {
    if (exponent == 0) {
        return IdentityMatrix(no_rows);
    }
    if (exponent % 2 != 0) {
        return A * repeated_squaring(A * A, (exponent - 1) / 2, no_rows);
    }
    return repeated_squaring(A * A, exponent / 2, no_rows);
}

Matrix direct_multiplication(Matrix A, int exponent) {
    Matrix C = A;
    for (int i = 0; i < exponent; i++) {
        A = A * C;
    }
    return A;
}

Matrix random_generator(int n) {
    Matrix A(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A.element(i, j) = 10 * (rand() / (double) RAND_MAX) - 5;
        }
    }
    return A;
}

void create_file(Matrix A, int no_rows) {
    vector <float> repeated_time;
    vector <float> direct_time;
    for (int i = 0; i <= 300; i++) {
        long time_before = clock();
        Matrix B = repeated_squaring(A, i, no_rows);
        long time_after = clock();
        float diff = ((float)time_after - (float)time_before);
        repeated_time.push_back(diff / CLOCKS_PER_SEC);
        time_before = clock();
        Matrix C = direct_multiplication(A, i);
        time_after = clock();
        diff = ((float)time_after - (float)time_before);
        direct_time.push_back(diff / CLOCKS_PER_SEC);
    }
    ofstream file("/Users/jzhou/Desktop/output.txt");
    for (int i = 0; i <= repeated_time.size(); i++) {
        file<< repeated_time[i] << ",";
        file<< direct_time[i] << endl;
    }
}

int main(int argc, const char * argv[]) {
    int no_row, exponent;
    sscanf (argv[1], "%d", &exponent);
    sscanf (argv[2], "%d", &no_row);
    cout << "The number of rows/columns in the square matrix is: " << no_row << endl;
    cout << "The exponent is: " << exponent << endl;
    Matrix A = random_generator(no_row);
    cout << "Original Matrix is: " << endl;
    cout << setw(5) << setprecision(1) << A << endl;
    long time_before = clock();
    Matrix B = repeated_squaring(A, exponent, no_row);
    long time_after = clock();
    float diff = ((float)time_after - (float)time_before);
    cout << "Repeated Squaring Result: " << endl;
    cout << setw(5) << setprecision(1) << B << endl;
    cout << "It took " << diff / CLOCKS_PER_SEC << " seconds to complete." << endl;
    time_before = clock();
    Matrix C = direct_multiplication(A, exponent);
    time_after = clock();
    diff = ((float)time_after - (float)time_before);
    cout << "Direct Multiplication Result: " << endl;
    //cout << setw(5) << setprecision(1) << C << endl;
    cout << "It took " << diff / CLOCKS_PER_SEC << " seconds to complete." << endl;
    
    create_file(A, no_row);
    return 0;
}

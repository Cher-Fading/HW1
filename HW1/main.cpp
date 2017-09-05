//
//  main.cpp
//  HW1
//
//  Created by Xiaoning Wang on 9/4/17.
//  Copyright © 2017 Cher Wang. All rights reserved.
//
//

#include <iostream>
#include <fstream>
#include <string>
#include "nr3.h"
#include "ludcmp.h"
#include "svd.h"
#include "gaussj.h"

int main(int argc, const char * argv[]) {
    MatDoub test;
    MatDoub A;
    VecDoub x;
    VecDoub b;
    VecDoub c;
    VecDoub x1;
    VecDoub x2;
    VecDoub r1;
    VecDoub r2;
    VecDoub x_numerical;
    VecDoub r_numerical;
    VecDoub error_numerical;
    
    
    //Homework 1 Question 1
    // Initialize vectors and matrix A
    Doub a[2]={0.1869169, 0.2155734};
    b = VecDoub(2,a);
    Doub a1[2]={1.0,-1.0};
    x = VecDoub(2,a1);
    Doub a2[2]={0.9999999,-1.0000001};
    x1 = VecDoub(2,a2);
    Doub a3[2]={0.4073666,-0.1945277};
    x2 = VecDoub(2,a3);
    A =MatDoub(2,2);
    A[0][0]=0.7073725;
    A[0][1]=0.5204556;
    A[1][0]=0.8158208;
    A[1][1]=0.6002474;
    
    //a)Compute residuals
    //Set up 2 vectors to contain residual vectors
    r1 = VecDoub(2);
    r2 = VecDoub(2);
    
    //Do matrix multiplication and vector subtraction
    for (int i = 0; i<A.nrows(); i++){
        for (int j = 0; j<A.ncols(); j++){
            r1[i]+=A[i][j]*x1[j];
            r2[i]+=A[i][j]*x2[j];
        }
        r1[i]=b[i]-r1[i];
        r2[i]=b[i]-r2[i];
    }
    //Print out result
    std::cout<<"1. a)"<<std::endl;
    std::cout<<"Residual of x_1 is ("<<r1[0]<<", "<<r1[1]<<")."<<std::endl;
    std::cout<<"Residual of x_2 is ("<<r2[0]<<", "<<r2[1]<<")."<<std::endl;
    Doub sr1 = sqrt(r1[0]*r1[0]+r1[1]*r1[1]);
    Doub sr2 = sqrt(r2[0]*r2[0]+r2[1]*r2[1]);
    std::cout<<"Residual magnitude of x_1 is r_1 = "<<sr1<<std::endl;
    std::cout<<"Residual magnitude of x_2 is r_2 = "<<sr2<<std::endl;;
    if (sr1 < sr2){
        std::cout << "The more accurate answer (x_1) does have smaller residual\n";
    }
    else {
        std::cout << "The more accurate answer (x_1) does NOT have smaller residual\n";
    }
    
    
    //b)Compute residuals and error of using LU decomposition
    //Use NR routine of LU decomposition to construct a LUdcmp struct, and then use its method solve to solve for the given vector b.
    std::cout<<std::endl<<"1. b)"<<std::endl;
    LUdcmp LUA = LUdcmp(A);
    x_numerical = VecDoub(2);
    LUA.solve(b, x_numerical);
    
    std::cout<<"x_numerical using LU decomposition routine is ("<<x_numerical[0]<<", "<<x_numerical[1]<<"). "<<std::endl;
    error_numerical  = VecDoub(2);
    for (int i = 0; i<error_numerical.size(); i++){
        error_numerical[i]=x_numerical[i]-x[i];
    }
    std::cout<<"x_numerical error using LU decomposition routine is ("<<error_numerical[0]<<", "<<error_numerical[1]<<"). "<<std::endl;
    
    r_numerical = VecDoub(2);
    for (int i = 0; i<A.nrows(); i++){
        for (int j = 0; j<A.ncols(); j++){
            r_numerical[i]+=A[i][j]*x_numerical[j];
        }
        r_numerical[i]=b[i]-r_numerical[i];
    }
    std::cout<<"x_numerical residual using LU decomposition routine is ("<<r_numerical[0]<<", "<<r_numerical[1]<<"). "<<std::endl;
    Doub sr_numerical = sqrt(r_numerical[0]*r_numerical[0]+r_numerical[1]*r_numerical[1]);
    Doub serror_numerical = sqrt(error_numerical[0]*error_numerical[0]+error_numerical[1]*error_numerical[1]);
    Int order = log10(sr_numerical/serror_numerical);
    std::cout<<"The solution residual is about an order of 10 to the exponent of "<<order<<" in comparison to the solution error."<<std::endl;
    
    //c) Compute condition number
    std::cout<<std::endl<<"1. c)"<<std::endl;
    SVD SVDA = SVD(A);
    Doub condA_inv = SVDA.inv_condition();
    std::cout<<"The condition number of A is "<<1./condA_inv<<std::endl;
    std::cout<<"It's reciprocal is "<<condA_inv<<std::endl;
    
    //2. Firstly load the matrix and vectors
    //To avoid memory leak, declare all variables for question 2 as new variables, to distinguish, add "_2" at the end of variable name.
    
    ifstream in;
    ifstream in2;
    ifstream in3;
    MatDoub A_2;
    VecDoub b_2;
    VecDoub c_2;
    VecDoub x_2;
    VecDoub x_2LU;
    VecDoub x_2imp;
    VecDoub r_2imp;
    Doub sr_2;
    Doub sr_2LU;
    Doub sr_2imp;
    
    A_2 = MatDoub(1600,1600);
    in.open("/Users/cherwang/Google Drive/2017 Fall/PHYS 4480/HW1/HW1/A.txt");
    if (!in){
        std::cout<<"error in opening A.txt"<<std::endl;
    }
    else{
        for (int i = 0; i < 1600; i++) {
            for (int j = 0; j < 1600; j++) {
                in >> A_2[i][j];
            }
        }
    }
    in.close();
    
    b_2 = VecDoub(1600);
    in2.open("/Users/cherwang/Google Drive/2017 Fall/PHYS 4480/HW1/HW1/b.txt");
    if (!in2){
        std::cout<<"error in opening b.txt"<<std::endl;
    }
    else{
        for (int i = 0; i < 1600; i++){
            in2 >> b_2[i];
        }
    }
    in2.close();
    
    c = VecDoub(1600);
    in3.open("/Users/cherwang/Google Drive/2017 Fall/PHYS 4480/HW1/HW1/c.txt");
    if (!in3){
        std::cout<<"error in opening c.txt"<<std::endl;
    }
    else{
        for (int i = 0; i < 1600; i++){
            in3 >> c[i];
        }
    }
    
    //a) Make a copy of loaded matrix A_2, and multiply it by b to get x
    std::cout<<std::endl<<"2. a)"<<std::endl;
    test = A_2;
    gaussj(test);
    x_2 = VecDoub(1600);
    sr_2 = 0;
    for (int i = 0; i < test.nrows(); i++){
        for (int j = 0; j < test.ncols(); j++){
            x_2[i]+=test[i][j]*b_2[j];
        }
        sr_2 += abs(x_2[i]-c[i]);
    }
    std::cout<< "The measure of error using matrix inverse is "<< sr_2 << std::endl;
    
    //b) Use LU decomposition to compute a solution and measure its error
    std::cout<<std::endl<<"2. b)"<<std::endl;
    LUdcmp LUA_2 = LUdcmp(A_2);
    x_2LU = VecDoub(1600);
    sr_2LU = 0;
    LUA_2.solve(b_2, x_2LU);
    for (int i = 0; i < x_2LU.size(); i++){
        sr_2LU += abs(x_2LU[i]-c[i]);
    }
    std::cout<< "The measure of error using LU decomposition is "<< sr_2LU << std::endl;
    
    //c) Use interactive improvement routine from NR codes and report on residual
    //Use a separate variable to store improved solution
    x_2imp = x_2LU;
    LUA_2.mprove(b_2, x_2imp);
    
    //initialize r_2imp
    r_2imp = VecDoub(16);
    sr_2imp = 0;
    std::cout<<"The residual of using improved LU decomposition is (";
    for (int i = 0; i < A_2.nrows(); i++){
        for (int j = 0; j < A_2.ncols(); j++){
            r_2imp[i]+=A_2[i][j]*x_2imp[j];
        }
        r_2imp[i] = b_2[i]-r_2imp[i];
        std::cout<<r_2imp[i];
        if ( i != A_2.nrows()-1) {
            std::cout<<", ";
        }
        sr_2imp += abs(x_2imp[i]-c[i]);
    }
    std::cout<<")."<<std::endl;
    std::cout<<"The measure of error using improved LU decomposition is "<< sr_2imp <<std::endl;
    
    SVD SVD_2 = SVD(A_2);
    Doub condA_2_inv = SVD_2.inv_condition();
    std::cout<<"The condition number of matrix A is" << 1./condA_2_inv<<std::endl;
    return 0;
}


//
//  main.cpp
//  HW1
//
//  Created by Cher Wang on 9/4/17.
//  Copyright © 2017 Cher Wang. All rights reserved.
//

#include <iostream>
#include "nr3.h"
#include "ludcmp.h"

int main(int argc, const char * argv[]) {
    //Homework 1 Question 1
    // Initialize vectors and matrix A
    Doub a[2]={0.1869169, 0.2155734};
    VecDoub b(2,a);
    Doub a1[2]={1.0,-1.0};
    VecDoub x(2,a1);
    Doub a2[2]={0.9999999,-1.0000001};
    VecDoub x1(2,a2);
    Doub a3[2]={0.4073666,-0.1945277};
    VecDoub x2(2,a3);
    MatDoub A(2,2);
    A[0][0]=0.7073725;
    A[0][1]=0.5204556;
    A[1][0]=0.8158208;
    A[1][1]=0.6002474;
    
    //a)Compute residuals
    //Set up 2 vectors to contain residual vectors
    VecDoub r1(2);
    VecDoub r2(2);
    
    //Do matrix multiplication and vector subtraction
    for (int i = 0; i<=1; i++){
        for (int j = 0; j<=1; j++){
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
    LUdcmp LUA = LUdcmp(A);
    VecDoub x_numerical(2);
    LUA.solve(b, x_numerical);
    
    return 0;
}

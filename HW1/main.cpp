//
//  main.cpp
//  HW1
//
//  Created by Cher Wang on 9/4/17.
//  Copyright © 2017 Cher Wang. All rights reserved.
//

#include <iostream>
#include "nr3.h"

int main(int argc, const char * argv[]) {
    // Initialize vector
    Doub a[2]={0.1869169, 0.2155734};
    VecDoub b(2,a);
    std::cout<<b[0]<<"b[0]"<<std::endl;
    Doub a1[2]={1.0,-1.0};
    VecDoub x(2,a1);
    Doub a2[2]={0.9999999,1.0000001};
    VecDoub x1(2,a2);
    Doub a3[2]={0.4073666,0.1945277};
    VecDoub x2(2,a3);
    
    std::cout<<x[0]<<"x[0]"<<std::endl;
    std::cout<<x1[0]<<"x1[0]"<<std::endl;
    std::cout<<x2[0]<<"x2[0]"<<std::endl;
    std::cout << "Hello, World!\n";
    return 0;
}

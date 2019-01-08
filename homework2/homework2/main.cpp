//
//  main.cpp
//  homework2
//
//  Created by jzhou on 9/13/18.
//  Copyright © 2018 jzhou. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "sudoku.h"

using namespace std;

int main (int argc, char * const argv[]) {
    Sudoku x;
    x.read_puzzle(argc, argv);
    x.print_puzzle();
    x.Solve(0,0);
    x.alternate_Solve(0, 0);
    //x.print_puzzle();
    return 0;
}

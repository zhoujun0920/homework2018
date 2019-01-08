//
//  main.cpp
//  homework1
//
//  Created by Jun Zhou on 8/28/18.
//  Copyright © 2018 jzhou. All rights reserved.
//

#include <iostream>
using namespace std;

int f(int n, int &m) {
    if (n < 10) {
        return m + 1;
    } else {
        m++;
        return (f(n / 10, m));
    }
}

void print_integer(int num) {
    if (num / 10) {
        print_integer(num / 10);
    }
    putchar(num % 10 + '0');
}

int main(int argc, const char * argv[]) {
    // insert code here...
    //print_integer(12345);
    cout << endl;
    return 0;
}

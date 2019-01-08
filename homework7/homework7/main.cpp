//
//  main.cpp
//  homework7
//
//  Created by Jun Zhou on 11/1/18.
//  Copyright © 2018 jzhou. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <random>
#include "hint.h"
using namespace std;

int main (int argc, char* argv[])
{
    I_have_nothing_apropos_for_this_class x;
    double alice_success_prob, bob_success_prob;
    sscanf (argv[1], "%lf", &alice_success_prob);
    sscanf (argv[2], "%lf", &bob_success_prob);

    cout << "Probability of success for Alice = " << alice_success_prob << endl;
    cout << "Probability of success for Bob = " << bob_success_prob << endl;

    x.set_probability(alice_success_prob, bob_success_prob);

    int optimal = x.search_result();
    if (optimal > 0)
        cout << "The optimal number of coin tosses in each game is " << optimal << endl;
    else {
        cout << "The optimal number of coin tosses in each game exceeds 100... Quitting" << endl;
    }
    x.get_prob_values(30);
}

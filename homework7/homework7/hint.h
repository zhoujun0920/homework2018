//
//  main.h
//  homework7
//
//  Created by Jun Zhou on 11/1/18.
//  Copyright © 2018 jzhou. All rights reserved.
//

#ifndef main_h
#define main_h

#include <cmath>
#include <random>
#include <iostream>
using namespace std;

class I_have_nothing_apropos_for_this_class
{
private:
    double alice_probability, bob_probability;
    default_random_engine generator;
    uniform_real_distribution<double> distribution;
    
    // private member function: uniform RV generator
    double get_uniform()
    {
        // write the appropriate code here
        return distribution(generator);
    }
    
    // private member function: nCi (i.e. n-take-i)
    int take(int n, int i)
    {
        // write a **RECURSIVE** implementation of n-take-i.
        // If you made it non-recurisive (i.e. n!/((n-i)!i!)) -- it
        // will take too long for large sizes
        if (i == 0 || n == i)
            return 1;
        else if (n > i)
            return take(n - 1, i) + take(n - 1, i - 1);
        return 0;
    }
    
    // this routine implements the probability that Alice has more
    // heads than Bob after n-many coin tosses
    double theoretical_value(double q, double p, int n)
    {
        // implement equation 1.1 of Addona-Wagon-Wilf paper
        double sum_prob = 0;
        for (int r = 0; r <= n; r++) {
            double sum_prob_bob_loss = take(n, r) * pow(bob_probability, r) * pow(1.0 - bob_probability, n - r);
            double sum_prob_alice_win = 0;
            for (int s = r + 1; s <= n; s++) {
                sum_prob_alice_win += take(n, s) * pow(alice_probability, s) * pow(1.0 - alice_probability, n - s);
            }
            sum_prob += sum_prob_bob_loss * sum_prob_alice_win;
        }
        return sum_prob;
    }
    
public:
    // public function:
    void set_probability(double alice_p, double bob_p)
    {
        alice_probability = alice_p;
        bob_probability = bob_p;
    }
    
    // probability of Alice winning the game.
    double simulated_value(int number_of_coin_tosses_in_each_game, int no_of_trials)
    {
        int no_of_wins_for_alice = 0;
        for (int i = 0; i < no_of_trials; i++)
        {
            int number_of_heads_for_alice = 0;
            int number_of_heads_for_bob = 0;
            for (int j = 0; j < number_of_coin_tosses_in_each_game; j++)
            {
                if (get_uniform() < alice_probability)
                    number_of_heads_for_alice++;
                if (get_uniform() < bob_probability)
                    number_of_heads_for_bob++;
            }
            if (number_of_heads_for_alice > number_of_heads_for_bob)
                no_of_wins_for_alice++;
        }
        return (((double) no_of_wins_for_alice)/((double) no_of_trials));
    }
    
    int search_result()
    {
        // implememt a discrete-search procedure for the optimal n-value.
        // start with n = 1 and find the discrete-value of n that has
        // the largest probability for Alice winning.  Why would this work?
        // See Theorem 2.2 of the paper for the reason!
        double largest_p = 0;
        int largest_n = 1;
        while (true) {
            double temp_p = theoretical_value(alice_probability, bob_probability, largest_n);
            if (temp_p >= largest_p) {
                largest_p = temp_p;
                largest_n++;
            } else {
                largest_n--;
                break;
            }
        }
        return largest_n;
    }
    
    void get_prob_values(int n) {
        cout << "Simulation" << ", " << "Theoretical" << endl;
        for (int i = 1; i <= n; i++) {
            cout << simulated_value(i, 500000) << ", " << theoretical_value(alice_probability, bob_probability, i) << endl;
        }
    }
};

#endif /* main_h */

//
//  main.cpp
//  homework10_2
//
//  Created by Jun Zhou on 12/2/18.
//  Copyright © 2018 jzhou. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <map>
using namespace std;

float up_factor, uptick_prob, downtick_prob, risk_free_rate, strike_price;
float initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;
map<pair<int, int>, float> lookup_table;

float max(float a, float b) {
    return (b < a )? a:b;
}

float american_call_option(int k, int i, float current_stock_price) {
    if (k == no_of_divisions)
        return max(0.0, (current_stock_price - strike_price));
    else {
        if (lookup_table.find(make_pair(k, i)) == lookup_table.end()) {
            float temp = max((current_stock_price - strike_price),
                             (uptick_prob * american_call_option(k + 1, i + 1, current_stock_price * up_factor) +
                              (1 - uptick_prob - downtick_prob) * american_call_option(k + 1, i, current_stock_price) +
                              downtick_prob * american_call_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
            lookup_table[make_pair(k, i)] = temp;
            return temp;
        }
        return lookup_table[make_pair(k, i)];
    }
}

float american_put_option(int k, int i, float current_stock_price) {
    if (k == no_of_divisions)
        return max(0.0, (strike_price - current_stock_price));
    else {
        if (lookup_table.find(make_pair(k, i)) == lookup_table.end()) {
            float temp = max((strike_price - current_stock_price),
                             (uptick_prob * american_put_option(k + 1, i + 1, current_stock_price * up_factor) +
                              (1 - uptick_prob - downtick_prob) * american_put_option(k + 1, i, current_stock_price) +
                              downtick_prob * american_put_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
            lookup_table[make_pair(k, i)] = temp;
            return temp;
        }
        return lookup_table[make_pair(k, i)];
    }
}

int main (int argc, char* argv[])
{
    sscanf (argv[1], "%f", &expiration_time);
    sscanf (argv[2], "%d", &no_of_divisions);
    sscanf (argv[3], "%f", &risk_free_rate);
    sscanf (argv[4], "%f", &volatility);
    sscanf (argv[5], "%f", &initial_stock_price);
    sscanf (argv[6], "%f", &strike_price);
    up_factor = exp(volatility * sqrt(2 * expiration_time / ((float) no_of_divisions)));
    R = exp(risk_free_rate * expiration_time / ((float) no_of_divisions));
    //uptick_prob = (R - (1 / up_factor)) / (up_factor - (1 / up_factor));
    uptick_prob = pow((sqrt(R) - 1.0 / sqrt(up_factor)) / (sqrt(up_factor) - 1.0 / sqrt(up_factor)), 2.0);
    downtick_prob = pow((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - 1.0 / sqrt(up_factor)), 2.0);
    cout << "(Memoized) Recursive Trinomial American-Asian Option Pricing" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "--------------------------------------" << endl;
    cout << "R = " << R << endl;
    cout << "Up Factor = " << up_factor << endl;
    cout << "Uptick Probability = " << uptick_prob << endl;
    cout << "Downtick Probability = " << downtick_prob << endl;
    cout << "Notick Probability = " << 1 - uptick_prob - downtick_prob << endl;
    cout << "--------------------------------------" << endl;
    double call_price = american_call_option(0, 0, initial_stock_price);
    cout << "Trinomial Price of an American Call Option = " << call_price << endl;
    cout << "--------------------------------------" << endl;
    lookup_table.clear();
    double put_price = american_put_option(0, 0, initial_stock_price);
    cout << "Trinomial Price of an American Put Option = " << put_price << endl;
    cout << "--------------------------------------" << endl;
    cout << "Verifying Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
    cout <<  initial_stock_price << " + " << put_price << " - " << call_price;
    cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
    cout << initial_stock_price + put_price - call_price << " = " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
    cout << "--------------------------------------" << endl;
}

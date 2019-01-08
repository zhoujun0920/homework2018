//
//  main.cpp
//  final_2
//
//  Created by Jun Zhou on 12/6/18.
//  Copyright © 2018 jzhou. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "normdist.h"
#define UNIT_STEP(x) ((x)>0?(1):(0))
using namespace std;

double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility, barrier_price;
int no_of_trials, no_of_discrete_barrier;

double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a=fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
};

double get_p(double final_stock_price) {
    if (initial_stock_price <= barrier_price || final_stock_price <= barrier_price) {
        return 1;
    }
    return exp(-(2.0 * log(initial_stock_price / barrier_price) *
                 log(final_stock_price / barrier_price)) / (pow(volatility, 2.0) * expiration_time));
}

double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}

int main (int argc, char* argv[])
{
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%lf", &risk_free_rate);
    sscanf (argv[3], "%lf", &volatility);
    sscanf (argv[4], "%lf", &initial_stock_price);
    sscanf (argv[5], "%lf", &strike_price);
    sscanf (argv[6], "%d", &no_of_trials);
    sscanf (argv[7], "%d", &no_of_discrete_barrier);
    sscanf (argv[8], "%lf", &barrier_price);
    
    // delta T
    double delta_expiration_time = expiration_time / no_of_discrete_barrier;
    
    // R and SD for delta T
    double R = (risk_free_rate - 0.5 * pow(volatility, 2)) * delta_expiration_time;
    double SD = volatility * sqrt(delta_expiration_time);
    
    cout << "--------------------------------" << endl;
    cout << "European Down-and-Out Discrete Barrier Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Duration T (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate r = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price S0 = " << initial_stock_price << endl;
    cout << "Strike Price K = " << strike_price << endl;
    cout << "Number of Trials n = " << no_of_trials << endl;
    cout << "Number of discrete barriers m = " << no_of_discrete_barrier << endl;
    cout << "Barrier Price B = " << barrier_price << endl;
    cout << "--------------------------------" << endl;
    
    double explicit_call_option_price = 0.0;
    double explicit_put_option_price = 0.0;
    
    double brownian_bride_correction_call_option_price = 0.0;
    double brownian_bride_correction_put_option_price = 0.0;
    
    for (int i = 0; i < no_of_trials; i++) {
        // generate unit-normals using Box-Muller Transform
        double x = get_uniform();
        double y = get_uniform();
        double a =  sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
        double b =  sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);
        
        // get four path for the price of one path
        double S1 = initial_stock_price * exp(R + SD * a);
        double S2 = initial_stock_price * exp(R - SD * a);
        double S3 = initial_stock_price * exp(R + SD * b);
        double S4 = initial_stock_price * exp(R - SD * b);
        
        // generate m-many price points, and knock out the price lower then barrier
        for (int j = 1; j < no_of_discrete_barrier; j++) {
            x = get_uniform();
            y = get_uniform();
            a =  sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
            b =  sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);
            
            double temp1 = S1 * exp(R + SD * a);
            double temp2 = S2 * exp(R - SD * a);
            double temp3 = S3 * exp(R + SD * b);
            double temp4 = S4 * exp(R - SD * b);
            
            S1 = temp1 <= barrier_price ? 0 : temp1;
            S2 = temp2 <= barrier_price ? 0 : temp2;
            S3 = temp3 <= barrier_price ? 0 : temp3;
            S4 = temp4 <= barrier_price ? 0 : temp4;
        }
        
        // calculatet the barrier option price summation
        explicit_call_option_price += (max(0.0, S1 - strike_price) +
                                       max(0.0, S2 - strike_price) +
                                       max(0.0, S3 - strike_price) +
                                       max(0.0, S4 - strike_price)) / 4.0;
        explicit_put_option_price += (max(0.0, strike_price - S1) +
                                      max(0.0, strike_price - S2) +
                                      max(0.0, strike_price - S3) +
                                      max(0.0, strike_price - S4)) / 4.0;
    }
    
    // get the R and SD for T
    R = (risk_free_rate - 0.5 * pow(volatility, 2)) * expiration_time;
    SD = volatility * sqrt(expiration_time);
    
    for (int i = 0; i < no_of_trials; i++) {
        // generate unit-normals using Box-Muller Transform
        double x = get_uniform();
        double y = get_uniform();
        double a =  sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
        double b =  sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);
        
        // get the ST
        double St1 = initial_stock_price * exp(R + SD * a);
        double St2 = initial_stock_price * exp(R - SD * a);
        double St3 = initial_stock_price * exp(R + SD * b);
        double St4 = initial_stock_price * exp(R - SD * b);
        
        double pd1 = 0.0;
        double pd2 = 0.0;
        double pd3 = 0.0;
        double pd4 = 0.0;
        
        // create the brownian bridge based on the equation
        for (int i = 1; i < no_of_discrete_barrier; i++) {
            double u1 = initial_stock_price + (delta_expiration_time * i - delta_expiration_time) /
            (expiration_time - delta_expiration_time) * (St1 - initial_stock_price);
            double u2 = initial_stock_price + (delta_expiration_time * i - delta_expiration_time) /
            (expiration_time - delta_expiration_time) * (St2 - initial_stock_price);
            double u3 = initial_stock_price + (delta_expiration_time * i - delta_expiration_time) /
            (expiration_time - delta_expiration_time) * (St3 - initial_stock_price);
            double u4 = initial_stock_price + (delta_expiration_time * i - delta_expiration_time) /
            (expiration_time - delta_expiration_time) * (St4 - initial_stock_price);
            double std = sqrt((delta_expiration_time * i - delta_expiration_time) *
                              (expiration_time - delta_expiration_time * i) / (expiration_time - delta_expiration_time));
            pd1 *= (1- N((barrier_price - u1) / std));
            pd2 *= (1- N((barrier_price - u2) / std));
            pd3 *= (1- N((barrier_price - u3) / std));
            pd4 *= (1- N((barrier_price - u4) / std));
        }
        
        
        // get the adjusted barrier option price summation
        brownian_bride_correction_call_option_price += (max(0.0, St1 - strike_price) * (1 - pd1) +
                                                        max(0.0, St2 - strike_price) * (1 - pd2) +
                                                        max(0.0, St3 - strike_price) * (1 - pd3) +
                                                        max(0.0, St4 - strike_price) * (1 - pd4)) / 4.0;
        brownian_bride_correction_put_option_price += (max(0.0, strike_price - St1) * (1 - pd1) +
                                                       max(0.0, strike_price - St2) * (1 - pd2) +
                                                       max(0.0, strike_price - St3) * (1 - pd3) +
                                                       max(0.0, strike_price - St4) * (1 - pd4)) / 4.0;
    }
    
    // get the discount price
    explicit_call_option_price = exp(-risk_free_rate * expiration_time) * (explicit_call_option_price / ((double) no_of_trials));
    explicit_put_option_price = exp(-risk_free_rate * expiration_time) * (explicit_put_option_price / ((double) no_of_trials));
    
    brownian_bride_correction_call_option_price = exp(-risk_free_rate * expiration_time) * (brownian_bride_correction_call_option_price / ((double) no_of_trials));
    brownian_bride_correction_put_option_price = exp(-risk_free_rate * expiration_time) * (brownian_bride_correction_put_option_price / ((double) no_of_trials));
    
    cout << "The average Call Price by explicit simulation = " << explicit_call_option_price << endl;
    cout << "The call price with Brownian-Bridge correction on the fianl price = " << brownian_bride_correction_call_option_price << endl;
    cout << "The average Put Price is " << explicit_put_option_price << endl;
    cout << "The put price with Brownian-Bridge correction on the fianl price = " << brownian_bride_correction_put_option_price << endl;
    cout << "--------------------------------" << endl;
}

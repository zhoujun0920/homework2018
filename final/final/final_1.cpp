//
//  main.cpp
//  final
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
int no_of_trials, no_of_sample_points;


double option_price_put_black_scholes(const double& S,      // spot price
                                      const double& K,      // Strike (exercise) price,
                                      const double& r,      // interest rate
                                      const double& sigma,  // volatility
                                      const double& time){
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time) {  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};

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

double option_price_delta_call_black_scholes(const double& S,     // spot price
                                             const double& K,     // Strike (exercise) price,
                                             const double& r,     // interest rate
                                             const double& sigma, // volatility
                                             const double& time){  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double delta = N(d1);
    return delta;
};

double option_price_delta_put_black_scholes(const double& S, // spot price
                                            const double& K, // Strike (exercise) price,
                                            const double& r,  // interest rate
                                            const double& sigma,
                                            const double& time) {
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double delta = -N(-d1);
    return delta;
}

// lession 7
float closed_form_down_and_out_european_call_option()
{
    // I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
    float K = (2*risk_free_rate)/(volatility*volatility);
    float A = option_price_call_black_scholes(initial_stock_price, strike_price,
                                              risk_free_rate, volatility, expiration_time);
    float B = (barrier_price*barrier_price)/initial_stock_price;
    float C = pow(initial_stock_price/barrier_price, -(K-1));
    float D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
    return (A - D*C);
}

// lession 7
float closed_form_down_and_in_european_put_option()
{
    // just making it easier by renaming the global variables locally
    float S = initial_stock_price;
    float r = risk_free_rate;
    float T = expiration_time;
    float sigma = volatility;
    float H = barrier_price;
    float X = strike_price;
    
    // Took these formulae from some online reference
    float lambda = (r+((sigma*sigma)/2))/(sigma*sigma);
    float temp = 2*lambda - 2.0;
    float x1 = (log(S/H)/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    float y = (log(H*H/(S*X))/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    float y1 = (log(H/S)/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    return (-S*N(-x1) + X*exp(-r*T)*N(-x1 + sigma*sqrt(T)) +
            S*pow(H/S, 2*lambda)*(N(y)-N(y1)) -
            X*exp(-r*T)*pow(H/S, temp)*(N(y-sigma*sqrt(T))-N(y1-sigma*sqrt(T))));
}

// lession 7
float closed_form_down_and_out_european_put_option()
{
    float vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price,
                                                       risk_free_rate, volatility, expiration_time);
    float put_down_in = closed_form_down_and_in_european_put_option();
    return (vanilla_put - put_down_in);
}

double max(double a, double b) {
    return (b < a )? a:b;
}

//Calculate the pc based on the formula
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
    sscanf (argv[7], "%d", &no_of_sample_points);
    sscanf (argv[8], "%lf", &barrier_price);
    
    
    //get delta T
    double delta_expiration_time = expiration_time / no_of_sample_points;
    
    //get R and SD for delta T
    double R = (risk_free_rate - 0.5 * pow(volatility, 2)) * delta_expiration_time;
    double SD = volatility * sqrt(delta_expiration_time);
    
    cout << "--------------------------------" << endl;
    cout << "European Down-and-Out Continuous Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Duration T (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate r = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price S0 = " << initial_stock_price << endl;
    cout << "Strike Price K = " << strike_price << endl;
    cout << "Number of Trials n = " << no_of_trials << endl;
    cout << "Number of sample-points m = " << no_of_sample_points << endl;
    cout << "Barrier Price B = " << barrier_price << endl;
    cout << "--------------------------------" << endl;
    
    double closed_form_down_and_out_european_call_option_price = 0.0;
    double closed_form_down_and_out_european_put_optionprice = 0.0;
    
    double p_adjustment_call_option_price = 0.0;
    double p_adjustment_put_option_price = 0.0;
    
    for (int i = 0; i < no_of_trials; i++) {
        // generate unit-normals using Box-Muller Transform
        double x = get_uniform();
        double y = get_uniform();
        double a =  sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
        double b =  sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);
        
        double S1 = initial_stock_price * exp(R + SD * a);
        double S2 = initial_stock_price * exp(R - SD * a);
        double S3 = initial_stock_price * exp(R + SD * b);
        double S4 = initial_stock_price * exp(R - SD * b);
        
        // get four paths for the price of on path and generate m many price points
        for (int j = 1; j < no_of_sample_points; j++) {
            x = get_uniform();
            y = get_uniform();
            a =  sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
            b =  sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);
        
            S1 = S1 * exp(R + SD * a);
            S2 = S2 * exp(R - SD * a);
            S3 = S3 * exp(R + SD * b);
            S4 = S4 * exp(R - SD * b);
        }
        
        // get Pc value
        double p1 = get_p(S1);
        double p2 = get_p(S2);
        double p3 = get_p(S3);
        double p4 = get_p(S4);
        
        
        // adjust the option price by Pc based on the equation
        p_adjustment_call_option_price += (max(0.0, S1 - strike_price) * (1 - p1) +
                              max(0.0, S2 - strike_price) * (1 - p2) +
                              max(0.0, S3 - strike_price) * (1 - p3) +
                              max(0.0, S4 - strike_price) * (1 - p4)) / 4.0;
        p_adjustment_put_option_price += (max(0.0, strike_price - S1) * (1 - p1) +
                             max(0.0, strike_price - S2) * (1 - p2) +
                             max(0.0, strike_price - S3) * (1 - p3) +
                             max(0.0, strike_price - S4) * (1 - p4)) / 4.0;
    }
    // get closed form down and out option price
    closed_form_down_and_out_european_call_option_price = closed_form_down_and_out_european_call_option();
    closed_form_down_and_out_european_put_optionprice = closed_form_down_and_out_european_put_option();
    
    // get the adjustment option price
    p_adjustment_call_option_price = exp(-risk_free_rate * expiration_time) * (p_adjustment_call_option_price / ((double) no_of_trials));
    p_adjustment_put_option_price = exp(-risk_free_rate * expiration_time) * (p_adjustment_put_option_price / ((double) no_of_trials));
    
    cout << "The average Call Price by explicit simulation = " << closed_form_down_and_out_european_call_option_price << endl;
    cout << "The call price using the (1 - p)-adjustment term = " << p_adjustment_call_option_price << endl;
    cout << "Call Price according to Black-Scholes = " <<
    option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
                                    volatility, expiration_time) << endl;
    cout << "The average Put Price is " << closed_form_down_and_out_european_put_optionprice << endl;
    cout << "The put price using the (1 - p)-adjustment term = " << p_adjustment_put_option_price << endl;
    cout << "Put Price according to Black-Scholes = " <<
    option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
                                   volatility, expiration_time) << endl;
    cout << "--------------------------------" << endl;
}

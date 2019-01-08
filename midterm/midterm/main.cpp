// Written by Prof. Sreenivas for IE523: Financial Computing

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include "lp_lib.h"

using namespace std;
# define ERROR() { fprintf(stderr, "Error\n"); exit(1); }

const double ERROR = 1e-10;
int number_of_cash_flows;
vector <double> price_list;
vector <int> maturity_list;
vector <double> yield_to_maturity;
vector <double> duration;
vector <double> convexity;
double debt_obligation_amount;
double time_when_debt_is_due;
vector <double> percentage_of_cash_flow_to_meet_debt_obligation;

double origin_function(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that computes f(r) in page 2 of lesson 3 of my notes
    double temp_sum = 0;
    for (int i = 0; i < cash_flow.size(); i++) {
        temp_sum += cash_flow.at(i) * pow(1 + rate, maturity - (i + 1));
    }
    double fx = price * pow(1 + rate, maturity) - temp_sum;
    return fx;
}

double derivative_function(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that computes f'(r) in the bottom of page 2 of lesson 3
    // of my notes
    double temp_sum = 0;
    for (int i = 0; i < cash_flow.size() - 1; i++) {
        temp_sum += cash_flow.at(i) * (maturity - (i + 1)) * pow(1 + rate, maturity - (i + 1) - 1);
    }
    double dev_fx = maturity * price * pow(1 + rate, maturity - 1) - temp_sum;
    return dev_fx;
}

double Newton_Raphson(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that finds the (only) +ve root of f(r) of page 2 of
    // lesson 3 using Newton-Raphson method
    double rate_old = rate;
    double rate_new = 0;
    double d = rate - rate_new;
    do {
        double fx = origin_function(cash_flow, price, maturity, rate_old);
        double dev_fx = derivative_function(cash_flow, price, maturity, rate_old);
        rate_new = rate_old - fx / dev_fx;
        d = abs(rate_old - rate_new);
        rate_old = rate_new;
    } while (d > ERROR);
    return rate_new;
}

double get_duration(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that computes the duration of a cash flow
    // based on the note
    double temp_sum = 0;
    for (int i = 0; i < cash_flow.size(); i++) {
        temp_sum += ((double)(i + 1) * cash_flow.at(i)) / pow(1 + rate, (double)(i + 1));
    }
    return temp_sum / price;
}

double get_convexity(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that computes the convexity of a cash flow
    // based on the note
    double temp_sum = 0;
    for (int i = 0; i < cash_flow.size(); i++) {
        temp_sum += ((double)(i + 1) * (double)(i + 1 + 1) * cash_flow.at(i)) / pow(1 + rate, (double)(i + 1 + 2));
    }
    return temp_sum / price;
}

double present_value_of_debt()
{
    // compute PV of future debt obligation
    // using the average-value-of-the-YTMs
    double temp_sum = 0;
    for (int i = 0; i < yield_to_maturity.size(); i++) {
        temp_sum += yield_to_maturity.at(i);
    }
    double average_ytm = temp_sum / (double) number_of_cash_flows; //average yield to maturity
    double pv = debt_obligation_amount / pow((average_ytm + 1), (double) time_when_debt_is_due); // discount value
    return pv;
}

void print_data(char *filename)
{
    cout << "Input File: " << filename << endl;
    cout << "We owe " << debt_obligation_amount << " in " << time_when_debt_is_due << " years" << endl;
    cout << "Number of Cash Flows: " << number_of_cash_flows << endl;
    for (int i = 0; i < number_of_cash_flows; i++)
    {
        cout << "---------------------------" << endl;
        cout << "Cash Flow #" << i+1 << endl;
        cout << "Price = " << price_list[i] << endl;
        cout << "Maturity = " << maturity_list[i] << endl;
        cout << "Yield to Maturity = " << yield_to_maturity[i] << endl;
        cout << "Duration = " << duration[i] << endl;
        cout << "Convexity = " << convexity[i] << endl;
        cout << "Percentage of Face Value that would meet the obligation = " <<
        percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
    }
    cout << "---------------------------" << endl;
}

void get_data(char* argv[])
{
    // write the code that reads the data from the file identified
    // on the command-line.
    ifstream input_file(argv[1]);
    if (input_file.is_open()) {
        input_file >> number_of_cash_flows; // read numbe of cash flow
        for (int i = 0; i < number_of_cash_flows; i++) {
            double price;
            input_file >> price; //read price
            price_list.push_back(price);
            int maturity;
            input_file >> maturity; //read maturity
            maturity_list.push_back(maturity);
            vector<double> cash_flow;
            for (int j = 0; j < maturity; j++) {
                double ci;
                input_file >> ci;
                cash_flow.push_back(ci);
            } // read cash flow
            double ytm = Newton_Raphson(cash_flow, price, maturity, 0.0); // calculate yield to maturity
            yield_to_maturity.push_back(ytm);
            double d = get_duration(cash_flow, price, maturity, ytm); // calculate duration
            duration.push_back(d);
            double c = get_convexity(cash_flow, price, maturity, ytm); // calculate convexity
            convexity.push_back(c);
        }
        input_file >> debt_obligation_amount; // read debt obligation amount
        input_file >> time_when_debt_is_due; // read due time
        for (int i = 0; i < number_of_cash_flows; i++) {
            double p = present_value_of_debt() / price_list.at(i); // calculate percentage
            percentage_of_cash_flow_to_meet_debt_obligation.push_back(p);
        }
    }
}

void get_optimal_portfolio()
{
    // write the lp_solve specific material that
    // computes the optimal_portfolio
    lprec *lp; // claim lp solver
    double solution[number_of_cash_flows]; // init solution array
    lp = make_lp(0, number_of_cash_flows); // init lp solver
    set_verbose(lp, 3);
    double eq0[number_of_cash_flows + 1]; // init first constraint about price and present value of debt
    eq0[0] = 0;
    for (int i = 1; i <= number_of_cash_flows + 1; i++){
        eq0[i] = price_list[i - 1];
    }
    if (!add_constraint(lp, eq0, EQ, present_value_of_debt())) // add constraint
        ERROR()
    double eq1[number_of_cash_flows + 1]; // init second constraint about duration and due time
    eq1[0] = 0;
    for (int i = 1; i < number_of_cash_flows + 1; i++) {
        eq1[i] = duration.at(i - 1) * (price_list[i - 1] / present_value_of_debt());
    }
    if (!add_constraint(lp, eq1, EQ, time_when_debt_is_due)) // add constraint
        ERROR();
    double obj[number_of_cash_flows + 1]; // init object, multiply -1 to maximize the object function
    obj[0] = 0;
    for (int i = 1; i < number_of_cash_flows + 1; i++) {
        obj[i] = -convexity[i - 1] * (price_list[i - 1] / present_value_of_debt());
    }
    if (!set_obj_fn(lp, obj)) // add object
        ERROR();
    printf("We can show the current problem with print_lp(lp)\n");
    print_lp(lp);
    int result = solve(lp); // solve lp
    if (result == 0) {
        cout << "Largest Convexity we can get is: " << -get_objective(lp) << endl; // print optimal values
        get_variables(lp, solution); // print optimizing value
        cout << "Optimal portfolio:" << endl;
        for (int i = 0; i < number_of_cash_flows; i++) {
            cout << "%Cash Flow:" << i+1 << "  ";
            cout << solution[i] << endl;
        }
        cout << endl << "That is, buy" << endl;
        for (int i = 0; i < number_of_cash_flows; i++)
        if (solution[i] !=0) {
            cout << "$" << price_list[i] * solution[i];
            cout << " of Cash Flow#" << i + 1 << endl;
        }
    } else {
        cout << "There is no portfolio that meets the duration constraint of ";
        cout << time_when_debt_is_due << " years" << endl;
    }
}

int main (int argc, char* argv[])
{
    if (argc == 1) {
        cout << "Input filename missing" << endl;
    }
    else
    {
        get_data(argv);
        
        print_data(argv[1]);
        
        get_optimal_portfolio();
    }
    return (0);
}

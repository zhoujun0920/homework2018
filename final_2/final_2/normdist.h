//
//  normdist.h
//  final_2
//
//  Created by Jun Zhou on 12/6/18.
//  Copyright © 2018 jzhou. All rights reserved.
//

#ifndef _NORMAL_DIST_H_
#define _NORMAL_DIST_H_

double n(const double& z);                          // normal distribution function
double N(const double& z);                          // cumulative probability of univariate normal
double N(const double& a, const double& b, const double& rho);// cumulative probability of bivariate normal
double N3(const double& h, const double& k, const double& j,
          const double& rho12, const double& rho13, const double& rho23); // trivariate
double random_normal(); // random numbers with mean zero and variance one
double random_uniform_0_1();
#endif

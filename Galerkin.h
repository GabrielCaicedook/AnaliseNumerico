//
//  ExactIntegration.hpp
//  AnaliseNumerico
//
//  Created by Gabriel Caicedo on 7/04/24.
//

#include <iostream>
#include <math.h>
#include <vector>

class ExactIntegration{
public:
    
    double DefinedIntegral(double x);
    ExactIntegration();
    
    double ValExact(double x1, double x2);
    
    std::vector<double> ValExact(double a, double b, int ref);
    
};

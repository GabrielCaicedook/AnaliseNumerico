
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;

class NumericalIntegration{
public:
    //0 Trapecios 1Simpsons 2...
    int ftype;
    std::function<double(double)> funcion;
   
    void SetFunction(std::function<double( double)> fun );
    
    NumericalIntegration();
    NumericalIntegration(int type);
    void SetType(int type);
    
    //Method 1
    double IntegrateTrapecio(double x1, double x2);
    std::vector<double> IntegrateTrapecio(double a, double b, int ref);
    
    //Method 2
    //IntegrateSimpson1 = Simpson1/3
    double IntegrateSimpson1(double x1, double x2);
    std::vector<double> IntegrateSimpson1(double a, double b, int ref);
    
    //Method 3
    //IntegrateSimpson2 = Simpson3/8
    double IntegrateSimpson2(double x1, double x2);
    std::vector<double> IntegrateSimpson2(double a, double b, int ref);
    
    //ExactIntegrate
    double DefinedIntegral(double x);
    //ExactIntegration();
    double ValExact(double x1, double x2);
    std::vector<double> ValExact(double a, double b, int ref);
    double IntegrateGauss(double x1, double x2, int np);
    std::vector<double> IntegrateGauss(double a, double b, int ref, int np);
    std::vector<std::pair<double, double> > Calc_points_weights(int num_points);
};

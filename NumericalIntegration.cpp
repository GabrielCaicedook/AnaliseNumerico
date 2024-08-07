
#include "NumericalIntegration.h"
#include <fstream>
#include <vector>
using namespace std;

NumericalIntegration::NumericalIntegration(){
    
}
NumericalIntegration::NumericalIntegration(int type){
    
}
void NumericalIntegration::SetType(int type){
    ftype = type;
}
void NumericalIntegration::SetFunction(std::function<double( double)> fun ){
    funcion = fun;
}
void NumericalIntegration::SetFunctionExact(std::function<double( double)> fun ){
    funcionExact = fun;
}

//EXACT INTEGRATION
//Implementation of definite integral
double NumericalIntegration::ValExact(double x1, double x2){
    double val1 = funcionExact(x1);
    double val2 = funcionExact(x2);
    double area = val2 - val1;
    if(area <0.0){
        std::cout<<"Area Negativa"<<std::endl;
    }
    return area;
    }
//Implementation according to the number of refinements
std::vector<double> NumericalIntegration::ValExact(double a, double b, int ref){
    int nel = pow(2, ref);
    std::vector<double> vectorAreas(nel,0);
    double h = (b-a)/nel;
    for (int i=0; i< nel; i++){
        double point1 = a + i*h;
        double point2 = a + (i+1)*h;
        double area = ValExact(point1, point2);
        vectorAreas[i]=area;
    }
    fsolExact = vectorAreas;
    return vectorAreas;
}

//TRAPEZIUM METHOD
double NumericalIntegration::IntegrateTrapecio(double x1, double x2){
    //x1 and x2 = limits of integration
    double val2 = funcion(x2);
    double val1 = funcion(x1);
    double area = 0.5*(x2-x1)*( val2+val1 );
    return area;
}
//Implementation according to the number of refinements
std::vector<double> NumericalIntegration::IntegrateTrapecio(double a, double b, int ref){
    double areatotal=0;
    int nel = pow(2, ref);
    std::vector<double> vectorAreas(nel,0);
    double h = (b-a)/nel;
    for (int i=0; i< nel; i++) {
        double point1 = a + i*h;
        double point2 = a + (i+1)*h;
        double area = IntegrateTrapecio(point1, point2);
        vectorAreas[i] = area;
        
    }
    fsol = vectorAreas;
    return vectorAreas;
}

// SIMPSON METHOD 1/3
double NumericalIntegration::IntegrateSimpson1(double x1,double x2){
    double area = ((x2-x1)/6)*(funcion(x1)+(4*funcion((x1+x2)/2))+funcion(x2));
    return area;
}

//Implementation according to the number of refinements
std::vector<double> NumericalIntegration::IntegrateSimpson1(double a, double b, int ref){
    
    int nel = pow(2, ref);
    std::vector<double> vectorAreas(nel,0);
    double h = (b-a)/nel;
    for (int i=0; i< nel; i++) {
        double point1 = a + i*h;
        double point2 = a + (i+1)*h;
        double area = IntegrateSimpson1(point1, point2);
        vectorAreas[i]=area;
    }
    fsol =vectorAreas;
    return vectorAreas;
}

//SIMPSON METHOD 3/8
double NumericalIntegration::IntegrateSimpson2(double x1, double x2){
    double area = ((x2-x1)/8)*(funcion(x1)+(3*funcion((2*x1+x2)/3))+(3*funcion((x1+(2*x2))/3))+funcion(x2));
    return area;
}

//Implementation according to the number of refinements
std::vector<double> NumericalIntegration::IntegrateSimpson2(double a, double b, int ref){
    int nel = pow(2, ref);
    std::vector<double> vectorAreas(nel,0);
    double h = (b-a)/nel;
    for (int i=0; i< nel; i++) {
        double point1 = a + i*h;
        double point2 = a + (i+1)*h;
        double area = IntegrateSimpson2(point1, point2);
        vectorAreas[i]=area;
    }
    fsol =vectorAreas;
    return vectorAreas;
}

//GAUSS METHOD
std::vector<std::pair<double, double> > NumericalIntegration::Calc_points_weights(int num_points) {
    // Table of points and weights for Gaussian squaring
    
    
    std::vector<std::pair<double, double> > points_weights;
    // Defining points and weights for Gaussian quadrature with num_points points
    switch (num_points) {
        case 1:
                  points_weights.resize(1);
                  points_weights[0] = std::make_pair(0.0, 2.0);
            
                  break;
              case 2:
                  points_weights.resize(2);
                  points_weights[0] = std::make_pair(-0.5773502691896257, 1.0);
                  points_weights[1] = std::make_pair(0.5773502691896257, 1.0);
                  break;
              case 3:
                  points_weights.resize(3);
                  points_weights[0] = std::make_pair(-0.7745966692414834, 0.5555555555555556);
                  points_weights[1] = std::make_pair(0.0, 0.8888888888888888);
                  points_weights[2] = std::make_pair(0.7745966692414834, 0.5555555555555556);
                  break;
              case 4:
                  points_weights.resize(4);
                  points_weights[0] = std::make_pair(-0.8611363115940526, 0.3478548451374539);
                  points_weights[1] = std::make_pair(-0.3399810435848563, 0.6521451548625461);
                  points_weights[2] = std::make_pair(0.3399810435848563, 0.6521451548625461);
                  points_weights[3] = std::make_pair(0.8611363115940526, 0.3478548451374539);
                  break;
              case 5:
                  points_weights.resize(5);
                  points_weights[0] = std::make_pair(-0.9061798459386640, 0.2369268850561891);
                  points_weights[1] = std::make_pair(-0.5384693101056831, 0.4786286704993665);
                  points_weights[2] = std::make_pair(0.0, 0.5688888888888889);
                  points_weights[3] = std::make_pair(0.5384693101056831, 0.4786286704993665);
                  points_weights[4] = std::make_pair(0.9061798459386640, 0.2369268850561891);
                  break;
              default:
                  std::cerr << "Número de puntos no soportado" << std::endl;
                  return points_weights;
    }

    return points_weights;
}

//Implemetation based on the number of points
double NumericalIntegration::IntegrateGauss(double x1, double x2, int np){
    
    std::vector<std::pair<double, double> > pointsw =Calc_points_weights(np);
    double area = 0;
    for (int i = 0; i< np; i++) {
        double xi =pointsw[i].first;
        double wi =pointsw[i].second;
        double xi_t = 0.5*(xi+1.0)*(x2 - x1) + x1;
        double wi_t = 0.5*(wi)*(x2 - x1);
        area += funcion(xi_t)*wi_t;
    }
    return area;
}

//Implementation according to the number of refinements

//(a, b = Integration limits); (ref = Domain refinements); (np = Number of points)
std::vector<double> NumericalIntegration::IntegrateGauss(double a, double b, int ref, int np){
    
    int nel = pow(2, ref);
    std::vector<double> vectorAreas(nel,0);
    double h = (b-a)/nel;
    for (int i=0; i< nel; i++) {
        double point1 = a + i*h;
        double point2 = a + (i+1)*h;
        double area = IntegrateGauss(point1, point2, np);
        vectorAreas[i]=area;
    }
    fsol =vectorAreas;
    return vectorAreas;
}
double NumericalIntegration::CalcError(){
    int nareas = fsol.size();
    double sum =0.0;
    ferror.resize(nareas);
    for (int ia = 0; ia< nareas; ia++) {
        ferror[ia] = std::abs(fsol[ia] );
        sum += ferror[ia] ;
    }
    return std::abs(sum - fsolExactval);
}

double NumericalIntegration::Hreturn(double a, double b, int ref){
        int nel = pow(2, ref);
        double h = (b-a)/nel;
    return h;
}

double NumericalIntegration::DefinedIntegral(double x){
    //Set the solution of the integral
    double fx = ((x*x*x)/3);
    return fx;
}

double NumericalIntegration::SumAreas(std::vector<double> areas){
    int nval = areas.size();
    double sum = 0.0;
    for (int iar=0; iar < nval; iar++) {
        sum = sum + areas[iar];
    }
    return sum;
}

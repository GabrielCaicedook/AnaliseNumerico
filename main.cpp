#include "NumericalIntegration.h"
#include "iostream"
#include "fstream"
#include <cmath>
using namespace std;
//#include "ExactIntegration.h"
//Incluir la clase "ExactIntegration.h"

int main(){
    
    //Definiendo funciones iniciales
    auto func1 = [](double x) -> double {
        if (0 <= x && x < 1 / M_PI) {
               return 2 * x + 5;
           } else if (1 / M_PI <= x && x < 2 / M_PI) {
               return (-5 * M_PI * M_PI * (x * x - 2 * x - 5) + 10 * M_PI * x + 2 * x) / (1 + 5 * M_PI * M_PI);
           } else if (2 / M_PI <= x && x <= 8 /M_PI) {
               return 2 * sin(2 * x) + (4 + 20 * M_PI * M_PI + 25 * M_PI * M_PI * M_PI) / (M_PI + 5 * M_PI * M_PI * M_PI) - 2 * sin(4 / M_PI);
           } else {
               // En caso de que x no esté en los rangos definidos, se podría retornar un valor especial o lanzar un error.
               // Aquí se retorna 0 por defecto, pero esto se puede ajustar según las necesidades.
               return 0;
           }
        return 0;
    };
    
    auto func2 = [](double x) -> double {
        if (0 <= x && x < 1 / M_PI) {
               return x*x + 5*x;
           } else if (1 / M_PI <= x && x < 2 / M_PI) {
               double pi = M_PI;
               return    (5 + pi * x * (3 * x + 5 * pi * (3 * x + pi * (15 - (-3 + x) * x)))) / (3 * (pi + 5 * pow(pi, 3)));
           } else if (2 / M_PI <= x && x <= 8 /M_PI) {
               const double pi = M_PI; // Definimos pi usando la constante M_PI de cmath

               double A = (5*pi)*(4+(5*pi))*x;
               double B = (1 + (5*pi*pi))*cos(4.0/pi);
               double bigterm1 = 3*pi*(-20 + A + B);
               double firstTerm = -12 + (pi*(25 + (12*x) + bigterm1));
               double C = -(3*pi)*(1 + (5*pi*pi));
               double D = pi*cos(2*x);
               double E = 2*(-2 + (pi*x))*sin(4/pi);
               double secondTerm = C*(D+E);
               double num = firstTerm + secondTerm;
               double den = 3*((pi*pi) + (5*(pi*pi*pi*pi)));
               return num/den;
 
               
           } else {
               // En caso de que x no esté en los rangos definidos, se podría retornar un valor especial o lanzar un error.
               // Aquí se retorna 0 por defecto, pero esto se puede ajustar según las necesidades.
               return 0;
           }
        return 0;
    };
//    auto antiderivadaTramo1 = [](double x) -> double {
//        return x*x + 5*x;
//    };
//
//    auto antiderivadaTramo2 = [](double x) -> double {
//
//        return (25 * M_PI * M_PI * x + x * x + 5 * M_PI * x * x + 5 * M_PI * M_PI * x * x - (5 * M_PI * M_PI * x * x * x) / 3) / (1 + 5 * M_PI * M_PI); // Reemplazar con la antiderivada real calculada
//    };
//
//    auto antiderivadaTramo3 = [](double x) -> double {
//
//        return ((4 + 20 * M_PI * M_PI + 25 * M_PI * M_PI * M_PI) * x) / (M_PI + 5 * M_PI * M_PI * M_PI) - cos(2 * x) - 2 * x * sin(4 / M_PI);
//
//    };
//
//    auto func2 = [&](double x) -> double {
//        if (0 <= x && x < 1 / M_PI) {
//            return antiderivadaTramo1(x);
//        } else if (1 / M_PI <= x && x < 2 / M_PI) {
//            // Se ajusta por la constante de integración, si es necesario
//            return antiderivadaTramo2(x) - antiderivadaTramo2(1 / M_PI);
//        } else if (2 / M_PI <= x && x <= 8 / M_PI) {
//            // Se ajusta por la constante de integración, si es necesario
//            return antiderivadaTramo3(x) - antiderivadaTramo3(2 / M_PI);
//        } else {
//            // Fuera del rango definido, se puede devolver un valor especial, error o continuar la última antiderivada
//            return 0; // Ajustar según sea necesario
//        }
//    };
//    auto fun2222 = [](double x) -> double {
//        return 2. * (1. - x) * x;
//    };
//    auto fun3 = [](double x) -> double {
//        return x*x;
//    };
//    auto fun2 = [](double x) -> double {
//
//        if (0 <= x && x < 1 / M_PI) {
//               return 2 * x + 5;
//           } else if (1 / M_PI <= x && x < 2 / M_PI) {
//               return (-5 * M_PI * M_PI * (x * x - 2 * x - 5) + 10 * M_PI * x + 2 * x) / (1 + 5 * M_PI * M_PI);
//           } else if (2 / M_PI <= x && x <= 8 /M_PI) {
//               return 2 * sin(2 * x) + (4 + 20 * M_PI * M_PI + 25 * M_PI * M_PI * M_PI) / (M_PI + 5 * M_PI * M_PI * M_PI) - 2 * sin(4 / M_PI);
//           } else {
//               // En caso de que x no esté en los rangos definidos, se podría retornar un valor especial o lanzar un error.
//               // Aquí se retorna 0 por defecto, pero esto se puede ajustar según las necesidades.
//               return 0;
//           }
//        return 0;
//    };
//
//    auto fun5 = [](double x) -> double {
//
//        return 0.0886227 * std::erf(10.0 * (-1.0 + x));
//    };
    
    //Creando archivo de .txt
    ofstream archivo1("Soltrapeciofunc3.txt");
    ofstream archivo2("SolExactfun3.txt");
//    ofstream archivo("Soltrapeciofunc3.txt");
    
    int np= 20;
    double a = 0;
    double b= 8/M_PI;
    double h = (b-a)/np;
    for (int ip=0; ip< np; ip++) {
        double xval = a + ip*h;
        double val = func2(xval);
        if(xval < 2.0/M_PI){
            val = -1.0;
        }
        archivo1 << xval <<" " <<val<< std::endl;

    }
    
    double test1 = cos(4);
    //Test Numerical Integration

    NumericalIntegration solve1;
    solve1.SetFunction(func1);
    for (int ip=2; ip<9; ip++) {
        auto result1 = solve1.IntegrateTrapecio(0,(8/M_PI),ip);
        ip=ip+1;
        for (int jp=0; jp<result1.size(); jp++) {
                    archivo1 << result1[jp] <<" ";
        }
        archivo1<<endl;
    }
         //Cerrar el archivo
    archivo1.close();
    
    NumericalIntegration solve2;
    solve2.SetFunction(func2);
    for (int ip=2; ip<9; ip++) {
        auto result2 = solve2.ValExact(0,(8/M_PI),ip);
        ip=ip+1;
        for (int jp=0; jp<result2.size(); jp++) {
                archivo2 << result2[jp] <<" ";
        }
        archivo2<<endl;
    }
         //Cerrar el archivo
    archivo2.close();
    
    
    
    
//    NumericalIntegration ejemplo1;
//    ejemplo1.SetFunction(fun1);
//    for (int ip=1; ip<5; ip++) {
//        auto result1 = ejemplo1.IntegrateGauss(2, 3, 3, ip);
//        for (int jp=0; jp<result1.size(); jp++) {
//                    archivo1 << result1[jp] <<" ";
//        }
//        archivo1<<endl;
//    }
//         //Cerrar el archivo
//    archivo1.close();
    
    
//    NumericalIntegration solfunc1;
//    solfunc1.SetFunction(fun1);
//    auto result1 = solfunc1.IntegrateGauss(0, 2, 0, 4);
//
//    int count = 1;
//        for (int i=0; i< result1.size(); i++) {
//            double val = result1[i];
//                archivo << i <<" "<< val << std::endl;
//                count++;
//            }
//    archivo<<endl;
//
        
    
    return 0;
}


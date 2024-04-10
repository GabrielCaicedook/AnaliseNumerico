#include "NumericalIntegration.h"
#include "iostream"
#include "fstream"
using namespace std;
//#include "ExactIntegration.h"
//Incluir la clase "ExactIntegration.h"

int main(){
    
    //Definiendo funciones iniciales
    auto fun1 = [](double x) -> double {
        return x*x;
    };
    auto fun2 = [](double x) -> double {
        return 2. * (1. - x) * x;
    };
    auto fun3 = [](double x) -> double {
        return x*x;
    };
    auto fun4 = [](double x) -> double {
        
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
    
    //Creando archivo de .txt
    ofstream archivo("ExactIntegrate.txt");
    
    int np= 100;
    double a = 0;
    double b= 8/M_PI;
    double h = (b-a)/np;
    for (int ip=0; ip< np; ip++) {
        double xval = a + ip*h;
        archivo << xval <<" " <<fun4(xval) << std::endl;
        
    }
    
    
    
    
    NumericalIntegration solfunc1;
    solfunc1.SetFunction(fun1);
    auto result1 = solfunc1.IntegrateGauss(0, 2, 0, 4);
   
    int count = 1;
        for (int i=0; i< result1.size(); i++) {
            double val = result1[i];
                archivo << i <<" "<< val << std::endl;
                count++;
            }
    archivo<<endl;
            
//    //Test Numerical Integration
//    NumericalIntegration ejemplo1;
//    ejemplo1.SetFunction(fun1);
//    for (int ip=1; ip<5; ip++) {
//        auto result1 = ejemplo1.IntegrateGauss(2, 3, 3, ip);
//        for (int jp=0; jp<result1.size(); jp++) {
//                    archivo << result1[jp] <<" ";
//        }
//        archivo<<endl;
//    }
        // Cerrar el archivo
        archivo.close();
        
    
    return 0;
}


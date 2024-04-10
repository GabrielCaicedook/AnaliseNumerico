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
    
    
    //Creando archivo de .txt
    ofstream archivo("ExactIntegrate.txt");
    
    
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


#include "NumericalIntegration.h"
#include "iostream"
#include "fstream"
#include <cmath>
#include <complex>

//#ifdef USING_BOOST
//#include <boost/math/complex.hpp>
#include <boost/math/special_functions/erf.hpp>

//#endif
void ResultsFun1();
void ResultsFun2();
void ResultsFun3();
void ResultsPlot();
using namespace std;
//#include "ExactIntegration.h"
//Incluir la clase "ExactIntegration.h"

int main(){
   
     ResultsFun1();
     ResultsFun2();
     ResultsFun3();
//     ResultsPlot();
    
    return 0;
}


void ResultsFun1(){
    
    auto func3 = [](double x) -> double {
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
    
    
    //FUNCION 2
    auto func2=[](double x) -> double {
        return x*std::sin(1/x);
    };
    
    
    //FUNC 1
    auto func1 = [](double x) -> double {
            double xi=0.1;
            return std::exp(-1.0*std::pow(x - 1, 2) / xi);
        };
    //func 4 solucion de la integral indefinida
    
    
    //Creando archivo de .txt
    ofstream archivo1("AllSols.txt");
    ofstream archivo2("SolExactfun3.txt");
    ofstream archivo3("f1Table1_FUNC1_TRAPE.txt");
    ofstream archivo4("f1Table1_FUNC1_S_1_3.txt");
    ofstream archivo5("f1Table1_FUNC1_S_8_3.txt");
    ofstream archivo6("f1Table1_FUNC1_GAUSS.txt");
    
    int dig = 5;
    archivo1.precision(dig);
    archivo2.precision(dig);
    archivo3.precision(dig);
    archivo4.precision(dig);
    archivo5.precision(dig);
    archivo6.precision(dig);

    
    double test1 = cos(4);
    //Test Numerical Integration

    NumericalIntegration solve1;
    NumericalIntegration solve2;
    NumericalIntegration solve3;
    NumericalIntegration solve4;
    double a= 0;
    double b = 2;
    solve1.a = a;
    solve1.b = b;
    solve2.a = a;
    solve2.b = b;
    solve3.a = a;
    solve3.b = b;
    solve4.a = a;
    solve4.b = b;
    
    std::vector<int> refs(9);
    refs[0]=0;
    refs[1]=1;
    refs[2]=2;
    refs[3]=3;
    refs[4]=4;
    refs[5]=5;
    refs[6]=6;
    refs[7]=7;
    refs[8]=8;
  
    
    double solexactafunc1 = 0.560495;
    double solexactafunc2 = -0.000898894;
    double solexactafunc3 = 11.6391;
    
    solve1.SetnReff(refs);
    solve1.SetFunction(func1);
    solve1.SetSolExact(solexactafunc1);
    solve1.SetType(0); //SolForTrapecios.
    
    solve2.SetType(1);
    solve2.SetnReff(refs);
    solve2.SetFunction(func1);
    solve2.SetSolExact(solexactafunc1);
    
    solve3.SetType(2);
    solve3.SetnReff(refs);
    solve3.SetSolExact(solexactafunc1);
    solve3.SetFunction(func1);
    
    solve4.SetType(3);
    solve4.SetnReff(refs);
    solve4.SetSolExact(solexactafunc1);
    solve4.SetFunction(func1);
    
    solve1.Run();
    solve2.Run();
    solve3.Run();
    solve4.Run();
    
    for (int it = 0 ; it< refs.size(); it++) {
        //archivo3 << solve1.fAllErrors[it] << " "<< solve1.fAllRates[it]<<std::endl;
        archivo3 << std::log(solve1.fAllhs[it])<<" "<<solve1.fAllSols[it] << " " << solve1.fAllErrors[it]  <<" "<< solve1.fAllRates[it]<<std::endl;
        //archivo3 << " " << solve1.fsol[it] << " "<< solve1.fsolExact[it]<<std::endl;
        //archivo3 << log(solve1.fAllhs[it]) << " " << log(solve1.fAllErrors[it])<<std::endl;
    }
   // solve1.SetType(1);
    for (int it = 0 ; it< refs.size(); it++) {
        archivo4<< std::log(solve2.fAllhs[it]) <<" "<< solve2.fAllSols[it] << " " << solve2.fAllErrors[it] <<" "<< solve2.fAllRates[it]<<std::endl;
    }
    //solve1.SetType(2);
    for (int it = 0 ; it< refs.size(); it++) {
        archivo5<< std::log(solve3.fAllhs[it]) <<" "<< solve3.fAllSols[it] << " " << solve3.fAllErrors[it] <<" "<< solve3.fAllRates[it]<<std::endl;
    }
    //solve1.SetType(3);
    for (int it = 0 ; it< refs.size(); it++) {
        archivo6<< std::log(solve4.fAllhs[it]) <<" "<< solve4.fAllSols[it] << " " << solve4.fAllErrors[it] <<" "<< solve4.fAllRates[it]<<std::endl;
    }
        
}
void ResultsFun2(){
    
    auto func3 = [](double x) -> double {
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
    
    
    //FUNCION 2
    auto func2=[](double x) -> double {
        return x*std::sin(1/x);
    };
    
    
    //FUNC 1
    auto func1 = [](double x) -> double {
            double xi=0.1;
            return std::exp(-1.0*std::pow(x - 1, 2) / xi);
        };
    //func 4 solucion de la integral indefinida
    
    
    //Creando archivo de .txt
    ofstream archivo1("AllSols.txt");
    ofstream archivo2("SolExactfun3.txt");
    ofstream archivo3("f2Table1_FUNC1_TRAPE.txt");
    ofstream archivo4("f2Table1_FUNC1_S_1_3.txt");
    ofstream archivo5("f2Table1_FUNC1_S_8_3.txt");
    ofstream archivo6("f2Table1_FUNC1_GAUSS.txt");
    
    int dig = 5;
    archivo1.precision(dig);
    archivo2.precision(dig);
    archivo3.precision(dig);
    archivo4.precision(dig);
    archivo5.precision(dig);
    archivo6.precision(dig);

    
    double test1 = cos(4);
    //Test Numerical Integration

    NumericalIntegration solve1;
    NumericalIntegration solve2;
    NumericalIntegration solve3;
    NumericalIntegration solve4;
    double a= 1/100.0;
    double b = 1/10.0;
    solve1.a = a;
    solve1.b = b;
    solve2.a = a;
    solve2.b = b;
    solve3.a = a;
    solve3.b = b;
    solve4.a = a;
    solve4.b = b;
    
    std::vector<int> refs(9);
    refs[0]=0;
    refs[1]=1;
    refs[2]=2;
    refs[3]=3;
    refs[4]=4;
    refs[5]=5;
    refs[6]=6;
    refs[7]=7;
    refs[8]=8;
  
    
    double solexactafunc1 = 0.560495;
    double solexactafunc2 = -0.000898894;
    double solexactafunc3 = 11.6391;
    
    solve1.SetnReff(refs);
    solve1.SetFunction(func2);
    solve1.SetSolExact(solexactafunc2);
    solve1.SetType(0); //SolForTrapecios.
    
    solve2.SetType(1);
    solve2.SetnReff(refs);
    solve2.SetFunction(func2);
    solve2.SetSolExact(solexactafunc2);
    
    solve3.SetType(2);
    solve3.SetnReff(refs);
    solve3.SetSolExact(solexactafunc2);
    solve3.SetFunction(func2);
    
    solve4.SetType(3);
    solve4.SetnReff(refs);
    solve4.SetSolExact(solexactafunc2);
    solve4.SetFunction(func2);
    
    solve1.Run();
    solve2.Run();
    solve3.Run();
    solve4.Run();
    
    for (int it = 0 ; it< refs.size(); it++) {
        //archivo3 << solve1.fAllErrors[it] << " "<< solve1.fAllRates[it]<<std::endl;
        archivo3 << std::log(solve1.fAllhs[it])<<" "<<solve1.fAllSols[it] << " " << solve1.fAllErrors[it]  <<" "<< solve1.fAllRates[it]<<std::endl;
        //archivo3 << " " << solve1.fsol[it] << " "<< solve1.fsolExact[it]<<std::endl;
        //archivo3 << log(solve1.fAllhs[it]) << " " << log(solve1.fAllErrors[it])<<std::endl;
    }
   // solve1.SetType(1);
    for (int it = 0 ; it< refs.size(); it++) {
        archivo4<< std::log(solve2.fAllhs[it]) <<" "<< solve2.fAllSols[it] << " " << solve2.fAllErrors[it] <<" "<< solve2.fAllRates[it]<<std::endl;
    }
    //solve1.SetType(2);
    for (int it = 0 ; it< refs.size(); it++) {
        archivo5<< std::log(solve3.fAllhs[it]) <<" "<< solve3.fAllSols[it] << " " << solve3.fAllErrors[it] <<" "<< solve3.fAllRates[it]<<std::endl;
    }
    //solve1.SetType(3);
    for (int it = 0 ; it< refs.size(); it++) {
        archivo6<< std::log(solve4.fAllhs[it]) <<" "<< solve4.fAllSols[it] << " " << solve4.fAllErrors[it] <<" "<< solve4.fAllRates[it]<<std::endl;
    }
        
}
void ResultsFun3(){
    
    auto func3 = [](double x) -> double {
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
    
    
    //FUNCION 2
    auto func2=[](double x) -> double {
        return x*std::sin(1/x);
    };
    
    
    //FUNC 1
    auto func1 = [](double x) -> double {
            double xi=0.1;
            return std::exp(-1.0*std::pow(x - 1, 2) / xi);
        };
    //func 4 solucion de la integral indefinida
    
    
    //Creando archivo de .txt
    ofstream archivo1("AllSols.txt");
    ofstream archivo2("SolExactfun3.txt");
    ofstream archivo3("f3Table1_FUNC1_TRAPE.txt");
    ofstream archivo4("f3Table1_FUNC1_S_1_3.txt");
    ofstream archivo5("f3Table1_FUNC1_S_8_3.txt");
    ofstream archivo6("f3Table1_FUNC1_GAUSS.txt");
    
    int dig = 5;
    archivo1.precision(dig);
    archivo2.precision(dig);
    archivo3.precision(dig);
    archivo4.precision(dig);
    archivo5.precision(dig);
    archivo6.precision(dig);

    
    double test1 = cos(4);
    //Test Numerical Integration

    NumericalIntegration solve1;
    NumericalIntegration solve2;
    NumericalIntegration solve3;
    NumericalIntegration solve4;
    double a= 0;
    double b = 8/M_PI;
    solve1.a = a;
    solve1.b = b;
    solve2.a = a;
    solve2.b = b;
    solve3.a = a;
    solve3.b = b;
    solve4.a = a;
    solve4.b = b;
    
    std::vector<int> refs(9);
    refs[0]=0;
    refs[1]=1;
    refs[2]=2;
    refs[3]=3;
    refs[4]=4;
    refs[5]=5;
    refs[6]=6;
    refs[7]=7;
    refs[8]=8;
  
    
    double solexactafunc1 = 0.560495;
    double solexactafunc2 = -0.000898894;
    double solexactafunc3 = 11.6391;
    
    solve1.SetnReff(refs);
    solve1.SetFunction(func3);
    solve1.SetSolExact(solexactafunc3);
    solve1.SetType(0); //SolForTrapecios.
    
    solve2.SetType(1);
    solve2.SetnReff(refs);
    solve2.SetFunction(func3);
    solve2.SetSolExact(solexactafunc3);
    
    solve3.SetType(2);
    solve3.SetnReff(refs);
    solve3.SetSolExact(solexactafunc3);
    solve3.SetFunction(func3);
    
    solve4.SetType(3);
    solve4.SetnReff(refs);
    solve4.SetSolExact(solexactafunc3);
    solve4.SetFunction(func3);
    
    solve1.Run();
    solve2.Run();
    solve3.Run();
    solve4.Run();
    
    for (int it = 0 ; it< refs.size(); it++) {
        //archivo3 << solve1.fAllErrors[it] << " "<< solve1.fAllRates[it]<<std::endl;
        archivo3 << std::log(solve1.fAllhs[it])<<" "<<solve1.fAllSols[it] << " " << solve1.fAllErrors[it]  <<" "<< solve1.fAllRates[it]<<std::endl;
        //archivo3 << " " << solve1.fsol[it] << " "<< solve1.fsolExact[it]<<std::endl;
        //archivo3 << log(solve1.fAllhs[it]) << " " << log(solve1.fAllErrors[it])<<std::endl;
    }
   // solve1.SetType(1);
    for (int it = 0 ; it< refs.size(); it++) {
        archivo4<< std::log(solve2.fAllhs[it]) <<" "<< solve2.fAllSols[it] << " " << solve2.fAllErrors[it] <<" "<< solve2.fAllRates[it]<<std::endl;
    }
    //solve1.SetType(2);
    for (int it = 0 ; it< refs.size(); it++) {
        archivo5<< std::log(solve3.fAllhs[it]) <<" "<< solve3.fAllSols[it] << " " << solve3.fAllErrors[it] <<" "<< solve3.fAllRates[it]<<std::endl;
    }
    //solve1.SetType(3);
    for (int it = 0 ; it< refs.size(); it++) {
        archivo6<< std::log(solve4.fAllhs[it]) <<" "<< solve4.fAllSols[it] << " " << solve4.fAllErrors[it] <<" "<< solve4.fAllRates[it]<<std::endl;
    }
        
}

void ResultsPlot(){
    
    auto func3 = [](double x) -> double {
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
    
    
    //FUNCION 2
    auto func2=[](double x) -> double {
        return x*std::sin(1/x);
    };
    
    
    //FUNC 1
    auto func1 = [](double x) -> double {
            double xi=0.1;
            return std::exp(-1.0*std::pow(x - 1, 2) / xi);
        };
    //func 4 solucion de la integral indefinida
    
    
    //Creando archivo de .txt
    ofstream archivo1("AllSols.txt");
    ofstream archivo2("SolExactfun3.txt");
    ofstream archivo3("f4Table1_FUNC1_TRAPE.txt");
    ofstream archivo4("f4Table1_FUNC1_S_1_3.txt");
    ofstream archivo5("f4Table1_FUNC1_S_8_3.txt");
    ofstream archivo6("P1_f4Table1_FUNC1_GAUSS.txt");
    
    int dig = 5;
    archivo1.precision(dig);
    archivo2.precision(dig);
    archivo3.precision(dig);
    archivo4.precision(dig);
    archivo5.precision(dig);
    archivo6.precision(dig);

    
    double test1 = cos(4);
    //Test Numerical Integration

    NumericalIntegration solve1;
    NumericalIntegration solve2;
    NumericalIntegration solve3;
    NumericalIntegration solve4;
    double a= 1/100.0;
    double b = 1/10.0;
    solve1.a = a;
    solve1.b = b;
    solve2.a = a;
    solve2.b = b;
    solve3.a = a;
    solve3.b = b;
    solve4.a = a;
    solve4.b = b;
    
    std::vector<int> refs(9);
    refs[0]=0;
    refs[1]=1;
    refs[2]=2;
    refs[3]=3;
    refs[4]=4;
    refs[5]=5;
    refs[6]=6;
    refs[7]=7;
    refs[8]=8;
  
    
    double solexactafunc1 = 0.560495;
    double solexactafunc2 = -0.000898894;
    double solexactafunc3 = 11.6391;
    
    solve1.SetnReff(refs);
    solve1.SetFunction(func2);
    solve1.SetSolExact(solexactafunc2);
    solve1.SetType(0); //SolForTrapecios.
    
    solve2.SetType(1);
    solve2.SetnReff(refs);
    solve2.SetFunction(func2);
    solve2.SetSolExact(solexactafunc2);
    
    solve3.SetType(2);
    solve3.SetnReff(refs);
    solve3.SetSolExact(solexactafunc2);
    solve3.SetFunction(func2);
    
    solve4.SetType(3);
    solve4.SetnReff(refs);
    solve4.SetSolExact(solexactafunc2);
    solve4.SetFunction(func2);
    
    solve1.Run();
    solve2.Run();
    solve3.Run();
    solve4.Run();
    
    for (int it = 0 ; it< refs.size(); it++) {
        //archivo3 << solve1.fAllErrors[it] << " "<< solve1.fAllRates[it]<<std::endl;
        archivo6 << solve1.fAllhs[it]<<" "<<solve1.fAllSols[it] << " " << solve1.fAllErrors[it]  <<" "<< solve1.fAllRates[it]<<std::endl;
        //archivo3 << " " << solve1.fsol[it] << " "<< solve1.fsolExact[it]<<std::endl;
        //archivo3 << log(solve1.fAllhs[it]) << " " << log(solve1.fAllErrors[it])<<std::endl;
    }
  
}


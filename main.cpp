#include "NumericalIntegration.h"
#include "iostream"
#include "fstream"
#include <cmath>
#include <complex>

//#ifdef USING_BOOST
//#include <boost/math/complex.hpp>
#include <boost/math/special_functions/erf.hpp>

//#endif

using namespace std;
//#include "ExactIntegration.h"
//Incluir la clase "ExactIntegration.h"

int main(){
    
//    auto func9 = [](double x) -> double {
//        double eps = 0.1;
//        // Utiliza la función error estándar con un número complejo
//          std::complex<double> z(0, x);
//          double test =0.1;
//          auto result = boost::math::erf(z);
//
//          // La parte imaginaria del resultado de erf(z) es equivalente a erfi(x)
//   //       return result.imag();
//    };
    
    
    
//#if USING_BOOST
//    std::cout<<"ok***********"<<std::endl;
//#endif
    //FUNCION 3
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
    //FUNCION 2
    auto func5=[](double x) -> double {
        return x*std::sin(1/x);
    };
    auto func6 =[](double x)-> double {
        return 0.5 * x * std::cos(1 / x) + 0.5 * x * x * std::sin(1 / x) + 0.5 * std::sinh(1/x);
        ;
    };
    
    //FUNC 1
    auto func3 = [](double x) -> double {
            double xi=0.1;
            return std::exp(-1.0*std::pow(x - 1, 2) / xi);
        };
    //func 4 solucion de la integral indefinida
    auto func4 = [](double x) -> double {
//        double eps = 0.1; // Suponiendo que eps es 0.1, ajusta según sea necesario.
//        double term = (-1 + x) / sqrt(eps);
//
//        // Usa std::complex en lugar de boost::math::complex
//        std::complex<double> z(0, term); // crea un número complejo con parte real 0 y parte imaginaria term
//
//        // Llama a erf con el número complejo y toma la parte imaginaria del resultado
//        double result = 0.5 * sqrt(eps) * sqrt(M_PI) * boost::math::erf(z).imag();
//
//        std::cout << "Result: " << result << std::endl;
    };
    
    
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
//    ofstream archivo3("List1_Func3.txt");
//    ofstream archivo("Soltrapeciofunc3.txt");
    
//    int np= 20;
//    double a = 0;
//    double b= 8/M_PI;
//    double h = (b-a)/np;
//    for (int ip=0; ip< np; ip++) {
//        double xval = a + ip*h;
//        double val = func2(xval);
//        if(xval < 2.0/M_PI){
//            val = -1.0;
//        }
//        archivo1 << xval <<" " <<val<< std::endl;
//
//    }
    
    double test1 = cos(4);
    //Test Numerical Integration

    NumericalIntegration solve1;
    NumericalIntegration solve2;
    NumericalIntegration solve3;
    NumericalIntegration solve4;
    double a= 1.0/100.0;
    double b = 1.0/10.0;
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
    
    solve1.SetnReff(refs);
    solve1.SetFunction(func5);
    solve1.SetSolExact(solexactafunc2);
    solve1.SetType(0); //SolForTrapecios.
    
    solve2.SetType(1);
    solve2.SetnReff(refs);
    solve2.SetFunction(func5);
    solve2.SetFunctionExact(func2);
    
    solve3.SetType(2);
    solve3.SetnReff(refs);
    solve3.SetSolExact(solexactafunc2);
    solve3.SetFunction(func5);
    
    solve4.SetType(3);
    solve4.SetnReff(refs);
    solve4.SetSolExact(solexactafunc2);
    solve4.SetFunction(func5);
    
    solve1.Run();
    solve2.Run();
    solve3.Run();
    solve4.Run();
    
    //archivo3 << solve1.fsol1<< std::endl;
    
    //solve1.SetType(1);
   
    
    //archivo1 << solve1.fAllhs[0] << " " <<"0"<<std::endl;
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
    
    int ok=0;
    return 0;
//    for (int ip=2; ip<9; ip++) {
//        auto result1 = solve1.IntegrateTrapecio(0,(8/M_PI),ip);
//        auto result2 = solve1.ValExact(0,(8/M_PI),ip);
//        solve1.CalcError();
//        ip=ip+1;
//        for (int jp=0; jp<result1.size(); jp++) {
//                    archivo1 << result1[jp] <<" ";
//        }
//        archivo1<<endl;
//    }
//         //Cerrar el archivo
//    archivo1.close();
//
//    NumericalIntegration solve2;
//    solve2.SetFunction(func2);
//    for (int ip=2; ip<9; ip++) {
//        auto result2 = solve2.ValExact(0,(8/M_PI),ip);
//        ip=ip+1;
//        for (int jp=0; jp<result2.size(); jp++) {
//                archivo2 << result2[jp] <<" ";
//        }
//        archivo2<<endl;
//    }
//         //Cerrar el archivo
//    archivo2.close();
    NumericalIntegration testhreturn;
    for (int i=2; i<10;i++){
        auto vec =testhreturn.Hreturn(0, 2, i);
        for (int j=0; j<9; j++){
           // archivo2<<vec[i]<<""<<std::endl;
            
        }
    }
    
    
    
    
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


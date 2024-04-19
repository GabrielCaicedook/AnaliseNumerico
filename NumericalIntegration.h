
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;

class NumericalIntegration{
public:
    //0 Trapecios 1Simpsons 2...
    int ftype;
    double fsolsum;
    double a=0;
    double b=0;
    double fsol1;
    std::function<double(double)> funcion;
    std::function<double(double)> funcionExact;
    std::vector<double> fsol;
    std::vector<double> fsolExact;
    std::vector<double> ferror;
    std::vector<double> fAllErrors;
    std::vector<double> fAllRates;
    std::vector<double> fAllSols;
    std::vector<double> fAllhs;
    std::vector<int> fnreffs;
    
    double fsolExactval = 0.0;
   
    void SetFunction(std::function<double( double)> fun );
    void SetFunctionExact(std::function<double( double)> fun );
    void SetnReff(std::vector<int> reffs){
        fnreffs = reffs;
    }
    void SetSolExact(double val){
        fsolExactval = val;
    }
    void Run(){
        fAllErrors.resize(fnreffs.size());
        fAllhs.resize(fnreffs.size());
        fAllSols.resize(fnreffs.size());
        
        std::ofstream solution("gaussolAprox.txt");
        std::ofstream solution2("gaussolExact.txt");
        for (int i = 0; i< fnreffs.size(); i++) {
            //if type == 0 -> Trapecios
            if (ftype == 0) {
               
                auto solAreas = IntegrateTrapecio(a, b, fnreffs[i]);
                fAllhs[i] = Hreturn(a, b, fnreffs[i]);
               // ValExact(a, b, fnreffs[i]);
                double sum = CalcError();
                fAllErrors[i] = sum;
                fAllSols[i] = SumAreas(solAreas);
                
            }
            if (ftype == 1) {
                auto solAreas = IntegrateSimpson1(a, b, fnreffs[i]);
                fAllhs[i] = Hreturn(a, b, fnreffs[i]);
             //   ValExact(a, b, fnreffs[i]);
                double sum = CalcError();
                fAllErrors[i] = sum;
                fAllSols[i] = SumAreas(solAreas);
               
            }
            if (ftype == 2) {
                auto solAreas = IntegrateSimpson2(a, b, fnreffs[i]);
                fAllhs[i] = Hreturn(a, b, fnreffs[i]);
              //  ValExact(a, b, fnreffs[i]);
                double sum = CalcError();
                fAllErrors[i] = sum;
                fAllSols[i] = SumAreas(solAreas);
              
            }
            if (ftype == 3) {
                
                
                auto solAreas = IntegrateGauss(a, b, fnreffs[i], 4
                                               );
                fAllhs[i] = Hreturn(a, b, fnreffs[i]);
             //   ValExact(a, b, fnreffs[i]);
//                for (int isol=0; isol< fsol.size(); isol++) {
//                    solution<< fsol[isol] << " ";
//                    solution2<< fsolExact[isol] << " ";
//                }
    
                double sum = CalcError();
                fAllErrors[i] = sum;
                fAllSols[i] = SumAreas(solAreas);
               
            }
            
            
        }
        CalcRates();
        
    }
    
    
    double CalcError();
    double CalcRates(){
        
        fAllRates.resize(fAllErrors.size() -1);
        for(int i=0; i< fAllRates.size(); i++ ){
            double e1 = fAllErrors[i];
            double e2 = fAllErrors[i+1];
            double num = log(e1/e2);
            //double num = log(fAllErrors[i];/fAllErrors[i+1]);
            double h1 = (1/pow(2,fnreffs[i]));
            double h2 = (1/pow(2,fnreffs[i+1]));
            double den = log(h1/h2);
            //double den = log((pow(2,fnreffs[i]))/(pow(2,fnreffs[i+1])));
            fAllRates[i] = num/den;
            
        }
    }
    
    NumericalIntegration();
    NumericalIntegration(int type);
    void SetType(int type);
    
    //HRETURN
    double Hreturn(double x1, double x2, int ref);
   // std::vector<double> Hreturn(double a, int b, int ref);
    
    
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
    double SumAreas(std::vector<double> areas);
};

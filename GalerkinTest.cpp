#include "NumericalIntegration.h"
#include "iostream"
#include "fstream"
#include <cmath>
#include <vector>
#include <complex>


//#ifdef USING_BOOST
//#include <boost/math/complex.hpp>
#include <boost/math/special_functions/erf.hpp>
void CreateSpace(double x, double y, int order, std::vector<double> &val, std::vector<double> &derivx, std::vector<double> &derivy);
void PrintMaTrix(std::vector<std::vector<double>> Mat);
using Matrix = std::vector<std::vector<double>>;
std::vector<double> solveSystem(const Matrix& A, const std::vector<double>& B);
std::vector<double> multiplyMatrixVector(const Matrix& matrix, const std::vector<double>& vec);
bool invertMatrix(const Matrix& A, Matrix& inverse);
void printMatrix(const Matrix& matrix);
Matrix createIdentityMatrix(int n);
std::vector<double> CalculateAlphas(int orderval);
double IntegrateSolution(int orderval, std::vector<double> alphas);
double exactSol(double x, double y);
int main(){
   
  
    
    //intsol retorna la norma del error
    std::ofstream file("normError.txt");
    std::ofstream file2("alphas.txt");
    for(int iorder=1; iorder<=5; iorder++){
        
        std::cout<<" ******************************************"<<std::endl;
        std::cout<<" Simulation for order: "<<iorder<<std::endl;
        std::vector<double> sol = CalculateAlphas(iorder);
        int nalpha = sol.size();
        for (int ialpha=0; ialpha<nalpha; ialpha++) {
            std::cout<<"alpha "<<ialpha<<" ="<<sol[ialpha]<<std::endl;
            file2<<sol[ialpha]<<" ";
        }
        file2<<std::endl;
        double normerror = IntegrateSolution(iorder, sol);
        file<<iorder<<" "<<normerror<<std::endl;
        std::cout<<" order: "<<iorder<<" normError: "<<normerror<<std::endl;
        std::cout<<" ******************************************"<<std::endl;
    }
   
    
    return 0;

}
std::vector<double> CalculateAlphas(int orderval){
    int order =orderval;
    std::vector<double> val;
    std::vector<double> derivx;
    std::vector<double> derivy;
    CreateSpace(0.5, 0.5, order, val, derivx, derivy);
    std::vector<double> Points={-0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
        -0.4333953941292472, -0.1488743389816312, 0.1488743389816312,
          0.4333953941292472, 0.6794095682990244, 0.8650633666889845,
          0.9739065285171717};
    std::vector<double> Pesos={
        0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529,
        0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881
    };
    double test1=5.0;
    

    
    int nrows = val.size();
    std::vector<std::vector<double>> A(nrows, std::vector<double>(nrows, 0.0));;
    std::vector<double> B(nrows, 0.0);
    //Create the Matrix
    
    for (int i=0; i<nrows; i++) {
        for (int j=0; j<nrows; j++) {
            double suma =0.0;
            for (int ipoint=0; ipoint<Points.size(); ipoint++) {
                double pointx = Points[ipoint];
                for (int jpoint=0; jpoint<Points.size(); jpoint++) {
                    double pointy = Points[jpoint];
                    double peso = Pesos[jpoint];
                    CreateSpace(pointx, pointy, order, val, derivx, derivy);
                    double fun = derivx[i]*derivx[j] + derivy[i]*derivy[j];
                    suma += Pesos[ipoint] * Pesos[jpoint] * fun;
                }
            }
            A[i][j]=suma;
        };
        double suma =0.0;
        double Rbc =0.0;
        double Tbc =0.0;
        for (int ipoint=0; ipoint<Points.size(); ipoint++) {
            double pointx = Points[ipoint];
            for (int jpoint=0; jpoint<Points.size(); jpoint++) {
                double pointy = Points[jpoint];
                double peso = Pesos[jpoint];
                CreateSpace(pointx, pointy, order, val, derivx, derivy);
                
                double term1 = (2 * (1 + pointx) * (-3 + 2 * pointy) * atan(1 - pointx)) / std::pow((2 + (-2 + pointy) * pointy), 2);
                double term2 = (2 * (-3 + 2 * pointx) * (1 + pointy) * atan(1 - pointy)) / std::pow((2 + (-2 + pointx) * pointx), 2);
                double fun = -1.0*term1 - term2;
                suma += Pesos[ipoint] * Pesos[jpoint] * fun*val[i];
           
            }
            //point x = point y
            CreateSpace(1.0, pointx, order, val, derivx, derivy);
            double pointxAux = pointx;
            double pointy =pointxAux;
            double fun = -2 *(1.0 + pointy) *atan(1 -pointy);
            Rbc += Pesos[ipoint] * fun * val[i];
            
            CreateSpace(pointx, 1.0, order, val, derivx, derivy);
            fun = -2*(1 +pointxAux) *atan(1 - pointxAux);
            Tbc += Pesos[ipoint] * fun*val[i];
            
        }
        B[i]= suma + Tbc + Rbc  ;
//        B[i]= suma   ;
      
    }
    
   
    std::vector<double> solution = solveSystem(A, B);
    return solution;
}

double IntegrateSolution(int orderval, std::vector<double> alphas){
    int order =orderval;
    std::vector<double> val;
    std::vector<double> derivx;
    std::vector<double> derivy;
    CreateSpace(0.5, 0.5, order, val, derivx, derivy);
    
    std::vector<double> Points={-0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
        -0.4333953941292472, -0.1488743389816312, 0.1488743389816312,
          0.4333953941292472, 0.6794095682990244, 0.8650633666889845,
          0.9739065285171717};
    std::vector<double> Pesos={
        0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529,
        0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881
    };
    double test1=5.0;
    

    
    int nrows = val.size();
    std::vector<std::vector<double>> A(nrows, std::vector<double>(nrows, 0.0));;
    std::vector<double> B(nrows, 0.0);
    //Create the Matrix
    
    
        double suma =0.0;
    double exactVal =0.0;
    for(int iphi=0; iphi< nrows; iphi++){
        for (int ipoint=0; ipoint<Points.size(); ipoint++) {
            double pointx = Points[ipoint];
            for (int jpoint=0; jpoint<Points.size(); jpoint++) {
                double pointy = Points[jpoint];
                double peso = Pesos[jpoint];
                CreateSpace(pointx, pointy, order, val, derivx, derivy);
                double fun = alphas[iphi]*val[iphi];
                suma += Pesos[ipoint] * Pesos[jpoint] * fun;
                if(iphi==0){
                    exactVal += Pesos[ipoint] * Pesos[jpoint] * exactSol(pointx, pointy);
                }
                
            }
        }
        
    }
    return sqrt((suma - exactVal)*(suma - exactVal));
}

void PrintMaTrix(std::vector<std::vector<double>> Mat){
   
    
    int nrows = Mat.size();
    for (int i=0; i<nrows; i++) {
        for (int j=0; j<nrows; j++) {
            std::cout<<" "<< Mat[i][j];
        };
        std::cout<<""<<std::endl;
    }
    

    
    return 0;
}

void CreateSpace(double x, double y, int order, std::vector<double> &val, std::vector<double> &derivx, std::vector<double> &derivy){
    if(order ==1){
        val.resize(1);
        derivx.resize(1);
        derivy.resize(1);
        val[0] = (1+x)*(1+y);
        derivx[0] = (1+y);
        derivy[0] = (1+x);
    }
    if(order ==2){
        val.resize(4);
        derivx.resize(4);
        derivy.resize(4);

        val[0] = (1 + x) *(1 + y);
        val[1] = (1 + x) *y* (1 + y);
        val[2] = x *(1 + x) *(1 + y);
        val[3] = x * (1 + x) *y* (1 + y);

        derivx[0] = (1+y);
        derivx[1] = y * (1 + y);
        derivx[2] = x * (1 + y) + (1 + x) *(1 + y);
        derivx[3] = x * y* (1 + y) + (1 + x) *y* (1 + y);

        derivy[0] = 1 + x;
        derivy[1] = (1 + x) *y + (1 + x)* (1 + y);
        derivy[2] = x *(1 + x);
        derivy[3] =  x *(1 + x) *y + x* (1 + x)* (1 + y);
    }

    if(order ==3){
        val.resize(9);
           derivx.resize(9);
           derivy.resize(9);

           // Llenado de los valores de las funciones
            val[0] = (1 + x) * (1 + y);
            val[1] = (1 + x) * y * (1 + y);
            val[2] = 0.5 * (1 + x) * (1 + y) * (-1 + 3 * y * y);
            val[3] = x * (1 + x) * (1 + y);
            val[4] = x * (1 + x) * y * (1 + y);
            val[5] = 0.5 * x * (1 + x) * (1 + y) * (-1 + 3 * y * y);
            val[6] = 0.5 * (1 + x) * (-1 + 3 * x * x) * (1 + y);
            val[7] = 0.5 * (1 + x) * (-1 + 3 * x * x) * y * (1 + y);
            val[8] = 0.25 * (1 + x) * (-1 + 3 * x * x) * (1 + y) * (-1 + 3 * y * y);

           // Llenado de las derivadas en x
            derivx[0] = 1 + y;
            derivx[1] = y * (1 + y);
            derivx[2] = 0.5 * (1 + y) * (-1 + 3 * y * y);
            derivx[3] = x * (1 + y) + (1 + x) * (1 + y);
            derivx[4] = x * y * (1 + y) + (1 + x) * y * (1 + y);
            derivx[5] = 0.5 * x * (1 + y) * (-1 + 3 * y * y) + 0.5 * (1 + x) * (1 + y) * (-1 + 3 * y * y);
            derivx[6] = 3 * x * (1 + x) * (1 + y) + 0.5 * (-1 + 3 * x * x) * (1 + y);
            derivx[7] = 3 * x * (1 + x) * y * (1 + y) + 0.5 * (-1 + 3 * x * x) * y * (1 + y);
            derivx[8] = (3.0 / 2.0) * x * (1 + x) * (1 + y) * (-1 + 3 * y * y) +
        (1.0 / 4.0) * (-1 + 3 * x * x) * (1 + y) * (-1 + 3 * y * y);

           // Llenado de las derivadas en y
            derivy[0] = 1 + x;
            derivy[1] = (1 + x) * y + (1 + x) * (1 + y);
            derivy[2] = 3 * (1 + x) * y * (1 + y) + 0.5 * (1 + x) * (-1 + 3 * y * y);
            derivy[3] = x * (1 + x);
            derivy[4] = x * (1 + x) * y + x * (1 + x) * (1 + y);
            derivy[5] = 3 * x * (1 + x) * y * (1 + y) + 0.5 * x * (1 + x) * (-1 + 3 * y * y);
            derivy[6] = 0.5 * (1 + x) * (-1 + 3 * x * x);
            derivy[7] = 0.5 * (1 + x) * (-1 + 3 * x * x) * y + 0.5 * (1 + x) * (-1 + 3 * x * x) * (1 + y);
            derivy[8] = 1.5 * (1 + x) * (-1 + 3 * x * x) * y * (1 + y) + 0.25 * (1 + x) * (-1 + 3 * x * x) * (-1 + 3 * y * y);
    }
    if(order ==4){

        val.resize(16);
            derivx.resize(16);
            derivy.resize(16);

            // Llenado de los valores de las funciones
            val[0] = (1 + x) * (1 + y);
            val[1] = (1 + x) * y * (1 + y);
            val[2] = 0.5 * (1 + x) * (1 + y) * (-1 + 3 * y * y);
            val[3] = 0.5 * (1 + x) * y * (1 + y) * (-3 + 5 * y * y);
            val[4] = x * (1 + x) * (1 + y);
            val[5] = x * (1 + x) * y * (1 + y);
            val[6] = 0.5 * x * (1 + x) * (1 + y) * (-1 + 3 * y * y);
            val[7] = 0.5 * x * (1 + x) * y * (1 + y) * (-3 + 5 * y * y);
            val[8] = 0.5 * (1 + x) * (-1 + 3 * x * x) * (1 + y);
            val[9] = 0.5 * (1 + x) * (-1 + 3 * x * x) * y * (1 + y);
            val[10] = 0.25 * (1 + x) * (-1 + 3 * x * x) * (1 + y) * (-1 + 3 * y * y);
            val[11] = 0.25 * (1 + x) * (-1 + 3 * x * x) * y * (1 + y) * (-3 + 5 * y * y);
            val[12] = 0.5 * x * (1 + x) * (-3 + 5 * x * x) * (1 + y);
            val[13] = 0.5 * x * (1 + x) * (-3 + 5 * x * x) * y * (1 + y);
            val[14] = 0.25 * x * (1 + x) * (-3 + 5 * x * x) * (1 + y) * (-1 + 3 * y * y);
            val[15] = 0.25 * x * (1 + x) * (-3 + 5 * x * x) * y * (1 + y) * (-3 + 5 * y * y);


            // Llenado de las derivadas en x
            derivx[0] = 1 + y;
            derivx[1] = y * (1 + y);
            derivx[2] = 0.5 * (1 + y) * (-1 + 3 * y * y);
            derivx[3] = 0.5 * y * (1 + y) * (-3 + 5 * y * y);
            derivx[4] = x * (1 + y) + (1 + x) * (1 + y);
            derivx[5] = x * y * (1 + y) + (1 + x) * y * (1 + y);
            derivx[6] = 0.5 * x * (1 + y) * (-1 + 3 * y * y) + 0.5 * (1 + x) * (1 + y) * (-1 + 3 * y * y);
            derivx[7] = 0.5 * x * y * (1 + y) * (-3 + 5 * y * y) + 0.5 * (1 + x) * y * (1 + y) * (-3 + 5 * y * y);
            derivx[8] = 3 * x * (1 + x) * (1 + y) + 0.5 * (-1 + 3 * x * x) * (1 + y);
            derivx[9] = 3 * x * (1 + x) * y * (1 + y) + 0.5 * (-1 + 3 * x * x) * y * (1 + y);
            derivx[10] = 1.5 * x * (1 + x) * (1 + y) * (-1 + 3 * y * y) + 0.25 * (-1 + 3 * x * x) * (1 + y) * (-1 + 3 * y * y);
            derivx[11] = 1.5 * x * (1 + x) * y * (1 + y) * (-3 + 5 * y * y) + 0.25 * (-1 + 3 * x * x) * y * (1 + y) * (-3 + 5 * y * y);
            derivx[12] = 5 * x * x * (1 + x) * (1 + y) + 0.5 * x * (-3 + 5 * x * x) * (1 + y) + 0.5 * (1 + x) * (-3 + 5 * x * x) * (1 + y);
            derivx[13] = 5 * x * x * (1 + x) * y * (1 + y) + 0.5 * x * (-3 + 5 * x * x) * y * (1 + y) + 0.5 * (1 + x) * (-3 + 5 * x * x) * y * (1 + y);
            derivx[14] = 2.5 * x * x * (1 + x) * (1 + y) * (-1 + 3 * y * y) + 0.25 * x * (-3 + 5 * x * x) * (1 + y) * (-1 + 3 * y * y) + 0.25 * (1 + x) * (-3 + 5 * x * x) * (1 + y) * (-1 + 3 * y * y);
            derivx[15] = 2.5 * x * x * (1 + x) * y * (1 + y) * (-3 + 5 * y * y) + 0.25 * x * (-3 + 5 * x * x) * y * (1 + y) * (-3 + 5 * y * y) + 0.25 * (1 + x) * (-3 + 5 * x * x) * y * (1 + y) * (-3 + 5 * y * y);


            // Llenado de las derivadas en y
            derivy[0] = 1 + x;
            derivy[1] = (1 + x) * y + (1 + x) * (1 + y);
            derivy[2] = 3 * (1 + x) * y * (1 + y) + 0.5 * (1 + x) * (-1 + 3 * y * y);
            derivy[3] = 5 * (1 + x) * y * y * (1 + y) + 0.5 * (1 + x) * y * (-3 + 5 * y * y) + 0.5 * (1 + x) * (1 + y) * (-3 + 5 * y * y);
            derivy[4] = x * (1 + x);
            derivy[5] = x * (1 + x) * y + x * (1 + x) * (1 + y);
            derivy[6] = 3 * x * (1 + x) * y * (1 + y) + 0.5 * x * (1 + x) * (-1 + 3 * y * y);
            derivy[7] = 5 * x * (1 + x) * y * y * (1 + y) + 0.5 * x * (1 + x) * y * (-3 + 5 * y * y) + 0.5 * x * (1 + x) * (1 + y) * (-3 + 5 * y * y);
            derivy[8] = 0.5 * (1 + x) * (-1 + 3 * x * x);
            derivy[9] = 0.5 * (1 + x) * (-1 + 3 * x * x) * y + 0.5 * (1 + x) * (-1 + 3 * x * x) * (1 + y);
            derivy[10] = 1.5 * (1 + x) * (-1 + 3 * x * x) * y * (1 + y) + 0.25 * (1 + x) * (-1 + 3 * x * x) * (-1 + 3 * y * y);
            derivy[11] = 2.5 * (1 + x) * (-1 + 3 * x * x) * y * y * (1 + y) + 0.25 * (1 + x) * (-1 + 3 * x * x) * y * (-3 + 5 * y * y) + 0.25 * (1 + x) * (-1 + 3 * x * x) * (1 + y) * (-3 + 5 * y * y);
            derivy[12] = 0.5 * x * (1 + x) * (-3 + 5 * x * x);
            derivy[13] = 0.5 * x * (1 + x) * (-3 + 5 * x * x) * y + 0.5 * x * (1 + x) * (-3 + 5 * x * x) * (1 + y);
            derivy[14] = 1.5 * x * (1 + x) * (-3 + 5 * x * x) * y * (1 + y) + 0.25 * x * (1 + x) * (-3 + 5 * x * x) * (-1 + 3 * y * y);
            derivy[15] = 2.5 * x * (1 + x) * (-3 + 5 * x * x) * y * y * (1 + y) + 0.25 * x * (1 + x) * (-3 + 5 * x * x) * y * (-3 + 5 * y * y) + 0.25 * x * (1 + x) * (-3 + 5 * x * x) * (1 + y) * (-3 + 5 * y * y);


    }
    if(order ==5){
        val.resize(25);
           derivx.resize(25);
           derivy.resize(25);

        // Llenar el vector con las posiciones de la matriz proporcionada
            val[0] = (1 + x) * (1 + y);
            val[1] = (1 + x) * y * (1 + y);
            val[2] = 0.5 * (1 + x) * (1 + y) * (-1 + 3 * y * y);
            val[3] = 0.5 * (1 + x) * y * (1 + y) * (-3 + 5 * y * y);
            val[4] = 0.125 * (1 + x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y);
            val[5] = x * (1 + x) * (1 + y);
            val[6] = x * (1 + x) * y * (1 + y);
            val[7] = 0.5 * x * (1 + x) * (1 + y) * (-1 + 3 * y * y);
            val[8] = 0.5 * x * (1 + x) * y * (1 + y) * (-3 + 5 * y * y);
            val[9] = 0.125 * x * (1 + x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y);
            val[10] = 0.5 * (1 + x) * (-1 + 3 * x * x) * (1 + y);
            val[11] = 0.5 * (1 + x) * (-1 + 3 * x * x) * y * (1 + y);
            val[12] = 0.25 * (1 + x) * (-1 + 3 * x * x) * (1 + y) * (-1 + 3 * y * y);
            val[13] = 0.25 * (1 + x) * (-1 + 3 * x * x) * y * (1 + y) * (-3 + 5 * y * y);
            val[14] = 0.0625 * (1 + x) * (-1 + 3 * x * x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y);
            val[15] = 0.5 * x * (1 + x) * (-3 + 5 * x * x) * (1 + y);
            val[16] = 0.5 * x * (1 + x) * (-3 + 5 * x * x) * y * (1 + y);
            val[17] = 0.25 * x * (1 + x) * (-3 + 5 * x * x) * (1 + y) * (-1 + 3 * y * y);
            val[18] = 0.25 * x * (1 + x) * (-3 + 5 * x * x) * y * (1 + y) * (-3 + 5 * y * y);
            val[19] = 0.0625 * x * (1 + x) * (-3 + 5 * x * x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y);
            val[20] = 0.125 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * (1 + y);
            val[21] = 0.125 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * y * (1 + y);
            val[22] = 0.0625 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * (1 + y) * (-1 + 3 * y * y);
            val[23] = 0.0625 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * y * (1 + y) * (-3 + 5 * y * y);
            val[24] = 0.015625 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y);



           // Llenado de las derivadas en x
            derivx[0] = 1 + y;
            derivx[1] = y * (1 + y);
            derivx[2] = 0.5 * (1 + y) * (-1 + 3 * y * y);
            derivx[3] = 0.5 * y * (1 + y) * (-3 + 5 * y * y);
            derivx[4] = 0.125 * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y);
            derivx[5] = x * (1 + y) + (1 + x) * (1 + y);
            derivx[6] = x * y * (1 + y) + (1 + x) * y * (1 + y);
            derivx[7] = 0.5 * x * (1 + y) * (-1 + 3 * y * y) + 0.5 * (1 + x) * (1 + y) * (-1 + 3 * y * y);
            derivx[8] = 0.5 * x * y * (1 + y) * (-3 + 5 * y * y) + 0.5 * (1 + x) * y * (1 + y) * (-3 + 5 * y * y);
            derivx[9] = 0.125 * x * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y) + 0.125 * (1 + x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y);
            derivx[10] = 3 * x * (1 + x) * (1 + y) + 0.5 * (-1 + 3 * x * x) * (1 + y);
            derivx[11] = 3 * x * (1 + x) * y * (1 + y) + 0.5 * (-1 + 3 * x * x) * y * (1 + y);
            derivx[12] = 1.5 * x * (1 + x) * (1 + y) * (-1 + 3 * y * y) + 0.25 * (-1 + 3 * x * x) * (1 + y) * (-1 + 3 * y * y);
            derivx[13] = 1.5 * x * (1 + x) * y * (1 + y) * (-3 + 5 * y * y) + 0.25 * (-1 + 3 * x * x) * y * (1 + y) * (-3 + 5 * y * y);
            derivx[14] = 0.375 * x * (1 + x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y) + 0.0625 * (-1 + 3 * x * x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y);
            derivx[15] = 5 * x * x * (1 + x) * (1 + y) + 0.5 * x * (-3 + 5 * x * x) * (1 + y) + 0.5 * (1 + x) * (-3 + 5 * x * x) * (1 + y);
            derivx[16] = 5 * x * x * (1 + x) * y * (1 + y) + 0.5 * x * (-3 + 5 * x * x) * y * (1 + y) + 0.5 * (1 + x) * (-3 + 5 * x * x) * y * (1 + y);
            derivx[17] = 2.5 * x * x * (1 + x) * (1 + y) * (-1 + 3 * y * y) + 0.25 * x * (-3 + 5 * x * x) * (1 + y) * (-1 + 3 * y * y) + 0.25 * (1 + x) * (-3 + 5 * x * x) * (1 + y) * (-1 + 3 * y * y);
            derivx[18] = 2.5 * x * x * (1 + x) * y * (1 + y) * (-3 + 5 * y * y) + 0.25 * x * (-3 + 5 * x * x) * y * (1 + y) * (-3 + 5 * y * y) + 0.25 * (1 + x) * (-3 + 5 * x * x) * y * (1 + y) * (-3 + 5 * y * y);
            derivx[19] = 0.625 * x * x * (1 + x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y) + 0.0625 * x * (-3 + 5 * x * x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y) + 0.0625 * (1 + x) * (-3 + 5 * x * x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y);
            derivx[20] = 0.125 * (1 + x) * (-60 * x + 140 * x * x * x) * (1 + y) + 0.125 * (3 - 30 * x * x + 35 * x * x * x * x) * (1 + y);
            derivx[21] = 0.125 * (1 + x) * (-60 * x + 140 * x * x * x) * y * (1 + y) + 0.125 * (3 - 30 * x * x + 35 * x * x * x * x) * y * (1 + y);
            derivx[22] = 0.0625 * (1 + x) * (-60 * x + 140 * x * x * x) * (1 + y) * (-1 + 3 * y * y) + 0.0625 * (3 - 30 * x * x + 35 * x * x * x * x) * (1 + y) * (-1 + 3 * y * y);
            derivx[23] = 0.0625 * (1 + x) * (-60 * x + 140 * x * x * x) * y * (1 + y) * (-3 + 5 * y * y) + 0.0625 * (3 - 30 * x * x + 35 * x * x * x * x) * y * (1 + y) * (-3 + 5 * y * y);
            derivx[24] = 0.015625 * (1 + x) * (-60 * x + 140 * x * x * x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y) + 0.015625 * (3 - 30 * x * x + 35 * x * x * x * x) * (1 + y) * (3 - 30 * y * y + 35 * y * y * y * y);

           // Llenado de las derivadas en y
            derivy[0] = 1 + x;
            derivy[1] = (1 + x) * y + (1 + x) * (1 + y);
            derivy[2] = 3 * (1 + x) * y * (1 + y) + 0.5 * (1 + x) * (-1 + 3 * y * y);
            derivy[3] = 5 * (1 + x) * y * y * (1 + y) + 0.5 * (1 + x) * y * (-3 + 5 * y * y) + 0.5 * (1 + x) * (1 + y) * (-3 + 5 * y * y);
            derivy[4] = 0.125 * (1 + x) * (1 + y) * (-60 * y + 140 * y * y * y) + 0.125 * (1 + x) * (3 - 30 * y * y + 35 * y * y * y * y);
            derivy[5] = x * (1 + x);
            derivy[6] = x * (1 + x) * y + x * (1 + x) * (1 + y);
            derivy[7] = 3 * x * (1 + x) * y * (1 + y) + 0.5 * x * (1 + x) * (-1 + 3 * y * y);
            derivy[8] = 5 * x * (1 + x) * y * y * (1 + y) + 0.5 * x * (1 + x) * y * (-3 + 5 * y * y) + 0.5 * x * (1 + x) * (1 + y) * (-3 + 5 * y * y);
            derivy[9] = 0.125 * x * (1 + x) * (1 + y) * (-60 * y + 140 * y * y * y) + 0.125 * x * (1 + x) * (3 - 30 * y * y + 35 * y * y * y * y);
            derivy[10] = 0.5 * (1 + x) * (-1 + 3 * x * x);
            derivy[11] = 0.5 * (1 + x) * (-1 + 3 * x * x) * y + 0.5 * (1 + x) * (-1 + 3 * x * x) * (1 + y);
            derivy[12] = 1.5 * (1 + x) * (-1 + 3 * x * x) * y * (1 + y) + 0.25 * (1 + x) * (-1 + 3 * x * x) * (-1 + 3 * y * y);
            derivy[13] = 2.5 * (1 + x) * (-1 + 3 * x * x) * y * y * (1 + y) + 0.25 * (1 + x) * (-1 + 3 * x * x) * y * (-3 + 5 * y * y) + 0.25 * (1 + x) * (-1 + 3 * x * x) * (1 + y) * (-3 + 5 * y * y);
            derivy[14] = 0.0625 * (1 + x) * (-1 + 3 * x * x) * (1 + y) * (-60 * y + 140 * y * y * y) + 0.0625 * (1 + x) * (-1 + 3 * x * x) * (3 - 30 * y * y + 35 * y * y * y * y);
            derivy[15] = 0.5 * x * (1 + x) * (-3 + 5 * x * x);
            derivy[16] = 0.5 * x * (1 + x) * (-3 + 5 * x * x) * y + 0.5 * x * (1 + x) * (-3 + 5 * x * x) * (1 + y);
            derivy[17] = 1.5 * x * (1 + x) * (-3 + 5 * x * x) * y * (1 + y) + 0.25 * x * (1 + x) * (-3 + 5 * x * x) * (-1 + 3 * y * y);
            derivy[18] = 2.5 * x * (1 + x) * (-3 + 5 * x * x) * y * y * (1 + y) + 0.25 * x * (1 + x) * (-3 + 5 * x * x) * y * (-3 + 5 * y * y) + 0.25 * x * (1 + x) * (-3 + 5 * x * x) * (1 + y) * (-3 + 5 * y * y);
            derivy[19] = 0.0625 * x * (1 + x) * (-3 + 5 * x * x) * (1 + y) * (-60 * y + 140 * y * y * y) + 0.0625 * x * (1 + x) * (-3 + 5 * x * x) * (3 - 30 * y * y + 35 * y * y * y * y);
            derivy[20] = 0.125 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x);
            derivy[21] = 0.125 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * y + 0.125 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * (1 + y);
            derivy[22] = 0.375 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * y * (1 + y) + 0.0625 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * (-1 + 3 * y * y);
            derivy[23] = 0.625 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * y * y * (1 + y) + 0.0625 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * y * (-3 + 5 * y * y) + 0.0625 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * (1 + y) * (-3 + 5 * y * y);
            derivy[24] = 0.015625 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * (1 + y) * (-60 * y + 140 * y * y * y) + 0.015625 * (1 + x) * (3 - 30 * x * x + 35 * x * x * x * x) * (3 - 30 * y * y + 35 * y * y * y * y);

    }
}


// Función para imprimir una matriz
void printMatrix(const Matrix& matrix) {
    for (const auto& row : matrix) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

// Función para crear una matriz identidad de tamaño n
Matrix createIdentityMatrix(int n) {
    Matrix identity(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        identity[i][i] = 1.0;
    }
    return identity;
}

// Función para realizar la eliminación Gaussiana para calcular la inversa
bool invertMatrix(const Matrix& A, Matrix& inverse) {
    int n = A.size();
    inverse = createIdentityMatrix(n);
    Matrix tempA = A; // Copiar A para no modificar la matriz original

    for (int i = 0; i < n; ++i) {
        // Encontrar el mayor valor en la columna i
        double maxEl = std::abs(tempA[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(tempA[k][i]) > maxEl) {
                maxEl = std::abs(tempA[k][i]);
                maxRow = k;
            }
        }

        // Intercambiar filas máximas con la fila actual en tempA e inverse
        for (int k = 0; k < n; ++k) {
            std::swap(tempA[maxRow][k], tempA[i][k]);
            std::swap(inverse[maxRow][k], inverse[i][k]);
        }

        // Hacer todas las filas debajo de esta cero en la columna actual
        double diagElement = tempA[i][i];
        if (std::abs(diagElement) < 1e-10) return false; // Singular matrix

        for (int k = 0; k < n; ++k) {
            tempA[i][k] /= diagElement;
            inverse[i][k] /= diagElement;
        }

        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double c = tempA[k][i];
                for (int j = 0; j < n; ++j) {
                    tempA[k][j] -= c * tempA[i][j];
                    inverse[k][j] -= c * inverse[i][j];
                }
            }
        }
    }
    return true;
}

// Función para multiplicar una matriz por un vector
std::vector<double> multiplyMatrixVector(const Matrix& matrix, const std::vector<double>& vec) {
    int n = matrix.size();
    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }

    return result;
}

// Función para calcular la inversa de A y multiplicarla por B
std::vector<double> solveSystem(const Matrix& A, const std::vector<double>& B) {
    Matrix inverse;
    if (!invertMatrix(A, inverse)) {
        throw std::runtime_error("La matriz es singular y no tiene inversa.");
    }
    return multiplyMatrixVector(inverse, B);
}
double exactSol(double x, double y){
  double sol =  (1 + x) *(1 + y) *atan(1 - x)* atan(1 - y);
    return sol;
}

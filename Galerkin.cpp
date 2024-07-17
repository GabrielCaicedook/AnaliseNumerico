//
//  ExactIntegration.cpp
//  AnaliseNumerico
//
//  Created by Gabriel Caicedo on 7/04/24.
//

#include "Galerkin.h"

ExactIntegration::ExactIntegration(){
}

double ExactIntegration::DefinedIntegral(double x){
    //Set the integral solve
    double y= (1/3)*pow(x, 3);
    
    return y;
}

double ExactIntegration::ValExact(double x1, double x2){
    double area = DefinedIntegral(x2)-DefinedIntegral(x1);
    return area;
    }

std::vector<double> ExactIntegration::ValExact(double a, double b, int ref){
    int nel = pow(2, ref);
    std::vector<double>vectorAreas(nel,0);
    double h = (b-a)/nel;
    for (int i=0; i<nel; i++){
        double point1 = a + i*h;
        double point2 = a + (i+1)*h;
        double area = ValExact(point1, point2);
        vectorAreas[i]=area;
    }
    return vectorAreas;
}


//CHECKEAR LÓGICA(ABAJO)
//Exact Integration Func

//Print solutions on .txt file


//    int count =1;
//    for (int i=0; i< TestExact.size(); i++) {
//        double val = TestExact[i];
//            archivo << i <<" "<< val << std::endl;
//            count++;
//        }
//
//    return;

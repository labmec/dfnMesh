//
//  InterpolPhil.cpp
//  Interpolator
//
//  Created by JOSE VILLEGAS on 22/5/18.
//  Copyright © 2018 JOSE VILLEGAS. All rights reserved.
//

#include "Interpol3RelPerm.h"

void Interpol3RelPerm::SetData(Matrix data, Kralpha alpa){
    switch (alpa) {
            
        case EKrw:
            Krw.SetData(data);
            break;
        case EKrow:
            Krow.SetData(data);
            break;
        case EKrg:
            Krg.SetData(data);
            break;
        case EKrog:
            Krog.SetData(data);
            break;
        default:
            break;
    }
}

double Interpol3RelPerm::Val(double Sw, double Sg){
 
    // calcula la saturacion inicial de petroleo
    int npointsw = Krw.GetData().Rows();
    int npointsg = Krw.GetData().Rows();
    double swc=0.0;
    for(int i=0; i<npointsw; i++){
        double sw_min = Krw.GetData().GetVal(i, 1);
        double sw_anal = Krw.GetData().GetVal(i+1, 1);
        if (sw_min != sw_anal) {
            swc = Krw.GetData().GetVal(i, 0);
            break;
        }
    }
    
    
   // double sorg = Krog.GetData().GetVal(npointsg-1, 0);
    
    //  Calcula las permeabilidades relativas para los sistemas krw-krow, krg-krog
    double krocw = Krow.Val(swc);
    double krw  = Krw.Val(Sw);
    double krow = Krow.Val(Sw);
    double krg  = Krg.Val(Sg);
    double krog = Krog.Val(Sg);
    
//    //Stone I
//    double alpha, bw, bg, Swe ;
//    alpha = 1.0 - (Sg / (1-swc - sorg));
//    
//    bw = (krow/krocw)/(1-Swe);
//    bg = (krow/krocw)/(1-Swe);
    
    
    
    
     // Stone II
    double val = (krocw)*((((krow/krocw)+krw)*((krog/krocw)+krg)) - (krw+krg));
    
    
    
    return val;
    
}
double Interpol3RelPerm::ValDeriv(double Sw, double Sg){
    
    int npoints = Krw.GetData().Rows();
    double swc=0.0;
    for(int i=0; i<npoints; i++){
        double sw_min = Krw.GetData().GetVal(i, 1);
        double sw_anal = Krw.GetData().GetVal(i+1, 1);
        if (sw_min != sw_anal) {
            swc = Krw.GetData().GetVal(i, 0);
            break;
        }
    }
    
    //calcula las derivadas
    std::pair<int, int> krocw = Krow.ValDeriv(swc);
    std::pair<int, int> krw  = Krw.ValDeriv(Sw);
    std::pair<int, int> krow = Krow.ValDeriv(Sw);
    std::pair<int, int> krg  = Krg.ValDeriv(Sg);
    std::pair<int, int> krog = Krog.ValDeriv(Sg);
    
    
    
}

void Interpol3RelPerm::ReadData(std::string data){
    std::ifstream file;
    file.open(data);
    
    int i=1;
    Matrix mdata;
    std::string line;
    int n_cols =0;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::istringstream issText(line);
        char l = line[0];
        int nchar = line.size();
        if(l != '/'){
            
            //calcula el numero de columnas
            if(n_cols==0){
                std::string w;
                for(int j=0; j< nchar; j++){
                    w="";
                    if (issText >> w);
                    if(w==""){
                        break;
                    }
                    n_cols ++;
                }
            }
            
            if (n_cols==2) {
                double a, b;
                mdata.Resize(i, 2);
                if(iss >> a >> b) ;
                mdata(i-1,0)=a;
                mdata(i-1,1)=b;
                i=i+1;
            }
            if (n_cols==3) {
                double a, b, c;
                mdata.Resize(i, 2);
                if(iss >> a >> b >> c) ;
                mdata(i-1,0)=a;
                mdata(i-1,1)=b;
                mdata(i-1,2)=c;
                i=i+1;
            }
        }
    }
    
    
    if(mdata.Rows()>0){
        std::cout<<"*************************"<<std::endl;
        std::cout<<"Reading file... ok!"<<std::endl;
        std::cout<<"*************************"<<std::endl;
       // SetData(mdata);
        mdata.Print(std::cout);
    }
    
}

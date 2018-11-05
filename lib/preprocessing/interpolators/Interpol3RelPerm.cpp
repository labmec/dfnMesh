//
//  InterpolPhil.cpp
//  Interpolator
//
//  Created by JOSE VILLEGAS on 22/5/18.
//  Copyright © 2018 JOSE VILLEGAS. All rights reserved.
//

#include "Interpol3RelPerm.h"



void Interpol3RelPerm::SetData(std::function<std::tuple<double, double>(double)> function, Kralpha alpa){
    switch (alpa) {
        case EKrw:
            Krw = function;
            break;
        case EKrow:
            Krow = function;
            break;
        case EKrg:
            Krg = function;
            break;
        case EKrog:
            Krog = function;
            break;
        default:
            break;
    }
}

double Interpol3RelPerm::Val(double Sw, double Sg){
    
    double krocw = std::get<0>(Krow(0.4)); //Saturacion inicial
    double krow = std::get<0>(Krow(Sw));
    double krw = std::get<0>(Krw(1.0 - Sw));
    double krg = std::get<0>(Krg(Sg));
    double krog = std::get<0>(Krog(1- Sg));//Saturacion inicial
    
    double intval =(krocw)*((((krow/krocw)+krw)*((krog/krocw)+krg)) - (krw + krg));
    
    return intval;
}

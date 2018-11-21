//
//  InterpolPhil.cpp
//  Interpolator
//
//  Created by JOSE VILLEGAS on 22/5/18.
//  Copyright © 2018 JOSE VILLEGAS. All rights reserved.

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

void Interpol3RelPerm::SetParam(REAL swc, REAL sorw, REAL sorg){
  
    REAL sum = swc + sorw + sorg;
    if(sum>1.0){DebugStop();}
    if(!(swc>0.0 && swc<1.0)){DebugStop();}
            fswc = swc;
    if(!(sorw>0.0 && sorw<1.0)){DebugStop();}
    fsorw = sorw;
    if(!(sorg>0.0 && sorg<1.0)){DebugStop();}
    fsorg = sorg;
    
}

void Interpol3RelPerm::SetKroModel(ModelInterpol model){
    fKroModel =model;
}

double Interpol3RelPerm::Val(double Sw, double Sg){
    
    int npointsw = Krw.GetData().Rows();
    int npointsg = Krg.GetData().Rows();
    
    if (fswc==-1){
        fswc = Krw.GetData().GetVal(0,0);
    }

    
    if(Sw<fswc){
        return 0.0;
    }
  
    double krocw = Krow.Val(fswc);
    double krw  = Krw.Val(Sw);
    double krow = Krow.Val(Sw);
    double krg  = Krg.Val(Sg);
    double krog = Krog.Val(Sg);
    
    //Stone I
    double alpha, bw, bg, So, Swe, Sge, Soe, Som ;
    
    //Normalizando las saturaciones:
    fsorg= 1.0 - (Krg.GetData().GetVal(npointsg-1, 0)) - fswc ;
    fsorw= 1.0 - Krw.GetData().GetVal(npointsw-1, 0);
    
    alpha = 1.0 - (Sg / (1- fswc - fsorg));
    
    Som = (alpha*fsorw) + ((1-alpha)*fsorg);
    Swe = (Sw - fswc)/(1.0 - fswc - Som);
    Sge = (Sg)/(1.0 - fswc - Som);
    So = 1.0 - Sw - Sg;
    Soe = (So - Som)/(1.0 - fswc - Som);

    if(So<Som){
        return 0.0;
    }
    
    bw = (krow/krocw)/(1.0 - Swe);
    bg = (krog/krocw)/(1.0 - Sge);
    
    
    double val =0.0;
    
    if (fKroModel==MStoneI) {
         val = krocw*Soe*bg*bw;
    }
    
    if (fKroModel==MStoneII) {
         val = (krocw)*((((krow/krocw)+krw)*((krog/krocw)+krg)) - (krw+krg));
    }

    
    if(val<0.0){
        val=0.0;
    }
    
    return val;

}
Matrix Interpol3RelPerm::Deriv(double Sw, double Sg){
    
    int npointsw = Krw.GetData().Rows();
    int npointsg = Krg.GetData().Rows();
    
    Matrix deriv(2,3,0.0);
    
    double valor = this->Val(Sw, Sg);
    if (valor==0){
        return deriv;
    }

    if (fswc==-1){
        fswc = Krw.GetData().GetVal(0,0);
    }
    if(Sw<fswc){
        return deriv;
    }
    
    double krocw = Krow.Val(fswc);
    
    std::pair<double, double> krw  = Krw.ValDeriv(Sw);
    std::pair<double, double> krow = Krow.ValDeriv(Sw);
    std::pair<double, double> krg  = Krg.ValDeriv(Sg);
    std::pair<double, double> krog = Krog.ValDeriv(Sg);

    double a, b, dKrwdSw, dKrgdSg, dKrodSw, dKrodSg;

    dKrwdSw = std::get<1>(Krw.ValDeriv(Sw));
    dKrgdSg = std::get<1>(Krg.ValDeriv(Sg));
    // calcula las derivadas con respecto a kro
    a = ((std::get<1>(krow))/krocw) + (std::get<1>(krw));
    b = ((std::get<0>(krog))/krocw) + (std::get<0>(krg));
    dKrodSw = krocw*((a*b) - (std::get<1>(krw)));
    
    a = ((std::get<0>(krow))/krocw) + (std::get<0>(krw));
    b = ((std::get<1>(krog))/krocw) + (std::get<1>(krg));
    dKrodSg = krocw*((a*b) - (std::get<1>(krg)));
    
    deriv(0,0)=dKrwdSw;
    deriv(0,1)=dKrodSw;
    deriv(1,1)=dKrodSg;
    deriv(1,2)=dKrgdSg;
    
    return deriv;
    
}

void Interpol3RelPerm::ReadData(std::string data, Kralpha alpa){
    switch (alpa) {
        case EKrw:
            Krw.ReadData(data);
            break;
        case EKrow:
            Krow.ReadData(data);
            break;
        case EKrg:
            Krg.ReadData(data);
            break;
        case EKrog:
            Krog.ReadData(data);
            break;
        default:
            break;
    }
}
std::function<double(double, double)> Interpol3RelPerm::GetFunction(){
    return [this](double Sw, double Sg){
        return this->Val(Sw,Sg);
    };
}
std::function<Matrix(double,double)> Interpol3RelPerm::GetFunctionDeriv(){
    return [this](double Sw, double Sg){
        return this->Deriv(Sw,Sg);
    };
}

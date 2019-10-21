// Derived class for building the asymptotic relation for radial p modes
// Created by Enrico Corsaro @ OACT - August 2019
// e-mail: emncorsaro@gmail.com
// Header file "NonRadialModesPatternModel.h"
// Implementations contained in "NonRadialModesPatternModel.cpp"


#ifndef NONRADIALMODESPATTERNMODEL_H
#define NONRADIALMODESPATTERNMODEL_H

#include <iostream>
#include "Model.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class NonRadialModesPatternModel : public Model
{
    public:
    
        NonRadialModesPatternModel(RefArrayXd const covariates, const int angularDegree, const double nuMax, const double DeltaNu, 
                                   const double inputEpsilon, const double inputAlpha, const double inputBeta); 
        ~NonRadialModesPatternModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);

    protected:


    private:

        int angularDegree;                  // The angular degree of the modes to fit
        double nuMax;                       // The frequency of maximum oscillation power
        double inputDeltaNu;                // The large frequency separation (here used as a fixed parameter)
        double inputEpsilon;                // The epsilon term of the asymptotic relation (here used as a fixed parameter)
        double inputAlpha;                  // The curvature term for DeltaNu
        double inputBeta;                   // The curvature term for the small frequency separation
}; 


#endif

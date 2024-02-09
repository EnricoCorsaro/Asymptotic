// Derived class for building the asymptotic relation for radial p modes
// Created by Enrico Corsaro @ OACT - August 2019
// e-mail: emncorsaro@gmail.com
// Header file "RadialModesPatternModel.h"
// Implementations contained in "RadialModesPatternModel.cpp"


#ifndef RADIALMODESPATTERNMODEL_H
#define RADIALMODESPATTERNMODEL_H

#include <iostream>
#include "Model.h"

using namespace std;
using Eigen::ArrayXd;
typedef Eigen::Ref<Eigen::ArrayXd> RefArrayXd;


class RadialModesPatternModel : public Model
{
    public:
    
        RadialModesPatternModel(RefArrayXd const covariates, const double nuMax, const double inputEpsilon, const double inputAlpha); 
        ~RadialModesPatternModel();

        virtual void predict(RefArrayXd predictions, RefArrayXd const modelParameters);
        virtual void computeVariance(RefArrayXd modelVariance, const RefArrayXd modelParameters){};

    protected:


    private:

        double nuMax;                       // The frequency of maximum oscillation power
        double inputEpsilon;                // The epsilon term of the asymptotic relation
        double inputAlpha;                  // The curvature term of the asymptotic relation

}; 


#endif

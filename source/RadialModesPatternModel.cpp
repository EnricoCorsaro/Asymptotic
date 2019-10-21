#include "RadialModesPatternModel.h"


// RadialModesPatternModel::RadialModesPatternModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//      nuMax:                  the frequency of maximum oscillation power (in microHz)
//      inputEpsilon:           the input value of the epsilon term of the asymptotic relation. If > 0, it is
//                              used as a fixed parameter.
//      inputAlpha:             the input value of the curvature term of the asymptotic relation. If > 0, it is
//                              used as a fixed parameter.
//

RadialModesPatternModel::RadialModesPatternModel(RefArrayXd const covariates, const double nuMax, 
                                                 const double inputEpsilon, const double inputAlpha)
: Model(covariates),
  nuMax(nuMax),
  inputEpsilon(inputEpsilon),
  inputAlpha(inputAlpha)
{
}










// RadialModesPatternModel::RadialModesPatternModel()
//
// PURPOSE: 
//      Destructor.
//

RadialModesPatternModel::~RadialModesPatternModel()
{
}










// RadialModesPatternModel::predict()
//
// PURPOSE:
//      Builds the predictions from a RadialModesPatternModel model to obtain
//      estimates of the asymptotic parameters for the radial modes asymptotic
//      pattern.
//
// INPUT:
//      predictions:        one-dimensional array to contain the predictions
//                          from the model
//      modelParameters:    one-dimensional array where each element
//                          contains the value of a free parameter of the model
//
// OUTPUT:
//      void
//
// NOTE:
//      The free parameters are to be given in the order
//      (i) DeltaNu
//      (ii) epsilon (phase term)
//      (iii) alpha (curvature term)
//

void RadialModesPatternModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    double DeltaNu = modelParameters(0);
    double epsilon;
    double alpha;
    
    if (inputEpsilon > 0)
    {
        epsilon = inputEpsilon;
    }
    else
    {
        epsilon = modelParameters(1);
    }
    
    if (inputAlpha != -99)
    {
        alpha = inputAlpha;
    }
    else
    {
        if (inputEpsilon > 0)
        {
            alpha = modelParameters(1);
        }
        else
        {
            alpha = modelParameters(2);
        }
    }

    predictions = DeltaNu*(covariates + epsilon + alpha/2.0*(covariates - nuMax/DeltaNu).square());
}

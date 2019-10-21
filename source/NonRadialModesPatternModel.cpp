#include "NonRadialModesPatternModel.h"


// NonRadialModesPatternModel::NonRadialModesPatternModel()
//
// PURPOSE: 
//      Constructor. Initializes model computation.
//
// INPUT:
//      covariates:             one-dimensional array containing the values
//                              of the independent variable.
//

NonRadialModesPatternModel::NonRadialModesPatternModel(RefArrayXd const covariates, const int angularDegree, const double nuMax, 
                                                       const double inputDeltaNu, const double inputEpsilon,
                                                       const double inputAlpha, const double inputBeta)
: Model(covariates),
  angularDegree(angularDegree),
  nuMax(nuMax),
  inputDeltaNu(inputDeltaNu),
  inputEpsilon(inputEpsilon),
  inputAlpha(inputAlpha),
  inputBeta(inputBeta)
{
}










// NonRadialModesPatternModel::NonRadialModesPatternModel()
//
// PURPOSE: 
//      Destructor.
//

NonRadialModesPatternModel::~NonRadialModesPatternModel()
{
}










// NonRadialModesPatternModel::predict()
//
// PURPOSE:
//      Builds the predictions from a NonRadialModesPatternModel model to obtain
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
//      (i) deltaNu0Degree (the small frequency separation for the given angular degree
//      (ii) alpha (curvature term for the large frequency separation)
//      (iii) beta (curvature term for small frequency separation)
//

void NonRadialModesPatternModel::predict(RefArrayXd predictions, RefArrayXd const modelParameters)
{
    double deltaNu0Degree = modelParameters(0);
    double alpha;
    double beta;

    if (inputAlpha != -99)
    {
        alpha = inputAlpha;
    }
    else
    {
        alpha = modelParameters(1);
    }

    if (inputBeta != -99)
    {
        beta = inputBeta;
    }
    else
    {
        if (inputAlpha != -99)
        {
            beta = modelParameters(1);
        }
        else
        {
            beta = modelParameters(2);
        }
    }

    predictions = inputDeltaNu*(covariates + inputEpsilon + angularDegree/2.0 + alpha/2.0*(covariates - nuMax/inputDeltaNu).square() - 
                  beta*(covariates - nuMax/inputDeltaNu)) - deltaNu0Degree;
}

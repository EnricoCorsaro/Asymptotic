// Main code for peak bagging by means of nested sampling analysis
// Created by Enrico Corsaro @ INAF-OACT - August 2019
// e-mail: emncorsaro@gmail.com
// Source code file "Asymptotic.cpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "Functions.h"
#include "File.h"
#include "MultiEllipsoidSampler.h"
#include "KmeansClusterer.h"
#include "EuclideanMetric.h"
#include "Prior.h"
#include "UniformPrior.h"
#include "NormalLikelihood.h"
#include "RadialModesPatternModel.h"
#include "NonRadialModesPatternModel.h"
#include "PowerlawReducer.h"
#include "Results.h"
#include "Ellipsoid.h"
#include "PrincipalComponentProjector.h"

int main(int argc, char *argv[])
{
    // Check number of arguments for main function
   
    if (argc != 8)
    {
        cerr << "Usage: ./asymptotic <Catalog ID> <Star ID> <output sub-directory> <run number> <input prior base filename> "
                "<input nuMax> <input angular degree>" << endl;
        exit(EXIT_FAILURE);
    }

    
    // ---------------------------
    // ----- Read input data -----
    // ---------------------------

    unsigned long Nrows;
    int Ncols;
    ArrayXXd data;
    string CatalogID(argv[1]);
    string StarID(argv[2]);
    string runNumber(argv[4]);
    string inputPriorBaseName(argv[5]);
    string inputNuMax(argv[6]);
    string inputDegree(argv[7]);
    double nuMax = stod(inputNuMax);
    int angularDegree = stoi(inputDegree);


    // Read the local path for the working session from an input ASCII file
    ifstream inputFile;
    File::openInputFile(inputFile, "localPath.txt");
    File::sniffFile(inputFile, Nrows, Ncols);
    vector<string> myLocalPath;
    myLocalPath = File::vectorStringFromFile(inputFile, Nrows);
    inputFile.close();


    // Set up some string paths used in the computation
    string outputSubDirName(argv[3]);
    string baseOutputDirName = myLocalPath[0] + "results/" + CatalogID + StarID + "/";
    string outputDirName = baseOutputDirName + outputSubDirName + "/";
    string outputPathPrefix = outputDirName + runNumber + "/asymptotic_";
    string baseInputDirName = baseOutputDirName + outputSubDirName + "/data/";
    string inputFileName = baseOutputDirName + outputSubDirName + "/data/" + runNumber + ".txt";

    cout << "------------------------------------------------------ " << endl;
    cout << " Performing asymptotic fit for l = " + inputDegree + " modes in " + CatalogID + StarID << endl;
    cout << "------------------------------------------------------ " << endl;


    // Read the input dataset
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nrows, Ncols);
    data = File::arrayXXdFromFile(inputFile, Nrows, Ncols);
    inputFile.close();


    // Creating frequency and PSD arrays
    ArrayXd covariates = data.col(0);
    ArrayXd observations = data.col(1);
    ArrayXd uncertainties = data.col(2);


    // -------------------------------------------------------
    // ----- First step. Set up all prior distributions ------
    // -------------------------------------------------------
    
    unsigned long Nparameters;              // Number of parameters for which prior distributions are defined

    // ---- Read prior hyper parameters for resolved modes -----
    inputFileName = outputDirName + inputPriorBaseName + "_" + runNumber + ".txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);
    ArrayXXd hyperParameters;
  
    if (Ncols == 1)
    {
        Ncols = 3;
        hyperParameters.conservativeResize(Nparameters,Ncols);
    }

    hyperParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();
    
    ArrayXd hyperParametersMinima = hyperParameters.col(0);
    ArrayXd hyperParametersMaxima = hyperParameters.col(1);
    // ---------------------------------------------------------

    int Ndimensions = Nparameters;              // Total number of dimensions of the peak bagging model

    if (angularDegree == 0)
    {
        if (Nparameters != 3)
        {
            cerr << "Wrong number of input prior hyper-parameters." << endl;
            cerr << "When performing an asymptotic fit for radial modes, three lines are " << endl;
            cerr << "expected from the input prior list (1 - DeltaNu, 2 - epsilon, 3 - alpha)." << endl;
        }
    }
    else
    {
        if (Nparameters != 5)
        {
            cerr << "Wrong number of input prior hyper-parameters." << endl;
            cerr << "When performing an asymptotic fit for radial modes, five lines are " << endl;
            cerr << "expected from the input prior list (1 - DeltaNu, 2 - epsilon, 3 - deltaNu0Degree, 4 - alpha, 5 - beta)." << endl;
        }
    }


    // Uniform Prior
    
    int NpriorTypes = 1;                                        // Total number of prior types included in the computation
    vector<Prior*> ptrPriors(NpriorTypes);
    
    double DeltaNu = 0;
    double epsilon = 0;
    double alpha = -99;
    double beta = -99;
    ArrayXd parametersMinima;                      // Minima for prior PDF
    ArrayXd parametersMaxima;                      // Maxima for prior PDF

    if (angularDegree == 0)
    {
        parametersMinima.resize(Ndimensions);
        parametersMaxima.resize(Ndimensions);
        parametersMinima << hyperParametersMinima;
        parametersMaxima << hyperParametersMaxima;
        
        if (hyperParametersMinima(1) == hyperParametersMaxima(1))
        {
            // In this case epsilon is a fixed value. Then reduce the number
            // of free parameters by 1 and remove the epsilon prior line from
            // the list of input prior parameters (only Deltanu and alpha are left).

            Ndimensions--;
            parametersMinima.resize(Ndimensions);
            parametersMaxima.resize(Ndimensions);
            parametersMinima << hyperParametersMinima(0), hyperParametersMinima(2);
            parametersMaxima << hyperParametersMaxima(0), hyperParametersMaxima(2); 
            epsilon = hyperParametersMinima(1);
        }

        if (hyperParametersMinima(2) == hyperParametersMaxima(2))
        {
            // In this case alpha is a fixed value. Then reduce the number
            // of free parameters by 1 and remove the alpha prior line from
            // the list of input prior parameters (only DeltaNu and epsilon, or only DeltaNu are left).

            Ndimensions--;
            parametersMinima.conservativeResize(Ndimensions);
            parametersMaxima.conservativeResize(Ndimensions);
            alpha = hyperParametersMinima(2);
        }
    }

    if (angularDegree == 1 || angularDegree == 2 || angularDegree == 3)
    {
        // Here DeltaNu, epsilon are fixed parameters. Then reduce the
        // number of free parameters by 2 and remove DeltaNu and epsilon
        // prior lines from the list of input prior parameters.

        Ndimensions = 3;
        parametersMinima.resize(Ndimensions);
        parametersMaxima.resize(Ndimensions);
        parametersMinima << hyperParametersMinima.segment(2,Ndimensions);
        parametersMaxima << hyperParametersMaxima.segment(2,Ndimensions); 
        DeltaNu = hyperParametersMinima(0);
        epsilon = hyperParametersMinima(1);


        // Check if alpha is fixed or not. If fixed, remove it from the prior
        // list and set its value to the one given in the prior file.

        if (hyperParametersMinima(3) == hyperParametersMaxima(3))
        {
            Ndimensions--;
            parametersMinima.resize(Ndimensions);
            parametersMaxima.resize(Ndimensions);
            parametersMinima << hyperParametersMinima(2), hyperParametersMinima(4);
            parametersMaxima << hyperParametersMaxima(2), hyperParametersMaxima(4);
            alpha = hyperParametersMinima(3);
        }


        // Check if beta is fixed or not. If fixed, remove it from the prior
        // list and set its value to the one given in the prior file.

        if (hyperParametersMinima(4) == hyperParametersMaxima(4))
        {
            Ndimensions--;
            parametersMinima.conservativeResize(Ndimensions);
            parametersMaxima.conservativeResize(Ndimensions);
            beta = hyperParametersMinima(4);
        }
    }

    UniformPrior uniformPrior(parametersMinima, parametersMaxima);
    ptrPriors[0] = &uniformPrior;
    
    string fullPathHyperParameters = outputPathPrefix + "hyperParametersUniform.txt";
    uniformPrior.writeHyperParametersToFile(fullPathHyperParameters);

    
    // -------------------------------------------------------------------
    // ---- Second step. Set up the asymptotic model ---------------------
    // -------------------------------------------------------------------
    
    Model *model = nullptr;

    if (angularDegree == 0)
    {
        model = new RadialModesPatternModel(covariates, nuMax, epsilon, alpha);
    }
   
    if (angularDegree > 0)
    {
        model = new NonRadialModesPatternModel(covariates, angularDegree, nuMax, DeltaNu, epsilon, alpha, beta);
    }
    

    // -----------------------------------------------------------------
    // ---- Third step. Set up the likelihood function to be used ------
    // -----------------------------------------------------------------
    
    NormalLikelihood likelihood(observations, uncertainties, *model);
    

    // -------------------------------------------------------------------------------
    // ----- Fourth step. Set up the X-means clusterer using an Euclidean metric -----
    // -------------------------------------------------------------------------------

    inputFileName = outputDirName + "Xmeans_configuringParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);

    if (Nparameters != 2)
    {
        cerr << "Wrong number of input parameters for X-means algorithm." << endl;
        exit(EXIT_FAILURE);
    }

    ArrayXd configuringParameters;
    configuringParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();
    
    int minNclusters = configuringParameters(0);
    int maxNclusters = configuringParameters(1);
    
    if ((minNclusters <= 0) || (maxNclusters <= 0) || (maxNclusters < minNclusters))
    {
        cerr << "Minimum or maximum number of clusters cannot be <= 0, and " << endl;
        cerr << "minimum number of clusters cannot be larger than maximum number of clusters." << endl;
        exit(EXIT_FAILURE);
    }
    
    int Ntrials = 10;
    double relTolerance = 0.01;     // k-means
  
    bool printNdimensions = false;
    PrincipalComponentProjector projector(printNdimensions);
    bool featureProjectionActivated = false;
    EuclideanMetric myMetric;
    
    KmeansClusterer clusterer(myMetric, projector, featureProjectionActivated, 
                           minNclusters, maxNclusters, Ntrials, relTolerance); 


    // ---------------------------------------------------------------------
    // ----- Fifth step. Configure and start nested sampling inference -----
    // ---------------------------------------------------------------------

    inputFileName = outputDirName + "NSMC_configuringParameters.txt";
    File::openInputFile(inputFile, inputFileName);
    File::sniffFile(inputFile, Nparameters, Ncols);
    configuringParameters.setZero();
    configuringParameters = File::arrayXXdFromFile(inputFile, Nparameters, Ncols);
    inputFile.close();

    if (Nparameters > 9 || Nparameters < 8)
    {
        cerr << "Wrong number of input parameters for NSMC algorithm." << endl;
        cerr << "There must be either 8 or 9 parameters. In this case the last parameter is always ignored." << endl; 
        exit(EXIT_FAILURE);
    }
    

    // Print results on the screen
    
    bool printOnTheScreen = true;                  


    // Initial number of live points
    
    int initialNobjects = configuringParameters(0);

    
    // Minimum number of live points 
    
    int minNobjects = configuringParameters(1);
    
    
    // Maximum number of attempts when trying to draw a new sampling point
    
    int maxNdrawAttempts = configuringParameters(2);

    
    // The first N iterations, we assume that there is only 1 cluster
    
    int NinitialIterationsWithoutClustering = configuringParameters(3);

    
    // Clustering is only happening every N iterations.
    
    int NiterationsWithSameClustering = configuringParameters(4);

    
    // Fraction by which each axis in an ellipsoid has to be enlarged
    // It can be a number >= 0, where 0 means no enlargement. configuringParameters(5)
    // Calibration from Corsaro et al. (2018)
    
    double initialEnlargementFraction = 0.369*pow(Ndimensions,0.574);    

    
    // Exponent for remaining prior mass in ellipsoid enlargement fraction.
    // It is a number between 0 and 1. The smaller the slower the shrinkage of the ellipsoids.
    
    double shrinkingRate = configuringParameters(6);        
                                                                                                                    
    
    // Termination factor for nested sampling process.                                 
    
    double terminationFactor = configuringParameters(7);

    
    // Total maximum number of nested iterations required to carry out the computation.
    // This is used only in the multi-modal approach.
    
    int maxNiterations = 0; 
    

    MultiEllipsoidSampler nestedSampler(printOnTheScreen, ptrPriors, likelihood, myMetric, clusterer, 
                                        initialNobjects, minNobjects, initialEnlargementFraction, shrinkingRate);
    
    double tolerance = 1.e2;
    double exponent = 0.4;
    PowerlawReducer livePointsReducer(nestedSampler, tolerance, exponent, terminationFactor);

    nestedSampler.run(livePointsReducer, NinitialIterationsWithoutClustering, NiterationsWithSameClustering, 
                      maxNdrawAttempts, terminationFactor, maxNiterations, outputPathPrefix);

    nestedSampler.outputFile << "# List of configuring parameters used for the ellipsoidal sampler and X-means" << endl;
    nestedSampler.outputFile << "# Row #1: Minimum Nclusters" << endl;
    nestedSampler.outputFile << "# Row #2: Maximum Nclusters" << endl;
    nestedSampler.outputFile << "# Row #3: Initial Enlargement Fraction" << endl;
    nestedSampler.outputFile << "# Row #4: Shrinking Rate" << endl;
    nestedSampler.outputFile << minNclusters << endl;
    nestedSampler.outputFile << maxNclusters << endl;
    nestedSampler.outputFile << initialEnlargementFraction << endl;
    nestedSampler.outputFile << shrinkingRate << endl;
    nestedSampler.outputFile << "# Other information on the run" << endl;
    nestedSampler.outputFile << "# Row #1: Local working path used" << endl;
    nestedSampler.outputFile << "# Row #2: Catalog and Star ID" << endl;
    nestedSampler.outputFile << "# Row #3: Run Directory" << endl;
    nestedSampler.outputFile << "# Row #4: Run Number" << endl;
    nestedSampler.outputFile << "# Row #5: nuMax (microHz)" << endl;
    nestedSampler.outputFile << "# Row #7: angular degree (either 0, 1, 2, or 3)" << endl;
    nestedSampler.outputFile << myLocalPath[0] << endl;
    nestedSampler.outputFile << CatalogID + StarID << endl;
    nestedSampler.outputFile << outputSubDirName << endl;
    nestedSampler.outputFile << runNumber << endl;
    nestedSampler.outputFile << nuMax << endl;
    nestedSampler.outputFile << angularDegree << endl;
    nestedSampler.outputFile.close();


    // -------------------------------------------------------
    // ----- Last step. Save the results in output files -----
    // -------------------------------------------------------
   
    Results results(nestedSampler);
    results.writeParametersToFile("parameter");
    results.writeLogLikelihoodToFile("logLikelihood.txt");
    results.writePosteriorProbabilityToFile("posteriorDistribution.txt");
    results.writeLogEvidenceToFile("logEvidence.txt");
    results.writeLogMeanLiveEvidenceToFile("logMeanLiveEvidence.txt");
    results.writeEvidenceInformationToFile("evidenceInformation.txt");
    

    // Print out parameter estimates only in the case of a uni-modal high-dimensional fit.

    double credibleLevel = 68.3;
    bool writeMarginalDistributionToFile = true;
    results.writeParametersSummaryToFile("parameterSummary.txt", credibleLevel, writeMarginalDistributionToFile);

    cout << "Process # " << runNumber << " under subdir: " + outputSubDirName + " has been completed." << endl;

    return EXIT_SUCCESS;
}

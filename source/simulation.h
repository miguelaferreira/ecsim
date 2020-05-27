#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <limits>
#include <vector>

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/Core>

using namespace std;

// minimal step sizes (except for step=0)
const double MIN_STEP_TIME = 1.0e-9; // step time > 1 ns seems a reasonable minimal time step (as in: no time step will ever be smaller than this)
const double MIN_STEP_POT = 1.0e-6; // step pot > 1 uV seems reasonable (as in: no potential step will ever be smaller than this)
const double MIN_RATE = 1.0e-6; // minimal rate > 1e-6 /s seems reasonable (as in: no rate will ever be smaller than this)
const double MIN_CONC = 1.0e-9; // minimal (initial!) concentration > 1 pM (because /1000) seems reasonable
const double MIN_AREA = 1.0e-18; // minimal electrode area > 1 nm2 seems reasonable

const double CONST_PI = 3.14159265358979323846;
const double CONST_F = 96485.33212;
const double CONST_R = 8.31446261815324;

#include "coefs_alpha_beta.h"
#include "electrodes.h"
#include "system.h"
#include "environment.h"
#include "experiment.h"
#include "matrixsystem.h"

uint64_t getPosixClockTime();

class Simulation
{
private:
    ostream &output_stream;
    MatrixSystem msys;
    
    double f; // F/RT
    double totalTime, totalTheta, deltaTheta = 0.2, sigma;

    double maxX, deltaX; // maxX = maximum distance from electrode [1], deltaX = size of smallest grid element [1]
    double maxT, deltaT; // maxT = total experiment time [1], deltaT = size of smallest time step [1]
    double paramGamma = 1.67, paramR0; // paramGamma = grid expansion factor, paramR0 is some electrode-geometry-dependent parameter
    double minF = 2.2, maxF = 6.2;
    double minLogRate = 3.0, maxLogRate = 7.0;
    double paramLambda; // dimensionless kinetic parameter of the system

    size_t numSpecies, numRedox, numReactions;
    size_t numGridPoints, numCurrentPoints = 5, numDerivPoints = 6; // number of points in grid, number of points for calculation of current, number of points for diffusion

    vector<double> gridCoordinate, paramGammai, paramGamma2i; // grid coordinates and expansion factor (and its ^2), pre-calculated for all points in grid
    vector<double> coeffAlpha, coeffBeta, coeffBeta0; // coefficients for diffusion (alpha and beta) and current (beta0)
    Eigen::MatrixXd independentTerms; // "independent terms" equals b, the right-hand side in Ax=b when solving for x [species_index, grid_location]
    Eigen::MatrixXd gridConcentration, gridConcentrationPrevious; // current concentrations in grid for all species [species_index, grid_location]
    
    // matrix that stores the 2nd order kinetic interactions (needed for Laasonen linearization)
    // secondOrderKinetics[] = {A, B, C, normalized_rate} for the 2nd order reaction A+B<->C
    vector<tuple<size_t, size_t, size_t, double>> secondOrderKinetics;
    vector<tuple<size_t, size_t, double>> firstOrderKinetics;
    vector<bool> speciesInRedox;
    vector<double> matrixB0DiagonalTerms;
    Eigen::MatrixXd currentContributionMatrix;
    Eigen::VectorXd currentContributionSpeciesFlux, currentContributionRedoxFlux;

public:
    double ipc, ipa, Epc, Epa; // for electroanalytical purposes (set during run())

    Environment env;
    Electrode el;
    Experiment exper;
    System sys;
    
    Simulation(ostream &_output) : output_stream( _output ) {}
    ~Simulation() {}

    void setGridSizing(double _gamma, double _minF, double _maxF, double _minlograte, double _maxlograte)
        { paramGamma = _gamma; minF = _minF; maxF = _maxF; minLogRate = _minlograte; maxLogRate = _maxlograte; }
    void setPotentialSizing(double _deltatheta) { deltaTheta = _deltatheta; }
    void setDifferentialOrders(size_t _numcurr, size_t _numderiv) { numCurrentPoints = _numcurr; numDerivPoints = _numderiv; }

    size_t run(vector<double>&, vector<double>&);
    
private:
    size_t runSimulation(vector<double>&, vector<double>&);
    void scanSegment(double, double, bool, vector<double>&, vector<double>&);
    void delaySegment(double, double);
    void solveSystem(double);
    void storeConcentrations();
    void invertMatrix();
    void updateIndependentTerms();
    void updateRedoxInMatrix(double);
    void updateKineticsInMatrix(bool);
    double calcCurrentFromFlux();

    void prepareSimulation();
    void initParametersEtc();
    void initVectorsEtc();
    void addRedoxToMatrix();
    void addBICoeffsToMatrix();
    void addKineticsToMatrix();
    void addHalfReactionKineticsToMatrix(Species*, Species*, Species*, Species*, double);
    void addFirstOrderKineticTerm(Species*, Species*, double);
    void addSecondOrderKineticTerm(Species*, Species*, Species*, double);

    double coeffMatrixN2(size_t, int, double);
};

#endif // SIMULATION_H

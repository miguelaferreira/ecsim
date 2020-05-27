#include "matrixsystem.h"

void MatrixSystem::initialize(size_t ns, size_t ng)
{
    numSpecies = ns;
    numGridPoints = ng;
    size_t numNonZeroesTotal = 15*numSpecies*numGridPoints;  // estimation (15 > 2*numDerivPoints)
    
    vecb = Eigen::VectorXd::Zero(numSpecies*numGridPoints);
    vecx = Eigen::VectorXd::Zero(numSpecies*numGridPoints);
    
    // initialize Triplet container
    tripletContainerMatrixA.clear();
    tripletContainerMatrixA.reserve(numNonZeroesTotal);

    // clear and resize matrix:
    matrixA.setZero(); // removes all non-zeros
    matrixA.resize(0,0);
    matrixA.data().squeeze(); // removes all items from matrix
    matrixA.resize(static_cast<Eigen::Index>(numSpecies*numGridPoints), static_cast<Eigen::Index>(numSpecies*numGridPoints));
    matrixA.reserve(static_cast<Eigen::Index>(numNonZeroesTotal));
    
    matrixPatternAnalyzed = false;
}

void MatrixSystem::solve(Eigen::MatrixXd *independentTerms, Eigen::MatrixXd *gridConcentration)
{
    // pattern needs to be analyzed once after matrix changed:
    if (!matrixPatternAnalyzed) { sparseMatrixSolver.analyzePattern(matrixA); matrixPatternAnalyzed = true; }
    
    sparseMatrixSolver.factorize(matrixA);

    for (size_t s = 0; s < numSpecies; s++)
        for (size_t x = 0; x < numGridPoints; x++)
            vecb[static_cast<Eigen::Index>( s*numGridPoints+x )] = (*independentTerms)(s, x);
    
    // this is where the magic happens:
    vecx = sparseMatrixSolver.solve(vecb);
    
    for (size_t s = 0; s < numSpecies; s++)
        for (size_t x = 0; x < numGridPoints; x++)
            (*gridConcentration)(s, x) = vecx[static_cast<Eigen::Index>( s*numGridPoints+x )];
}

//void MatrixSystem::addToMatrix(size_t row, size_t col, double value)
void MatrixSystem::addToMatrix(size_t s1, size_t g1, size_t s2, size_t g2, double value)
{
    tripletContainerMatrixA.emplace_back( Triplet(static_cast<int>(s1*numGridPoints+g1),
                                                  static_cast<int>(s2*numGridPoints+g2), value) );
}

void MatrixSystem::createMatrix()
{
    // set initial values in matrix:
    matrixA.setFromTriplets(tripletContainerMatrixA.begin(), tripletContainerMatrixA.end());
}

void MatrixSystem::changeValue(size_t s1, size_t g1, size_t s2, size_t g2, double value)
{
    matrixA.coeffRef(static_cast<int>(s1*numGridPoints+g1),
                     static_cast<int>(s2*numGridPoints+g2)) = value;
}

void MatrixSystem::addToValue(size_t s1, size_t g1, size_t s2, size_t g2, double value)
{
    matrixA.coeffRef(static_cast<int>(s1*numGridPoints+g1),
                     static_cast<int>(s2*numGridPoints+g2)) += value;
}
 

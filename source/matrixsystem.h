#ifndef MATRIXSYSTEM_H
#define MATRIXSYSTEM_H

#include <vector>

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/Core>

using namespace std;

using Triplet = Eigen::Triplet<double>; //triplets (row, column, value) to fill sparse matrix
using Solver  = Eigen::SparseLU<Eigen::SparseMatrix<double>>;  //Turns out that SparseLU is actually fast enough, but only when running as Release :)

class MatrixSystem
{
private:
    vector<Triplet> tripletContainerMatrixA; // this container holds (row, col, value) triplets before matrixA is assembled
    Eigen::SparseMatrix<double> matrixA; // matrix to be inverted: Ax=b corresponds to matrixA * vecx = vecb
    bool matrixPatternAnalyzed;
    Solver sparseMatrixSolver; // solver that inverts the sparse matrixA using LU factorization
    Eigen::VectorXd vecb, vecx; // vectors, defininition see above
    
    size_t numSpecies, numGridPoints; // heavily used, as this is the number of species * the number of grid points (effectively, to store in a matrix in a vector)
public:
    void initialize(size_t, size_t);
    void solve(Eigen::MatrixXd*, Eigen::MatrixXd*);
    void addToMatrix(size_t, size_t, size_t, size_t, double); // add to tripletcontainer
    void createMatrix(); // create matrix from tripletcontainer
    void changeValue(size_t, size_t, size_t, size_t, double); // change directly in matrix:   matrix.coeffref(x,y) = val
    void addToValue(size_t, size_t, size_t, size_t, double); // sum value directly in matrix: matrix.coeffref(x,y) += val
};

#endif // MATRIXSYSTEM_H

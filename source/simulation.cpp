#include <iostream>
#include "simulation.h"

uint64_t getPosixClockTime()
{
    struct timespec ts = {};
    return (clock_gettime (CLOCK_MONOTONIC, &ts) == 0) ? static_cast<uint64_t>(ts.tv_sec * 1e9 + ts.tv_nsec) : 0;
}

/*===============================================================================================
 * HIGHEST LEVEL CLASS METHODS
 *=============================================================================================*/

size_t Simulation::run(vector<double> &current, vector<double> &potential)
{
    // prepare simulation:
    prepareSimulation();

    // run simulation
    size_t numDataPoints = runSimulation(current, potential);
    
    return numDataPoints;
}

size_t Simulation::runSimulation(vector<double> &current, vector<double> &potential)
{
    // initialize electroanalytical parameters:
    ipc = 0.0;
    ipa = 0.0;
    Epc = 0.0;
    Epa = 0.0;

    uint64_t startTime, endTime; // timing variables [ns]
    startTime = getPosixClockTime();
    
    // conditioning step:
    delaySegment(exper.conditioningPotential, exper.conditioningTime);
    
    // equilibration step:
    delaySegment(exper.initialPotential, exper.equilibrationTime);
    
    // scan cycles:
    double segmentStartPotential; // start potential of scan segment
    bool recordCurrent; // record current in segment?
    for (int cycle = 0; cycle < exper.numCycles; cycle++) {
        recordCurrent = (cycle == exper.numCycles - 1); // warning: hard-coded option to record current only of last cycle

        segmentStartPotential = exper.initialPotential;
        for (auto vp: exper.vertexPotentials) {
            // scan up to vertex:
            scanSegment(segmentStartPotential, vp, recordCurrent, current, potential);
            segmentStartPotential = vp;
            // vertex delay:
            delaySegment(vp, exper.vertexDelay);
        }
        // scan from last vertex to final potential:
        scanSegment(segmentStartPotential, exper.finalPotential, recordCurrent, current, potential);
    }
    
    endTime = getPosixClockTime();
    output_stream << "Est. simulation time (" << current.size() << " steps &amp; " << numGridPoints << " grid points) = " << static_cast<double>(endTime - startTime)/1.0e9 << " s<br />" << endl;

    return current.size();
}

void Simulation::delaySegment(double potential, double time)
{
    if (time > MIN_STEP_TIME)
    {
        int num_points = static_cast<int>( (time/ (totalTime*deltaTheta/totalTheta) ));
        
        for (int d = 0; d < num_points; d++)
            solveSystem(f*potential);
    }
}

void Simulation::scanSegment(double potential_start, double potential_stop, bool record_current, vector<double> &current, vector<double> &potential)
{
    double curr;
    double theta_start = f*potential_start;
    double theta_stop = f*potential_stop;

    double sign = 0.0;
    if (theta_start < theta_stop) sign = 1.0;
    else if (theta_start > theta_stop) sign = -1.0;

    int num_points = static_cast<int>(abs(theta_start-theta_stop)/deltaTheta);
    double theta = theta_start;
    for (int d = 0; d < num_points; d++)
    {
        // do magic:
        solveSystem(theta);

        if (record_current) {
            // get current:
            curr = calcCurrentFromFlux();

            // update electroanalytical metrics:
            if (curr > ipa) { ipa = curr; Epa = theta/f; }
            if (curr < ipc) { ipc = curr; Epc = theta/f; }

            // add data points to current & potential vectors:
            potential.push_back(theta/f);
            current.push_back(curr);
        }

        // increase potential:
        theta += sign*deltaTheta;
    }
}

/*===============================================================================================
 * THIS IS THE SIMULATION
 *=============================================================================================*/

void Simulation::solveSystem(double theta)
{
    storeConcentrations();
    updateIndependentTerms();
    updateRedoxInMatrix(theta); // update potentials in matrix
    updateKineticsInMatrix(true);
    msys.solve(&independentTerms, &gridConcentration);
    updateKineticsInMatrix(false);
}

void Simulation::storeConcentrations()
{
    gridConcentrationPrevious = gridConcentration;
}

void Simulation::updateKineticsInMatrix(bool add)
{
    double rate;
    size_t spec1, spec2, spec3;

    for (auto kinTerm: secondOrderKinetics)
    {
        spec1 = get<0>(kinTerm);
        spec2 = get<1>(kinTerm);
        spec3 = get<2>(kinTerm);

        for (size_t x = 1; x < numGridPoints; x++)
        {
            rate = get<3>(kinTerm) * deltaX*deltaX*paramGamma2i[x];
            
            if (add)
            {
                msys.addToValue( spec3, x, spec1, x, rate * gridConcentration(spec2, x) );
                msys.addToValue( spec3, x, spec2, x, rate * gridConcentration(spec1, x) );
                
                independentTerms(spec3, x) += rate * gridConcentration(spec1, x) * gridConcentration(spec2, x);
            }
            else
            {
                msys.addToValue( spec3, x, spec1, x, -rate * gridConcentrationPrevious(spec2, x) );
                msys.addToValue( spec3, x, spec2, x, -rate * gridConcentrationPrevious(spec1, x) );
            }
        }
    }
}

void Simulation::updateIndependentTerms()
{
    size_t EndNormal = numGridPoints - numDerivPoints + 2;

    //for(int s = 0; s < sys.num_species; s++)
    for (auto spec: sys.vecSpecies)
    {
        //Zero at matrix rows storing surface conditions
        independentTerms(spec->getIndex(), 0) = 0.0;

        for (size_t x = 1; x < numGridPoints; x++)
        {
            /* independentTerms is the RHS ('b') of the matrix equation A*x=b, where 'x' is the concentration at T+dT (the 'next' concentration)
             * the system of equations is ('a' are coefficients, 'x(i,s)' the next concentrations at grid point i for species s):
             * a(-1)*x(i-1,s) + a(0)*x(i,s) + a(+1)*x(i+1,s) + a(+2)*x(i+2,s) + (etc.) = b(i,s)
             */
            independentTerms(spec->getIndex(), x) = -gridConcentration(spec->getIndex(), x) * paramGamma2i[x]/paramLambda;
            // however, when i>numGridPoints, then x(i+i) reduces to the (never changing) bulk concentration and we add them into 'b', because they are now known:
            if (x >= EndNormal)
            {
                for (int jj = 0; jj <= static_cast<int>(x - EndNormal); jj++)
                {
                    independentTerms(spec->getIndex(), x) -= coeffMatrixN2(x, static_cast<int>(numDerivPoints)-jj-2, spec->getDiffNorm()) * spec->getConcNormEquil();
                }
            }
        }
    }
}

void Simulation::updateRedoxInMatrix(double theta)
{
    double p, Kred, Kox;

    // reset diagonal terms:
    fill(matrixB0DiagonalTerms.begin(), matrixB0DiagonalTerms.end(), 0.0);

    for (auto redox: sys.vecRedox)
    {
        p = static_cast<double>(redox->getNe()) * (theta - f * redox->getE0()); // normalized potential

        Kred = deltaX * redox->getKeNorm() * exp(-redox->getAlpha() * p); // B-V kinetics
        Kox = deltaX * redox->getKeNorm() * exp((1.0 - redox->getAlpha()) * p); // B-V kinetics

        // add terms to matrix per redox reaction:
        msys.changeValue( redox->getSpecOx()->getIndex(), 0,
                          redox->getSpecRed()->getIndex(), 0,
                          -Kox / redox->getSpecOx()->getDiffNorm() );   // B0
        msys.changeValue( redox->getSpecRed()->getIndex(), 0,
                          redox->getSpecOx()->getIndex(), 0, 
                          -Kred / redox->getSpecRed()->getDiffNorm() ); // B0

        // add to diagonal terms:
        matrixB0DiagonalTerms[redox->getSpecOx()->getIndex()] += Kred / redox->getSpecOx()->getDiffNorm();
        matrixB0DiagonalTerms[redox->getSpecRed()->getIndex()] += Kox / redox->getSpecRed()->getDiffNorm();
    }

    // add diagonal terms to matrix:
    for (size_t s = 0; s < numSpecies; s++)
    {
        if (speciesInRedox[s]) // active redox species
        {
            msys.changeValue( s, 0, s, 0,  1.0 + matrixB0DiagonalTerms[s] ); // B0
            msys.changeValue( s, 0, s, 1, -1.0 ); // B1
        }
    }
}

double Simulation::calcCurrentFromFlux()
{
    double totalflux, speciesflux, currentFromFlux;
    currentFromFlux = el.epsilon * CONST_F * sys.getDiffMax() * sys.getConcMax();

    // determine flux per species:
    for (auto spec: sys.vecSpecies)
    {
        speciesflux = 0;
        for(size_t kk = 0; kk < numCurrentPoints; kk++)
        {
            speciesflux += coeffBeta0[kk] * gridConcentration(spec->getIndex(), kk);
        }
        currentContributionSpeciesFlux[static_cast<Eigen::Index>(spec->getIndex())] = speciesflux * spec->getDiffNorm();
    }

    /* SOLVE FOR THE FLUX PER REDOX REACTION BY LEAST SQUARES
     * - The method is described in Compton (pages 117 through 119) but is not general
     * - I have generalized the method as follows, and it works for all mechanisms I have tried:
     *   - The nett flux per redox reaction flux_redox = k_red*[ox] - k_ox*[red]
     *   - The total flux at the electrode (from which the current can be calculated) can be expressed as:
     *     total_flux = - SUM(flux_redox, for each redox reaction)
     *   - The flux of each species at the electrode is:
     *     +redox_flux if the species is ox in that redox reaction
     *     -redox_flux if the species is red in that redox reaction
     *     0 if it doesn't participate in that redox reaction
     *     --> I have put all these "current contributions" in currentContributionMatrix(species, redox) where each value is -1.0/0.0/1.0
     *   - The flux of each species at the electrode is ALSO:
     *     - Calculated by the "current function", using beta0 values (see for loop above)
     *     --> These species fluxes are stored in currentContributionSpeciesFlux(species)
     *   - The nett flux per redox reaction is then determined by solving the (overdetermined!) matrix equation:
     *     F*fr=fs where F = currentContributionMatrix, fr = flux per redox reaction, fs = species flux
     *   - The total flux is then determined by summing the elements in fr
     *   :)
     */
    // determine nett flux per redox reaction
    currentContributionRedoxFlux = (currentContributionMatrix.transpose() * currentContributionMatrix).ldlt().solve(
                currentContributionMatrix.transpose() * currentContributionSpeciesFlux );

    // sum to get total flux:
    totalflux = 0.0;
    for (auto redox: sys.vecRedox)
    {
        // this doesn't work well for n > 1 --> why not?
        totalflux -= currentContributionRedoxFlux[static_cast<Eigen::Index>(redox->getIndex())]
                     * static_cast<double>(redox->getNe());
    }
    totalflux /= deltaX; // totalflux := dC/dX

    // calculate current from flux:
    return totalflux * currentFromFlux;
}

/*===============================================================================================
 * BEFORE SIMULATION
 *=============================================================================================*/

void Simulation::prepareSimulation()
{
    // initialize all simulation/grid parameters:
    initParametersEtc();
    // initialize all vectors and matrices:
    initVectorsEtc();

    // set boundary conditions in matrix:
    addRedoxToMatrix();
    // fill matrix with backwards implicit (Laasonen) coefficients:
    addBICoeffsToMatrix();
    // add (1st and 2nd order) kinetics to matrix:
    addKineticsToMatrix();
    // create matrix:
    msys.createMatrix();
}

void Simulation::initParametersEtc()
{
    // some parameters
    sys.finalize(el.epsilon);
    f = CONST_F / (CONST_R * env.temperature);
    totalTime = exper.totalTime();
    totalTheta = f * totalTime * exper.scanRate;
    paramR0 = el.epsilon / sqrt(sys.getDiffMax()*totalTime);

    /* system dimensioning according to Compton (Understanding Voltammetry)
     * deltaTheta is a fixed value, taken from Compton
     * to increase accuracy, we need to decrease deltaX, which doesn't need a lot of extra time to solve!!!
     * the factor F in deltaX = pow(10.0, -F-0.5*log10(sigma)); should e.g. go from 2.2 to 4.2
     */
    double F,logRate = log10(sys.getRateMax());
    if (logRate < minLogRate) { F = minF; }
    else if (logRate <= maxLogRate) { F = minF + (logRate - minLogRate) * (maxF - minF) / (maxLogRate - minLogRate); }
    else { F = maxF; }

    sigma = el.epsilon * el.epsilon / sys.getDiffMax() * f * exper.scanRate; // dimensionless scan rate
    maxT = totalTheta / sigma;
    deltaT = deltaTheta / sigma;
    maxX = 6.0*sqrt(maxT);
    //deltaX = pow(10.0, -deltaXSizingFactor-0.5*log10(sigma));
    deltaX = pow(10.0, -F) / sqrt(sigma); // Compton (page 64) leads to: dX = 10^(-2.2-0.5*sigma)
    paramLambda = deltaT / (deltaX * deltaX);

    // print scaling parameters:
    //output_stream << "Diff coeff [m2/s]: max = " << sys.getDiffMax() << endl;
    //output_stream << "Conc [mol/m3]: max = " << sys.getConcMax() << endl;
    //output_stream << "Electrode [m]: epsilon = " << el.epsilon << ", area = " << el.electrodeArea << endl;
    output_stream << "Initial concentrations after equilibration: (";
    for (auto spec: sys.vecSpecies)
        output_stream << "[" << spec->getName() << "] = " << spec->getConcNormEquil() << ((spec!=sys.vecSpecies.back())?", ":"");
    output_stream << ")<br />" << endl;
    output_stream << "V<sub>step</sub> [V] = " << deltaTheta/f << ", t<sub>max</sub> [s] = " << totalTime << ", " << endl;
    output_stream << "r<sub>max</sub> [1/s] = " << sys.getRateMax() << ", F = " << F << "<br />" << endl;
    output_stream << "&sigma; = " << sigma << ", &theta;<sub>max</sub> = " << totalTheta << ", &Delta;&theta; = " << deltaTheta << ", " << endl;
    output_stream << "X<sub>max</sub> = " << maxX << ", &Delta;X = " << deltaX << ", T<sub>max</sub> = " << maxT << ", &Delta;T = " << deltaT << "<br />" << endl;

    // create expanding grid and set dimensions:
    numGridPoints = 1;
    do { numGridPoints++; } while (deltaX < maxX * (paramGamma - 1.0) / ( pow(paramGamma, numGridPoints-1) - 1.0 ));
    numCurrentPoints = min(numCurrentPoints, numGridPoints); // truncate number of current points if it is larger than number of grid points
    numSpecies = sys.vecSpecies.size();
    numRedox = sys.vecRedox.size();
    numReactions = sys.vecReactions.size();
    
    // output parametrisation of grid:
    output_stream << "Number of grid points = " << numGridPoints << "<br />" << endl;
    output_stream << "Number of current points = " << numCurrentPoints << endl;
    output_stream << "Number of deriv points = " << numDerivPoints << endl;
}

void Simulation::initVectorsEtc()
{
    // create matrix system (num_rows)
    msys.initialize(numSpecies, numGridPoints);
    
    // initialize independent terms and b & x vectors:
    independentTerms = Eigen::MatrixXd::Zero(numSpecies, numGridPoints);
    gridConcentration = Eigen::MatrixXd(numSpecies, numGridPoints);
    gridConcentrationPrevious = Eigen::MatrixXd::Zero(numSpecies, numGridPoints);
    
    // initialize distance in grid, expansion factor gamma and concentration vectors:
    gridCoordinate.resize(numGridPoints);
    paramGammai.resize(numGridPoints);
    paramGamma2i.resize(numGridPoints);
    for (size_t x = 0; x < numGridPoints; x++)
    {
        paramGammai[x] = pow(paramGamma, x);
        paramGamma2i[x] = paramGammai[x]*paramGammai[x];
        gridCoordinate[x] = (x == 0) ? 0.0 : gridCoordinate[x-1] + deltaX*paramGammai[x-1];

        for (auto spec: sys.vecSpecies)
            gridConcentration(spec->getIndex(), x) = spec->getConcNormEquil();
    }
    
    // obtain flux condition coefficients:
    coeffBeta0.resize(numCurrentPoints);
    for(size_t i = 0; i < numCurrentPoints; i++)
        coeffBeta0[i] = Beta_N_1(static_cast<int>(numCurrentPoints), static_cast<int>(i), paramGamma);
    
    // obtain diffusion coefficients:
    coeffAlpha.resize(numDerivPoints);
    coeffBeta.resize(numDerivPoints);
    for(size_t d = 0; d < numDerivPoints; d++)
    {
        coeffAlpha[d] = Alpha_N_2(static_cast<int>(numDerivPoints), static_cast<int>(d)-1, paramGamma);
        coeffBeta[d]  =  Beta_N_2(static_cast<int>(numDerivPoints), static_cast<int>(d)-1, paramGamma);
    }

    // initilize redox/current vectors & matrix:
    speciesInRedox.clear();
    speciesInRedox.resize(numSpecies, false);
    matrixB0DiagonalTerms.clear();
    matrixB0DiagonalTerms.resize(numSpecies, 0.0);
    currentContributionMatrix = Eigen::MatrixXd::Zero(numSpecies, numSpecies);
    currentContributionSpeciesFlux = Eigen::VectorXd::Zero(numSpecies);
    currentContributionRedoxFlux = Eigen::VectorXd::Zero(numRedox);
}

void Simulation::addKineticsToMatrix()
{
    for (auto rxn: sys.vecReactions) {
        addHalfReactionKineticsToMatrix(rxn->getSpecLHS1(), rxn->getSpecLHS2(),
                                        rxn->getSpecRHS1(), rxn->getSpecRHS2(), rxn->getKfNorm()); // FORWARD REACTIONS
        addHalfReactionKineticsToMatrix(rxn->getSpecRHS1(), rxn->getSpecRHS2(),
                                        rxn->getSpecLHS1(), rxn->getSpecLHS2(), rxn->getKbNorm()); // BACKWARD REACTIONS
    }
}

void Simulation::addHalfReactionKineticsToMatrix(Species *f1, Species *f2, Species *b1, Species *b2, double normrate)
{
    if (f1 != nullptr && f2 != nullptr)
    {
        // second-order reaction
        addSecondOrderKineticTerm(f1, f2, f1, -normrate);
        addSecondOrderKineticTerm(f1, f2, f2, -normrate);
        addSecondOrderKineticTerm(f1, f2, b1, normrate);
        addSecondOrderKineticTerm(f1, f2, b2, normrate);
    }
    else if (f1 != nullptr || f2 != nullptr)
    {
        // first-order reaction
        Species *f = (f1 != nullptr) ? f1 : f2;
        addFirstOrderKineticTerm(f, f, -normrate);
        addFirstOrderKineticTerm(f, b1, normrate);
        addFirstOrderKineticTerm(f, b2, normrate);
    }
    else
    {
        // reaction has no products or reactants; is that a problem?
    }
}

void Simulation::addFirstOrderKineticTerm(Species *spec1, Species *spec2, double normrate)
{
    if (abs(normrate) > MIN_RATE && spec2 != nullptr)
        for(size_t x = 1; x < numGridPoints; x++)
            msys.addToMatrix(spec2->getIndex(), x, spec1->getIndex(), x, normrate*deltaX*deltaX*paramGamma2i[x]);
}

void Simulation::addSecondOrderKineticTerm(Species *spec1, Species *spec2, Species *spec3, double normrate)
{
    if (abs(normrate) > MIN_RATE && spec3 != nullptr)
        secondOrderKinetics.push_back( make_tuple(spec1->getIndex(), spec2->getIndex(), spec3->getIndex(), normrate) );
}

void Simulation::addRedoxToMatrix()
{
    double p, Kred, Kox;

    for (auto redox: sys.vecRedox)
    {
        p = static_cast<double>(redox->getNe()) * f * (exper.initialPotential - redox->getE0()); // normalized potential

        Kred = deltaX * redox->getKeNorm() * exp(      -redox->getAlpha()  * p); // K_red B-V kinetics ("Kf")
        Kox  = deltaX * redox->getKeNorm() * exp((1.0 - redox->getAlpha()) * p); // K_ox  B-V kinetics ("Kb")

        // add terms per redox reaction:
        msys.addToMatrix(redox->getSpecOx()->getIndex(),  0,
                         redox->getSpecRed()->getIndex(), 0,
                         -Kox / redox->getSpecOx()->getDiffNorm());  //B0
        msys.addToMatrix(redox->getSpecRed()->getIndex(), 0,
                         redox->getSpecOx()->getIndex(),  0,
                         -Kred / redox->getSpecRed()->getDiffNorm()); //B0
        
        // keep track of which species partake in redox steps:
        speciesInRedox[redox->getSpecOx()->getIndex()] = true;
        speciesInRedox[redox->getSpecRed()->getIndex()] = true;
        
        // add to diagonal terms:
        matrixB0DiagonalTerms[redox->getSpecOx()->getIndex()] += Kred / redox->getSpecOx()->getDiffNorm();
        matrixB0DiagonalTerms[redox->getSpecRed()->getIndex()] += Kox / redox->getSpecRed()->getDiffNorm();

        // set current contribution matrix:
        currentContributionMatrix.coeffRef(static_cast<Eigen::Index>(redox->getSpecOx()->getIndex()), static_cast<Eigen::Index>(redox->getIndex())) = 1.0;
        currentContributionMatrix.coeffRef(static_cast<Eigen::Index>(redox->getSpecRed()->getIndex()), static_cast<Eigen::Index>(redox->getIndex())) = -1.0;
    }

    // add zero flux condition for all species NOT in redox step:
    for (size_t s = 0; s < numSpecies; s++)
    {
        if (!speciesInRedox[s]) // species does not participate in any redox step
        {
            for (size_t x = 0; x < numCurrentPoints; x++)
            {
                // instead of zero flux condition function:
                msys.addToMatrix(s, 0, s, x, coeffBeta0[x]/deltaX);
            }
        }
        else // active redox species
        {
            msys.addToMatrix(s, 0, s, 0,  1.0 + matrixB0DiagonalTerms[s]);  // B0
            msys.addToMatrix(s, 0, s, 1, -1.0);  // B1
        }
    }
}

void Simulation::addBICoeffsToMatrix() // Backwards Implicit coefficients
{
    size_t relidx_max;
    for (auto spec: sys.vecSpecies)
    {
        for (size_t x = 1; x < numGridPoints; x++) // since x == 0 corresponds to the boundary condition
        {
            if (x < numGridPoints - numDerivPoints + 2)
                relidx_max = numDerivPoints-1;
            else
                relidx_max = numGridPoints-x;
            
            for (size_t relidx = 0; relidx < relidx_max+1; relidx++)
            {
                msys.addToMatrix(spec->getIndex(), x, spec->getIndex(), x+relidx-1,
                                 coeffMatrixN2(x, static_cast<int>(relidx)-1, spec->getDiffNorm()));
            }
        }
    }
}

double Simulation::coeffMatrixN2(size_t x, int relative_position, double diffnorm)
{
    size_t coeffidx = static_cast<size_t>(relative_position+1); // relative_position >= -1
    
    double coeff = coeffAlpha[coeffidx] + el.electrodeGeometryFactor*deltaX*paramGammai[x]/(paramR0+gridCoordinate[x])*coeffBeta[coeffidx];
    coeff *= diffnorm; // Compton page 88/89
    coeff -= (relative_position == 0) ? paramGamma2i[x]/paramLambda : 0.0;
    
    return coeff;
}



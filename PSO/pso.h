#ifndef PSO_H
#define PSO_H
#include <iostream>
#include <random>
#include <ctime>
using namespace std;

std::random_device rd;        // Will be used to obtain a seed for the random number engine
std::mt19937 generator(rd()); // Standard mersenne_twister_engine seeded with rd()
class PSO
{
public:
    // Constructor
    PSO(double (*function)(const double *, const int),
        int numDim,
        double xMax,
        double xMin,
        double k,
        double c1,
        double c2,
        int numParticle = 20,
        double wMax = 0.9,
        double wMin = 0.4)
        : _function(function),
          _numDim(numDim),
          _xMax(xMax),
          _xMin(xMin),
          _numParticle(numParticle),
          _maxIter(numDim * 10000),
          _wMax(wMax),
          _wMin(wMin),
          _c1(c1),
          _c2(c2),
          _k(k)
    {
        // init particles
        for (int i = 0; i < _numParticle; i++)
        {
            X[i] = new double[_numDim];
            pBest[i] = new double[_numDim];
            V[i] = new double[_numDim];
            for (int j = 0; j < _numDim; j++)
            {
                // random particle position from uniform distribution
                uniform_real_distribution<double> unif_X(_xMin, _xMax);
                X[i][j] = unif_X(generator);
                pBest[i][j] = X[i][j];
                pBestFitness[i] = _function(X[i], _numDim);

                // random velocity for each particle from uniform distribution
                uniform_real_distribution<double>
                    unif_V(_vMin, _vMax);
                V[i][j] = unif_V(generator);
                // V[i][j] = 0;
            }
        }
        gBest = X[0];
        gBestFitness = pBestFitness[0];
        for (int i = 1; i < _numParticle; i++)
        {
            if (pBestFitness[i] < gBestFitness)
            {
                gBest = X[i];
                gBestFitness = pBestFitness[i];
            }
        }
    }
    void update()
    {
        // start iteration
        for (int t = 0; t < _maxIter; t++)
        {
            for (int i = 0; i < _numParticle; i++)
            {
                // calculate w (linear decrease)
                double w = _wMax - t * (_wMax - _wMin) / _maxIter;
                // double w = 0.8;
                // update velocity
                double *nextV = new double[_numDim];
                uniform_real_distribution<double> unif_R1(0, 1);
                uniform_real_distribution<double> unif_R2(0, 1);
                double r1 = unif_R1(generator);
                double r2 = unif_R2(generator);

                for (int j = 0; j < _numDim; j++)
                {
                    nextV[j] = w * V[i][j] + _c1 * (pBest[i][j] - X[i][j]) * r1 + _c2 * (gBest[j] - X[i][j]) * r2;

                    if (nextV[j] > _vMax)
                        nextV[j] = _vMax;
                    else if (nextV[j] < _vMin)
                        nextV[j] = _vMin;
                    V[i][j] = nextV[j];

                    X[i][j] += nextV[j];
                    if (X[i][j] > _xMax)
                        X[i][j] = _xMax;
                    else if (X[i][j] < _xMin)
                        X[i][j] = _xMin;
                }
                delete[] nextV;
            }
            for (int i = 0; i < _numParticle; i++)
            {
                // Calculate particle fitness
                double currentFitness = _function(X[i], _numDim);

                // update Personal Best
                if (currentFitness <= pBestFitness[i])
                {
                    pBest[i] = X[i];
                    pBestFitness[i] = currentFitness;
                }
                // update Global Best
                if (currentFitness <= gBestFitness)
                {
                    gBest = X[i];
                    gBestFitness = currentFitness;
                }
            }
        }
    }
    double get_gBestFitness()
    {
        return gBestFitness;
    }
    double *get_gBest()
    {
        return gBest;
    }

private:
    double (*_function)(const double *, const int);  // Pointer to test function
    int _numDim;                                     // Dimension
    double _xMax;                                    // Max of x
    double _xMin;                                    // Min of x
    int _numParticle;                                // Number of particles
    int _maxIter;                                    // Max iterations number
    double _wMax;                                    // Max of inertia weight
    double _wMin;                                    // Min inertia weight
    double _c1;                                      // Cognitive parameter
    double _c2;                                      // Social parameter
    double _k;                                       // Constant k
    double _vMin = -_k * (_xMax - _xMin) / 2;        // Min of Velocity
    double _vMax = _k * (_xMax - _xMin) / 2;         // Max of Velocity
    double **X = new double *[_numParticle];         // Particle X[_numParticle][_numDim]
    double **pBest = new double *[_numParticle];     // Personal Best pBest[_numParticle][_numDim]
    double *pBestFitness = new double[_numParticle]; // Personal Best fitness value
    double *gBest = new double[_numDim];             // Global Best gBest[_numDim]
    double gBestFitness;                             // Global Best fitness value
    double **V = new double *[_numParticle];         // Velocity [_numParticle][_numDim]
};

double Ackley(const double *x, const int d)
{
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < d; ++i)
    {
        sum1 += x[i] * x[i];
        sum2 += cos(2.0 * M_PI * x[i]);
    }
    return (-20.0) * exp(-0.2 * sqrt(sum1 / double(d))) - exp(sum2 / double(d)) + 20.0 + exp(1.0);
}

double Rastrigin(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 0; i < d; ++i)
    {
        sum1 += (x[i] * x[i]) - (10.0 * cos(2.0 * M_PI * x[i]));
    }
    return sum1 + 10.0 * double(d);
}

double HappyCat(const double *x, const int d)
{
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < d; ++i)
    {
        sum1 += x[i] * x[i];
        sum2 += x[i];
    }
    return pow(fabs(sum1 - double(d)), 0.25) + (0.5 * sum1 + sum2) / double(d) + 0.5;
}

double Rosenbrock(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 0; i < d - 1; ++i)
    {
        sum1 += 100.0 * (x[i + 1] - (x[i] * x[i])) * (x[i + 1] - (x[i] * x[i])) + ((x[i] - 1.0) * (x[i] - 1.0));
    }
    return sum1;
}

double Zakharov(const double *x, const int d)
{
    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < d; ++i)
    {
        sum1 += x[i] * x[i];
        sum2 += 0.5 * (i + 1) * x[i];
    }
    return sum1 + pow(sum2, 2) + pow(sum2, 4);
}

double Michalewicz(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 0; i < d; ++i)
    {
        sum1 += sin(x[i]) * pow(sin((double(i + 1) * x[i] * x[i]) / M_PI), 20.0);
    }
    return sum1 * (-1);
}

#endif
#ifndef PSO_H
#define PSO_H
#include <iostream>
#include <random>
#include <ctime>
#include <iomanip>

using namespace std;

std::random_device rd;
std::mt19937 generator(rd());

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
        int numParticle,
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
        // 1. initialzation
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
                // init velocity  0
                V[i][j] = 0;
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
            // 2. transition
            for (int i = 0; i < _numParticle; i++)
            {
                // calculate w (linear decrease)
                double w = _wMax - t * (_wMax - _wMin) / _maxIter;
                // double w = 0.5;

                // update velocity
                double *nextV = new double[_numDim];
                uniform_real_distribution<double> unif(0, 1);
                double r1 = unif(generator);
                double r2 = unif(generator);
                for (int j = 0; j < _numDim; j++)
                {
                    V[i][j] = w * V[i][j] + _c1 * (pBest[i][j] - X[i][j]) * r1 + _c2 * (gBest[j] - X[i][j]) * r2;

                    if (V[i][j] > _vMax)
                        V[i][j] = _vMax;
                    else if (V[i][j] < _vMin)
                        V[i][j] = _vMin;

                    X[i][j] += V[i][j];
                    if (X[i][j] > _xMax)
                        X[i][j] = _xMax;
                    else if (X[i][j] < _xMin)
                        X[i][j] = _xMin;
                }
                delete[] nextV;
            }

            // new idea
            if (t > 0)
            {
                double improvement = prevGBestFitness - gBestFitness;
                if (improvement == 0)
                    localopt_counter++;
                if (localopt_counter >= _threshold)
                {
                    // cout << "局部最佳解 :" << gBestFitness << endl;
                    double *shadowParticle = new double[_numDim];

                    // 局部搜索
                    for (int i = 0; i < 10; ++i)
                    {
                        int index = rand() % _numParticle; // 隨機選擇一個粒子

                        shadowParticle = X[index];
                        for (int j = 0; j < _numDim / 2; ++j)
                        {
                            uniform_real_distribution<double> unif(0, 1);
                            uniform_int_distribution<int> unif_int(0, _numDim - 1);

                            int dim = unif_int(generator);    // 隨機選擇一個維度
                            double randVal = unif(generator); // 生成一個隨機數
                            shadowParticle[dim] += randVal;   // 在選擇的維度上加上隨機數
                            // 計算新粒子的適應度
                            double newFitness = _function(shadowParticle, _numDim);
                            // 將新粒子的適應度與原始粒子的適應度進行比較，如果新粒子的適應度更佳，則取代原始粒子
                            if (newFitness < pBestFitness[index])
                            {
                                pBestFitness[index] = newFitness;
                                for (int j = 0; j < _numDim; ++j)
                                {
                                    pBest[index][j] = shadowParticle[j];
                                }
                                // cout << "success" << endl;
                            }
                        }
                    }

                    localopt_counter = 0;
                }
            }

            // 更新 prevGBestFitness
            prevGBestFitness = gBestFitness;

            // 3. evalution
            for (int i = 0; i < _numParticle; i++)
            {
                // Calculate particle fitness
                double currentFitness = _function(X[i], _numDim);

                // update Personal Best
                if (currentFitness <= pBestFitness[i])
                {
                    for (int j = 0; j < _numDim; j++)
                    {
                        pBest[i][j] = X[i][j];
                    }
                    pBestFitness[i] = currentFitness;
                }
                // update Global Best
                if (currentFitness <= gBestFitness)
                {
                    for (int j = 0; j < _numDim; j++)
                    {
                        gBest[j] = X[i][j];
                    }
                    gBestFitness = currentFitness;
                }
            }
            gBest_list[t] = gBestFitness;
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
    double *get_gBest_list()
    {
        return gBest_list;
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
    double *gBest_list = new double[_maxIter];       // Record GBest
    double gBestFitness;                             // Global Best fitness value
    double **V = new double *[_numParticle];         // Velocity [_numParticle][_numDim]
    double prevGBestFitness;
    int localopt_counter = 0;
    int _threshold = 100;
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

double Schwefel(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 0; i < d; ++i)
    {
        sum1 += x[i] * sin(sqrt(fabs(x[i])));
    }
    return 418.9829 * double(d) - sum1;
}

double BentCigar(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 1; i < d; ++i)
    {
        sum1 += x[i] * x[i];
    }
    return x[0] * x[0] + pow(10.0, 6) * sum1;
}

double DropWave(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 0; i < d; ++i)
    {
        sum1 += x[i] * x[i];
    }
    return 1.0 - ((1.0 + cos(12.0 * sqrt(sum1))) / (0.5 * sum1 + 2.0));
}

double Step(const double *x, const int d)
{
    double sum1 = 0.0;
    for (int i = 0; i < d; ++i)
    {
        sum1 += floor(x[i] + 0.5) * floor(x[i] + 0.5);
    }
    return sum1;
}

#endif
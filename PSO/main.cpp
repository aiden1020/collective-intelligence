#include "pso.h"
#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <function_type> <run> <dimension> <k> <c1> <c2>" << endl;
        return 1;
    }

    string function_type = argv[1];
    double run = atof(argv[2]);
    int dimension = atoi(argv[3]);
    double k = atof(argv[4]);
    double c1 = atof(argv[5]);
    double c2 = atof(argv[6]);

    // Define null function pointer
    typedef double (*FunctionPtr)(const double *, const int);
    FunctionPtr function = nullptr;

    double upper_bound, lower_bound;
    if (function_type == "Ackley")
    {
        function = Ackley;
        upper_bound = 32.768;
        lower_bound = -32.768;
    }
    else if (function_type == "Rastrigin")
    {
        function = Rastrigin;
        upper_bound = 5.12;
        lower_bound = -5.12;
    }
    else if (function_type == "HappyCat")
    {
        function = HappyCat;
        upper_bound = 20.0;
        lower_bound = -20.0;
    }
    else if (function_type == "Rosenbrock")
    {
        function = Rosenbrock;
        upper_bound = 10.0;
        lower_bound = -10.0;
    }
    else if (function_type == "Zakharov")
    {
        function = Zakharov;
        upper_bound = 10.0;
        lower_bound = -10.0;
    }
    else if (function_type == "Michalewicz")
    {
        function = Michalewicz;
        upper_bound = M_PI;
        lower_bound = 0.0;
    }
    else
    {
        cerr << "Unknown function type: " << function_type << endl;
        return 1;
    }
    double avg_gBest = 0;
    if (function != nullptr)
    {
        for (int i = 0; i < run; i++)
        {
            PSO pso(function, dimension, upper_bound, lower_bound, k, c1, c2);
            pso.update();
            double *gBest_ = pso.get_gBest();
            double gBestFitness_ = pso.get_gBestFitness();
            cout << "Run" << i + 1 << " Global best Fitness :" << endl;
            cout << gBestFitness_ << endl;
            cout << "gBest : " << endl;
            for (int j = 0; j < dimension; j++)
            {
                cout << gBest_[j] << " ";
            }
            cout << endl
                 << "-------------------------------"
                 << endl;
            avg_gBest += gBestFitness_;
        }
        avg_gBest /= run;
        cout << endl
             << "Avg Global Best " << avg_gBest << endl;
    }

    return 0;
}

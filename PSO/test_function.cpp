#include <iostream>
#include <random>
#include <ctime>
using namespace std;
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
int main()
{
    double x[2] = {0.0, 0.0};
    cout << Ackley(x, 2);
    return 0;
}

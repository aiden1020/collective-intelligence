#include "pso.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream> // 用於 ostringstream

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <function_type> <run> <dimension> <k> <c1> <c2> <numParticle>" << endl;
        return 1;
    }
    cout << scientific << setprecision(3);
    string function_type = argv[1];
    double run = atof(argv[2]);
    int dimension = atoi(argv[3]);
    double k = atof(argv[4]);
    double c1 = atof(argv[5]);
    double c2 = atof(argv[6]);
    int numParticle = atoi(argv[7]);
    // 使用 ostringstream 來格式化文件名
    ostringstream filename_stream;
    // filename_stream << "result/50P/" << function_type << "_" << dimension << "D.txt";
    filename_stream << "result/coverage/" << function_type << "_" << dimension << "D.txt";

    string filename = filename_stream.str();
    ofstream outputFile(filename);
    if (!outputFile.is_open())
    {
        cerr << "Error opening file!" << endl;
        return 1;
    }
    // 續寫模式打開 run_time.txt
    ofstream runTimeFile("run_time.txt", ios::app);
    if (!runTimeFile.is_open())
    {
        cerr << "Error opening run_time.txt file!" << endl;
        return 1;
    }

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
    else if (function_type == "Schwefel")
    {
        function = Schwefel;
        upper_bound = 500.0;
        lower_bound = -500.0;
    }
    else if (function_type == "BentCigar")
    {
        function = BentCigar;
        upper_bound = 100;
        lower_bound = -100;
    }
    else if (function_type == "DropWave")
    {
        function = DropWave;
        upper_bound = 5.12;
        lower_bound = -5.12;
    }
    else if (function_type == "Step")
    {
        function = Step;
        upper_bound = 100;
        lower_bound = -100;
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
            double START, END;
            START = clock();
            PSO pso(function, dimension, upper_bound, lower_bound, k, c1, c2, numParticle);
            pso.update();
            END = clock();
            double *gBest_ = pso.get_gBest();
            double gBestFitness_ = pso.get_gBestFitness();
            double *gBest_list_ = pso.get_gBest_list();
            cout << "Run" << i + 1 << " Global best Fitness :" << endl;
            cout << gBestFitness_ << endl;
            cout << "gBest : " << endl;
            for (int j = 0; j < dimension; j++)
            {
                cout << gBest_[j] << " ";
            }
            cout << endl
                 << "run time : " << (END - START) / CLOCKS_PER_SEC << " S" << endl;
            cout << endl
                 << "-------------------------------"
                 << endl;
            avg_gBest += gBestFitness_;
            // output file
            for (int j = 0; j < dimension * 10000; j++)
            {
                outputFile << scientific << setprecision(3) << gBest_list_[j] << endl;
            }
            // outputFile << scientific << setprecision(3) << gBestFitness_ << endl;

            // 將 function_type 和 runtime 寫入 run_time.txt
            double runtime = (END - START) / CLOCKS_PER_SEC;
            runTimeFile << function_type << " , " << runtime << endl;
        }
        avg_gBest /= run;
        cout << endl
             << "Avg Global Best " << avg_gBest << endl;
    }

    runTimeFile.close(); // 關閉 run_time.txt 文件
    return 0;
}

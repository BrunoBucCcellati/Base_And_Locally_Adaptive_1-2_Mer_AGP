#include <random>
#define M_PI 3.14159265358979323846
using namespace std;
template <class Value>
short Sign(Value Val)
{
    if (Val == 0.)  return 0;
    if (Val > 0.)  return 1;
    else return -1;
}
extern "C" __declspec(dllexport) double HillFunc(double x)
{
    const int HH_COUNT = 14;
    double a = 0.0;
    double b = 1.0;
    //
    default_random_engine generator;
    uniform_real_distribution<double> distribution(a, b);
    double A[20] = {};
    double B[20] = {};
    //
    for (int i = 0; i < 20; ++i)
    {
        A[i] = -1.1 + distribution(generator) * 2.0;
        B[i] = -1.1 + distribution(generator) * 2.0;
    }
    //
    double res = A[0];
    for (int i = 1; i < HH_COUNT; i++)
    {
        res += A[i] * sin(i * 2 * M_PI * x) + B[i] * cos(i * 2 * M_PI * x);
    }
    return res;
}
extern "C" __declspec(dllexport) double ShekelFunc(double x)
{
    const int SH_COUNT = 10;
    double a = 0.0;
    double b = 1.0;
    //
    default_random_engine generator;
    uniform_real_distribution<double> distribution(a, b);
    double A[20] = {};
    double B[20] = {};
    double K[20] = {};
    //
    for (int i = 0; i < 20; ++i)
    {
        A[i] = 10.0 * distribution(generator);
        B[i] = 1.0 + 0.2 * distribution(generator) * 2.0;
        K[i] = 5.0 + 20 * distribution(generator);
    }
    //
    double res = 0.0;
    for (int i = 1; i < SH_COUNT; i++)
    {
        res -= 1.0 / ((K[i] * pow((x - A[i]), 2.0) + B[i]));
    }
    return res;
}
extern "C" __declspec(dllexport) double GrishaginFunc(double x1, double x2)
{
    const int GR_COUNT = 8;
    double a = -1.0;
    double b = 1.0;
    //
    default_random_engine generator;
    uniform_real_distribution<double> distribution(a, b);
    double a_x1x2[20][20] = {};
    double b_x1x2[20][20] = {};
    double A[20][20] = {};
    double B[20][20] = {};
    double C[20][20] = {};
    double D[20][20] = {};
    //
    double part1 = 0.0, part2 = 0.0;
    //
    for (int i = 0; i < 20; ++i)
    {
        for (int j = 0; j < 20; ++j)
        {
            a_x1x2[i][j] = sin(M_PI * i * x1) * sin(M_PI * j * x2);
            b_x1x2[i][j] = cos(M_PI * i * x1) * cos(M_PI * j * x2);
            A[i][j] = distribution(generator);
            B[i][j] = distribution(generator);
            C[i][j] = distribution(generator);
            D[i][j] = distribution(generator);
            //
            if (i >= 1 && i <= 7 && j >= 1 && j <= 7)
            {
                part1 += A[i][j] * a_x1x2[i][j] + B[i][j] * b_x1x2[i][j];
                part2 += C[i][j] * a_x1x2[i][j] - D[i][j] * b_x1x2[i][j];
            }
        }
    }
    //
    double res = -pow((pow(part1, 2) + pow(part2, 2)), 0.5);
    return res;
}
extern "C" __declspec(dllexport) double Characteristic(double _m, double x1, double x2, double y1, double y2, unsigned short _N)
{
    double delta_x = pow(abs(x2 - x1), (1 / double(_N)));
    return ((delta_x) + pow((y2 - y1), 2) / (pow(_m, 2) * (delta_x)) - 2 * (y2 + y1) / _m);
}
extern "C" __declspec(dllexport) double Shag(double _m, double x1, double x2, double y1, double y2, unsigned short _N)
{
    return ((x1 + x2) / 2) - Sign(y2 - y1) * pow((abs(y2 - y1) / (2 * _m)), _N);
}
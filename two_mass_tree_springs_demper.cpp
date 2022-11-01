#include <iostream>
#include <tuple>
#include <vector>
#include <math.h>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//----------- system dynamic parameters --------------------

float m1 = 1.0;
float m2 = 1.0;
float k1 = 3.0;
float k2 = 3.0;
float k3 = 3.0;
float b1 = 0.5;
float b2 = 0.5;
float b3 = 0.5;

float dt = 0.001;

//-----------------------------------------------------------
// x1_dot
float function1(float x1, float x2, float x1_dot, float x2_dot)
{

    return x1_dot;
}

//-----------------------------------------------------------
// x1_dot_dot
float function3(float x1, float x2, float x1_dot, float x2_dot)
{
    // m1 x1'' = −k1 (x1 − R1) + k2 (x2 − x1 − w1 − R2)

    float x1_dot_dot = (-(b1 + b2) * x1_dot + b2 * x2_dot - (k1 + k2) * x1 + k2 * x2) / m1;
    return x1_dot_dot;
}

//-----------------------------------------------------------
// x2_dot
float function2(float x1, float x2, float x1_dot, float x2_dot)
{

    return x2_dot;
}

//-----------------------------------------------------------
// x2_dot_dot
float function4(float x1, float x2, float x1_dot, float x2_dot)
{

    float x2_dot_dot = (b2 * x1_dot - (b2 + b3) * x2_dot + k2 * x1 - (k2 + k3) * x2) / m2;

    return x2_dot_dot;
}

//-----------------------------------------------------------

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> methodRK4_TwoMass()
{

    std::vector<float> diffEq1;
    std::vector<float> diffEq2;
    std::vector<float> diffEq3;
    std::vector<float> diffEq4;

    std::vector<float> time;

    // init values
    float x1 = 3.2; //
    float x2 = 0.0; //
    float x3 = 6.4; // theta1_dot
    float x4 = 0;   // theta2_dot
    float t = 0.0;  // init time

    diffEq1.push_back(x1);
    diffEq2.push_back(x2);
    diffEq3.push_back(x3);
    diffEq4.push_back(x4);
    time.push_back(t);

    for (int ii = 0; ii < 50000; ii++)
    {
        t = t + dt;
        float k11 = function1(x1, x2, x3, x4);
        float k12 = function2(x1, x2, x3, x4);
        float k13 = function3(x1, x2, x3, x4);
        float k14 = function4(x1, x2, x3, x4);

        float k21 = function1(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k22 = function2(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k23 = function3(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k24 = function4(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);

        float k31 = function1(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k32 = function2(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k33 = function3(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k34 = function4(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);

        float k41 = function1(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k42 = function2(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k43 = function3(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k44 = function4(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);

        x1 = x1 + dt / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        x2 = x2 + dt / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);
        x3 = x3 + dt / 6.0 * (k13 + 2 * k23 + 2 * k33 + k43);
        x4 = x4 + dt / 6.0 * (k14 + 2 * k24 + 2 * k34 + k44);

        diffEq1.push_back(x1);
        diffEq2.push_back(x2);
        diffEq3.push_back(x3);
        diffEq4.push_back(x4);
        time.push_back(t);
    }

    return std::make_tuple(diffEq1, diffEq2, diffEq3, diffEq4, time);
}

//---------------------------------------------------------------------------------------------------------

void plot2D(std::tuple<std::vector<float>, std::vector<float>> data1)
{

    std::vector<float> xX1 = std::get<0>(data1);
    std::vector<float> yY1 = std::get<1>(data1);

    plt::plot(xX1, yY1);
    plt::xlabel("pos[m]");
    plt::ylabel("speed[m/s]");
    plt::show();
}

//---------------------------------------------------------------

void plot2D2D(std::tuple<std::vector<float>, std::vector<float>> data1, std::tuple<std::vector<float>, std::vector<float>> data2)
{

    std::vector<float> xX1 = std::get<0>(data1);
    std::vector<float> yY1 = std::get<1>(data1);

    std::vector<float> xX2 = std::get<0>(data2);
    std::vector<float> yY2 = std::get<1>(data2);

    plt::plot(xX1, yY1);
    plt::plot(xX2, yY2);
    plt::xlabel("time[s]");
    plt::ylabel("A[m]");
    plt::show();
}

//---------------------------------------------------------------
int main()
{

    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> dyn2mass = methodRK4_TwoMass();

    std::tuple<std::vector<float>, std::vector<float>> mass12 = std::make_tuple(std::get<0>(dyn2mass), std::get<2>(dyn2mass));

    std::tuple<std::vector<float>, std::vector<float>> mass1t = std::make_tuple(std::get<4>(dyn2mass), std::get<0>(dyn2mass));
    std::tuple<std::vector<float>, std::vector<float>> mass2t = std::make_tuple(std::get<4>(dyn2mass), std::get<2>(dyn2mass));

    plot2D(mass12);
    plot2D2D(mass1t, mass2t);
}

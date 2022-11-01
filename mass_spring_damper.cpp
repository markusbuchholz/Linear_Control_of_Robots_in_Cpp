#include <iostream>
#include <tuple>
#include <vector>
#include <math.h>
#include <cmath>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//----------- system dynamic parameters --------------------

float k = 1;   // spring constant
float b = 0.3; // damping coefficient
float m = 1;   // mass
float kp = 2.0;
float kv = 0.5;
float kpri = k + kp;
float bpri = b + kv;

float dt = 0.001;

//-----------------------------------------------------------
// x_dot
float function1(float x, float x_dot)
{

    return x_dot;
}

//-----------------------------------------------------------
// x_dot_dot
float function2(float x, float x_dot)
{

    float x_dot_dot = (-bpri * x_dot - kpri * x) / m;

    return x_dot_dot;
}

//-----------------------------------------------------------

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>> RK4MotionSpringDamperMass()
{

    std::vector<float> diffEq1;
    std::vector<float> diffEq2;
    std::vector<float> time;

    // init values
    float x1 = 0.3; // x
    float x2 = 0.0; // x_dot
    float t = 0.0;  // init time

    diffEq1.push_back(x1);
    diffEq2.push_back(x2);
    time.push_back(t);

    for (int ii = 0; ii < 30000; ii++)
    {
        t = t + dt;
        float k11 = function1(x1, x2);
        float k12 = function2(x1, x2);

        float k21 = function1(x1 + dt / 2 * k11, x2 + dt / 2 * k12);
        float k22 = function2(x1 + dt / 2 * k11, x2 + dt / 2 * k12);

        float k31 = function1(x1 + dt / 2 * k21, x2 + dt / 2 * k22);
        float k32 = function2(x1 + dt / 2 * k21, x2 + dt / 2 * k22);

        float k41 = function1(x1 + dt * k31, x2 + dt * k32);
        float k42 = function2(x1 + dt * k31, x2 + dt * k32);

        x1 = x1 + dt / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        x2 = x2 + dt / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);

        diffEq1.push_back(x1);
        diffEq2.push_back(x2);
        time.push_back(t);
    }

    return std::make_tuple(diffEq1, diffEq2, time);
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

    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>> motion = RK4MotionSpringDamperMass();

    std::tuple<std::vector<float>, std::vector<float>> motionX = std::make_tuple(std::get<2>(motion), std::get<0>(motion));
    std::tuple<std::vector<float>, std::vector<float>> motionV = std::make_tuple(std::get<2>(motion), std::get<1>(motion));
    std::tuple<std::vector<float>, std::vector<float>> motionXV = std::make_tuple(std::get<0>(motion), std::get<1>(motion));
    plot2D(motionX);
    plot2D(motionV);
    plot2D(motionXV);
    plot2D2D(motionX, motionV);
}

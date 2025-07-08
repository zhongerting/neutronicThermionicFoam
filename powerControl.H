#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <typeinfo>

class powerControl
{
    public:
    bool isFirst;

    public:
    Foam::scalar power;
    Foam::scalar powerOld;
    Foam::scalar powerOldOld;

    Foam::scalar fissPower;
    Foam::scalar decayPower;

    public:
    Foam::scalar powerRelative;
    Foam::scalar powerRelativeOld;

    Foam::scalar powerInitial;
    Foam::scalar powerRelativeInitial;

    std::vector<Foam::scalar> timeInput;
    std::vector<Foam::scalar> dataInput;

    public:
    Foam::scalar fuelReac;
    Foam::scalar emitterReac;
    Foam::scalar collectorReac;
    Foam::scalar moderatorReac;

    Foam::scalar reactFeedback;

    public:
    bool isPower;
    bool isTotalReactivity;
    bool isInputReactivity;

    public:
    powerControl()
    {
    }
    ~powerControl()
    {
    }

    // 读取输入方式
    public:
    std::string getPower(std::string mode)
    {
        if (mode == "power")
        {
            isPower = true;
            isTotalReactivity = false;
            isInputReactivity = false;
        }
        else if (mode == "totalReactivity")
        {
            isPower = false;
            isTotalReactivity = true;
            isInputReactivity = false;
        }
        else if (mode == "inputReactivity")
        {
            isPower = false;
            isTotalReactivity = false;
            isInputReactivity = true;
        }

        return mode;
    }

    // 插值函数，假设均为线性变化
    Foam::scalar linearInterpolation(Foam::scalar time)
    {
        if (time <= timeInput[0])
        {
            Info << "time < timeInput[0]" << endl;
            return dataInput[0];
        }
        else if (time > timeInput[int(timeInput.size()) - 1])
        {
            Info << "time > timeInput[int(timeInput.size()) - 1]" << endl;
            return dataInput[int(timeInput.size()) - 1];
        }
        else
        {
            for (int i = 0; i < int(timeInput.size()) - 1; ++i)
            {
                if (time > timeInput[i] && time <= timeInput[i + 1])
                {
                    return dataInput[i] + (dataInput[i + 1] - dataInput[i]) / (timeInput[i + 1] - timeInput[i]) * (time - timeInput[i]);
                }
            }
        }
        return -1.0;
    }

    // 反应性反馈：燃料反应性反馈
    Foam::scalar reactivityFeedbackFuel(scalar T)
    {
        return 0.001360811 - 6.47927757e-6 * T + 2.321231e-9 * Foam::pow(T, 2) - 3.52e-13 * Foam::pow(T, 3);
    }
    // 反应性反馈：发射极反馈
    Foam::scalar reactivityFeedbackEmitter(scalar T)
    {
        return 1.e-4 * (3.46455e-6 * Foam::pow(T, 2) - 0.03232167 * T + 0.74202216);
    }
    // 反应性反馈：接收极反馈
    Foam::scalar reactivityFeedbackCollector(scalar T)
    {
        return 1.e-4 * (1.7094427e-6 * Foam::pow(T, 2) - 0.017584414 * T + 0.7409044027);
    }
    // 反应性反馈：慢化剂反馈
    Foam::scalar reactivityFeedbackModerator(scalar T)
    {
        return 1.0 * (-1.7251979e-8 * Foam::pow(T, 2) + 0.000059009 * T - 0.0148233537);
    }

    // 反应性反馈计算
    void reactivityCalc(std::vector<Foam::scalar> T)
    {
        fuelReac = reactivityFeedbackFuel(T[0]);
        emitterReac = reactivityFeedbackEmitter(T[1]);
        collectorReac = reactivityFeedbackCollector(T[2]);
        moderatorReac = reactivityFeedbackModerator(T[3]);
        reactFeedback = fuelReac + emitterReac + collectorReac + moderatorReac;
    }
};

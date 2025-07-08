PC1->powerOldOld = PC1->powerOld;
PC1->powerOld = PC1->power;

// 每个时间步都进行若干步中子计算
Foam::scalar timeNow = std::stod(runTime.timeName());
Foam::scalar timeNext = timeNow + runTime.deltaTValue();

Foam::scalar inputReact = 0.0;

// 功率控制
if (PC1->isPower)
{
    PC1->power = PC1->linearInterpolation(timeNow);
    Info << PC1->power << nl;
    PC1->powerRelativeOld = PC1->powerRelative;
    PC1->powerRelative = PC1->power / PC1->powerInitial;
    PC1->isFirst = false;
}
else if (PC1->isTotalReactivity)
{
    inputReact = PC1->linearInterpolation(timeNow);
    ode.react = inputReact;
    Info << ode.react << nl;
    odeSolver->solve(timeNow, timeNext, y, 0, dxEst);
    PC1->fissPower = y[0];
    PC1->decayPower = 0.0;
    for (int i = 0; i < int(ode.decayFra.size()); i++)
    {
        PC1->decayPower += y[7 + i] * ode.decayCon[i];
    }
    PC1->power = PC1->fissPower + PC1->decayPower;
    PC1->powerRelativeOld = PC1->powerRelative;
    PC1->powerRelative = PC1->power / PC1->powerInitial;
    PC1->isFirst = false;
}
else if (PC1->isInputReactivity)
{
    // 计算平均温度
    #include "calculateTemperature.H"
    // 计算总反应性
    ode.react = PC1->fuelReac + PC1->emitterReac + PC1->collectorReac + PC1->moderatorReac;
    if (PC1->isFirst)
    {
        Info << "initial reactivity feedback is " << ode.react << nl;
        Foam::sleep(3);
    }    
    inputReact = PC1->linearInterpolation(timeNow);
    ode.react += inputReact;
    // 计算功率
    odeSolver->solve(timeNow, timeNext, y, 0, dxEst);
    PC1->fissPower = y[0];
    PC1->decayPower = 0.0;
    for (int i = 0; i < int(ode.decayFra.size()); i++)
    {
        PC1->decayPower += y[7 + i] * ode.decayCon[i];
    }
    PC1->power = PC1->fissPower + PC1->decayPower;
    PC1->powerRelativeOld = PC1->powerRelative;
    PC1->powerRelative = PC1->power / PC1->powerInitial;
    PC1->isFirst = false;
}

// 与电输出同步输出文件内容
if (runTime.value() > neutronTimeWrite)
{
    outputPowerFile << runTime.value() << tab << PC1->power <<  tab << PC1->fissPower <<  tab << PC1->decayPower <<  tab;
    if (!PC1->isPower)
    {
        outputPowerFile << ode.react << tab << inputReact << tab;
        if (PC1->isInputReactivity)
        {
            outputPowerFile << PC1->fuelReac << tab << PC1->emitterReac << tab << PC1->collectorReac << tab << PC1->moderatorReac << tab;
        }
    }
    outputPowerFile << endl;

    neutronTimeWrite += neutronTimeWriteInterval;
}

forAll(solidRegions, regionI)
{
    if (solidRegions[regionI].name().find("fuel") != std::string::npos) // 说明是燃料
    {
        powerFactorSolid[regionI] = PC1->powerRelative;
    }
}

// FatalError << "test This" << exit(FatalError);

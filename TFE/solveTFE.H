// 若存在TEC元件，规定间隔时间计算一次
if(runTime.value() > TFEcalcTime)
{
    Info << "Time = " << runTime.value() << nl;

    if (!emitterIndex.empty())
    {
        // 串联电路计算
        C1->circuitTECsCalc();
        Info << "Calculation of TECs is over!" << nl;
        // std::getchar();
        // 并联电路计算
        C1->parallelTECsCalc();
    }

    // 对电子流动造成的导热量进行估计和赋值（对应发射极外表面）
    for (int i = 0; i < int(emitterIndex.size()); ++i)
    {
        label regionI = emitterIndex[i];

        volScalarField& heatElectronE = addHeatFluxSolid[regionI];

        for (int j = 0; j < int(outerFaceId[regionI].size()); ++j)
        {
            for (int k = 0; k < int(outerFaceId[regionI][j].size()); ++k)
            {
                heatElectronE.boundaryFieldRef()[outerBoundIndex[regionI]][outerFaceId[regionI][j][k]]
                    = C1->TECqeList[i][j];
            }
        }
    }
    
    // 对电子流动造成的导热量进行估计和赋值（对应接收极内表面）
    for (int i = 0; i < int(collectorIndex.size()); ++i)
    {
        label regionI = collectorIndex[i];

        volScalarField& heatElectronC = addHeatFluxSolid[regionI];

        for (int j = 0; j < int(innerFaceId[regionI].size()); ++j)
        {
            for (int k = 0; k < int(innerFaceId[regionI][j].size()); ++k)
            {
                heatElectronC.boundaryFieldRef()[innerBoundIndex[regionI]][innerFaceId[regionI][j][k]]
                    = -C1->TECqeList[i][j];
            }
        }
    }

    std::vector<Foam::scalar> jouleHeatSingleEmitterListSeries(C1->TECs.size(), 0.0);
    std::vector<Foam::scalar> jouleHeatSingleCollectorListSeries(C1->TECs.size(), 0.0);
    std::vector<Foam::scalar> jouleHeatSingleEmitterListParallel(C1->TECsParallel.size(), 0.0);
    std::vector<Foam::scalar> jouleHeatSingleCollectorListParallel(C1->TECsParallel.size(), 0.0);
    // 赋值每根元件发射极焦耳热功率分布，作为后续处理的源项
    for (int i = 0; i < int(emitterIndex.size()); ++i)
    {
        // 每根元件轴向焦耳热（分元件首先进行计算）
        std::vector<Foam::scalar> jouleHeatEmitterList(int(emitterCellId[i].size()), 0.0);
        // 并联TECs
        if (i < int(C1->TECsParallel.size() + 0.1))
        {
            // 若此路没有断
            if (C1->isNormalParallel[i])
            {
                // 计算每根元件焦耳热（轴向分布）
                for (int j = 0; j < int(C1->TECsParallel[i]->J.size()); ++j)
                {
                    jouleHeatEmitterList[j] = Foam::pow(C1->TECsParallel[i]->IEsecSingle[j], 2) * C1->TECsParallel[i]->resistanceE[j];
                }
            }
            else
            {
                // 若此路断
                for (int j = 0; j < int(C1->TECsParallel[i]->J.size()); ++j)
                {
                    jouleHeatEmitterList[j] = 0.0;
                }
            }            
        }
        // 串联TECs
        else
        {
            if (C1->isNormalSeries)
            {
                // 计算每根元件焦耳热（轴向分布）
                for (int j = 0; j < int(C1->TECs[i - 3]->J.size()); ++j)
                {
                    jouleHeatEmitterList[j] = Foam::pow(C1->TECs[i - 3]->IEsecSingle[j], 2) * C1->TECs[i - 3]->resistanceE[j];
                }
            }
            else
            {
                // 若此路断
                for (int j = 0; j < int(C1->TECs[i - 3]->J.size()); ++j)
                {
                    jouleHeatEmitterList[j] = 0.0;
                }
            }
        }

        // 创建焦耳热功率场引用
        volScalarField& jouleHeatEmitter = jouleHeatEmitterSolid[emitterIndex[i]];
        // 赋值焦耳热功率场分布
        for (int j = 0; j < int(emitterCellId[i].size()); ++j)
        {
            for (int k = 0; k < int(emitterCellId[i][j].size()); ++k)
            {
                jouleHeatEmitter[emitterCellId[i][j][k]] = jouleHeatEmitterList[j] / emitterCellId[i][j].size();
            }
        }

        // 计算每根元件发射极焦耳热总和，用于后续计算
        // 并联部分
        if (i < int(C1->TECsParallel.size() + 0.1))
        {
            for (int j = 0; j < int(C1->TECsParallel[i]->J.size()); ++j)
            {
                jouleHeatSingleEmitterListParallel[i] += jouleHeatEmitterList[j];
            }
        }
        // 串联部分
        else
        {
            for (int j = 0; j < int(C1->TECs[i - 3]->J.size()); ++j)
            {
                jouleHeatSingleEmitterListSeries[i - 3] += jouleHeatEmitterList[j];
            }
        }
    }
    

    // 赋值每根元件接收极焦耳热功率分布，作为后续处理的源项
    for (int i = 0; i < int(collectorIndex.size()); ++i)
    {
        // 每根元件轴向焦耳热（分元件首先进行计算）
        std::vector<Foam::scalar> jouleHeatCollectorList(int(collectorCellId[i].size()), 0.0);
        // 并联TECs
        if (i < int(C1->TECsParallel.size() + 0.1))
        {
            if (C1->isNormalParallel[i])
            {
                // 若此路没有断
                for (int j = 0; j < int(C1->TECsParallel[i]->J.size()); ++j)
                {
                    jouleHeatCollectorList[j] = Foam::pow(C1->TECsParallel[i]->ICsecSingle[j], 2) * C1->TECsParallel[i]->resistanceC[j];
                }
            }
            else
            {
                // 若此路断
                for (int j = 0; j < int(C1->TECsParallel[i]->J.size()); ++j)
                {
                    jouleHeatCollectorList[j] = 0.0;
                }
            }
        }
        // 串联TECs
        else
        {
            if (C1->isNormalSeries)
            {
                // 计算每根元件焦耳热（轴向分布）
                for (int j = 0; j < int(C1->TECs[i - 3]->J.size()); ++j)
                {
                    jouleHeatCollectorList[j] = Foam::pow(C1->TECs[i - 3]->ICsecSingle[j], 2) * C1->TECs[i - 3]->resistanceC[j];
                }
            }
            else
            {
                // 计算每根元件焦耳热（轴向分布）
                for (int j = 0; j < int(C1->TECs[i - 3]->J.size()); ++j)
                {
                    jouleHeatCollectorList[j] = 0.;
                }
            }
        }

        // 创建焦耳热功率场引用
        volScalarField& jouleHeatCollector = jouleHeatCollectorSolid[collectorIndex[i]];
        // 赋值焦耳热功率场分布
        for (int j = 0; j < int(collectorCellId[i].size()); ++j)
        {
            for (int k = 0; k < int(collectorCellId[i][j].size()); ++k)
            {
                jouleHeatCollector[collectorCellId[i][j][k]] = jouleHeatCollectorList[j] / collectorCellId[i][j].size();
            }
        }   

        // 计算每根元件接收极焦耳热总和，用于后续计算
        // 并联部分
        if (i < int(C1->TECsParallel.size() + 0.1))
        {
            for (int j = 0; j < int(C1->TECsParallel[i]->J.size()); ++j)
            {
                jouleHeatSingleCollectorListParallel[i] += jouleHeatCollectorList[j];
            }
        }
        // 串联部分
        else
        {
            for (int j = 0; j < int(C1->TECs[i - 3]->J.size()); ++j)
            {
                jouleHeatSingleCollectorListSeries[i - 3] += jouleHeatCollectorList[j];
            }
        }
    }

    // Info << "test1" << nl;

    double Iall = C1->Iout;
    double Uall = C1->Uout;
    // 总 串联 电输出功率
    Foam::scalar PoutSeries = Iall * Uall;
    std::vector<Foam::scalar> PoutListSeries(int(C1->TECs.size()), 0.0);
    // 计算每个串联元件的电功率(按照两端电压加权分配)
    // 除分配电功率外，还需要考虑焦耳热功率
    Foam::scalar UTFEall = 0.;
    if (C1->isNormalSeries)
    {
        for (int i = 0; i < int(C1->TECs.size()); ++i)
        {
            UTFEall += 0.5 * (C1->TECs[i]->terminalPointUE1 + C1->TECs[i]->terminalPointUE2
                - C1->TECs[i]->terminalPointUC1 - C1->TECs[i]->terminalPointUC2);
        }
        // Info << "UTFEall = " << UTFEall << endl;
        for (int i = 0; i < int(PoutListSeries.size()); ++i)
        {
            PoutListSeries[i] = 0.5 * (C1->TECs[i]->terminalPointUE1 + C1->TECs[i]->terminalPointUE2
                - C1->TECs[i]->terminalPointUC1 - C1->TECs[i]->terminalPointUC2) / UTFEall * PoutSeries;
    
            PoutListSeries[i] += (jouleHeatSingleEmitterListSeries[i] + jouleHeatSingleCollectorListSeries[i]);
    
            // Info << "\tPoutListSeries[" << i << "] = " << PoutListSeries[i] << 
            //     "\tjouleHeatSingleEmitterListSeries[" << i << "] = " << jouleHeatSingleEmitterListSeries[i] << 
            //     "\tjouleHeatSingleCollectorListSeries[" << i << "] = " << jouleHeatSingleCollectorListSeries[i] << endl;
        }
    }
    else
    {
        for (int i = 0; i < int(PoutListSeries.size()); ++i)
        {    
            PoutListSeries[i] = 0.;
        }
    }

    // Info << endl << endl;

    // Foam::sleep(5);

    std::vector<Foam::scalar> PoutListParallel(int(C1->TECsParallel.size()), 0.0);
    for (int i = 0; i < int(PoutListParallel.size()); ++i)
    {
        if (C1->isNormalParallel[i])
        {
            PoutListParallel[i] = C1->parallelVoltage * C1->TECsParallel[i]->I;
            PoutListParallel[i] += (jouleHeatSingleEmitterListParallel[i] + jouleHeatSingleCollectorListParallel[i]);
        }
        else
        {
            PoutListParallel[i] = 0.;
        }
    }

    // Foam::sleep(5);

    // 赋值每根元件电功率分布，作为后续处理的源项
    for (int i = 0; i < int(collectorIndex.size()); ++i)
    {       
        // 电流密度轴向分配
        // 计算总电流密度
        Foam::scalar Jsum = 0.;
        if (i < int(C1->TECsParallel.size() + 0.1))
        {
            if (C1->isNormalParallel[i])
            {
                for (int j = 0; j < int(C1->TECsParallel[i]->J.size()); ++j)
                {
                    Jsum += C1->TECsParallel[i]->J[j];
                }
            }
            else
            {
                Jsum = 1.;
            }
            
        }
        else
        {
            if (C1->isNormalSeries)
            {
                for (int j = 0; j < int(C1->TECs[i - 3]->J.size()); ++j)
                {
                    Jsum += C1->TECs[i - 3]->J[j];
                }
            }
            else
            {
                Jsum = 1.;
            }
        }

        std::vector<Foam::scalar> Plist(int(collectorCellId[i].size()), 0.0);
        // 并联TECs
        if (i < int(C1->TECsParallel.size() + 0.1))
        {
            // 读取单根功率
            Foam::scalar Psingle = PoutListParallel[i];
             // 将总功率按照电流密度J轴向分配
             for (int j = 0; j < int(C1->TECsParallel[i]->J.size()); ++j)
            {
                Plist[j] = Psingle * C1->TECsParallel[i]->J[j] / Jsum;
            }
        }
        // 串联TECs
        else
        {
            // 读取单根功率
            Foam::scalar Psingle = PoutListSeries[i - 3];
             // 将总功率按照电流密度J轴向分配
             for (int j = 0; j < int(C1->TECs[i - 3]->J.size()); ++j)
            {
                Plist[j] = Psingle * C1->TECs[i - 3]->J[j] / Jsum;
            }
        }
        // 创建电功率场引用
        volScalarField& electricPowerC = electricPowerSolid[collectorIndex[i]];
        // 赋值功率场分布
        for (int j = 0; j < int(collectorCellId[i].size()); ++j)
        {
            for (int k = 0; k < int(collectorCellId[i][j].size()); ++k)
            {
                electricPowerC[collectorCellId[i][j][k]] = Plist[j] / collectorCellId[i][j].size();
            }
        }
    }

    outputFile << runTime.value() << tab << Iall << tab << Uall << endl;
    for (int i = 0; i < int(C1->TECs.size()); ++i){
        const double U1 = C1->TECs[i]->wireU[0] - C1->TECs[i]->wireU[2];
        outputFile << tab << U1;
    }
    outputFile << endl;

    for (int i = 0; i < int(C1->TECs.size()); ++i){
        const double U2 = C1->TECs[i]->wireU[1] - C1->TECs[i]->wireU[3];
        outputFile << tab << U2;
    }
    outputFile << endl;

    for (int i = 0; i < int(C1->TECs.size()); ++i){
        const double I1 = C1->Ileft[i];
        outputFile << tab << I1;
    }
    outputFile << endl;

    Info << "outputFile.size() = " << outputFiles.size() << endl;
    for(int i = 0; i < int(outputFiles.size()); ++i){
        Info << "i = " << i << "\tnumberEmitter = " << numberEmitter[i] << endl;
        outputFiles[i]() << runTime.value() << tab;
        // outputFiles[i]() << tab << C1->deltaU1[i] << tab << C1->deltaU2[i] << endl;
        // 并联需要额外输出电流
        if (i < int(C1->TECsParallel.size() + 0.1))
        {
            outputFiles[i]() << C1->TECsParallel[i]->I << tab << C1->parallelVoltage << endl;
        }
        // 串联输出的是两端电势差，用作对比
        else
        {
            outputFiles[i]() << tab << C1->deltaU1[i - 3] << tab << C1->deltaU2[i - 3] << endl;
            for(int j =  0; j < int(C1->TECs[i - 3]->J.size()); ++j){
                outputFiles[i]() << C1->TECs[i - 3]->J[j] << " ";
            }
            outputFiles[i]() << endl;
            for(int j = 0; j < int(C1->TECs[i - 3]->J.size()); ++j){
                outputFiles[i]() << C1->TECs[i - 3]->UE[j] << " ";
            }
            outputFiles[i]() << endl;
            for(int j = 0; j < int(C1->TECs[i - 3]->J.size()); ++j){
                outputFiles[i]() << C1->TECs[i - 3]->UC[j] << " ";
            }
            outputFiles[i]() << endl;
        }
    }

    TFEcalcTime += timeStep;

    Info << "TFEcalcTime = " << TFEcalcTime << nl;
    // FatalError << "Set for testing" << exit(FatalError);
}


// std::vector<std::vector<double>> testInput;
// for(int i = 0; i < int(emitterIndex.size()); ++i)
// if(int(emitterIndex.size()) == 1){
//     testInput.push_back(TElist[0]);
//     testInput.push_back(TClist[0]);
//     std::vector<double> dlC(int(TElist[0].size()), 0.00754);
//     std::vector<double> dlE(dlC);
//     testInput.push_back(dlC);
//     testInput.push_back(dlE);
//     std::vector<double> cross;
//     cross.push_back(6.6669999999999997e-05);
//     cross.push_back(0.00010786000000000000);
//     testInput.push_back(cross);
//     std::vector<double> sideAreaE(int(TElist[0].size()), 0.00092855424159680002 / 2.);
//     std::vector<double> sideAreaC(int(TElist[0].size()), 0.00097592945800480005 / 2.);
//     testInput.push_back(sideAreaE);
//     testInput.push_back(sideAreaC);
//     std::vector<double> REC;
//     // REC.push_back(re0Solid[0]);
//     // REC.push_back(re1Solid[0]);
//     // REC.push_back(rc0Solid[0]);
//     // REC.push_back(rc1Solid[0]);
//     REC.push_back(0.00097 * 1.8);
//     REC.push_back(0.00064 * 1.8);
//     REC.push_back(0.00021 * 1.8);
//     REC.push_back(0.00038 * 1.8);
//     testInput.push_back(REC);
//     std::vector<double> other = { currentSolid[emitterIndex[0]], 0.5 };
//     testInput.push_back(other);
//     std::vector<double> tcs(int(TElist[0].size()), 600.);
//     testInput.push_back(tcs);
//     TFE* testT = new TFE(testInput);
//     testT->initial();
//     testT->calc();
//     currentSolid[emitterIndex[0]] = testT->I;
//     voltageSolid[emitterIndex[0]] = testT->U;
//     Info<< testT->I << "\t\t" << testT->U << nl;
// }
// // else{
// //     FatalError << "Not added here" << exit(FatalError);
// // }

#include <unistd.h>

// Type of the solid regions
List<scalar> electrodSolid(solidRegions.size());

// Electric parameters of TFEs
List<scalar> currentSolid(solidRegions.size());
List<scalar> voltageSolid(solidRegions.size());
List<scalar> powerSolid(solidRegions.size());
List<scalar> re0Solid(solidRegions.size());
List<scalar> re1Solid(solidRegions.size());
List<scalar> rc0Solid(solidRegions.size());
List<scalar> rc1Solid(solidRegions.size());
List<scalar> cesiumTemperatureSolid(solidRegions.size());

// Electric parameters of nodes
PtrList<volScalarField> currentDensitySolid(solidRegions.size());
PtrList<volScalarField> electricPotentialSolid(solidRegions.size());
PtrList<volScalarField> thermoElectricSolid(solidRegions.size());
PtrList<volScalarField> electricPowerSolid(solidRegions.size());
PtrList<volScalarField> jouleHeatEmitterSolid(solidRegions.size());
PtrList<volScalarField> jouleHeatCollectorSolid(solidRegions.size());
PtrList<volScalarField> jouleHeatCollectorForEmitterSolid(solidRegions.size()); // 目前没有使用 

// 几何参数注册
PtrList<IOdictionary> geometricParameterSolid(solidRegions.size());

forAll(solidRegions, i){
    Info<< "    Reading the type of Geometry\n" << endl;

    geometricParameterSolid.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "geometricParameter",
                runTime.constant(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    IOdictionary& geometricParameter = geometricParameterSolid[i];

    // IOdictionary geometricParameter
    // (
    //     IOobject
    //     (
    //         "geometricParameter",
    //         runTime.constant(),
    //         solidRegions[i],
    //         IOobject::MUST_READ,
    //         IOobject::NO_WRITE
    //     )
    // );

    // Read the set eletrodSolid;
    // > 0 for emitter;
    // < 0 for collector;
    // Number means the ID (order of connection) of TECs
    // 0 for other units
    Info<< "    Adding to electrodSolid\n" << endl;
    dimensionedScalar electrodTemp
    (
        "electrod",
        dimless,
        geometricParameter.lookup("electrod")
    );
    electrodSolid[i] = electrodTemp.value();

    // Info << electrodSolid[i] << nl;

    // Foam::sleep(3);

    if(electrodSolid[i] != 0){
        Info<< "    Adding to currentSolid\n" << endl;
        dimensionedScalar currentTemp
        (
            "current",
            dimless,
            geometricParameter.lookup("current")
        );
        currentSolid[i] = currentTemp.value();

        Info<< "    Adding to voltageSolid\n" << endl;
        dimensionedScalar voltageTemp
        (
            "voltage",
            dimless,
            geometricParameter.lookup("voltage")
        );
        voltageSolid[i] = voltageTemp.value();

        Info<< "    Adding to powerSolid\n" << endl;
        dimensionedScalar powerTemp
        (
            "power",
            dimless,
            geometricParameter.lookup("power")
        );
        powerSolid[i] = powerTemp.value();

        Info<< "    Adding to re0Solid\n" << endl;
        dimensionedScalar re0Temp
        (
            "re0",
            dimless,
            geometricParameter.lookup("re0")
        );
        re0Solid[i] = re0Temp.value();

        Info<< "    Adding to re1Solid\n" << endl;
        dimensionedScalar re1Temp
        (
            "re1",
            dimless,
            geometricParameter.lookup("re1")
        );
        re1Solid[i] = re1Temp.value();

        Info<< "    Adding to rc0Solid\n" << endl;
        dimensionedScalar rc0Temp
        (
            "rc0",
            dimless,
            geometricParameter.lookup("rc0")
        );
        rc0Solid[i] = rc0Temp.value();

        Info<< "    Adding to rc1Solid\n" << endl;
        dimensionedScalar rc1Temp
        (
            "rc1",
            dimless,
            geometricParameter.lookup("rc1")
        );
        rc1Solid[i] = rc1Temp.value();

        Info<< "    Adding to cesiumTemperatureSolid\n" << endl;
        dimensionedScalar cesiumTemperatureTemp
        (
            "cesiumTemperature",
            dimless,
            geometricParameter.lookup("cesiumTemperature")
        );
        cesiumTemperatureSolid[i] = cesiumTemperatureTemp.value();

        Info<< "    Adding to currentDensitySolid\n" << endl;
        currentDensitySolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "currentDensity",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar(dimCurrent, 0.0)
            )
        );

        Info<< "    Adding to electricPotentialSolid\n" << endl;
        electricPotentialSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "electricPotential",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar
                (
                    dimPower / dimCurrent,
                    0.0
                )
            )
        );
        
        Info<< "    Adding to thermoElectricSolid\n" << endl;
        thermoElectricSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "thermoElectric",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar(dimPower, 0.0)
            )
        );

        // 注册电功率场
        Info<< "    Adding to electricPowerSolid\n" << endl;
        electricPowerSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "electricPower",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar(dimPower, 0.0)
            )
        );

        // 注册发射极焦耳热场
        Info<< "    Adding to jouleHeatEmitterSolid\n" << endl;
        jouleHeatEmitterSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "jouleHeatEmitter",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar(dimPower, 0.0)
            )
        );

        // 注册接收极焦耳热场
        Info<< "    Adding to jouleHeatCollectorSolid\n" << endl;
        jouleHeatCollectorSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "jouleHeatCollector",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar(dimPower, 0.0)
            )
        );

        // 注册发射极焦耳热场（用于赋给接收极，减去发热部分）（目前没有使用）
        Info<< "    Adding to jouleHeatCollectorForEmitterSolid\n" << endl;
        jouleHeatCollectorForEmitterSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "jouleHeatCollectorForEmitter",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar(dimPower, 0.0)
            )
        );
    }
}


// Construct the index of emitters and collectors
std::vector<int> numberEmitter;
std::vector<int> numberCollector;
std::vector<int> emitterIndex;
std::vector<int> collectorIndex;
forAll(solidRegions, i)
{
    // Judged as emitter
    if (electrodSolid[i] > 0.01)
    {
        emitterIndex.push_back(i);
        // Info<< "i = " << i << ", region name is "
        //     << solidRegions[i].name() << tab << emitterIndex.back() << endl;
    }
    if (electrodSolid[i] < -0.01)
    {
        collectorIndex.push_back(i);
        // Info<< "i = " << i << ", region name is "
        //     << solidRegions[i].name() << tab << collectorIndex.back()  << endl;
    }
}

// Construct  面索引
std::vector<std::vector<std::vector<int>>> emitterFaceId(int(emitterIndex.size()));
std::vector<std::vector<std::vector<int>>> collectorFaceId(int(collectorIndex.size()));
// 体索引----发射极+接收极
std::vector<std::vector<std::vector<int>>> emitterCellId(int(emitterIndex.size()));
std::vector<std::vector<std::vector<int>>> collectorCellId(int(collectorIndex.size()));


std::vector<int> emitterCoupledInterfaceIndex;
std::vector<int> collectorCoupledInterfaceIndex;

if(!(collectorIndex.empty() && emitterIndex.empty()))
{
    if(collectorIndex.size() != emitterIndex.size())
    {
        FatalError << "The number of emitter and receiver stages does not match" << exit(FatalError);
    }

    // 改为在constant/.../geometricParameter中读取连接编号
    // Acquire every emitter's TFE id
    for(int i = 0; i < int(emitterIndex.size()); ++i)
    {
        // Info << solidRegions[emitterIndex[i]].name() << tab << extractSingleIntFromString(solidRegions[emitterIndex[i]].name()) << nl;
        numberEmitter.push_back(int(0.5 + electrodSolid[emitterIndex[i]]));
        // Info << "Emitter, numberEmitter = " << numberEmitter.back() << tab << solidRegions[emitterIndex[i]].name() << nl;
        // numberEmitter.push_back(extractSingleIntFromString(solidRegions[emitterIndex[i]].name()));
        // Info << numberEmitter.back() << nl;
    }
    // Foam::sleep(60);
    // Acquire every collector's TFE id
    for(int i = 0; i < int(collectorIndex.size()); ++i)
    {
        // Info << solidRegions[collectorIndex[i]].name() << tab << extractSingleIntFromString(solidRegions[collectorIndex[i]].name()) << nl;
        numberCollector.push_back(int(-0.5 + electrodSolid[collectorIndex[i]]));
        // Info << "Collector, numberCollector = " << numberCollector.back() << tab << solidRegions[collectorIndex[i]].name() << nl;
        // numberCollector.push_back(extractSingleIntFromString(solidRegions[collectorIndex[i]].name()));
    }
    // Info << "numberEmitter = " << numberEmitter.size() << tab;
    // Info << "emitterIndex = " << emitterIndex.size() << nl;
    // Foam::sleep(60);

    // Sort emitter's id
    bool isNeg = false;
    for(int num : numberEmitter)
    {
        // Info << "Emitter, num = " << num << nl;
        if(num == -1)
        {
            isNeg = true;
            break;
        }
    }
    if(isNeg)
    {  // Only one TFE is included
        if(numberEmitter.size() != 1)
        {
            FatalError << "Incorrect TFE number" << exit(FatalError);
        }
    }
    else
    {
        sort_with_correspondence(numberEmitter, emitterIndex);
        // for(int i = 0; i < int(emitterIndex.size()); ++i){
        //     Info << numberEmitter[i] << tab << emitterIndex[i] << endl;
        // }
    }

    // Sort collector's id
    isNeg = false;
    for(int num : numberCollector)
    {
        if(num == -1)
        {
            isNeg = true;
            break;
        }
    }
    if(isNeg)
    {
        if(numberEmitter.size() > 1)
        {
            FatalError << "Incorrect TFE number" << exit(FatalError);
        }
    }
    else 
    {
        sort_with_correspondence(numberCollector, collectorIndex);
        // for(int i = 0; i < int(collectorIndex.size()); ++i){
        //     Info << numberCollector[i] << tab << collectorIndex[i] << endl;
        // }
        std::reverse(collectorIndex.begin(), collectorIndex.end());
    }

    // Create interface index of emitter and collector
    for(int i = 0; i < int(emitterIndex.size()); ++i)
    {
        std::string coupledNameE = solidRegions[emitterIndex[i]].name() + "_to_" + 
            solidRegions[collectorIndex[i]].name();
        // Info<< "region name = " << solidRegions[emitterIndex[i]].name() << tab << "coupledNameE = " << coupledNameE << endl;
        const fvBoundaryMesh& boundaryE = solidRegions[emitterIndex[i]].boundary();
        forAll(boundaryE, patchI)
        {
            const fvPatch& patch = boundaryE[patchI];
            std::string patchName = std::string(patch.name());
            // Info<< "patchI = " << patchI << ", " << "patchName = " << patchName << endl;
            if (patchName == coupledNameE)
            {
                emitterCoupledInterfaceIndex.push_back(patchI);
                // Info << solidRegions[emitterIndex[i]].name() << ": " << patchName << endl;
            }
        }
        std::string coupledNameC = solidRegions[collectorIndex[i]].name() + "_to_" + 
            solidRegions[emitterIndex[i]].name();
        const fvBoundaryMesh& boundaryC = solidRegions[collectorIndex[i]].boundary();
        // Info<< "region name = " << solidRegions[collectorIndex[i]].name() << tab << "coupledNameE = " << coupledNameC << endl;
        forAll(boundaryC, patchI)
        {
            const fvPatch& patch = boundaryC[patchI];
            std::string patchName = std::string(patch.name());
            //Info<< "patchI = " << patchI << ", " << "patchName = " << patchName << endl;
            if (patchName == coupledNameC)
            {
                collectorCoupledInterfaceIndex.push_back(patchI);
                Info << solidRegions[collectorIndex[i]].name() << ": " << patchName << endl;
            }
        }
    }

    if(collectorCoupledInterfaceIndex.empty() || emitterCoupledInterfaceIndex.empty())
    {
        FatalError << "No coupling boundary was found" << exit(FatalError);
    }

    // 记录下指定区域和指定边界的索引及对应的排序
    for(int i = 0; i < int(emitterIndex.size()); ++i)
    {
        std::vector<int> faceZIndexE;
        std::vector<double> faceZE;
        const vectorField& CfE = solidRegions[emitterIndex[i]].boundary()[emitterCoupledInterfaceIndex[i]].Cf();
        forAll(CfE, faceI)
        {
            faceZIndexE.push_back(faceI);
            faceZE.push_back(CfE[faceI].z());
        }
        emitterFaceId[i] = processVectors(faceZE, faceZIndexE);
        // for(int j = 0; j < int(emitterFaceId[i].size()); ++j){
        //     for(int k = 0; k < int(emitterFaceId[i][j].size()); ++k){
        //         const label faceI = emitterFaceId[i][j][k];
        //         Info << faceI << "\t" << CfE[faceI].z() << nl;
        //     }
        //     Info << nl;taogh
        // }

        std::vector<int> faceZIndexC;
        std::vector<double> faceZC;
        const vectorField& CfC = solidRegions[collectorIndex[i]].boundary()[collectorCoupledInterfaceIndex[i]].Cf();
        forAll(CfC, faceI){
            faceZIndexC.push_back(faceI);
            faceZC.push_back(CfC[faceI].z());
        }
        collectorFaceId[i] = processVectors(faceZC, faceZIndexC);
        // for(int j = 0; j < int(collectorFaceId[i].size()); ++j){
        //     for(int k = 0; k < int(collectorFaceId[i][j].size()); ++k){
        //         const label faceI = collectorFaceId[i][j][k];
        //         Info << faceI << "\t" << CfC[faceI].z() << nl;
        //     }
        //     Info << nl;
        // }
    }

    // 记录下指定区域的体网格索引及排序--发射极
    for (int i = 0; i < int(emitterIndex.size()); ++i)
    {
        std::vector<int> cellZIndexE;
        std::vector<double> cellZE;
        const vectorField& cellE = solidRegions[emitterIndex[i]].C();
        forAll(cellE, cellI)
        {
            cellZIndexE.push_back(cellI);
            cellZE.push_back(cellE[cellI].z());
        }
        emitterCellId[i] = processVectors(cellZE, cellZIndexE);
    }

    // 记录下指定区域的体网格索引及排序--接收极
    for (int i = 0; i < int(collectorIndex.size()); ++i)
    {
        std::vector<int> cellZIndexC;
        std::vector<double> cellZC;
        const vectorField& cellC = solidRegions[collectorIndex[i]].C();
        forAll(cellC, cellI)
        {
            cellZIndexC.push_back(cellI);
            cellZC.push_back(cellC[cellI].z());
        }
        collectorCellId[i] = processVectors(cellZC, cellZIndexC);
    }

    // FatalError << "Not added here" << exit(FatalError);
}
else
{
    Info<< "No emitter/collector found!" << endl;
}

Info << "Create Temperature lists of emitter and collector" << endl;
std::vector<std::vector<double>> TElist(int(emitterIndex.size()));
std::vector<std::vector<double>> TClist(int(collectorFaceId.size()));

for(int i = 0; i < int(TElist.size()); ++i)
{
    TElist[i].resize(int(emitterFaceId[0].size()));
}
for(int i = 0; i < int(TClist.size()); ++i)
{
    TClist[i].resize(int(collectorFaceId[0].size()));
}

// Info<<int(collectorFaceId[0].size()) << "\t\t" << int(emitterFaceId[0].size()) <<endl;

Info << "Create list of TECs" << endl;
// 创建TEC元件列表并进行初始化
std::vector<std::singleThermionicEnergyConversion*> thermionicEnergyConversions;
std::vector<std::singleThermionicEnergyConversion*> thermionicEnergyConversionsParallel;

// 初始化电路中每个TEC元件
// 目前的设置是：
// 前三个（2-4）是并联元件
// 后面的34根（5-38）是串联元件
for(int i = 0; i < int(emitterIndex.size()); ++i)
{
    std::vector<std::vector<double>> input;
    std::vector<double> TE(int(emitterFaceId[i].size()), 743.);
    std::vector<double> TC(int(emitterFaceId[i].size()), 743.);
    std::vector<double> dlE(int(emitterFaceId[i].size()), 0.377 / int(emitterFaceId[i].size()));
    std::vector<double> dlC(dlE);
    std::vector<double> cross = { 6.6669999999999997e-05, 0.00010786000000000000 };
    std::vector<double> sideAreaE(int(emitterFaceId[i].size()), 0.0232 / int(emitterFaceId[i].size()));
    std::vector<double> sideAreaC(int(emitterFaceId[i].size()), 0.0244 / int(emitterFaceId[i].size()));
    std::vector<double> REC;
    REC.push_back(re0Solid[emitterIndex[i]]);
    REC.push_back(re1Solid[emitterIndex[i]]);
    REC.push_back(rc0Solid[emitterIndex[i]]);
    REC.push_back(rc1Solid[emitterIndex[i]]);
    std::vector<double> other = {0.9, 0.5};
    std::vector<double> tcs(int(emitterFaceId[i].size()), cesiumTemperatureSolid[emitterIndex[i]]);
    // Info << "tcs = " << cesiumTemperatureSolid[emitterIndex[i]] << endl;
    // Foam::sleep(1);
    std::vector<double> vinput(int(emitterFaceId[i].size()), 0.9);
    std::vector<double> itat = {200.};
    std::vector<double> uwire = {0.8, 0.8, 0., 0.};

    input.push_back(TE);
    input.push_back(TC);
    input.push_back(dlE);
    input.push_back(dlC);
    input.push_back(cross);
    input.push_back(sideAreaE);
    input.push_back(sideAreaC);
    input.push_back(REC);
    input.push_back(other);
    input.push_back(tcs);
    input.push_back(vinput);
    input.push_back(itat);
    input.push_back(uwire);

    std::singleThermionicEnergyConversion* S1 = new std::singleThermionicEnergyConversion(input);

    // Info << S1->crossAreaC << tab << S1->crossAreaE << endl;
    // Info << S1->rhoC[0] << tab << S1->rhoE[0] << tab << S1->dlC[0] << tab << S1->dlE[0] << nl;

    S1->initial();

    if (numberEmitter[i] <= 4.5){
        thermionicEnergyConversionsParallel.push_back(S1);
    }
    else{
        thermionicEnergyConversions.push_back(S1);
    }
}

// 初始化电路
Info << "Initialize the circuit" << endl;
// Info << "emitterIndex[0] is " << solidRegions[emitterIndex[0]].name() << nl;
// Info << "emitterIndex[3] is " << solidRegions[emitterIndex[0]].name() << nl;
// Foam::sleep(60);
std::circuitTECs* C1 = new std::circuitTECs();
C1->TECs = thermionicEnergyConversions;
C1->TECsParallel = thermionicEnergyConversionsParallel;
C1->isFixedU = true;
C1->isFixedR = false;
C1->isFixedI = false;
C1->Utarget = voltageSolid[emitterIndex[3]];
C1->Uout = voltageSolid[emitterIndex[3]];
C1->parallelVoltage = voltageSolid[emitterIndex[0]];
// C1->Utarget = 30.00;
// C1->Uout = 30.00;
// 电流计算初值
C1->Iout = 198.63;
C1->isFirst = true;
// 暂定电阻，按需修改
C1->Rload = 0.145;

// Create the outputFile for electric
Info << "Create the output Files" << endl;
OFstream outputFile(runTime.rootPath()/runTime.caseName()/"IU.dat");
outputFile << "Time" << tab << "I" <<  tab << "U" << endl;

std::vector<autoPtr<OFstream>> outputFiles;
for (int i = 0; i < int(numberEmitter.size()); ++i)
{
    std::string fileName = "JandUEUC" + std::to_string(i) + ".dat";
    outputFiles.emplace_back
    (
        new OFstream(runTime.rootPath()/runTime.caseName()/fileName)
    );

    // 写入文件头
    outputFiles[i]() << "Time" << tab;
    outputFiles[i]() << solidRegions[collectorIndex[i]].name() << endl;
}

scalar time0 = runTime.value();
// TFE计算时间间隔为2.0s
scalar timeStep = 2.0;
scalar TFEcalcTime = time0;

// std::vector<double> 

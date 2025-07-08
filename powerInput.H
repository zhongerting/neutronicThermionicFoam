powerControl* PC1 = new powerControl();

std::ifstream file("inputControl.DAT");

if (!file.is_open()) 
{
    FatalError << "'inputControl.DAT' can not open" << exit(FatalError);
}

// 读取第一行内容（功率控制方式）
std::string controlLine;
std::getline(file, controlLine);
Info << "控制方式为: " << controlLine << nl;

// 读取第二行内容（表格列数）
int columnNum;
std::string columnNumStr;
std::getline(file, columnNumStr);
columnNum = std::stoi(columnNumStr);
Info << "表格列数为: " << columnNum << nl;

// 读取第三行内容（时间序列，以空格分割）
std::string timeLine;
std::getline(file, timeLine);
std::istringstream timeLineStream(timeLine);
std::vector<Foam::scalar> timeLineArray;
Foam::scalar num;
while (timeLineStream >> num)
{
    timeLineArray.push_back(num);
}
Info << "时间序列为: ";
for (Foam::scalar n : timeLineArray)
{
    Info << n << tab;
}
Info << nl;

// 读取第四行内容（输入序列，以空格分隔）
std::string inputLine;
std::getline(file, inputLine);
std::istringstream inputLineStream(inputLine);
std::vector<Foam::scalar> inputLineArray;
while (inputLineStream >> num)
{
    inputLineArray.push_back(num);
}
Info << "输入序列为: ";
for (Foam::scalar n : inputLineArray)
{
    Info << n << tab;
}
Info << nl;

// 创建常微分方程并初始化
neutronODE ode;
dictionary dict;
dict.add("solver", args[1]);
// Create the selected ODE system solver
autoPtr<ODESolver> odeSolver = ODESolver::New(ode, dict);

// 创建初始值并进行初始化（本求解器适用于初始功率为115kW的所有工况）
scalar initialPower = 108.062e3;
scalarField yStart(ode.nEqns());
yStart[0] = initialPower * (1.0 - ode.decayTotal);
yStart[1] = initialPower * ode.betaFra[0] / ode.neulife / ode.betaCon[0];
yStart[2] = initialPower * ode.betaFra[1] / ode.neulife / ode.betaCon[1];
yStart[3] = initialPower * ode.betaFra[2] / ode.neulife / ode.betaCon[2];
yStart[4] = initialPower * ode.betaFra[3] / ode.neulife / ode.betaCon[3];
yStart[5] = initialPower * ode.betaFra[4] / ode.neulife / ode.betaCon[4];
yStart[6] = initialPower * ode.betaFra[5] / ode.neulife / ode.betaCon[5];
yStart[7] = initialPower * ode.decayFra[0] / ode.decayCon[0];
yStart[8] = initialPower * ode.decayFra[1] / ode.decayCon[1];
yStart[9] = initialPower * ode.decayFra[2] / ode.decayCon[2];
yStart[10] = initialPower * ode.decayFra[3] / ode.decayCon[3];
scalarField y(yStart);
scalarField dyStart(ode.nEqns());
scalar dxEst = 0.5;

// 给PC1中的变量赋值
std::string modePC = PC1->getPower(controlLine);
PC1->timeInput = timeLineArray;
PC1->dataInput = inputLineArray;
PC1->powerRelative = 1.0;
PC1->powerRelativeOld = PC1->powerRelative;

PC1->isFirst = true;

// 如果采用功率控制
if (PC1->isPower)
{
    // 初始化最初功率
    Info << "使用功率控制" << nl;
    PC1->power = PC1->dataInput[0];
    PC1->powerOld = PC1->power;
    PC1->powerOldOld = PC1->powerOld;

    PC1->powerInitial = PC1->power;
}

// 如果采用反应性控制，则创建常微分方程
if (!PC1->isPower)
{
    if (PC1->isTotalReactivity)
    {
        Info << "使用反应性控制（总反应性）" << nl;
    }
    if (PC1->isInputReactivity)
    {
        Info << "使用反应性控制（外界反应性）" << nl;
    }

    // 初始化计算稳态结果作为输入
    ode.react = 0.0;
    odeSolver->solve(0.0, 1.0, y, 0, dxEst);

    Info << "y: " << y << nl;

    // 初始化功率
    PC1->fissPower = y[0];
    PC1->decayPower = 0.0;
    for (int i = 0; i < int(ode.decayFra.size()); i++)
    {
        PC1->decayPower += y[7 + i] * ode.decayCon[i];
    }
    PC1->power = PC1->fissPower + PC1->decayPower;
    PC1->powerOld = PC1->power;
    PC1->powerOldOld = PC1->powerOld;

    PC1->powerInitial = PC1->power;

    Info << PC1->fissPower << tab << PC1->decayPower << tab << PC1->power << nl;
}

// 创建输出文件
Info << "Create the power output Files" << endl;
OFstream outputPowerFile(runTime.rootPath()/runTime.caseName()/"POWEROUTPUT.dat");
outputPowerFile << "Time" << tab << "power" <<  tab << "fission power" <<  tab << "decay power" <<  tab << "react"
 <<  tab << "input react" <<  tab << "fuel react" <<  tab << "emitter react" <<  tab << "collector react" 
  <<  tab << "moderator react" << endl;

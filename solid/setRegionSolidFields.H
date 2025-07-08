const fvMesh& mesh = solidRegions[i];

solidThermo& thermo = thermoSolid[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

volScalarField& e = thermo.he();

const Foam::fvModels& fvModels = fvModelsSolid[i];
Foam::fvConstraints& fvConstraints = fvConstraintsSolid[i];

#include "checkRadiationModel.H"

solidNoLoopControl& pimple = pimples.solid(i);

// 声明自定义边界温度场的引用
volScalarField& boundaryT = boundTemperatureSolid[i];

// 找到内壁面温度分布的引用
if (innerBoundIndex[i] >= 0)
{
    const scalarField& patchTInner = thermoSolid[i].T().boundaryField()[innerBoundIndex[i]];
    for(int j = 0; j < int(innerFaceId[i].size()); ++j)
    {
        // 计算内表面平均轴向温度分布，并进行赋值
        double TempInner = 0.;
        for (int k = 0; k < int(innerFaceId[i][j].size()); ++k)
        {
            TempInner += patchTInner[innerFaceId[i][j][k]];
        }
        for (int k = 0; k < int(innerFaceId[i][j].size()); ++k){
            boundaryT.boundaryFieldRef()[innerBoundIndex[i]][innerFaceId[i][j][k]] = TempInner / int(innerFaceId[i][j].size());
        }
    }
}

// 找到外壁面温度分布的引用
if (outerBoundIndex[i] >= 0)
{
    const scalarField& patchTOuter = thermoSolid[i].T().boundaryField()[outerBoundIndex[i]];
    for(int j = 0; j < int(outerFaceId[i].size()); ++j)
    {
        // 计算外表面平均轴向温度分布，并进行赋值
        double TempOuter = 0.;
        for (int k = 0; k < int(outerFaceId[i][j].size()); ++k)
        {
            TempOuter += patchTOuter[outerFaceId[i][j][k]];
        }
        for (int k = 0; k < int(outerFaceId[i][j].size()); ++k){
            boundaryT.boundaryFieldRef()[outerBoundIndex[i]][outerFaceId[i][j][k]] = TempOuter / int(outerFaceId[i][j].size());
        }
    }
}


Foam::scalar tFuelAve = 0.;
Foam::scalar VFuelAve = 0.;
Foam::scalar tEmitterAve = 0.;
Foam::scalar VEmitterAve = 0.;
Foam::scalar tCollectorAve = 0.;
Foam::scalar VCollectorAve = 0.;
Foam::scalar tModeratorAve = 0.;
Foam::scalar VModeratorAve = 0.;

forAll(solidRegions, regionI)
{
    if (solidRegions[regionI].name().find("fuel") != std::string::npos) // 说明是燃料
    {
        volScalarField& T = thermoSolid[regionI].T();
        const DimensionedField<scalar, volMesh>& V = solidRegions[regionI].V();
        forAll(T, cellI)
        {
            tFuelAve += T[cellI] * V[cellI];
            VFuelAve += V[cellI]; 
        }
    }
    else if (solidRegions[regionI].name().find("emitter") != std::string::npos) // 说明是燃料
    {
        volScalarField& T = thermoSolid[regionI].T();
        const DimensionedField<scalar, volMesh>& V = solidRegions[regionI].V();
        forAll(T, cellI)
        {
            tEmitterAve += T[cellI] * V[cellI];
            VEmitterAve += V[cellI];
        }
    }
    else if (solidRegions[regionI].name().find("collector") != std::string::npos) // 说明是燃料
    {
        volScalarField& T = thermoSolid[regionI].T();
        const DimensionedField<scalar, volMesh>& V = solidRegions[regionI].V();
        forAll(T, cellI)
        {
            tCollectorAve += T[cellI] * V[cellI];
            VCollectorAve += V[cellI];
        }
    }
    else if (solidRegions[regionI].name().find("moderator_moderator") != std::string::npos) // 说明是燃料
    {
        volScalarField& T = thermoSolid[regionI].T();
        const DimensionedField<scalar, volMesh>& V = solidRegions[regionI].V();
        forAll(T, cellI)
        {
            tModeratorAve += T[cellI] * V[cellI];
            VModeratorAve += V[cellI];
        }
    }
    else
    {
        continue;
    }
}

tFuelAve /= VFuelAve;
tEmitterAve /= VEmitterAve;
tCollectorAve /= VCollectorAve;
tModeratorAve /= VModeratorAve;

ode.aveTemperature[0] = tFuelAve;
ode.aveTemperature[1] = tEmitterAve;
ode.aveTemperature[2] = tCollectorAve;
ode.aveTemperature[3] = tModeratorAve;

// 获取反应性反馈
PC1->reactivityCalc(ode.aveTemperature);

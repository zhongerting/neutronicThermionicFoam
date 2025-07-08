/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtMultiRegionFoam

Description
    Solver for steady or transient fluid flow and solid heat conduction, with
    conjugate heat transfer between regions, buoyancy effects, turbulence,
    reactions and radiation modelling.

\*---------------------------------------------------------------------------*/
#include "thermionicConvertion.H"
#include "TFEsingle.H"
#include "extractSingleIntFromString.H"
#include "sort_with_correspondence.H"
#include "processVectors.H"
#include "thermionicEmission.H"
#include "singleThermionicEnergyConversion.H"
#include "circuitTECs.H"

#include "fvCFD.H"
#include "dynamicMomentumTransportModel.H"
#include "fluidReactionThermophysicalTransportModel.H"
#include "fluidReactionThermo.H"
#include "combustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "coordinateSystem.H"
#include "pimpleMultiRegionControl.H"
#include "pressureReference.H"
#include "hydrostaticInitialisation.H"

#include "neutronODE.H"
#include "powerControl.H"

#include <omp.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // 验证并行是否生效
    #pragma omp single
    {
        Info << "Hello, omp_get_max_threads = " << omp_get_max_threads() << endl;
    }

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    // 创建常微分方程求解类
    argList::validArgs.append("ODESolver");

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    pimpleMultiRegionControl pimples(fluidRegions, solidRegions);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createFluidPressureControls.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    // 读入输入功率求解方法，并进行初始化
    #include "powerInput.H"

    while (pimples.run(runTime))
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // 功率计算, 计算当前时间功率
        #include "powerCalc.H"

        // --- PIMPLE loop
        while (pimples.loop())
        {
            // #pragma omp parallel for
            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "solveFluid.H"
            }

            // #pragma omp parallel for
            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "solveSolid.H"
            }

            // 计算电路相关参数
            #include "setEmitterAndCollector.H"
            #include "solveTFE.H"   // 如果不需要电路，可以注释掉，只用注释这一个！！！！
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

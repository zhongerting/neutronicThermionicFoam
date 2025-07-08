// 注册功率因子参数（创建场进行计算）
PtrList<volScalarField> powerFactorSolid(solidRegions.size());
forAll(solidRegions, i){
    Info << "    Adding powerFactorSolid\n" << endl;
    powerFactorSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "powerFactor",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensioned<scalar>(dimless, 1.0)
        )
    );
}

scalar time00 = runTime.value();
scalar neutronTimeWriteInterval = 1.0;
scalar neutronTimeWrite = time00;


// PtrList<solidThermo> thermoSolid(solidRegions.size());



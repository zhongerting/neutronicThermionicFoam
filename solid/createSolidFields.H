// Initialise solid field pointer lists
PtrList<solidThermo> thermoSolid(solidRegions.size());
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size());
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size());

// PtrList<volScalarField> innerFlagSolid(solidRegions.size());
// PtrList<volScalarField> outerFlagSolid(solidRegions.size());
// 边界温度，用于存储边界上轴向各个位置平均温度
PtrList<volScalarField> boundTemperatureSolid(solidRegions.size());
// 边界标记，读入变量，用于标记内/外边界
PtrList<volScalarField> boundFlagSolid(solidRegions.size());
// 附加热量，包含电子导热部分（专用于发射极和接收极）
PtrList<volScalarField> addHeatFluxSolid(solidRegions.size());
// PtrList<volScalarField> outerTemperatureSolid(solidRegions.size());

List<label> innerBoundIndex(solidRegions.size(), -1);
List<label> outerBoundIndex(solidRegions.size(), -1);

// Populate solid field pointer lists
// #pragma omp parallel for
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoSolid\n" << endl;
    thermoSolid.set(i, solidThermo::New(solidRegions[i]));

    Info<< "    Adding to fvModelsSolid\n" << endl;
    fvModelsSolid.set
    (
        i,
        new Foam::fvModels(solidRegions[i])
    );

    Info<< "    Adding fvConstraintsSolid\n" << endl;
    fvConstraintsSolid.set
    (
        i,
        new fvConstraints(solidRegions[i])
    );

    Info<< "    Adding boundTemperatureSolid\n" << endl;
    boundTemperatureSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "boundTemperature",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensioned<scalar>(dimless, 743.)
        )
    );

    Info<< "    Adding boundFlagSolid\n" << endl;
    boundFlagSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "boundFlag",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    Info<< "    Adding addHeatFluxSolid\n" << endl;
    addHeatFluxSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "addHeatFlux",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar(dimless, 0.0)
        )
    );

    // Info<< "    Adding outerTemperatureSolid\n" << endl;
    // outerTemperatureSolid.set
    // (
    //     i,
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             "outerTemperature",
    //             runTime.timeName(),
    //             solidRegions[i],
    //             IOobject::MUST_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         solidRegions[i],
    //         dimensionedScalar(dimless, 0.0)
    //     )
    // );

    // Info<< "    Adding innerFlagSolid\n" << endl;
    // innerFlagSolid.set
    // (
    //     i,
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             "innerFlag",
    //             runTime.timeName(),
    //             solidRegions[i],
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         solidRegions[i],
    //         dimensionedScalar(dimless, 0.0)
    //     )
    // );

    // Info<< "    Adding outerFlagSolid\n" << endl;
    // outerFlagSolid.set
    // (
    //     i,
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             "outerFlag",
    //             runTime.timeName(),
    //             solidRegions[i],
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         solidRegions[i],
    //         dimensionedScalar(dimless, 0.0)
    //     )
    // );
}

Info << "all Fields are created" << nl;

// 标记内外表面索引
forAll(solidRegions, i)
{
    Info << "region name is " << solidRegions[i].name() << endl;
    const volScalarField& boundT = boundFlagSolid[i];
    forAll(boundT.boundaryField(), patchI)
    {
        Info << "patchI = " << patchI << tab;
        if (abs(boundT.boundaryField()[patchI][0] - -1.) < 1.e-4)
        {
            innerBoundIndex[i] = patchI;
            Info << "value = " << boundT.boundaryField()[patchI][0] << tab << "inner boundary name is " << solidRegions[i].boundary()[patchI].name() << endl;
        }
        else if (abs(boundT.boundaryField()[patchI][0] - 1.) < 1.e-4)
        {
            outerBoundIndex[i] = patchI;
            Info << "value = " << boundT.boundaryField()[patchI][0] << tab << "outer boundary name is " << solidRegions[i].boundary()[patchI].name() << endl;
        }
        else{
            Info << "value = " << boundT.boundaryField()[patchI][0] << tab << "normal boundary name is " << solidRegions[i].boundary()[patchI].name() << endl;
        }
    }
    Info << "innerBoundIndex is " << innerBoundIndex << tab << "outerBoundIndex is " << outerBoundIndex << endl;

    // const fvBoundaryMesh& boundTemp = solidRegions[i].boundary();
    // forAll (boundTemp, patchI)
    // {

    //     const fvPatch& patch = boundTemp[patchI];
        
    // }
    // forAll(solidRegions[i], patchI)
    // {
    //     const fvPatch& patch = solidRegions[i].boundary[patchI];
    //     Info << "Processing patch: " << patch.name() << endl;

    //     fixedValueFvPatchScalarField& fvPatchField = 
    //         refCast<fixedValueFvPatchScalarField>(patch.lookupPatchField<scalar>("boundTemperature"));
    //     scalar valueTemp = fvPatchField.value()[0];

    //     if (abs(valueTemp - -1.) < 1.e-5)
    //     {
    //         innerBoundIndex[i] = patchI;
    //     }
    //     else if (abs(valueTemp - 1.) < 1.e-5)
    //     {
    //         outerBoundIndex[i] = patchI;
    //     }
    // }
    // {
    //     // const fvPatch& patch = boundaryE[patchI];
    //     // std::string patchName = std::string(patch.name());
    //     // if ()
    //     // {
    //     //     innerBoundIndex = patchI;
    //     // }
    //     // else if ()
    //     // {
    //     //     outerBoundIndex = patchI;
    //     // }
    // }
}


std::vector<std::vector<std::vector<int>>> innerFaceId(int(solidRegions.size()));
std::vector<std::vector<std::vector<int>>> outerFaceId(int(solidRegions.size()));

// 标记轴向索引
forAll(solidRegions, i)
{
    Info << "region name is " << solidRegions[i].name() << endl;
    // 内表面轴向索引
    if (innerBoundIndex[i] >= 0.)
    {
        std::vector<int> faceZIndexInner;
        std::vector<double> faceZInner;
        const vectorField& CfInner = solidRegions[i].boundary()[innerBoundIndex[i]].Cf();
        forAll(CfInner, faceI){
            faceZIndexInner.push_back(faceI);
            faceZInner.push_back(CfInner[faceI].z());
        }
        innerFaceId[i] = processVectors(faceZInner, faceZIndexInner);
        // for(int j = 0; j < int(innerFaceId[i].size()); ++j){
        //     for(int k = 0; k < int(innerFaceId[i][j].size()); ++k){
        //         const label faceI = innerFaceId[i][j][k];
        //         Info << faceI << "\t" << CfInner[faceI].z() << nl;
        //     }
        //     Info << nl;
        // }
    }
    
    if (outerBoundIndex[i] >= 0)
    {
        // 外表面轴向索引
        std::vector<int> faceZIndexOuter;
        std::vector<double> faceZOuter;
        const vectorField& CfOuter = solidRegions[i].boundary()[outerBoundIndex[i]].Cf();
        forAll(CfOuter, faceI){
            faceZIndexOuter.push_back(faceI);
            faceZOuter.push_back(CfOuter[faceI].z());
        }
        outerFaceId[i] = processVectors(faceZOuter, faceZIndexOuter);
        // for(int j = 0; j < int(outerFaceId[i].size()); ++j){
        //     for(int k = 0; k < int(outerFaceId[i][j].size()); ++k){
        //         const label faceI = outerFaceId[i][j][k];
        //         Info << faceI << "\t" << CfOuter[faceI].z() << nl;
        //     }
        //     Info << nl;
        // }
    }
    
}

// FatalError << "Error test!" << exit(FatalError);

#include "createTFEFields.H"
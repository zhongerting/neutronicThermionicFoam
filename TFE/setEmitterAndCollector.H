// 获取发射极和接收极表面温度

// 第i个TEC装置
for(int i = 0; i < int(emitterIndex.size()); ++i)
{
    const volScalarField& TE = thermoSolid[emitterIndex[i]].T();
    // Info<< "name = " << TE.mesh().name() << nl;
    const scalarField& patchTE = TE.boundaryField()[emitterCoupledInterfaceIndex[i]];
    for(int j = 0; j < int(emitterFaceId[i].size()); ++j){
        double Temp = 0.;
        // Info << "j = " << j << nl;
        for(int k = 0; k < int(emitterFaceId[i][j].size()); ++k){
            // Info << "T = " << patchTE[emitterFaceId[i][j][k]] << nl;
            Temp += patchTE[emitterFaceId[i][j][k]];
        }
        // Info << nl;
        TElist[i][j] = Temp / emitterFaceId[i][j].size();
    }
    const volScalarField& TC = thermoSolid[collectorIndex[i]].T();
    // Info<< "name = " << TC.mesh().name() << nl;
    const scalarField& patchTC = TC.boundaryField()[collectorCoupledInterfaceIndex[i]];
    for(int j = 0; j < int(collectorFaceId[i].size()); ++j){
        scalar Temp = 0.;
        for(int k = 0; k < int(collectorFaceId[i][j].size()); ++k){
            Temp += patchTC[collectorFaceId[i][j][k]];
        }
        TClist[i][j] = Temp / collectorFaceId[i][j].size();
    }
    // for(int j = 0; j < int(collectorFaceId[i].size()); ++j){
    //     Info<< TElist[i][j] << "\t\t" << TClist[i][j] << nl;
    // }
    // const fvBoundaryMesh& boundary = mesh.boundary();
    // const fvPatch& patch = boundary[emitterCoupledInterfaceIndex[i]];
    // const pointField& points = patch.patch().localPoints();
    // for(int j = 0; j < int(emitterFaceId[i].size()); ++j){
    //     for(int k = 0; k < int(emitterFaceId[i][j].size()); ++k){
    //         Info<< "faceID = " << emitterFaceId[i][j][k] <<
    //         "\t\tz = " << points[emitterFaceId[i][j][k]].z()
    //         << "\t\tT = " << patchT[emitterFaceId[i][j][k]]
    //         <<nl;
    //     }
    //     Info << nl;
    // }

    Info << nl;
    
}

// FatalError << "Error test!" << exit(FatalError);

// 对每个TEC元件进行赋值
// 先进行并联赋值
// 再进行串联赋值
for (int i = 0; i < int(C1->TECsParallel.size()); ++i)
{
    // Info << "i = " << i << nl;
    for(int j = 0; j < int(TElist[i].size()); ++j)
    {
        C1->TECsParallel[i]->Temitter[j] = TElist[i][j];
        // Info << C1->TECsParallel[i]->Temitter[j] << tab;
    }
    // Info << nl;
    for(int j = 0; j < int(TClist[i].size()); ++j)
    {
        C1->TECsParallel[i]->Tcollector[j] = TClist[i][j];
        // Info << C1->TECsParallel[i]->Tcollector[j] << tab;
    }
    // Info << nl;
}
// Foam::sleep(60);
for(int i = 0; i < int(C1->TECs.size()); ++i){
    for(int j = 0; j < int(TElist[i].size()); ++j){
        C1->TECs[i]->Temitter[j] = TElist[i + int(C1->TECsParallel.size())][j];
        //Info << C1->TECs[i]->Temitter[j] << ", ";
    }
    //Info << endl;
    for(int j = 0; j < int(TClist[i].size()); ++j){
        C1->TECs[i]->Tcollector[j] = TClist[i + int(C1->TECsParallel.size())][j];
        //Info << C1->TECs[i]->Tcollector[j] << ", ";
    }
    //Info << endl;
}
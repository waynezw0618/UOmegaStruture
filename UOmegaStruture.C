/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the teRMS of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    vorticity

Description
    Calculates and writes the correlation between U and Vorticity.

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include <vector>
#define Nx 40//192
#define Ny0  25//24
#define Ny1 25//36
#define Ny2 25
#define Ny3 8
#define Nz  30//160
#define Nblocks 4


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int getGlobalID(std::vector<int> localID)
{
    if (localID[3]==0) {
        return (Nx*Ny0)*localID[2]+Nx*localID[1]+localID[0];
    }
    else if(localID[3]==1){
        return (Nx*Ny1)*localID[2]+Nx*localID[1]+localID[0]+Nx*Ny0*Nz;
    }
    else if(localID[3]==2){
        return (Nx*Ny2)*localID[2]+Nx*localID[1]+localID[0]+Nx*(Ny0+Ny1)*Nz;
    }
    else if(localID[3]==3){
        return (Nx*Ny3)*localID[2]+Nx*localID[1]+localID[0]+Nx*(Ny0+Ny1+Ny2)*Nz;
    }
    else{
        cout << "out of range, when calculate GlobalID for point i=["
           << localID[0] <<"],j=["<< localID[1]<<"], and k=["
           << localID[2]<<"], and blk_ID=["<< localID[3]<<"]"
           << std::endl;
        return -1;
    }
}

std::vector<int> getLocalID(int ID)
{
    std::vector<int> localID;

    if(ID<Nx*Ny0*Nz){
        localID.push_back((ID%(Nx*Ny0))%Nx);
        localID.push_back((ID%(Nx*Ny0))/Nx);
        localID.push_back( ID/(Nx*Ny0));
        localID.push_back(0);
    }
    else if(ID>Nx*Ny0*Nz && ID<Nx*(Ny0+Ny1)*Nz ){
        localID.push_back(((ID-Nx*Ny0*Nz)%(Nx*Ny1))%Nx);
        localID.push_back(((ID-Nx*Ny0*Nz)%(Nx*Ny1))/Nx + Ny0);
        localID.push_back( (ID-Nx*Ny0*Nz)/(Nx*Ny1));
        localID.push_back(1);
        
    }
    else if(ID>Nx*(Ny0+Ny1)*Nz && ID<Nx*(Ny0+Ny1+Ny2)*Nz){
        localID.push_back(((ID-Nx*(Ny0+Ny1)*Nz)%(Nx*Ny2))%Nx);
        localID.push_back(((ID-Nx*(Ny0+Ny1)*Nz)%(Nx*Ny2))/Nx + Ny0 + Ny1);
        localID.push_back( (ID-Nx*(Ny0+Ny1)*Nz)/(Nx*Ny2));
        localID.push_back(2);
    }
    else if(ID>Nx*(Ny0+Ny1+Ny2)*Nz && ID<Nx*(Ny0+Ny1+Ny2+Ny3)*Nz){
        localID.push_back(((ID-Nx*(Ny0+Ny1+Ny2)*Nz)%(Nx*Ny3))%Nx);
        localID.push_back(((ID-Nx*(Ny0+Ny1+Ny2)*Nz)%(Nx*Ny3))/Nx + Ny0 + Ny1 + Ny2);
        localID.push_back( (ID-Nx*(Ny0+Ny1+Ny2)*Nz)/(Nx*Ny3));
        localID.push_back(3);
    }
    else{
        cout << "out of range, when calculate LocalID for global" <<ID << std::endl;
    }
    
    return localID;
}

int getsamRefGloablID(std::vector<int> samPtLocalID,  std::vector<int> orgPtLocalID, std::vector<int> refLocalID)
{
    std::vector<int> samRefLocalID(4);
    //get idx in x
    samRefLocalID[0]  = samPtLocalID[0] + (refLocalID[0] - orgPtLocalID[0]);
    if (samRefLocalID[0]<0) {
        samRefLocalID[0] = samRefLocalID[0] + Nx;
    } else if(samRefLocalID[0] >Nx){
        samRefLocalID[0] = samRefLocalID[0] - Nx;
    }
    
    //get idx in z
    samRefLocalID[2]  = samPtLocalID[2] + (refLocalID[2] - orgPtLocalID[2]);
    if (samRefLocalID[2]<0) {
        samRefLocalID[2] = samRefLocalID[2] + Nz;
    } else if(samRefLocalID[2] >Nz){
        samRefLocalID[2] = samRefLocalID[2] - Nz;
    }
    
    if ( (samRefLocalID[0]<0 || samRefLocalID[0]>Nx )
         || (samRefLocalID[2]<0 || samRefLocalID[2]>Nz)
        ) {
        cout << "out of range, local point is not in the proper position" << std::endl;
    }
    
    samRefLocalID[1] = refLocalID[1];
    samRefLocalID[3] = refLocalID[3];
    
    return getGlobalID(samRefLocalID);
}

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");
    
    Info<<"reading the 6 basic data for analysis."<<endl;
    
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField UMeanMap
    (
        IOobject
        (
            "UMeanMap",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volSymmTensorField URMSMap
    (
        IOobject
        (
            "URMSMap",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField vorticity
    (
        IOobject
        (
            "vorticity",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField vorticityMeanMap
    (
        IOobject
        (
            "vorticityMeanMap",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volSymmTensorField vorticityRMSMap
    (
        IOobject
        (
            "vorticityRMSMap",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );


    volVectorField Uturb(U-UMeanMap);
    volVectorField Wturb(vorticity-vorticityMeanMap);


    
    Info << "Reading reference points" << endl;
    IOdictionary refPointsProperties
    (
        IOobject
        (
            "refPointsProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
         )
    );
    
    List<vector> refpts = List<vector>(refPointsProperties.lookup("refPoints"));
    
    forAll(refpts,ref_idx)
    {
        int gID_refpts = mesh.findCell(refpts[ref_idx]);
        std::vector<int> localID_refpts=getLocalID(gID_refpts);
        Info << "Calculate R for the point (" <<refpts[ref_idx]<<"),and it is localted at"
             << mesh.C()[gID_refpts] << nl
             << "it's global ID is " << gID_refpts << " and its local id is i=["
             << localID_refpts[0] <<"],j=["<< localID_refpts[1]<<"], and k=["
             << localID_refpts[2]<<"], and blk_ID=["<< localID_refpts[3]<<"]"
             << endl;
        
        volVectorField RUW(
                                IOobject
                                (
                                    "RUW"+name(ref_idx),
                                    runTime.timeName(),
                                    mesh,
                                    IOobject::NO_READ,
                                    IOobject::AUTO_WRITE
                                ),
                                mesh,
                                dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0, 0, 0), vector::zero)
                        );


        
        forAll(U,cellI)
        {
            std::vector<int> localID_pts=getLocalID(cellI);
            //loop over layer
            for (int i=0; i<Nx; i++) {
                for (int k=0; k<Nz; k++) {
                        std::vector<int> localID_samPt(4,0);
                        localID_samPt[0]=i;
                        localID_samPt[1]=localID_pts[1];
                        localID_samPt[2]=k;
                        localID_samPt[3]=localID_pts[3];
                    Info << "a"<< endl;
                    Info << "i=" << localID_samPt[0] << ",j="<<localID_samPt[1]<<", k=" <<localID_samPt[2]<<",nBK="<<localID_samPt[3]<<endl;
                        int gID_samPt    = getGlobalID(localID_samPt);
                    Info << "b" << endl;
                        int gID_samRefPt = getsamRefGloablID(localID_samPt,localID_pts,localID_refpts);
                        RUW.component(vector::X)()[cellI] = RUW.component(vector::X)()[cellI]
                                                            + Uturb.component(vector::X)()[gID_samRefPt]*Wturb.component(vector::X)()[gID_samPt];
                    // Foam::sqrt(URMSMap.component(tensor::XX)()[gID_samRefPt]*vorticityRMSMap.component(tensor::XX)()[gID_samPt]);
                        RUW.component(vector::Y)()[cellI] = RUW.component(vector::Y)()[cellI]
                                                            + Uturb.component(vector::X)()[gID_samRefPt]*Wturb.component(vector::Y)()[gID_samPt];
                                                           // Foam::sqrt(URMSMap.component(tensor::XX)()[gID_samRefPt]*vorticityRMSMap.component(tensor::YY)()[gID_samPt]);
                        RUW.component(vector::Z)()[cellI] = RUW.component(tensor::XX)()[cellI]
                                                            + Uturb.component(vector::X)()[gID_samRefPt]*Wturb.component(vector::Z)()[gID_samPt];
                                                           //Foam::sqrt(URMSMap.component(tensor::XX)()[gID_samRefPt]*vorticityRMSMap.component(tensor::ZZ)()[gID_samPt]);
                }
            }
            RUW[cellI]=RUW[cellI]/(Nx*Nz);
            
        }
        
        if (writeResults)
        {
            RUW.write();
        }
        
    }
    
    Info<< "\nEnd\n" << endl;
        
}


// ************************************************************************* //

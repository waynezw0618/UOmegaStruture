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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

    volVectorField URMSMap
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

    volVectorField vorticityRMSMap
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
    volScalarField Q
    (
        IOobject
        (
            "Q",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField QMeanMap
    (
        IOobject
        (
            "QMeanMap",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField QRMSMap
    (
        IOobject
        (
            "QRMSMap",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField UPrime(U-UMeanMap);
    volVectorField vorticityPrime(vorticity-vorticityMeanMap);
    volScalarField QPrime(Q-QMeanMap);


/*    timeSelector::addOptions();
    argList::validArgs.append('Yr');
    scalar Yr(args.additionalArgs([0]);
*/
    dimensionedScalar Xr("Xr",dimless,6.283145314); //reference Y coordinate
    dimensionedScalar Xcoor("Xcoor",dimless,0.0); //reference Y coordinate
    dimensionedScalar xTol("xTol",dimless,1e-3); //reference Y coordinate
    
    dimensionedScalar Yr("Zr",dimless,1.7445e-02); //reference Y coordinate
    dimensionedScalar Ycoor("Zcoor",dimless,0.0); //reference Y coordinate
    dimensionedScalar yTol("zTol",dimless,1e-3); //reference Y coordinate

    dimensionedScalar Zr("Yr",dimless,3.14); //reference Y coordinate
    dimensionedScalar Zcoor("Ycoor",dimless,0.0); //reference Y coordinate
    dimensionedScalar zTol("yTol",dimless,1e-3); //reference Y coordinate
    scalar k(0); //reference Y coordinate
/*    dimensionedScalar temp1("temp1",dimless,0); //reference Y coordinate
    dimensionedScalar temp2("temp2",dimless,0); //reference Y coordinate
    dimensionedScalar RUU11("RUU11",dimless,0); //self-correation coefficient of U and U
    dimensionedScalar RUU22("RUU22",dimless,0); //self-correation coefficient of V and V
    dimensionedScalar RUU33("RUU33",dimless,0); //self-correation coefficient of W and W
*/

    scalar temp1(0); //reference Y coordinate
    scalar temp2(0); //reference Y coordinate
    scalar RUU11(0); //self-correation coefficient of U and U
    scalar RUU22(0); //self-correation coefficient of V and V
    scalar RUU33(0); //self-correation coefficient of W and W
    scalar small_scalar(SMALL);

    volScalarField RUW13
    (
        IOobject
        (
            "RUW13",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimless,0)
    );
    volScalarField RUW23
    (
        IOobject
        (
            "RUW23",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
                mesh,
                dimensionedScalar("zero",dimless,0)
    );
    volScalarField RUW33
    (
        IOobject
        (
            "RUW33",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
                mesh,
                dimensionedScalar("zero",dimless,0)
    );

    volScalarField RU1Q
    (
        IOobject
        (
            "RU1Q",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
                mesh,
                dimensionedScalar("zero",dimless,0)
    );
    volScalarField RU2Q
    (
        IOobject
        (
            "RU2Q",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
                mesh,
                dimensionedScalar("zero",dimless,0)
    );

    volScalarField RU3Q
    (
        IOobject
        (
            "RU3Q",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
                mesh,
                dimensionedScalar("zero",dimless,0)
    );

    volScalarField small_
    (
        IOobject
        (
            "small",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
                mesh,
                dimensionedScalar("zero",dimensionSet(0, 0, -1, 0, 0),SMALL)
    );

    volScalarField small2_
    (
        IOobject
        (
            "small2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
                mesh,
                dimensionedScalar("zero",dimensionSet(0, 0, -2, 0, 0),SMALL)
    );
Info<<"loop starting !"<<endl;
    forAll(U, cellI)
    {
        Xcoor=mesh.C()[cellI].component(vector::X);
        Ycoor=mesh.C()[cellI].component(vector::Y);
        Zcoor=mesh.C([cellI].component(vector::Z);
        if(Foam::mag(Xcoor-Xr)<yTol){
            if (Foam::mag(Zcoor-Zr)<yTol) {
                //poit 1
                if (Foam::mag(Ycoor-Yr)<yTol) {
                    temp1=UPrime.component(vector::X)()[cellI];
                    temp2=URMSMap.component(vector::X)()[cellI];
                    //Info <<"before RUU11" << endl;
                    RUU11=RUU11+temp1*temp1/(temp2*temp2+small_scalar);
                    
                    temp1=UPrime.component(vector::Y)()[cellI];
                    temp2=URMSMap.component(vector::Y)()[cellI];
                    //Info <<"before RUU22" << endl;
                    RUU22=RUU22+temp1*temp1/(temp2*temp2+small_scalar);
                    
                    temp1=UPrime.component(vector::Z)()[cellI];
                    temp2=URMSMap.component(vector::Z)()[cellI];
                    //Info <<"before RUU33" << endl;
                    RUU33=RUU33+temp1*temp1/(temp2*temp2+small_scalar);
                    
                    temp1=UPrime.component(vector::X)()[cellI];
                    temp2=URMSMap.component(vector::X)()[cellI];
                    //Info <<"before RUW13" << endl;
                    RUW13=RUW13+temp1*vorticityPrime.component(vector::Z)()/(temp2*vorticityRMSMap.component(vector::Z)()+small_);
                    
                    temp1=UPrime.component(vector::Y)()[cellI];
                    temp2=URMSMap.component(vector::Y)()[cellI];
                    //Info <<"before RUW23"<< endl;
                    RUW23=RUW23+temp1*vorticityPrime.component(vector::Z)()/(temp2*vorticityRMSMap.component(vector::Z)()+small_);
                    
                    temp1=UPrime.component(vector::Z)()[cellI];
                    temp2=URMSMap.component(vector::Z)()[cellI];
                    //Info <<"before RUW33"<< endl;
                    RUW33=RUW33+temp1*vorticityPrime.component(vector::Z)()/(temp2*vorticityRMSMap.component(vector::Z)()+small_);
                    
                    //Info <<"before RUQ"<< endl;
                    RU1Q=RU1Q+UPrime.component(vector::X)()[cellI]*QPrime/(URMSMap.component(vector::X)()[cellI]*QRMSMap+small2_);
                    RU2Q=RU2Q+UPrime.component(vector::Y)()[cellI]*QPrime/(URMSMap.component(vector::Y)()[cellI]*QRMSMap+small2_);
                    RU3Q=RU3Q+UPrime.component(vector::Z)()[cellI]*QPrime/(URMSMap.component(vector::Z)()[cellI]*QRMSMap+small2_);
                    
                    k=k+1;
                    Info << "the " << k<<"th iteration in the loop" << endl;
                }
            }
        }
    }
    RUU11=RUU11/k;
    RUU22=RUU22/k;
    RUU33=RUU33/k;

    RUW13=RUW13/k;
    RUW23=RUW23/k;
    RUW33=RUW33/k;
    
    RU1Q=RU1Q/k;
    RU2Q=RU2Q/k;
    RU3Q=RU3Q/k;

            Info<<"RUU11 is "<<RUU11<<endl;
            Info<<"RUU22 is "<<RUU22<<endl;
            Info<<"RUU33 is "<<RUU33<<endl;

    Info<<"Formulating the Uprime field"<<endl;

        if (writeResults)
        {
            //RUU11.write();
            //RUU22.write();
            //RUU33.write();
            Info<<"RUU11 is "<<RUU11<<endl;
            Info<<"RUU22 is "<<RUU22<<endl;
            Info<<"RUU33 is "<<RUU33<<endl;

            RUW13.write();
            RUW23.write();
            RUW33.write();

            RU1Q.write();
            RU2Q.write();
            RU3Q.write();
  
      }
        else
        {
            Info<< "    No U" << endl;
        }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //

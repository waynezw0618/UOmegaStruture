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

/*    volVectorField UMeanMap
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


    volVectorField UPrime(U-UMeanMap);
    volVectorField vorticityPrime(vorticity-vorticityMeanMap);

*/
    forAll(U, cellI)
    {
       
        Info<<"In cell ID ["<< cellI<<"], cell center is" << mesh.C()[cellI]<<endl;
   /*  Info<<"Formulating the Uprime field"<<endl;

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
    */
        Info<< "\nEnd\n" << endl;
        
}


// ************************************************************************* //

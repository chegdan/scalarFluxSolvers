/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    scalarTransportFoam

Description
    Solves a transport equation for a passive scalar

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
  
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
/*#   include "readTimeControls.H"//added
#   include "CourantNo.H"//added

//Info<<"Setting initial delta time"<<endl;
#   include "setInitialDeltaT.H"//added--only need to set timestep once
#   include "showCoNum.H"//output the Courant number after timestep change
*/
#   include "readSIMPLEControls.H"//added--reads in tSchmidt to see if turbulent schmidt relation should be used

//treat relationship to R as a turbulent diffusivity

    volSymmTensorField Dt = ((k/(Sct*epsilon)))*R;

//correct negative values in Dt to positive
forAll(Dt, cellI){

	Dt[cellI].xx() = ( Dt[cellI].xx() < 0 ) ? -Dt[cellI].xx() : Dt[cellI].xx() ;
	Dt[cellI].xy() = ( Dt[cellI].xy() < 0 ) ? -Dt[cellI].xy() : Dt[cellI].xy() ;
	Dt[cellI].xz() = ( Dt[cellI].xz() < 0 ) ? -Dt[cellI].xz() : Dt[cellI].xz() ;
	Dt[cellI].yy() = ( Dt[cellI].yy() < 0 ) ? -Dt[cellI].yy() : Dt[cellI].yy() ;
	Dt[cellI].yz() = ( Dt[cellI].yz() < 0 ) ? -Dt[cellI].yz() : Dt[cellI].yz() ;
	Dt[cellI].zz() = ( Dt[cellI].zz() < 0 ) ? -Dt[cellI].zz() : Dt[cellI].zz() ;

}

    //volVectorField gradC = fvc::grad(C);

    //Dt.write();//must write the Dturbulent field if changed by ScNo.H
    phi.write();//write the phi field in initial directory


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {

	fvScalarMatrix CEqn
	(
           fvm::div(phi, C)
	 + fvm::SuSp(-fvc::div(phi), C)//added for boundedness from post (http://www.cfd-online.com/Forums/openfoam/64602-origin-fvm-sp-fvc-div-phi_-epsilon_-kepsilon-eqn.html)
	 - fvm::laplacian(D, C)
         - fvm::laplacian(Dt, C)
	);

	CEqn.relax();
	
	solve(CEqn);






/*            solve
            (
                fvm::ddt(C)
              + fvm::div(phi, C)
	      + fvm::SuSp(-fvc::div(phi), C)//added for boundedness from post (http://www.cfd-online.com/Forums/openfoam/64602-origin-fvm-sp-fvc-div-phi_-epsilon_-kepsilon-eqn.html)
	      - fvm::laplacian(D, C)
              - fvm::laplacian(Dt, C)

            );
*/
        }

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

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
    Solves a transient transport equation for a passive scalar in a laminar flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
  
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readTimeControls.H"//added
#   include "CourantNo.H"//added
#   include "setInitialDeltaT.H"//added--only need to set timestep once
#   include "showCoNum.H"//output the Courant number after timestep change

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

	//Runge-Kutta Constant Declaration and initialization
	volScalarField k1 = C;
	volScalarField k2 = C;
	volScalarField k3 = C;
	volScalarField k4 = C;

	volScalarField C2 = C;
	volScalarField C3 = C;
	volScalarField C4 = C;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

	k1 == runTime.deltaT()*( fvc::laplacian(D, C) - fvc::div(phi, C) + fvc::div(phi)*C);
	
	C2 == C + 0.5 * k1;

	k2 == runTime.deltaT()*( fvc::laplacian(D, C2) - fvc::div(phi, C2) + fvc::div(phi)*C2);

	C3 == C + 0.5 * k2;

	k3 == runTime.deltaT()*( fvc::laplacian(D, C3) - fvc::div(phi, C3) + fvc::div(phi)*C3);

	C4 == C + k3;

	k4 == runTime.deltaT()*( fvc::laplacian(D, C4) - fvc::div(phi, C4) + fvc::div(phi)*C4);

	C == C + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

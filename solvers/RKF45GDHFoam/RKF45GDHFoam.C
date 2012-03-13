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
    Explicitly solves a transient transport equation for a passive scalar in a turbulent flow

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

    Dt = nut/Sct;

    Dt.write();//must write the Dturbulent field if changed by ScNo.H
    phi.write();//write the phi field in initial directory

    scalar eps = 1;
    dictionary simple = mesh.solutionDict().subDict("SIMPLE");
    simple.readIfPresent("eps", eps);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating implicit turbulent scalar transport\n" << endl;

	//Runge-Kutta Constant Declaration and initialization
	volScalarField k1 = C;
	volScalarField k2 = C;
	volScalarField k3 = C;
	volScalarField k4 = C;
	volScalarField k5 = C;
	volScalarField k6 = C;

	volScalarField C2 = C;
	volScalarField C3 = C;
	volScalarField C4 = C;
	volScalarField C5 = C;
	volScalarField C6 = C;

	volScalarField Y = C;
	volScalarField Z = C;
	
	volScalarField diff = C;

	dimensionedScalar s 
	( 
	"s", 
	dimensionSet(0,0,0,0,0,0,0), 
	scalar(1.0) 
	); 

	dimensionedScalar h = runTime.deltaT();


	dimensionedScalar one 
	( 
	"one", 
	dimensionSet(0,0,-0.25,0,0,0,0), 
	scalar(1.0) 
	); 


    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

	k1 == h*( fvc::laplacian(D, C) + fvc::laplacian(Dt, C) - fvc::div(phi, C) + fvc::div(phi)*C);
	
	C2 == C + 0.25 * k1;

	k2 == h*( fvc::laplacian(D, C2) + fvc::laplacian(Dt, C2) - fvc::div(phi, C2) + fvc::div(phi)*C2);

	C3 == C + (3.0/32.0) * k1 + (9.0/32.0) * k2;

	k3 == h*( fvc::laplacian(D, C3) + fvc::laplacian(Dt, C3) - fvc::div(phi, C3) + fvc::div(phi)*C3);

	C4 == C + (1932.0/2197.0) * k1 + (-7200.0/2197.0) * k2 + (7297.0/2197.0) * k3 ;

	k4 == h*( fvc::laplacian(D, C4) + fvc::laplacian(Dt, C4) - fvc::div(phi, C4) + fvc::div(phi)*C4);

	C5 == C + (439.0/216.0) * k1 + (-8.0) * k2 + (3600.0/513.0) * k3 + (-845.0/4104.0) * k4;
 
	k5 == h*( fvc::laplacian(D, C5) + fvc::laplacian(Dt, C5) - fvc::div(phi, C5) + fvc::div(phi)*C5);
	
	C6 == C + (-8.0/27.0) * k1 + (2.0) * k2 + (-3544.0/2565.0) * k3 + (1859.0/4104.0) * k4 + (-11.0/40.0)*k5;

	k6 == h*( fvc::laplacian(D, C6) + fvc::laplacian(Dt, C6) - fvc::div(phi, C6) + fvc::div(phi)*C6);

	Y == C + (25.0/216.0)*k1 + (1408.0/2565.0) *k3 + (2197.0/4104.0)  *k4 + (-1.0/5.0) *k5;
	Z == C + (16.0/135.0)*k1 + (6656.0/12825.0)*k3 + (28561.0/56430.0)*k4 + (-9.0/50.0)*k5 + (2.0/55.0)*k6;

	
	diff == sqrt(pow(Z - Y,2));
	//Info << "max: " << gMax(diff.internalField()) << endl;

	//dimensionedScalar maxDiff = gMax((Z-Y));
	s = one*pow(h/(2.0*gMax(diff.internalField())),0.25);

	C == Y;

	h = s*h;	

	runTime.deltaT() = h;

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
	    << " Step size = "<< h.value() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

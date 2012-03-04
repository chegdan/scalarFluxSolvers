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
    Solves a transient transport equation for a passive scalar in a turbulent flow

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating implicit scalar transport\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

	C == runTime.deltaT()*( fvc::laplacian(D, C) + fvc::laplacian(Dt, C) - fvc::div(phi, C) + fvc::div(phi)*C) + C;

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

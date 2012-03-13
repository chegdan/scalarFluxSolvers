
if (mesh.nInternalFaces())
{
    surfaceScalarField SfUfbyDelta = 
        mesh.surfaceInterpolation::deltaCoeffs()*mag(phi);

    CoNum = max(SfUfbyDelta/mesh.magSf())
        .value()*runTime.deltaT().value();

    meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
        .value()*runTime.deltaT().value();

    velMag = max(mag(phi)/mesh.magSf()).value();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum 
    << " velocity magnitude: " << velMag 
    << " timestep size "<<runTime.deltaT().value()<<" [sec] "<< endl;

// ************************************************************************* //

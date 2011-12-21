// check convergence

if (maxResidual < convergenceCriterion)
{
    Info<< "reached convergence criterion: " << convergenceCriterion << endl;
    runTime.writeAndEnd();
    U.write();
    phi.write();
    Dt.write();
    nut.write();
    Info<< "latestTime = " << runTime.timeName() << endl;
}


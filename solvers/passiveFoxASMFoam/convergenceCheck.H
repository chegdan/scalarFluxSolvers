// check convergence

if (maxResidual < convergenceCriterion)
{
    Info<< "reached convergence criterion: " << convergenceCriterion << endl;
    runTime.writeAndEnd();

    U.write();
    phi.write();
    k.write();
    R.write();
    epsilon.write();

    Info<< "latestTime = " << runTime.timeName() << endl;
}


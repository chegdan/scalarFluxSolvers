// initialize values for convergence checks

    scalar eqnResidual = 1, maxResidual = 0;
    scalar convergenceCriterion = 0;
    scalar eps = 1;

    simple.readIfPresent("convergence", convergenceCriterion);
    simple.readIfPresent("eps", eps);


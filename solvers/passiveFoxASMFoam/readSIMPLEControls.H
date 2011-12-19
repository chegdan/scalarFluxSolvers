    dictionary simple = mesh.solutionDict().subDict("SIMPLE");

    int nNonOrthCorr = 0;
    if (simple.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(simple.lookup("nNonOrthogonalCorrectors"));
    }

    int tSchmidt = 0;
    if (simple.found("tSchmidt"))
    {
	tSchmidt = readInt(simple.lookup("tSchmidt"));
    }

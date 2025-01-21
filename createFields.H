Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

dimensionedScalar rho		// MHD part
(
    "rho",
    dimDensity,
    transportProperties.lookup("rho")
);

dimensionedScalar sigma		// MHD part
(
    "sigma",
    transportProperties.lookup("sigma")
);

dimensionedVector B		// MHD part
(
    "B",
    transportProperties.lookup("B")
);

volScalarField PotE
(
    IOobject
    (
        "PotE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField lorentz
(
    IOobject
    (
        "lorentz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedVector(dimensionSet(1, -2, -2, 0, 0 ,0 ,0), Zero)
);

volVectorField jcenter
(
    IOobject
    (
        "jcenter",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedVector(dimensionSet(0, -2, 0, 0, 0 ,1 ,0), Zero)
);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

label PotERefCell = 0;
scalar PotERefValue = 0.0;
setRefCell(PotE, mesh.solutionDict().subDict("PotE"), PotERefCell, PotERefValue);
mesh.setFluxRequired(PotE.name());

const dictionary& PotEpiso = mesh.solutionDict().subDict("PotE");
int PotEcorr = PotEpiso.lookupOrDefault<int>("PotEcorrs", 1);
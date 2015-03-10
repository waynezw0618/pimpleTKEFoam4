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

Info<< "Reading field UMeanMap\n" << endl;
volVectorField UMeanMap
(
    IOobject
    (
        "UMeanMap",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField Uturb(U-UMeanMap);

Info<< "Reading field URMSMap\n" << endl;
volVectorField URMSMap
(
    IOobject
    (
        "URMSMap",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
     ),
    mesh
 );
Info<< "Calculating vorticity \n" << endl;

volVectorField vorticity
(
   IOobject
   (
        "vorticity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
   ),
   fvc::curl(U)
);

Info<< "Reading field vorticityMeanMap\n" << endl;

volVectorField vorticityMeanMap
(
    IOobject
    (
        "vorticityMeanMap",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField Wturb(vorticity-vorticityMeanMap);


Info<< "Reading field vorticityRMSMap\n" << endl;
volVectorField vorticityRMSMap
(
    IOobject
    (
        "vorticityRMSMap",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
     ),
    mesh
);

#include "TKE.H"
#include "Enstrophy.H"

volScalarField Lambda2
(
    IOobject
    (
        "Lambda2",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -2, 0, 0, 0, 0), 0.0)
);


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);


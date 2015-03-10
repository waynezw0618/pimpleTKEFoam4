Info << "creating variable for Enstrophy analysis" <<endl;

volTensorField ShearStress
(
    IOobject
    (
        "ShearStress",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    0.5*(fvc::grad(U)+fvc::grad(U)().T())
);

volTensorField WpWp
(
    IOobject
    (
        "WpWP",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Wturb*Wturb
);

volScalarField Pens21
(
    IOobject
    (
        "Pens21",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -2, 0, 0, 0, 0), 0.0)
);

volScalarField Pens22
(
    IOobject
    (
        "Pens22",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -2, 0, 0, 0, 0), 0.0)
 );

volScalarField Pens23
(
    IOobject
    (
        "Pens23",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -2, 0, 0, 0, 0), 0.0)
);

volScalarField Pens3
(
    IOobject
    (
        "Pens3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("zero", dimensionSet(0, 0, -3, 0, 0, 0, 0), 0.0)
);

volTensorField Pens4
(
    IOobject
    (
        "Pens4",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb*Wturb
);

volVectorField Tens
(
    IOobject
    (
        "Tens",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb*magSqr(Wturb)
);

volScalarField DisEnstropy
(
    IOobject
    (
        "DisEnstropy",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    fvc::grad(Wturb)&&fvc::grad(Wturb)().T()
);

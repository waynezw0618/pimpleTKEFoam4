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

volScalarField Pens31
(
    IOobject
    (
        "Pens31",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("zero", dimensionSet(0, 0, -3, 0, 0, 0, 0), 0.0)
);

volScalarField Pens32
(
    IOobject
    (
        "Pens32",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -3, 0, 0, 0, 0), 0.0)
);

volScalarField Pens33
(
    IOobject
    (
        "Pens33",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -3, 0, 0, 0, 0), 0.0)
);

volTensorField Pens4  //term (4)
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

volScalarField Tens1   //term 5 (1)
(
    IOobject
    (
        "Tens1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Foam::magSqr(Wturb.component(vector::X))*Uturb.component(vector::Y)
);

volScalarField Tens2    //term 5 (2)
(
    IOobject
    (
        "Tens2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Foam::magSqr(Wturb.component(vector::Y))*Uturb.component(vector::Y)
);

volScalarField Tens3 //term 5 (3)
(
    IOobject
    (
        "Tens3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Foam::magSqr(Wturb.component(vector::Z))*Uturb.component(vector::Y)
 );

volScalarField DisEnstropy1
(
    IOobject
    (
        "DisEnstropy1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Foam::magSqr(fvc::grad(Wturb.component(vector::X)))
);

volScalarField DisEnstropy2
(
    IOobject
    (
        "DisEnstropy2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Foam::magSqr(fvc::grad(Wturb.component(vector::Y)))
);


volScalarField DisEnstropy3
(
    IOobject
    (
        "DisEnstropy3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
  ),
    Foam::magSqr(fvc::grad(Wturb.component(vector::Z)))
);



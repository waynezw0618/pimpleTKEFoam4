Info << "creating variable for TKE analysis" <<endl;

volVectorField pU
(
    IOobject
    (
        "pU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p*U
);

volScalarField Uxxx
(
    IOobject
    (
        "Uxxx",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::X)()*U.component(vector::X)()*U.component(vector::X)()
);

volScalarField Uxxy
(
    IOobject
    (
        "Uxxy",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    U.component(vector::X)()*U.component(vector::X)()*U.component(vector::Y)()
);

volScalarField Uxxz
(
    IOobject
    (
        "Uxxz",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::X)()*U.component(vector::X)()*U.component(vector::Z)()
);

volScalarField Uyyy
(
    IOobject
    (
        "Uyyy",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    U.component(vector::Y)()*U.component(vector::Y)()*U.component(vector::Y)()
);

volScalarField Uyyx
(
    IOobject
    (
        "Uyyx",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    U.component(vector::Y)()*U.component(vector::Y)()*U.component(vector::X)()
 );

volScalarField Uyyz
(
    IOobject
    (
        "Uyyz",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    U.component(vector::Y)()*U.component(vector::Y)()*U.component(vector::Z)()
);

volScalarField Uzzz
(
    IOobject
    (
        "Uzzz",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    U.component(vector::Z)()*U.component(vector::Z)()*U.component(vector::Z)()
);

volScalarField Uzzx
(
    IOobject
    (
        "Uzzx",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::Z)()*U.component(vector::Z)()*U.component(vector::X)()
);

volScalarField Uzzy
(
    IOobject
    (
        "Uzzy",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    U.component(vector::Z)()*U.component(vector::Z)()*U.component(vector::Y)()
);

//epsilon
volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(Uturb)&&fvc::grad(Uturb)().T()
);

volScalarField Trans
(
    IOobject
    (
        "Trans",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb.component(vector::X)()*Uturb.component(vector::X)()*Uturb.component(vector::Y)()\
    +Uturb.component(vector::Y)()*Uturb.component(vector::Y)()*Uturb.component(vector::Y)()\
    +Uturb.component(vector::Z)()*Uturb.component(vector::Z)()*Uturb.component(vector::Y)()
);

//epsilon
volScalarField epsilonXX
(
    IOobject
    (
        "epsilonXX",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::XX) \
    +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::XY) \
    +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::XZ)
);

volScalarField epsilonXY
(
    IOobject
    (
        "epsilonXY",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::YX) \
    +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::YY) \
    +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::YZ)
);

volScalarField epsilonXZ
(
    IOobject
    (
        "epsilonXZ",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::ZX) \
    +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::ZY) \
    +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::ZZ)
 );

volScalarField epsilonYY
(
    IOobject
    (
        "epsilonYY",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    fvc::grad(Uturb)().component(tensor::YX)*fvc::grad(Uturb)().component(tensor::YX) \
    +fvc::grad(Uturb)().component(tensor::YY)*fvc::grad(Uturb)().component(tensor::YY) \
    +fvc::grad(Uturb)().component(tensor::YZ)*fvc::grad(Uturb)().component(tensor::YZ)
);

volScalarField epsilonYZ
(
    IOobject
    (
        "epsilonYZ",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    fvc::grad(Uturb)().component(tensor::YX)*fvc::grad(Uturb)().component(tensor::ZX) \
    +fvc::grad(Uturb)().component(tensor::YY)*fvc::grad(Uturb)().component(tensor::ZY) \
    +fvc::grad(Uturb)().component(tensor::YZ)*fvc::grad(Uturb)().component(tensor::ZZ)
);

volScalarField epsilonZZ
(
    IOobject
    (
        "epsilonZZ",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    fvc::grad(Uturb)().component(tensor::ZX)*fvc::grad(Uturb)().component(tensor::ZX) \
    +fvc::grad(Uturb)().component(tensor::ZY)*fvc::grad(Uturb)().component(tensor::ZY) \
    +fvc::grad(Uturb)().component(tensor::ZZ)*fvc::grad(Uturb)().component(tensor::ZZ)
);



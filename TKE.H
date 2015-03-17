Info << "creating variable for TKE analysis" <<endl;

volVectorField pU   //Pressure diffusion term (3)
(
    IOobject
    (
        "pU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    pturb*Uturb
);

volTensorField UpUp  // for diffusion term (4)
(
    IOobject
    (
        "UpUp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb*Uturb
 );

volTensorField Trans  // three dimensioned for Reynolds stress analysis (2)
(
    IOobject
    (
        "Trans",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    UpUp*Uturb.component(vector::Y)
);

volTensorField PreStrain  // pressure strain term. (6)
(
    IOobject
    (
        "PreStrain",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.5*pturb*(fvc::grad(Uturb) + fvc::grad(Uturb).T())
 );

/*volScalarField Uxxx
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
*/
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

//epsilon(5)
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


